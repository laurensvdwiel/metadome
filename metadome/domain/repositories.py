import logging
import traceback
from collections import defaultdict

from sqlalchemy.orm.exc import NoResultFound, MultipleResultsFound
from sqlalchemy.sql.functions import func
from sqlalchemy.exc import ResourceClosedError as AlchemyResourceClosedError
from sqlalchemy.exc import OperationalError as AlchemyOperationalError
from Bio.Data.IUPACData import protein_letters_1to3

from metadome.database import db
from metadome.domain.models.mapping import Mapping
from metadome.domain.models.gene import Gene
from metadome.domain.models.meta_domain_mapping import MetaDomainMapping
from metadome.domain.models.meta_domain_position import MetaDomainPosition
from metadome.domain.models.protein import Protein
from metadome.domain.models.interpro import Interpro
from metadome.domain.error import RecoverableError

from flask import current_app

_log = logging.getLogger(__name__)

class RepositoryException(Exception):
    pass

class MalformedAARegionException(Exception):
    pass

class RepositoryCacheException(Exception):
    """Raised when the repository cache layer cannot be used safely."""
    pass

def _cleanup_failed_session(_session):
    try:
        _session.rollback()
    except Exception:
        _log.exception("Session rollback failed after DB operational error")

    try:
        bind = _session.get_bind()
        if bind is not None:
            invalidate = getattr(bind, "invalidate", None)
            if callable(invalidate):
                invalidate()
            else:
                dispose = getattr(bind, "dispose", None)
                if callable(dispose):
                    dispose()
    except Exception:
        _log.exception("Connection/engine cleanup failed after DB operational error")

class RepositoryCacheHelper:
    @staticmethod
    def _get_cache():
        cache = getattr(current_app, "cache", None)
        if cache is None:
            raise RepositoryCacheException(
                "Flask cache is not configured: current_app.cache is missing."
            )
        return cache

    @staticmethod
    def _resolve_timeout(timeout: int | None) -> int:
        if timeout is not None:
            return timeout

        resolved = (
            current_app.config.get("REPOSITORY_CACHE_TIMEOUT")
            or current_app.config.get("CACHE_DEFAULT_TIMEOUT")
        )
        if resolved is None:
            raise RepositoryCacheException(
                "No cache timeout configured. Set REPOSITORY_CACHE_TIMEOUT "
                "(preferred) or CACHE_DEFAULT_TIMEOUT."
            )
        return int(resolved)

    @staticmethod
    def cache_get(cache_key: str) -> tuple[bool, object | None]:
        """
        Returns (hit, value).

        Raises RepositoryCacheException on cache backend errors or misconfiguration.
        """
        cache = RepositoryCacheHelper._get_cache()
        try:
            has_method = getattr(cache, "has", None)
            if callable(has_method):
                if not cache.has(cache_key):
                    return False, None
                return True, cache.get(cache_key)

            # Fallback: cannot distinguish cached None from miss (treat None as miss)
            value = cache.get(cache_key)
            if value is None:
                return False, None
            return True, value
        except Exception as e:
            raise RepositoryCacheException(f"Cache GET failed for key={cache_key}") from e

    @staticmethod
    def cache_set(cache_key: str, value: object, timeout: int | None = None) -> None:
        """
        Cache value with timeout.

        Raises RepositoryCacheException on cache backend errors or misconfiguration.
        """
        cache = RepositoryCacheHelper._get_cache()
        resolved_timeout = RepositoryCacheHelper._resolve_timeout(timeout)

        try:
            cache.set(cache_key, value, timeout=resolved_timeout)
        except Exception as e:
            raise RepositoryCacheException(
                f"Cache SET failed for key={cache_key}, timeout={resolved_timeout}"
            ) from e

class MetaDomainRepository:

    @staticmethod
    def get_meta_domain_annotation_for_variants(_variant_dicts):
        """Retrieve mapping and metadomain information for a batch of variants."""
        _session = db._make_scoped_session(options={})

        def _normalize_gene_id(gene_id):
            return gene_id.split('.', 1)[0] if gene_id else gene_id

        def _normalize_chr(chromosome):
            return chromosome.lower() if chromosome else chromosome

        def _reverse_complement(base):
            mapping = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
            return mapping.get(base.upper()) if base else None

        try:
            normalized_variants = []
            positions = set()
            chromosomes = set()
            genome_build_prefixes = set()

            for variant in _variant_dicts:
                normalized_gene_id = _normalize_gene_id(variant['gene_id'])
                normalized_variant = {
                    'chr': variant['chr'],
                    'pos': int(variant['pos']),
                    'ref': variant.get('ref'),
                    'gene_id': normalized_gene_id,
                    'genome_build': variant['genome_build'],
                }
                normalized_variants.append(normalized_variant)
                positions.add(normalized_variant['pos'])
                chromosomes.add(_normalize_chr(normalized_variant['chr']))
                genome_build_prefixes.add(normalized_variant['genome_build'][:6].upper())

            if len(normalized_variants) == 0:
                return {}

            results = _session.query(
                Mapping.chromosome,
                Mapping.chromosome_position,
                Mapping.base_pair,
                Gene.gencode_gene_id,
                Gene.genome_build,
                Gene.strand,
                Mapping.id,
                MetaDomainPosition.ext_db_id,
                MetaDomainPosition.consensus_position
            ) \
                .join(Gene, Gene.id == Mapping.gene_id) \
                .outerjoin(MetaDomainMapping, MetaDomainMapping.mapping_id == Mapping.id) \
                .outerjoin(MetaDomainPosition, MetaDomainPosition.id == MetaDomainMapping.meta_domain_position_id) \
                .filter(
                    Mapping.chromosome_position.in_(positions),
                    func.lower(Mapping.chromosome).in_(chromosomes),
                    func.upper(func.substr(Gene.genome_build, 1, 6)).in_(genome_build_prefixes)
                ) \
                .all()

            grouped_results = {}
            for variant in normalized_variants:
                variant_key = (
                    variant['chr'],
                    variant['pos'],
                    variant['ref'],
                    variant['gene_id'],
                    variant['genome_build'],
                )
                grouped_results[variant_key] = {
                    'MetaDomainPositions': '',
                    'MetaDomainStatus': 'no_mapping',
                    'RefMatchStatus': 'not_checked',
                }

            for variant in normalized_variants:
                variant_key = (
                    variant['chr'],
                    variant['pos'],
                    variant['ref'],
                    variant['gene_id'],
                    variant['genome_build'],
                )

                matching_rows = []
                for row in results:
                    db_gene_id = row.gencode_gene_id or ''
                    db_genome_build = row.genome_build or ''

                    if (row.chromosome or '').lower() != variant['chr'].lower():
                        continue
                    if row.chromosome_position != variant['pos']:
                        continue
                    if not db_gene_id.startswith(variant['gene_id']):
                        continue
                    if not db_genome_build.upper().startswith(variant['genome_build'][:6].upper()):
                        continue

                    matching_rows.append(row)

                if len(matching_rows) == 0:
                    continue

                grouped_results[variant_key]['MetaDomainStatus'] = 'mapping_no_metadomain'

                ref_statuses = set()
                md_positions = set()

                for row in matching_rows:
                    if variant['ref'] and row.base_pair:
                        if variant['ref'].upper() == row.base_pair.upper():
                            ref_statuses.add('direct_match')
                        elif _reverse_complement(row.base_pair) == variant['ref'].upper():
                            ref_statuses.add('reverse_complement_match')
                        else:
                            ref_statuses.add('mismatch')

                    if row.ext_db_id is not None and row.consensus_position is not None:
                        md_positions.add((row.ext_db_id, row.consensus_position))

                if len(md_positions) > 0:
                    grouped_results[variant_key]['MetaDomainStatus'] = 'metadomain_found'
                    grouped_results[variant_key]['MetaDomainPositions'] = ';'.join(
                        f'{ext_db_id}:{consensus_position}'
                        for ext_db_id, consensus_position in sorted(md_positions)
                    )

                if 'direct_match' in ref_statuses:
                    grouped_results[variant_key]['RefMatchStatus'] = 'direct_match'
                elif 'reverse_complement_match' in ref_statuses:
                    grouped_results[variant_key]['RefMatchStatus'] = 'reverse_complement_match'
                elif 'mismatch' in ref_statuses:
                    grouped_results[variant_key]['RefMatchStatus'] = 'mismatch'

            return grouped_results
        except (AlchemyResourceClosedError, AlchemyOperationalError) as e:
            _cleanup_failed_session(_session)
            raise RecoverableError(str(e))
        except:
            _log.error(traceback.format_exc())
            raise
        finally:
            _session.remove()

    @staticmethod
    def retrieve_all_mappings_for_meta_domain(_ext_db_id, _genome_build):
        """Retrieves all Mappings, Genes, MetaDomainMapping, MetaDomainPosition, and Protein for
         a given metadomain ext_db_id and genome_build"""
        # Open as session
        _session = db._make_scoped_session(options={})

        try:
            results = _session.query(Mapping, Gene, MetaDomainMapping, MetaDomainPosition, Protein)\
                    .join(Gene, Mapping.gene_id == Gene.id)\
                    .join(Protein, Mapping.protein_id == Protein.id)\
                    .join(MetaDomainMapping, MetaDomainMapping.mapping_id == Mapping.id)\
                    .join(MetaDomainPosition, MetaDomainPosition.id == MetaDomainMapping.meta_domain_position_id)\
                    .filter(\
                        MetaDomainPosition.ext_db_id == _ext_db_id,\
                        Gene.genome_build.like(_genome_build + '%')\
                    ).all()

            return results
        except (AlchemyResourceClosedError, AlchemyOperationalError) as e:
            _cleanup_failed_session(_session)
            raise RecoverableError(str(e))
        except:
            _log.error(traceback.format_exc())
            raise
        finally:
            # Close this session, thus all items are cleared and memory usage is kept at a minimum
            _session.remove()

class GeneRepository:

    @staticmethod
    def retrieve_gene_names_for_multiple_transcript_ids(_transcript_ids):
        """Retrieves all gene names for a given set of gencode transcripts 
        based on multiple Gene objects as {gencode_transcription_id: gene_name}"""
        # Open as session
        _session = db._make_scoped_session(options={})

        try:
            _gene_name_per_gencode_transcription_id = {}
            for gene in _session.query(Gene).filter(Gene.gencode_transcription_id.in_(_transcript_ids)).all():
                _gene_name_per_gencode_transcription_id[gene.gencode_transcription_id] = gene.gene_name
            return _gene_name_per_gencode_transcription_id
        except (AlchemyResourceClosedError, AlchemyOperationalError) as e:
            _cleanup_failed_session(_session)
            raise RecoverableError(str(e))
        except:
            _log.error(traceback.format_exc())
            raise
        finally:
            # Close this session, thus all items are cleared and memory usage is kept at a minimum
            _session.remove()

    @staticmethod
    def retrieve_genes_by_ids(_gene_ids):
        """Retrieve multiple genes by IDs as a dict {id: Gene}"""
        _session = db._make_scoped_session(options={})

        try:
            if not _gene_ids:
                return {}
            genes = _session.query(Gene).filter(Gene.id.in_(_gene_ids)).all()
            return {g.id: g for g in genes}
        except (AlchemyResourceClosedError, AlchemyOperationalError) as e:
            _cleanup_failed_session(_session)
            raise RecoverableError(str(e))
        except:
            _log.error(traceback.format_exc())
            raise
        finally:
            _session.remove()

    @staticmethod
    def retrieve_transcript_id_for_multiple_protein_ids(_protein_ids):
        """Retrieves all gencode transcripts for multiple protein ids as {gene_id: gencode_transcription_id}"""
        # Open as session
        _session = db._make_scoped_session(options={})

        try:
            _gencode_transcription_id_per_gene_id = {}
            for gene in _session.query(Gene).filter(Gene.protein_id.in_(_protein_ids)).all():
                _gencode_transcription_id_per_gene_id[gene.id] = gene.gencode_transcription_id
            return _gencode_transcription_id_per_gene_id
        except (AlchemyResourceClosedError, AlchemyOperationalError) as e:
            _cleanup_failed_session(_session)
            raise RecoverableError(str(e))
        except:
            _log.error(traceback.format_exc())
            raise
        finally:
            # Close this session, thus all items are cleared and memory usage is kept at a minimum
            _session.remove()

    @staticmethod
    def retrieve_all_transcript_ids_with_mappings():
        """Retrieves all transcripts for which there are mappings"""
        # Open as session
        _session = db._make_scoped_session(options={})

        try :
            return [transcript for transcript in _session.query(Gene.gencode_transcription_id).filter(Gene.protein_id != None).all()]
        except (AlchemyResourceClosedError, AlchemyOperationalError) as e:
            _cleanup_failed_session(_session)
            raise RecoverableError(str(e))
        except:
            _log.error(traceback.format_exc())
            raise
        finally:
            # Close this session, thus all items are cleared and memory usage is kept at a minimum
            _session.remove()

    @staticmethod
    def retrieve_all_transcript_ids_with_mappings_for_genome_build(genome_build):
        """Retrieves all transcripts with mappings for a single genome build.

        genome_build is the full build string as stored in the database (e.g. 'GRCh38').
        """
        # Open as session
        _session = db._make_scoped_session(options={})

        try:
            return [
                transcript
                for transcript in _session.query(Gene.gencode_transcription_id)
                .filter(Gene.protein_id != None,
                        func.lower(Gene.genome_build) == genome_build.lower())
                .all()
            ]
        except (AlchemyResourceClosedError, AlchemyOperationalError) as e:
            _cleanup_failed_session(_session)
            raise RecoverableError(str(e))
        except:
            _log.error(traceback.format_exc())
            raise
        finally:
            # Close this session, thus all items are cleared and memory usage is kept at a minimum
            _session.remove()

    @staticmethod
    def retrieve_all_gene_names_from_db():
        """Retrieves all gene names present in the database"""
        cache_key = "GeneRepository:gene_names:distinct"

        try:
            hit, cached = RepositoryCacheHelper.cache_get(cache_key)
        except RepositoryCacheException:
            _log.exception("Repository cache failure in retrieve_all_gene_names_from_db (key=%s)", cache_key)
            hit, cached = (False, None)

        if hit:
            return cached

        # Open as session
        _session = db._make_scoped_session(options={})

        try:
            result = [
                gene_name[0]
                for gene_name in _session.query(Gene.gene_name).distinct(Gene.gene_name).all()
            ]

            try:
                RepositoryCacheHelper.cache_set(cache_key, result)
            except RepositoryCacheException:
                _log.exception("Repository cache failure in retrieve_all_gene_names_from_db (set key=%s)", cache_key)

            return result
        except (AlchemyResourceClosedError, AlchemyOperationalError) as e:
            _cleanup_failed_session(_session)
            raise RecoverableError(str(e))
        except:
            _log.error(traceback.format_exc())
            raise
        finally:
            # Close this session, thus all items are cleared and memory usage is kept at a minimum
            _session.remove()

    @staticmethod
    def retrieve_all_genome_builds_from_db():
        """Retrieves all genome builds present in the database"""
        cache_key = "GeneRepository:genome_builds:distinct"

        try:
            hit, cached = RepositoryCacheHelper.cache_get(cache_key)
        except RepositoryCacheException:
            _log.exception("Repository cache failure in retrieve_all_genome_builds_from_db (key=%s)", cache_key)
            hit, cached = (False, None)

        if hit:
            return cached

        # Open as session
        _session = db._make_scoped_session(options={})

        try:
            result = [
                genome_build[0]
                for genome_build in _session.query(Gene.genome_build).distinct(Gene.genome_build).all()
            ]

            try:
                RepositoryCacheHelper.cache_set(cache_key, result)
            except RepositoryCacheException:
                _log.exception("Repository cache failure in retrieve_all_genome_builds_from_db (set key=%s)", cache_key)

            return result
        except (AlchemyResourceClosedError, AlchemyOperationalError) as e:
            _cleanup_failed_session(_session)
            raise RecoverableError(str(e))
        except:
            _log.error(traceback.format_exc())
            raise
        finally:
            # Close this session, thus all items are cleared and memory usage is kept at a minimum
            _session.remove()

    @staticmethod
    def retrieve_all_transcript_ids(genome_build, gene_name):
        """Retrieves all transcript ids for a gene name as plain serializable dicts"""
        cache_key = f"GeneRepository:transcript_ids:{genome_build.lower()}:{gene_name.lower()}"

        try:
            hit, cached = RepositoryCacheHelper.cache_get(cache_key)
        except RepositoryCacheException:
            _log.exception(
                "Repository cache failure in retrieve_all_transcript_ids (key=%s)",
                cache_key
            )
            hit, cached = (False, None)

        if hit:
            return cached

        # Open as session
        _session = db._make_scoped_session(options={})

        try:
            result = []
            for transcript in _session.query(Gene).filter(
                func.lower(Gene.gene_name) == gene_name.lower(),
                func.lower(Gene.genome_build) == genome_build.lower()
            ).all():
                result.append({
                    'gene_name': transcript.gene_name,
                    'sequence_length': transcript.sequence_length,
                    'gencode_transcription_id': transcript.gencode_transcription_id,
                    'refseq_transcript_id': transcript.refseq_transcript_id,
                    'mane_transcript_type': transcript.mane_transcript_type,
                    'protein_id': transcript.protein_id,
                })

            try:
                RepositoryCacheHelper.cache_set(cache_key, result)
            except RepositoryCacheException:
                _log.exception(
                    "Repository cache failure in retrieve_all_transcript_ids (set key=%s)",
                    cache_key
                )

            return result
        except (AlchemyResourceClosedError, AlchemyOperationalError) as e:
            _cleanup_failed_session(_session)
            raise RecoverableError(str(e))
        except:
            _log.error(traceback.format_exc())
            raise
        finally:
            # Close this session, thus all items are cleared and memory usage is kept at a minimum
            _session.remove()

    @staticmethod
    def retrieve_gene(transcription_id, genome_build):
        """Retrieves the gene object for a given transcript id and genome_build"""
        # Open as session
        _session = db._make_scoped_session(options={})
        try:
            return _session.query(Gene).filter(Gene.gencode_transcription_id == transcription_id, func.lower(Gene.genome_build) == genome_build.lower()).one()
        except (AlchemyResourceClosedError, AlchemyOperationalError) as e:
            _cleanup_failed_session(_session)
            raise RecoverableError(str(e))
        except MultipleResultsFound as e:
            error_message = "GeneRepository.retrieve_gene(transcription_id, genome_build): Multiple results found while expecting uniqueness for transcription_id '"+str(transcription_id)+"', genome_build='"+genome_build+"'. "+str(e)
            _log.error(error_message)
            raise RepositoryException(error_message)
        except NoResultFound as e:
            error_message = "GeneRepository.retrieve_gene(transcription_id, genome_build): Expected results but found none for transcription_id '"+str(transcription_id)+"', genome_build='"+genome_build+"'. "+str(e)
            _log.error(error_message)
            raise RepositoryException(error_message)
        except Exception as e:
            error_message = "GeneRepository.retrieve_gene(transcription_id, genome_build): Unexpected exception for transcription_id '"+str(transcription_id)+"', genome_build='"+genome_build+"'. "+str(e)
            _log.error(error_message)
            raise RepositoryException(error_message)
        finally:
            # Close this session, thus all items are cleared and memory usage is kept at a minimum
            _session.remove()

class InterproRepository:

    @staticmethod
    def get_all_Pfam_identifiers_suitable_for_metadomains():
        """Retrieves all pfam identifiers that occur at least two times in the database"""
        for domain_entry in Interpro.query.with_entities(Interpro.ext_db_id, func.count(Interpro.id)).filter(Interpro.ext_db_id.like('PF%')).group_by(Interpro.ext_db_id).having(func.count(Interpro.id)>2).distinct(Interpro.ext_db_id):
            yield domain_entry.ext_db_id


    @staticmethod
    def get_all_Pfam_identifiers():
        """Retrieves all pfam identifiers present in the database"""
        for domain_entry in Interpro.query.filter(Interpro.ext_db_id.like('PF%')).distinct(Interpro.ext_db_id):
            yield domain_entry.ext_db_id

    @staticmethod
    def get_domains_for_ext_domain_id(ext_domain_id):
        """Retrieves all interpro entries of the corresponding ext_db_id"""
        # Open as session
        _session = db._make_scoped_session(options={})

        try:
            return [interpro_domain for interpro_domain in _session.query(Interpro).filter(Interpro.ext_db_id == ext_domain_id).all()]
        except (AlchemyResourceClosedError, AlchemyOperationalError) as e:
            _cleanup_failed_session(_session)
            raise RecoverableError(str(e))
        except:
            _log.error(traceback.format_exc())
            raise
        finally:
            # Close this session, thus all items are cleared and memory usage is kept at a minimum
            _session.remove()

    @staticmethod
    def get_domains_for_protein(protein_id):
        """Retrieves all interpro entries for a given protein_id"""
        # Open as session
        _session = db._make_scoped_session(options={})

        try:
            return [interpro_domain for interpro_domain in _session.query(Interpro).filter(Interpro.protein_id == protein_id).all()]
        except (AlchemyResourceClosedError, AlchemyOperationalError) as e:
            _cleanup_failed_session(_session)
            raise RecoverableError(str(e))
        except:
            _log.error(traceback.format_exc())
            raise
        finally:
            # Close this session, thus all items are cleared and memory usage is kept at a minimum
            _session.remove()

class ProteinRepository:

    @staticmethod
    def retrieve_protein_ac_for_multiple_protein_ids(_protein_ids):
        """Retrieves all uniprot accession codes for multiple Protein objects as {protein_id: uniprot_ac}"""
        # Open as session
        _session = db._make_scoped_session(options={})

        try:
            _protein_ac_per_protein_id = {}
            for protein in _session.query(Protein).filter(Protein.id.in_(_protein_ids)).all():
                _protein_ac_per_protein_id[protein.id] = protein.uniprot_ac
            return _protein_ac_per_protein_id
        except (AlchemyResourceClosedError, AlchemyOperationalError) as e:
            _cleanup_failed_session(_session)
            raise RecoverableError(str(e))
        except:
            _log.error(traceback.format_exc())
            raise
        finally:
            # Close this session, thus all items are cleared and memory usage is kept at a minimum
            _session.remove()

    @staticmethod
    def retrieve_protein_id_for_multiple_protein_acs(_protein_acs):
        """Retrieves all protein ids for multiple Protein objects as {protein_ac: uniprot_id}"""
        # Open as session
        _session = db._make_scoped_session(options={})

        try:
            _protein_id_per_protein_ac = {}
            for protein in _session.query(Protein).filter(Protein.uniprot_ac.in_(_protein_acs)).all():
                _protein_id_per_protein_ac[protein.uniprot_ac] = protein.id
            return _protein_id_per_protein_ac
        except (AlchemyResourceClosedError, AlchemyOperationalError) as e:
            _cleanup_failed_session(_session)
            raise RecoverableError(str(e))
        except:
            _log.error(traceback.format_exc())
            raise
        finally:
            # Close this session, thus all items are cleared and memory usage is kept at a minimum
            _session.remove()

    @staticmethod
    def retrieve_protein(protein_id):
        """Retrieves the protein object for a given protein id"""
        # Open as session
        _session = db._make_scoped_session(options={})

        try:
            return _session.query(Protein).filter(Protein.id == protein_id).one()
        except (AlchemyResourceClosedError, AlchemyOperationalError) as e:
            _cleanup_failed_session(_session)
            raise RecoverableError(str(e))
        except MultipleResultsFound as e:
            _log.error("ProteinRepository.retrieve_protein(protein_id): Multiple results found while expecting uniqueness for protein_id '"+str(protein_id)+"'. "+str(e))
            return None
        except NoResultFound as  e:
            _log.error("ProteinRepository.retrieve_protein(protein_id): Expected results but found none for protein_id '"+str(protein_id)+"'. "+str(e))
            return None
        except:
            _log.error(traceback.format_exc())
            raise
        finally:
            # Close this session, thus all items are cleared and memory usage is kept at a minimum
            _session.remove()

    @staticmethod
    def retrieve_proteins_by_ids(_protein_ids):
        """Retrieve multiple proteins by IDs as a dict {id: Protein}"""
        _session = db._make_scoped_session(options={})

        try:
            if not _protein_ids:
                return {}
            proteins = _session.query(Protein).filter(Protein.id.in_(_protein_ids)).all()
            return {p.id: p for p in proteins}
        except (AlchemyResourceClosedError, AlchemyOperationalError) as e:
            _cleanup_failed_session(_session)
            raise RecoverableError(str(e))
        except:
            _log.error(traceback.format_exc())
            raise
        finally:
            _session.remove()

class MappingRepository:

    @staticmethod
    def _normalize_chromosome_input(chromosome: str) -> str:
        c = (chromosome or "").strip()
        if not c:
            return c
        if c.lower().startswith("chr"):
            suffix = c[3:]
        else:
            suffix = c
        return "chr" + suffix.upper()

    @staticmethod
    def _format_positions_to_ranges(chromosome: str, positions: list[int]) -> str:
        """Formats [1,2,3,7,8,10] -> chr1:1-3, 7-8, 10"""
        if not positions:
            return f"{chromosome}:"
        positions = sorted(set(positions))

        ranges = []
        start = positions[0]
        prev = positions[0]

        for p in positions[1:]:
            if p == prev + 1:
                prev = p
                continue
            if start == prev:
                ranges.append(str(start))
            else:
                ranges.append(f"{start}-{prev}")
            start = p
            prev = p

        if start == prev:
            ranges.append(str(start))
        else:
            ranges.append(f"{start}-{prev}")

        return f"{chromosome}:{', '.join(ranges)}"

    @staticmethod
    def lookup_mapping_hits_for_position(chromosome: str, position: int, genome_build: str | None = None):
        """
        Lightweight mapping-first lookup for position.
        Returns minimal mapping hit dicts for enrichment in a second step.
        """
        normalized_chr = MappingRepository._normalize_chromosome_input(chromosome)
        _session = db._make_scoped_session(options={})

        try:
            query = _session.query(
                Mapping.id,
                Mapping.gene_id,
                Mapping.protein_id,
                Mapping.chromosome,
                Mapping.chromosome_position,
                Mapping.uniprot_position,
                Mapping.uniprot_residue,
                Mapping.base_pair,
                Mapping.codon,
                Mapping.codon_base_pair_position,
                Mapping.strand
            ).filter(
                func.lower(Mapping.chromosome) == normalized_chr.lower(),
                Mapping.chromosome_position == position
            ).order_by(Mapping.id)

            if genome_build:
                query = query.join(Gene, Mapping.gene_id == Gene.id).filter(
                    func.lower(Gene.genome_build) == genome_build.lower())

            rows = query.all()

            return [
                {
                    "mapping_id": r.id,
                    "gene_id": r.gene_id,
                    "protein_id": r.protein_id,
                    "chromosome": r.chromosome,
                    "chromosome_position": r.chromosome_position,
                    "uniprot_position": r.uniprot_position,
                    "uniprot_residue": r.uniprot_residue,
                    "base_pair": r.base_pair or "",
                    "codon": r.codon or "",
                    "codon_position": (
                                r.codon_base_pair_position + 1) if r.codon_base_pair_position is not None else None,
                    "strand": r.strand.value if r.strand is not None else "",
                }
                for r in rows
                if r.uniprot_position is not None
            ]
        except (AlchemyResourceClosedError, AlchemyOperationalError) as e:
            _cleanup_failed_session(_session)
            raise RecoverableError(str(e))
        except:
            _log.error(traceback.format_exc())
            raise
        finally:
            _session.remove()

    @staticmethod
    def lookup_position_results(chromosome: str, position: int, genome_build: str | None = None):
        """
        Returns UI-ready rows for position lookup table.
        Sorted with MANE_Select first.
        """
        hits = MappingRepository.lookup_mapping_hits_for_position(
            chromosome=chromosome,
            position=position,
            genome_build=genome_build
        )
        if not hits:
            return []

        gene_ids = {h["gene_id"] for h in hits if h["gene_id"] is not None}
        protein_ids = {h["protein_id"] for h in hits if h["protein_id"] is not None}
        mapping_ids = {h["mapping_id"] for h in hits if h["mapping_id"] is not None}

        genes_by_id = GeneRepository.retrieve_genes_by_ids(gene_ids)
        proteins_by_id = ProteinRepository.retrieve_proteins_by_ids(protein_ids)

        mapping_pfam_domains = defaultdict(list)
        if mapping_ids:
            _session = db._make_scoped_session(options={})
            try:
                pfam_rows = _session.query(
                    MetaDomainMapping.mapping_id,
                    Interpro.ext_db_id,
                    Interpro.region_name
                ).join(
                    Interpro, Interpro.id == MetaDomainMapping.interpro_id
                ).filter(
                    MetaDomainMapping.mapping_id.in_(mapping_ids)
                ).all()

                seen_per_mapping = defaultdict(set)
                for row in pfam_rows:
                    if not row.ext_db_id:
                        continue

                    if row.ext_db_id in seen_per_mapping[row.mapping_id]:
                        continue

                    seen_per_mapping[row.mapping_id].add(row.ext_db_id)
                    mapping_pfam_domains[row.mapping_id].append({
                        "name": row.region_name or row.ext_db_id,
                        "ext_db_id": row.ext_db_id,
                    })
            except (AlchemyResourceClosedError, AlchemyOperationalError) as e:
                _cleanup_failed_session(_session)
                raise RecoverableError(str(e))
            except:
                _log.error(traceback.format_exc())
                raise
            finally:
                _session.remove()

        grouped = defaultdict(lambda: {
            "positions": [],
            "genome_build": "",
            "chromosome": "",
            "gene_name": "",
            "strand": "",
            "strand_symbol": "",
            "gencode_transcript": "",
            "refseq_transcript": "",
            "mane_transcript_type": "",
            "uniprot_ac": "",
            "codon": "",
            "codon_position": None,
            "amino_acid": "",
            "protein_position": None,
            "pfam_domains": []
        })

        for h in hits:
            gene = genes_by_id.get(h["gene_id"])
            if gene is None:
                continue

            protein = proteins_by_id.get(h["protein_id"])
            protein_position_1_based = h["uniprot_position"] + 1
            strand_value = str(h.get("strand") or "").strip()

            key = (
                gene.id,
                h["protein_id"],
                gene.gencode_transcription_id,
                protein_position_1_based,
                h["codon"],
                h["codon_position"],
                strand_value,
                h["uniprot_residue"]
            )

            item = grouped[key]
            item["positions"].append(h["chromosome_position"])
            item["genome_build"] = gene.genome_build
            item["chromosome"] = h["chromosome"]
            item["gene_name"] = gene.gene_name
            item["strand"] = strand_value if strand_value in {"+", "-"} else ""
            item["base_pair"] = h.get("base_pair", "")
            item["gencode_transcript"] = gene.gencode_transcription_id
            item["refseq_transcript"] = gene.refseq_transcript_id or ""
            item["mane_transcript_type"] = gene.mane_transcript_type or ""
            item["uniprot_ac"] = protein.uniprot_ac if protein is not None else ""
            item["codon"] = h["codon"] or ""
            item["codon_position"] = h.get("codon_position")
            item["amino_acid"] = protein_letters_1to3.get(h["uniprot_residue"], "") if h["uniprot_residue"] else ""
            item["protein_position"] = protein_position_1_based

            for pfam_domain in mapping_pfam_domains.get(h["mapping_id"], []):
                if pfam_domain["ext_db_id"] not in {d["ext_db_id"] for d in item["pfam_domains"]}:
                    item["pfam_domains"].append(pfam_domain)

        results = []
        for _, item in grouped.items():
            results.append({
                "gene_name": item["gene_name"],
                "genome_build": item["genome_build"],
                "chromosome": item["chromosome"],
                "genomic_position": position,
                "base_pair": item["base_pair"],
                "strand": item["strand"],
                "gencode_transcript": item["gencode_transcript"],
                "refseq_transcript": item["refseq_transcript"],
                "mane_transcript_type": item["mane_transcript_type"],
                "uniprot_ac": item["uniprot_ac"],
                "codon": item["codon"],
                "codon_position": item["codon_position"],
                "amino_acid": item["amino_acid"],
                "protein_position": item["protein_position"],
                "pfam_domains": sorted(item["pfam_domains"], key=lambda d: d["ext_db_id"]),
            })

        results.sort(key=lambda r: (
            0 if r["mane_transcript_type"] == "MANE_Select" else 1,
            r["gene_name"],
            r["gencode_transcript"],
            r["protein_position"] if r["protein_position"] is not None else 10 ** 9
        ))

        return results

    @staticmethod
    def get_mappings_for_multiple_protein_ids(_protein_ids):
        """Retrieves all mappings for a multiple Protein objects as {protein_id: [ Mapping ]}"""
        # Open as session
        _session = db._make_scoped_session(options={})

        try:
            _mappings_per_protein = {}
            for mapping in _session.query(Mapping).filter(Mapping.protein_id.in_(_protein_ids)).all():
                if not mapping.protein_id in _mappings_per_protein:
                    _mappings_per_protein[mapping.protein_id] = []

                _mappings_per_protein[mapping.protein_id].append(mapping)
            return _mappings_per_protein
        except (AlchemyResourceClosedError, AlchemyOperationalError) as e:
            _cleanup_failed_session(_session)
            raise RecoverableError(str(e))
        except:
            _log.error(traceback.format_exc())
            raise
        finally:
            # Close this session, thus all items are cleared and memory usage is kept at a minimum
            _session.remove()

    @staticmethod
    def get_mappings_for_protein(_protein):
        """Retrieves all mappings for a Protein object"""
        # Open as session
        _session = db._make_scoped_session(options={})

        try:
            return [x for x in _session.query(Mapping).filter(Mapping.protein_id == _protein.id).all()]
        except (AlchemyResourceClosedError, AlchemyOperationalError) as e:
            _cleanup_failed_session(_session)
            raise RecoverableError(str(e))
        except:
            _log.error(traceback.format_exc())
            raise
        finally:
            # Close this session, thus all items are cleared and memory usage is kept at a minimum
            _session.remove()


    @staticmethod
    def get_mappings_for_gene(_gene):
        """Retrieves all mappings for a Gene object, restricted to the gene's selected protein if present"""
        # Open as session
        _session = db._make_scoped_session(options={})

        try:
            query = _session.query(Mapping).filter(Mapping.gene_id == _gene.id)

            if _gene.protein_id is not None:
                query = query.filter(Mapping.protein_id == _gene.protein_id)

            return [x for x in query.all()]
        except (AlchemyResourceClosedError, AlchemyOperationalError) as e:
            _cleanup_failed_session(_session)
            raise RecoverableError(str(e))
        except:
            _log.error(traceback.format_exc())
            raise
        finally:
            # Close this session, thus all items are cleared and memory usage is kept at a minimum
            _session.remove()

class SequenceRepository:

    @staticmethod
    def get_aa_sequence(mappings, skip_asterix_at_end=False):
        """For a list of mappings, returns the amino acid sequence based on the uniprot positions
        If an asterix is expected at the end (e.g. a stop codon) there is the posibility to skip that"""
        _aa_sequence = ""
        mappings = {x.uniprot_position:x.uniprot_residue for x in mappings}
        for key in sorted(mappings, key=lambda x: (x is None, x)):
            # type check this is all the same protein and gene

            # type check if there are any gaps

            if skip_asterix_at_end and key is None:
                continue
            _aa_sequence+= mappings[key]
        return _aa_sequence

    @staticmethod
    def get_aa_region(sequence, region_start, region_stop):
        """For a given sequence, returns the sub-sequence
        based on the region_start and region stop"""
        # type check region start and region stop
        if region_start < 0 or region_stop > len(sequence) or region_stop < region_start:
            raise MalformedAARegionException("For sequence of length '"+str(len(sequence))+
                                              "': received a faulty  attempted to build a gene region where_region_stop < _region_start")


        # Return the sub sequence
        return sequence[region_start-1:region_stop]

    @staticmethod
    def get_cDNA_sequence(mappings):
        """For a list of mappings, returns the cDNA sequence based on the cDNA positions"""
        _cDNA_sequence = ""
        mappings = {x.cDNA_position:x.base_pair for x in mappings}
        for key in sorted(mappings.keys()):
            _cDNA_sequence+= mappings[key]
        return _cDNA_sequence
