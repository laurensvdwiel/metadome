import logging
import traceback

from sqlalchemy.orm.exc import NoResultFound, MultipleResultsFound
from sqlalchemy.sql.functions import func
from sqlalchemy.sql.expression import and_, distinct
from sqlalchemy.exc import ResourceClosedError as AlchemyResourceClosedError
from sqlalchemy.exc import OperationalError as AlchemyOperationalError
from psycopg2 import OperationalError as PsycopOperationalError

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
        except (AlchemyResourceClosedError, AlchemyOperationalError, PsycopOperationalError) as e:
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
        except (AlchemyResourceClosedError, AlchemyOperationalError, PsycopOperationalError) as e:
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
        except (AlchemyResourceClosedError, AlchemyOperationalError, PsycopOperationalError) as e:
            raise RecoverableError(str(e))
        except:
            _log.error(traceback.format_exc())
            raise
        finally:
            # Close this session, thus all items are cleared and memory usage is kept at a minimum
            _session.remove()
    
    @staticmethod
    def retrieve_transcript_id_for_multiple_gene_ids(_gene_ids):
        """Retrieves all gencode transcripts for multiple Gene objects as {gene_id: gencode_transcription_id}"""
        # Open as session
        _session = db._make_scoped_session(options={})

        try:
            _gencode_transcription_id_per_gene_id = {}
            for gene in _session.query(Gene).filter(Gene.id.in_(_gene_ids)).all():
                _gencode_transcription_id_per_gene_id[gene.id] = gene.gencode_transcription_id
            return _gencode_transcription_id_per_gene_id
        except (AlchemyResourceClosedError, AlchemyOperationalError, PsycopOperationalError) as e:
            raise RecoverableError(str(e))
        except:
            _log.error(traceback.format_exc())
            raise
        finally:
            # Close this session, thus all items are cleared and memory usage is kept at a minimum
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
        except (AlchemyResourceClosedError, AlchemyOperationalError, PsycopOperationalError) as e:
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
        except (AlchemyResourceClosedError, AlchemyOperationalError, PsycopOperationalError) as e:
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
                gene_name
                for gene_name in _session.query(Gene.gene_name).distinct(Gene.gene_name).all()
            ]

            try:
                RepositoryCacheHelper.cache_set(cache_key, result)
            except RepositoryCacheException:
                _log.exception("Repository cache failure in retrieve_all_gene_names_from_db (set key=%s)", cache_key)

            return result
        except (AlchemyResourceClosedError, AlchemyOperationalError, PsycopOperationalError) as e:
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
                genome_build
                for genome_build in _session.query(Gene.genome_build).distinct(Gene.genome_build).all()
            ]

            try:
                RepositoryCacheHelper.cache_set(cache_key, result)
            except RepositoryCacheException:
                _log.exception("Repository cache failure in retrieve_all_genome_builds_from_db (set key=%s)", cache_key)

            return result
        except (AlchemyResourceClosedError, AlchemyOperationalError, PsycopOperationalError) as e:
            raise RecoverableError(str(e))
        except:
            _log.error(traceback.format_exc())
            raise
        finally:
            # Close this session, thus all items are cleared and memory usage is kept at a minimum
            _session.remove()

    @staticmethod
    def retrieve_all_transcript_ids(genome_build, gene_name):
        """Retrieves all transcript ids for a gene name"""
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
            result = [
                transcript
                for transcript in _session.query(Gene).filter(
                    func.lower(Gene.gene_name) == gene_name.lower(),
                    func.lower(Gene.genome_build) == genome_build.lower()
                ).all()
            ]

            try:
                RepositoryCacheHelper.cache_set(cache_key, result)
            except RepositoryCacheException:
                _log.exception(
                    "Repository cache failure in retrieve_all_transcript_ids (set key=%s)",
                    cache_key
                )

            return result
        except (AlchemyResourceClosedError, AlchemyOperationalError, PsycopOperationalError) as e:
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
        except (AlchemyResourceClosedError, AlchemyOperationalError, PsycopOperationalError) as e:
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
        except (AlchemyResourceClosedError, AlchemyOperationalError, PsycopOperationalError) as e:
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
        except (AlchemyResourceClosedError, AlchemyOperationalError, PsycopOperationalError) as e:
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
        except (AlchemyResourceClosedError, AlchemyOperationalError, PsycopOperationalError) as e:
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
        except (AlchemyResourceClosedError, AlchemyOperationalError, PsycopOperationalError) as e:
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
        except (AlchemyResourceClosedError, AlchemyOperationalError, PsycopOperationalError) as e:
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

class MappingRepository:

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
        except (AlchemyResourceClosedError, AlchemyOperationalError, PsycopOperationalError) as e:
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
        except (AlchemyResourceClosedError, AlchemyOperationalError, PsycopOperationalError) as e:
            raise RecoverableError(str(e))
        except:
            _log.error(traceback.format_exc())
            raise
        finally:
            # Close this session, thus all items are cleared and memory usage is kept at a minimum
            _session.remove()


    @staticmethod
    def get_mappings_for_gene(_gene):
        """Retrieves all mappings for a Gene object"""
        # Open as session
        _session = db._make_scoped_session(options={})

        try:
            return [x for x in _session.query(Mapping).filter(Mapping.gene_id == _gene.id).all()]
        except (AlchemyResourceClosedError, AlchemyOperationalError, PsycopOperationalError) as e:
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
