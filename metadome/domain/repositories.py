import logging
import traceback

from sqlalchemy.orm.exc import NoResultFound, MultipleResultsFound
from sqlalchemy.sql.functions import func
from metadome.default_settings import GENE_NAMES_FILE
from sqlalchemy.sql.expression import and_, distinct
from sqlalchemy.exc import ResourceClosedError as AlchemyResourceClosedError
from sqlalchemy.exc import OperationalError as AlchemyOperationalError
from psycopg2 import OperationalError as PsycopOperationalError

from metadome.database import db
from metadome.domain.models.mapping import Mapping
from metadome.domain.models.gene import Gene
from metadome.domain.models.protein import Protein
from metadome.domain.models.interpro import Interpro
from metadome.domain.error import RecoverableError

_log = logging.getLogger(__name__)

class RepositoryException(Exception):
    pass

class MalformedAARegionException(Exception):
    pass


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
    def retrieve_all_gene_names_from_file():
        """Retrieves all gene names present in the static file"""
        try:
            with open(GENE_NAMES_FILE, 'r') as gene_names_file:
                return gene_names_file.read().splitlines()
        except OSError:
            return []

    @staticmethod
    def retrieve_all_gene_names_from_db():
        """Retrieves all gene names present in the database"""
        # Open as session
        _session = db._make_scoped_session(options={})

        try:
            return [gene_name for gene_name in _session.query(Gene.gene_name).distinct(Gene.gene_name).all()]
        except (AlchemyResourceClosedError, AlchemyOperationalError, PsycopOperationalError) as e:
            raise RecoverableError(str(e))
        except:
            _log.error(traceback.format_exc())
            raise
        finally:
            # Close this session, thus all items are cleared and memory usage is kept at a minimum
            _session.remove()

    @staticmethod
    def retrieve_all_transcript_ids(gene_name):
        """Retrieves all transcript ids for a gene name"""
        # Open as session
        _session = db._make_scoped_session(options={})

        try:
            return [transcript for transcript in _session.query(Gene).filter(func.lower(Gene.gene_name) == gene_name.lower()).all()]
        except (AlchemyResourceClosedError, AlchemyOperationalError, PsycopOperationalError) as e:
            raise RecoverableError(str(e))
        except:
            _log.error(traceback.format_exc())
            raise
        finally:
            # Close this session, thus all items are cleared and memory usage is kept at a minimum
            _session.remove()

    @staticmethod
    def retrieve_gene(transcription_id):
        """Retrieves the gene object for a given transcript id"""
        # Open as session
        _session = db._make_scoped_session(options={})
        try:
            return _session.query(Gene).filter(Gene.gencode_transcription_id == transcription_id).one()
        except (AlchemyResourceClosedError, AlchemyOperationalError, PsycopOperationalError) as e:
            raise RecoverableError(str(e))
        except MultipleResultsFound as e:
            error_message = "GeneRepository.retrieve_gene(transcription_id): Multiple results found while expecting uniqueness for transcription_id '"+str(transcription_id)+"'. "+str(e)
            _log.error(error_message)
            raise RepositoryException(error_message)
        except NoResultFound as e:
            error_message = "GeneRepository.retrieve_gene(transcription_id): Expected results but found none for transcription_id '"+str(transcription_id)+"'. "+str(e)
            _log.error(error_message)
            raise RepositoryException(error_message)
        except Exception as e:
            error_message = "GeneRepository.retrieve_gene(transcription_id): Unexpected exception for transcription_id '"+str(transcription_id)+"'. "+str(e)
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
