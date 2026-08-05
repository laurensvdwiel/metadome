from metadome.domain.data_generation.mapping.meta_domain_mapping import retrieve_pfam_aligned_codons
from metadome.domain.models.entities.gene_region import GenomeBuild
from metadome.domain.services.annotation.codon_annotation import annotate_ClinVar_SNVs_for_codons,\
    annotate_gnomAD_SNVs_for_codons
from metadome.domain.models.entities.single_nucleotide_variant import SingleNucleotideVariant, VariantSource
from metadome.domain.models.entities.codon import Codon
from metadome.default_settings import METADOMAIN_DIR,\
    METADOMAIN_MAPPING_FILE_NAME, METADOMAIN_DETAILS_FILE_NAME,\
    METADOMAIN_CACHE_MAXSIZE, METADOMAIN_SNV_ANNOTATION_FILE_NAME_PER_SOURCE

import pandas as pd
import numpy as np
import json
import os

import logging

_log = logging.getLogger(__name__)

from collections import OrderedDict

# Per-process LRU cache of built MetaDomain objects, keyed by (domain_id, genome_build.value).
# MetaDomain is read-only after construction, so reuse across transcripts is safe.
_METADOMAIN_CACHE = OrderedDict()

class MetaDomainException(Exception):
    pass

class UnsupportedMetaDomainIdentifier(MetaDomainException):
    pass

class ConsensusPositionOutOfBounds(MetaDomainException):
    pass

class NotInMetaDomain(MetaDomainException):
    pass

class MetaDomain(object):
    """
    MetaDomain Model Entity
    Used for representation of meta domains
    
    Variables
    name                                    description
    domain_id                               str the id / accession code of this domain
    genome_build                            GenomeBuild the genome build for which this metadomain is constructed
    consensus_length                        int length of the domain consensus
    consensus_positions                     set(int) of all consensus positions for this domain
    n_proteins                              int number of unique proteins containing this domain
    n_instances                             int number of unique instances containing this domain
    n_transcripts                           int number of unique transcripts containing this domain
    meta_domain_mapping                     pandas.DataFrame containing all codons annotated with corresponding consensus position
    meta_domain_annotation_per_source       pandas.DataFrame containing all SNVs with corresponding consensus position split by source
    """

    def get_annotation(self, variant_source):
        """The SNV annotation DataFrame for one variant source, read from disk on first use."""
        if variant_source in self._annotation_per_source:
            return self._annotation_per_source[variant_source]

        annotation_file = METADOMAIN_DIR + self.genome_build.value + '/' + self.domain_id + '/' + \
                          METADOMAIN_SNV_ANNOTATION_FILE_NAME_PER_SOURCE[variant_source.value]

        if not os.path.exists(annotation_file):
            raise MetaDomainException("No '{}' SNV annotation for domain '{}', expected it at '{}'".format(
                variant_source, self.domain_id, annotation_file))

        _log.info("Reading '{}'".format(annotation_file))
        try:
            self._annotation_per_source[variant_source] = pd.read_csv(annotation_file)
        except pd.errors.EmptyDataError:
            # A source with no SNVs for this domain is written as an empty file
            self._annotation_per_source[variant_source] = pd.DataFrame()

        return self._annotation_per_source[variant_source]

    def get_annotated_SNVs_for_consensus_position(self, consensus_position, variant_sources=None):
        """Retrieves SNVs for this consensus position as:
        {SingleNucleotideVariant.unique_var_str_representation():  dict()}"""
        snvs = dict()
        
        if consensus_position <= 0:
            raise ConsensusPositionOutOfBounds("The provided consensus position ('{}') is below zero, this position foes not exist".format(consensus_position))
        if consensus_position > self.consensus_length:
            raise ConsensusPositionOutOfBounds("The provided consensus position ('{}') is above the maximum consensus length ('{}'), this position foes not exist".format(consensus_position, self.consensus_length))
        if variant_sources is None:
            variant_sources = list(VariantSource)

        for variant_source in variant_sources:
            annotation = self.get_annotation(variant_source)
            if 'consensus_pos' not in annotation.columns:
                continue  # A source with no annotated SNVs is a legitimately empty DataFrame (no columns) — skip it.

            # Retrieve all codons aligned to the consensus position
            aligned_to_position = annotation[annotation.consensus_pos == consensus_position].to_dict('records')

            for snv in aligned_to_position:
                # aggregate duplicate chromosomal regions
                if not snv['unique_snv_str_representation'] in snvs.keys():
                    snvs[snv['unique_snv_str_representation']] = []

                # add the codon to the dictionary
                snvs[snv['unique_snv_str_representation']].append(snv)
                
        # return the codons that correspond to this position
        return snvs
    
    def get_consensus_positions_for_uniprot_position(self, uniprot_ac, uniprot_position):
        """Retrieves the consensus positions for this MetaDomain
        based on the uniprot ac and position"""
        consensus_positions = []
        # Retrieve all codons aligned to the consensus position
        aligned_to_position = self.meta_domain_mapping[(self.meta_domain_mapping.uniprot_ac == uniprot_ac) &\
                                                       (self.meta_domain_mapping.uniprot_position == uniprot_position)]
        
        # check if there are any matches
        if len(aligned_to_position) > 0:
            # check how many matches and type check if all positions are the same
            unique_consensus_positions = pd.unique(aligned_to_position.consensus_pos)
            
            if len(unique_consensus_positions) > 1:
                _log.warning("There are more than one consensus positions assigned ('"+str(unique_consensus_positions)+"') to the protein '"+str(uniprot_ac)+"' for position '"+str(uniprot_position)+"'")
            
            consensus_positions = []
            for c_pos in unique_consensus_positions:
                consensus_positions.append(int(c_pos))
        else:
            _log.info("No alignment for domain '"+str(self.domain_id)+"' for uniprot_ac '"+str(uniprot_ac)+"' on position '"+str(uniprot_position)+"'" )
        
        return consensus_positions
    
    def get_codon_for_transcript_and_position(self, transcript_id, protein_position):
        """Construct the codon for a provided position"""
        # Retrieve all codons aligned to the consensus position
        aligned_to_position = self.meta_domain_mapping[(self.meta_domain_mapping.gencode_transcription_id == transcript_id) & (self.meta_domain_mapping.uniprot_position == protein_position)].to_dict('records')
        
        if len(aligned_to_position) == 0:
            raise NotInMetaDomain("No codons found to be aligned for metadomain '"+str(self.domain_id)+"' for transcript '"+str(transcript_id)+"' at position '"+str(protein_position)+"'")
        else:
            return Codon.initializeFromDict(aligned_to_position[0])
        
    
    def get_codons_aligned_to_consensus_position(self, consensus_position):
        """Retrieves codons for this consensus position as:
        {Codon.unique_str_representation(): Codon}"""
        codons = dict()

        if consensus_position <= 0:
            raise ConsensusPositionOutOfBounds("The provided consensus position ('{}') is below zero, this position foes not exist".format(consensus_position))
        if consensus_position > self.consensus_length:
            raise ConsensusPositionOutOfBounds("The provided consensus position ('{}') is above the maximum consensus length ('{}'), this position foes not exist".format(consensus_position, self.consensus_length))

        # Retrieve all codons aligned to the consensus position
        aligned_to_position = self.meta_domain_mapping[self.meta_domain_mapping.consensus_pos == consensus_position].to_dict('records')
        
        # first check if the consensus position is present in the mappings_per_consensus_pos
        if len(aligned_to_position) >0:
            for codon_dict in aligned_to_position:
                # initialize a codon from the dataframe row
                codon = Codon.initializeFromDict(codon_dict)
                
                # aggregate duplicate chromosomal regions
                if not codon.unique_str_representation() in codons.keys():
                    codons[codon.unique_str_representation()] = []

                # add the codon to the dictionary
                codons[codon.unique_str_representation()].append(codon)
                
        # return the codons that correspond to this position
        return codons
    
    def get_alignment_depth_for_consensus_position(self, consensus_position):
        """Retrieves the number of aligned codons for this consensus position"""
        if consensus_position <= 0:
            raise ConsensusPositionOutOfBounds("The provided consensus position ('{}') is below zero, this position foes not exist".format(consensus_position))
        if consensus_position > self.consensus_length:
            raise ConsensusPositionOutOfBounds("The provided consensus position ('{}') is above the maximum consensus length ('{}'), this position foes not exist".format(consensus_position, self.consensus_length))
        
        # Retrieve all codons aligned to the consensus position
        aligned_to_position = self.meta_domain_mapping[self.meta_domain_mapping.consensus_pos == consensus_position].to_dict('records')

        unique_keys = [Codon.initializeFromDict(codon_dict).unique_str_representation() for codon_dict in aligned_to_position]
        return len(np.unique(unique_keys))
    
    def get_max_alignment_depth(self):
        alignment_depths = [ self.get_alignment_depth_for_consensus_position(consensus_position) for consensus_position in self.consensus_positions]
        return int(np.max(alignment_depths))
    
    def annotate_metadomain(self, reannotate=False):
        """Ensure the per-source SNV annotations for this meta domain exist on disk;
        reading them back is deferred until get_annotation()"""
        meta_domain_dir = METADOMAIN_DIR + self.genome_build.value + '/' + self.domain_id
        annotation_files = {VariantSource(source): meta_domain_dir + '/' + file_name
                              for source, file_name in METADOMAIN_SNV_ANNOTATION_FILE_NAME_PER_SOURCE.items()}

        if not reannotate:
            if all(os.path.exists(f) for f in annotation_files.values()):
                _log.info('Previously annotated MetaDomain available for domain id: '+str(self.domain_id))
                return

        _log.info('Start annotation of MetaDomain for domain id: '+str(self.domain_id))

        annotation_per_source = {source: [] for source in annotation_files.keys()}
        for consensus_position in self.consensus_positions:
            meta_codons = self.get_codons_aligned_to_consensus_position(consensus_position)

            # Annotate ClinVar and gnomAD SNVs
            for unique_str_repr in meta_codons.keys():
                for snv in annotate_ClinVar_SNVs_for_codons(meta_codons[unique_str_repr], self.genome_build):
                    snv['consensus_pos'] = consensus_position
                    annotation_per_source[VariantSource.clinvar].append(snv)
                for snv in annotate_gnomAD_SNVs_for_codons(meta_codons[unique_str_repr], self.genome_build):
                    snv['consensus_pos'] = consensus_position
                    annotation_per_source[VariantSource.gnomad].append(snv)

        # save each annotation to disk
        for source, annotation in annotation_per_source.items():
            annotation = pd.DataFrame(annotation)
            _tmp = annotation_files[source] + '.tmp.' + str(os.getpid())
            annotation.to_csv(_tmp)
            os.replace(_tmp, annotation_files[source])
            self._annotation_per_source[source] = annotation

        _log.info('Finished annotation of MetaDomain for domain id: '+str(self.domain_id))
    
    def __init__(self, domain_id, genome_build, consensus_length, consensus_positions, n_instances, meta_domain_mapping, meta_domain_annotation_per_source=None):
        self.domain_id = domain_id
        self.genome_build = genome_build
        self.consensus_length = consensus_length
        self.consensus_positions = consensus_positions
        self.n_instances = n_instances
        self.meta_domain_mapping = meta_domain_mapping
        self._annotation_per_source = dict(meta_domain_annotation_per_source or {})
        
        # derive from meta_domain_mapping
        self.n_proteins = len(pd.unique(self.meta_domain_mapping.uniprot_ac))
        self.n_transcripts = len(pd.unique(self.meta_domain_mapping.gencode_transcription_id))
        
    @classmethod
    def initializeFromDomainID(cls, domain_id, genome_build, recreate=False):
        _log.info('Start initialization of MetaDomain for domain id: '+str(domain_id))
        
        # Set values needed for construction of this class
        consensus_length = 0
        meta_domain_mapping = []
        consensus_positions = set()

        # check if genome build is GenomeBuild - an enum defined in gene_region
        if not (genome_build is GenomeBuild.GRCh37 or genome_build is GenomeBuild.GRCh38):
            raise MetaDomainException("Expected genome_build to be of type GenomeBuild, instead received '"+str(type(genome_build))+"'")

        cache_key = (domain_id, genome_build.value)
        if not recreate and cache_key in _METADOMAIN_CACHE:
            _METADOMAIN_CACHE.move_to_end(cache_key)
            return _METADOMAIN_CACHE[cache_key]

        # Double check this conserns a Pfam domain
        if domain_id.startswith('PF'):
            # check if a Meta Domain is already mapped
            meta_domain_dir = METADOMAIN_DIR+genome_build.value+'/'+domain_id
            meta_domain_details_file = meta_domain_dir+'/'+METADOMAIN_DETAILS_FILE_NAME
            meta_domain_mapping_file = meta_domain_dir+'/'+METADOMAIN_MAPPING_FILE_NAME
            
            # first check if the metadomain dir exist
            os.makedirs(meta_domain_dir, exist_ok=True)
            
            # Check if the mapping has previously been build already
            if os.path.exists(meta_domain_mapping_file) and os.path.exists(meta_domain_details_file) and not recreate:
                # The mapping exists, load it
                _log.info('Loading previously build creation of MetaDomain for domain id: '+str(domain_id)+' genome build: '+str(genome_build.value))
                # Read the files
                _log.info("Reading '{}'".format(meta_domain_mapping_file))
                meta_domain_mapping = pd.read_csv(meta_domain_mapping_file)
                _log.info("Reading '{}'".format(meta_domain_details_file))
                with open(meta_domain_details_file) as f:
                    meta_domain_details = json.load(f)
                    
                consensus_length = meta_domain_details['consensus_length']
                n_instances = meta_domain_details['n_instances']
                # Populate consensus_positions from the loaded mapping
                consensus_positions = set(meta_domain_mapping['consensus_pos'].unique())
            else:
                # The mapping does not exists yet, we need to create it
                _log.info('Start creation of MetaDomain for domain id: '+str(domain_id)+' genome build: '+str(genome_build.value))
               
                # create the meta domain mapping alignment
                meta_codons_per_consensus_pos, consensus_length, n_instances = retrieve_pfam_aligned_codons(domain_id, genome_build.value)
                
                # create the meta_domain_details
                meta_domain_details = {}
                meta_domain_details['consensus_length'] = consensus_length
                meta_domain_details['n_instances'] = n_instances
                
                # create the dataframe context for this meta_domain
                for consensus_pos in meta_codons_per_consensus_pos.keys():
                    consensus_positions.add(int(consensus_pos))
                    for codon in meta_codons_per_consensus_pos[consensus_pos]:
                        _meta_codon = codon.toDict()
                        _meta_codon['consensus_pos'] = consensus_pos
                        _meta_codon['domain_id'] = domain_id
                        
                        meta_domain_mapping.append(_meta_codon)

                # convert meta_domain_mapping to a pandas Dataframe
                meta_domain_mapping = pd.DataFrame(meta_domain_mapping)

                # A domain without aligned codons cannot become a metadomain; bail out before writing
                # an empty mapping file that would break every later read of this domain
                if meta_domain_mapping.empty:
                    raise UnsupportedMetaDomainIdentifier("No aligned codons found for domain '" + str(domain_id) + "' for genome build '" + str(genome_build.value) + "'")
                
                ## Save the results to disk
                # save meta_domain_details
                _tmp = meta_domain_details_file + '.tmp.' + str(os.getpid())
                with open(_tmp, 'w') as f:
                    json.dump(meta_domain_details, f)
                os.replace(_tmp, meta_domain_details_file)

                # save meta_domain_mapping to disk
                _tmp = meta_domain_mapping_file + '.tmp.' + str(os.getpid())
                meta_domain_mapping.to_csv(_tmp, index=False)
                os.replace(_tmp, meta_domain_mapping_file)
        else:
            raise UnsupportedMetaDomainIdentifier("Expected a Pfam domain, instead the identifier '"+str(domain_id)+"' was received")
        
        # Attempt to create the object
        meta_domain = cls(domain_id, genome_build, consensus_length, consensus_positions, n_instances, meta_domain_mapping)
        
        # Annotate this meta domain
        meta_domain.annotate_metadomain()

        # Cache for reuse across transcripts in this process (LRU-bounded).
        _METADOMAIN_CACHE[cache_key] = meta_domain
        _METADOMAIN_CACHE.move_to_end(cache_key)
        while len(_METADOMAIN_CACHE) > METADOMAIN_CACHE_MAXSIZE:
            _METADOMAIN_CACHE.popitem(last=False)
        
        # return the object
        return meta_domain
    
    def __repr__(self):
        return "<MetaDomain(domain_id='%s', consensus_length='%s', n_proteins='%s', n_instances='%s')>" % (
                            self.domain_id, self.consensus_length, self.n_proteins, self.n_instances)
        