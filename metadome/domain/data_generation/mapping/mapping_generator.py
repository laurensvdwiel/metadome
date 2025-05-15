from metadome.default_settings import UNIPROT_SPROT_SPECIES_FILTER
from metadome.domain.wrappers.gencode import retrieveGeneTranslations_gencode,\
    retrieveNucleotideSequence_gencode, retrieveCodingGenomicLocations_gencode,\
    retrieveStrandDirection_gencode, retrieveMRNAValidatedTranslations_gencode,\
    NoGeneTranslationsFoundException, TranscriptionNotContainingCDS,\
    TranscriptionNotEncodingForTranslation, TranscriptionStrandMismatchException,\
    MissMatchTranscriptIDToMatchingTranscript
from metadome.domain.wrappers.uniprot import retrieveIdenticalUniprotMatch,\
    NoUniProtACFoundException
from metadome.domain.data_generation.mapping.Gene2ProteinMapping import createMappingOfGeneTranscriptionToTranslationToProtein
from metadome.domain.models.gene import Gene
from metadome.domain.models.protein import Protein
from metadome.domain.models.interpro import Interpro

import logging

_log = logging.getLogger(__name__)
        
def generate_gene_to_swissprot_mapping(gene_name):
    """
    Given a gene_name, this method generates a mapping between swissprot and 
    every GENCODE Basic protein-coding translation for that gene
    """
    _log.info("Starting swissprot mapping for gene '"+gene_name+"'")
    
    # The result of this method; a dictionary containing all items that will be added to the database
    to_be_added_db_entries = {"genes":dict(), "proteins":dict(), "chromosome_positions":dict(), "mappings":dict()}    
    
    # retrieve all translations for the gene
    matching_translations = retrieveGeneTranslations_gencode(gene_name)
    
    # filter out translations that are not validated on the mRNA level
    mrna_translations =  [a['transcription_id'] for a in retrieveMRNAValidatedTranslations_gencode(gene_name)]
    mrna_filtered_translations = [translation for translation in matching_translations if translation['transcription-id'] in mrna_translations]
    n_translations = len(matching_translations)
    n_translations_after_mrna_filter  = len(mrna_filtered_translations)
    if n_translations > n_translations_after_mrna_filter:
        _log.debug("Filtered out '"+str(n_translations - n_translations_after_mrna_filter)+"' from the original '"+str(n_translations)+"' matching translation(s) for gene "+gene_name+", by focussing on mRNA level validation")
    
    # filter out any non-coding genes, by excluding any sequences that do not start with a methionine ('M')
    matching_coding_translations = [translation for translation in mrna_filtered_translations if translation['sequence'].startswith('M')]
    n_translations_after_coding_filter  = len(matching_coding_translations)
    if n_translations_after_mrna_filter > n_translations_after_coding_filter:
        _log.debug("Filtered out '"+str(n_translations - n_translations_after_coding_filter)+"' from the original '"+str(n_translations)+"' matching translation(s) for gene "+gene_name+", by checking if the sequence starts with a Methionine")
    
    # retrieve the lengths for the sequences
    sequences_lengths = [len(s['sequence']) for s in matching_coding_translations]
    _log.info("Found '"+str(len(matching_coding_translations))+"' matching protein coding translation(s) for gene "+gene_name+", with lengths: "+str(sequences_lengths))

    # start creation of the mapping between the gene and swissprot
    for matching_coding_translation in matching_coding_translations:
        # retrieve the nucleotide and coding sequence information
        try:
            gene_transcription = retrieveNucleotideSequence_gencode(matching_coding_translation)
            gene_transcription['CDS_annotation'] = retrieveCodingGenomicLocations_gencode(matching_coding_translation)
            matching_coding_translation['strand'] = retrieveStrandDirection_gencode(gene_transcription['CDS_annotation'])            
        except (NoGeneTranslationsFoundException, MissMatchTranscriptIDToMatchingTranscript,
                TranscriptionNotContainingCDS, TranscriptionStrandMismatchException,
                TranscriptionNotEncodingForTranslation) as e:
            _log.error(e)
            continue
        
        # Create the gene transcription entry
        to_be_added_db_entries["genes"][matching_coding_translation['transcription-id']] =  Gene(_strand=matching_coding_translation['strand'],
            _gene_name = matching_coding_translation['gene-name'],
            _gencode_transcription_id = matching_coding_translation['transcription-id'],
            _gencode_translation_name = matching_coding_translation['translation-name'],
            _gencode_gene_id = matching_coding_translation['gene_name-id'],
            _havana_gene_id = matching_coding_translation['Havana-gene_name-id'],
            _havana_translation_id = matching_coding_translation['Havana-translation-id'],
            _sequence_length = matching_coding_translation['sequence-length'])
        
        # retrieve uniprot match
        try:
            uniprot = retrieveIdenticalUniprotMatch(matching_coding_translation, species_filter=UNIPROT_SPROT_SPECIES_FILTER)
        except (NoUniProtACFoundException) as e:
            _log.error("For gene '"+str(gene_name)+
                                    "' with translation '"+str(matching_coding_translation['translation-name'])+
                                    "', no match could be made with uniprot due to: "+str(e))
            continue
        
        # Create the uniprot_result database entry
        to_be_added_db_entries["proteins"][matching_coding_translation['transcription-id']] = Protein(_uniprot_ac = uniprot['uniprot_ac'],
                _uniprot_name = uniprot['uniprot_name'],
                _source = uniprot['database'])
        
        # create the mapping between both sequences
        mappings = createMappingOfGeneTranscriptionToTranslationToProtein(gene_transcription, matching_coding_translation, uniprot)
        to_be_added_db_entries["mappings"][matching_coding_translation['transcription-id']] = mappings
        
        _log.info("For gene '"+str(gene_name)+
                                    "' the translation '"+str(matching_coding_translation['translation-name'])+
                                    "' is matched and mapped with '"+str(uniprot['uniprot_ac'])+"'")
    
    return to_be_added_db_entries