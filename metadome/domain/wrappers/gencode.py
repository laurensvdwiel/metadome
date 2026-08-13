# TODO(dead-data-generation): unreachable at runtime.
#  Reached only through database_creation and install.py, which no longer exists.
import logging
from metadome.default_settings import GENCODE_HG_TRANSLATION_FILE,\
    GENCODE_SWISSPROT_FILE, GENCODE_HG_TRANSCRIPTION_FILE,\
    GENCODE_HG_ANNOTATION_FILE_GFF3, GENCODE_BASIC_FILE, GENCODE_REFSEQ_FILE
from metadome.domain.parsers import gff3
from Bio.Seq import translate
import urllib

_log = logging.getLogger(__name__)

class NoGeneTranslationsFoundException(Exception):
    pass

class NoSwissProtEntryFoundException(Exception):
    pass

class MissMatchTranscriptIDToMatchingTranscript(Exception):
    pass

class TranscriptionNotEncodingForTranslation(Exception):
    pass

class TranscriptionNotContainingCDS(Exception):
    pass

class TranscriptionStrandMismatchException(Exception):
    pass

def retrieveAllCharacterPositionsFromString(string_to_check, char_to_check):
    return [pos for pos, char in enumerate(string_to_check) if char == char_to_check]

def retrieveLongestTranslation(translations):
    """Retrieves the longest translation based on the sequence-length, within a list of matching transcripts"""
    longest_translation = None
    if len(translations) > 0:
        _log.debug("Finding longest translation")
        for translation in translations:
            _log.debug("Translation '"+str(translation['translation-name'])+"' has length "+str(translation['sequence-length']))
            if longest_translation is None or longest_translation['sequence-length']<translation['sequence-length']:
                longest_translation = translation
        
        _log.info("Using translation '"+str(longest_translation['translation-name'])+"', since it has the largest length")
    else:
        raise NoGeneTranslationsFoundException("Found zero matching transcripts")
        
    return longest_translation

def retrieveSwissProtIDs(gencode_ids):
    """Retrieve the SwissProt accession code and id based on multiple GencodeIDs"""
    results = {}
    with open(GENCODE_SWISSPROT_FILE) as gencode_swissprot:
        # read the lines in the file
        lines = gencode_swissprot.readlines()
        for line in lines:
            # check for each of the gencode_ids
            for gencode_id in gencode_ids:
                # check if the unique identifier is on the current line
                if gencode_id in line:
                    if not(gencode_id in results.keys()):
                        #Add the result to hits
                        tokens = line.split('\t')
                        
                        # Only add the translation to the translation list if the gene_name exactly matches the one we are looking for
                        if gencode_id == tokens[0]:                            
                            # add the results
                            result = {}
                            result['ac'] = tokens[1].strip()
                            result['swissprot_id'] = tokens[2].strip()
                            
                            # add the result to the results using the gencode id as an identifier
                            results[gencode_id] = result
                        else:
                            _log.warning("When searching for swissprot id using Gencode ID '"+gencode_id+"', found entry '"+tokens[0]+"', which only partially matched the id")
                    else:
                        _log.warning("For Gencode ID '"+gencode_id+"', found an additional swissprot entry: '"+line+"', we already found '"+str(results[gencode_id])+"'")
    if len(results) == 0:
        raise NoSwissProtEntryFoundException("For Gencode ID '"+gencode_id+"', found no matching swissprot entry")
    
    return results

def retrieveSwissProtID(gencode_id):
    """Retrieve the SwissProt accession code and id based on the GencodeID"""
    
    result = {}
    with open(GENCODE_SWISSPROT_FILE) as gencode_swissprot:
        # read the lines in the file
        lines = gencode_swissprot.readlines()
        for line in lines:
            # check if the unique identifier is on the current line
            if gencode_id in line:
                if len(result) == 0:
                    #Add the result to hits
                    tokens = line.split('\t')
                    
                    # Only add the translation to the translation list if the gene_name exactly matches the one we are looking for
                    if gencode_id == tokens[0]:
                        # add the results
                        result['ac'] = tokens[1].strip()
                        result['swissprot_id'] = tokens[2].strip()
                    else:
                        _log.warning("When searching for swissprot id using Gencode ID '"+gencode_id+"', found entry '"+tokens[0]+"', which only partially matched the id")
                else:
                    _log.warning("For Gencode ID '"+gencode_id+"', found an additional swissprot entry: '"+line+"', we already found '"+str(result)+"'")
    if len(result) == 0:
        raise NoSwissProtEntryFoundException("For Gencode ID '"+gencode_id+"', found no matching swissprot entry")
    
    return result

def retrieveStrandDirection_gencode(CDS_annotations):
    """Retrieves the strand direction from cds annotations
    To ensure that all strands are the same across the provided CDS annotations
    an error will be raised when this is not the case."""
    strand = None
    for cd in CDS_annotations:
        if strand is None: strand = cd.strand
        if cd.strand != strand: raise TranscriptionStrandMismatchException('Strands did not match for ')
    
    return strand    

def retrieveNucleotideSequence_gencode(translation):
    """Retrieve the gencode translations for a given gene_name based on the gene_name name.
    Gencode's gene transcription format is:
    Coding transcripts include transcripts with the following biotypes:
    protein_coding, nonsense_mediated_decay, non_stop_decay, polymorphic_pseudogene, IG_*_gene
     and TR_*_gene from Gencode v10 onwards (in previous versions they include only
     protein-coding and nonsense_mediated_decay biotypes).
     The pc_transcripts.fa file excludes all cases of transcripts with an internal stop codon."""
    
    matching_transcript = {}
    
    transcription_id = translation['transcription-id']
    
    _log.debug("Retrieving nucleotide sequence for transcription id")
    with open(GENCODE_HG_TRANSCRIPTION_FILE) as hg_trancription:
        # read the lines in the file
        lines = hg_trancription.readlines()
        for i in range(len(lines)):
            # retrieve current line
            line = lines[i]
            
            # check if the line is a fasta header, otherwise skip it 
            # Check if the string of the gene_name we are trying to find occurs in the line
            if line.startswith('>') and transcription_id in line:
                ## Parse the Line
                # Remove the fasta header syntax
                line = line[1:]
                if line.endswith('\n'):
                    # remove newline statement
                    line = line[:len(line)-1]
                
                # convert the fasta header to tokens
                tokens = line.split('|')
                
                # retrieve the various fields from the tokens
                gene_name_id = tokens[1] # == translation['gene_name-id']
                havana_gene_name_id = tokens[2] # == translation['Havana-gene_name-id']
                havana_translation_id = tokens[3] # == translation['Havana-translation-id']
                translation_name = tokens[4] # == translation['translation-name']
                gene_name = tokens[5] # == translation['gene-name']
                nucleotide_sequence_length = tokens[6]
                
                # perform checks that we are considering the correct sequence
                if transcription_id == tokens[0] \
                    and gene_name_id == translation['gene_name-id'] \
                    and havana_gene_name_id == translation['Havana-gene_name-id'] \
                    and havana_translation_id == translation['Havana-translation-id'] \
                    and translation_name == translation['translation-name'] \
                    and gene_name == translation['gene-name']:
                    
                    utr5=None
                    utr3=None
                    CDS=None
                    
                    for token in tokens[7:]:
                        if token.startswith('UTR3'):utr3=token
                        if token.startswith('CDS'):CDS=token
                        if token.startswith('UTR5'): utr5 = token
                            
                    if(utr3 is None):
                        _log.error("Non-fatal error: no UTR3 found for transcript '"+transcription_id+"'")
                    if(utr5 is None):
                        _log.error("Non-fatal error: no UTR5 found for transcript '"+transcription_id+"'")
                    if(CDS is None):
                        raise TranscriptionNotContainingCDS("No CDS found in transcript '"+transcription_id+"'")
                    
                    sequence = lines[i+1]
                    if sequence.endswith('\n'):
                        sequence = sequence[:len(sequence)-1]
                
                    # append transcript elements to the object
                    matching_transcript['sequence'] = sequence
                    matching_transcript['sequence_length'] = nucleotide_sequence_length
                    matching_transcript['UTR5'] = utr5
                    matching_transcript['UTR3'] = utr3
                    matching_transcript['CDS'] = CDS
                    
                    # check if the CDS is properly formatted
                    if CDS.startswith('CDS:'):
                        # retrieve the ranges of the Coding DNA Sequence
                        cds_range = [int(a) for a in CDS[4:].split('-')]
                        
                        # correct the first number for index mismatch
                        cds_range[0] = cds_range[0]-1
                        
                        # retrieve the sequence
                        matching_transcript['coding-sequence'] = sequence[cds_range[0]:cds_range[1]]
                        
                        # check if it can be correctly translated to the already found translation
                        translation_check = translate(matching_transcript['coding-sequence'])
                        if not(translation_check == translation['sequence']+'*'):
                            #translation may contain a Selenocysteine (DNA:'TGA', RNA:'UGA') or a Pyrrolysine (DNA:'TAG', RNA:'UAG')
                            U_pos = retrieveAllCharacterPositionsFromString(translation['sequence'], 'U')
                            O_pos = retrieveAllCharacterPositionsFromString(translation['sequence'], 'O')
                            Stop_pos = retrieveAllCharacterPositionsFromString(translation_check, '*')
                            
                            if  len(Stop_pos) < 1:
                                raise TranscriptionNotEncodingForTranslation("No stop codons present in the transcript '"+transcription_id+"'")
                            
                            if len(Stop_pos) < (len(U_pos)+len(O_pos)):
                                raise TranscriptionNotEncodingForTranslation("When translating the transcript '"+transcription_id+"' the translation did not match the nucleotide translation from BioPython")
                            
                            # filter out the U and O positions from the measured stop positions
                            true_stop_codons = [pos for pos in Stop_pos if not(pos in U_pos) and not(pos in O_pos)]
                            
                            if len(true_stop_codons) > 1:
                                raise TranscriptionNotEncodingForTranslation("After filtering out Selenocysteine and Pyrrolysine from the sequence, there were still other stop codons present before the final stop codon for the transcript '"+transcription_id+"'")
                            if len(true_stop_codons) < 1:
                                raise TranscriptionNotEncodingForTranslation("After filtering out Selenocysteine and Pyrrolysine from the sequence, there were No stop codons left for the transcript '"+transcription_id+"'")
                            
                            if len(U_pos) > 0 and len(O_pos) > 0:
                                # sequence contains both a Selenocysteine and a Pyrrolysine according to translation
                                _log.info("Both a Selenocysteine and a Pyrrolysine present in transcript '"+transcription_id+"'")
                            elif len(U_pos) > 0:
                                # sequence contains a Selenocysteine
                                _log.info("A Selenocysteine is present in transcript '"+transcription_id+"'")
                            elif len(O_pos) > 0:
                                # sequence contains a Pyrrolysine
                                _log.info("A Pyrrolysine is present in transcript '"+transcription_id+"'")                                                                
                    else:
                        raise MissMatchTranscriptIDToMatchingTranscript("Transcription ID "+transcription_id+" passed checks, but its CDS did not start with 'CDS:'")
                    
                    return matching_transcript
                else:
                    raise MissMatchTranscriptIDToMatchingTranscript("A missmatch occurred during transcript evaluation on: "+\
                    "transcription_id "+str(transcription_id == tokens[0])+\
                    " gene_name_id "+str(gene_name_id == translation['gene_name-id'])+\
                    " havana_gene_name_id "+str(havana_gene_name_id == translation['Havana-gene_name-id'])+\
                    " havana_translation_id "+str(havana_translation_id == translation['Havana-translation-id'])+\
                    " translation_name "+str(translation_name == translation['translation-name'])+\
                    " gene_name" +str(gene_name == translation['gene-name']))
                
                
    return matching_transcript
    
def retrieveCodingGenomicLocations_gencode(transcription):    
    transcription_records = []
    # Extract records that have the specific transcription id from the Gencode database 
    for record in gff3.parseGFF3(GENCODE_HG_ANNOTATION_FILE_GFF3, transcription['transcription-id']):
        transcription_records.append(record)
        
    _log.debug('Found'+str(len(transcription_records))+' records in the GFF file')
    
    # retrieve the Coding Domain Sequences
    CDS = [record for record in transcription_records if record.type=='CDS']
    
    _log.debug('Found'+str(len(CDS))+' CDS parts in the GFF file')
    
    return CDS

def retrieveMRNAValidatedTranslations_gencode(gene_name):
    """Retrieve the gencode basic translations for a given gene_name based on the gene_name name.
    Gencode basic set is validated on the mRNA level.
    Gencode's Basic gene translation format is:
    #bin|name|chrom|strand|txStart|txEnd|cdsStart|cdsEnd|exonCount|exonStarts|exonEnds|score|
     name2|cdsStartStat|cdsEndStat|exonFrames
    With name
     """
    # ensure we are dealing with the gene_name in upper characters
    if not(gene_name.isupper()):
        _log.warning("Provided gene_name '"+gene_name+"' is not provided in upper case")
    
    _log.info("Starting search for mRNA validated matching translations")
    
    with open(GENCODE_BASIC_FILE) as infile:
        for line in infile:
            if line.startswith("#"): continue
            # Check if the string of the gene_name we are trying to find occurs in the line
            if gene_name in line:
                parts = line.strip().split("\t")
                #If this fails, the file format is not standard-compatible
                assert len(parts) == 16
                #Normalize data
                gencode_basic_translation = {
                    "#bin" : None if parts[0] == "" else int(parts[0]),
                    "transcription_id" : None if parts[1] == "" else urllib.request.unquote(parts[1]),
                    "chrom" : None if parts[2] == "" else urllib.request.unquote(parts[2]),
                    "strand" : None if parts[3] == "" else urllib.request.unquote(parts[3]),
                    "txStart" : None if parts[4] == "" else int(parts[4]),
                    "txEnd" : None if parts[5] == "" else int(parts[5]),
                    "cdsStart" : None if parts[6] == "" else int(parts[6]),
                    "cdsEnd" : None if parts[7] == "" else int(parts[7]),
                    "exonCount" : None if parts[8] == "" else int(parts[8]),
                    "exonStarts" : None if parts[9] == "" else urllib.request.unquote(parts[9]),
                    "exonEnds" : None if parts[10] == "" else urllib.request.unquote(parts[10]),
                    "score" : None if parts[11] == "" else int(parts[11]),
                    "gene_name" : None if parts[12] == "" else urllib.request.unquote(parts[12]),
                    "cdsStartStat" : None if parts[13] == "" else urllib.request.unquote(parts[13]),
                    "cdsEndStat" : None if parts[14] == "" else urllib.request.unquote(parts[14]),
                    "exonFrames" : None if parts[15] == "" else urllib.request.unquote(parts[15]),
                }
                #yoeld the record
                yield gencode_basic_translation    

def retrieveGeneTranslations_gencode(gene_name):
    """Retrieve the gencode translations for a given gene_name based on the gene_name name.
    Gencode's gene translation format is:
    transcription-id|gene_name-id|Havana-gene_name-id (if the gene_name contains manually annotated
     transcripts, '-' otherwise)|Havana-translation-id (if this translation was
     manually annotated, '-' otherwise)|translation-name|gene_name-name|sequence-length"""
    # ensure we are dealing with the gene_name in upper characters
    if not(gene_name.isupper()):
        _log.warning("Provided gene_name '"+gene_name+"', this may lead to duplicate gene names")
    
    _log.info("Starting search for matching translations")
    # Look up the best matching translation
    matching_translations = []
    with open(GENCODE_HG_TRANSLATION_FILE) as hg_translation:
        # read the lines in the file
        lines = hg_translation.readlines()
        for i in range(len(lines)):
            # retrieve current line
            line = lines[i]
            
            # check if the line is a fasta header, otherwise skip it
            if line.startswith('>'):
                # Check if the string of the gene_name we are trying to find occurs in the line
                if gene_name in line:
                    ## Parse the Line
                    # Remove the fasta header syntax
                    line = line[1:]
                    if line.endswith('\n'):
                        # remove newline statement
                        line = line[:len(line)-1]
                    
                    tokens = line.split('|')
                    
                    # Only add the translation to the translation list if the gene_name exactly matches the one we are looking for
                    if gene_name == tokens[5]:
                        sequence = lines[i+1]
                        if sequence.endswith('\n'):
                            sequence = sequence[:len(sequence)-1]
                        
                        translation = {'transcription-id': tokens[0],
                                      'gene_name-id': tokens[1],
                                      'Havana-gene_name-id': tokens[2],
                                      'Havana-translation-id': tokens[3],
                                      'translation-name':tokens[4],
                                      'gene-name':tokens[5],
                                      'sequence-length':int(tokens[6]),
                                      'sequence':sequence}
                        
                        matching_translations.append(translation)
                        
    return matching_translations

def retrieve_all_protein_coding_gene_names():
    """Retrieves the protein-coding gene names for the Gencode dataset"""
    gene_names = set()
    with open(GENCODE_HG_TRANSLATION_FILE) as hg_translation:
        # read the lines in the file
        lines = hg_translation.readlines()
        for i in range(len(lines)):
            # retrieve current line
            line = lines[i]
            
            # check if the line is a fasta header, otherwise skip it
            if not line.startswith('>'):
                continue
            
            tokens = line.strip().split('|')
            
            if len(tokens) > 0:
                # add the gene name to the set
                gene_names.add(tokens[5])
                
    return list(gene_names)

def retrieve_refseq_identifiers_for_transcript(gencode_id):
    """Retrieves the refseq identifiers for a Gencode transcript"""
    result = {}
    result['NP'] = []
    result['NM'] = []
    result['NR'] = []
    with open(GENCODE_REFSEQ_FILE) as gencode_refseq:
        # read the lines in the file
        lines = gencode_refseq.readlines()
        for line in lines:
            # check if the unique identifier is on the current line
            if gencode_id in line:
                #Add the result to hits
                tokens = line.split('\t')
                
                # Only add the translation to the translation list if the gene_name exactly matches the one we are looking for
                if gencode_id == tokens[0]:
                    # add the results
                    for token in tokens[1:]:
                        token = token.strip()
                        if token.startswith('NP'):
                            result['NP'].append(token)
                        elif token.startswith('NM'):
                            result['NM'].append(token)
                        elif token.startswith('NR'):
                            result['NR'].append(token)
                        elif len(token) == 0:
                            continue
                        else:
                            _log.warning("When retrieving matching RefSeq ids for "+gencode_id+" unexpected token: "+token)    
    return result