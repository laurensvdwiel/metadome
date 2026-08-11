# TODO(dead-data-generation): unreachable at runtime.
#  Reached only via mapping_generator -> database_creation -> install.py, which no longer exists.
import logging
from metadome.domain.wrappers.blast import run_blast, interpret_blast_as_uniprot
from metadome.domain.wrappers.gencode import retrieveSwissProtIDs,\
    NoSwissProtEntryFoundException
from metadome.default_settings import UNIPROT_SPROT_CANONICAL,\
    UNIPROT_MAX_BLAST_RESULTS, UNIPROT_SPROT_CANONICAL_AND_ISOFORM,\
    UNIPROT_SPROT_ISOFORM

_log = logging.getLogger(__name__)


class NoUniProtACFoundException(Exception):
    pass

class NoneExistingUniprotDatabaseProvidedException(Exception):
    pass

class MissMatchingSwissProtEntriesFoundException(Exception):
    pass

def retrieveIdenticalUniprotMatch(geneTranslation, uniprot_database="sp", species_filter=None):
    """Retrieves the top swiss-/Uniprot match for the provided sequence"""
    # Which database to blast to
    if uniprot_database=="sp":
        _log.debug('Querying the gene translation to Uniprot_Sprot database')
        blast_db=UNIPROT_SPROT_CANONICAL_AND_ISOFORM
    else:
        raise NoneExistingUniprotDatabaseProvidedException("Possible options are 'sp' for SwissProt and 'tr' for Uniprot_Trembl, provided was '"+str(uniprot_database)+"'")
    
    # retrieve the results
    sequence_results = retrieveMatchingUniprotSequences(geneTranslation, blast_db, species_filter)
    _log.debug("Retrieved '"+str(len(sequence_results))+"' results from blasting to Uniprot")
    
    # retrieve the best matching result
    top_result = None
    for result in sequence_results:
        _result_sequence = getUniprotSequence(result['accession_code'], blast_db)
        if _result_sequence == geneTranslation['sequence']:
            _log.debug('input and blast result sequences are identical')
            top_result = result
            break
    
    if top_result is None:
        raise NoUniProtACFoundException("No identical uniprot sequence found for "+geneTranslation['gene-name']+" after blast search")
    
    # construct the result    
    uniprot_result = {
        "sequence": getUniprotSequence(top_result['accession_code'], blast_db),
        "database": "swissprot" if top_result['database_id']=='sp' else "uniprot_trembl",
        "uniprot_name": top_result['entry_name'],
        "uniprot_ac": top_result['accession_code'],
        "blast_result":top_result,
    }
    
    return uniprot_result
    

def getUniprotSequence(accession_code , blast_db=UNIPROT_SPROT_CANONICAL_AND_ISOFORM):
    if blast_db != UNIPROT_SPROT_CANONICAL_AND_ISOFORM and blast_db != UNIPROT_SPROT_CANONICAL and blast_db != UNIPROT_SPROT_ISOFORM:
        raise NotImplementedError('Only swissprot is supported for now, uniprot needs to be handled by a database structure')
    
    sequence = ""
    with open(blast_db) as fasta_file:
        fasta_lines = fasta_file.read().splitlines()
        n_total_lines = len(fasta_lines)
        for i in range(n_total_lines):
            if not(fasta_lines[i].startswith('>')):continue
            if accession_code in fasta_lines[i]:
                for j in range(i+1, n_total_lines):
                    if fasta_lines[j].startswith('>'):
                        return sequence
                    if j < n_total_lines:
                        sequence = sequence+fasta_lines[j]
                    else:
                        return sequence
                
    return sequence    

def retrieveMatchingUniprotSequences(geneTranslation, blast_db=UNIPROT_SPROT_CANONICAL_AND_ISOFORM, species_filter=None):
    """Blasts the translated gene sequence to the uniprot/swissprot fasta database and returns 
     the blast results as a list of dictionaries with keys: database_id, accession_code,
     entry_name, qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart,
     send, evalue, bitscore"""
    sequence_results = []
    blast_results = run_blast(sequence=geneTranslation['sequence'], blastdb=blast_db, max_target_seqs=UNIPROT_MAX_BLAST_RESULTS)
    
    if not len(blast_results) == 0 :
        sequence_results = [interpret_blast_as_uniprot(blast_result) for blast_result in blast_results]
    
    if not species_filter is None:
        if blast_db == UNIPROT_SPROT_CANONICAL or blast_db == UNIPROT_SPROT_ISOFORM or blast_db == UNIPROT_SPROT_CANONICAL_AND_ISOFORM:
            sequence_results = [entry for entry in sequence_results if entry['entry_name'].split("_")[1] == species_filter]
        else:
            raise NotImplementedError("There is no species filtering available at the moment for uniprot (non-swissprot) entries")
    
    if len(sequence_results) == 0 :
        raise NoUniProtACFoundException("No uniprot sequence found for "+geneTranslation['gene-name'])
    
    
    return sequence_results
