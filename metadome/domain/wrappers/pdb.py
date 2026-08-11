# TODO(dead-data-generation): unreachable at runtime.
#  Reached only via mapping_generator -> database_creation -> install.py, which no longer exists.
import logging
from metadome.domain.wrappers.blast import run_blast, interpretation_to_string,\
    interpret_blast_as_pdb
from metadome.default_settings import PDB_SEQRES_FASTA, PDB_STRUCTURE_DIR,\
    MINIMAL_BLASTPE_VALUE, MINIMAL_TRANSLATION_TO_STRUCTURE_PIDENT_VALUE,\
    MINIMAL_XRAY_STRUCTURE_RESOLUTION
import warnings
import numpy as np
import gzip
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import is_aa
from Bio.SeqUtils import seq1
from Bio.PDB.PDBExceptions import PDBConstructionWarning

_log = logging.getLogger(__name__)

class NoPDBStructureFoundException(Exception):
    pass

class NoPDBFileFoundException(Exception):
    pass

class AtomicSequenceIDMappingFailedException(Exception):
    pass

def retrieveStructureFromPDB(pdb_id):
    """Retrieves and reads the structure from a PDB id, returns the structure as a 
     BioPython PDB Parsed structure"""
    try:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', PDBConstructionWarning)
            pdb_file = gzip.open(PDB_STRUCTURE_DIR+ "xx/".replace("xx",pdb_id[1:3])+"pdbxxxx.ent.gz".replace('xxxx', pdb_id), 'rt', encoding='utf-8')
            structure = PDBParser().get_structure(pdb_id, pdb_file)
            return structure
    except FileNotFoundError as e:
        raise NoPDBFileFoundException("No PDB file found for pdb_id '"+pdb_id+"', with FileNotFoundError: "+str(e))    

def retrieveAtomicStructure(pdb_sequence):
    """Retrieves the atomic structure for a single chain of a PDB file, based on the measured structure"""
    pdb_structure = retrieveStructureFromPDB(pdb_sequence['pdb_id'])
    
    return {residue.get_id()[1]:residue.get_resname() for residue in pdb_structure[0][pdb_sequence['chain_id']].get_residues() if is_aa(residue)}         
    
def retrieveAtomicStructureMapping(pdb_sequence, translation_to_structure_mapping):
    """Retrieves the mapping to atoms in a PDB file, based on the measured structure"""
    measured_structure = retrieveAtomicStructure(pdb_sequence)
    
    seq_ids = [x for x in sorted(measured_structure.keys())]
    
    # aligned_sequence to atomic structure mapping
    atomic_structure_mapping = {}
    
    index = 0
    for i in range(len(translation_to_structure_mapping['secondary_sequence'])):
        if translation_to_structure_mapping['secondary_sequence'][i] != '-':
            if seq1(measured_structure[seq_ids[index]]) != translation_to_structure_mapping['secondary_sequence'][i]:
                raise AtomicSequenceIDMappingFailedException("Alternative mapping for atomic sequence to atomic seqids failed for pdb structure "+pdb_sequence['pdb_id']+", chain "+pdb_sequence['chain_id'])
            atomic_structure_mapping[i] = seq_ids[index]
            index+=1
            
    return atomic_structure_mapping

def retrieveAtomicStructureSequence(pdb_sequence):
    """Retrieves the sequence of a PDB file, based on the measured structure"""
    measured_structure = retrieveAtomicStructure(pdb_sequence)
    
    seq_ids = [x for x in measured_structure.keys()]
    measure_structure_sequence = ""
    for seq_id in range(np.min(seq_ids), np.max(seq_ids)+1):
        if seq_id in measured_structure.keys():
            measure_structure_sequence += seq1(measured_structure[seq_id])
        else:
            measure_structure_sequence += "-"
    
    return measure_structure_sequence
    
def retrieveMatchingSeqresSequences(geneTranslation):
    """Blasts the translated gene sequence to the PDB fasta file and returns the blast results
    as a list of dictionaries with keys: pdb_id, chain_id, qseqid, sseqid, pident, length, 
     mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore"""
    seqres_sequence_results = []
    blast_results = run_blast(geneTranslation['sequence'], PDB_SEQRES_FASTA)
    
    if len(blast_results) == 0 :
        raise NoPDBStructureFoundException("No structures found for "+geneTranslation['gene-name'])
    else:
        seqres_sequence_results = [interpret_blast_as_pdb(blast_result) for blast_result in blast_results]
    
    return seqres_sequence_results

def cleanBlastResults(blastResults, length_of_translation):
    """Cleans up the blast results, ensures only valid blast results remain based on the treshold settings
    e.g.: minimal p-value, length of structure vs length of translation, percentage of identity"""
    # save the amount of blast results
    n_original_results = len(blastResults)
    if n_original_results == 0:
        raise NoPDBStructureFoundException("Blast to seqres provided no results: No PDB structure found")
    
    # ensure we only meet the minimal e value requirement
    filteredResults = [result for result in blastResults if result['evalue'] < MINIMAL_BLASTPE_VALUE]
    
    # ensure we only work with results that have high enough identity
    filteredResults = [result for result in filteredResults if result['pident']>= MINIMAL_TRANSLATION_TO_STRUCTURE_PIDENT_VALUE]

    # Add structural information to the filtered blast results
    appendStructureDataToBlastResults(filteredResults)
    
    # Only keep high quality resolution files
    highQualtyFilteredResults = []
    for result in filteredResults:
        if result['structure_method'] == 'x-ray diffraction':
            if not result['resolution'] is None:
                if float(result['resolution']) <= MINIMAL_XRAY_STRUCTURE_RESOLUTION:
                    highQualtyFilteredResults.append(result)
            else:
                _log.warning("Structure from '"+str(result['pdb_id'])+"' was specified as a x-ray diffraction method, but did not have a resolution")
        else:
            highQualtyFilteredResults.append(result)
    n_filtered_results = len(highQualtyFilteredResults)

    # check if there are any results left
    if n_filtered_results > 0 :
        _log.debug('From '+str(n_original_results)+' original blast results, '+str(n_filtered_results)+' remained after clean up, from executing Blast to PDB seqres Fasta')
        return highQualtyFilteredResults
    else:
        raise NoPDBStructureFoundException("Blast to seqres provided "+str(n_original_results)+" results, but not enough identity was found for the seqres sequence. With top structure: "+interpretation_to_string(blastResults[0]))


def appendStructureDataToBlastResults(blastResults):
    """Appends additional information from the pdb files given a list of blastresults"""
    for result in blastResults:
        try:
            result_structure = retrieveStructureFromPDB(result['pdb_id'])
            result['release_date'] = result_structure.header['release_date']
            result['resolution'] = result_structure.header['resolution']
            result['structure_method'] = result_structure.header['structure_method']
        except (NoPDBFileFoundException) as e:
            _log.error(e)
            result['release_date'] = None
            result['resolution'] = None
            result['structure_method'] = None
        
def retrieveBestMatchingSeqresSequence(blastResults):
    """Retrieves the top structure_result from blast results of the translated gene sequence 
    to the PDB fasta file. First selects the best structure as the one with the percentage
    of identical residues, then trims this number down by checking best resolution.
    If no resoluton is present, an NRM or cryoEM structure is chosen instead    
    """
    bestSeqres = {}
    
    # save the amount of blast results
    n_original_results = len(blastResults)
    
    # sort by the most identical percentage of residues
    sortedBlastResults = sorted(blastResults, key=lambda result: result['pident'], reverse=True)
    
    # choose the best structure
    for result in sortedBlastResults:
        if bestSeqres == {}: 
            bestSeqres = result; continue
        
        # if the percentage of identical residues are equal
        if bestSeqres['pident'] == result['pident']:
            # ensure we use the best resolution
            if result['resolution'] is None: 
                continue
            if bestSeqres['resolution'] is None: 
                bestSeqres = result; continue
            if float(bestSeqres['resolution']) > float(result['resolution']):
                bestSeqres = result; continue

    _log.debug('From '+str(n_original_results)+' cleaned blast results picked the top structure_result: '+interpretation_to_string(bestSeqres))
    
    if bestSeqres == {}:
        raise NoPDBStructureFoundException("In method retrieveBestMatchingSeqresSequence: no PDB structure was found that met the criteria for 'best' structure")
    
    return bestSeqres

def getRefseqForPDBfile(pdb_id, chain_id):
    sequence = ""
    with open(PDB_SEQRES_FASTA) as a:
        pdbfind = a.read().splitlines()
    
        n_pdb_seqres_lines = len(pdbfind)
        for i in range(n_pdb_seqres_lines):
            if pdbfind[i].startswith('>'+pdb_id+'_'+chain_id):
                if i+1 < n_pdb_seqres_lines:
                    sequence = pdbfind[i+1]
                    return sequence
                else:
                    return sequence
                
    return sequence