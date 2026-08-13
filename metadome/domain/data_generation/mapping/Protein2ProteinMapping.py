import logging

from metadome.domain.wrappers.clustal import clustalw_pairwiseAlignment

_log = logging.getLogger(__name__)


# TODO(dead-data-generation): createMappingOfAASequenceToAASequence is unreachable;
#  its only importer is Gene2ProteinMapping. The clustal import at the top of this file exists
#  solely for this method.
def createMappingOfAASequenceToAASequence(primary_sequence, secondary_sequence):
    """Annotates a blast result with the atomic sequence and a mapping of the translated
     sequence based on the gene to the atomic sequence (e.g. the measured structure)"""
    # align the two sequences
    seq1, seq2 = clustalw_pairwiseAlignment(primary_sequence, secondary_sequence)
    # Create a mapping of the sequence 
    mapping_seq1_seq2 = createAlignedSequenceMapping(seq1, seq2)
    
    return {'mapping':mapping_seq1_seq2, 'primary_sequence':seq1, 'secondary_sequence':seq2}

# TODO(dead-data-generation): map_single_residue is unreachable;
#  its only importer is Gene2ProteinMapping. The clustal import at the top of this file exists
#  solely for this method.
def map_single_residue(aligned_mapping, cur_protein_position, alternate_position_mapping=None):
    """Maps a single residue for a previously mapped sequence at a given position"""
    # check if the number is in the mapping
    if cur_protein_position in aligned_mapping['mapping'].keys():
        # add the position
        mapped_position = aligned_mapping['mapping'][cur_protein_position]
        residue_at_mapped_position = aligned_mapping['secondary_sequence'][mapped_position]
        
        # if an alternative mapping is provided
        if not alternate_position_mapping is None:
            # remap the sequence position to the mapped position provided
            mapped_position = alternate_position_mapping[mapped_position]
            
    else:
        mapped_position ='-'
        residue_at_mapped_position = '-'
    
    return mapped_position, residue_at_mapped_position

def createAlignedSequenceMapping(seq1, seq2, strict=True):
    """Creates a position mapping from seq1 to seq2 for two aligned sequences.
    If strict is True, non-matching residues will be ignored in the mapping
    If strict is False only gaps (i.e. '-') will be ignored."""    
    try:
        assert len(seq1) == len(seq2)
    except AssertionError as e:
        _log.error("{}".format(e.output))
        return {}
    
    if strict:
        return {i:i for i, v_i in enumerate(seq1) if v_i == seq2[i] and v_i != '-'}
    else:
        return {i:i for i, v_i in enumerate(seq1) if v_i != '-' and seq2[i] != '-'}
