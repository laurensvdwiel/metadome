# TODO(dead-data-generation): unreachable at runtime.
#  The Pfam HMM binaries are used only by generate_pfam_alignments,
#  called from meta_domain_creation and so from install.py, which no
#  longer exists. The live meta-domain path is retrieve_pfam_aligned_codons,
#  which reads the database.
import logging
from metadome.default_settings import PFAM_HMM_DAT, HMMFETCH_EXECUTABLE, PFAM_HMM,\
    HMMALIGN_EXECUTABLE, HMMEMIT_EXECUTABLE, HMMSTAT_EXECUTABLE,\
    HMMLOGO_EXECUTABLE, PFAM_ALIGNMENT_DIR
from metadome.domain.parsers.fasta import unwrap_fasta_alignment
from builtins import FileNotFoundError
import os
import tempfile
import gzip
import subprocess
import re
import errno

_log = logging.getLogger(__name__)

class FoundNoPfamHMMException(Exception):
    pass

class FoundMoreThanOnePfamHMMException(Exception):
    pass

def retrieve_PFAM_ID_by_AC(pfam_ac):
    """Retrieves the PFAM ID w.r.t. the provided PFAM Accession Code
    The PFAM ID is retrieved from the PFAM HMM metadata file, thus
    ensuring a HMM seed is present for the provided AC"""
    openFunc = gzip.open if PFAM_HMM_DAT.endswith(".gz") else open
    with openFunc(PFAM_HMM_DAT) as infile:
        lines = infile.readlines()
        for line_index in range(len(lines)):
            line = lines[line_index].decode()
            if line.startswith("#=GF AC   "+pfam_ac):
                AC = line[10:].split('\n')[0]
                ID =  lines[line_index-1].decode()[10:].split('\n')[0]
                yield ID, AC

def retrieve_PFAM_full_alignment_by_AC(pfam_ac):
    """
    retrieves and parses the full alignment for a 
     Pfam full alignment file in the following format:
    'AC' : pfam accession code
    'alignments' : all the alignments for this Pfam
        accession code,
    'consensus' : 
        {
          'identifier' : consensus identifier,
          'sequence' : consensus sequence
        },
    SS_cons (optional):
        Mark up of sequence of known structure for 
        secondary structure.
    seq_cons, consensus sequence (i.e 60% or above,
        of the amino acids in this column belong to
        this class of residue):
        class           key  residues
        A               A    A
        C               C    C
        D               D    D
        E               E    E
        F               F    F
        G               G    G
        H               H    H
        I               I    I
        K               K    K
        L               L    L
        M               M    M
        N               N    N
        P               P    P
        Q               Q    Q
        R               R    R
        S               S    S
        T               T    T
        V               V    V
        W               W    W
        Y               Y    Y
        alcohol         o    S,T
        aliphatic       l    I,L,V
        any             .    A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y
        aromatic        a    F,H,W,Y
        charged         c    D,E,H,K,R
        hydrophobic     h    A,C,F,G,H,I,K,L,M,R,T,V,W,Y
        negative        -    D,E
        polar           p    C,D,E,H,K,N,Q,R,S,T
        positive        +    H,K,R
        small           s    A,C,D,G,N,P,S,T,V
        tiny            u    A,G,S
        turnlike        t    A,C,D,E,G,H,K,N,Q,R,S,T
    """
    
    # Construct the filepath of the alignment file
    alignment_file = PFAM_ALIGNMENT_DIR+pfam_ac+"/full"
    
    # check if the path exists
    if not os.path.exists(alignment_file):
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), alignment_file)
    
    # check if there is an ID associated with the provided accession code
    pfam_hmm_meta_data = retrieve_PFAM_ID_by_AC(pfam_ac)
    found_ids = 0
    for pfam_id, _ in pfam_hmm_meta_data:
        found_ids+=1
    
    if found_ids > 1:
        raise FoundMoreThanOnePfamHMMException("Found more than one Pfam ids that match the Pfam ac '"+pfam_ac+"' when searching for matching HMMER HMM's")
    if found_ids == 0:
        raise FoundNoPfamHMMException("Found no matching Pfam ids that for Pfam ac '"+pfam_ac+"' when searching for matching HMMER HMM's")
    
    # create a temporary HMM file
    tmp_hmm_file = tempfile.NamedTemporaryFile(suffix=".hmm", delete=False)
    
    # fetch the HMM
    fetch_args = [HMMFETCH_EXECUTABLE, "-o", tmp_hmm_file.name, PFAM_HMM, pfam_id]
    try:
        subprocess.call(fetch_args)
    except subprocess.CalledProcessError as e:
        _log.error("{}".format(e.output))
    consensus_identifier, consensus_sequence = retrieve_consensus_sequence(tmp_hmm_file.name)
    
    pfam_alignment_output = {'AC':pfam_ac,
                            'alignments':[], 
                            'consensus':{'identifier':consensus_identifier,
                                        'sequence':consensus_sequence},
                            'SS_cons':None, 
                            'seq_cons':None}
    
    
    # interpret the full alignments input for this domain
    n_sequences = 0
    try:
        with open(alignment_file) as a:
            # retrieve the aligned sequences from the output file
            Pfam_alignments = a.readlines()
            for line in Pfam_alignments:
                if not line.startswith("#") and not line.startswith("//"):
                    # An alignment row.
                    alignment_row = [al for al in line.strip().split(" ") if len(al)>0]
                    alignment = alignment_row[1]
                    alignment_consensus = convert_pfam_fasta_alignment_to_strict_consensus_sequence(alignment)
                    sequence_origin = alignment_row[0].split("/")
                    sequence_identifier =  sequence_origin[0]
                    sequence_region = sequence_origin[1].split('-')
                    sequence_start = int(sequence_region[0])
                    sequence_stop = int(sequence_region[1])
                    
                    # add the alignment to the output
                    pfam_alignment_output['alignments'].append({'seq_id':sequence_identifier, 
                                                                'seq_start':sequence_start, 
                                                                'seq_stop':sequence_stop, 
                                                                'alignment':alignment,
                                                                'alignment_consensus':alignment_consensus})
                elif line.startswith("#=GC SS_cons"):
                    pfam_alignment_output['SS_cons'] = [al for al in line.strip().split(" ") if len(al)>0][1]
                elif line.startswith("#=GC seq_cons"):
                    pfam_alignment_output['seq_cons'] = [al for al in line.strip().split(" ") if len(al)>0][1]
                elif line.startswith("#=GF SQ"):
                    # Number of sequences, start of alignment.
                    n_sequences = int([al for al in line.strip().split(" ") if len(al)>0][2])
    except IOError as e:
        _log.error("{}".format(e.output))
    
    if n_sequences != len(pfam_alignment_output['alignments']):
        _log.error("Number of sequences in full alignment for Pfam "+pfam_ac+" ('"+n_sequences+"') does not match the actual number of sequences '"+len(pfam_alignment_output['alignments'])+"'")

    # remove the temporary files
    os.remove(tmp_hmm_file.name)
    
    # return the output
    return pfam_alignment_output
        
        
def retrieve_consensus_sequence(hmm_file):
    """For a given HMM file ('hmm_file'), this method returns the 
    majority-rule consensus sequence that and the consensus identifier"""
    # create the temporary file wherein the consensus will be saved
    tmp_consensus_sequence_file = tempfile.NamedTemporaryFile(suffix=".fasta", delete=False)
    consensus_identifier = None
    consensus_sequence = None
    try:
        # check if the provided file exists
        if os.path.exists(hmm_file):    
            
            # emit the consensus sequence of the HMM
            emit_args = [HMMEMIT_EXECUTABLE, "-o", tmp_consensus_sequence_file.name, "-C", hmm_file]
            try:
                subprocess.call(emit_args)
            except subprocess.CalledProcessError as e:
                _log.error("{}".format(e.output))
            
            # Process the output of the hmmemit call
            if os.path.exists(tmp_consensus_sequence_file.name):
                with open(tmp_consensus_sequence_file.name) as a:
                    aln = unwrap_fasta_alignment(a.read().splitlines())
                
                    # Perform checks on the consensus fasta file:
                    if not(len(aln) == 2) or not(aln[0].startswith(">")):
                        raise Exception("something went wrong when finding the consensus")
                        
                    consensus_identifier = aln[0][1:]
                    consensus_sequence = aln[1]
                    
    except IOError as e:
        _log.error("{}".format(e.output))        
    os.remove(tmp_consensus_sequence_file.name)
    
    if consensus_sequence is None:
        raise Exception("No consensus found")
    
    # return the identifier and consensus
    return consensus_identifier, consensus_sequence

def report_hmm_logo_for_pfam(pfam_ac):
    """Creates HMM logo from a pfam HMM based on the pfam accession
    code and generates dictionaries based on the output of a hmmlogo
    with the following keys:
     
      pos
        The position of the residue within the consensus sequence,
        starting from 1.
      total_relent
        the relative entropy at this position
      A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y
        relent for the specific amino acid residue at this position
    """
    # check if there is an ID associated with the provided accession code
    pfam_hmm_meta_data = retrieve_PFAM_ID_by_AC(pfam_ac)
    found_ids = 0
    for pfam_id, _ in pfam_hmm_meta_data:
        found_ids+=1
    
    if found_ids > 1:
        raise FoundMoreThanOnePfamHMMException("Found more than one Pfam ids that match the Pfam ac '"+pfam_ac+"' when searching for matching HMMER HMM's")
    if found_ids == 0:
        raise FoundNoPfamHMMException("Found no matching Pfam ids that for Pfam ac '"+pfam_ac+"' when searching for matching HMMER HMM's")
    
    # create a temporary HMM file
    tmp_hmm_file = tempfile.NamedTemporaryFile(suffix=".hmm", delete=False)
    
    # fetch the HMM
    fetch_args = [HMMFETCH_EXECUTABLE, "-o", tmp_hmm_file.name, PFAM_HMM, pfam_id]
    try:
        subprocess.call(fetch_args)
    except subprocess.CalledProcessError as e:
        _log.error("{}".format(e.output))
    
    # Create temp file for storing the hmm logo
    tmp_hmm_logo_file = tempfile.NamedTemporaryFile(suffix=".hmm_hmmlogo", delete=False)
    
    # use the HMM as a hmmalign using the found HMM model    
    hmmlogo_args = [HMMLOGO_EXECUTABLE, "--height_relent_all", "--no_indel", tmp_hmm_file.name]
    try:
        subprocess.call(hmmlogo_args, stdout=tmp_hmm_logo_file)
    except subprocess.CalledProcessError as e:
        _log.error("{}".format(e.output))
    
    # interpret the hmm logo
    hmmlogo_output = []
    try:
        # retrieve the aligned sequences from the output file
        if os.path.exists(tmp_hmm_logo_file.name):
            with open(tmp_hmm_logo_file.name) as a:
                hmm_logo_lines = a.readlines()
                for line in hmm_logo_lines:
                    if not line[0].isdigit(): continue # empty strings are false in python
                    tokens = line.split()
                    hmmlogo_entry = {\
                                     "pos":int(tokens[0].split(':')[0]),\
                                     "A":float(tokens[1]),\
                                     "C":float(tokens[2]),\
                                     "D":float(tokens[3]),\
                                     "E":float(tokens[4]),\
                                     "F":float(tokens[5]),\
                                     "G":float(tokens[6]),\
                                     "H":float(tokens[7]),\
                                     "I":float(tokens[8]),\
                                     "K":float(tokens[9]),\
                                     "L":float(tokens[10]),\
                                     "M":float(tokens[11]),\
                                     "N":float(tokens[12]),\
                                     "P":float(tokens[13]),\
                                     "Q":float(tokens[14]),\
                                     "R":float(tokens[15]),\
                                     "S":float(tokens[16]),\
                                     "T":float(tokens[17]),\
                                     "V":float(tokens[18]),\
                                     "W":float(tokens[19]),\
                                     "Y":float(tokens[20]),\
                                     "total_relent":float(tokens[22].split(')')[0]),\
                                    }
                    hmmlogo_output.append(hmmlogo_entry)
    except IOError as e:
        _log.error("{}".format(e.output))
    
    # remove the temporary files
    os.remove(tmp_hmm_file.name)
    os.remove(tmp_hmm_logo_file.name)
    
    # return the hmmstat
    return hmmlogo_output
    
def report_hmm_stat_for_pfam(pfam_ac):
    """Retrieves a pfam HMM based on pfam accession code and generates 
    dictionaries based on the output of a hmmstat report with the 
    following keys:
     
      idx
        The index of this profile, numbering each on in the file
        starting from 1.
      name
        The name of the profile.
      accession
        The optional accession of the profile, or "-" if there is none.
      nseq
        The number of sequences that the profile was estimated from.
      eff_nseq
        The effective number of sequences that the profile was estimated
        from, after HMMER applied an effective sequence number calculation
        such as the default entropy weighting.
      M
        The length of the model in consensus residues (match states).
      relent
        Mean relative entropy per match state, in bits. This is the
        expected (mean) score per consensus position. This is what the
        default entropy-weighting method for effective sequence number
        estimation focuses on, so for default HMMER3 models, you expect
        this value to reflect the default target for entropy-weighting.
      info
        Mean information content per match state, in bits. Probably not
        useful. Information content is a slightly different calculation
        than relative entropy.
      p_relE
        Mean positional relative entropy, in bits. This is a fancier
        version of the per-match-state relative entropy, taking into
        account the transition (insertion/deletion) probabilities; it
        may be a more accurate estimation of the average score
        contributed per model consensus position.
      compKL
        Kullback-Leibler distance between the model's overall average
        residue composition and the default background frequency
        distribution. The higher this number, the more biased the
        residue composition of the profile is. Highly biased profiles
        can slow the HMMER3 acceleration pipeline, by causing  too  many
        nonhomologous sequences to pass the filters.
    """
    # check if there is an ID associated with the provided accession code
    pfam_hmm_meta_data = retrieve_PFAM_ID_by_AC(pfam_ac)
    found_ids = 0
    for pfam_id, _ in pfam_hmm_meta_data:
        found_ids+=1
    
    if found_ids > 1:
        raise FoundMoreThanOnePfamHMMException("Found more than one Pfam ids that match the Pfam ac '"+pfam_ac+"' when searching for matching HMMER HMM's")
    if found_ids == 0:
        raise FoundNoPfamHMMException("Found no matching Pfam ids that for Pfam ac '"+pfam_ac+"' when searching for matching HMMER HMM's")
    
    # create a temporary HMM file
    tmp_hmm_file = tempfile.NamedTemporaryFile(suffix=".hmm", delete=False)
    
    # fetch the HMM
    fetch_args = [HMMFETCH_EXECUTABLE, "-o", tmp_hmm_file.name, PFAM_HMM, pfam_id]
    try:
        subprocess.call(fetch_args)
    except subprocess.CalledProcessError as e:
        _log.error("{}".format(e.output))
        
    
    # Create temp file for storing the statistics
    tmp_hmm_stat_file = tempfile.NamedTemporaryFile(suffix=".hmm_hmmstats", delete=False)
    
    # use the HMM as a hmmalign using the found HMM model    
    hmmstat_args = [HMMSTAT_EXECUTABLE, tmp_hmm_file.name]
    try:
        subprocess.call(hmmstat_args, stdout=tmp_hmm_stat_file)
    except subprocess.CalledProcessError as e:
        _log.error("{}".format(e.output))
    
    # interpret the statistics from thePfam's HMM
    hmmstat_output = []
    try:
        # retrieve the aligned sequences from the output file
        if os.path.exists(tmp_hmm_stat_file.name):
            with open(tmp_hmm_stat_file.name) as a:
                Pfam_alignments = a.readlines()
                for line in Pfam_alignments:
                    if line.startswith('#'): continue
                    if not line.strip(): continue # empty strings are false in python
                    tokens = line.split()
                    hmmstat_entry = {\
                                     "idx":int(tokens[0]),\
                                     "name":tokens[1],\
                                     "accession":tokens[2],\
                                     "nseq":int(tokens[3]),\
                                     "eff_nseq":float(tokens[4]),\
                                     "M":int(tokens[5]),\
                                     "relent":float(tokens[6]),\
                                     "info":float(tokens[7]),\
                                     "p_relE":float(tokens[8]),\
                                     "compKL":float(tokens[9]),\
                                    }
                    assert len( hmmstat_entry ) == len(tokens)
                    hmmstat_output.append(hmmstat_entry)
    except IOError as e:
        _log.error("{}".format(e.output))
    
    # remove the temporary files
    os.remove(tmp_hmm_file.name)
    os.remove(tmp_hmm_stat_file.name)
    
    # return the hmmstat
    return hmmstat_output

def create_alignment_of_sequences_according_to_PFAM_HMM(sequences, pfam_ac, target_directory, target_file_alignments):
    """ aligns homologues pfam domains contained in 'sequences' based 
    on the Pfam-HMM specified by the 'pfam_ac'.
    
    The assumption of this method is that the 'sequences' are found to
    be all part of the domain specified by 'pfam_ac' (_src Pfam accession
    code ex.: 'PF00102')
    
    Output: _src dictionary containing:
    KEY - VALUE DESCRIPTION
    'alignments' - _src python list type containing python dictionaries
                    corresponding to the amount of sequences provided
                    by 'sequences'. The key 'seq_nr' herein corresponds
                    to the index of the sequence in 'sequences' and the
                    key 'alignment' corresponds to the alignment created
                    in the context of the Pfam HMM w.r.t. the other
                    sequences in 'sequences'
   'consensus' - _src python dictionary containing the consensus_identifier,
                    consensus_sequence and alignments with corresponding
                    keys
    'PP_cons' - consensus posterior probability annotation for the entire
                    column. It’s calculated simply as the arithmetic mean
                    of the per-residue posterior probabilities in that
                    column.
    'RF' - reference coordinate annotation, with an x marking each column
                    that the profile considered to be consensus.
    
    """
    # check if there is an ID associated with the provided accession code
    pfam_hmm_meta_data = retrieve_PFAM_ID_by_AC(pfam_ac)
    found_ids = 0
    for pfam_id, pfam_ac_hmm in pfam_hmm_meta_data:
        found_ids+=1
    
    if found_ids > 1:
        raise FoundMoreThanOnePfamHMMException("Found more than one Pfam ids that match the Pfam ac '"+pfam_ac+"' when searching for matching HMMER HMM's")
    if found_ids == 0:
        raise FoundNoPfamHMMException("Found no matching Pfam ids that for Pfam ac '"+pfam_ac+"' when searching for matching HMMER HMM's")
    
    # create _src temporary HMM file
    tmp_hmm_file = tempfile.NamedTemporaryFile(suffix=".hmm", delete=False)
    
    # fetch the HMM
    fetch_args = [HMMFETCH_EXECUTABLE, "-o", tmp_hmm_file.name, PFAM_HMM, pfam_id]
    try:
        subprocess.call(fetch_args)
    except subprocess.CalledProcessError as e:
        _log.error("{}".format(e.output))
    
    # get the consensus sequence
    consensus_identifier, consensus_sequence = retrieve_consensus_sequence(tmp_hmm_file.name)
    
    # create an fasta file from the various sequences
    with tempfile.NamedTemporaryFile(suffix=".fasta", delete=False) as tmp_sequences_file:
        tmp_sequences_file.write(('\n'.join(['>'+str(consensus_identifier), consensus_sequence])).encode(encoding='utf_8', errors='strict'))
        for sequence in sequences:
            tmp_sequences_file.write(('\n'.join(['>'+str(sequence['uniprot_ac'])+'/'+str(sequence['start'])+'-'+str(sequence['stop']), sequence['sequence']])).encode(encoding='utf_8', errors='strict'))
    
    # Create temp file for storing the alignment
    tmp_hmm_alignment_file = tmp_hmm_file.name+"_alignment"
    
    # use the HMM as _src hmmalign using the found HMM model    
    hmmalign_args = [HMMALIGN_EXECUTABLE, "-o", tmp_hmm_alignment_file, "--outformat", "Pfam", tmp_hmm_file.name, tmp_sequences_file.name]
    try:
        subprocess.call(hmmalign_args)
    except subprocess.CalledProcessError as e:
        _log.error("{}".format(e.output))
    
    try:
        # retrieve the aligned sequences from the output file
        if os.path.exists(tmp_hmm_alignment_file):
            # first create the metadomain dir if it does not yet exist
            if not os.path.isdir(target_directory):
                _log.info('Directory '+target_directory+' did not exist yet, creating ...')
                os.mkdir(target_directory)
                
            # first create the specific metadomain dir if it does not yet exist
            if not os.path.isdir(target_directory+pfam_ac):
                _log.info('Directory '+target_directory+pfam_ac+' did not exist yet, creating ...')
                os.mkdir(target_directory+pfam_ac)
                
            with open(tmp_hmm_alignment_file) as _src:
                with open(target_directory+pfam_ac+'/'+target_file_alignments, 'wt') as _dst:
                    Pfam_alignments = _src.readlines()
                    for line in Pfam_alignments:
                        # write the line
                        _dst.write(line)
                        if line.startswith('# STOCKHOLM 1.0'):
                            # this is just after the start of the file, appending comments
                            _dst.write('#=GF ID '+consensus_identifier+'\n')
                            _dst.write('#=GF AC '+pfam_ac+'\n')
                            _dst.write('#=GF DC This alignment file only contains Pfam domains found in the human species'+'\n')
                            _dst.write('#=GF CC consensus_sequence:'+consensus_sequence+'\n')
    except IOError as e:
        _log.error("{}".format(e.output))
        
    # remove the temporary files
    os.remove(tmp_hmm_file.name)
    os.remove(tmp_sequences_file.name)
    os.remove(tmp_hmm_alignment_file)

def interpret_hmm_alignment_file(hmm_alignment_file):
    # interpret the alignments made by Pfam's HMM
    hmmalign_output = {}
    try:
        # retrieve the aligned sequences from the output file
        if os.path.exists(hmm_alignment_file):
            hmmalign_output['alignments'] = []
            hmmalign_output['consensus'] = {'identifier':None, 'alignment':None, 'sequence':None}
            hmmalign_output['PP_cons'] = None
            hmmalign_output['RF'] = None
            
            with open(hmm_alignment_file) as a:
                Pfam_alignments = a.readlines()
                processed_comments = False
                for line in Pfam_alignments:
                    if processed_comments:
                        if line.startswith(hmmalign_output['consensus']['identifier']):
                            # this is the consensus sequence
                            alignment = [al for al in line.strip().split(" ") if len(al)>0]
                            hmmalign_output['consensus']['alignment'] = alignment[1]
                        elif line.startswith("#=GC PP_cons"):
                            # handle the consensus posterior probability annotation for the entire column
                            pp_cons = [al.strip() for al in line.strip().split("#=GC PP_cons") if len(al)>0]
                            hmmalign_output['PP_cons'] = pp_cons[0]
                        elif line.startswith("#=GC RF"):
                            # handle the reference coordinate annotation
                            rf = [al.strip() for al in line.strip().split("#=GC RF") if len(al)>0]
                            hmmalign_output['RF'] = rf[0]
                        else:
                            alignment = [al for al in line.strip().split(" ") if len(al)>0]
                            seq_identifier = alignment[0]
                            if '/' in seq_identifier and '-' in seq_identifier:
                                # this is one of the alignments
                                uniprot_ac_and_positions = seq_identifier.split('/')
                                uniprot_ac = uniprot_ac_and_positions[0]
                                positions = uniprot_ac_and_positions[1].split('-')
                                start_pos = positions[0]
                                stop_pos = positions[1] 
                                
                                hmmalign_output['alignments'].append({'seq_id': seq_identifier,'uniprot_ac':uniprot_ac, 'start_pos':start_pos, 'stop_pos':stop_pos, 'alignment':alignment[1]})
                            
                    else:
                        if line.startswith('#=GF ID '):
                            hmmalign_output['consensus']['identifier'] = line.strip().split('#=GF ID ')[1]
                        elif line.startswith('#=GF CC consensus_sequence:'):
                            hmmalign_output['consensus']['sequence'] = line.strip().split('#=GF CC consensus_sequence:')[1]
                        elif line == '\n':
                            processed_comments = True
    except IOError as e:
        _log.error("{}".format(e.output))
        
    
    # return the alignments
    return hmmalign_output

def map_sequence_to_aligned_sequence(original_sequence, original_sequence_aligned):
    
    # ensure we are dealing with the same sequence
    assert original_sequence == re.sub('-', '', original_sequence_aligned)

    mapping_sequence_to_aligned_sequence = {}
    domain_sequence_position = 0
    for alignment_i, value_i in enumerate(original_sequence_aligned):
        if value_i != '-':
            mapping_sequence_to_aligned_sequence[alignment_i] = domain_sequence_position
            
            # double check the sequence mapping
            assert value_i == original_sequence[domain_sequence_position]
            domain_sequence_position += 1
            
    return mapping_sequence_to_aligned_sequence

def convert_pfam_fasta_alignment_to_original_aligned_sequence(alignment_sequence):
    """Removes any '.', and converts lower case characters to upper case in the alignment
    sequence"""
    return re.sub(r'[.]', '-', alignment_sequence).upper()


def convert_pfam_fasta_alignment_to_strict_fasta(alignment_sequence):
    """Converts any '.' characters or lower case characters to the fasta
    gap character '-'.
    From the HMMER manual: a ’-’ character means a deletion relative to
    the consensus. In an insert column, residues are lower case, and a
    ’.’ is padding."""
    return re.sub(r'[a-z,.]', '-', alignment_sequence)

def convert_pfam_fasta_alignment_to_strict_sequence(alignment_sequence):
    """Removes any '.', '-' or lower case characters in the alignment
    sequence in order to obtain the strict fasta sequence"""
    return re.sub(r'[a-z,.-]', '', alignment_sequence)

def convert_pfam_fasta_alignment_to_strict_consensus_sequence(alignment_sequence):
    """Removes any '.', '-' or lower case characters in the alignment
    sequence in order to obtain the strict fasta sequence"""
    return re.sub(r'[a-z,.]', '', alignment_sequence)