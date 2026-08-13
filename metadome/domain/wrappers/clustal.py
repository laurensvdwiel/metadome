# TODO(dead-data-generation): unreachable at runtime.
#  Reached only via mapping_generator -> database_creation -> install.py, which no longer exists.
import logging
import tempfile
import subprocess
import os
from Bio.Align.Applications._Clustalw import ClustalwCommandline
from metadome.default_settings import CLUSTALW_EXECUTABLE
from metadome.domain.parsers.fasta import unwrap_fasta_alignment

_log = logging.getLogger(__name__)

def clustalw_pairwiseAlignment(seq1, seq2):
    """Creates a pairwise alignment for two given sequences.
    Returns the aligned sequences"""
    # create a temporary file
    tmp_file = tempfile.NamedTemporaryFile(suffix=".fasta", delete=False)
    with tmp_file as f:
        f.write(('\n'.join(['>1', seq1, '>2', seq2])).encode(encoding='utf_8', errors='strict'))
    
    aln_seq1 = ""
    aln_seq2 = ""
    
    # construct the clustalw command
    cline = ClustalwCommandline(CLUSTALW_EXECUTABLE, infile=tmp_file.name, outfile=tmp_file.name+"_al", output="FASTA", outorder="input")
    _log.debug("Running command: "+str(cline))    
    try:        
        # run clustalw command
        stdout, stderr = cline()
        
        # log any output
        if stdout:
            _log.debug(stdout)
        # log any errors
        if stderr:
            _log.error(stderr)
            
        try:
            # retrieve the aligned sequences from the output file
            if os.path.exists(tmp_file.name + '_al'):
                with open(tmp_file.name + '_al') as a:
                    aln = unwrap_fasta_alignment(a.read().splitlines())
                aln_seq1 = aln[1]
                aln_seq2 = aln[3]
        except IOError as e:
            _log.error("{}".format(e.output))
    
    except subprocess.CalledProcessError as e:
        _log.error("{}".format(e.output))
        
    # remove any created files
    if os.path.exists(tmp_file.name + '_al'):
        os.remove(tmp_file.name + '_al')
    if os.path.exists(tmp_file.name):
        os.remove(tmp_file.name)
    
    # return the aligned sequences
    return aln_seq1, aln_seq2