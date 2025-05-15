import logging
import tempfile
import subprocess
import os
from metadome.default_settings import BLASTP_EXECUTABLE

_log = logging.getLogger(__name__)

# The supported format specifiers are:
#          qseqid means Query Seq-id
#             qgi means Query GI
#            qacc means Query accesion
#         qaccver means Query accesion.version
#            qlen means Query sequence length
#          sseqid means Subject Seq-id
#       sallseqid means All subject Seq-id(s), separated by a ';'
#             sgi means Subject GI
#          sallgi means All subject GIs
#            sacc means Subject accession
#         saccver means Subject accession.version
#         sallacc means All subject accessions
#            slen means Subject sequence length
#          qstart means Start of alignment in query
#            qend means End of alignment in query
#          sstart means Start of alignment in subject
#            send means End of alignment in subject
#            qseq means Aligned part of query sequence
#            sseq means Aligned part of subject sequence
#          evalue means Expect value
#        bitscore means Bit score
#           score means Raw score
#          length means Alignment length
#          pident means Percentage of identical matches
#          nident means Number of identical matches
#        mismatch means Number of mismatches
#        positive means Number of positive-scoring matches
#         gapopen means Number of gap openings
#            gaps means Total number of gaps
#            ppos means Percentage of positive-scoring matches
#          frames means Query and subject frames separated by a '/'
#          qframe means Query frame
#          sframe means Subject frame
#            btop means Blast traceback operations (BTOP)
#         staxids means unique Subject Taxonomy ID(s), separated by a ';'
#                       (in numerical order)
#       sscinames means unique Subject Scientific Name(s), separated by a
#       scomnames means unique Subject Common Name(s), separated by a ';'
#      sblastnames means unique Subject Blast Name(s), separated by a ';'
#                       (in alphabetical order)
#      sskingdoms means unique Subject Super Kingdom(s), separated by a ';
#                       (in alphabetical order)
#          stitle means Subject Title
#      salltitles means All Subject Title(s), separated by a '<>'
#         sstrand means Subject Strand
#           qcovs means Query Coverage Per Subject
#         qcovhsp means Query Coverage Per HSP
# When not provided, the default value is:
# 'qseqid sseqid pident length mismatch gapopen qstart qend sstart send
# evalue bitscore', which is equivalent to the keyword 'std'

def run_blast(sequence, blastdb, max_target_seqs=None):
    """ Default output:
    qseqid - Query Seq-id,
    sseqid - Subject Seq-id,
    lident - length of identical matches,
    pident - Percentage of identical matches,
    nident - Number of identical matches,
    length - Alignment length,
    mismatch - Number of mismatches,
    gapopen - Number of gap openings,
    qstart - Start of alignment in query,
    qend - End of alignment in query,
    sstart - Start of alignment in subject,
    send - End of alignment in subject,
    evalue - Expect value, The E-value is a parameter that describes the number of hits one can expect to see by chance when searching a database of a particular size. It decreases exponentially with the score (S) that is assigned to a match between two sequences. Essentially, the E-value describes the random background noise that exists for matches between sequences. For example, an E-value of 1 assigned to a hit can be interpreted as meaning that in a database of the current size, one might expect to see one match with a similar score simply by chance. This means that the lower the E-value, or the closer it is to 0, the higher is the significance of the match. However, it is important to note that searches with short sequences can be virtually identical and have relatively high E-value. This is because the calculation of the E-value also takes into account the length of the query sequence. This is because shorter sequences have a high probability of occurring in the database purely by chance. 
    bitscore - The bit score gives an indication of how good the alignment is; the higher the score, the better the alignment. In general terms, this score is calculated from a formula that takes into account the alignment of similar or identical residues, as well as any gaps introduced to align the sequences.
    """
    
    with tempfile.NamedTemporaryFile(suffix=".fasta", delete=False) as tmp_file:
        tmp_file.write(sequence.encode('utf-8'))
    out_blast = tmp_file.name + '.blastp'
    
    if max_target_seqs is None:
        args = [BLASTP_EXECUTABLE, "-query", tmp_file.name, "-evalue", "1e-5",
            "-num_threads", "15", "-db", blastdb,
            "-out", out_blast, '-outfmt', '10 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore nident']
    else:
        args = [BLASTP_EXECUTABLE, "-query", tmp_file.name, "-evalue", "1e-5",
            "-num_threads", "15", "-db", blastdb,
            "-out", out_blast, '-outfmt', '10 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore nident',
            "-max_target_seqs", str(max_target_seqs)]
    try:
        subprocess.call(args)
    except subprocess.CalledProcessError as e:
        _log.error("{}".format(e.output))

    if os.path.exists(out_blast):
        with open(out_blast) as a:
            output = a.read().splitlines()
        logging.debug(output)
        if output:
            logging.debug('output True')
        os.remove(out_blast)
    else:
        output = []
    os.remove(tmp_file.name)
    
    return output

def interpret_blast_as_uniprot(line):
    """ return a blast result as a dictionary with
    key - value
    database_id - sp for swissprot, 
    accession_code - accession code,
    entry_name - entry name
    qseqid - Query Seq-id,
    sseqid - Subject Seq-id,
    lident - length of identical matches,
    pident - Percentage of identical matches,
    nident - Number of identical matches,
    length - Alignment length,
    mismatch - Number of mismatches,
    gapopen - Number of gap openings,
    qstart - Start of alignment in query,
    qend - End of alignment in query,
    sstart - Start of alignment in subject,
    send - End of alignment in subject,
    evalue - Expect value,
    bitscore - Bit score    
    """
    result = interpret_blast(line)
    
    # convert the id to a database id, an accession code and an entry name
    result['database_id'] = result['sseqid'].split('|')[0]
    result['accession_code'] = result['sseqid'].split('|')[1]
    result['entry_name'] = result['sseqid'].split('|')[2]
    
    
    return result

def interpret_blast_as_pdb(line):
    """ returns a blast result as a dictionary with
    key - value
    pdb_id - The PDB id
    chain_id - The chain id
    qseqid - Query Seq-id,
    sseqid - Subject Seq-id,
    lident - length of identical matches,
    pident - Percentage of identical matches,
    nident - Number of identical matches,
    length - Alignment length,
    mismatch - Number of mismatches,
    gapopen - Number of gap openings,
    qstart - Start of alignment in query,
    qend - End of alignment in query,
    sstart - Start of alignment in subject,
    send - End of alignment in subject,
    evalue - Expect value,
    bitscore - Bit score    
    """
    result = interpret_blast(line)
    
    # convert the id to a pdb_id and a chain id
    result['pdb_id'] = result['sseqid'].split('_')[0]
    result['chain_id'] = result['sseqid'].split('_')[1]
    
    return result
    

def interpret_blast(line):
    """ returns a blast result as a dictionary with
    key - value
    qseqid - Query Seq-id,
    sseqid - Subject Seq-id,
    lident - length of identical matches,
    pident - Percentage of identical matches,
    nident - Number of identical matches,
    length - Alignment length,
    mismatch - Number of mismatches,
    gapopen - Number of gap openings,
    qstart - Start of alignment in query,
    qend - End of alignment in query,
    sstart - Start of alignment in subject,
    send - End of alignment in subject,
    evalue - Expect value,
    bitscore - Bit score    
    """

    result_line = line.split(',')

    result = {'qseqid': result_line[0],
              'sseqid': result_line[1],
              'pident': float(result_line[2]),
              'length': int(result_line[3]),
              'mismatch': int(result_line[4]),
              'gapopen': int(result_line[5]),
              'qstart': int(result_line[6]),
              'qend': int(result_line[7]),
              'sstart': int(result_line[8]),
              'send': int(result_line[9]),
              'evalue': float(result_line[10]),
              'bitscore': float(result_line[11]),
              'nident': int(result_line[12]),
              }
    
    result['lident'] = (result['send']-result['sstart'])*(result['pident']/100)
    
    return result

def interpretation_to_string(interpretation):
    return 'ID:'+str(interpretation['sseqid'])+', identity '+str(interpretation['pident'])+'%, identity length '+str(interpretation['lident'])+' total length:'+str(interpretation['length'])+', bit score:'+str(interpretation['bitscore'])+', E-value:'+str(interpretation['evalue'])