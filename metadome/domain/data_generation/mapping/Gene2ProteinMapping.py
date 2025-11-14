from metadome.domain.data_generation.mapping.Protein2ProteinMapping import \
    createMappingOfAASequenceToAASequence, map_single_residue
from metadome.domain.models.mapping import Mapping

import logging
from metadome.domain.models.gene import Strand

_log = logging.getLogger(__name__)

def createMappingOfGeneTranscriptionToTranslationToProtein(gene_transcription, matching_coding_translation, uniprot):
    # Retrieve the amino acid sequence according to the translation
    gene_protein_translation_sequence = matching_coding_translation['sequence']
    # Retrieve the canonical uniprot sequence
    canonical_protein_sequence = uniprot["sequence"]
    # nucleotide sequence containing of the CDS
    coding_sequence = gene_transcription['coding-sequence']
    # CDS for nucleotide positions
    cds =  gene_transcription['CDS_annotation']
        
    # check if there are multiple chromosomal positions
    chr_in_cds = list(set([cd.seqid for cd in cds]))
    if len(chr_in_cds) > 1:
        _log.warning('Found multiple chromosomes in coding sequence: '+str(chr_in_cds)+', using only '+str(chr_in_cds[0])+', this may indicate pseudoautosomal genes')
    
    # ensure we have the stop codon at the end of the translation sequence
    gene_protein_translation_sequence+='*'
    
    # align cDNA sequence with uniprot canonical sequence
    translation_to_uniprot_mapping  = createMappingOfAASequenceToAASequence(gene_protein_translation_sequence, canonical_protein_sequence)
    
    # Retrieve the codons from the cDNA translation 
    translationCodons = []
    for i in range(0, len(gene_protein_translation_sequence)): translationCodons.append(coding_sequence[i * 3 : i * 3 + 3])
    
    # Create mapping between Gene and cDNA
    cDNA_pos = 0
    currentChr = ''
    to_be_added_mapping = []
    for cd in cds:
        if currentChr == '': currentChr = cd.seqid # set it as the first
        elif not(cd.seqid == currentChr):
            _log.warning('Coding sequence with '+str(cd.seqid)+' does not match on chromosome '+str(currentChr)+', skipping this coding sequence...')
            continue # skip this Chromosome
        if cd.strand == '-':
            custom_range = range(cd.end, cd.start-1, -1)
            strand = Strand.minus
        else:
            custom_range = range(cd.start, cd.end+1)
            strand = Strand.plus
        for i in custom_range:
            # create mapping object for this position
            base_pair = coding_sequence[cDNA_pos]
            
            # increment cDNA position
            cDNA_pos = cDNA_pos+1
            aa_pos = int((cDNA_pos-1) / 3)
            
            # add codon to the mapping
            codon = translationCodons[aa_pos]                
            # add codon base pair number to the mapping
            codon_base_pair_position = ((cDNA_pos-1)%3)
            # Add residue from translation
            amino_acid_residue = gene_protein_translation_sequence[aa_pos]
            # Add information for uniprot
            if amino_acid_residue == '*':
                uniprot_residue = '*'
                uniprot_position = None
                aa_pos = None
            else:
                uniprot_position, uniprot_residue = map_single_residue(translation_to_uniprot_mapping, aa_pos)
            
            # create the mapping entry
            to_be_added_mapping.append(Mapping(
                base_pair = base_pair,
                cDNA_position = cDNA_pos,
                strand = strand,
                codon = codon,
                codon_base_pair_position = codon_base_pair_position,
                uniprot_position = uniprot_position,
                uniprot_residue = uniprot_residue,
                chromosome = str(cd.seqid),
                chromosome_position = i,
                ))
                
    return to_be_added_mapping

def extract_pdb_from_gene_region(gene_mapping, gene_region):
    """Expects a gene_region, as created by 'extract_gene_region' with 
    key/values: {'chr':chromosome, 'regions':chromosome_ranges}
    Returns: 
    {'PDB_ID': 
        {'CHAIN_ID': 
            {UNIPROT_POS: 
                {'pos': 'PDB_STRUCTURE_POS', 
                'residue': 'PDB_STRUCTURE_RESIDUE'
                }
            },
            ...
        },
        ...
    },..
    """
    pdb_structures = dict()
    not_pdb_keys = ['translation_codon_allele_pos', 'uniprot_residue', 'cDNA', 'uniprot', 'translation_codon', 'translation_residue', 'allele']
    for chr_range in gene_region['regions']:
        for i in range(chr_range[0], chr_range[1]+1):
            chr_pos = gene_region['chr']+":"+str(i)
            for key in gene_mapping['GenomeMapping']['Genome'][chr_pos].keys():
                uniprot_pos = gene_mapping['GenomeMapping']['Genome'][chr_pos]['uniprot']
                # focus only on pdb structure entries
                if key not in not_pdb_keys:
                    key_tokenized = key.split('_')
                    pdb_id = key_tokenized[0]
                    chain_id = key_tokenized[1]
    
                    if pdb_id not in pdb_structures.keys():
                        pdb_structures[pdb_id] = dict()
    
                    if chain_id not in pdb_structures[pdb_id].keys():
                        pdb_structures[pdb_id][chain_id] = dict()
                        
                    if uniprot_pos not in pdb_structures[pdb_id][chain_id].keys():
                        pdb_structures[pdb_id][chain_id][uniprot_pos] = dict()
    
                    if 'residue' in key:
                        pdb_structures[pdb_id][chain_id][uniprot_pos]['residue']= gene_mapping['GenomeMapping']['Genome'][chr_pos][key]
                        # check if the residue is identical to the swissprot residue
                        if gene_mapping['GenomeMapping']['Genome'][chr_pos][key] == gene_mapping['GenomeMapping']['Genome'][chr_pos]['uniprot_residue']:
                            pdb_structures[pdb_id][chain_id][uniprot_pos]['identical_sprot_residue'] = True
                        else:
                            pdb_structures[pdb_id][chain_id][uniprot_pos]['identical_sprot_residue'] = False
                    else:
                        pdb_structures[pdb_id][chain_id][uniprot_pos]['pos']= gene_mapping['GenomeMapping']['Genome'][chr_pos][key]
                        
    return pdb_structures
