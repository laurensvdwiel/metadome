from metadome.domain.metrics.codon_statistics import codon_background_rates
from metadome.domain.models.entities.single_nucleotide_variant import SingleNucleotideVariant

class ExternalREFAlleleNotEqualsTranscriptionException(Exception):
    pass

def retrieve_background_variant_counts(mappings_per_chromosome):
    # the value that is to be returned
    variant_type_counts = dict()
    
    for chrom_pos in mappings_per_chromosome.keys():
        codon = mappings_per_chromosome[chrom_pos]['base_pair_representation']
        residue_position = mappings_per_chromosome[chrom_pos]['uniprot_position']
    
        if not residue_position in variant_type_counts.keys():
            variant_type_counts[residue_position] = dict()
            
            # Add the background mutation type rates for this codon        
            variant_type_counts[residue_position]['background_missense'] = codon_background_rates[codon]['missense']
            variant_type_counts[residue_position]['background_nonsense'] = codon_background_rates[codon]['nonsense']
            variant_type_counts[residue_position]['background_synonymous'] = codon_background_rates[codon]['synonymous']
            
            # Add the background mutation type rates for this codon        
            variant_type_counts[residue_position]['missense'] = 0
            variant_type_counts[residue_position]['nonsense'] = 0
            variant_type_counts[residue_position]['synonymous'] = 0
    
    # return the variant type counts
    return variant_type_counts
            
def retrieve_variant_type_counts(mappings_per_chromosome, annotated_region):
    """Given an annotated region, returns a per position count"""
    # First retrieve the background for the region
    variant_type_counts = retrieve_background_variant_counts(mappings_per_chromosome)
    
    for chrom_pos in annotated_region.keys():
        codon = mappings_per_chromosome[chrom_pos]['base_pair_representation']
        codon_pos = mappings_per_chromosome[chrom_pos]['codon_base_pair_position']
        residue_position = mappings_per_chromosome[chrom_pos]['uniprot_position']
        
        for annotation_entry in annotated_region[chrom_pos]:
            if mappings_per_chromosome[chrom_pos]['base_pair'] != annotation_entry['REF']:
                raise ExternalREFAlleleNotEqualsTranscriptionException("For gene id '"+str(mappings_per_chromosome[chrom_pos]['gencode_transcription_id'])+
                                                   "', at chrom_pos '"+str(chrom_pos)+"' : analysis of protein region could"+
                                                   " not be made due to: gene_region.retrieve_mappings_per_chromosome()[chrom_pos].base_pair"+
                                                   " != annotation_entry['REF']")
            
            # interpret the variant type
            alt_codon = SingleNucleotideVariant.interpret_alt_codon(ref_basepair_representation=codon, var_codon_pos=codon_pos, alt_nucleotide=annotation_entry['ALT'])
            variant_type = SingleNucleotideVariant.interpret_variant_type_from_codon_basepair_representations(ref_basepair_representation=codon, alt_basepair_representation=alt_codon).value
            
            # add variant_type to variant_type counts
            variant_type_counts[residue_position][variant_type] += 1
                
    return variant_type_counts