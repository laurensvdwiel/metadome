from metadome.domain.models.entities.gene_region import GenomeBuild
from metadome.domain.services.annotation.annotation import annotateSNVs
from metadome.domain.models.entities.single_nucleotide_variant import SingleNucleotideVariant
from metadome.domain.services.annotation.gene_region_annotators import GRCh37_annotateTranscriptWithGnomADData, \
    GRCh37_annotateTranscriptWithClinvarData, GRCh38_annotateTranscriptWithClinvarData, GRCh38_annotateTranscriptWithGnomADData

def annotate_ClinVar_SNVs_for_codons(codons, genome_build):
    """Annotate provided codons with ClinVar SNVs"""
    if genome_build == GenomeBuild.GRCh37:
        clinvar_annotation_function = GRCh37_annotateTranscriptWithClinvarData
    elif genome_build == GenomeBuild.GRCh38:
        clinvar_annotation_function = GRCh38_annotateTranscriptWithClinvarData
    else:
        raise Exception("Could not determine ClinVar annotation function based on genome build '" + str(
            genome_build.value) + "'")
    return annotate_SNVs_for_codons(annotateTranscriptFunction=clinvar_annotation_function, codons=codons, variant_source='ClinVar')

def annotate_gnomAD_SNVs_for_codons(codons, genome_build):
    """Annotate provided codons with gnomAD SNVs"""
    if genome_build == GenomeBuild.GRCh37:
        gnomad_annotation_function = GRCh37_annotateTranscriptWithGnomADData
    elif genome_build == GenomeBuild.GRCh38:
        gnomad_annotation_function = GRCh38_annotateTranscriptWithGnomADData
    else:
        raise Exception("Could not determine gnomAD annotation function based on genome build '" + str(
            genome_build.value) + "'")
    return annotate_SNVs_for_codons(annotateTranscriptFunction=gnomad_annotation_function, codons=codons, variant_source='gnomAD')

def annotate_SNVs_for_codons(annotateTranscriptFunction, codons, variant_source):
    """Annotate provided codons with SNVs from the provided variant source and transcript annotation function"""
    # the list that will be returned
    SNVs = []
    
    # aggregate the identical codons
    meta_codons = {}
    for codon in codons:
        if not codon.unique_str_representation() in meta_codons.keys():
            meta_codons[codon.unique_str_representation()] = []
        
        meta_codons[codon.unique_str_representation()].append(codon)
        
    # iterate over meta_codons and add to metadom_entry
    for meta_codon_repr in meta_codons.keys():
        # Check if we are dealing with the gene and protein_pos of interest
        # just take the first
        meta_codon = meta_codons[meta_codon_repr][0]
    
        # annotate the variant
        variant_annotation = annotateSNVs(annotateTranscriptFunction,
                                                 mappings_per_chr_pos=meta_codon.retrieve_mappings_per_chromosome(),
                                                 strand=meta_codon.strand, 
                                                 chromosome=meta_codon.chr,
                                                 regions=meta_codon.regions)
    
        for chrom_pos in variant_annotation.keys():
            for variant in variant_annotation[chrom_pos]:
                for _codon in meta_codons[meta_codon_repr]:
                    # create a variant object for this variant
                    SNV = SingleNucleotideVariant.initializeFromVariant(_codon=_codon, _chr_position=chrom_pos, _alt_nucleotide=variant['ALT'], _variant_source=variant_source)
                    
                    # create the variant entry
                    variant_entry = extract_variant_source_info(SNV.toDict(), variant, variant_source)
                    variant_entry['unique_snv_str_representation'] = SNV.unique_snv_str_representation()
                                        
                    SNVs.append(variant_entry)
    
    return SNVs

def extract_variant_source_info(variant_entry, variant_annotation, variant_source):
    if variant_source == 'gnomAD':
        return extract_gnomAD_variant_info(variant_entry, variant_annotation)
    elif variant_source == 'ClinVar':
        return extract_ClinVar_variant_info(variant_entry, variant_annotation)
    else:
        raise NotImplementedError("Variant source '"+str(variant_source)+"' is not supported for additional info")

def extract_gnomAD_variant_info(variant_entry, variant_annotation):
    # append gnomAD specific information
    variant_entry['allele_number'] = variant_annotation['AN']
    variant_entry['allele_count'] = variant_annotation['AC']
    
    return variant_entry

def extract_ClinVar_variant_info(variant_entry, variant_annotation):    
    # append ClinVar specific information
    variant_entry['clinvar_ID'] = variant_annotation['ID']
    
    return variant_entry
