from metadome.default_settings import GNOMAD_GRCH37_VCF_FILE, GNOMAD_GRCH38_VCF_FILE, GNOMAD_ACCEPTED_FILTERS,\
 CLINVAR_CONSIDERED_CLINSIG, CLINVAR_GRCH37_VCF_FILE, CLINVAR_GRCH38_VCF_FILE
from metadome.domain.parsers.tabix import tabix_query, variant_coordinate_system

def GRCh37_annotateTranscriptWithClinvarData(chromosome, regions):
    """
    Annotates variants found within the ClinVar dataset.
    Specific clinical signifcance level of interest is
    filtered (See INFO field CLNSIG)
    
    Additional fields:
    ID            ClinVar Variation ID
    
    ClinVar specific INFO fields:
    AF_ESP              allele frequencies from GO-ESP
    AF_EXAC             allele frequencies from ExAC
    AF_TGP              allele frequencies from TGP
    ALLELEID            The ClinVar Allele ID
    CLNDN               ClinVar's preferred disease name for the concept specified by disease identifiers in CLNDISDB
    CLNDNINCL           For included Variant : ClinVar's preferred disease name for the concept specified by disease identifiers in CLNDISDB
    CLNDISDB            Tag-value pairs of disease database name and identifier, e.g. OMIM:NNNNNN
    CLNDISDBINCL        For included Variant: Tag-value pairs of disease database name and identifier, e.g. OMIM:NNNNNN
    CLNHGVS             Top-level (primary assembly, alt, or patch) HGVS expression.
    CLNREVSTAT          ClinVar review status for the Variation ID
    CLNSIG              Clinical significance for this single variant
    CLNSIGINCL          Clinical significance for a haplotype or genotype that includes this variant. Reported as pairs of VariationID:clinical significance.
    CLNVC               Variant type
    CLNVCSO             Sequence Ontology id for variant type
    CLNVI               The variant's clinical sources reported as tag-value pairs of database and variant identifier
    DBVARID             nsv accessions from dbVar for the variant
    GENEINFO            Gene(s) for the variant reported as gene symbol:gene id. The gene symbol and id are delimited by a colon (:) and each pair is delimited by a vertical bar (|)
    MC                  comma separated list of molecular consequence in the form of Sequence Ontology ID|molecular_consequence
    ORIGIN              Allele origin. One or more of the following values may be added: 0 - unknown; 1 - germline; 2 - somatic; 4 - inherited; 8 - paternal; 16 - maternal; 32 - de-novo; 64 - biparental; 128 - uniparental; 256 - not-tested; 512 - tested-inconclusive; 1073741824 - other
    RS                  dbSNP ID (i.e. rs number)
    SSR                 Variant Suspect Reason Codes. One or more of the following values may be added: 0 - unspecified, 1 - Paralog, 2 - byEST, 4 - oldAlign, 8 - Para_EST, 16 - 1kg_failed, 1024 - other
    """
    for gene_sub_region in regions:
        for tabix_record in tabix_query(CLINVAR_GRCH37_VCF_FILE, chromosome[3:], gene_sub_region[0], gene_sub_region[1], variant_coordinate_system.one_based):
            for _, item in enumerate(tabix_record.ALT):
                clinvar_info = {"AF_ESP": None if "AF_ESP" not in tabix_record.INFO.keys() else tabix_record.INFO["AF_ESP"],
                            "AF_EXAC": None if "AF_EXAC" not in tabix_record.INFO.keys() else tabix_record.INFO["AF_EXAC"],
                            "AF_TGP": None if "AF_TGP" not in tabix_record.INFO.keys() else tabix_record.INFO["AF_TGP"],
                            "CLNDNINCL": None if "CLNDNINCL" not in tabix_record.INFO.keys() else tabix_record.INFO["CLNDNINCL"],
                            "CLNDISDBINCL": None if "CLNDISDBINCL" not in tabix_record.INFO.keys() else tabix_record.INFO["CLNDISDBINCL"],
                            "CLNSIGINCL": None if "CLNSIGINCL" not in tabix_record.INFO.keys() else tabix_record.INFO["CLNSIGINCL"],
                            "DBVARID": None if "DBVARID" not in tabix_record.INFO.keys() else tabix_record.INFO["DBVARID"],
                            "SSR": None if "SSR" not in tabix_record.INFO.keys() else tabix_record.INFO["SSR"],
                            "CLNVI": None if "CLNVI" not in tabix_record.INFO.keys() else tabix_record.INFO["CLNVI"],
                            "MC": None if "MC" not in tabix_record.INFO.keys() else tabix_record.INFO["MC"],
                            "RS": None if "RS" not in tabix_record.INFO.keys() else tabix_record.INFO["RS"],
                            "ORIGIN": None if "ORIGIN" not in tabix_record.INFO.keys() else tabix_record.INFO["ORIGIN"],
                            "CLNDN": None if "CLNDN" not in tabix_record.INFO.keys() else tabix_record.INFO["CLNDN"],
                            "CLNDISDB": None if "CLNDISDB" not in tabix_record.INFO.keys() else tabix_record.INFO["CLNDISDB"], 
                            "CLNSIG": None if "CLNSIG" not in tabix_record.INFO.keys() else tabix_record.INFO["CLNSIG"],
                            "GENEINFO": None if "GENEINFO" not in tabix_record.INFO.keys() else tabix_record.INFO["GENEINFO"],
                            "CLNREVSTAT": None if "CLNREVSTAT" not in tabix_record.INFO.keys() else tabix_record.INFO["CLNREVSTAT"],
                            "ALLELEID": None if "ALLELEID" not in tabix_record.INFO.keys() else tabix_record.INFO["ALLELEID"],
                            "CLNHGVS": None if "CLNHGVS" not in tabix_record.INFO.keys() else tabix_record.INFO["CLNHGVS"],
                            "CLNVC": None if "CLNVC" not in tabix_record.INFO.keys() else tabix_record.INFO["CLNVC"],
                            "CLNVCSO": None if "CLNVCSO" not in tabix_record.INFO.keys() else tabix_record.INFO["CLNVCSO"],
                            "ID":tabix_record.ID}
            
                clinvar_record = {'CHROM': tabix_record.CHROM, 'POS': tabix_record.POS, 'REF':tabix_record.REF, 'ALT':item, 'INFO':clinvar_info}
                
                if not clinvar_record['INFO']['CLNSIG'] is None and\
                        len( clinvar_record['INFO']['CLNSIG'] ) == 1 and\
                        clinvar_record['INFO']['CLNSIG'][0] in CLINVAR_CONSIDERED_CLINSIG and \
                        clinvar_record['ALT'] != 'N' and \
                        clinvar_record['REF'] != 'N':
                    yield clinvar_record


def GRCh38_annotateTranscriptWithClinvarData(chromosome, regions):
    """
    Annotates variants found within the ClinVar dataset.
    Specific clinical signifcance level of interest is
    filtered (See INFO field CLNSIG)

    Additional fields:
    ID            ClinVar Variation ID

    ClinVar specific INFO fields:
    AF_ESP              allele frequencies from GO-ESP
    AF_EXAC             allele frequencies from ExAC
    AF_TGP              allele frequencies from TGP
    ALLELEID            The ClinVar Allele ID
    CLNDN               ClinVar's preferred disease name for the concept specified by disease identifiers in CLNDISDB
    CLNDNINCL           For included Variant : ClinVar's preferred disease name for the concept specified by disease identifiers in CLNDISDB
    CLNDISDB            Tag-value pairs of disease database name and identifier, e.g. OMIM:NNNNNN
    CLNDISDBINCL        For included Variant: Tag-value pairs of disease database name and identifier, e.g. OMIM:NNNNNN
    CLNHGVS             Top-level (primary assembly, alt, or patch) HGVS expression.
    CLNREVSTAT          ClinVar review status for the Variation ID
    CLNSIG              Clinical significance for this single variant
    CLNSIGINCL          Clinical significance for a haplotype or genotype that includes this variant. Reported as pairs of VariationID:clinical significance.
    CLNVC               Variant type
    CLNVCSO             Sequence Ontology id for variant type
    CLNVI               The variant's clinical sources reported as tag-value pairs of database and variant identifier
    DBVARID             nsv accessions from dbVar for the variant
    GENEINFO            Gene(s) for the variant reported as gene symbol:gene id. The gene symbol and id are delimited by a colon (:) and each pair is delimited by a vertical bar (|)
    MC                  comma separated list of molecular consequence in the form of Sequence Ontology ID|molecular_consequence
    ORIGIN              Allele origin. One or more of the following values may be added: 0 - unknown; 1 - germline; 2 - somatic; 4 - inherited; 8 - paternal; 16 - maternal; 32 - de-novo; 64 - biparental; 128 - uniparental; 256 - not-tested; 512 - tested-inconclusive; 1073741824 - other
    RS                  dbSNP ID (i.e. rs number)
    SSR                 Variant Suspect Reason Codes. One or more of the following values may be added: 0 - unspecified, 1 - Paralog, 2 - byEST, 4 - oldAlign, 8 - Para_EST, 16 - 1kg_failed, 1024 - other
    """
    for gene_sub_region in regions:
        for tabix_record in tabix_query(CLINVAR_GRCH38_VCF_FILE, chromosome[3:], gene_sub_region[0], gene_sub_region[1],
                                        variant_coordinate_system.one_based):
            for _, item in enumerate(tabix_record.ALT):
                clinvar_info = {"AF_ESP": None if "AF_ESP" not in tabix_record.INFO.keys() else tabix_record.INFO["AF_ESP"],
                            "AF_EXAC": None if "AF_EXAC" not in tabix_record.INFO.keys() else tabix_record.INFO["AF_EXAC"],
                            "AF_TGP": None if "AF_TGP" not in tabix_record.INFO.keys() else tabix_record.INFO["AF_TGP"],
                            "CLNDNINCL": None if "CLNDNINCL" not in tabix_record.INFO.keys() else tabix_record.INFO["CLNDNINCL"],
                            "CLNDISDBINCL": None if "CLNDISDBINCL" not in tabix_record.INFO.keys() else tabix_record.INFO["CLNDISDBINCL"],
                            "CLNSIGINCL": None if "CLNSIGINCL" not in tabix_record.INFO.keys() else tabix_record.INFO["CLNSIGINCL"],
                            "DBVARID": None if "DBVARID" not in tabix_record.INFO.keys() else tabix_record.INFO["DBVARID"],
                            "SSR": None if "SSR" not in tabix_record.INFO.keys() else tabix_record.INFO["SSR"],
                            "CLNVI": None if "CLNVI" not in tabix_record.INFO.keys() else tabix_record.INFO["CLNVI"],
                            "MC": None if "MC" not in tabix_record.INFO.keys() else tabix_record.INFO["MC"],
                            "RS": None if "RS" not in tabix_record.INFO.keys() else tabix_record.INFO["RS"],
                            "ORIGIN": None if "ORIGIN" not in tabix_record.INFO.keys() else tabix_record.INFO["ORIGIN"],
                            "CLNDN": None if "CLNDN" not in tabix_record.INFO.keys() else tabix_record.INFO["CLNDN"],
                            "CLNDISDB": None if "CLNDISDB" not in tabix_record.INFO.keys() else tabix_record.INFO["CLNDISDB"],
                            "CLNSIG": None if "CLNSIG" not in tabix_record.INFO.keys() else tabix_record.INFO["CLNSIG"],
                            "GENEINFO": None if "GENEINFO" not in tabix_record.INFO.keys() else tabix_record.INFO["GENEINFO"],
                            "CLNREVSTAT": None if "CLNREVSTAT" not in tabix_record.INFO.keys() else tabix_record.INFO["CLNREVSTAT"],
                            "ALLELEID": None if "ALLELEID" not in tabix_record.INFO.keys() else tabix_record.INFO["ALLELEID"],
                            "CLNHGVS": None if "CLNHGVS" not in tabix_record.INFO.keys() else tabix_record.INFO["CLNHGVS"],
                            "CLNVC": None if "CLNVC" not in tabix_record.INFO.keys() else tabix_record.INFO["CLNVC"],
                            "CLNVCSO": None if "CLNVCSO" not in tabix_record.INFO.keys() else tabix_record.INFO["CLNVCSO"],
                            "ID":tabix_record.ID}

                clinvar_record = {'CHROM': tabix_record.CHROM, 'POS': tabix_record.POS, 'REF': tabix_record.REF,
                                  'ALT': item, 'INFO': clinvar_info}

                if not clinvar_record['INFO']['CLNSIG'] is None and \
                        len(clinvar_record['INFO']['CLNSIG']) == 1 and \
                        clinvar_record['INFO']['CLNSIG'][0] in CLINVAR_CONSIDERED_CLINSIG and \
                        clinvar_record['ALT'] != 'N' and \
                        clinvar_record['REF'] != 'N':
                    yield clinvar_record

def GRCh37_annotateTranscriptWithGnomADData(chromosome, regions):
    """
    Annotates variants found within the gnomAD dataset for GRCh37 with specific FILTER settings.
        PASS : passed all variant filters imposed, all such variants are considered real variants.
    """
    for gene_sub_region in regions:
        for tabix_record in tabix_query(GNOMAD_GRCH37_VCF_FILE, chromosome[3:], gene_sub_region[0], gene_sub_region[1], variant_coordinate_system.one_based):
            for i, item in enumerate(tabix_record.ALT):
                gnomad_filter = ''
                if len(tabix_record.FILTER) == 0:
                    gnomad_filter = 'PASS'
                else:
                    gnomad_filter = tabix_record.FILTER[0]
                
                if gnomad_filter in GNOMAD_ACCEPTED_FILTERS:
                    gnomad_record = {'CHROM': tabix_record.CHROM, 'POS': tabix_record.POS, 'FILTER':gnomad_filter, 'REF':tabix_record.REF, 'ALT':item, 'INFO':{'AC':tabix_record.INFO['AC'][i], 'AF':tabix_record.INFO['AF'][i], 'AN':tabix_record.INFO['AN']}}
                    yield gnomad_record

def GRCh38_annotateTranscriptWithGnomADData(chromosome, regions):
    """
    Annotates variants found within the gnomAD dataset for GRCh38 with specific FILTER settings.
        PASS : passed all variant filters imposed, all such variants are considered real variants.
    """
    for gene_sub_region in regions:
        for tabix_record in tabix_query(GNOMAD_GRCH38_VCF_FILE, chromosome, gene_sub_region[0],
                                        gene_sub_region[1], variant_coordinate_system.one_based):
            for i, item in enumerate(tabix_record.ALT):
                gnomad_filter = ''
                if len(tabix_record.FILTER) == 0:
                    gnomad_filter = 'PASS'
                else:
                    gnomad_filter = tabix_record.FILTER[0]

                if gnomad_filter in GNOMAD_ACCEPTED_FILTERS:
                    gnomad_record = {'CHROM': tabix_record.CHROM[3:], 'POS': tabix_record.POS, 'FILTER': gnomad_filter,
                                     'REF': tabix_record.REF, 'ALT': item,
                                     'INFO': {'AC': tabix_record.INFO['AC_exomes'][i], 'AF': tabix_record.INFO['AF_exomes'][i],
                                              'AN': tabix_record.INFO['AN_exomes']}}
                    yield gnomad_record