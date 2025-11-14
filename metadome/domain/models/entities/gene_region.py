from metadome.domain.repositories import MappingRepository, ProteinRepository, InterproRepository
from metadome.domain.services.helper_functions import convertListOfIntegerToRanges
from metadome.domain.models.entities.codon import Codon
import enum

class GenomeBuild(enum.Enum):
    GRCh37 = 'GRCh37'
    GRCh38 = 'GRCh38'

class FailedToConstructGeneRegion(Exception):
    pass

class RegioncDNALengthDoesNotEqualProteinLengthException(Exception):
    pass

class MalformedGeneRegionException(Exception):
    pass

class GeneRegion(object):
    """
    GeneRegion Model Entity
    Used for representation of (partial) gene regions for annotation purposes
    
    Variables
    name                      description
    chr                       str containing the chromosome number in form 'chr...'
    gene_name                 str containing the gene name of the region
    gencode_transcription_id  str containing the transcription id
    gencode_version           str containing the gencode version
    gencode_basic             bool indicating if this is a gencode basic transcript
    genome_build_full         str containing the genome build
    genome_build              GenomeBuild enum containing the genome build
    refseq_transcript_id      str containing the refseq transcript id
    mane_transcript_type      str containing the mane transcript type
    strand                    str containing the strand
    protein_region_start      int of start of the region 
    protein_region_stop       int of stop of the region
    protein_region_length     int representing the region length based on the amino acid sequence
    cDNA_region_length        int representing the region length based on the cDNA sequence
    protein_id                int identifier corresponding to the protein id in the database
    uniprot_ac                str representing the region's uniprot accession code
    uniprot_name              str representing the region's uniprot name
    interpro_domains          list containing interpro domains (models.interpro.Interpro)
    regions                   list of regions
    mappings_per_chromosome   dictionary of mappings per chromosome positions; {POS: models.mapping.Mapping}
    """
    
    def retrieve_specific_domains_in_gene(self, domain_id):
        """retrieves specific domains in a gene based on the 
        provided domain id, sorted by the start position"""
        return sorted([domain for domain in self.interpro_domains if domain.ext_db_id == domain_id], key=lambda k: k.uniprot_start)
    
    def retrieve_codon_for_protein_position(self, protein_pos):
        """Returns the codon for this gene region at the protein position"""
        _codon = None
        # retrieve te mappings
        _mappings = []
        for cDNA_pos in self.protein_pos_to_cDNA[protein_pos]:
            _mappings.append(self.mappings_per_cDNA[cDNA_pos])
        # create the codon
        _codon = Codon.initializeFromMapping(_mappings, self.gencode_transcription_id, self.uniprot_ac)
        return _codon
    
    def retrieve_mappings_per_chromosome(self):
        """Returns the mappings for this gene region per chromosome position"""
        mappings_per_chromosome = dict()
        for chromosome_position in self.chromosome_pos_to_cDNA.keys():            
            mappings_per_chromosome[chromosome_position] = {}
            mappings_per_chromosome[chromosome_position]['base_pair_representation'] = self.mappings_per_cDNA[self.chromosome_pos_to_cDNA[chromosome_position]].codon
            mappings_per_chromosome[chromosome_position]['codon_base_pair_position'] = self.mappings_per_cDNA[self.chromosome_pos_to_cDNA[chromosome_position]].codon_base_pair_position
            mappings_per_chromosome[chromosome_position]['uniprot_position'] = self.mappings_per_cDNA[self.chromosome_pos_to_cDNA[chromosome_position]].uniprot_position
            mappings_per_chromosome[chromosome_position]['base_pair'] = self.mappings_per_cDNA[self.chromosome_pos_to_cDNA[chromosome_position]].base_pair
            mappings_per_chromosome[chromosome_position]['gencode_transcription_id'] = self.gencode_transcription_id
        return mappings_per_chromosome
    
    def get_domains_for_position(self, position):
        """Checks if there are any protein domains present at the given position"""
        return [domain for domain in self.interpro_domains if domain.uniprot_start-1 <= position <= domain.uniprot_stop]
    
    
    def __init__(self, _gene):
        self.chr = str()
        self.gene_name = str()
        self.gene_id = str()
        self.gencode_transcription_id = str()
        self.gencode_version = str()
        self.gencode_basic = bool()
        self.genome_build_full = str()
        self.genome_build = None
        self.refseq_transcript_id = str()
        self.mane_transcript_type = str()
        self.strand = str()
        self.protein_region_start = int()
        self.protein_region_stop = int()
        self.protein_region_length = int()
        self.cDNA_region_length = int()
        self.protein_id = int()
        self.uniprot_ac = str()
        self.uniprot_name = str()
        self.interpro_domains = []
        self.regions = []
        self.mappings_per_cDNA = dict()
        self.chromosome_pos_to_cDNA = dict()
        self.protein_pos_to_cDNA = dict()
        self.initialize_from_gene(_gene)
    
    def initialize_from_gene(self, _gene):
        # First set the fields that can directly be set from the _gene
        self.gene_name = _gene.gene_name
        self.gene_id = _gene.id
        self.protein_region_start = 0
        self.protein_region_stop = _gene.sequence_length            
        self.gencode_transcription_id = _gene.gencode_transcription_id
        self.gencode_version = _gene.gencode_version
        self.gencode_basic = _gene.gencode_basic
        self.genome_build_full = _gene.genome_build
        self.refseq_transcript_id = _gene.refseq_transcript_id if not _gene.refseq_transcript_id is None else str()
        self.mane_transcript_type = _gene.mane_transcript_type if not _gene.mane_transcript_type is None else str()
        self.protein_region_length = self.protein_region_stop - self.protein_region_start
        self.strand = _gene.strand

        # set the genome build enum
        if self.genome_build_full.upper()[:6] == 'GRCH37':
            self.genome_build = GenomeBuild.GRCh37
        elif self.genome_build_full.upper()[:6] == 'GRCH38':
            self.genome_build = GenomeBuild.GRCh38
        else:
            raise Exception("No valid genome build found for gene '"+str(self.gene_name)+"', found '"+str(self.genome_build_full)+"'")
        
        # retrieve the protein id
        _protein = ProteinRepository.retrieve_protein(_gene.protein_id)
        
        # set the uniprot ac and name
        self.protein_id = _protein.id
        self.uniprot_ac = _protein.uniprot_ac
        self.uniprot_name = _protein.uniprot_name
        self.interpro_domains = InterproRepository.get_domains_for_protein(self.protein_id)
        
        # Retrieve mappings from the database
        self.mappings_per_cDNA = {x.cDNA_position:x for x in MappingRepository.get_mappings_for_gene(_gene)}

        # Check f everything went fine so far
        if len(self.mappings_per_cDNA) == 0:
            raise FailedToConstructGeneRegion("For transcript '"+str(self.gencode_transcription_id)+
                                              "': No proper mapping was found in the database,"+
                                              " failed to initialize GeneRegion object")
        
        # iterate over all cDNA positions
        for cDNA_position in sorted( self.mappings_per_cDNA.keys(), key=lambda x: int(x) ):
            _mapping = self.mappings_per_cDNA[cDNA_position]
            
            # type check protein id
            if self.protein_id !=  _mapping.protein_id:
                raise MalformedGeneRegionException("For transcript '"+str(self.gencode_transcription_id)+
                                                   "': Found mappings for a single transcript aligned to multiple different proteins")


            # test if we already set the chromosome
            if self.chr == str():
                self.chr = _mapping.chromosome
            else:
                # type check chromosome
                if self.chr !=  _mapping.chromosome:
                    raise MalformedGeneRegionException("For transcript '"+str(self.gencode_transcription_id)+
                                                       "': Found alignments to multiple chromosomes")
            
            # ensure we have alignment here
            if _mapping.uniprot_residue == '-' or _mapping.uniprot_residue == '*':
                continue
            
            # Add the other mapping of other positions to this cDDA
            self.chromosome_pos_to_cDNA[_mapping.chromosome_position] = int(cDNA_position)
            if not _mapping.uniprot_position in self.protein_pos_to_cDNA.keys():
                self.protein_pos_to_cDNA[_mapping.uniprot_position] = []
            self.protein_pos_to_cDNA[_mapping.uniprot_position].append(int(_mapping.cDNA_position))
        
        # ensure that the region is fully covered
        if(self.protein_region_length != len(self.protein_pos_to_cDNA.keys())):
            raise RegioncDNALengthDoesNotEqualProteinLengthException("For transcript '"+str(self.gencode_transcription_id)+
                                                                     "': Attempted to recover chromosomal region, but the region was not"+
                                                                     " perfectly aligned, expected AA sequence length to be '"+str(self.protein_region_length)+"'"+
                                                                     ", but was '"+str(len(self.protein_pos_to_cDNA.keys()))+"'")
        
        # now add the cDNA region length
        self.cDNA_region_length = len(self.mappings_per_cDNA.keys())

        # convert the chromosome to ranges
        self.regions = list(convertListOfIntegerToRanges(sorted(self.chromosome_pos_to_cDNA.keys())))

    def __repr__(self):
        return "<GeneRegion(chr='%s', gene_name='%s', gencode_transcription_id='%s', protein_region_length='%s', cDNA_region_length='%s', strand='%s', number of chromosomal regions='%s')>" % (
                            self.chr, self.gene_name, self.gencode_transcription_id, self.protein_region_length, self.cDNA_region_length, self.strand, len(self.regions))