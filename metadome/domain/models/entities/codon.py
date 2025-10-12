from metadome.domain.services.helper_functions import convertListOfIntegerToRanges, list_of_stringified_of_ranges
from metadome.domain.models.gene import Strand
from Bio.Data.IUPACData import protein_letters_1to3

class MalformedCodonException(Exception):
    pass

class Codon(object):
    """
    Codon Model Entity
    Used for a codon-representation of models.mapping.Mappings that make up a codon
    
    Variables
    name                                   description
    gencode_transcription_id               str the transcription id associated to this codon
    uniprot_ac                             str the uniprot protein ac associated to this codon
    strand                                 Enum the strand represented as models.gene.Strand
    base_pair_representation               str the base pair representation of this codon (e.g. ATG)
    amino_acid_residue                     str the amino acid residue of this codon
    amino_acid_position                    int the amino acid position of this codon in the protein or gene translation
    chr                                    str the chromosome
    regions                                tuple the range of chromosomal positions of this codon
    chromosome_position_base_pair_one      int the position corresponding to the first base pair in the chromosome
    chromosome_position_base_pair_two      int the position corresponding to the second base pair in the chromosome
    chromosome_position_base_pair_three    int the position corresponding to the third base pair in the chromosome
    cDNA_position_one                      int the position corresponding to the first base pair in the cDNA
    cDNA_position_two                      int the position corresponding to the second base pair in the cDNA
    cDNA_position_three                    int the position corresponding to the third base pair in the cDNA
    """
    @staticmethod
    def one_to_three_letter_amino_acid_residue(amino_acid_residue):
        """Returns a three letter representation of the provided amino acid residue"""
        # Check if this is a Pyrrolysine
        if amino_acid_residue == 'O':
            return 'Pyl'
        # Check if this is a Selenocysteine
        if amino_acid_residue == 'U':
            return 'Sec'
        # Return one of the 20 amino acid residues
        return protein_letters_1to3[amino_acid_residue];

    def unique_str_representation(self):
        return str(self.chr)+":"+str(self.regions)+"::("+str(self.strand)+")"
    
    def retrieve_mappings_per_chromosome(self):
        """Returns the mappings for this codon per chromosome position"""
        mappings_per_chromosome = dict()
        
        # first base pair
        mappings_per_chromosome[self.chromosome_position_base_pair_one] = {}
        mappings_per_chromosome[self.chromosome_position_base_pair_one]['base_pair_representation'] = self.base_pair_representation
        mappings_per_chromosome[self.chromosome_position_base_pair_one]['codon_base_pair_position'] = 0
        mappings_per_chromosome[self.chromosome_position_base_pair_one]['amino_acid_position'] = self.amino_acid_position
        mappings_per_chromosome[self.chromosome_position_base_pair_one]['base_pair'] = self.base_pair_representation[0]
        mappings_per_chromosome[self.chromosome_position_base_pair_one]['gencode_transcription_id'] = self.gencode_transcription_id

        # second base pair
        mappings_per_chromosome[self.chromosome_position_base_pair_two] = {}
        mappings_per_chromosome[self.chromosome_position_base_pair_two]['base_pair_representation'] = self.base_pair_representation
        mappings_per_chromosome[self.chromosome_position_base_pair_two]['codon_base_pair_position'] = 1
        mappings_per_chromosome[self.chromosome_position_base_pair_two]['amino_acid_position'] = self.amino_acid_position
        mappings_per_chromosome[self.chromosome_position_base_pair_two]['base_pair'] = self.base_pair_representation[1]
        mappings_per_chromosome[self.chromosome_position_base_pair_two]['gencode_transcription_id'] = self.gencode_transcription_id

        # third base pair
        mappings_per_chromosome[self.chromosome_position_base_pair_three] = {}
        mappings_per_chromosome[self.chromosome_position_base_pair_three]['base_pair_representation'] = self.base_pair_representation
        mappings_per_chromosome[self.chromosome_position_base_pair_three]['codon_base_pair_position'] = 2
        mappings_per_chromosome[self.chromosome_position_base_pair_three]['amino_acid_position'] = self.amino_acid_position
        mappings_per_chromosome[self.chromosome_position_base_pair_three]['base_pair'] = self.base_pair_representation[2]
        mappings_per_chromosome[self.chromosome_position_base_pair_three]['gencode_transcription_id'] = self.gencode_transcription_id

        return mappings_per_chromosome
    
    def three_letter_amino_acid_residue(self):
        """Returns a three letter representation of the amino acid residue for this codon"""
        return Codon.one_to_three_letter_amino_acid_residue(self.amino_acid_residue)
    
    def pretty_print_cDNA_region(self):
        return "c."+str(self.cDNA_position_one)+"-"+str(self.cDNA_position_three)
    
    def pretty_print_chr_region(self):
        _stringified_list = list_of_stringified_of_ranges(self.regions)
        return "".join("g."+_stringified_list[i] if i+1 == len(_stringified_list) else "g."+_stringified_list[i]+", " for i in range(len(_stringified_list)))

    def __init__(self, _gencode_transcription_id, _uniprot_ac, 
                             _strand, _base_pair_representation, 
                             _amino_acid_residue, _amino_acid_position, 
                             _chr, _chromosome_position_base_pair_one, 
                             _chromosome_position_base_pair_two, 
                             _chromosome_position_base_pair_three,
                             _cDNA_position_one, _cDNA_position_two,
                             _cDNA_position_three):
        
        self.gencode_transcription_id = _gencode_transcription_id
        self.uniprot_ac = _uniprot_ac
        self.strand = str()
        if _strand == '-':
            self.strand = Strand.minus
        elif _strand == '+':
            self.strand = Strand.plus
        else:
            raise MalformedCodonException('Strand cold not be converted to domain.models.Gene.Strand(Enum) for provided: '+str(_strand))
        self.base_pair_representation = _base_pair_representation
        self.amino_acid_residue = _amino_acid_residue
        self.amino_acid_position = _amino_acid_position
        self.chr = _chr
        self.chromosome_position_base_pair_one = _chromosome_position_base_pair_one
        self.chromosome_position_base_pair_two = _chromosome_position_base_pair_two
        self.chromosome_position_base_pair_three = _chromosome_position_base_pair_three
        self.cDNA_position_one = _cDNA_position_one
        self.cDNA_position_two = _cDNA_position_two
        self.cDNA_position_three = _cDNA_position_three
        
        # Set the regions the same way as we set the regions in a gene_region
        self.regions = list(convertListOfIntegerToRanges([self.chromosome_position_base_pair_one, self.chromosome_position_base_pair_two, self.chromosome_position_base_pair_three]))
    
    @classmethod
    def initializeFromMapping(cls, _mappings, _gencode_transcription_id, _uniprot_ac):
        """Initializes a Codon object via a list of three mapping 
        objects that represent the same codon"""
        # check the mappings cover exactly 3 mappings (thus represent a codon)
        if len(_mappings) != 3:
            raise MalformedCodonException("Malformed codon mapping: Expected exactly 3 mappings but got '"+str(len(_mappings))+"'.")
        
        # sort the mappings
        _mappings = sorted(_mappings, key=lambda k: k.codon_base_pair_position)
        _mapping_ids = [_mappings[0].id, _mappings[1].id, _mappings[2].id, ]
        
        # double check the mapping for gene ids
        if not(_mappings[0].gene_id == _mappings[1].gene_id == _mappings[2].gene_id):
            raise MalformedCodonException("Malformed codon mapping: The gene ids are not the same for mapping ids: '"+str(_mapping_ids)+"'.")
        
        # double check the mapping for gene ids
        if not(_mappings[0].protein_id == _mappings[1].protein_id == _mappings[2].protein_id):
            raise MalformedCodonException("Malformed codon mapping: The protein ids are not the same for mapping ids: '"+str(_mapping_ids)+"'.")
        
        # double check the mapping for amino_acid_residue
        if not(_mappings[0].amino_acid_residue == _mappings[1].amino_acid_residue == _mappings[2].amino_acid_residue):
            raise MalformedCodonException("Malformed codon mapping: The encoded amino acid residue does not correspond for mapping ids: '"+str(_mapping_ids)+"'.")
        
        # double check the mapping for amino_acid_position
        if not(_mappings[0].amino_acid_position == _mappings[1].amino_acid_position == _mappings[2].amino_acid_position):
            raise MalformedCodonException("Malformed codon mapping: The encoded amino acid position does not correspond for mapping ids: '"+str(_mapping_ids)+"'.")
        
        # double check the mapping for chromosome
        if not(_mappings[0].chromosome == _mappings[1].chromosome == _mappings[2].chromosome):
            raise MalformedCodonException("Malformed codon mapping: The chromosome does not correspond for mapping ids: '"+str(_mapping_ids)+"'.")
        
        # double check the mapping for strand
        if not(_mappings[0].strand == _mappings[1].strand == _mappings[2].strand):
            raise MalformedCodonException("Malformed codon mapping: The strand does not correspond for mapping ids: '"+str(_mapping_ids)+"'.")
        
        # double check the cDNA positions for the mappins
        if not(_mappings[0].cDNA_position < _mappings[1].cDNA_position < _mappings[2].cDNA_position):
            raise MalformedCodonException("Malformed codon mapping: The cDNA positions do not agree follow in order for mapping ids: '"+str(_mapping_ids)+"'.")
        
        # set the checked values
        _strand = _mappings[0].strand.value
        _amino_acid_residue = _mappings[0].amino_acid_residue
        _amino_acid_position = _mappings[0].amino_acid_position
        if _mappings[0].chromosome.lower().startswith('chr'):
            _chr = 'chr'+_mappings[0].chromosome[3:]
        else:
            _chr = 'chr'+_mappings[0].chromosome
        _cDNA_position_one = _mappings[0].cDNA_position
        _cDNA_position_two = _mappings[1].cDNA_position
        _cDNA_position_three = _mappings[2].cDNA_position
        
        # Add the chromosome positions for each basepair
        _chromosome_position_base_pair_one = _mappings[0].chromosome_position
        _chromosome_position_base_pair_two = _mappings[1].chromosome_position
        _chromosome_position_base_pair_three = _mappings[2].chromosome_position
        
        # check if the codon corresponds to the base pairs
        _base_pair_representation = _mappings[0].base_pair+_mappings[1].base_pair+_mappings[2].base_pair
        if not(_mappings[0].codon == _mappings[1].codon == _mappings[2].codon == _base_pair_representation):
            raise MalformedCodonException("Malformed codon mapping: The codons and base pairs of the mapping pair is malformed for mapping ids: '"+str(_mapping_ids)+"'.")
        
        _codon = cls(_gencode_transcription_id=_gencode_transcription_id,
                      _uniprot_ac=_uniprot_ac, _strand=_strand, 
                      _base_pair_representation=_base_pair_representation,
                      _amino_acid_residue=_amino_acid_residue, 
                      _amino_acid_position=_amino_acid_position, 
                      _chr=_chr, _chromosome_position_base_pair_one=_chromosome_position_base_pair_one, 
                      _chromosome_position_base_pair_two=_chromosome_position_base_pair_two, 
                      _chromosome_position_base_pair_three=_chromosome_position_base_pair_three,
                      _cDNA_position_one=_cDNA_position_one,
                      _cDNA_position_two=_cDNA_position_two,
                      _cDNA_position_three=_cDNA_position_three)
        
        return _codon
    
    @classmethod
    def initializeFromDict(cls, _d):
        """ Initializes a Codon object via a Dict
        Expected to contain at least the fields
        as given in the Codon.toDict() method
        """
        try:
            _gencode_transcription_id = _d['gencode_transcription_id']
            _uniprot_ac = _d['uniprot_ac']
            _strand = _d['strand']
            _base_pair_representation = _d['base_pair_representation']
            _amino_acid_residue = _d['amino_acid_residue']
            _amino_acid_position = _d['amino_acid_position']
            if _d['chr'].lower().startswith('chr'):
                _chr = 'chr' + _d['chr'][3:]
            else:
                _chr = 'chr'+_d['chr']
            _chromosome_position_base_pair_one = _d['chromosome_position_base_pair_one']
            _chromosome_position_base_pair_two = _d['chromosome_position_base_pair_two']
            _chromosome_position_base_pair_three = _d['chromosome_position_base_pair_three']
            _cDNA_position_one = _d['cDNA_position_one']
            _cDNA_position_two = _d['cDNA_position_two']
            _cDNA_position_three = _d['cDNA_position_three']
            
            _codon = cls(_gencode_transcription_id=_gencode_transcription_id,
                      _uniprot_ac=_uniprot_ac, _strand=_strand, 
                      _base_pair_representation=_base_pair_representation,
                      _amino_acid_residue=_amino_acid_residue, 
                      _amino_acid_position=_amino_acid_position, 
                      _chr=_chr, _chromosome_position_base_pair_one=_chromosome_position_base_pair_one, 
                      _chromosome_position_base_pair_two=_chromosome_position_base_pair_two, 
                      _chromosome_position_base_pair_three=_chromosome_position_base_pair_three,
                      _cDNA_position_one=_cDNA_position_one,
                      _cDNA_position_two=_cDNA_position_two,
                      _cDNA_position_three=_cDNA_position_three)
        
            return _codon
        except KeyError as e:
            raise MalformedCodonException("Malformed codon from dict: KeyError with message: "+str(e))
    
    def toDict(self):
        """Casts the basic values of a codon to a dictionary"""
        _d = {}
        _d['gencode_transcription_id'] = self.gencode_transcription_id
        _d['uniprot_ac'] = self.uniprot_ac
        _d['strand'] = self.strand.value
        _d['base_pair_representation'] = self.base_pair_representation
        _d['amino_acid_residue'] = self.amino_acid_residue
        _d['amino_acid_position'] = self.amino_acid_position
        _d['chr'] = self.chr
        _d['chromosome_position_base_pair_one'] = self.chromosome_position_base_pair_one
        _d['chromosome_position_base_pair_two'] = self.chromosome_position_base_pair_two
        _d['chromosome_position_base_pair_three'] = self.chromosome_position_base_pair_three
        _d['cDNA_position_one'] = self.cDNA_position_one
        _d['cDNA_position_two'] = self.cDNA_position_two
        _d['cDNA_position_three'] = self.cDNA_position_three
        
        return _d
    
    def toCodonJson(self):
        json_entry = {}
        
        # Add positional information
        json_entry['strand'] = self.strand.value
        json_entry['protein_pos'] = self.amino_acid_position
        json_entry['cdna_pos'] = self.pretty_print_cDNA_region()
        json_entry['chr'] = self.chr
        json_entry['chr_positions'] = self.pretty_print_chr_region()
        
        # Add residue and nucleotide information
        json_entry['ref_aa'] = self.amino_acid_residue
        json_entry['ref_aa_triplet'] = self.three_letter_amino_acid_residue()
        json_entry['ref_codon'] = self.base_pair_representation
        
        return json_entry
    
    def __repr__(self):
        return "<Codon(representation='%s', amino_acid_residue='%s', chr='%s', chr_positions='%s', strand='%s')>" % (
                            self.base_pair_representation, self.amino_acid_residue, self.chr, str(self.regions), self.strand )
