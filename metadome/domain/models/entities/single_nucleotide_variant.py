from metadome.domain.models.entities.codon import Codon, MalformedCodonException
from Bio.Seq import translate
import enum

class MalformedVariantException(Exception):
    pass

class VariantType(enum.Enum):
    missense = 'missense'
    nonsense = 'nonsense'
    synonymous = 'synonymous'

class SingleNucleotideVariant(Codon):
    """
    SingleNucleotideVariant Model Entity
    Used for a single nucleotide variant representation of models.entities.Codon
    
    Extends the Codon class with the following variables
    name                                   description
    var_codon_position                     int the position in the codon corresponding to the variant (0,1,2)
    ref_nucleotide                         str the reference nucleotide at the var_codon_position
    alt_nucleotide                         str the alternative nucleotide at the var_codon_position
    alt_amino_acid_residue                 str the alternative amino acid residue of this variant
    variant_type                           Enum the variant type (i.e.: missense, synonymous, nonsense)
    alt_base_pair_representation           str the alternative representation of this codon
    variant_source                         str the source of this variant
    """
    
    @staticmethod
    def interpret_alt_codon(ref_basepair_representation, var_codon_pos, alt_nucleotide):
        """Given the reference base pair representation of the codon, the variant 
        position in the codon and the alternative nucleotide, returns the 
        alternative base pair representation of the new codon"""
        alt_basepair_representation =""
        for i in range(len(ref_basepair_representation)):
            if i == var_codon_pos:
                alt_basepair_representation+= alt_nucleotide
            else:
                alt_basepair_representation+= ref_basepair_representation[i]
        
        return alt_basepair_representation
    
    @staticmethod
    def interpret_variant_type_from_codon_basepair_representations(ref_basepair_representation, alt_basepair_representation):
        """Interprets the models.entities.SingleNucleotideVariant.VariantType
        from translating the ref and alt basepair representations of the codon"""
        # translate the residues
        ref_residue = translate(ref_basepair_representation)
        alt_residue = translate(alt_basepair_representation)
        # interpret the variant type
        return SingleNucleotideVariant.interpret_variant_type_from_residues(ref_residue, alt_residue)
    
    @staticmethod
    def interpret_variant_type_from_residues(ref_residue, alt_residue):
        """Interprets the models.entities.SingleNucleotideVariant.VariantType
        from translating the ref and alt amino acid residues"""
        if alt_residue == '*':
            return VariantType.nonsense
        elif alt_residue != ref_residue:
            return VariantType.missense
        else:
            return VariantType.synonymous

    def get_genomic_position_of_change(self):
        """Retrieves the genomic variant position"""
        if self.var_codon_position == 0:
            return self.chromosome_position_base_pair_one
        elif self.var_codon_position == 1:
            return self.chromosome_position_base_pair_two
        elif self.var_codon_position == 2:
            return self.chromosome_position_base_pair_three
        else:
            raise MalformedVariantException("The var codon positions is not set correctly, expected in range 0-2m but was '"+str(self.var_codon_position)+"'")
        
    def alt_three_letter_amino_acid_residue(self):
        """Returns a three letter representation of the amino acid residue for this codon"""
        if self.alt_amino_acid_residue == '*':
            return self.alt_amino_acid_residue
        else:
            return Codon.one_to_three_letter_amino_acid_residue(self.alt_amino_acid_residue)
    
    def unique_snv_str_representation(self):
        return str(self.variant_source)+":"+self.unique_str_representation()+\
            ";"+str(self.base_pair_representation)+">"+str(self.alt_base_pair_representation)+\
            ";"+str(self.amino_acid_residue)+">"+str(self.alt_amino_acid_residue)

    def __init__(self, _gencode_transcription_id, _uniprot_ac, 
                             _strand, _base_pair_representation, 
                             _amino_acid_residue, _uniprot_position,
                             _chr, _exon_number_base_pair_one,
                             _exon_number_base_pair_two,
                             _exon_number_base_pair_three,
                             _chromosome_position_base_pair_one,
                             _chromosome_position_base_pair_two, 
                             _chromosome_position_base_pair_three,
                             _cDNA_position_one, _cDNA_position_two,
                             _cDNA_position_three, 
                             _variant_type,
                             _alt_amino_acid_residue,
                             _ref_nucleotide,
                             _alt_nucleotide,
                             _var_codon_position,
                             _variant_source):
        # Init the parent object
        super().__init__(_gencode_transcription_id=_gencode_transcription_id,
                         _uniprot_ac=_uniprot_ac, _strand=_strand, 
                         _base_pair_representation=_base_pair_representation, 
                         _amino_acid_residue=_amino_acid_residue, 
                         _uniprot_position=_uniprot_position,
                         _chr=_chr,
                         _exon_number_base_pair_one=_exon_number_base_pair_one,
                         _exon_number_base_pair_two=_exon_number_base_pair_two,
                         _exon_number_base_pair_three=_exon_number_base_pair_three,
                         _chromosome_position_base_pair_one=_chromosome_position_base_pair_one, 
                         _chromosome_position_base_pair_two=_chromosome_position_base_pair_two, 
                         _chromosome_position_base_pair_three=_chromosome_position_base_pair_three,
                         _cDNA_position_one=_cDNA_position_one, 
                         _cDNA_position_two=_cDNA_position_two,
                         _cDNA_position_three=_cDNA_position_three)
        # Init the rest for SingleNucleotideVariant
        self.ref_nucleotide = str()
        self.alt_nucleotide = str()
        self.var_codon_position = int()
        self.variant_type = str()
        self.alt_amino_acid_residue = _alt_amino_acid_residue
        self.alt_base_pair_representation = str()
        self.variant_source = _variant_source
        # start the type and rule checks
        if _ref_nucleotide == _alt_nucleotide:
            raise MalformedVariantException("No SNV could be made: found identical ref and alt nucleotides")
        if len(_ref_nucleotide) == len(_alt_nucleotide) == 1:
            self.ref_nucleotide = _ref_nucleotide
            self.alt_nucleotide = _alt_nucleotide
        else:
            raise MalformedVariantException("No SNV could be made: The ref '"+str(_ref_nucleotide)+"' and alt '"+str(_alt_nucleotide)+"' nucleotides should both have total length '1'")
        
        # double check that the codon position is indeed in range 0,1,2
        if 0 <= _var_codon_position <= 2:
            self.var_codon_position = _var_codon_position
        else:
            raise MalformedVariantException("No SNV could be made: _var_codon_position, indicating the variant in the basepair position in the codon should be 0, 1, or 2 but was '"+str(_var_codon_position)+"'")
        
        # check if the var position correspond to the ref nucleotide
        if self.base_pair_representation[self.var_codon_position] != self.ref_nucleotide:
            raise MalformedVariantException("No SNV could be made: the ref nucleoide '"+str(self.ref_nucleotide)+"' does not correspond to the var_codon_position '"+str(self.var_codon_position)+"' in the base_pair_representation '"+str(self.base_pair_representation)+"'")
        
        # set the alt base pair representation
        self.alt_base_pair_representation = SingleNucleotideVariant.interpret_alt_codon(self.base_pair_representation, self.var_codon_position, self.alt_nucleotide)        
         
        # check that the variant type is of missense, nonsense, synonymous
        if _variant_type == 'missense':
            self.variant_type = VariantType.missense
        elif _variant_type == 'synonymous':
            self.variant_type = VariantType.synonymous
        elif _variant_type == 'nonsense':
            self.variant_type = VariantType.nonsense
        else:
            raise MalformedVariantException('No SNV could be made: Variant type cold not be converted to domain.models.entities.SingleNucleotideVariant.VariantType(Enum) for provided: '+str(_variant_type))
    
    def toDict(self):
        # create the dictionary based on the codon
        _d = Codon.toDict(self)
        
        # Add the variables of SNV
        _d['ref_nucleotide'] = self.ref_nucleotide
        _d['alt_nucleotide'] = self.alt_nucleotide
        _d['var_codon_position'] = self.var_codon_position
        _d['variant_type'] = self.variant_type.value
        _d['alt_amino_acid_residue'] = self.alt_amino_acid_residue
        _d['variant_source'] = self.variant_source

        return _d
    
    @classmethod
    def initializeFromDict(cls, _d):
        try:
            # Double check the _codon is properly formatted
            _codon = Codon.initializeFromDict(_d)
            
            # initialize the other expected values from 
            _ref_nucleotide = _d['ref_nucleotide']
            _alt_nucleotide = _d['alt_nucleotide']
            _var_codon_position = _d['var_codon_position']
            _variant_type = _d['variant_type']
            _alt_amino_acid_residue = _d['alt_amino_acid_residue']
            _variant_source = _d['variant_source']
            
            SNV = cls(_gencode_transcription_id=_codon.gencode_transcription_id,
                                    _uniprot_ac=_codon.uniprot_ac, _strand=_codon.strand.value,
                                    _base_pair_representation=_codon.base_pair_representation, 
                                    _amino_acid_residue=_codon.amino_acid_residue, 
                                    _uniprot_position=_codon.uniprot_position,
                                    _chr=_codon.chr,
                                    _exon_number_base_pair_one=_codon.exon_number_base_pair_one,
                                    _exon_number_base_pair_two=_codon.exon_number_base_pair_two,
                                    _exon_number_base_pair_three=_codon.exon_number_base_pair_three,
                                    _chromosome_position_base_pair_one=_codon.chromosome_position_base_pair_one, 
                                    _chromosome_position_base_pair_two=_codon.chromosome_position_base_pair_two, 
                                    _chromosome_position_base_pair_three=_codon.chromosome_position_base_pair_three,
                                    _cDNA_position_one=_codon.cDNA_position_one, 
                                    _cDNA_position_two=_codon.cDNA_position_two, 
                                    _cDNA_position_three=_codon.cDNA_position_three,
                                    _variant_type=_variant_type, 
                                    _alt_amino_acid_residue=_alt_amino_acid_residue, 
                                    _ref_nucleotide=_ref_nucleotide, 
                                    _alt_nucleotide=_alt_nucleotide,
                                    _var_codon_position=_var_codon_position,
                                    _variant_source=_variant_source)
            
            return SNV
        except MalformedCodonException as e:
            raise MalformedCodonException("No SNV could be made: Malformed codon from dict: KeyError with message: "+str(e))
        except KeyError as e:
            raise MalformedVariantException("No SNV could be made: Malformed variant from dict: KeyError with message: "+str(e))
    
    @classmethod
    def initializeFromMapping(cls, _mappings, _gencode_transcription_id, _uniprot_ac):
        raise NotImplementedError("The function "'initializeFromMapping'" is not supported for SingleNucleotideVariant class")

    @classmethod
    def initializeFromVariant(cls, _codon, _chr_position, _alt_nucleotide, _variant_source):
        """Interprets the variant as a models.entities.SingleNucleotideVariant"""
        # retrieve information needed from the codon
        _variant_position_information = {}
        try:
            _variant_position_information = _codon.retrieve_mappings_per_chromosome()[_chr_position]
        except KeyError as e:
            raise MalformedVariantException("No SNV could be made: chromosome position '"+str(_chr_position)+"' does not exist in provided codon '"+str(_codon)+"' with error:"+str(e))
        
        # Retrieve the position and base pair of the provided codon
        _var_codon_position = _variant_position_information['codon_base_pair_position']
        _ref_nucleotide = _variant_position_information['base_pair']
        
        # interpret the variant in the context of the provided codon
        _alt_codon = SingleNucleotideVariant.interpret_alt_codon(_codon.base_pair_representation, _var_codon_position, _alt_nucleotide)
        _alt_amino_acid_residue = translate(_alt_codon)
        _variant_type = SingleNucleotideVariant.interpret_variant_type_from_residues(_codon.amino_acid_residue, _alt_amino_acid_residue)
        
        # create the SNV
        return cls(_gencode_transcription_id=_codon.gencode_transcription_id,
                                _uniprot_ac=_codon.uniprot_ac, _strand=_codon.strand.value,
                                _base_pair_representation=_codon.base_pair_representation, 
                                _amino_acid_residue=_codon.amino_acid_residue, 
                                _uniprot_position=_codon.uniprot_position,
                                _chr=_codon.chr,
                                _exon_number_base_pair_one=_codon.exon_number_base_pair_one,
                                _exon_number_base_pair_two=_codon.exon_number_base_pair_two,
                                _exon_number_base_pair_three=_codon.exon_number_base_pair_three,
                                _chromosome_position_base_pair_one=_codon.chromosome_position_base_pair_one, 
                                _chromosome_position_base_pair_two=_codon.chromosome_position_base_pair_two, 
                                _chromosome_position_base_pair_three=_codon.chromosome_position_base_pair_three,
                                _cDNA_position_one=_codon.cDNA_position_one, 
                                _cDNA_position_two=_codon.cDNA_position_two, 
                                _cDNA_position_three=_codon.cDNA_position_three,
                                _variant_type=_variant_type.value, 
                                _alt_amino_acid_residue=_alt_amino_acid_residue, 
                                _ref_nucleotide=_ref_nucleotide, 
                                _alt_nucleotide=_alt_nucleotide,
                                _var_codon_position=_var_codon_position,
                                _variant_source=_variant_source)
    
    def toJson(self):
        json_entry = {}
        json_entry["alt"] = self.alt_nucleotide
        json_entry["alt_aa"] = self.alt_amino_acid_residue
        json_entry["alt_aa_triplet"] = self.alt_three_letter_amino_acid_residue()
        json_entry["alt_codon"] = self.alt_base_pair_representation
        json_entry["pos"] = self.get_genomic_position_of_change()
        json_entry["ref"] = self.ref_nucleotide
        json_entry["type"] = self.variant_type.value
        return json_entry
    
    def toClinVarJson(self, ClinVar_id):
        json_entry = self.toJson()
        # append ClinVar specific information
        json_entry["clinvar_ID"] = ClinVar_id
        return json_entry
    
    def toGnommADJson(self, allele_count, allele_number):
        json_entry = self.toJson()
        # append gnomAD specific information
        json_entry['allele_number'] = allele_number
        json_entry['allele_count'] = allele_count
        return json_entry
    
    def __repr__(self):
        return "<SingleNucleotideVariant(chr='%s', chr_positions='%s', strand='%s', ref_codon='%s', ref_aa='%s',)>" % (
                            self.chr, str(self.regions), self.strand, self.base_pair_representation, self.amino_acid_residue)
