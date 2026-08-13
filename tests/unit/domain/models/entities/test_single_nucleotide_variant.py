import unittest
from metadome.domain.models.entities.codon import MalformedCodonException, Codon
from metadome.domain.models.entities.single_nucleotide_variant import SingleNucleotideVariant,\
    MalformedVariantException, VariantType
from builtins import NotImplementedError

class mock_Codon(Codon):
    @classmethod
    def mock_Methionine(cls):
        _d = {
            'chr': 'chr1',
            'gencode_transcription_id': 'test_transcript',
            'cDNA_position_three': 3,
            'uniprot_position': 125,
            'chromosome_position_base_pair_three': 3,
            'cDNA_position_one': 1,
            'strand': '+',
            'uniprot_ac': 'test_protein_ac',
            'base_pair_representation': 'ATG',
            'chromosome_position_base_pair_two': 2,
            'cDNA_position_two': 2,
            'chromosome_position_base_pair_one': 1,
            'amino_acid_residue': 'M',
            'exon_number_base_pair_one': 1,
            'exon_number_base_pair_two': 1,
            'exon_number_base_pair_three': 1,
        }
        return super(mock_Codon, cls).initializeFromDict(_d)
    
class Test_SingleNucleotideVariant(unittest.TestCase):
    
    def test_interpret_variant_type_from_residues(self):
        self.assertTrue(SingleNucleotideVariant.interpret_variant_type_from_residues('M', 'I') == VariantType.missense)
        self.assertTrue(SingleNucleotideVariant.interpret_variant_type_from_residues('M', 'M') == VariantType.synonymous)
        self.assertTrue(SingleNucleotideVariant.interpret_variant_type_from_residues('M', '*') == VariantType.nonsense)    

    def test_interpret_variant_type_from_codon_basepair_representations(self):
        self.assertTrue(SingleNucleotideVariant.interpret_variant_type_from_codon_basepair_representations('ATG', 'ATA') == VariantType.missense)
        self.assertTrue(SingleNucleotideVariant.interpret_variant_type_from_codon_basepair_representations('ATA', 'ATC') == VariantType.synonymous)
        self.assertTrue(SingleNucleotideVariant.interpret_variant_type_from_codon_basepair_representations('TAC', 'TAA') == VariantType.nonsense)
    
    def test_interpret_alt_codon(self):
        self.assertTrue(SingleNucleotideVariant.interpret_alt_codon('ATG', 0, 'C') == 'CTG')
        self.assertTrue(SingleNucleotideVariant.interpret_alt_codon('ATG', 1, 'C') == 'ACG')
        self.assertTrue(SingleNucleotideVariant.interpret_alt_codon('ATG', 2, 'C') == 'ATC')

    def test_initializations(self):
        _codon = mock_Codon.mock_Methionine()
        _variant_type = 'missense'
        _alt_amino_acid_residue='I'
        _ref_nucleotide='G'
        _alt_nucleotide='A'
        _var_codon_position=2
        _chromosome_position=3
        _variant_source = 'test_data'
        # should succeed
        _var_from_var = SingleNucleotideVariant.initializeFromVariant(_codon, _chromosome_position, _alt_nucleotide, _variant_source)
        _var_from_init = SingleNucleotideVariant(_gencode_transcription_id=_codon.gencode_transcription_id,
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
                                _ref_nucleotide=_ref_nucleotide, _alt_nucleotide=_alt_nucleotide,
                                _var_codon_position=_var_codon_position, _variant_source=_variant_source)
        
        # Check if the init went okay
        self.assertTrue(_var_from_init.alt_amino_acid_residue == _var_from_var.alt_amino_acid_residue == _alt_amino_acid_residue)
        self.assertTrue(_var_from_init.ref_nucleotide == _var_from_var.ref_nucleotide == _ref_nucleotide)
        self.assertTrue(_var_from_init.alt_nucleotide == _var_from_var.alt_nucleotide == _alt_nucleotide)
        self.assertTrue(_var_from_init.variant_type.value == _var_from_var.variant_type.value == _variant_type)
        self.assertTrue(_var_from_init.var_codon_position == _var_from_var.var_codon_position == _var_codon_position)
        self.assertTrue(_var_from_init.alt_base_pair_representation == _var_from_var.alt_base_pair_representation == 'ATA')
        self.assertTrue(_var_from_init.variant_source == _var_from_var.variant_source == _variant_source)
        
        # Create a dictionary from the variant
        _d = _var_from_init.toDict()
        _var_from_dict = SingleNucleotideVariant.initializeFromDict(_d)
        
        # Check if the conversion went okay and it is the same as init
        self.assertTrue(_var_from_dict.alt_amino_acid_residue == _alt_amino_acid_residue)
        self.assertTrue(_var_from_dict.ref_nucleotide == _ref_nucleotide)
        self.assertTrue(_var_from_dict.alt_nucleotide == _alt_nucleotide)
        self.assertTrue(_var_from_dict.variant_type.value == _variant_type)
        self.assertTrue(_var_from_dict.var_codon_position == _var_codon_position)
        self.assertTrue(_var_from_dict.alt_base_pair_representation == 'ATA')
        self.assertTrue(_var_from_dict.variant_source == _variant_source)
        
    def test_initializeFromVariantFailures(self):
        _codon = mock_Codon.mock_Methionine()
        _alt_nucleotide='A'
        _chromosome_position=4
        _variant_source = 'test_data'
        
        with self.assertRaises(MalformedVariantException):
            SingleNucleotideVariant.initializeFromVariant(_codon, _chromosome_position, _alt_nucleotide, _variant_source)
        
        
    def test_initializeFromDictFailures(self):
        _codon = mock_Codon.mock_Methionine()
        
        # check if a malformed codon exception is raised
        with self.assertRaises(MalformedCodonException):
            SingleNucleotideVariant.initializeFromDict({})
            
        # check if a malformed variant is raised with just the codon values
        with self.assertRaises(MalformedVariantException):
            SingleNucleotideVariant.initializeFromDict(_codon.toDict())

    def test_initialization_fails(self):
        _codon = mock_Codon.mock_Methionine()

        #raise MalformedVariantException("No SNV could be made: found identical ref and alt nucleotides")
        with self.assertRaises(MalformedVariantException):
            SingleNucleotideVariant(_gencode_transcription_id=_codon.gencode_transcription_id,
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
                                    _variant_type='missense',
                                    _alt_amino_acid_residue='I',
                                    _ref_nucleotide='G', _alt_nucleotide='G',
                                    _var_codon_position=2, _variant_source='test_data')

        #raise MalformedVariantException("No SNV could be made: The ref '"+str(_ref_nucleotide)+"' and alt '"+str(_alt_nucleotide)+"' nucleotides should both have total length '1'")
        with self.assertRaises(MalformedVariantException):
            _test = SingleNucleotideVariant(_gencode_transcription_id=_codon.gencode_transcription_id,
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
                                    _variant_type='missense',
                                    _alt_amino_acid_residue='I',
                                    _ref_nucleotide='TG', _alt_nucleotide='A',
                                    _var_codon_position=2, _variant_source='test_data')
        with self.assertRaises(MalformedVariantException):
            _test = SingleNucleotideVariant(_gencode_transcription_id=_codon.gencode_transcription_id,
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
                                    _variant_type='missense',
                                    _alt_amino_acid_residue='I',
                                    _ref_nucleotide='G', _alt_nucleotide='TA',
                                    _var_codon_position=2, _variant_source='test_data')

        #raise MalformedVariantException("No SNV could be made: _var_codon_position, indicating the variant in the basepair position in the codon should be 0, 1, or 2 but was '"+str(_var_codon_position)+"'")
        with self.assertRaises(MalformedVariantException):
            _test = SingleNucleotideVariant(_gencode_transcription_id=_codon.gencode_transcription_id,
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
                                    _variant_type='missense',
                                    _alt_amino_acid_residue='I',
                                    _ref_nucleotide='G', _alt_nucleotide='A',
                                    _var_codon_position=4, _variant_source='test_data')

        #raise MalformedVariantException("No SNV could be made: the ref nucleoide '"+str(self.ref_nucleotide)+"' does not correspond to the var_codon_position '"+str(self.var_codon_position)+"' in the base_pair_representation '"+str(self.base_pair_representation)+"'")
        with self.assertRaises(MalformedVariantException):
            _test = SingleNucleotideVariant(_gencode_transcription_id=_codon.gencode_transcription_id,
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
                                    _variant_type='missense',
                                    _alt_amino_acid_residue='I',
                                    _ref_nucleotide='T', _alt_nucleotide='A',
                                    _var_codon_position=2, _variant_source='test_data')

        # raise MalformedVariantException('No SNV could be made: Variant type cold not be converted to domain.models.entities.SingleNucleotideVariant.VariantType(Enum) for provided: '+str(_variant_type))
        with self.assertRaises(MalformedVariantException):
            _test = SingleNucleotideVariant(_gencode_transcription_id=_codon.gencode_transcription_id,
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
                                    _variant_type='not a variant type',
                                    _alt_amino_acid_residue='I',
                                    _ref_nucleotide='G', _alt_nucleotide='A',
                                    _var_codon_position=2, _variant_source='test_data')

    def test_initializeFromMappings(self):
        with self.assertRaises(NotImplementedError):
            SingleNucleotideVariant.initializeFromMapping(None, None, None)

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.test_initializeFromVariant']
    unittest.main()