import unittest

from metadome.domain.models.entities.codon import Codon, MalformedCodonException
from metadome.domain.models.gene import Strand

class mock_Mapping(object):
    strand = Strand.plus
    uniprot_position = 125
    chromosome = 'chr1'
    gene_id = 1
    protein_id = 1
    
    def __init__(self, id, base_pair, codon, codon_base_pair_position,
                 cDNA_position, amino_acid_residue, chromosome_position):
        self.id = id
        self.base_pair = base_pair
        self.codon = codon
        self.codon_base_pair_position = codon_base_pair_position
        self.cDNA_position = cDNA_position
        self.amino_acid_residue = amino_acid_residue
        self.chromosome_position = chromosome_position

class Test_codon(unittest.TestCase):
    
    def test_init_from_dict_fail(self):
        _d = {}
        
        with self.assertRaises(MalformedCodonException):
            Codon.initializeFromDict(_d)
        
    def test_initializations(self):
        # init test variables
        _transcript =  'test_transcript'
        _protein_ac =  'test_test_protein_ac'
        _codon_repr = 'CTT'
        _residue = 'L'
        _strand = '+'
        _chromosome_position_base_pair_one = 231
        _chromosome_position_base_pair_two = 232
        _chromosome_position_base_pair_three = 233
        _cDNA_position_one = 333
        _cDNA_position_two = 334
        _cDNA_position_three = 335
        
        # init mappings
        _mappings = []
        _mappings.append(mock_Mapping(id=1, base_pair='C', codon=_codon_repr, codon_base_pair_position=0, amino_acid_residue=_residue, cDNA_position=_cDNA_position_one, chromosome_position=_chromosome_position_base_pair_one))
        _mappings.append(mock_Mapping(id=2, base_pair='T', codon=_codon_repr, codon_base_pair_position=1, amino_acid_residue=_residue, cDNA_position=_cDNA_position_two, chromosome_position=_chromosome_position_base_pair_two))
        _mappings.append(mock_Mapping(id=3, base_pair='T', codon=_codon_repr, codon_base_pair_position=2, amino_acid_residue=_residue, cDNA_position=_cDNA_position_three, chromosome_position=_chromosome_position_base_pair_three))
         
        # Create the Codon from the mapping
        _codon_from_mapping = Codon.initializeFromMapping(_mappings=_mappings, _gencode_transcription_id=_transcript, _uniprot_ac=_protein_ac)
        _codon_from_init = Codon(_gencode_transcription_id=_transcript, _uniprot_ac=_protein_ac,
                                 _strand=_strand, _base_pair_representation=_codon_repr,
                                 _amino_acid_residue=_residue, _uniprot_position=mock_Mapping.uniprot_position,
                                 _chr=mock_Mapping.chromosome, _chromosome_position_base_pair_one=_chromosome_position_base_pair_one,
                                 _chromosome_position_base_pair_two=_chromosome_position_base_pair_two,
                                 _chromosome_position_base_pair_three=_chromosome_position_base_pair_three,
                                 _cDNA_position_one=_cDNA_position_one, _cDNA_position_two=_cDNA_position_two,
                                 _cDNA_position_three=_cDNA_position_three)
        
        # Check if the conversion went okay and it is the same as init
        self.assertTrue(_codon_from_init.three_letter_amino_acid_residue() == _codon_from_mapping.three_letter_amino_acid_residue() == 'Leu')
        self.assertTrue(_codon_from_init.gencode_transcription_id == _codon_from_mapping.gencode_transcription_id == _transcript)
        self.assertTrue(_codon_from_init.uniprot_ac == _codon_from_mapping.uniprot_ac == _protein_ac)
        self.assertTrue(_codon_from_init.strand == _codon_from_mapping.strand == mock_Mapping.strand)
        self.assertTrue(_codon_from_init.base_pair_representation == _codon_from_mapping.base_pair_representation == _codon_repr) 
        self.assertTrue(_codon_from_init.amino_acid_residue == _codon_from_mapping.amino_acid_residue == _residue)
        self.assertTrue(_codon_from_init.uniprot_position == _codon_from_mapping.uniprot_position == mock_Mapping.uniprot_position)
        self.assertTrue(_codon_from_init.chr == _codon_from_mapping.chr == mock_Mapping.chromosome)
        self.assertTrue(_codon_from_init.chromosome_position_base_pair_one == _codon_from_mapping.chromosome_position_base_pair_one == _chromosome_position_base_pair_one)
        self.assertTrue(_codon_from_init.chromosome_position_base_pair_two == _codon_from_mapping.chromosome_position_base_pair_two == _chromosome_position_base_pair_two)
        self.assertTrue(_codon_from_init.chromosome_position_base_pair_three == _codon_from_mapping.chromosome_position_base_pair_three == _chromosome_position_base_pair_three)
        self.assertTrue(_codon_from_init.cDNA_position_one == _codon_from_mapping.cDNA_position_one == _cDNA_position_one)
        self.assertTrue(_codon_from_init.cDNA_position_two == _codon_from_mapping.cDNA_position_two == _cDNA_position_two)
        self.assertTrue(_codon_from_init.cDNA_position_three == _codon_from_mapping.cDNA_position_three == _cDNA_position_three)
        
        # Create a dictionary from the codon
        _d = _codon_from_mapping.toDict()
        _codon_from_dict = Codon.initializeFromDict(_d)
        
        # Check if the conversion went okay and it is the same as init
        self.assertTrue(_codon_from_dict.three_letter_amino_acid_residue() == _codon_from_mapping.three_letter_amino_acid_residue())
        self.assertTrue(_codon_from_dict.gencode_transcription_id == _codon_from_mapping.gencode_transcription_id)
        self.assertTrue(_codon_from_dict.uniprot_ac == _codon_from_mapping.uniprot_ac)
        self.assertTrue(_codon_from_dict.strand == _codon_from_mapping.strand)
        self.assertTrue(_codon_from_dict.base_pair_representation == _codon_from_mapping.base_pair_representation) 
        self.assertTrue(_codon_from_dict.amino_acid_residue == _codon_from_mapping.amino_acid_residue)
        self.assertTrue(_codon_from_dict.uniprot_position == _codon_from_mapping.uniprot_position)
        self.assertTrue(_codon_from_dict.chr == _codon_from_mapping.chr)
        self.assertTrue(_codon_from_dict.chromosome_position_base_pair_one == _codon_from_mapping.chromosome_position_base_pair_one)
        self.assertTrue(_codon_from_dict.chromosome_position_base_pair_two == _codon_from_mapping.chromosome_position_base_pair_two)
        self.assertTrue(_codon_from_dict.chromosome_position_base_pair_three == _codon_from_mapping.chromosome_position_base_pair_three)
        self.assertTrue(_codon_from_dict.cDNA_position_one == _codon_from_mapping.cDNA_position_one)
        self.assertTrue(_codon_from_dict.cDNA_position_two == _codon_from_mapping.cDNA_position_two)
        self.assertTrue(_codon_from_dict.cDNA_position_three == _codon_from_mapping.cDNA_position_three)        
        
    def test_three_letter_amino_acid_residue_selenocysteine(self):
        _mappings = []
        _mappings.append(mock_Mapping(id=1, base_pair='T', codon='TGA', codon_base_pair_position=0, amino_acid_residue='U', cDNA_position=1, chromosome_position=1))
        _mappings.append(mock_Mapping(id=2, base_pair='G', codon='TGA', codon_base_pair_position=1, amino_acid_residue='U', cDNA_position=2, chromosome_position=2))
        _mappings.append(mock_Mapping(id=3, base_pair='A', codon='TGA', codon_base_pair_position=2, amino_acid_residue='U', cDNA_position=3, chromosome_position=3))
         
        # Create the Codon
        codon = Codon.initializeFromMapping(_mappings, 'test_transcript', 'test_protein_ac')
         
        self.assertTrue(codon.three_letter_amino_acid_residue() == 'Sec')
        
    def test_three_letter_amino_acid_residue_pyrrolysine(self):

        _mappings = []
        _mappings.append(mock_Mapping(id=1, base_pair='T', codon='TAG', codon_base_pair_position=0, amino_acid_residue='U', cDNA_position=1, chromosome_position=1))
        _mappings.append(mock_Mapping(id=2, base_pair='A', codon='TAG', codon_base_pair_position=1, amino_acid_residue='U', cDNA_position=2, chromosome_position=2))
        _mappings.append(mock_Mapping(id=3, base_pair='G', codon='TAG', codon_base_pair_position=2, amino_acid_residue='U', cDNA_position=3, chromosome_position=3))
         
        # Create the Codon
        codon = Codon.initializeFromMapping(_mappings, 'test_transcript', 'test_protein_ac')
         
        self.assertTrue(codon.three_letter_amino_acid_residue() == 'Sec')

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()