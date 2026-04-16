import unittest
from builtins import NotImplementedError

from metadome.domain.models.entities.codon import MalformedCodonException, Codon
from metadome.domain.models.entities.single_nucleotide_variant import SingleNucleotideVariant, \
    MalformedVariantException, VariantType


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
        _alt_amino_acid_residue = 'I'
        _ref_nucleotide = 'G'
        _alt_nucleotide = 'A'
        _var_codon_position = 2
        _chromosome_position = 3
        _variant_source = 'test_data'

        _var_from_var = SingleNucleotideVariant.initializeFromVariant(_codon, _chromosome_position, _alt_nucleotide, _variant_source)
        _var_from_init = SingleNucleotideVariant(
            _gencode_transcription_id=_codon.gencode_transcription_id,
            _uniprot_ac=_codon.uniprot_ac,
            _strand=_codon.strand.value,
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
            _variant_source=_variant_source
        )

        self.assertTrue(_var_from_init.alt_amino_acid_residue == _var_from_var.alt_amino_acid_residue == _alt_amino_acid_residue)
        self.assertTrue(_var_from_init.ref_nucleotide == _var_from_var.ref_nucleotide == _ref_nucleotide)
        self.assertTrue(_var_from_init.alt_nucleotide == _var_from_var.alt_nucleotide == _alt_nucleotide)
        self.assertTrue(_var_from_init.variant_type.value == _var_from_var.variant_type.value == _variant_type)
        self.assertTrue(_var_from_init.var_codon_position == _var_from_var.var_codon_position == _var_codon_position)
        self.assertTrue(_var_from_init.alt_base_pair_representation == _var_from_var.alt_base_pair_representation == 'ATA')
        self.assertTrue(_var_from_init.variant_source == _var_from_var.variant_source == _variant_source)

        _d = _var_from_init.toDict()
        _var_from_dict = SingleNucleotideVariant.initializeFromDict(_d)

        self.assertTrue(_var_from_dict.alt_amino_acid_residue == _alt_amino_acid_residue)
        self.assertTrue(_var_from_dict.ref_nucleotide == _ref_nucleotide)
        self.assertTrue(_var_from_dict.alt_nucleotide == _alt_nucleotide)
        self.assertTrue(_var_from_dict.variant_type.value == _variant_type)
        self.assertTrue(_var_from_dict.var_codon_position == _var_codon_position)
        self.assertTrue(_var_from_dict.alt_base_pair_representation == 'ATA')
        self.assertTrue(_var_from_dict.variant_source == _variant_source)

    def test_initializeFromVariantFailures(self):
        _codon = mock_Codon.mock_Methionine()
        _alt_nucleotide = 'A'
        _chromosome_position = 4
        _variant_source = 'test_data'

        with self.assertRaises(MalformedVariantException):
            SingleNucleotideVariant.initializeFromVariant(_codon, _chromosome_position, _alt_nucleotide, _variant_source)

    def test_initializeFromDictFailures(self):
        _codon = mock_Codon.mock_Methionine()

        with self.assertRaises(MalformedCodonException):
            SingleNucleotideVariant.initializeFromDict({})

        with self.assertRaises(MalformedVariantException):
            SingleNucleotideVariant.initializeFromDict(_codon.toDict())

    def test_initialization_fails(self):
        _codon = mock_Codon.mock_Methionine()

        with self.assertRaises(MalformedVariantException):
            SingleNucleotideVariant(
                _gencode_transcription_id=_codon.gencode_transcription_id,
                _uniprot_ac=_codon.uniprot_ac,
                _strand=_codon.strand.value,
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
                _ref_nucleotide='G',
                _alt_nucleotide='G',
                _var_codon_position=2,
                _variant_source='test_data'
            )

        with self.assertRaises(MalformedVariantException):
            SingleNucleotideVariant(
                _gencode_transcription_id=_codon.gencode_transcription_id,
                _uniprot_ac=_codon.uniprot_ac,
                _strand=_codon.strand.value,
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
                _ref_nucleotide='TG',
                _alt_nucleotide='A',
                _var_codon_position=2,
                _variant_source='test_data'
            )

        with self.assertRaises(MalformedVariantException):
            SingleNucleotideVariant(
                _gencode_transcription_id=_codon.gencode_transcription_id,
                _uniprot_ac=_codon.uniprot_ac,
                _strand=_codon.strand.value,
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
                _ref_nucleotide='G',
                _alt_nucleotide='TA',
                _var_codon_position=2,
                _variant_source='test_data'
            )

        with self.assertRaises(MalformedVariantException):
            SingleNucleotideVariant(
                _gencode_transcription_id=_codon.gencode_transcription_id,
                _uniprot_ac=_codon.uniprot_ac,
                _strand=_codon.strand.value,
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
                _ref_nucleotide='G',
                _alt_nucleotide='A',
                _var_codon_position=4,
                _variant_source='test_data'
            )

        with self.assertRaises(MalformedVariantException):
            SingleNucleotideVariant(
                _gencode_transcription_id=_codon.gencode_transcription_id,
                _uniprot_ac=_codon.uniprot_ac,
                _strand=_codon.strand.value,
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
                _ref_nucleotide='T',
                _alt_nucleotide='A',
                _var_codon_position=2,
                _variant_source='test_data'
            )

        with self.assertRaises(MalformedVariantException):
            SingleNucleotideVariant(
                _gencode_transcription_id=_codon.gencode_transcription_id,
                _uniprot_ac=_codon.uniprot_ac,
                _strand=_codon.strand.value,
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
                _ref_nucleotide='G',
                _alt_nucleotide='A',
                _var_codon_position=2,
                _variant_source='test_data'
            )

    def test_initializeFromMappings(self):
        with self.assertRaises(NotImplementedError):
            SingleNucleotideVariant.initializeFromMapping(None, None, None)


if __name__ == "__main__":
    unittest.main()