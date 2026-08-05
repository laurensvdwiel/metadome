import unittest

import pandas as pd

from metadome.domain.models.entities.meta_domain import (
    MetaDomain,
    UnsupportedMetaDomainIdentifier,
    ConsensusPositionOutOfBounds,
)
from metadome.domain.models.entities.gene_region import GenomeBuild
from metadome.domain.models.entities.single_nucleotide_variant import VariantSource


class mock_MetaDomain(MetaDomain):
    @classmethod
    def mock_PF00907_metadomain_first_three_consensus_positions(cls):
        def codon_row(consensus_pos, idx, uniprot_ac, uniprot_position, transcript_id):
            base_offset = consensus_pos * 100000 + idx * 10
            return {
                "gencode_transcription_id": transcript_id,
                "uniprot_ac": uniprot_ac,
                "strand": "+",
                "base_pair_representation": "ACC",
                "amino_acid_residue": "T",
                "uniprot_position": uniprot_position,
                "chr": "chr15",
                "exon_number_base_pair_one": 1,
                "exon_number_base_pair_two": 1,
                "exon_number_base_pair_three": 1,
                "chromosome_position_base_pair_one": base_offset + 1,
                "chromosome_position_base_pair_two": base_offset + 2,
                "chromosome_position_base_pair_three": base_offset + 3,
                "cDNA_position_one": base_offset + 4,
                "cDNA_position_two": base_offset + 5,
                "cDNA_position_three": base_offset + 6,
                "consensus_pos": consensus_pos,
                "domain_id": "PF00907",
            }

        rows = []

        rows.append(codon_row(1, 0, "O43435", 111, "ENST00000329705.7"))
        rows.append(codon_row(2, 0, "O43435", 112, "ENST00000329705.7"))
        rows.append(codon_row(3, 0, "O43435", 113, "ENST00000329705.7"))

        rows.append(codon_row(1, 1, "Q8IWI9-4", 76, "ENST00000219905.7"))
        rows.append(codon_row(2, 1, "Q8IWI9-4", 77, "ENST00000219905.7"))
        rows.append(codon_row(3, 1, "Q8IWI9-4", 78, "ENST00000219905.7"))

        for idx in range(2, 13):
            rows.append(codon_row(1, idx, f"P00001-{idx}", 100 + idx, f"ENST00000{idx:05d}.1"))

        for idx in range(2, 15):
            rows.append(codon_row(2, idx, f"P00002-{idx}", 200 + idx, f"ENST00001{idx:05d}.1"))

        for idx in range(2, 17):
            rows.append(codon_row(3, idx, f"P00003-{idx}", 300 + idx, f"ENST00002{idx:05d}.1"))

        meta_domain_mapping = pd.DataFrame(rows)
        meta_domain_annotation_per_source = {VariantSource.clinvar: pd.DataFrame(), VariantSource.gnomad: pd.DataFrame()}

        return cls(
            domain_id="PF00907",
            genome_build=GenomeBuild.GRCh37,
            consensus_length=3,
            consensus_positions={1, 2, 3},
            n_instances=3,
            meta_domain_mapping=meta_domain_mapping,
            meta_domain_annotation_per_source=meta_domain_annotation_per_source,
        )


class Test_meta_domain(unittest.TestCase):
    def test_get_alignment_depth_for_consensus_position(self):
        mock_metadom = mock_MetaDomain.mock_PF00907_metadomain_first_three_consensus_positions()

        with self.assertRaises(ConsensusPositionOutOfBounds):
            mock_metadom.get_alignment_depth_for_consensus_position(4)

        with self.assertRaises(ConsensusPositionOutOfBounds):
            mock_metadom.get_alignment_depth_for_consensus_position(0)

        self.assertEqual(mock_metadom.get_alignment_depth_for_consensus_position(1), 13)
        self.assertEqual(mock_metadom.get_alignment_depth_for_consensus_position(2), 15)
        self.assertEqual(mock_metadom.get_alignment_depth_for_consensus_position(3), 17)

    def test_get_max_alignment_depth(self):
        mock_metadom = mock_MetaDomain.mock_PF00907_metadomain_first_three_consensus_positions()
        self.assertEqual(mock_metadom.get_max_alignment_depth(), 17)

    def test_invalid_domain_id(self):
        with self.assertRaises(UnsupportedMetaDomainIdentifier):
            MetaDomain.initializeFromDomainID("TEST", GenomeBuild.GRCh37)

    def test_get_codons_aligned_to_consensus_position(self):
        mock_metadom = mock_MetaDomain.mock_PF00907_metadomain_first_three_consensus_positions()

        self.assertEqual(len(mock_metadom.get_codons_aligned_to_consensus_position(1)), 13)
        self.assertEqual(len(mock_metadom.get_codons_aligned_to_consensus_position(2)), 15)
        self.assertEqual(len(mock_metadom.get_codons_aligned_to_consensus_position(3)), 17)

        with self.assertRaises(ConsensusPositionOutOfBounds):
            mock_metadom.get_codons_aligned_to_consensus_position(4)

        with self.assertRaises(ConsensusPositionOutOfBounds):
            mock_metadom.get_codons_aligned_to_consensus_position(0)

    def test_get_consensus_position_for_uniprot_position(self):
        mock_metadom = mock_MetaDomain.mock_PF00907_metadomain_first_three_consensus_positions()

        self.assertEqual(
            mock_metadom.get_consensus_positions_for_uniprot_position(uniprot_ac="O43435", uniprot_position=111),
            [1],
        )
        self.assertEqual(
            mock_metadom.get_consensus_positions_for_uniprot_position(uniprot_ac="O43435", uniprot_position=112),
            [2],
        )
        self.assertEqual(
            mock_metadom.get_consensus_positions_for_uniprot_position(uniprot_ac="O43435", uniprot_position=113),
            [3],
        )

        self.assertEqual(
            mock_metadom.get_consensus_positions_for_uniprot_position(uniprot_ac="Q8IWI9-4", uniprot_position=76),
            [1],
        )
        self.assertEqual(
            mock_metadom.get_consensus_positions_for_uniprot_position(uniprot_ac="Q8IWI9-4", uniprot_position=77),
            [2],
        )
        self.assertEqual(
            mock_metadom.get_consensus_positions_for_uniprot_position(uniprot_ac="Q8IWI9-4", uniprot_position=78),
            [3],
        )

        self.assertEqual(
            len(mock_metadom.get_consensus_positions_for_uniprot_position(uniprot_ac="TESTFAIL", uniprot_position=9999)),
            0,
        )

        malformed_row = {
            "gencode_transcription_id": "ENST00000219905.7",
            "uniprot_ac": "Q8IWI9-4",
            "strand": "+",
            "base_pair_representation": "ACC",
            "amino_acid_residue": "T",
            "uniprot_position": 77,
            "chr": "chr15",
            "exon_number_base_pair_one": 1,
            "exon_number_base_pair_two": 1,
            "exon_number_base_pair_three": 1,
            "chromosome_position_base_pair_one": 41961324,
            "chromosome_position_base_pair_two": 41961325,
            "chromosome_position_base_pair_three": 41961326,
            "cDNA_position_one": 232,
            "cDNA_position_two": 233,
            "cDNA_position_three": 234,
            "consensus_pos": 1,
            "domain_id": "PF00907",
        }
        mock_metadom.meta_domain_mapping = pd.concat(
            [mock_metadom.meta_domain_mapping, pd.DataFrame([malformed_row])],
            ignore_index=True,
        )

        self.assertEqual(
            mock_metadom.get_consensus_positions_for_uniprot_position(uniprot_ac="Q8IWI9-4", uniprot_position=77),
            [2, 1],
        )


if __name__ == "__main__":
    unittest.main()