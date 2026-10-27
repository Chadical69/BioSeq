import unittest

from dna_analysis import (
    analyze_sequence,
    complementary_strand,
    normalize_sequence,
    reverse_complement,
    sequence_statistics,
    transcribe_dna,
    translate_sequence,
    validate_sequence,
)


class TestBioSeq(unittest.TestCase):
    def test_normalize_fasta(self):
        self.assertEqual(normalize_sequence(">demo\natg c\nGTA\n"), "ATGCGTA")

    def test_validation(self):
        self.assertEqual(validate_sequence("ATGCN"), [])
        self.assertTrue(validate_sequence("ATGX"))
        self.assertTrue(validate_sequence(""))

    def test_statistics(self):
        stats = sequence_statistics("AATTGGCCN")
        self.assertEqual(stats["sequence_length"], 9)
        self.assertEqual(stats["canonical_bases"], 8)
        self.assertEqual(stats["ambiguous_bases"], 1)
        self.assertAlmostEqual(stats["gc_content"], 50.0)
        self.assertAlmostEqual(stats["at_content"], 50.0)

    def test_complements(self):
        self.assertEqual(complementary_strand("ATGC"), "TACG")
        self.assertEqual(reverse_complement("ATGC"), "GCAT")
        self.assertEqual(complementary_strand("RYSWKMBDHVN"), "YRSWMKVHDBN")

    def test_transcription(self):
        self.assertEqual(transcribe_dna("ATGC"), "AUGC")
        self.assertEqual(transcribe_dna("TACG", "Template strand (3′→5′)"), "AUGC")

    def test_translation(self):
        result = translate_sequence("ATGGCC TAA".replace(" ", ""), 1)
        self.assertEqual(result["protein_sequence"], "MA*")
        self.assertEqual(result["start_codon_positions"], [1])
        self.assertEqual(result["stop_codon_positions"], [7])

    def test_analyze_rejects_invalid(self):
        with self.assertRaises(ValueError):
            analyze_sequence("ATGX")


if __name__ == "__main__":
    unittest.main()
