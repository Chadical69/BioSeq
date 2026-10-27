import json
import unittest

from dna_analysis import analyze_sequence
from export_utils import build_csv_report, build_fasta, build_json_report, build_text_report


class TestExports(unittest.TestCase):
    def setUp(self):
        self.result = analyze_sequence(">demo\nATGGCTTAA")

    def test_csv_sections(self):
        csv_report = build_csv_report(self.result)
        self.assertIn("gc_content_percent", csv_report)
        self.assertIn("codon_frequency", csv_report)

    def test_text_report(self):
        report = build_text_report(self.result)
        self.assertIn("BIOSEQ", report)
        self.assertIn("Protein sequence: MA*", report)

    def test_fasta_format(self):
        fasta = build_fasta("ATGGCTTAA")
        self.assertTrue(fasta.startswith(">BioSeq_sequence length=9\n"))
        self.assertIn("ATGGCTTAA", fasta)

    def test_json_report(self):
        payload = json.loads(build_json_report(self.result))
        self.assertEqual(payload["sequence"], "ATGGCTTAA")


if __name__ == "__main__":
    unittest.main()
