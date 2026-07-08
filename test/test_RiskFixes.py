import os
import tempfile
import unittest
from unittest.mock import patch

import Report
from Utils import FetchUtil, SeqUtil


class TestRiskFixes(unittest.TestCase):
    @patch("Report.subprocess.run")
    def test_consense_uses_subprocess_without_shell(self, mock_run):
        with tempfile.TemporaryDirectory() as temp_dir:
            input_path = os.path.join(temp_dir, "input.txt")
            with open(input_path, "w") as in_file:
                in_file.write("dummy")

            expected_tree_path = os.path.splitext(input_path)[0] + ".tre"
            with open(expected_tree_path, "w") as tree_file:
                tree_file.write("(A,B);")

            result = Report.consense(input_path)
            self.assertEqual(result, "(A,B);")
            mock_run.assert_called_once()
            _, kwargs = mock_run.call_args
            self.assertEqual(kwargs["check"], True)
            self.assertEqual(kwargs["text"], True)
            self.assertIn(input_path, kwargs["input"])

    @patch("Utils.SeqUtil.clustal_align")
    def test_remove_gaps_nexus_does_not_mutate_dict_during_iteration(self, mock_clustal_align):
        with tempfile.TemporaryDirectory() as temp_dir:
            in_path = os.path.join(temp_dir, "sample.nex")
            out_path = os.path.join(temp_dir, "out.nex")
            with open(in_path, "w") as handle:
                handle.write(
                    "#NEXUS\n"
                    "begin data;\n"
                    "dimensions ntax=3 nchar=4;\n"
                    "format datatype=protein interleave=no gap=-;\n"
                    "matrix\n"
                    "A A--A\n"
                    "B A--A\n"
                    "C AAAA\n"
                    ";\n"
                    "end;\n"
                )

            SeqUtil.splice_align(in_path, out_path)
            mock_clustal_align.assert_called_once()

    @patch("Utils.FetchUtil.fetch_protein")
    def test_fetch_organism_handles_definition_tokens_without_closing_bracket(self, mock_fetch_protein):
        mock_fetch_protein.return_value = [
            "DEFINITION\n",
            "some protein [Organism only\n",
            "ACCESSION\n",
        ]

        result = FetchUtil.fetch_organism("ACC")
        self.assertIsNone(result)


if __name__ == "__main__":
    unittest.main()
