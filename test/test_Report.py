import os
import tempfile
import unittest
from unittest.mock import patch

import Report
from helpers.constants import AlignsPath, BayesPath, DataPath, MLPath, ProtPath, ReportsPath


def _write_file(path_obj, content):
    with open(str(path_obj), "w") as handle:
        handle.write(content)


class TestReport(unittest.TestCase):
    def setUp(self):
        self.original_cwd = os.getcwd()
        self.temp_dir = tempfile.TemporaryDirectory()
        os.chdir(self.temp_dir.name)

    def tearDown(self):
        os.chdir(self.original_cwd)
        self.temp_dir.cleanup()

    @patch("Report.consense", return_value="(MLTREE);")
    @patch("Report.FetchUtil.fetch_definition", side_effect=lambda acc: f"definition-for-{acc}")
    @patch("Report.FetchUtil.fetch_organism", return_value=["Query organism", "Bacteria"])
    def test_generate_report_prefers_con_file_and_writes_outputs(
        self, _mock_fetch_organism, _mock_fetch_definition, mock_consense
    ):
        dom_name = "all-demo"
        models = {"JTT+G": ["0.1", "0"]}

        _write_file(
            DataPath(f"{dom_name}.fas"),
            "Organism Alpha: Bacteria ACC1 1/1 1/1\nOrganism Beta: Archaea ACC2 1/1 1/1\n",
        )
        _write_file(
            AlignsPath(f"{dom_name}.best.nex"),
            "#NEXUS\nbegin data;\ndimensions ntax=2 nchar=4;\nformat datatype=protein interleave=no gap=-;\n"
            "matrix\nOrgA AAAA\nOrgB AAAA\n;\nend;\n",
        )
        _write_file(
            ProtPath(f"{dom_name}.pro"),
            "Header\nBest model according to BIC\nskip1\nskip2\nskip3\nskip4\nJTT+G 100 10 0.9 -1000\nWAG 250 20 0.1 -1200\n",
        )

        prefix = f"{dom_name}-JTT-ori"
        _write_file(MLPath(f"{prefix}_phyml_boot_trees.txt"), "tree1\n")

        _write_file(
            BayesPath(f"{dom_name}-bayes.nxs.con"),
            "translate\n1 OrgA,\n2 OrgB,\n;\n"
            "tree con_50 = [&R] (1[&prob=1.0],2[&prob=1.0]);\n",
        )
        _write_file(
            BayesPath(f"{dom_name}-bayes.nxs.con.tre"),
            "translate\n1 WrongA,\n2 WrongB,\n;\n"
            "tree con_50 = [&R] (1[&prob=1.0],2[&prob=1.0]);\n",
        )

        Report.generate_report("demo", "QUERY123", models, "all")

        report_path = str(ReportsPath(f"Report-{dom_name}.txt"))
        trees_path = str(ReportsPath(f"{dom_name}-trees.tre"))

        self.assertTrue(os.path.exists(report_path))
        self.assertTrue(os.path.exists(trees_path))

        with open(report_path) as report_file:
            report_text = report_file.read()

        self.assertIn("Reciprocal Best BLAST Results", report_text)
        self.assertIn("Query organism", report_text)
        self.assertIn("Model          deltaBIC*", report_text)
        self.assertIn("JTT+G 100", report_text)
        self.assertIn("Tree found by PhyML using the JTT model", report_text)
        self.assertIn("OrgA[&prob=1.0]", report_text)
        self.assertNotIn("WrongA[&prob=1.0]", report_text)

        mock_consense.assert_called_once_with(str(MLPath(f"{prefix}_phyml_boot_trees.txt")))

    @patch("Report.consense", return_value="(MLTREE);")
    @patch("Report.FetchUtil.fetch_definition", side_effect=lambda acc: f"definition-for-{acc}")
    @patch("Report.FetchUtil.fetch_organism", return_value=["Query organism", "Bacteria"])
    def test_generate_report_falls_back_to_con_tre_when_con_missing(
        self, _mock_fetch_organism, _mock_fetch_definition, _mock_consense
    ):
        dom_name = "bac-demo2"
        models = {"WAG": ["0", "0"]}

        _write_file(DataPath(f"{dom_name}.fas"), "Organism Gamma: Bacteria ACC9 1/1 1/1\n")
        _write_file(
            AlignsPath(f"{dom_name}.best.nex"),
            "#NEXUS\nbegin data;\ndimensions ntax=1 nchar=4;\nformat datatype=protein interleave=no gap=-;\n"
            "matrix\nOrgC AAAA\n;\nend;\n",
        )
        _write_file(
            ProtPath(f"{dom_name}.pro"),
            "Best model according to BIC\nskip1\nskip2\nskip3\nskip4\nWAG 100 10 0.9 -1000\nWAG2 250 20 0.1 -1200\n",
        )

        _write_file(
            BayesPath(f"{dom_name}-bayes.nxs.con.tre"),
            "translate\n1 OrgC,\n;\n"
            "tree con_50 = [&R] (1[&prob=1.0]);\n",
        )

        Report.generate_report("demo2", "QUERY999", models, "bac")

        report_path = str(ReportsPath(f"Report-{dom_name}.txt"))
        self.assertTrue(os.path.exists(report_path))

        with open(report_path) as report_file:
            report_text = report_file.read()

        self.assertIn("Tree found by MrBayes using the best model", report_text)
        self.assertIn("OrgC[&prob=1.0]", report_text)


if __name__ == "__main__":
    unittest.main()
