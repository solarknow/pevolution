import os
import unittest

from Bio import Entrez

import Reciprocal
from helpers.constants import DictsPath, OrthosPath, XMLPath


class TestReciprocal(unittest.TestCase):
    def setUp(self):
        Entrez.email = "example@gmail.com"
        self.test_accession = "4557757"
        self.expected_organism = "Homo sapiens"
        self.mock_db_name = (
            os.sep.join([os.getcwd(), "test_files", "mock_db"])
            if os.getcwd().split(os.sep)[-1] == "test"
            else os.sep.join([os.getcwd(), "test", "test_files", "mock_db"])
        )

    def tearDown(self):
        if os.path.exists(str(DictsPath(self.test_accession))):
            os.remove(str(OrthosPath(self.test_accession + ".fasta")))
            os.remove(str(XMLPath(self.test_accession + "_".join(self.expected_organism.split()) * 2 + ".xml")))
            os.remove(str(DictsPath(self.test_accession)))

    def test_best_reciprocal_blast_returns_expected_results(self):
        results = Reciprocal.best_reciprocal_blast(self.expected_organism, self.test_accession, db=self.mock_db_name)
        self.assertEqual(1, len(results))


if __name__ == "__main__":
    unittest.main()
