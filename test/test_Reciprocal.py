import os
import shutil
import unittest

from Bio import Entrez

import Reciprocal
from helpers.constants import XMLPath, OrthosPath

ortho_dir = OrthosPath.dir
xml_dir = XMLPath.dir

class TestReciprocal(unittest.TestCase):
    def setUp(self):
        if not os.path.isdir(xml_dir):
            os.mkdir(xml_dir)
        if not os.path.isdir(ortho_dir):
            os.mkdir(ortho_dir)
        Entrez.email = 'example@gmail.com'
        self.test_accession = '4557757'
        self.expected_organism = "Homo sapiens"

    def tearDown(self):
        shutil.rmtree(ortho_dir)
        shutil.rmtree(xml_dir)

    def test_best_reciprocal_blast_returns_expected_results(self):
        results = Reciprocal.best_reciprocal_blast(self.expected_organism,self.test_accession)
        self.assertEqual(1, len(results))


if __name__ == '__main__':
    unittest.main()
