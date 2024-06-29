import unittest

from Bio import Entrez

import Reciprocal


class TestReciprocal(unittest.TestCase):
    def setUp(self):
        Entrez.email = 'example@gmail.com'
        self.test_accession = '4557757'
        self.expected_organism = "Homo sapiens"

    def test_best_reciprocal_blast_returns_expected_results(self):
        results = Reciprocal.best_reciprocal_blast(self.expected_organism,self.test_accession)
        self.assertEqual(1, len(results))


if __name__ == '__main__':
    unittest.main()
