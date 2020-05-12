import unittest

import FetchUtil
import SeqUtil
from Bio import Entrez


class TestSeqUtil(unittest.TestCase):
    def setUp(self):
        Entrez.email = 'example@gmail.com'
        self.test_accession = 'AJA33470.1'
        self.expected_organism = ['Vibrio sp. AN61', 'Bacteria']
        self.expected_path = 'Orthos/AJA33470.1.fasta'
        self.expected_definition = 'MreB, partial [Vibrio sp. AN61].'

    def test_findlonglen_returns_longest_length(self):
        test_dict = {
            'length1': 15,
            'length10': 150,
            'length100': 1,
        }
        self.assertEqual(len('length100'), SeqUtil.findlonglen(test_dict))

    def test_fetch_protein_returns_list(self):
        test_proteins = FetchUtil.fetch_protein(self.test_accession)
        self.assertIsInstance(test_proteins, list)

    def test_fetch_organism_returns_expected_list(self):
        test_organism = FetchUtil.fetch_organism(self.test_accession)
        self.assertListEqual(test_organism, self.expected_organism)

    def test_fetch_fasta_returns_path(self):
        test_path = FetchUtil.fetch_fasta(self.test_accession)
        self.assertEqual(test_path, self.expected_path)

    def test_fetch_definition_returns_definition(self):
        test_definition = FetchUtil.fetch_definition(self.test_accession)
        self.assertEqual(test_definition, self.expected_definition)


if __name__ == '__main__':
    unittest.main()
