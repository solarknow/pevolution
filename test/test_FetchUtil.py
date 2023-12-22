import os
import unittest
import FetchUtil
from Bio import Entrez
from glob import glob


class TestFetchUtil(unittest.TestCase):
    def setUp(self):
        Entrez.email = 'example@gmail.com'
        self.test_accession = 'AJA33470.1'
        self.expected_organism = ['Vibrio sp. AN61', 'Bacteria']
        self.expected_path = 'Orthos' + os.sep + 'AJA33470.1.fasta'
        self.expected_definition = 'MreB, partial [Vibrio sp. AN61].'

    def test_set_email_sets_email(self):
        FetchUtil.set_email('mdsarwade@gmail.com')
        self.assertEqual(Entrez.email, 'mdsarwade@gmail.com')

    def test_fetch_protein_returns_list(self):
        test_proteins = FetchUtil.fetch_protein(self.test_accession)
        self.assertIsInstance(test_proteins, list)
        self.assertGreater(len(test_proteins), 0)

    def test_fetch_organism_returns_expected_list(self):
        test_organism = FetchUtil.fetch_organism(self.test_accession)
        self.assertListEqual(test_organism, self.expected_organism)

    def test_fetch_fasta_writes_file(self):
        FetchUtil.fetch_fasta(self.test_accession)
        self.assertIn(self.expected_path, glob('Orthos' + os.sep + '*'))

    def test_fetch_definition_returns_definition(self):
        test_definition = FetchUtil.fetch_definition(self.test_accession)
        self.assertEqual(test_definition, self.expected_definition)


if __name__ == '__main__':
    unittest.main()
