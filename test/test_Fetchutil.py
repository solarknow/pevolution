import unittest
import Fetchutil
from Bio import Entrez

class TestFetchutil(unittest.TestCase):
    def setUp(self):
        Entrez.email = 'example@gmail.com'
        self.test_accession = 'AJA33470.1'
        self.expected_organism = ['Vibrio sp. AN61', 'Bacteria']
        self.expected_path = 'Orthos/AJA33470.1.fasta'
        self.expected_definition = 'MreB, partial [Vibrio sp. AN61].'
        self.expected_gi='734850974'
        
    def test_set_email_sets_email(self):
        Fetchutil.set_email('mdsarwade@gmail.com')
        self.assertEqual(Entrez.email, 'mdsarwade@gmail.com')
    
    def test_fetch_protein_returns_list(self):
        test_proteins = Fetchutil.fetch_protein(self.test_accession)
        self.assertIsInstance(test_proteins, list)
    
    def test_fetch_organism_returns_expected_list(self):
        test_organism = Fetchutil.fetch_organism(self.test_accession)
        self.assertAlmostEqual(test_organism, self.expected_organism)
        
    def test_fetch_fasta_returns_path(self):
        test_path = Fetchutil.fetch_fasta(self.test_accession)
        self.assertEqual(test_path, self.expected_path)
    
    def test_fetch_definition_returns_definition(self):
        test_definition = Fetchutil.fetch_definition(self.test_accession)
        self.assertEqual(test_definition, self.expected_definition)

    def test_toGI_returns_GI_number(self):
        test_gi = Fetchutil.toGI(self.test_accession)
        self.assertEqual(test_gi, self.expected_gi)
    
if __name__ == '__main__':
    unittest.main()