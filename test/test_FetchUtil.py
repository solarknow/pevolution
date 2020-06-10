import os
import unittest
from pathlib import Path

import FetchUtil
from Bio import Entrez
from glob import glob


class TestFetchUtil(unittest.TestCase):
    def setUp(self):
        self.test_email = 'example@gmail.com'
        self.test_accession = 'AJA33470'
        self.expected_organism = 'Vibrio sp. AN61'
        self.expected_domain = 'Bacteria'
        self.expected_path = 'Orthos' + os.sep + 'AJA33470.fasta'
        self.expected_definition = 'MreB, partial [Vibrio sp. AN61].'
        self.expected_seq = 'ANTLIYVKGQGIVLDEPSVVANRQDRVGSAKSVAAVGHAAKQMLGRTPGNISAIRPMKDG\n' + \
                            'VIADFYVTEKMLQHFIKQVHDNSILKPSPRVLVCVPCGSTQVERRAIRESALGAGAREVY\n' + \
                            'LIDEPMAAAIGAGLRVSEPTGSMVVDIGGGTTEVAVISLNGVVYSSSVRIGGDRFDEAVI\n' + \
                            'NYVRRNYGSLIGEATAEKIKHEIGSAYPGDEVQEIEVRGRNLAEGVPRSFSLNSNEILEA\n' + \
                            'LQEPLSGIVSAVMVALEQCPPELASDISENGMVLTGGGALLKDLDRLLMEETGIPVVIAE\n' + \
                            'DPLTCV\n'
        self.expected_xml_path = 'XML' + os.sep + 'test.xml'

    def test_set_email_sets_email(self):
        FetchUtil.set_email(self.test_email)
        self.assertEqual(Entrez.email, 'example@gmail.com')

    def test_fetch_protein_returns_expected_protein(self):
        FetchUtil.set_email(self.test_email)
        test_protein = FetchUtil.fetch_protein(self.test_accession)
        self.assertTrue(test_protein.accession == self.test_accession)
        self.assertTrue(test_protein.binomial == self.expected_organism)
        self.assertTrue(test_protein.domain == self.expected_domain)
        self.assertTrue(test_protein.definition == self.expected_definition)
        self.assertTrue(test_protein.seq == self.expected_seq)

    def test_manual_protein_is_fetched_protein(self):
        FetchUtil.set_email(self.test_email)
        expected_protein = FetchUtil.fetch_protein(self.test_accession)
        test_protein = FetchUtil.Protein.manual(self.expected_seq, self.expected_definition, self.expected_domain,
                                                self.expected_organism, self.test_accession)
        self.assertTrue(expected_protein == test_protein, str(expected_protein) + ' is not ' + str(test_protein))

    def test_fetch_fasta_writes_file(self):
        FetchUtil.set_email(self.test_email)
        test_protein = FetchUtil.fetch_protein(self.test_accession)
        FetchUtil.write_fasta(test_protein)
        self.assertIn(self.expected_path, glob('Orthos' + os.sep + '*'))

    def test_remote_blast_writes(self):
        FetchUtil.remote_blast(self.expected_path, 1e-10, self.expected_xml_path, self.expected_organism)
        self.assertTrue(os.path.exists(self.expected_xml_path) and Path(self.expected_xml_path).stat().st_size > 0)


if __name__ == '__main__':
    unittest.main()
