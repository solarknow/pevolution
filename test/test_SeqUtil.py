import json
import os
import unittest

from Bio import Entrez

import SeqUtil
from constants import DICTS_PATH, ORTHOS_PATH


class TestSeqUtil(unittest.TestCase):
    def setUp(self):
        Entrez.email = 'example@gmail.com'
        self.test_dict = {
            'length1': 15,
            'length10': 150,
            'length100': 1,
        }
        self.test_seq = '''>Hvol1
MSEAAESGDAPAGREIWIEKYRPQTFDDVYGQDDIVERLRSYIERDDLPHLLFAGPAGVG
KTTSATAIARAIYGDDWRGNFLELNASDERGIDVVRDRIKNFARSSFGGHDYRVIFLDEA
DSLTNDAQSALRRTMEQFSDNTRFILSCNYSSKIIDPIQSRCAVFRFSPLGDDAIAEQVR
DIAAAEDIEVTEDGLDALVYAAGGDMRRAINSLQAAATTGEVVDEEAVYMITSTARPEDI
EEMVRAAIDGEFTAARKQLETLIVDTGMAGGDIIDQLHRSVWEFDLDERDAVRLMERIGE
ADYRISEGANEQVQLEALLASLALSQN
'''
        self.expected_organism = ['Vibrio sp. AN61', 'Bacteria']
        self.expected_path = ORTHOS_PATH + 'TVT88608.fasta'
        self.expected_header = '>Hvol\n'
        self.test_dict_path = DICTS_PATH + 'test'

    def test_longest_key_length_returns_longest_length(self):
        self.assertEqual(len('length100'), SeqUtil.longest_key_length(self.test_dict))

    def test_dict_extract_returns_expected_dictionary(self):
        if os.path.exists(self.test_dict_path):
            self.assertEqual(SeqUtil.dict_extract(self.test_dict_path), self.test_dict)
        else:
            os.makedirs(DICTS_PATH, exist_ok=True)
            with open(self.test_dict_path, 'w') as out:
                out.write(json.dumps(self.test_dict) + '\n')
            self.assertEqual(SeqUtil.dict_extract(self.test_dict_path), self.test_dict)

    def test_shorten_sequence_names_returns_expected_header(self):
        test_path = self.expected_path.split('.')[0]
        os.remove(test_path)
        SeqUtil.shorten_sequence_names(self.expected_path)
        with open(self.expected_path.split('.')[0]) as o:
            firstline = o.readline()
        self.assertEqual(firstline, self.expected_header)


if __name__ == '__main__':
    unittest.main()
