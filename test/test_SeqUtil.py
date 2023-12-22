import unittest

import SeqUtil
from Bio import Entrez


class TestSeqUtil(unittest.TestCase):
    def setUp(self):
        Entrez.email = 'example@gmail.com'
        self.test_accession = 'AJA33470.1'
        self.expected_organism = ['Vibrio sp. AN61', 'Bacteria']
        self.expected_path = 'Orthos/AJA33470.1.fasta'
        self.expected_definition = 'MreB, partial [Vibrio sp. AN61].'

    def test_find_longest_key_returns_longest_length(self):
        test_dict = {
            'length1': 15,
            'length10': 150,
            'length100': 1,
        }
        self.assertEqual(len('length100'), SeqUtil.find_longest_key(test_dict))

    def test_find_longest_key_returns_first_longest_length(self):
        test_dict = {
            'length1': 15,
            'length10': 150,
            'length100': 1,
            'length101': 2,
        }
        self.assertEqual(len('length100'), SeqUtil.find_longest_key(test_dict))


if __name__ == '__main__':
    unittest.main()
