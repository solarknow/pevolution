import ast
import os
import shutil
import unittest

from Bio import Entrez

import SeqUtil

temp_dir = os.getcwd() + os.sep + 'Orthos'


class TestSeqUtil(unittest.TestCase):
    def setUp(self):
        if not os.path.isdir(temp_dir):
            os.mkdir(temp_dir)
        Entrez.email = 'example@gmail.com'
        self.test_accession = 'AJA33470.1'
        self.test2_accession = 'WP_105493738.1'
        self.expected_organism = ['Vibrio sp. AN61', 'Bacteria']
        self.expected_path = 'Orthos/AJA33470.1.fasta'
        self.expected_definition = 'MreB, partial [Vibrio sp. AN61].'
        self.test_path = os.getcwd() + os.sep + 'test_files' + os.sep

    def tearDown(self):
        shutil.rmtree(temp_dir)

    def test_findlonglen_returns_longest_length(self):
        test_dict = {
            'length1': 15,
            'length10': 150,
            'length100': 1,
        }
        self.assertEqual(len('length100'), SeqUtil.findlonglen(test_dict))

    def test_clusttofasta_returns_fasta(self):
        sample_fasta = self.test_path + 'sample.fasta'
        sample_clust = self.test_path + 'sample.aln'
        calculated_fasta = temp_dir + os.sep + 'out'
        SeqUtil.clusttofasta(sample_clust, calculated_fasta)
        self.assertEqual(SeqUtil.count_fasta_seqs(sample_fasta), SeqUtil.count_fasta_seqs(calculated_fasta))

    def test_dict_extract_parses_dictionary_from_file(self):
        sample_dict = {"name": "solarknow", 'age': 30, "alive": True}
        file_path = temp_dir+os.sep+'dict_obj'
        with open(file_path, 'w') as write_file:
            write_file.write(str(sample_dict))
        extracted = SeqUtil.dict_extract(file_path)
        self.assertDictEqual(sample_dict, extracted)

    def test_rename_fasta_renames_organism_correctly(self):
        org_fasta = self.test_path + "sample_organism.fasta"
        target_species = {'Ecol', 'Vsp.'}
        SeqUtil.rename_seqs(org_fasta, temp_dir + os.sep + 'sample_organism')
        with open(temp_dir + os.sep + 'sample_organism') as so:
            generated_species = {s[1:].strip() for s in so.readlines() if s.startswith('>')}
        self.assertSetEqual(target_species, generated_species)

    def test_removes_gaps_removes_gapped_positions(self):
        original_fasta = self.test_path + "sample.nexus"
        target_dict_file = self.test_path + "sample_degapped"
        with open(target_dict_file) as target:
            target_dict = ast.literal_eval(target.read())
        calculated_dict = SeqUtil.remove_gaps_nexus(original_fasta)
        self.assertDictEqual(target_dict, calculated_dict)

if __name__ == '__main__':
    unittest.main()
