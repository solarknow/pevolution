import os
import unittest

from Bio import Entrez

from Utils import FileUtil
from helpers.constants import DataPath

class TestFileUtil(unittest.TestCase):
    def setUp(self):
        Entrez.email = 'example@gmail.com'
        self.accession_dict = {
            'Gemmata obscuriglobus': ['WP_109571177.1', '1/1', '11/18'],
            'Prosthecobacter dejongeii': ['WP_246431095.1','1/1','5/18'],
            'Verrucomicrobium spinosum': ['WP_009961088.1','1/2','3/5'],
            'Agrobacterium tumefaciens': ['NTE00583.1','1/29','13/18'],
            'Klebsiella multispecies': ['WP_032692722.1','1/500','11/18']
        }
        self.expected_accessions = ['AAA17374.1', 'NP_000240.1','BAG35497.1','AAT44531.1','BAD96530.1']
        self.test_path = os.getcwd() + os.sep + 'test_files' + os.sep

    def test_parse_accession_numbers_from_XML(self):
        xml_path=self.test_path + 'sample.xml'
        self.assertEqual(FileUtil.XML_parse_and_extract_accession_numbers(xml_path), self.expected_accessions)

    def test_merging_multiple_fastas(self):
        FileUtil.merge_domain_fastas('bac_test.fas', self.accession_dict)
        self.assertEqual(os.path.getsize(str(DataPath('bac_test.fas'))), os.path.getsize(self.test_path+'sample.fas'))