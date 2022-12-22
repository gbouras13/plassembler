'''
Unit tests for plassembler.
Usage: python -m unittest -v tests/plassembler_test.py 
'''

# import
import unittest
import os
import logging

# import  modules
import plassemblerModules 

# move to folder with mock files. 
# First try Github structure, then try pulled repository structure
try:
    os.chdir('/plassembler/tests/test_data/')
except FileNotFoundError:
    try:
        os.chdir('test_data/')
    except FileNotFoundError:
        os.chdir('tests/test_data/')


class test_contig_count(unittest.TestCase):
    """ Test for the contig_count"""
    @classmethod
    def setUpClass(cls):
        cls.logger = logging.getLogger('test_logger.log')
        cls.logger.setLevel(logging.INFO)

    def test_contig_count_case_one(self):
        out_dir = 'case_one/output'

        expected_return = 1

        return_object = plassemblerModules.contig_count(out_dir)

        self.assertEqual(return_object, expected_return)

class test_extract_chromosome(unittest.TestCase):
    """ Tests for extract_chromosome"""
    @classmethod
    def setUpClass(cls):
        cls.logger = logging.getLogger('test_logger.log')
        cls.logger.setLevel(logging.INFO)

    def test_get_extract_chromosome_only_chromosome_case_one(self):
        out_dir = 'case_one/output'

        expected_return  = True

        chrom_length = 2400000

        no_plasmid_flag = True

        return_object = plassemblerModules.extract_chromosome(out_dir, chrom_length, no_plasmid_flag)

        self.assertEqual(return_object, expected_return)

    def test_get_extract_chromosome_case_three(self):
        out_dir = 'case_three/output'

        expected_return  = True

        chrom_length = 2400000

        no_plasmid_flag = True

        return_object = plassemblerModules.extract_chromosome(out_dir, chrom_length, no_plasmid_flag)

        self.assertEqual(return_object, expected_return)

    def test_get_extract_chromosome_no_assembly_insufficient_depth(self):
        out_dir = 'insufficient_depth'

        expected_return  = False

        chrom_length = 2400000

        no_plasmid_flag = True

        return_object = plassemblerModules.extract_chromosome(out_dir, chrom_length, no_plasmid_flag)

        self.assertEqual(return_object, expected_return)
    

class test_validate_fastq(unittest.TestCase):
    """ Test for validate_fastq"""
    @classmethod
    def setUpClass(cls):
        cls.logger = logging.getLogger('test_logger.log')
        cls.logger.setLevel(logging.INFO)

    def test_validate_fastq_no_gzip(self):
        fastq = 'test.fastq'

        expected_return = False

        return_object = plassemblerModules.validate_fastq(fastq)

        self.assertEqual(return_object, expected_return)
    
    def test_validate_fastq_gzip(self):
        fastq = 'test_gzip.fastq.gz'

        expected_return = True

        return_object = plassemblerModules.validate_fastq(fastq)

        self.assertEqual(return_object, expected_return)


class test_get_contig_lengths(unittest.TestCase):
    """ Test for get_contig_lengths"""
    @classmethod
    def setUpClass(cls):
        cls.logger = logging.getLogger('test_logger.log')
        cls.logger.setLevel(logging.INFO)

    def test_get_contig_lengths_case_one(self):
        out_dir = 'case_one/output'

        expected_return  = {
                "chromosome": 2834724,
                "1": 2473
                }

        return_object = plassemblerModules.get_contig_lengths(out_dir)

        self.assertEqual(return_object, expected_return)

    def test_get_contig_lengths_case_three(self):
        out_dir = 'case_three/output'

        expected_return  = {
                "chromosome": 2857100,
                "1": 29025
                }

        return_object = plassemblerModules.get_contig_lengths(out_dir)

        self.assertEqual(return_object, expected_return)

    
class test_get_contig_circularity(unittest.TestCase):
    """ Test for get_contig_circularity"""
    @classmethod
    def setUpClass(cls):
        cls.logger = logging.getLogger('test_logger.log')
        cls.logger.setLevel(logging.INFO)

    def test_get_contig_circularity(self):
        out_dir = 'case_one/output'

        expected_return  = {
                "chromosome": 'circular',
                "1": 'circular'
                }

        return_object = plassemblerModules.get_contig_circularity(out_dir)

        self.assertEqual(return_object, expected_return)









if __name__ == '__main__':
    unittest.main()
