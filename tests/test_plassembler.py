'''
Unit tests for plassembler.

Usage: pytest

'''

# import
import unittest
import os
from pathlib import Path


# import functions
from src import input_commands


# move to folder with mock files. First try Github structure, then try pulled repository structure

test_data = Path("tests/test_data")
val_data =  os.path.join(test_data, 'validation')



class TestValidateFasta(unittest.TestCase):
    """ Tests of Fasta validation functions"""

    def test_non_fasta_input(self):
        with self.assertRaises(SystemExit):
            non_fasta_file =  os.path.join(val_data, 'test_not_fasta.txt')
            input_commands.validate_fasta(non_fasta_file)


# class test_extract_chromosome(unittest.TestCase):
#     """ Tests for extract_chromosome"""
#     @classmethod
#     def setUpClass(cls):
#         cls.logger = logging.getLogger('test_logger.log')
#         cls.logger.setLevel(logging.INFO)

#     def test_get_extract_chromosome_only_chromosome_case_one(self):
#         out_dir = 'case_one/output'

#         expected_return  = True

#         chrom_length = 2400000

#         no_plasmid_flag = True

#         return_object = plassemblerModules.extract_chromosome(out_dir, chrom_length, no_plasmid_flag)

#         self.assertEqual(return_object, expected_return)

#     def test_get_extract_chromosome_case_three(self):
#         out_dir = 'case_three/output'

#         expected_return  = True

#         chrom_length = 2400000

#         no_plasmid_flag = True

#         return_object = plassemblerModules.extract_chromosome(out_dir, chrom_length, no_plasmid_flag)

#         self.assertEqual(return_object, expected_return)

#     def test_get_extract_chromosome_no_assembly_insufficient_depth(self):
#         out_dir = 'insufficient_depth'

#         expected_return  = False

#         chrom_length = 2400000

#         no_plasmid_flag = True

#         return_object = plassemblerModules.extract_chromosome(out_dir, chrom_length, no_plasmid_flag)

#         self.assertEqual(return_object, expected_return)
    

# class test_validate_fastq(unittest.TestCase):
#     """ Test for validate_fastq"""
#     @classmethod
#     def setUpClass(cls):
#         cls.logger = logging.getLogger('test_logger.log')
#         cls.logger.setLevel(logging.INFO)

#     def test_validate_fastq_no_gzip(self):
#         fastq = 'test.fastq'

#         expected_return = False

#         return_object = plassemblerModules.validate_fastq(fastq)

#         self.assertEqual(return_object, expected_return)
    
#     def test_validate_fastq_gzip(self):
#         fastq = 'test_gzip.fastq.gz'

#         expected_return = True

#         return_object = plassemblerModules.validate_fastq(fastq)

#         self.assertEqual(return_object, expected_return)


# class test_get_contig_lengths(unittest.TestCase):
#     """ Test for get_contig_lengths"""
#     @classmethod
#     def setUpClass(cls):
#         cls.logger = logging.getLogger('test_logger.log')
#         cls.logger.setLevel(logging.INFO)

#     def test_get_contig_lengths_case_one(self):
#         out_dir = 'case_one/output'

#         expected_return  = {
#                 "chromosome": 2834724,
#                 "1": 2473
#                 }

#         return_object = plassemblerModules.get_contig_lengths(out_dir)

#         self.assertEqual(return_object, expected_return)

#     def test_get_contig_lengths_case_three(self):
#         out_dir = 'case_three/output'

#         expected_return  = {
#                 "chromosome": 2857100,
#                 "1": 29025
#                 }

#         return_object = plassemblerModules.get_contig_lengths(out_dir)

#         self.assertEqual(return_object, expected_return)

    
# class test_get_contig_circularity(unittest.TestCase):
#     """ Test for get_contig_circularity"""
#     @classmethod
#     def setUpClass(cls):
#         cls.logger = logging.getLogger('test_logger.log')
#         cls.logger.setLevel(logging.INFO)

#     def test_get_contig_circularity(self):
#         out_dir = 'case_one/output'

#         expected_return  = {
#                 "chromosome": 'circular',
#                 "1": 'circular'
#                 }

#         return_object = plassemblerModules.get_contig_circularity(out_dir)

#         self.assertEqual(return_object, expected_return)

# # tests the run_mash.py 

# class test_get_contig_count_run_mash(unittest.TestCase):
#     """ Test for the get_contig_count from run_mash.sh"""
#     @classmethod
#     def setUpClass(cls):
#         cls.logger = logging.getLogger('test_logger.log')
#         cls.logger.setLevel(logging.INFO)

#     def test_get_contig_count_run_mash(self):
#         plasmid_fasta = 'case_one/output/unicycler_output/assembly.fasta'

#         expected_return = 1

#         return_object = plassemblerModules.get_contig_count(plasmid_fasta)

#         self.assertEqual(return_object, expected_return)








# if __name__ == '__main__':
#     unittest.main()
