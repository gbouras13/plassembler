"""
Unit tests for plassembler.

Usage: pytest

"""

# import
import unittest
import os
from pathlib import Path
import pytest
import logging


# import functions
from src import input_commands
from src import concat
from src import depth


# move to folder with mock files. First try Github structure, then try pulled repository structure

test_data = Path("tests/test_data")
val_data = os.path.join(test_data, "validation")
fake_out_dir = os.path.join(test_data, "fake_out_dir")
bad_dir = os.path.join(test_data, "bad_dir")

# make fake tempdir for testing
@pytest.fixture(scope="session")
def tmp_dir(tmpdir_factory):
    return tmpdir_factory.mktemp("tmp")


def make_logger():
    logger = logging.getLogger()
    return logger


class TestValidateFasta(unittest.TestCase):
    """Tests of Fasta validation functions"""

    # tests non-fasta input
    def test_non_fasta_input(self):
        with self.assertRaises(SystemExit):
            non_fasta_file = os.path.join(val_data, "test_not_fasta.txt")
            input_commands.validate_fasta(non_fasta_file)

    # tests legit input
    def test_fasta_input(self):
        fasta_file = os.path.join(val_data, "test.fasta")
        # dummy to make sure the function works
        tmp = 1
        input_commands.validate_fasta(fasta_file)
        self.assertEqual(tmp, 1)

    # assembled multifasta chrom
    def test_validate_fastas_assembled_mode_multiFASTA(self):
        with self.assertRaises(SystemExit):
            input_plasmids = os.path.join(val_data, "test.fasta")
            input_chromosome = os.path.join(val_data, "test_multi.fasta")
            input_commands.validate_fastas_assembled_mode(
                input_chromosome, input_plasmids
            )

    # assembled single chrom works fine
    def test_validate_fastas_assembled_mode_single(self):
        input_plasmids = os.path.join(val_data, "test.fasta")
        input_chromosome = os.path.join(val_data, "test.fasta")
        tmp = 1
        input_commands.validate_fastas_assembled_mode(input_chromosome, input_plasmids)
        self.assertEqual(tmp, 1)


class test_validate_fastq(unittest.TestCase):
    """Test for validate_fastq"""

    # no gzip
    def test_validate_fastq_no_gzip(self):
        test_fastq = os.path.join(val_data, "test.fastq")
        # zipped
        expected_return = False
        return_object = input_commands.validate_fastq(test_fastq)
        self.assertEqual(return_object, expected_return)

    # gzip
    def test_validate_fastq_gzip(self):
        test_fastq = os.path.join(val_data, "test_2.fastq.gz")
        # zipped
        expected_return = True
        return_object = input_commands.validate_fastq(test_fastq)
        self.assertEqual(return_object, expected_return)

    # fastq
    def test_validate_fastqs_fasta_as_fastq(self):
        with self.assertRaises(ValueError):
            fasta = os.path.join(val_data, "test.fasta")
            input_commands.validate_fastq(fasta)

    # assembled fastq
    def test_validate_fastqs_assembled_mode_fasta_as_fastq(self):
        with self.assertRaises(ValueError):
            fasta = os.path.join(val_data, "test.fasta")
            s1 = os.path.join(val_data, "test_2.fastq.gz")
            s2 = os.path.join(val_data, "test.fastq")
            input_commands.validate_fastqs_assembled_mode(fasta, s1, s2)


class test_validate_pacbio_model(unittest.TestCase):
    """Test for validate_pacbio_model"""

    # bad pacbio model
    def test_validate_pacbio_model_bad(self):
        with self.assertRaises(SystemExit):
            pacbio_model = "not_a_model"
            logger = make_logger()
            input_commands.validate_pacbio_model(pacbio_model, logger)


class test_concat(unittest.TestCase):
    """Test for concat.py"""

    # concat single
    def test_concat_single_file_fasta(self):
        with self.assertRaises(ValueError):
            f1 = os.path.join(val_data, "test.fastq")
            fasta = os.path.join(val_data, "test.fasta")
            out_f = os.path.join(val_data, "concat.fastq")
            concat.concatenate_single(f1, fasta, out_f)

    # concat single good
    def test_concat_single_file_good(self):
        f1 = os.path.join(val_data, "test.fastq")
        fasta = os.path.join(val_data, "test.fastq")
        out_f = os.path.join(val_data, "concat.fastq")
        expected_return = True
        concat.concatenate_single(f1, fasta, out_f)
        self.assertEqual(expected_return, True)

    # bad pacbio model
    def test_concatenate_short_fastqs_bad_dir(self):
        with self.assertRaises(SystemExit):
            concat.concatenate_short_fastqs(val_data)

class test_depth(unittest.TestCase):
    """Test for concat.py"""

    # concat single
    def test_concatenate_chrom_plasmids_wrong_dir(self):
        with self.assertRaises(SystemExit):
            logger = make_logger()
            depth.concatenate_chrom_plasmids(val_data, logger)

    # good
    def test_get_contig_lengths_good_dir(self):
        expected_return = True
        depth.get_contig_lengths(val_data)
        self.assertEqual(expected_return, True)
    # bad dir
    def test_get_contig_lengths_bad_dir(self):
        with self.assertRaises(FileNotFoundError):
            depth.get_contig_lengths(bad_dir)

    def test_get_contig_circularity_good_dir(self):
        expected_return = True
        depth.get_contig_circularity(val_data)
        self.assertEqual(expected_return, True)

    # bad dir
    def test_get_contig_circularity_bad_dir(self):
        with self.assertRaises(FileNotFoundError):
            depth.get_contig_circularity(bad_dir)









    # # concat single good
    # def test_concat_single_file_good(self):
    #     f1 = os.path.join(val_data, "test.fastq")
    #     fasta = os.path.join(val_data, "test.fastq")
    #     out_f = os.path.join(val_data, "concat.fastq")
    #     expected_return = True
    #     concat.concatenate_single(f1, fasta, out_f)
    #     self.assertEqual(expected_return, True)

    # # bad pacbio model
    # def test_concatenate_short_fastqs_bad_dir(self):
    #     with self.assertRaises(SystemExit):
    #         concat.concatenate_short_fastqs(val_data)



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
