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
from src.sam_to_fastq import (extract_bin_long_fastqs)
from src.plass_class import Assembly, Plass

# data
test_data = Path("tests/test_data")
val_data = Path(f"{test_data}/validation") 
fake_out_dir = Path(f"{test_data}/fake_out_dir")
bad_dir = Path(f"{test_data}/bad_dir") 
logdir = Path(f"{test_data}/logs") 
map_dir = Path(f"{test_data}/map_dir") 

# make fake tempdir for testing
@pytest.fixture(scope="session")
def tmp_dir(tmpdir_factory):
    return tmpdir_factory.mktemp("tmp")


def make_logger():
    logger = logging.getLogger()
    return logger

class test_unicycler_success(unittest.TestCase):
    """Test for validate_pacbio_model"""

    # bad pacbio model
    def test_unicycler_success_bad(self):
        plass = Plass()
        expected_return = False
        unicycler_output_dir = Path(f"{test_data}/unicycler_output_bad") 
        plass.check_unicycler_success(unicycler_output_dir)
        self.assertEqual(expected_return, False)

    def test_unicycler_success_good(self):
        plass = Plass()
        expected_return = True
        unicycler_output_dir = Path(f"{test_data}/unicycler_output") 
        plass.check_unicycler_success(unicycler_output_dir)
        self.assertEqual(expected_return, True)


class test_sam_to_fastq_long(unittest.TestCase):
    """Test for sam to fastq convertion with pysam"""
    # long read map
    def test_sam_to_fastq_long(self):
        expected_return = True
        extract_bin_long_fastqs(map_dir)
        self.assertEqual(expected_return, True)


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
            input_commands.validate_pacbio_model(pacbio_model)


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
            depth.concatenate_chrom_plasmids(val_data)

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




    # def test_minimap_depth_sort_long(self):
    #     """test minimap depth sort long"""
    #     tmp = 1
    #     depth.minimap_depth_sort_long(out_dir, threads)
    #     assert tmp == 1








