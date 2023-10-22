"""
Unit tests for plassembler.

Usage: pytest

"""

import os
import shutil
import subprocess as sp
import sys
import unittest
from pathlib import Path

import pytest
from loguru import logger

from src.plassembler import begin_plassembler, end_plassembler
from src.plassembler.utils.concat import (
    concatenate_short_fastqs,
    concatenate_single_fasta,
    concatenate_single_fastq,
)
from src.plassembler.utils.depth import (
    concatenate_chrom_plasmids,
    get_contig_circularity,
    get_contig_lengths,
    get_depths_from_bam,
)

# import functions
from src.plassembler.utils.input_commands import (
    check_dependencies,
    validate_fasta,
    validate_fastas_assembled_mode,
    validate_fastq,
    validate_fastqs_assembled_mode,
    validate_pacbio_model,
)
from src.plassembler.utils.plass_class import Plass
from src.plassembler.utils.qc import copy_sr_fastq_file

# data
test_data = Path("tests/test_data")
val_data = Path(f"{test_data}/validation")
fake_out_dir = Path(f"{test_data}/fake_out_dir")
bad_dir = Path(f"{test_data}/bad_dir")
logdir = Path(f"{test_data}/logs")
map_dir = Path(f"{test_data}/map_dir")
assembly_class = Path(f"{test_data}/assembly_class")


# make fake tempdir for testing
@pytest.fixture(scope="session")
def tmp_dir(tmpdir_factory):
    return tmpdir_factory.mktemp("tmp")


def remove_directory(dir_path):
    if os.path.exists(dir_path):
        shutil.rmtree(dir_path)


# to ensure sys exit on logger error
logger.add(lambda _: sys.exit(1), level="ERROR")


class test_plassembler_misc(unittest.TestCase):
    """Tests for miscellaneous plassembler functions"""

    def test_end_plassembler(self):
        end_plassembler(23)

    def test_begin_plassembler_force(self):
        outdir: Path = f"{fake_out_dir}/test_out"
        begin_plassembler(outdir, True)
        remove_directory(outdir)

    def test_begin_plassembler_no_force_pass(self):
        outdir: Path = f"{fake_out_dir}/test_out"
        begin_plassembler(outdir, False)
        remove_directory(outdir)

    def test_begin_plassembler_no_force_fail(self):
        with self.assertRaises(SystemExit):
            begin_plassembler(fake_out_dir, False)


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


class TestInputCommands(unittest.TestCase):
    """Tests input commands"""

    # tests non-fasta input
    def test_non_fasta_input(self):
        with self.assertRaises(SystemExit):
            non_fasta_file = os.path.join(val_data, "test_not_fasta.txt")
            validate_fasta(non_fasta_file)

    # tests legit input
    def test_fasta_input(self):
        fasta_file = os.path.join(val_data, "test.fasta")
        # dummy to make sure the function works
        tmp = 1
        validate_fasta(fasta_file)
        self.assertEqual(tmp, 1)

    # assembled multifasta chrom
    def test_validate_fastas_assembled_mode_multiFASTA(self):
        with self.assertRaises(SystemExit):
            input_plasmids: Path = Path(val_data) / "test.fasta"
            input_chromosome: Path = Path(val_data) / "test_multi.fasta"
            no_copy_numbers = False
            validate_fastas_assembled_mode(input_chromosome, input_plasmids, no_copy_numbers)

    # assembled single chrom works fine
    def test_validate_fastas_assembled_mode_single(self):
        input_chromosome: Path = Path(val_data) / "test.fasta"
        input_plasmids: Path = Path(val_data) / "test.fasta"
        no_copy_numbers = False
        tmp = 1
        validate_fastas_assembled_mode(input_chromosome, input_plasmids, no_copy_numbers)
        self.assertEqual(tmp, 1)

    # no gzip
    def test_validate_fastq_no_gzip(self):
        test_fastq = os.path.join(val_data, "test.fastq")
        # zipped
        expected_return = False
        return_object = validate_fastq(test_fastq)
        self.assertEqual(return_object, expected_return)

    # gzip
    def test_validate_fastq_gzip(self):
        test_fastq = os.path.join(val_data, "test_2.fastq.gz")
        # zipped
        expected_return = True
        return_object = validate_fastq(test_fastq)
        self.assertEqual(return_object, expected_return)

    # fastq
    def test_validate_fastqs_fasta_as_fastq(self):
        with self.assertRaises(ValueError):
            fasta = os.path.join(val_data, "test.fasta")
            validate_fastq(fasta)

    # assembled fastq
    def test_validate_fastqs_assembled_mode_fasta_as_fastq(self):
        with self.assertRaises(ValueError):
            fasta = os.path.join(val_data, "test.fasta")
            s1 = os.path.join(val_data, "test_2.fastq.gz")
            s2 = os.path.join(val_data, "test.fastq")
            validate_fastqs_assembled_mode(fasta, s1, s2)

    # bad pacbio model
    def test_validate_pacbio_model_bad(self):
        with self.assertRaises(SystemExit):
            pacbio_model = "not_a_model"
            validate_pacbio_model(pacbio_model)

    # bad pacbio model
    def test_deps(self):
        expected_return = True
        check_dependencies()
        self.assertEqual(expected_return, True)


class test_concat(unittest.TestCase):
    """Test for concat.py"""

    # concat single
    def test_concat_single_fastq_bad(self):
        with self.assertRaises(ValueError):
            fasta1: Path = Path(val_data) / "test.fasta"
            fasta2: Path = Path(val_data) / "test.fasta"
            out_f: Path = Path(val_data) / "concat.fastq"
            concatenate_single_fastq(fasta1, fasta2, out_f)

    # concat single good
    def test_concatenate_single_fastq_good(self):
        f1: Path = Path(val_data) / "test.fastq"
        f2: Path = Path(val_data) / "test.fastq"
        out_f = os.path.join(val_data, "fastq")
        expected_return = True
        concatenate_single_fastq(f1, f2, out_f)
        self.assertEqual(expected_return, True)

    # concat single good
    def test_concatenate_single_fasta_good(self):
        f1: Path = Path(val_data) / "test.fasta"
        f2: Path = Path(val_data) / "test.fasta"
        out_f = os.path.join(val_data, "concat.fasta")
        expected_return = True
        concatenate_single_fasta(f1, f2, out_f)
        self.assertEqual(expected_return, True)

    def test_concatenate_short_fastqs_bad_dir(self):
        with self.assertRaises(SystemExit):
            concatenate_short_fastqs(val_data)


class test_qc(unittest.TestCase):
    """Test for qc.py"""

    def test_copy_sr_fastq_file_not_fastq(self):
        with self.assertRaises(SystemExit):
            infile: Path = Path(val_data) / "test.fasta"
            outfile: Path = Path(val_data) / "test2.fasta"
            copy_sr_fastq_file(infile, outfile)


class test_depth(unittest.TestCase):
    """Test for depth.py"""

    # concat single
    def test_concatenate_chrom_plasmids_wrong_dir(self):
        with self.assertRaises(SystemExit):
            concatenate_chrom_plasmids(val_data)

    # good
    def test_get_contig_lengths_good_dir(self):
        expected_return = True
        fasta: Path = Path(f"{val_data}/combined.fasta")
        get_contig_lengths(fasta)
        self.assertEqual(expected_return, True)

    # bad dir
    def test_get_contig_lengths_bad_dir(self):
        fasta: Path = Path(f"{bad_dir}/combined.fasta")
        with self.assertRaises(FileNotFoundError):
            get_contig_lengths(fasta)

    def test_get_contig_circularity_good_dir(self):
        expected_return = True
        fasta: Path = Path(f"{val_data}/combined.fasta")
        get_contig_circularity(fasta)
        self.assertEqual(expected_return, True)

    # bad dir
    def test_get_contig_circularity_bad_dir(self):
        fasta: Path = Path(f"{bad_dir}/combined.fasta")
        with self.assertRaises(FileNotFoundError):
            get_contig_circularity(fasta)

    def test_get_get_depths_from_bam_unsorted_error(self):
        bam_file: Path = Path(f"{map_dir}/short_read.bam")
        fasta: Path = Path(f"{map_dir}/combined.fasta")
        contig_lengths = get_contig_lengths(fasta)
        with self.assertRaises(sp.CalledProcessError):
            get_depths_from_bam(bam_file, contig_lengths=contig_lengths)
