"""
Unit tests for plassembler.

Usage: pytest

"""

import os
import shutil
import tempfile

# import
import unittest
from pathlib import Path
from unittest.mock import patch

import pytest

# import functions
from src.plassembler.utils.assembly import run_flye, run_raven
from src.plassembler.utils.bam import bam_to_fastq_short, sam_to_bam, split_bams
from src.plassembler.utils.cleanup import remove_directory, remove_file
from src.plassembler.utils.external_tools import ExternalTool
from src.plassembler.utils.mapping import minimap_long_reads, minimap_short_reads
from src.plassembler.utils.qc import chopper, fastp
from src.plassembler.utils.run_mash import get_contig_count, mash_sketch, run_mash
from src.plassembler.utils.run_unicycler import run_unicycler
from src.plassembler.utils.sam_to_fastq import (
    extract_long_fastqs_fast,
    extract_long_fastqs_slow_keep_fastqs,
)

test_data = Path("tests/test_data")
val_data = Path(f"{test_data}/validation")
fake_out_dir = Path(f"{test_data}/fake_out_dir")
bad_dir = Path(f"{test_data}/bad_dir")
logdir = Path(f"{test_data}/logs")
map_dir = Path(f"{test_data}/map_dir")
mash_dir = Path(f"{test_data}/mash_dir")
plassembler_db_dir = Path(f"{test_data}/Plassembler_Test_DB")


# make fake tempdir for testing
@pytest.fixture(scope="session")
def tmp_dir(tmpdir_factory):
    return tmpdir_factory.mktemp("tmp")


@pytest.mark.slow
class test_mash(unittest.TestCase):
    """Tests for run_mash.py"""

    # sam to bam
    def test_mash_sketch(self):
        expected_return = True
        with tempfile.TemporaryDirectory() as tmp:
            workdir = Path(tmp) / "mash_dir"
            shutil.copytree(mash_dir, workdir)
            fasta = workdir / "unicycler_plasmids.fasta"
            mash_sketch(workdir, fasta, logdir)
        self.assertEqual(expected_return, True)

    def test_run_mash(self):
        expected_return = True
        with tempfile.TemporaryDirectory() as tmp:
            workdir = Path(tmp) / "mash_dir"
            shutil.copytree(mash_dir, workdir)
            run_mash(workdir, plassembler_db_dir, logdir)
        self.assertEqual(expected_return, True)

    def test_get_contig_count(self):
        fasta: Path = Path(f"{mash_dir}/unicycler_plasmids.fasta")
        count = get_contig_count(fasta)
        self.assertEqual(count, 1)


@pytest.mark.slow
class test_bam(unittest.TestCase):
    """Tests for bam.py"""

    # sam to bam
    def test_sam_to_bam(self):
        expected_return = True
        threads = 1
        with tempfile.TemporaryDirectory() as tmp:
            workdir = Path(tmp) / "map_dir"
            shutil.copytree(map_dir, workdir)
            samfile = workdir / "sam_to_bam/test.sam"
            bamfile = workdir / "sam_to_bam/test.bam"
            sam_to_bam(samfile, bamfile, threads, logdir)
        self.assertEqual(expected_return, True)

    def test_split(self):
        expected_return = True
        with tempfile.TemporaryDirectory() as tmp:
            workdir = Path(tmp) / "map_dir"
            shutil.copytree(map_dir, workdir)
            split_bams(workdir, threads=1, logdir=logdir)
        self.assertEqual(expected_return, True)

    def test_bam_to_fastq_short(self):
        expected_return = True
        with tempfile.TemporaryDirectory() as tmp:
            workdir = Path(tmp) / "map_dir"
            shutil.copytree(map_dir, workdir)
            bam_to_fastq_short(workdir, threads=1, logdir=logdir)
        self.assertEqual(expected_return, True)


@pytest.mark.slow
class test_sam_to_fastq(unittest.TestCase):
    """Tests for sam_to_fastq.py"""

    # sam to bam
    def test_extract_long_fastqs_slow_keep_fastqs(self):
        expected_return = True
        with tempfile.TemporaryDirectory() as tmp:
            workdir = Path(tmp) / "map_dir"
            shutil.copytree(map_dir, workdir)
            samfile = workdir / "long_read.sam"
            plasmidfastq = workdir / "sam_to_bam/plasmid_long.fastq"
            outdir = workdir / "sam_to_bam"
            extract_long_fastqs_slow_keep_fastqs(outdir, samfile, plasmidfastq)
        self.assertEqual(expected_return, True)

    def test_extract_long_fastqs_fast(self):
        expected_return = True
        threads = 4
        with tempfile.TemporaryDirectory() as tmp:
            workdir = Path(tmp) / "map_dir"
            shutil.copytree(map_dir, workdir)
            samfile = workdir / "long_read.sam"
            plasmidfastq = workdir / "sam_to_bam/plasmid_long.fastq"
            extract_long_fastqs_fast(samfile, plasmidfastq, threads)
        self.assertEqual(expected_return, True)


@pytest.mark.slow
class test_mapping(unittest.TestCase):
    """Test for mapping"""

    # long read map
    def test_minimap_long_reads(self):
        expected_return = True
        pacbio_model = ""
        with tempfile.TemporaryDirectory() as tmp:
            workdir = Path(tmp) / "map_dir"
            shutil.copytree(map_dir, workdir)
            input_long_reads = workdir / "chopper_long_reads.fastq.gz"
            fasta = workdir / "flye_renamed.fasta"
            samfile = workdir / "test.sam"
            threads = 1
            minimap_long_reads(
                input_long_reads, fasta, samfile, threads, pacbio_model, logdir
            )
        self.assertEqual(expected_return, True)

        # short read map

    def test_minimap_short_reads(self):
        expected_return = True
        with tempfile.TemporaryDirectory() as tmp:
            workdir = Path(tmp) / "map_dir"
            shutil.copytree(map_dir, workdir)
            r1 = workdir / "trimmed_R1.fastq"
            r2 = workdir / "trimmed_R2.fastq"
            fasta = workdir / "flye_renamed.fasta"
            samfile = workdir / "test.sam"
            threads = 1
            minimap_short_reads(r1, r2, fasta, samfile, threads, logdir)
        self.assertEqual(expected_return, True)


@pytest.mark.slow
class test_qc_gzip(unittest.TestCase):
    """Test for qc"""

    # chopper
    def test_chopper_gzip(self):
        expected_return = True
        input_long_reads = os.path.join(test_data, "test_long.fastq.gz")
        chopper(
            input_long_reads, fake_out_dir, "500", "9", True, "1", logdir
        )  # True for gunzip
        remove_file(os.path.join(fake_out_dir, "chopper_long_reads.fastq.gz"))
        self.assertEqual(expected_return, True)

    def test_chopper_not_gzip(self):
        expected_return = True
        input_long_reads = os.path.join(test_data, "test_long.fastq")
        chopper(
            input_long_reads, fake_out_dir, "500", "9", False, "1", logdir
        )  # fasle for gunzip
        remove_file(os.path.join(fake_out_dir, "chopper_long_reads.fastq.gz"))
        self.assertEqual(expected_return, True)

    def test_fastp_gzip(self):
        expected_return = True
        short_one = Path(f"{test_data}/C11_subsetsim_R1.fastq.gz")
        short_two = Path(f"{test_data}/C11_subsetsim_R2.fastq.gz")
        fastp(short_one, short_two, fake_out_dir, logdir)
        remove_file(os.path.join(fake_out_dir, "trimmed_R1.fastq"))
        remove_file(os.path.join(fake_out_dir, "trimmed_R2.fastq"))
        remove_file("fastp.html")
        remove_file("fastp.json")
        self.assertEqual(expected_return, True)

    def test_fastp_nozip(self):
        expected_return = True
        short_one = Path(f"{test_data}/C11_subsetsim_R1.fastq")
        short_two = Path(f"{test_data}/C11_subsetsim_R2.fastq")
        fastp(short_one, short_two, fake_out_dir, logdir)
        remove_file(os.path.join(fake_out_dir, "trimmed_R1.fastq"))
        remove_file(os.path.join(fake_out_dir, "trimmed_R2.fastq"))
        remove_file("fastp.html")
        remove_file("fastp.json")
        self.assertEqual(expected_return, True)


@pytest.mark.slow
class test_assemblers(unittest.TestCase):
    """Test for assembles"""

    def test_flye(self):
        expected_return = True
        # C11 sim reads
        run_flye(test_data, 8, raw_flag=False, pacbio_model="nothing", logdir=logdir)
        shutil.rmtree(os.path.join(test_data, "00-assembly"))
        shutil.rmtree(os.path.join(test_data, "10-consensus"))
        shutil.rmtree(os.path.join(test_data, "20-repeat"))
        shutil.rmtree(os.path.join(test_data, "30-contigger"))
        shutil.rmtree(os.path.join(test_data, "40-polishing"))
        remove_file(os.path.join(test_data, "assembly.fasta"))
        remove_file(os.path.join(test_data, "assembly_info.txt"))
        remove_file(os.path.join(test_data, "assembly_graph.gfa"))
        remove_file(os.path.join(test_data, "assembly_graph.gv"))
        remove_file(os.path.join(test_data, "flye.log"))
        self.assertEqual(expected_return, True)

    def test_raven(self):
        expected_return = True
        # C11 sim reads
        run_raven(test_data, 1, logdir=logdir)
        remove_file(os.path.join(test_data, "assembly.fasta"))
        remove_file(os.path.join(test_data, "assembly_graph.gfa"))
        remove_file(os.path.join(test_data, "params.json"))
        remove_file("raven.cereal")
        self.assertEqual(expected_return, True)

    def test_unicycler_good(self):
        expected_return = True
        # C11 sim reads
        short_one = Path(f"{test_data}/short_read_concat_good_R1.fastq")
        short_two = Path(f"{test_data}/short_read_concat_good_R2.fastq")
        longreads = Path(f"{test_data}/plasmid_long_good.fastq")
        unicycler_output_dir = Path(f"{test_data}/unicycler_output")
        threads = 1
        run_unicycler(
            threads,
            logdir,
            short_one,
            short_two,
            longreads,
            unicycler_output_dir,
            unicycler_options=None,
            spades_options=None,
        )
        remove_directory(unicycler_output_dir)
        self.assertEqual(expected_return, True)

    def test_unicycler_bad(self):
        expected_return = True
        # C11 sim reads
        short_one = Path(f"{test_data}/C11_subsetsim_R1.fastq")
        short_two = Path(f"{test_data}/C11_subsetsim_R2.fastq")
        longreads = Path(f"{test_data}/plasmid_long_good.fastq")
        unicycler_output_dir = Path(f"{test_data}/unicycler_output_bad")
        threads = 1
        run_unicycler(
            threads,
            logdir,
            short_one,
            short_two,
            longreads,
            unicycler_output_dir,
            unicycler_options=None,
            spades_options=None,
        )
        remove_directory(unicycler_output_dir)
        self.assertEqual(expected_return, True)


class TestExternalTools:
    @patch.object(
        ExternalTool,
        ExternalTool._build_command.__name__,
        return_value=["mocked", "command", "arg"],
    )
    @patch.object(Path, Path.mkdir.__name__)
    def test___constructor(self, mkdir_mock, build_command_mock):
        logdir = Path("logs")

        external_tool = ExternalTool(
            "tool", "input", "output", "params", logdir, "outfile"
        )

        assert external_tool.command == ["mocked", "command", "arg"]
        assert external_tool.command_as_str == "mocked command arg"
        assert (
            external_tool.out_log
            == "logs/tool_c238863b32d18040bbf255fa3bf0dc91e9afa268335b56f51abe3c6d1fd83261.out"
        )
        assert (
            external_tool.err_log
            == "logs/tool_c238863b32d18040bbf255fa3bf0dc91e9afa268335b56f51abe3c6d1fd83261.err"
        )

        build_command_mock.assert_called_once_with("tool", "input", "output", "params")
        mkdir_mock.assert_called_once_with(parents=True, exist_ok=True)

    def test___build_command___simple_command(self):
        expected_escaped_command = ["tool", "param1", "param2", "-o", "out", "-i", "in"]
        actual_escaped_command = ExternalTool._build_command(
            "tool", "-i in", "-o out", "param1 param2"
        )
        assert expected_escaped_command == actual_escaped_command

    def test___build_command___single_quote_escaped(self):
        expected_escaped_command = [
            "tool",
            "params",
            "with",
            "escaped arg",
            "-o",
            "escaped out",
            "-i",
            "escaped in",
        ]
        actual_escaped_command = ExternalTool._build_command(
            "tool", "-i 'escaped in'", "-o 'escaped out'", "params with 'escaped arg'"
        )
        assert expected_escaped_command == actual_escaped_command

    def test___build_command___double_quote_escaped(self):
        expected_escaped_command = [
            "tool",
            "params",
            "with",
            "escaped arg",
            "-o",
            "escaped out",
            "-i",
            "escaped in",
        ]
        actual_escaped_command = ExternalTool._build_command(
            "tool", '-i "escaped in"', '-o "escaped out"', 'params with "escaped arg"'
        )
        assert expected_escaped_command == actual_escaped_command
