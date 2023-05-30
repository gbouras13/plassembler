"""
Unit tests for plassembler.

Usage: pytest

"""

# import
import os
import shutil
import subprocess
import sys
import unittest
from pathlib import Path

import pytest
from loguru import logger

# import functions

# data
test_data = Path("tests/test_data")
end_to_end = Path(f"{test_data}/end_to_end")
plassembler_db_dir = Path(f"{test_data}/Plassembler_Test_DB")


# make fake tempdir for testing
@pytest.fixture(scope="session")
def tmp_dir(tmpdir_factory):
    return tmpdir_factory.mktemp("tmp")


# to ensure sys exit on logger error
logger.add(lambda _: sys.exit(1), level="ERROR")
# import functions


def remove_directory(dir_path):
    if os.path.exists(dir_path):
        shutil.rmtree(dir_path)


def exec_command(cmnd, stdout=subprocess.PIPE, stderr=subprocess.PIPE):
    """executes shell command and returns stdout if completes exit code 0
    Parameters
    ----------
    cmnd : str
      shell command to be executed
    stdout, stderr : streams
      Default value (PIPE) intercepts process output, setting to None
      blocks this."""

    proc = subprocess.Popen(cmnd, shell=True, stdout=stdout, stderr=stderr)
    out, err = proc.communicate()
    if proc.returncode != 0:
        raise RuntimeError(f"FAILED: {cmnd}\n{err}")
    return out.decode("utf8") if out is not None else None


def test_plassembler_download():
    """test plassembler download"""
    fake_db = f"{end_to_end}/db"
    cmd = f"plassembler download -d {fake_db}"
    exec_command(cmd)
    remove_directory(fake_db)


def test_citation():
    """test plassembler citation"""
    cmd = "plassembler citation"
    exec_command(cmd)


# test running end to end
# uncomment for mac running to check

# 70kbp, 44kbp and 9kbp plasmid reads are from
# the 70kbp is a fake chromosome


# def test_plassembler(tmp_dir):
#     """test plassembler run"""
#     longreads: Path = f"{end_to_end}/input_fastq.gz"
#     s1: Path = f"{end_to_end}/input_R1.fastq.gz"
#     s2: Path = f"{end_to_end}/input_R2.fastq.gz"
#     chromosome = 50000
#     outdir: Path = f"{end_to_end}/test_out"
#     cmd = f"plassembler run -l {longreads} -c {chromosome} -1 {s1} -2 {s2} -d {plassembler_db_dir} -o {outdir}  -t 8 -f"
#     exec_command(cmd)
#     remove_directory(outdir)


# def test_plassembler_long(tmp_dir):
#     """test plassembler long"""
#     longreads: Path = f"{end_to_end}/input_fastq.gz"
#     chromosome = 50000
#     outdir: Path = f"{end_to_end}/test_out"
#     cmd = f"plassembler long -l {longreads} -c {chromosome} -d {plassembler_db_dir} -o {outdir}  -t 8 -f"
#     exec_command(cmd)
#     remove_directory(outdir)


# def test_plassembler_assembled(tmp_dir):
#     """test plassembler assembled"""
#     longreads: Path = f"{end_to_end}/input_fastq.gz"
#     s1: Path = f"{end_to_end}/input_R1.fastq.gz"
#     s2: Path = f"{end_to_end}/input_R2.fastq.gz"
#     outdir: Path = f"{end_to_end}/test_out"
#     chromosome = 50000
#     input_plasmids = f"{end_to_end}/test_plasmids.fasta"
#     input_chromosome = f"{end_to_end}/test_chromosome.fasta"
#     cmd = f"plassembler assembled -l {longreads} -c {chromosome} -1 {s1} -2 {s2} -d {plassembler_db_dir} -o {outdir} --input_plasmids {input_plasmids} --input_chromosome {input_chromosome}  -t 8 -f"
#     exec_command(cmd)
#     remove_directory(outdir)


#################################
#     # on my mac the test works fine, but for CI Unicycler install is crook
#     # therefore just error out
#################################
# something to sort out with unicycler conda recipe maybe


class TestPlassemblerEndToEnd(unittest.TestCase):
    """Tests of plassembler end to end functions"""

    def test_plassembler(self):
        """test plassembler run"""
        with self.assertRaises(SystemExit):
            longreads: Path = f"{end_to_end}/input_fastq.gz"
            s1: Path = f"{end_to_end}/input_R1.fastq.gz"
            s2: Path = f"{end_to_end}/input_R2.fastq.gz"
            chromosome = 50000
            outdir: Path = f"{end_to_end}/test_out"
            cmd = f"plassembler run -l {longreads} -c {chromosome} -1 {s1} -2 {s2} -d {plassembler_db_dir} -o {outdir}  -t 8 -f"
            exec_command(cmd)
            remove_directory(outdir)

    def test_plassembler_long(self):
        """test plassembler long"""
        with self.assertRaises(SystemExit):
            longreads: Path = f"{end_to_end}/input_fastq.gz"
            chromosome = 50000
            outdir: Path = f"{end_to_end}/test_out"
            cmd = f"plassembler long -l {longreads} -c {chromosome} -d {plassembler_db_dir} -o {outdir}  -t 8 -f"
            exec_command(cmd)
            remove_directory(outdir)

    def test_plassembler_assembled(self):
        """test plassembler assembled"""
        with self.assertRaises(SystemExit):
            longreads: Path = f"{end_to_end}/input_fastq.gz"
            s1: Path = f"{end_to_end}/input_R1.fastq.gz"
            s2: Path = f"{end_to_end}/input_R2.fastq.gz"
            outdir: Path = f"{end_to_end}/test_out"
            chromosome = 50000
            input_plasmids = f"{end_to_end}/test_plasmids.fasta"
            input_chromosome = f"{end_to_end}/test_chromosome.fasta"
            cmd = f"plassembler assembled -l {longreads} -c {chromosome} -1 {s1} -2 {s2} -d {plassembler_db_dir} -o {outdir} --input_plasmids {input_plasmids} --input_chromosome {input_chromosome}  -t 8 -f"
            exec_command(cmd)
            remove_directory(outdir)
