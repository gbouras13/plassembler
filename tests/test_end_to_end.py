"""
Unit tests for plassembler.

Usage: pytest

"""

# import
import os
from pathlib import Path
import pytest
from loguru import logger
import sys
import shutil
import subprocess

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


#### 70kbp, 44kbp and 9kbp plasmid reads are from
### the 70kbp is a fake chromosome

def test_plassembler(tmp_dir):
    """test plassembler run"""
    longreads: Path = f"{end_to_end}/input_fastq.gz"
    s1: Path = f"{end_to_end}/input_R1.fastq.gz"
    s2: Path = f"{end_to_end}/input_R2.fastq.gz"
    chromosome = 50000
    outdir: Path = f"{end_to_end}/test_out"
    cmd = f"plassembler run -l {longreads} -c {chromosome} -1 {s1} -2 {s2} -d {plassembler_db_dir} -o {outdir}  -t 8 -f"
    exec_command(cmd)
    remove_directory(outdir)

def test_plassembler_long(tmp_dir):
    """test plassembler long"""
    longreads: Path = f"{end_to_end}/input_fastq.gz"
    chromosome = 50000
    outdir: Path = f"{end_to_end}/test_out"
    cmd = f"plassembler long -l {longreads} -c {chromosome} -d {plassembler_db_dir} -o {outdir}  -t 8 -f"
    exec_command(cmd)
    remove_directory(outdir)

def test_plassembler_assembled(tmp_dir):
    """test plassembler assembled"""
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