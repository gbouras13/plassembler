"""
Unit tests for plassembler.

Usage: pytest

"""

# import
import unittest
import os
from pathlib import Path
import pytest
import shutil


# import functions
from src import input_commands
from src import concat
from src import depth
from src import assembly
from src.cleanup import remove_file
from src.qc import (chopper, fastp)
from src.mapping import (minimap_long_reads, minimap_short_reads)
from src.bam import (sam_to_bam_short)


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

class test_bam(unittest.TestCase):
    """Tests for bam.py"""
    # sam to bam
    def test_sam_to_bam_short(self):
        expected_return = True
        sam_to_bam_short(map_dir, threads, logdir)
        self.assertEqual(expected_return, True)





class test_mapping(unittest.TestCase):
    """Test for mapping"""
    # long read map
    def test_minimap_long_reads(self):
        expected_return = True
        pacbio_model = ""
        minimap_long_reads(map_dir, 1, pacbio_model, logdir)
        self.assertEqual(expected_return, True)

        # short read map
    def test_minimap_short_reads(self):
        expected_return = True
        minimap_short_reads(map_dir, 1, logdir)
        self.assertEqual(expected_return, True)




class test_qc_gzip(unittest.TestCase):
    """Test for qc"""
    # chopper
    def test_chopper_gzip(self):
        expected_return = True
        input_long_reads = os.path.join(test_data, "test_long.fastq.gz")
        chopper(input_long_reads, fake_out_dir, "500", "9", True, "1", logdir) # True for gunzip
        remove_file(os.path.join(fake_out_dir, "chopper_long_reads.fastq.gz"))
        self.assertEqual(expected_return, True)

    def test_chopper_not_gzip(self):
        expected_return = True
        input_long_reads = os.path.join(test_data, "test_long.fastq")
        chopper(input_long_reads, fake_out_dir, "500", "9", False, "1", logdir) # fasle for gunzip
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



# class test_assemblers(unittest.TestCase):
#     """Test for assembles"""

#     def test_flye(self):
#         expected_return = True
#         # C11 sim reads
#         assembly.run_flye(test_data, 8, raw_flag = False, pacbio_model = "nothing", logdir = logdir)
#         shutil.rmtree(os.path.join(test_data, "00-assembly"))
#         shutil.rmtree(os.path.join(test_data, "10-consensus"))
#         shutil.rmtree(os.path.join(test_data, "20-repeat"))
#         shutil.rmtree(os.path.join(test_data, "30-contigger"))
#         shutil.rmtree(os.path.join(test_data, "40-polishing"))
#         remove_file(os.path.join(test_data, "assembly.fasta"))
#         remove_file(os.path.join(test_data, "assembly_info.txt"))
#         remove_file(os.path.join(test_data, "assembly_graph.gfa"))
#         remove_file(os.path.join(test_data, "assembly_graph.gv"))
#         remove_file(os.path.join(test_data, "flye.log"))
#         self.assertEqual(expected_return, True)


    def test_raven(self):
        expected_return = True
        # C11 sim reads
        assembly.run_raven(test_data, 1,  logdir = logdir)
        remove_file(os.path.join(test_data, "assembly.fasta"))
        remove_file(os.path.join(test_data, "assembly_graph.gfa"))
        remove_file(os.path.join(test_data, "params.json"))
        remove_file("raven.cereal")
        self.assertEqual(expected_return, True)
        


 







    # def test_minimap_depth_sort_long(self):
    #     """test minimap depth sort long"""
    #     tmp = 1
    #     depth.minimap_depth_sort_long(out_dir, threads)
    #     assert tmp == 1








