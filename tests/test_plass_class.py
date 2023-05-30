"""
Unit tests for plassembler.

Usage: pytest

"""

# import
import unittest
from pathlib import Path
from loguru import logger



# import functions

from src.plassembler.utils.plass_class import Assembly, Plass
from src.plassembler.utils.cleanup import remove_file


# data
test_data = Path("tests/test_data")
plass_class_dir = Path(f"{test_data}/plass_class")
plass_class_depth_dir = Path(f"{plass_class_dir}/depth")
assembly_class_dir = Path(f"{test_data}/assembly_class")
assembly_depth_dir = Path(f"{assembly_class_dir}/depth")
plassembler_db_dir = Path(f"{test_data}/Plassembler_Test_DB")
logdir = Path(f"{test_data}/logs")


class test_plass_class(unittest.TestCase):
    """Tests for Plass class"""

    def test_get_contig_count(self):
        expected_return = True
        plass = Plass()
        plass.outdir = plass_class_dir
        plass.get_contig_count()
        # should be 2
        self.assertEqual(plass.contig_count, 2)

    def test_identify_chromosome_process_raven(self):
        plass = Plass()
        plass.outdir = plass_class_dir
        chrom_len = 10000000
        plass.identify_chromosome_process_raven(chrom_len)
        # should be False, too big chrom len
        self.assertEqual(plass.chromosome_flag, False)

    def test_identify_chromosome_process_flye(self):
        plass = Plass()
        plass.outdir = plass_class_dir
        chrom_len = 10000000
        plass.identify_chromosome_process_flye(chrom_len)
        # should be False, no chrom
        self.assertEqual(plass.chromosome_flag, False)

    def test_check_unicycler_success(self):
        plass = Plass()
        plass.check_unicycler_success(plass_class_dir)
        # should be True, as assembly.fasta exists
        self.assertEqual(plass.unicycler_success, True)

    def test_check_get_depth(self):
        expected = True
        plass = Plass()
        pacbio_model = "nothing"
        threads = 1
        # set to the depth dir for  intermediate files
        plass.outdir = plass_class_depth_dir
        plass.get_depth(logdir, pacbio_model, threads)
        remove_file(Path(f"{plass_class_depth_dir}/combined_long.sam"))
        remove_file(Path(f"{plass_class_depth_dir}/combined_short.sam"))
        remove_file(Path(f"{plass_class_depth_dir}/combined_sorted_long.bam"))
        remove_file(Path(f"{plass_class_depth_dir}/combined_sorted_short.bam"))
        remove_file(Path(f"{plass_class_depth_dir}/combined.fasta"))
        self.assertEqual(expected, True)

    def test_check_get_depth_long(self):
        expected = True
        plass = Plass()
        pacbio_model = "nothing"
        threads = 1
        # set to the depth dir for  intermediate files
        plass.outdir = plass_class_depth_dir
        plass.get_depth_long(logdir, pacbio_model, threads)
        remove_file(Path(f"{plass_class_depth_dir}/combined_long.sam"))
        remove_file(Path(f"{plass_class_depth_dir}/combined_sorted_long.bam"))
        self.assertEqual(expected, True)

    def test_process_mash_tsv(self):
        expected = True
        plass = Plass()
        # should be nothing due to fake db
        plass.outdir = plass_class_depth_dir
        plass.process_mash_tsv(plassembler_db_dir)
        self.assertEqual(expected, True)

    def test_combine_tsv(self):
        expected = True
        plass = Plass()
        plass.outdir = plass_class_depth_dir
        prefix = "plassembler"
        pacbio_model = "nothing"
        threads = 1
        plass.get_depth(logdir, pacbio_model, threads)
        plass.process_mash_tsv(plassembler_db_dir)
        plass.combine_depth_mash_tsvs(prefix)
        remove_file(Path(f"{plass_class_depth_dir}/combined_long.sam"))
        remove_file(Path(f"{plass_class_depth_dir}/combined_short.sam"))
        remove_file(Path(f"{plass_class_depth_dir}/combined_sorted_long.bam"))
        remove_file(Path(f"{plass_class_depth_dir}/combined_sorted_short.bam"))
        remove_file(Path(f"{plass_class_depth_dir}/{prefix}_summary.tsv"))
        self.assertEqual(expected, True)


class test_assembly_class(unittest.TestCase):
    """Tests for Assembly class"""

    def test_combine_input_fastas_good(self):
        expected_return = True
        assembly = Assembly()
        assembly.outdir = assembly_class_dir
        chrom_fasta: Path = Path(f"{assembly_class_dir}/chromosome.fasta")
        plasmid_fasta: Path = Path(f"{assembly_class_dir}/plasmid.fasta")
        assembly.combine_input_fastas(chrom_fasta, plasmid_fasta)
        self.assertEqual(expected_return, True)

    def test_check_get_depth(self):
        expected = True
        assembly = Assembly()
        pacbio_model = "nothing"
        threads = 1
        # set to the depth dir for  intermediate files
        assembly.outdir = assembly_depth_dir
        assembly.get_depth(logdir, pacbio_model, threads)
        remove_file(Path(f"{assembly_depth_dir}/combined_long.sam"))
        remove_file(Path(f"{assembly_depth_dir}/combined_short.sam"))
        remove_file(Path(f"{assembly_depth_dir}/combined_sorted_long.bam"))
        remove_file(Path(f"{assembly_depth_dir}/combined_sorted_short.bam"))
        self.assertEqual(expected, True)

    def test_process_mash_tsv(self):
        expected = True
        assembly = Assembly()
        plasmid_fasta = Path(f"{assembly_depth_dir}/plasmids.fasta")
        # should be nothing due to fake db
        assembly.outdir = assembly_depth_dir
        assembly.process_mash_tsv(plassembler_db_dir, plasmid_fasta)
        self.assertEqual(expected, True)

    def test_combine_tsv(self):
        expected = True
        assembly = Assembly()
        assembly.outdir = assembly_depth_dir
        prefix = "plassembler"
        pacbio_model = "nothing"
        threads = 1
        plasmid_fasta = Path(f"{assembly_depth_dir}/plasmids.fasta")
        assembly.get_depth(logdir, pacbio_model, threads)
        assembly.process_mash_tsv(plassembler_db_dir, plasmid_fasta)
        assembly.combine_depth_mash_tsvs(prefix)
        remove_file(Path(f"{assembly_depth_dir}/combined_long.sam"))
        remove_file(Path(f"{assembly_depth_dir}/combined_short.bam"))
        remove_file(Path(f"{assembly_depth_dir}/combined_sorted_long.bam"))
        remove_file(Path(f"{assembly_depth_dir}/combined_sorted_short.bam"))
        remove_file(Path(f"{assembly_depth_dir}/{prefix}_summary.tsv"))
        self.assertEqual(expected, True)
