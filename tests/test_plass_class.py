"""
Unit tests for plassembler.

Usage: pytest

"""

# import
import unittest
import os
from pathlib import Path
import pytest
from loguru import logger
import sys
import subprocess as sp


# import functions

from src.plass_class import Assembly, Plass


# data
test_data = Path("tests/test_data")
plass_class_dir = Path(f"{test_data}/plass_class") 
assembly_class_dir = Path(f"{test_data}/assembly_class") 


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


class test_assembly_class(unittest.TestCase):
    """Tests for Assembly class"""

    def test_combine_input_fastas_good(self):
        expected_return = True
        assembly = Assembly()
        assembly.outdir = assembly_class_dir
        chrom_fasta: Path = Path(f"{assembly_class_dir}/chromosome.fasta") 
        plasmid_fasta: Path =  Path(f"{assembly_class_dir}/plasmid.fasta") 
        assembly.combine_input_fastas(chrom_fasta, plasmid_fasta)
        self.assertEqual(expected_return, True)
