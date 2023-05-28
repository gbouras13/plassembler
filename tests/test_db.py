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
import shutil


# import functions
from src.db import (get_database_zenodo, check_db_installation, instantiate_db_dir)

# data
test_data = Path("tests/test_data")
db_path = Path(f"{test_data}/Plassembler_Test_DB") 
val_data = Path(f"{test_data}/validation") 
tmp_db_path = Path(f"{test_data}/Plassembler_Test_DB_test") 

# make fake tempdir for testing
@pytest.fixture(scope="session")
def tmp_dir(tmpdir_factory):
    return tmpdir_factory.mktemp("tmp")

# to ensure sys exit on logger error
logger.add(lambda _: sys.exit(1), level="ERROR")


class test_install(unittest.TestCase):
    """Test for db"""

    def test_check_db_installation_good(self):
        expected_return = False
        check_db_installation(db_path)
        self.assertEqual(expected_return, False)

    def test_instantiate_db_good(self):
        expected_return = True
        instantiate_db_dir(db_path)
        self.assertEqual(expected_return, False)

    def test_check_db_installation_bad(self):
        with self.assertRaises(SystemExit):
            check_db_installation(val_data)

    # def test_get_database_zenodo(self):
    #     expected_return = False
    #     get_database_zenodo(tmp_db_path)
    #     # remove it after downloading
    #     shutil.rmtree(tmp_db_path)
    #     self.assertEqual(expected_return, False)

