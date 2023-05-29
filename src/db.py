#!/usr/bin/env python3
import os
import subprocess as sp
from pathlib import Path
from loguru import logger
import shutil
import urllib
import hashlib
import tarfile


# define the 2 output files in the mash database
MASH_DB_NAMES = ["plsdb.msh", "plsdb.tsv"]


def instantiate_db_dir(db_dir: Path):
    """instatiate database dir
        :param db_dir: database directory
    :return:
    """
    if os.path.isdir(db_dir) == False:
        os.mkdir(db_dir)


def check_db_installation(db_dir: Path, install_flag: bool):
    """checks database is installed correctly
    :param db_dir: database directory
    :return: install_flag: whether this is being called to download or just check the db
    """
    downloaded_flag = True
    # Mash files
    for file_name in MASH_DB_NAMES:
        file_path: Path  = db_dir/f"{file_name}"
        if file_path.exists() == False:
            if install_flag == True:
                logger.info(f"Database directory is missing {file_path}. Plassembler database will be downloaded.")
                downloaded_flag = False
            else:
                logger.error(f"Database directory is missing {file_path}. Plassembler database needs to be downloaded using the plassembler download command.")
    return downloaded_flag



def get_database_zenodo(db_dir: Path):
    logger.info("Downloading Plassembler Database.")
    tarball = "plsdb_110222_plassembler_v0.1.4_databases.tar.gz"
    tar_path = os.path.join(db_dir, tarball)
    url = "https://zenodo.org/record/7499200/files/plsdb_110222_plassembler_v0.1.4_databases.tar.gz"
    try:
        # remvoe the directory
        if os.path.exists(db_dir):
            shutil.rmtree(db_dir)

        # make db dir
        if not os.path.exists(db_dir):
            os.mkdir(db_dir)
        # download the tarball
        urllib.request.urlretrieve(url=url, filename=tar_path)

        # untar tarball into database directory
        with tarfile.open(tar_path, "r:gz") as tar:
            tar.extractall(path=db_dir)

        # move files around
        shutil.move(os.path.join(db_dir, "plsdb_110222_plassembler_v0.1.4_databases", "plsdb.msh"), db_dir)
        shutil.move(os.path.join(db_dir, "plsdb_110222_plassembler_v0.1.4_databases", "plsdb.tsv"), db_dir)
        shutil.rmtree(os.path.join(db_dir, "plsdb_110222_plassembler_v0.1.4_databases"))

        # remove tarball
        if os.path.exists(tar_path):
            os.remove(tar_path)
    except:
        logger.error(
            "Plassembler Database Install Failed. \n Please try again or use the manual option detailed at https://github.com/gbouras13/plassembler.git \n to download the database from https://zenodo.org/record/7499200/files/plsdb_110222_plassembler_v0.1.4_databases.tar.gz"
        )




