#!/usr/bin/env python3
import os
import shutil
import tarfile
import urllib
from pathlib import Path

from loguru import logger


def check_db_installation(db_dir: Path, install_flag: bool):
    """checks database is installed correctly
    :param db_dir: database directory
    :param install_flag: whether to check or install database
    """
    # Mash files

    mash_db_names = ["plsdb.msh", "plsdb.tsv"]

    f1: Path = db_dir / f"{mash_db_names[0]}"
    f2: Path = db_dir / f"{mash_db_names[1]}"

    if f1.exists() and f2.exists():
        logger.info("PLSDB Database at {database} has already been downloaded")
    else:
        for file_name in mash_db_names:
            file_path: Path = db_dir / f"{file_name}"
            if file_path.exists() is False:
                if install_flag is True:
                    logger.info(
                        f"Database directory is missing file {file_path}. Plassembler Database will be downloaded."
                    )
                    get_database_zenodo(db_dir)
                    break
                else:
                    logger.error(
                        f"Database directory is missing {file_path}. Plassembler database needs to be downloaded using the plassembler download command."
                    )


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
        shutil.move(
            os.path.join(
                db_dir, "plsdb_110222_plassembler_v0.1.4_databases", "plsdb.msh"
            ),
            db_dir,
        )
        shutil.move(
            os.path.join(
                db_dir, "plsdb_110222_plassembler_v0.1.4_databases", "plsdb.tsv"
            ),
            db_dir,
        )
        shutil.rmtree(os.path.join(db_dir, "plsdb_110222_plassembler_v0.1.4_databases"))

        # remove tarball
        if os.path.exists(tar_path):
            os.remove(tar_path)
    except Exception:
        logger.error(
            "Plassembler Database Install Failed. \n Please try again or use the manual option detailed at https://github.com/gbouras13/plassembler.git \n to download the database from https://zenodo.org/record/7499200/files/plsdb_110222_plassembler_v0.1.4_databases.tar.gz"
        )
