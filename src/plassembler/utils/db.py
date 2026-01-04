#!/usr/bin/env python3
"""
taken from pharokka and therefore from bakta
"""

import hashlib
import os
import shutil
import tarfile
from pathlib import Path

import requests
from alive_progress import alive_bar
from loguru import logger

from plassembler.utils.cleanup import remove_directory
from plassembler.utils.util import get_version


def check_db_installation(db_dir: Path, force: bool, install_flag: bool):
    """checks database is installed correctly
    :param db_dir: database directory
    :param force: boolean if force is true
    :param install_flag: whether to check or install database
    """
    # Mash files

    mash_db_names = ["plsdb_2023_11_03_v2.msh", "plsdb_2023_11_03_v2.tsv"]

    f1: Path = db_dir / f"{mash_db_names[0]}"
    f2: Path = db_dir / f"{mash_db_names[1]}"

    if install_flag is True:
        if force is True:
            if os.path.isdir(db_dir) is True:
                logger.info(f"Removing the directory {db_dir} as --force was specified")
                shutil.rmtree(db_dir)
            else:
                logger.info(
                    f"--force was specified even though the directory {db_dir} does not already exist. Continuing"
                )
        else:
            if os.path.isdir(db_dir) is True:
                logger.error(
                    f"Directory {db_dir} already exists and force was not specified. Please specify -f or --force to overwrite {db_dir}"
                )

    # instantiate outdir
    if os.path.isdir(db_dir) is False:
        os.mkdir(db_dir)

    if f1.exists() and f2.exists():
        logger.info(f"PLSDB Database mash sketch at {f1} exists.")
        logger.info(f"PLSDB Database tsv metadata file at {f2} exists.")
        logger.info(f"PLSDB Database at {db_dir} has already been downloaded")
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
    logger.info(f"Downloading Plassembler Database to the directory {db_dir}")
    tarball = "201123_plassembler_v1.5.0_databases.tar.gz"
    tar_path = Path(f"{db_dir}/{tarball}")
    db_url = "https://zenodo.org/record/10158040/files/201123_plassembler_v1.5.0_databases.tar.gz"
    requiredmd5 = "3a24bacc05bb857dc044fc6662b58db7"

    version = get_version().strip()

    headers = {
        "User-Agent": f"plassembler/{version} (contact: george.bouras@adelaide.edu.au)"
    }

    try:
        with tar_path.open("wb") as fh_out, requests.get(
            db_url, stream=True, headers=headers
        ) as resp:
            total_length = resp.headers.get("content-length")
            if total_length is not None:  # content length header is set
                total_length = int(total_length)
            with alive_bar(total=total_length, scale="SI") as bar:
                for data in resp.iter_content(chunk_size=1024 * 1024):
                    fh_out.write(data)
                    bar(count=len(data))
    except IOError:
        logger.error(
            f"ERROR: Could not download file from Zenodo! url={db_url}, path={tar_path}"
        )

    md5_sum = calc_md5_sum(tar_path)

    if md5_sum == requiredmd5:
        logger.info(f"Database file download OK: {md5_sum}")
    else:
        logger.error(
            f"Error: corrupt database file! MD5 should be '{requiredmd5}' but is '{md5_sum}'"
        )

    logger.info(f"Extracting DB tarball: file={tar_path}, output={db_dir}")
    untar(tar_path, db_dir)
    tar_path.unlink()
    logger.info(f"Plassembler Database download into {db_dir} successful.")


def calc_md5_sum(tarball_path: Path, buffer_size: int = 1024 * 1024) -> str:
    """
    gets md5 of a file
    """
    md5 = hashlib.md5()
    with tarball_path.open("rb") as fh:
        data = fh.read(buffer_size)
        while data:
            md5.update(data)
            data = fh.read(buffer_size)
    return md5.hexdigest()


def untar(tarball_path: Path, output_path: Path):
    """
    untars a file
    """
    try:
        with tarball_path.open("rb") as fh_in, tarfile.open(
            fileobj=fh_in, mode="r:gz"
        ) as tar_file:
            tar_file.extractall(path=str(output_path))

        # get untarred directory
        untarpath = os.path.join(output_path, "201123_plassembler_v1.5.0_databases")

        # Get a list of all files in the source directory
        files_to_move = [
            f
            for f in os.listdir(untarpath)
            if os.path.isfile(os.path.join(untarpath, f))
        ]

        # Move each file to the destination directory
        for file_name in files_to_move:
            source_path = os.path.join(untarpath, file_name)
            destination_path = os.path.join(output_path, file_name)
            shutil.move(source_path, destination_path)
        # remove the directort
        remove_directory(untarpath)

    except OSError:
        logger.error(f"Could not extract {tarball_path} to {output_path}")
