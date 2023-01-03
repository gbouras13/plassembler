#!/usr/bin/env python3
import argparse
import os
from argparse import RawTextHelpFormatter
import subprocess as sp
import sys


# define the 2 output files in the mash database
MASH_DB_NAMES = ['plsdb.msh',
'plsdb.tsv']

def get_db_input():
    """gets input for install_database.py
    :return: args
    """
    parser = argparse.ArgumentParser(description='script to download required PLSDB databases', formatter_class=RawTextHelpFormatter)
    parser.add_argument('-d', '--database', action="store", required = True, help='Database Directory.')
    args = parser.parse_args()
    return args


def instantiate_db_dir(db_dir):
    """ instatiate database dir
	:param db_dir: database directory
    :return: 
    """
    if os.path.isdir(db_dir) == False:
        os.mkdir(db_dir)
 
def check_db_installation(db_dir):
    """ checks database is installed correctly
	:param db_dir: database directory
    :return: downloaded_flag: boolean whether database is correctly downloaded.
    """
    downloaded_flag = True
    # Mash files
    for file_name in MASH_DB_NAMES:
        path = os.path.join(db_dir, file_name)
        if os.path.isfile(path) == False:
            print("Databases are missing. Plassembler Database Needs to be Downloaded.")
            downloaded_flag = False
            break
    return downloaded_flag
    

def get_database_zenodo(db_dir):
    print("Downloading Plassembler Database.")
    tarball = 'plsdb_110222_plassembler_v0.1.4_databases.tar.gz'
    url = "https://zenodo.org/record/7499200/files/plsdb_110222_plassembler_v0.1.4_databases.tar.gz"
    try:
        # remvoe the directory
        sp.call(["rm", "-rf", os.path.join(db_dir)])
        # make db dir
        sp.call(["mkdir", "-p", os.path.join(db_dir)])
        # download the tarball
        sp.call(["curl", url, "-o", os.path.join(db_dir,tarball)])
        # untar tarball into database directory
        sp.call(["tar", "-xzf", os.path.join(db_dir, tarball), "-C", db_dir, "--strip-components=1"])
        # remove tarball
        sp.call(["rm","-f", os.path.join(db_dir,tarball)])
    except:
        sys.stderr.write("Error: Plassembler Database Install Failed. \n Please try again or use the manual option detailed at https://github.com/gbouras13/plassembler.git \n to download the database from https://zenodo.org/record/7499200/files/plsdb_110222_plassembler_v0.1.4_databases.tar.gzz")  
        return 0


if __name__ == "__main__":
    args = get_db_input()
    instantiate_db_dir(args.database)
    downloaded_flag = check_db_installation(args.database)
    if downloaded_flag == True:
        print("PLSDB Database has already been Downloaded and Checked.")
    else:
        get_database_zenodo(args.database)

