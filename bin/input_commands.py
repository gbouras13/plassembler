import argparse
import os
import sys
import gzip
from argparse import RawTextHelpFormatter
from Bio import SeqIO
import shutil
import subprocess as sp
from version import __version__

v = __version__

### GLOBAL VARIABLES

def get_input():
	"""gets input for plassembler
    :return: args
    """
	parser = argparse.ArgumentParser(description='plassembler: accurate extra-chromosomal plasmid assembler pipeline for haploid bacterial genomes.', formatter_class=RawTextHelpFormatter)
	parser.add_argument('-d', '--database', action="store", help='Directory of PLSDB database downloaded using install_database.py.',  required=True)
	parser.add_argument('-l', '--longreads', action="store", help='Fastq File of ONT Long Reads. Required',  required=True)
	parser.add_argument('-1', '--short_one', action="store", help='R1 short read fastq file. Required.',  required=True)
	parser.add_argument('-2', '--short_two', action="store", help='R2 short read fastq file. Required.',  required=True)
	parser.add_argument('-c', '--chromosome', action="store", help='Approximate chromosome length of bacteria. Defaults to 2500000.',  default=2500000)
	parser.add_argument('-o', '--outdir', action="store", help='Directory to write the output to. Defaults to output/', default=os.path.join(os.getcwd(), "output/") )
	parser.add_argument('-m', '--min_length', action="store", help='minimum length for long reads for nanofilt. Defaults to 500.',  default='500')
	parser.add_argument('-t', '--threads', help="Number of threads for flye and unicycler. Defaults to 1.", action="store", default = str(1))
	parser.add_argument('-f', '--force', help="Overwrites the output directory.", action="store_true" )
	parser.add_argument('-r', '--raw_flag', help="Use --nano-raw for Flye Guppy FAST reads. \nBy default, Flye will assume SUP or HAC reads and use --nano-hq", action="store_true" )
	parser.add_argument('-p', '--prefix', action="store", help='Prefix for output files. This is not required',  default='Default')
	parser.add_argument('-q', '--min_quality', action="store", help='minimum quality of long reads for nanofilt. Defaults to 9.',  default=str(9))
	parser.add_argument('-V', '--version', action='version',help='show plassembler version and exit.', version=v)
	args = parser.parse_args()

	return args

def instantiate_dirs(output_dir, force):
	"""checks that the output directory doesn't already exist, overwrites if forced
	:param output_dir: output directory path
    :param force: flag to overwrite the output directory
    :return: output_dir
    """
	# remove outdir on force
	if force == True:
		if os.path.isdir(output_dir) == True:
			shutil.rmtree(output_dir)
		else:
			print("\n--force was specified even though the outdir does not already exist. Continuing \n")
	else:
		if os.path.isdir(output_dir) == True:
			sys.exit("\nOutput directory already exists and force was not specified. Please specify -f or --force to overwrite the output directory. \n")  
	# instantiate outdir
	if os.path.isdir(output_dir) == False:
		os.mkdir(output_dir)
	return output_dir

def validate_fastq(file):
	"""Checks the input fastq is really a fastq
	:param file: fastq file
    :return: zipped - Boolean whether the input fastq is gzipped.
    """
	# to get extension
	filename, file_extension = os.path.splitext(file)
	# flag for whether file is zipped
	zipped = True
	if file_extension == ".gz":
	# if gzipped 
		with gzip.open(file, "rt") as handle:
			fastq = SeqIO.parse(handle, "fastq")
			if any(fastq):
				print("FASTQ " + file + " checked")
			else:
				sys.exit("Error: Input file is not in the FASTQ format.\n")  
	else:
		zipped = False
		with open(file, "r") as handle:
			fastq = SeqIO.parse(handle, "fastq")
			if any(fastq):
				print("FASTQ " +file + " checked")
			else:
				sys.exit("Error: Input file is not in the FASTQ format.\n") 
	return zipped

def check_dependencies(logger):
	"""Checks the version of Unicycler, spades and Flye
    :return:
    """
	# Flye
	process = sp.Popen(["flye", "--version"], stdout=sp.PIPE, stderr=sp.STDOUT) 
	flye_out, _ = process.communicate()
	flye_out = flye_out.decode().strip()
	flye_major_version = int(flye_out.split('.')[0])
	flye_minor_version = int(flye_out.split('.')[1])
	flye_minorest_version = flye_out.split('.')[2]

	print("Flye version found is v" + str(flye_major_version) +"." + str(flye_minor_version) +"."+flye_minorest_version + ".")
	logger.info("Flye version found is v" + str(flye_major_version) +"." + str(flye_minor_version) +"."+flye_minorest_version +".")

	if flye_major_version != 2:
		sys.exit("Flye is too old - please reinstall plassembler.")
	if flye_minor_version != 9:
		sys.exit("Flye is too old - please reinstall plassembler.")

	print("Flye version is ok.")
	logger.info("Flye version is ok.")

	# unicycler
	process = sp.Popen(["unicycler", "--version"], stdout=sp.PIPE, stderr=sp.STDOUT) 
	unicycler_out, _ = process.communicate()
	unicycler_out = unicycler_out.decode()
	unicycler_version = unicycler_out.split(' ')[1]
	# get rid of the "v"
	unicycler_version = unicycler_version[1:]

	unicycler_major_version = int(unicycler_version.split('.')[0])
	unicycler_minor_version = int(unicycler_version.split('.')[1])
	unicycler_minorest_version = int(unicycler_version.split('.')[2])

	print("Unicycler version found is v" + str(unicycler_major_version) +"." + str(unicycler_minor_version) +"."+str(unicycler_minorest_version)+".")
	logger.info("Unicycler version found is v" + str(unicycler_major_version) +"." + str(unicycler_minor_version) +"."+str(unicycler_minorest_version)+".")

	if unicycler_minor_version < 4 :
		sys.exit("Unicycler is too old - please reinstall plassembler, see instructions at https://github.com/gbouras13/plassembler.")
	elif unicycler_minor_version == 4 and unicycler_minorest_version < 8:
		sys.exit("Unicycler is too old - please reinstall plassembler, see instructions at https://github.com/gbouras13/plassembler.")
	elif unicycler_minor_version == 4 and unicycler_minorest_version >= 8:
		print("Unicycler version is older than v0.5.0 - plassembler will continue but please consider installing Unicycler v0.5.0. See instructions as https://github.com/gbouras13/plassembler.")
		logger.info("Unicycler version is older than v0.5.0 - plassembler will continue but please consider installing Unicycler v0.5.0. See instructions as https://github.com/gbouras13/plassembler.")
	else:
		print("Unicycler version is ok.")
		logger.info("Unicycler version is ok.")










	# # to get extension
	# filename, file_extension = os.path.splitext(file)
	# # flag for whether file is zipped
	# zipped = True
	# if file_extension == ".gz":
	# # if gzipped 
	# 	with gzip.open(file, "rt") as handle:
	# 		fastq = SeqIO.parse(handle, "fastq")
	# 		if any(fastq):
	# 			print("FASTQ " + file + " checked")
	# 		else:
	# 			sys.exit("Error: Input file is not in the FASTQ format.\n")  
	# else:
	# 	zipped = False
	# 	with open(file, "r") as handle:
	# 		fastq = SeqIO.parse(handle, "fastq")
	# 		if any(fastq):
	# 			print("FASTQ " +file + " checked")
	# 		else:
	# 			sys.exit("Error: Input file is not in the FASTQ format.\n") 
	# return zipped


