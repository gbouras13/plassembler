import argparse
import os
import sys
import gzip
from argparse import RawTextHelpFormatter
from Bio import SeqIO
import shutil
from version import __version__

v = __version__

### GLOBAL VARIABLES

def get_input():
	parser = argparse.ArgumentParser(description='plassembler: accurate extra-chromosomal plasmid assembler pipeline for haploid bacterial genomes.', formatter_class=RawTextHelpFormatter)
	parser.add_argument('-l', '--longreads', action="store", help='Fastq File of ONT Long Reads. Required',  required=True)
	parser.add_argument('-o', '--outdir', action="store", help='Directory to write the output to.', default=os.path.join(os.getcwd(), "output/") )
	parser.add_argument('-s1', '--short_one', action="store", help='R1 short read fastq file. Required.',  required=True)
	parser.add_argument('-s2', '--short_two', action="store", help='R2 short read fastq file. Required.',  required=True)
	parser.add_argument('-m', '--min_length', action="store", help='minimum length for long reads for nanofilt. Defaults to 1000.',  default='1000')
	parser.add_argument('-t', '--threads', help="Number of threads for flye and unicycler. Defaults to 8.", action="store", default = str(8))
	parser.add_argument('-f', '--force', help="Overwrites the output directory.", action="store_true" )
	parser.add_argument('-p', '--prefix', action="store", help='Prefix for output files. This is not required',  default='Default')
	parser.add_argument('-c', '--chromosome', action="store", help='Approximate chromosome length of bacteria',  default=2500000)
	parser.add_argument('-q', '--min_quality', action="store", help='minimum quality of long reads for nanofilt. Defaults to 9.',  default=str(9))
	parser.add_argument('-V', '--version', action='version', version=v)
	args = parser.parse_args()

	return args

def instantiate_dirs(output_dir, force):
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
	# to get extension
	filename, file_extension = os.path.splitext(file)
	# flag for whether file is zipped
	zipped = True
	if file_extension == ".gz":
	# if gzipped 
		with gzip.open(file, "rt") as handle:
			fastq = SeqIO.parse(handle, "fastq")
			if any(fastq):
				print("FASTQ " +file + " checked")
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


