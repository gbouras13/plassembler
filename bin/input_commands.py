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
	parser = argparse.ArgumentParser(description='plassembler: accurate plasmid assembler pipeline for bacterial genomes.', formatter_class=RawTextHelpFormatter)
	parser.add_argument('-l', '--longreads', action="store", help='Fastq File of Illumina Long Reads.',  required=True)
	parser.add_argument('-o', '--outdir', action="store", help='Directory to write the output to.', default=os.path.join(os.getcwd(), "output/") )
	parser.add_argument('-s1', '--short_one', action="store", help='R1 short read fastq file.',  required=True)
	parser.add_argument('-s2', '--short_two', action="store", help='R2 short read fastq file.',  required=True)
	parser.add_argument('-t', '--threads', help="Number of threads for flye and unicycler. Defaults to 8.", action="store", default = str(8))
	parser.add_argument('-f', '--force', help="Overwrites the output directory.", action="store_true" )
	parser.add_argument('-p', '--prefix', action="store", help='Prefix for output files. This is not required',  default='Default')
	parser.add_argument('-c', '--chromosome', action="store", help='Approximate chromosome length of bacteria',  default=2500000)
	parser.add_argument('-V', '--version', action='version', version=v)
	args = parser.parse_args()

	return args

def instantiate_dirs(output_dir, force):
	# remove outdir on force
	if force == True:
		if os.path.isdir(output_dir) == True:
			#shutil.rmtree(output_dir)
			print('placeholder')
		else:
			print("\n--force was specified even though the outdir does not already exist. Continuing \n")
	else:
		if os.path.isdir(output_dir) == True:
			sys.exit("\nOutput directory already exists and force was not specified. Please specify -f or --force to overwrite the output directory. \n")  
	# instantiate outdir
	if os.path.isdir(output_dir) == False:
		os.mkdir(output_dir)
	return output_dir

def validate_fastq(filename):
	with gzip.open(filename, "rt") as handle:
		fastq = SeqIO.parse(handle, "fastq")
		if any(fastq):
			print("FASTQ " +filename + " checked")
		else:
			sys.exit("Error: Input file is not in the FASTQ format.\n")  


