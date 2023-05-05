import argparse
import os
import sys
import gzip
from argparse import RawTextHelpFormatter
from Bio import SeqIO
import shutil
import subprocess as sp
from version import __version__
import log


v = __version__

### GLOBAL VARIABLES

def get_input():
	"""gets input for plassembler
    :return: args
    """
	parser = argparse.ArgumentParser(description='plassembler: automated bacterial plasmid assembly tool.', formatter_class=RawTextHelpFormatter)
	parser.add_argument('-d', '--database', action="store", help='Directory of PLSDB database downloaded using install_database.py.',  required=True)
	parser.add_argument('-l', '--longreads', action="store", help='Fastq file of long reads.',  default='nothing')
	parser.add_argument('-1', '--short_one', action="store", help='R1 short read fastq file.',  default='nothing')
	parser.add_argument('-2', '--short_two', action="store", help='R2 short read fastq file.',  default='nothing')
	parser.add_argument('-c', '--chromosome', action="store", help='Approximate lower-bound chromosome length of bacteria. \nDefaults to 2500000.',  default=2500000)
	parser.add_argument('-o', '--outdir', action="store", help='Directory to write the output to. Defaults to output/', default=os.path.join(os.getcwd(), "output/") )
	parser.add_argument('-m', '--min_length', action="store", help='minimum length for filtering long reads with chopper. Defaults to 500.',  default='500')
	parser.add_argument('-q', '--min_quality', action="store", help='minimum quality for filtering long reads with chopper. Defaults to 9.',  default=str(9))
	parser.add_argument('-t', '--threads', help="Number of threads. Defaults to 1.", action="store", default = str(1))
	parser.add_argument('-f', '--force', help="Overwrites the output directory.", action="store_true" )
	parser.add_argument('-r', '--raw_flag', help="Use --nano-raw for Flye. \n Designed for Guppy fast configuration reads. \nBy default, Flye will assume SUP or HAC reads and use --nano-hq", action="store_true" )
	parser.add_argument('-p', '--prefix', action="store", help='Prefix for output files. This is not required',  default='Default')
	parser.add_argument('-s', '--subsample_depth',  help='Subsample long-read depth as an integer. \nUsed combined with the coverage of the chromosome length provided with -c. \nDefaults to 30.',  default=str(30) )
	parser.add_argument('-k', '--kmer_mode',  help='Very high quality Nanopore R10.4 and above reads. \nNo short reads required. Experimental for now.', action="store_true" )
	parser.add_argument('--pacbio_model',  help='Pacbio Flye model. Must be pacbio-raw, pacbio-corr or pacbio-hifi. \nUse pacbio-raw for PacBio regular CLR reads (<20 percent error), \npacbio-corr for PacBio reads that were corrected with other methods (<3 percent error) \nor pacbio-hifi for PacBio HiFi reads (<1 percent error).', action="store", default='nothing')
	parser.add_argument('--no_subsample',  help='Turns off long-read sub-sampling. \nRecommended if long-read sets have low N50s/N90s, \nor are of a difficult-to-assemble species with lots of repeats.', action="store_true" )
	parser.add_argument('--keep_fastqs',  help='Whether you want to keep fastq files containing putative plasmid reads.', action="store_true")
	parser.add_argument('--keep_chromosome',  help='Whether you want to keep the unpolished Flye chromosome assembly.', action="store_true")
	parser.add_argument('-a', '--assembled_mode',  help='Activates assembled mode..', action="store_true")
	parser.add_argument('--input_chromosome',  help='Input FASTA file consisting of already assembled chromosome with assembled mode. \nMust be 1 complete contig.', action="store", default='nothing')
	parser.add_argument('--input_plasmids',  help='Input FASTA file consisting of already assembled plasmids with assembled mode. \nRequires FASTQ file input (short only, long only or long + short).', action="store", default='nothing')
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
				sys.exit("Error: Input " + file + "is not in the FASTQ format.\n")  
	else:
		zipped = False
		with open(file, "r") as handle:
			fastq = SeqIO.parse(handle, "fastq")
			if any(fastq):
				print("FASTQ " +file + " checked")
			else:
				sys.exit("Error: Input " + file + "is not in the FASTQ format.\n") 
	return zipped

def validate_fasta(filename):
	"""Checks the input insta is really a fasta
	:param file: fasta file
    :return: 
    """
	with open(filename, "r") as handle:
		fasta = SeqIO.parse(handle, "fasta")
		if any(fasta):
			print("FASTA checked")
		else:
			sys.exit("Error: Input file is not in the FASTA format.\n")  


def validate_fastas_assembled_mode(input_chromosome, input_plasmids):
	"""Checks the input insta is really a fasta
	:param file: fasta file
    :return: 
    """
    # chromosome
	validate_fasta(input_chromosome)

	with open(input_chromosome, "r") as fasta:
		# count contigs
		records = list(SeqIO.parse(fasta, "fasta"))
		num_contigs = len(records)
		if num_contigs > 1:
			sys.exit("Error: There are multiple contigs in your chromosome FASTA. Please input a completed chromosome.\n") 

	# plasmids
	validate_fasta(input_plasmids)



def validate_fastqs_assembled_mode(longreads, short_one, short_two):
	"""Checks the input instq are really fastqs
	:param longreads: long read file
	:param short_one: short_one read file
	:param short_two: short_two read file
    :return: 
    """

    # long
	long_flag = False
	long_gzipped = False
	if longreads != "nothing":
		print("You have input long read FASTQs for depth calculation.")
		long_gzipped = validate_fastq(longreads)
		long_flag = True
	
	# short
	short_flag = False
	if short_one != "nothing" and short_two != "nothing":
		print("You have input paired short read FASTQs for depth calculation.")
		s1_gzipped = validate_fastq(short_one)
		s2_gzipped = validate_fastq(short_two)
		if s1_gzipped != s2_gzipped:
			sys.exit("R1 and R2 files are inconsistenly compressed. Please check the compression format and try again.")
		short_flag = True
	
	if short_flag == False and long_flag == False:
		sys.exit("No valid long read or paired short read FASTQs were input. Please check your input and try again.")

	if (short_one != "nothing" and short_two == "nothing") or (short_one == "nothing" and short_two != "nothing") :
		sys.exit("Only 1 short read file was found. Please check your input and try again.")

	return (short_flag, long_flag, long_gzipped)




def check_dependencies(logger):
	"""Checks the version of Unicycler, spades and Flye
    :return:
    """
	# Flye
	try:
		process = sp.Popen(["flye", "--version"], stdout=sp.PIPE, stderr=sp.STDOUT) 
		flye_out, _ = process.communicate()
		flye_out = flye_out.decode().strip()
		flye_major_version = int(flye_out.split('.')[0])
		flye_minor_version = int(flye_out.split('.')[1])
		flye_minorest_version = flye_out.split('.')[2]
	except:
		sys.exit("Flye not found. Please reinstall Plassembler.")

	message = "Flye version found is v" + str(flye_major_version) +"." + str(flye_minor_version) +"."+flye_minorest_version + "."
	log.write_message(message, logger)

	if flye_major_version != 2:
		sys.exit("Flye is too old - please reinstall plassembler, see instructions at https://github.com/gbouras13/plassembler.")
	if flye_minor_version < 9:
		sys.exit("Flye is too old - please reinstall plassembler, see instructions at https://github.com/gbouras13/plassembler.")

	message = "Flye version is ok."
	log.write_message(message, logger)

	# unicycler
	try:
		process = sp.Popen(["unicycler", "--version"], stdout=sp.PIPE, stderr=sp.STDOUT) 
		unicycler_out, _ = process.communicate()
		unicycler_out = unicycler_out.decode()
		unicycler_version = unicycler_out.split(' ')[1]
		# get rid of the "v"
		unicycler_version = unicycler_version[1:]

		unicycler_major_version = int(unicycler_version.split('.')[0])
		unicycler_minor_version = int(unicycler_version.split('.')[1])
		unicycler_minorest_version = int(unicycler_version.split('.')[2])
	except:
		sys.exit("Unicycler not found. Please reinstall Plassembler, see instructions at https://github.com/gbouras13/plassembler.")

	message = "Unicycler version found is v" + str(unicycler_major_version) +"." + str(unicycler_minor_version) +"."+str(unicycler_minorest_version)+"."
	log.write_message(message, logger)

	if unicycler_minor_version < 4 :
		sys.exit("Unicycler is too old - please reinstall Plassembler, see instructions at https://github.com/gbouras13/plassembler.")
	elif unicycler_minor_version == 4 and unicycler_minorest_version < 8:
		sys.exit("Unicycler is too old - please reinstall Plassembler, see instructions at https://github.com/gbouras13/plassembler.")
	elif unicycler_minor_version == 4 and unicycler_minorest_version >= 8:
		message = "Unicycler version is older than v0.5.0 - Plassembler will continue but please consider installing Unicycler v0.5.0. See instructions at https://github.com/gbouras13/plassembler."
		log.write_message(message, logger)
	else:
		message = "Unicycler version is ok."
		log.write_message(message, logger)

#spades
	try:
		process = sp.Popen(["spades.py", "--version"], stdout=sp.PIPE, stderr=sp.PIPE)  
		spades_out, _ = process.communicate()
		spades_out = spades_out.decode()
		spades_version = spades_out.split(' ')[3]
		spades_version = spades_version.split("\n")[0]
		message ="SPAdes " + str(spades_version) + " found."
		log.write_message(message, logger)
	except:
		sys.exit("SPAdes not found.\n")  

#samtools
	try:
		process = sp.Popen(["samtools", "--version"], stdout=sp.PIPE, stderr=sp.PIPE) 
		samtools_out, _ = process.communicate()
		samtools_out = samtools_out.decode()
		samtools_version = samtools_out.split("\n")[0].split(' ')[1] # get second line, and then second component of line
		message ="Samtools v" + str(samtools_version) + " found."
		log.write_message(message, logger)
	except:
		sys.exit("Samtools not found.\n")  

#minimap2
	try:
		process = sp.Popen(["minimap2", "--version"], stdout=sp.PIPE, stderr=sp.PIPE) 
		minimap2_out, _ = process.communicate()
		minimap2_version = minimap2_out.decode()
		minimap2_version = minimap2_version.split("\n")[0]
		message ="minimap2 v" + str(minimap2_version) + " found."
		log.write_message(message, logger)
	except:
		sys.exit("minimap2 not found.\n")  

#fastp
	try:
		process = sp.Popen(["fastp", "--version"], stdout=sp.PIPE, stderr=sp.PIPE) 
		_, fastp_out = process.communicate()
		fastp_version = fastp_out.decode()
		fastp_version = fastp_version.split("\n")[0].split(' ')[1]
		message ="fastp v" + str(fastp_version) + " found."
		log.write_message(message, logger)
	except:
		sys.exit("fastp not found.\n")  

#chopper
	try:
		process = sp.Popen(["chopper", "--version"], stdout=sp.PIPE, stderr=sp.PIPE) 
		chopper_out, _ = process.communicate()
		chopper_version = chopper_out.decode()
		chopper_version = chopper_version.split("\n")[0].split(' ')[1]
		message ="chopper v" + str(chopper_version) + " found."
		log.write_message(message, logger)
	except:
		sys.exit("chopper not found.\n")  

#seqkit
	try:
		process = sp.Popen(["seqkit", "version"], stdout=sp.PIPE, stderr=sp.PIPE) 
		seqkit_out, _ = process.communicate()
		seqkit_version = seqkit_out.decode()
		seqkit_version = seqkit_version.split("\n")[0].split(' ')[1]
		message ="seqkit " + str(seqkit_version) + " found."
		log.write_message(message, logger)
	except:
		sys.exit("seqkit not found.\n")  

#mash
	try:
		process = sp.Popen(["mash", "version"], stdout=sp.PIPE, stderr=sp.PIPE) 
		mash_out, _ = process.communicate()
		mash_out = mash_out.decode()
		version_line = []
		for line in mash_out.split("\n"):
			if "version" in line:
				version_line.append(line)
		mash_version = version_line[0].split(' ')[2]
		message = "mash v" + str(mash_version) + " found." 
		log.write_message(message, logger)
	except:
		sys.exit("mash not found.\n")  

#rasusa
	try:
		process = sp.Popen(["rasusa", "--version"], stdout=sp.PIPE, stderr=sp.PIPE) 
		rasusa_out, _ = process.communicate()
		rasusa_version = rasusa_out.decode()
		rasusa_version = rasusa_version.split("\n")[0].split(' ')[1]
		message ="rasusa v" + str(rasusa_version) + " found."
		log.write_message(message, logger)
	except:
		sys.exit("rasusa not found.\n")  
	
	# all dependencies found
	print("All dependencies found.")
	logger.info("All dependencies found.")



def validate_pacbio_model(pacbio_model, logger):
	"""Checks the input insta is really a fasta
	:param file: fasta file
    :return: 
    """
	
	message = "You have specified using a pacbio model for Flye with --pacbio_model. Checking the input."
	log.write_message(message, logger)

	if pacbio_model == "pacbio-raw":
		message = "You have selected pacbio-raw designed for PacBio regular CLR reads (<20% error)."
		log.write_message(message, logger)
		pacbio_model = "--pacbio-raw"
	elif pacbio_model == "pacbio-corr":
		message = "You have selected pacbio-corr designed for PacBio reads that were corrected with other methods (<3% error)."
		log.write_message(message, logger)
		pacbio_model = "--pacbio-corr"
	elif pacbio_model == "pacbio-hifi":
		message = "You have selected pacbio-hifi designed for PacBio HiFi reads (<1% error)."
		log.write_message(message, logger)
		pacbio_model = "--pacbio-hifi"
	else:
		sys.exit('You pacbio model was not pacbio-raw, pacbio-corr or pacbio-hifi. Please check your input and run plassembler again.')
	
	return pacbio_model
	
