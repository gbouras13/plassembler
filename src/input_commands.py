import os
import sys
import gzip
from Bio import SeqIO
import shutil
import subprocess as sp
from loguru import logger


### GLOBAL VARIABLES

def instantiate_dirs(output_dir, force):
    """checks that the output directory doesn't already exist, overwrites if forced
        :param output_dir: output directory path
    :param force: flag to overwrite the output directory
    :return: output_dir
    """
    # remove outdir on force
    if force == True:
        if os.path.isdir(output_dir) == True:
            #shutil.rmtree(output_dir)
            print('l')
        else:
            logger.info(f"--force was specified even though the outdir does not already exist. Continuing ")
    else:
        if os.path.isdir(output_dir) == True:
            sys.exit(
                "\nOutput directory already exists and force was not specified. Please specify -f or --force to overwrite the output directory. \n"
            )
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
                logger.info(f"FASTQ {file} checked")
            else:
                logger.error(f"Input file {file} is not in the FASTQ format.")
    else:
        zipped = False
        with open(file, "r") as handle:
            fastq = SeqIO.parse(handle, "fastq")
            if any(fastq):
                logger.info(f"FASTQ {file} checked")
            else:
                logger.error(f"Input file {file} is not in the FASTQ format.")
    return zipped


def validate_fasta(filename):
    """Checks the input insta is really a fasta
        :param file: fasta file
    :return:
    """
    with open(filename, "r") as handle:
        fasta = SeqIO.parse(handle, "fasta")
        if any(fasta):
            logger.info(f"FASTA {filename} checked")
        else:
            logger.error(f"Input file {filename} is not in the FASTA format.")


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
            logger.error(f"Error: There are multiple contigs in your chromosome FASTA {input_chromosome}. Please input a completed chromosome..")

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
        logger.info("You have input long read FASTQs for depth calculation.")
        long_gzipped = validate_fastq(longreads)
        long_flag = True

    # short
    short_flag = False
    if short_one != "nothing" and short_two != "nothing":
        logger.info("You have input paired short read FASTQs for depth calculation.")
        s1_gzipped = validate_fastq(short_one)
        s2_gzipped = validate_fastq(short_two)
        if s1_gzipped != s2_gzipped:
            logger.error(
                "R1 and R2 files are inconsistenly compressed. Please check the compression format and try again."
            )
        short_flag = True

    if short_flag == False and long_flag == False:
        logger.error(
            "No valid long read or paired short read FASTQs were input. Please check your input and try again."
        )

    if (short_one != "nothing" and short_two == "nothing") or (
        short_one == "nothing" and short_two != "nothing"
    ):
        logger.error(
            "Only 1 short read file was found. Please check your input and try again."
        )

    return (short_flag, long_flag, long_gzipped)


def check_dependencies():
    """Checks the version of Unicycler, spades and Raven
    :return:
    """

    # Flye
    try:
        process = sp.Popen(["flye", "--version"], stdout=sp.PIPE, stderr=sp.STDOUT)
        flye_out, _ = process.communicate()
        flye_out = flye_out.decode().strip()
        flye_major_version = int(flye_out.split(".")[0])
        flye_minor_version = int(flye_out.split(".")[1])
        flye_minorest_version = flye_out.split(".")[2]
    except:
        logger.error("Flye not found. Please reinstall Plassembler.")

    message = (
        "Flye version found is v"
        + str(flye_major_version)
        + "."
        + str(flye_minor_version)
        + "."
        + flye_minorest_version
        + "."
    )
    logger.info(message)

    if flye_major_version != 2:
        message = "Flye is too old - please reinstall plassembler, see instructions at https://github.com/gbouras13/plassembler."
        logger.error(message)
    if flye_minor_version < 9:
        message = "Flye is too old - please reinstall plassembler, see instructions at https://github.com/gbouras13/plassembler."
        logger.error(message)

    message = "Flye version is ok."
    logger.info(message)

    # raven
    try:
        process = sp.Popen(["raven", "--version"], stdout=sp.PIPE, stderr=sp.PIPE)
        raven_out, _ = process.communicate()
        raven_version = raven_out.decode()
        raven_version = raven_version.split("\n")[0]
        message = "Raven v" + str(raven_version) + " found."
        logger.info(message)
        message = "Raven version is ok."
        logger.info(message)
    except:
        logger.error("Raven not found")

    # unicycler
    try:
        process = sp.Popen(["unicycler", "--version"], stdout=sp.PIPE, stderr=sp.STDOUT)
        unicycler_out, _ = process.communicate()
        unicycler_out = unicycler_out.decode()
        unicycler_version = unicycler_out.split(" ")[1]
        # get rid of the "v"
        unicycler_version = unicycler_version[1:]

        unicycler_major_version = int(unicycler_version.split(".")[0])
        unicycler_minor_version = int(unicycler_version.split(".")[1])
        unicycler_minorest_version = int(unicycler_version.split(".")[2])
    except:
        message = "Unicycler not found. Please reinstall Plassembler, see instructions at https://github.com/gbouras13/plassembler."
        logger.error(message)

    message = (
        "Unicycler version found is v"
        + str(unicycler_major_version)
        + "."
        + str(unicycler_minor_version)
        + "."
        + str(unicycler_minorest_version)
        + "."
    )
    logger.info(message)

    if unicycler_minor_version < 4:
        message = "Unicycler is too old - please reinstall Plassembler, see instructions at https://github.com/gbouras13/plassembler."
        logger.error(message)
    elif unicycler_minor_version == 4 and unicycler_minorest_version < 8:
        message = "Unicycler is too old - please reinstall Plassembler, see instructions at https://github.com/gbouras13/plassembler."
        logger.error(message)
    elif unicycler_minor_version == 4 and unicycler_minorest_version >= 8:
        message = "Unicycler version is older than v0.5.0 - Plassembler will continue but please consider installing Unicycler v0.5.0. See instructions at https://github.com/gbouras13/plassembler."
        logger.info(message)
    else:
        message = "Unicycler version is ok."
        logger.info(message)

    # spades
    try:
        process = sp.Popen(["spades.py", "--version"], stdout=sp.PIPE, stderr=sp.PIPE)
        spades_out, _ = process.communicate()
        spades_out = spades_out.decode()
        spades_version = spades_out.split(" ")[3]
        spades_version = spades_version.split("\n")[0]
        message = "SPAdes " + str(spades_version) + " found."
        logger.info(message)
    except:
        logger.error("SPAdes not found.")

    # samtools
    try:
        process = sp.Popen(["samtools", "--version"], stdout=sp.PIPE, stderr=sp.PIPE)
        samtools_out, _ = process.communicate()
        samtools_out = samtools_out.decode()
        samtools_version = samtools_out.split("\n")[0].split(" ")[
            1
        ]  # get second line, and then second component of line
        message = "Samtools v" + str(samtools_version) + " found."
        logger.info(message)
    except:
        logger.error("Samtools not found.")
        

    # minimap2
    try:
        process = sp.Popen(["minimap2", "--version"], stdout=sp.PIPE, stderr=sp.PIPE)
        minimap2_out, _ = process.communicate()
        minimap2_version = minimap2_out.decode()
        minimap2_version = minimap2_version.split("\n")[0]
        message = "minimap2 v" + str(minimap2_version) + " found."
        logger.info(message)
    except:
        logger.error("minimap2 not found.")

    # fastp
    try:
        process = sp.Popen(["fastp", "--version"], stdout=sp.PIPE, stderr=sp.PIPE)
        _, fastp_out = process.communicate()
        fastp_version = fastp_out.decode()
        fastp_version = fastp_version.split("\n")[0].split(" ")[1]
        message = "fastp v" + str(fastp_version) + " found."
        logger.info(message)
    except:
        logger.error("fastp not found.")
        

    # chopper
    try:
        process = sp.Popen(["chopper", "--version"], stdout=sp.PIPE, stderr=sp.PIPE)
        chopper_out, _ = process.communicate()
        chopper_version = chopper_out.decode()
        chopper_version = chopper_version.split("\n")[0].split(" ")[1]
        message = "chopper v" + str(chopper_version) + " found."
        logger.info(message)
    except:
        logger.error("chopper not found.")

    # mash
    try:
        process = sp.Popen(["mash", "version"], stdout=sp.PIPE, stderr=sp.PIPE)
        mash_out, _ = process.communicate()
        mash_out = mash_out.decode()
        version_line = []
        for line in mash_out.split("\n"):
            if "version" in line:
                version_line.append(line)
        mash_version = version_line[0].split(" ")[2]
        message = "mash v" + str(mash_version) + " found."
        logger.info(message)
    except:
        logger.error("mash not found")

    # all dependencies found
    logger.info("All dependencies found.")


def validate_pacbio_model(pacbio_model):
    """Checks the input insta is really a fasta
        :param file: fasta file
    :return:
    """

    message = "You have specified using a pacbio model for Flye with --pacbio_model. Checking the input."
    logger.info(message)

    if pacbio_model == "pacbio-raw":
        message = "You have selected pacbio-raw designed for PacBio regular CLR reads (<20% error)."
        logger.info(message)
        pacbio_model = "--pacbio-raw"
    elif pacbio_model == "pacbio-corr":
        message = "You have selected pacbio-corr designed for PacBio reads that were corrected with other methods (<3% error)."
        logger.info(message)
        pacbio_model = "--pacbio-corr"
    elif pacbio_model == "pacbio-hifi":
        message = (
            "You have selected pacbio-hifi designed for PacBio HiFi reads (<1% error)."
        )
        logger.info(message)
        pacbio_model = "--pacbio-hifi"
    else:
        message = "You pacbio model was not pacbio-raw, pacbio-corr or pacbio-hifi. Please check your input and run plassembler again."
        logger.error(message)
        

    return pacbio_model
