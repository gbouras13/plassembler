import os
import sys
import subprocess as sp
import logging
import shutil

from src import log


def chopper(input_long_reads, out_dir, min_length, min_quality, gzip_flag, threads):
    """Filters long reads using chopper

    :param input_long_reads: input ONT reads file
    :param out_dir: output directory
    :param min_length: minimum length for long reads - defaults to 1000
    :param min_quality:  minimum quality for long reads - defaults to 8
    :param gzip_flag: whether or not the long reads are gzipped
    :return:
    """
    filtered_long_reads = os.path.join(out_dir, "chopper_long_reads.fastq.gz")
    f = open(filtered_long_reads, "w")
    if gzip_flag == True:
        try:
            unzip = sp.Popen(["gunzip", "-c", input_long_reads], stdout=sp.PIPE)
            chopper = sp.Popen(
                [
                    "chopper",
                    "-q",
                    min_quality,
                    "--threads",
                    threads,
                    "-l",
                    min_length,
                    "--headcrop",
                    "25",
                    "--tailcrop",
                    "25",
                ],
                stdin=unzip.stdout,
                stdout=sp.PIPE,
            )
            gzip = sp.Popen(["gzip"], stdin=chopper.stdout, stdout=f, stderr=sp.PIPE)
            output = gzip.communicate()[0]
        except:
            sys.exit("Error with chopper\n")
    else:
        try:
            cat = sp.Popen(["cat", input_long_reads], stdout=sp.PIPE)
            chopper = sp.Popen(
                [
                    "chopper",
                    "-q",
                    min_quality,
                    "--threads",
                    threads,
                    "-l",
                    min_length,
                    "--headcrop",
                    "25",
                    "--tailcrop",
                    "25",
                ],
                stdin=cat.stdout,
                stdout=sp.PIPE,
            )
            gzip = sp.Popen(["gzip"], stdin=chopper.stdout, stdout=f, stderr=sp.PIPE)
            output = gzip.communicate()[0]
        except:
            sys.exit("Error with chopper\n")


def rasusa(out_dir, subsample_flag, subsample_depth, chromosome_length, logger):
    """Subsets long reads using rasusa for faster Flye

    :param input_long_reads: input ONT reads file
    :param out_dir: output directory
    :param subsample_flag: whether user turns on subsampling (will just return chopper reads if False)
    :param subsample_depth:  depth of subsampling (defaults to 30x).
    :param chromosome_length: int chromosome length input by user
    :param logger: logger
    :return:
    """

    chopper_long_reads = os.path.join(out_dir, "chopper_long_reads.fastq.gz")
    subset_long_reads = os.path.join(out_dir, "final_filtered_long_reads.fastq.gz")

    # if chosen not to subset then just take chopper output
    if subsample_flag == False:
        shutil.copy2(chopper_long_reads, subset_long_reads)
    # use rasusa to subset to
    else:
        f = open(subset_long_reads, "w")
        try:
            rasusa = sp.run(
                [
                    "rasusa",
                    "-i",
                    chopper_long_reads,
                    "-s",
                    "13",
                    "--coverage",
                    str(subsample_depth),
                    "--genome-size",
                    str(chromosome_length),
                    "-O",
                    "g",
                ],
                stdout=f,
                stderr=sp.PIPE,
            )
            logger.log(logging.INFO, rasusa.stderr)
        except:
            sys.exit("Error with Rasusa\n")


def trim_short_read(short_one, short_two, out_dir, logger):
    """Trims short reads using fastp

    :param short_one:  R1 short read file
    :param short_two:  R2 short read file
    :param out_dir: output directory
    :param logger: logger
    :return:
    """
    out_one = os.path.join(out_dir, "trimmed_R1.fastq")
    out_two = os.path.join(out_dir, "trimmed_R2.fastq")
    try:
        fastp = sp.Popen(
            [
                "fastp",
                "--in1",
                short_one,
                "--in2",
                short_two,
                "--out1",
                out_one,
                "--out2",
                out_two,
            ],
            stdout=sp.PIPE,
            stderr=sp.PIPE,
        )
        log.write_to_log(fastp.stderr, logger)
    except:
        sys.exit("Error with Fastp\n")
