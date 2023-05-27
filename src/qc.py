import os
import sys
import subprocess as sp
from src.external_tools import ExternalTool
from loguru import logger
from pathlib import Path



def chopper(input_long_reads, out_dir, min_length, min_quality, gzip_flag, threads, logdir):
    """Filters long reads using chopper

    :param input_long_reads: input ONT reads file
    :param out_dir: output directory
    :param min_length: minimum length for long reads - defaults to 1000
    :param min_quality:  minimum quality for long reads - defaults to 8
    :param gzip_flag: whether or not the long reads are gzipped
    :param logdir
    :return:
    """
    filtered_long_reads = os.path.join(out_dir, "chopper_long_reads.fastq.gz")
    logger.info(f"Started running chopper")
    logdir.mkdir(parents=True, exist_ok=True)
    tool = "chopper"
    tool_name = Path(tool).name
    logfile_prefix: Path = logdir / f"{tool_name}"
    err_log = open(f"{logfile_prefix}.err", "w")
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
                stderr=err_log
            )
            gzip = sp.Popen(["gzip"], stdin=chopper.stdout, stdout=f)
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
                stderr=err_log
            )
            gzip = sp.Popen(["gzip"], stdin=chopper.stdout, stdout=f)
            output = gzip.communicate()[0]
        except:
            sys.exit("Error with chopper\n")
    logger.info(f"Finised running chopper")



def fastp(short_one, short_two, out_dir, logdir):
    """Trims short reads using fastp

    :param short_one:  R1 short read file
    :param short_two:  R2 short read file
    :param out_dir: output directory
    :param logger: logger
    :return:
    """
    out_one = os.path.join(out_dir, "trimmed_R1.fastq")
    out_two = os.path.join(out_dir, "trimmed_R2.fastq")

    fastp = ExternalTool(
        tool="fastp",
        input=f"--in1 {short_one} --in2 {short_two}",
        output=f"--out1 {out_one} --out2 {out_two}",
        params=f"",
        logdir=logdir,
        outfile=""
    )

    ExternalTool.run_tool(fastp, to_stdout=False)
