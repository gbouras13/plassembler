import os
import sys
import subprocess as sp
import logging



def nanofilt(input_long_reads, out_dir, min_length, min_quality, gzip_flag):
    """Filters long reads using nanofilt

    :param input_long_reads: input ONT reads file
    :param out_dir: output directory
    :param min_length: minimum length for long reads - defaults to 1000
    :param min_quality:  minimum quality for long reads - defaults to 8
    :param gzip_flag: whether or not the long reads are gzipped
    :return:
    """
    filtered_long_reads = os.path.join(out_dir, "filtered_long_reads.fastq.gz")
    f = open(filtered_long_reads, "w")
    if gzip_flag == True:
        try:
            unzip = sp.Popen(["gunzip", "-c", input_long_reads ], stdout=sp.PIPE) 
            nanofilt = sp.Popen(["NanoFilt", "-q", min_quality, "-l", min_length, "--headcrop", "50"  ], stdin=unzip.stdout, stdout=sp.PIPE ) 
            gzip = sp.Popen(["gzip" ], stdin=nanofilt.stdout,stdout=f,stderr=sp.PIPE ) 
            output = gzip.communicate()[0]
        except:
            sys.exit("Error with nanofilt\n")  
    else:
        try:
            cat = sp.Popen(["cat", input_long_reads ], stdout=sp.PIPE) 
            nanofilt = sp.Popen(["NanoFilt", "-q", min_quality, "-l", min_length, "--headcrop", "50"  ], stdin=cat.stdout, stdout=sp.PIPE ) 
            gzip = sp.Popen(["gzip" ], stdin=nanofilt.stdout,stdout=f,stderr=sp.PIPE ) 
            output = gzip.communicate()[0]
        except:
            sys.exit("Error with nanofilt\n")  

def trim_short_read(short_one, short_two, out_dir,  logger):
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
        fastp = sp.Popen(["fastp", "--in1", short_one, "--in2", short_two, "--out1", out_one, "--out2", out_two ], stdout=sp.PIPE, stderr=sp.PIPE) 
        logging.write_to_log(fastp.stderr, logger)
    except:
        sys.exit("Error with Fastp\n")  