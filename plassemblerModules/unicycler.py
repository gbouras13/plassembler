import sys
import subprocess as sp
import plassemblerModules


### unicycler and deduplicating

def run_unicycler(short_only, threads, logger, short_one, short_two, long, unicycler_output_dir):
    """ runs Unicycler
    :param short_only: boolean flag whether or not this is short read only
	:param short_one: R1 short read fastq
    :param short_two: R2 short read fastq
    :param long: long read fastq
    :param unicycler_output_dir: unicycler Output Directory
    :param threads: threads
    :param logger: logger
    :return: 
    """
    try:
        if short_only == False:
            unicycler = sp.Popen(["unicycler", "-1", short_one, "-2", short_two, "-l", long, "-t", threads, "-o",  unicycler_output_dir ], stdout=sp.PIPE, stderr=sp.PIPE) 
        else: 
            unicycler = sp.Popen(["unicycler", "-1", short_one, "-2", short_two, "-t", threads, "-o",  unicycler_output_dir ], stdout=sp.PIPE, stderr=sp.PIPE) 
        plassemblerModules.write_to_log(unicycler.stdout, logger)
    except:
        sys.exit("Error with Unicycler.\n")  