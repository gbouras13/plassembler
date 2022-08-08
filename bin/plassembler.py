#!/usr/bin/env python3
import input_commands
import processes
import os
import logging
import time
import datetime
import depth


if __name__ == "__main__":

    print("Starting plassembler.")

    # get start time
    start_time = time.time()

    # getting time for log file 

    time_for_log = datetime.datetime.now().strftime("%m%d%Y_%H%M%S")

    args = input_commands.get_input()
        # set the prefix
    if args.prefix == "Default":
        prefix = "plassembler"
    else:
        prefix = args.prefix
    
    out_dir = input_commands.instantiate_dirs(args.outdir, args.force) # incase there is already an outdir

    LOG_FILE = os.path.join(args.outdir, prefix + "_" + str(time_for_log) + ".log")
    logger = logging.getLogger()
    logging.basicConfig(level=logging.INFO,filename=LOG_FILE,format='%(asctime)s - %(levelname)s - %(message)s')
    logger.info("Starting plassembler")
    print("Checking input fastqs.")
    logger.info("Checking input fastqs")

    # instantiation/checking fastq 
    long_zipped = input_commands.validate_fastq(args.longreads)
    s1_zipped = input_commands.validate_fastq(args.short_one)
    s2_zipped = input_commands.validate_fastq(args.short_two)

    print("Filtering and Trimming long reads.")
    logger.info("Filtering and Trimming long reads.")

    processes.nanofilt(args.longreads, out_dir, args.min_length, args.min_quality, long_zipped)

    print("Running Flye.")
    logger.info("Running Flye")
    processes.run_flye( out_dir, args.threads, logger)

    print("Counting Contigs.")
    logger.info("Counting Contigs")
    contig_count = processes.contig_count(args.outdir)
    # flag for creating output if the program fails 
    fail = False

    if contig_count == 1:
        logger.info("Only one contig was assembled. There are no plasmids.")
        print("Only one contig was assembled. There are no plasmids.")
        fail = True 
        processes.move_and_copy_files(out_dir, prefix, fail)
        processes.remove_intermediate_files(out_dir)
    else:
        print("Extracting Chromosome.")
        logger.info("Extracting Chromosome.")
        chromosome_cirularised_flag = processes.extract_chromosome(args.outdir, args.chromosome)
        if chromosome_cirularised_flag == False:
            print('Insufficient long read depth for chromosome to circularise. Increasing sequencing depth is recommended.')
            logger.info("Insufficient long read depth for chromosome to circularise. Increasing sequencing depth is recommended.")
            fail = True 
            processes.move_and_copy_files(out_dir, prefix, fail)
            processes.remove_intermediate_files(out_dir)
        else:
            print('Trimming short reads.')
            logger.info("Trimming short reads.")
            processes.trim_short_read(args.short_one, args.short_two, out_dir,  logger)
            ##### modules
            # assembly plasmids
            processes.plasmid_assembly(out_dir, args.threads,logger)
            # get double mapped regions
            processes.double_mapping_analysis(out_dir, args.threads,logger)
            # get copy number 
            print('Calculating Plasmid Copy Numbers.')
            logger.info("Calculating Plasmid Copy Numbers.")
            depth.get_depth(out_dir, logger,  args.threads)
            processes.move_and_copy_files(out_dir, prefix, fail)
            processes.remove_intermediate_files(out_dir)


    # Determine elapsed time
    elapsed_time = time.time() - start_time
    elapsed_time = round(elapsed_time, 2)

    # Show elapsed time for the process
    logger.info("plassembler has finished")
    logger.info("Elapsed time: "+str(elapsed_time)+" seconds")

    print("plassembler has finished")
    print("Elapsed time: "+str(elapsed_time)+" seconds")

    




