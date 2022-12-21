import sys
import os
import logging
import time
import datetime
import plassemblerModules

def main(argv):

    # get start time
    start_time = time.time()

    # getting time for log file 

    time_for_log = datetime.datetime.now().strftime("%m%d%Y_%H%M%S")

    args = plassemblerModules.get_input()

    # set the prefix
    if args.prefix == "Default":
        prefix = "plassembler"
    else:
        prefix = args.prefix
    
    # instiate the output directory
    out_dir = plassemblerModules.instantiate_dirs(args.outdir, args.force) # incase there is already an out_dir
    #out_dir = args.outdir

    # beginning logging
    LOG_FILE = os.path.join(out_dir, prefix + "_" + str(time_for_log) + ".log")
    logger = logging.getLogger()
    logging.basicConfig(level=logging.INFO,filename=LOG_FILE,format='%(asctime)s - %(levelname)s - %(message)s')
    print("Starting plassembler.")
    logger.info("Starting plassembler")
    print("Checking input fastqs.")
    logger.info("Checking input fastqs")

    # checking fastq 
    long_zipped = plassemblerModules.validate_fastq(args.longreads)
    s1_zipped = plassemblerModules.validate_fastq(args.short_one)
    s2_zipped = plassemblerModules.validate_fastq(args.short_two)

    # filtering long readfastq
    print("Filtering long reads.")
    logger.info("Filtering long reads.")
    plassemblerModules.nanofilt(args.longreads, out_dir, args.min_length, args.min_quality, long_zipped)

    # running Flye
    print("Running Flye.")
    logger.info("Running Flye")
    plassemblerModules.run_flye( out_dir, args.threads,args.raw_flag, logger)

    # count contigs
    print("Counting Contigs.")
    logger.info("Counting Contigs")
    contig_count = plassemblerModules.contig_count(out_dir)

    # flag for extracting the plasmid chromosomes
    # no_plasmids_flag = false means there are plasmids
    no_plasmids_flag = False

    ####################################################################
    # Case 1: where there is only 1 contig -> no plasmids in the long read only assembly
    ####################################################################

    if contig_count == 1:
        logger.info("Only one contig was assembled with Flye.")
        print("Only one contig was assembled with Flye.")
        print('Plassembler will now try to leverage short reads to assemble plasmids.')
        logger.info("Plassembler will now try to leverage short reads to assemble plasmids.")

        # no_plasmids_flag = True as no plasmids
        no_plasmids_flag = True
        chromosome_flag = plassemblerModules.extract_chromosome(out_dir, args.chromosome, no_plasmids_flag)
        
        print('Trimming short reads.')
        logger.info("Trimming short reads.")
        plassemblerModules.trim_short_read(args.short_one, args.short_two, out_dir,  logger)

        print('Recovering possible plasmids from short reads.')
        logger.info("Recovering possible plasmids from short reads.")
        successful_unicycler_recovery = plassemblerModules.case_one(out_dir, args.threads, logger)

        # if unicycler worked, calculate the plasmid copy numbers
        if successful_unicycler_recovery == True:
            print('Unicycler found plasmids. Calculating Plasmid Copy Numbers.')
            logger.info("Unicycler found plasmids. Calculating Plasmid Copy Numbers.")
            plassemblerModules.get_depth(out_dir, logger,  args.threads, prefix)
            #plassemblerModules.move_and_copy_files(out_dir, prefix, successful_unicycler_recovery)
            #plassemblerModules.remove_intermediate_files(out_dir)
        else: # successful recovery false, just touch the files empty for downstream (snakemake)
            print('placeholder')
            #plassemblerModules.move_and_copy_files(out_dir, prefix, successful_unicycler_recovery)
            #plassemblerModules.remove_intermediate_files(out_dir)

    # where more than 1 contig was assembled
    else:
        logger.info("More than one contig was assembled with Flye.")
        print("More than one contig was assembled with Flye.")
        print("Extracting Chromosome.")
        logger.info("Extracting Chromosome.")
        chromosome_flag = plassemblerModules.extract_chromosome(out_dir, args.chromosome, no_plasmids_flag)
        ####################################################################
        # Case 2 - where no chromosome was identified (below read length) - need more long reads - exit plassembler
        ####################################################################
        if chromosome_flag == False:
            print('Insufficient long read depth for Flye assembly to generate chromosome. Please check you -c or --chromosome value. Increasing sequencing depth is recommended.')
            logger.info("Insufficient long read depth for Flye assembly to generate chromosome. Please check you -c or --chromosome value. Increasing sequencing depth is recommended.")
            #plassemblerModules.move_and_copy_files(out_dir, prefix, chromosome_flag)
            #plassemblerModules.remove_intermediate_files(out_dir)
        ####################################################################
        # Case 3 - where a chromosome and plasmids were identified in the Flye assembly
        ####################################################################
        else:
            print('Trimming short reads.')
            logger.info("Trimming short reads.")
            plassemblerModules.trim_short_read(args.short_one, args.short_two, out_dir,  logger)
            ##### modules
            # assembly plasmids
            plassemblerModules.case_three(out_dir, args.threads,logger)
            # get copy number 
            print('Calculating Plasmid Copy Numbers.')
            logger.info("Calculating Plasmid Copy Numbers.")
            # assumes unicycler works - write a test if not
            plassemblerModules.get_depth(out_dir, logger,  args.threads, prefix)
            #plassemblerModules.move_and_copy_files(out_dir, prefix, chromosome_flag)
            #plassemblerModules.remove_intermediate_files(out_dir)


    # Determine elapsed time
    elapsed_time = time.time() - start_time
    elapsed_time = round(elapsed_time, 2)

    # Show elapsed time for the process
    logger.info("plassembler has finished")
    logger.info("Elapsed time: "+str(elapsed_time)+" seconds")

    print("plassembler has finished")
    print("Elapsed time: "+str(elapsed_time)+" seconds")


def run():
    main(sys.argv)


