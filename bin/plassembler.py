#!/usr/bin/env python3
import os
import logging
import time
import datetime
import input_commands
import qc
import run_flye
import extract
import case_one
import case_one_kmer
import case_three
import case_three_kmer
import depth
import extract
import cleanup
import run_mash
import install_database
import assembly
import sys

from version import __version__

v = __version__

if __name__ == "__main__":

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

    
    # instiate the output directory
    out_dir = input_commands.instantiate_dirs(args.outdir, args.force) # incase there is already an out_dir

    # beginning logging
    LOG_FILE = os.path.join(out_dir, prefix + "_" + str(time_for_log) + ".log")
    logger = logging.getLogger()
    logging.basicConfig(level=logging.INFO,filename=LOG_FILE,format='%(asctime)s - %(levelname)s - %(message)s')
    print("Starting plassembler v" + v)
    logger.info("Starting plassembler v" + v )

    # add the inputs to the log
    logging.info("Input args: %r", args)

    # check deps 
    print("Checking dependencies.")
    logger.info("Checking dependencies.")
    input_commands.check_dependencies(logger)


    # check the mash database is installed
    print("Checking database installation.")
    logger.info("Checking database installation.")
    database_installed = install_database.check_db_installation(args.database)
    if database_installed == True:
        print("Database successfully checked.")
        logger.info("Database successfully checked.")
    else:
        sys.exit("\nPlease run install_database.py \n") 

#############################
######### assembled_mode == true 
#############################

    if args.assembled_mode ==True:
        print("You have chosen to specify an input assembly FASTA file containing plasmids to calculate depth and PLSDB type. No assembly will be conducted.")
        logger.info("You have chosen to specify an input assembly FASTA file containing plasmids to calculate depth and PLSDB type. No assembly will be conducted.")
        print("###########################\nAssembled Mode Activated\n###########################")
        logger.info("###########################\nAssembled Mode Activated\n###########################")
        # validation
        print("Checking input FASTA.")
        logger.info("Checking input FASTA.")
        input_commands.validate_fasta(args.input)
        print("Checking input long fastqs.")
        logger.info("Checking input long fastqs.")
        long_zipped = input_commands.validate_fastq(args.longreads)
        print("Filtering long reads.")
        logger.info("Filtering long reads.")
        qc.nanofilt(args.longreads, out_dir, args.min_length, args.min_quality, long_zipped)
        # flag for whether there are short reads
        short_reads = False
        if args.short_one != "nothing" and args.short_two != "nothing":
            print("Checking input short read fastqs.")
            logger.info("Checking input short read fastqs")
            short_reads = True
            s1_zipped = input_commands.validate_fastq(args.short_one)
            s2_zipped = input_commands.validate_fastq(args.short_two)
            print('Trimming short reads.')
            logger.info("Trimming short reads.")
            qc.trim_short_read(args.short_one, args.short_two, out_dir,  logger)
        else:
            print("You have not input any short read fastqs. Proceeding with long only.")
            logger.info("You have not input any short read fastqs. Proceeding with long only.")
        
        print("Calculating Depths.")
        logger.info("Calculating Depths.")

        if short_reads == True:
            assembly.get_depth_assembly(args.input, out_dir, logger,  args.threads, prefix)
        else:
            assembly.get_depth_assembly_long_only(args.input, out_dir, logger,  args.threads, prefix)

        # run mash
        print('Calculating mash distances to PLSDB.')
        logger.info('Calculating mash distances to PLSDB.')
        assembly.mash_sketch_assembly(out_dir, args.input, logger)
        assembly.run_mash(out_dir, os.path.join(out_dir, "input.fasta.msh"), args.database, logger)
        mash_empty = assembly.process_mash_tsv_assembly(out_dir, args.database, prefix, args.input)
        # rename contigs and update copy bumber with plsdb
        assembly.rename_contigs_assembly(out_dir, args.input, prefix, short_reads) 
        cleanup.update_copy_number_summary_plsdb(out_dir, prefix, mash_empty)
        assembly.remove_intermediate_files(out_dir)

#############################
# not in assembled_mode#
#############################
    else:

        print("Checking input fastqs.")
        logger.info("Checking input fastqs")

        if args.kmer_mode == False:
            if args.short_one == 'nothing':
                logger.info("ERROR: You have running hybrid mode and have forgotten to specify short reads fastq files. Please try again and specify these with -1 and -2.")
                sys.exit("ERROR: You have running hybrid mode and have forgotten to specify short reads fastq files. Please try again and specify these with -1 and -2.")
        else:
            logger.info("You have chosen --kmer_mode with long reads only. Ignoring any short reads.")
            print("You have chosen --kmer_mode with long reads only. Ignoring any short reads.")


        # checking fastq 
        long_zipped = input_commands.validate_fastq(args.longreads)
        if args.kmer_mode == False:
            s1_zipped = input_commands.validate_fastq(args.short_one)
            s2_zipped = input_commands.validate_fastq(args.short_two)

        # filtering long readfastq
        print("Filtering long reads.")
        logger.info("Filtering long reads.")
        qc.nanofilt(args.longreads, out_dir, args.min_length, args.min_quality, long_zipped)

        # running Flye
        print("Running Flye.")
        logger.info("Running Flye")
        run_flye.run_flye( out_dir, args.threads,args.raw_flag, logger)

        # count contigs
        print("Counting Contigs.")
        logger.info("Counting Contigs")
        contig_count = run_flye.contig_count(out_dir)

        # flag for extracting the plasmid chromosomes
        # no_plasmids_flag = false means there are plasmids
        no_plasmids_flag = False

        ####################################################################
        # Case 1: where there is only 1 contig -> no plasmids in the long read only assembly
        ####################################################################

        if contig_count == 1:
            logger.info("Only one contig was assembled with Flye.")
            print("Only one contig was assembled with Flye.")

            if args.kmer_mode == False:
                print('Plassembler will now try to use short reads to find possible plasmids.')
                logger.info("Plassembler will now try to use short reads to find possible plasmids.")
            else:
                print('Plassembler will now try to use Unicycler to find possible plasmids that Flye may have missed in the long reads.')
                logger.info('Plassembler will now try to use Unicycler to find possible plasmids that Flye may have missed in the long reads.')


            # no_plasmids_flag = True as no plasmids
            no_plasmids_flag = True
            chromosome_flag = extract.extract_chromosome(out_dir, args.chromosome, no_plasmids_flag)
            
            if args.kmer_mode == False:
                print('Trimming short reads.')
                logger.info("Trimming short reads.")
                qc.trim_short_read(args.short_one, args.short_two, out_dir,  logger)

                print('Recovering possible plasmids from short reads.')
                logger.info("Recovering possible plasmids from short reads.")
                successful_unicycler_recovery = case_one.case_one(out_dir, args.threads, logger)
            else:
                print('Recovering possible plasmids using Unicycler.')
                logger.info("Recovering possible plasmids using Unicycler.")
                successful_unicycler_recovery = case_one_kmer.case_one_kmer(out_dir, args.threads, logger)


            # if unicycler successfully finished, calculate the plasmid copy numbers
            if successful_unicycler_recovery == True:

                print('Unicycler identified plasmids. Calculating Plasmid Copy Numbers.')
                logger.info("Unicycler identified plasmids. Calculating Plasmid Copy Numbers.")
                if args.kmer_mode == False:
                    depth.get_depth(out_dir, logger,  args.threads, prefix)
                else:
                    depth.get_depth_kmer(out_dir, logger,  args.threads, prefix)

                # run mash
                print('Calculating mash distances to PLSDB.')
                logger.info('Calculating mash distances to PLSDB.')
                run_mash.mash_sketch(out_dir, logger)
                run_mash.run_mash(out_dir, os.path.join(out_dir,"unicycler_output", "assembly.fasta.msh"), args.database, logger)
                mash_empty = run_mash.process_mash_tsv(out_dir, args.database, prefix)
                # rename contigs and update copy bumber with plsdb
                cleanup.rename_contigs(out_dir, prefix)
                cleanup.update_copy_number_summary_plsdb(out_dir, prefix, mash_empty)

                cleanup.move_and_copy_files(out_dir, prefix, successful_unicycler_recovery)
                cleanup.remove_intermediate_files(out_dir)
            ####################################################################
            # Case 4: where there are truly no plasmids
            ####################################################################
            else: # unicycler did not successfully finish, just touch the files empty for downstream (snakemake)
                print('No plasmids found.')
                logger.info("No plasmids found.")
                cleanup.move_and_copy_files(out_dir, prefix, successful_unicycler_recovery)
                cleanup.remove_intermediate_files(out_dir)

        # where more than 1 contig was assembled
        else:
            logger.info("More than one contig was assembled with Flye.")
            print("More than one contig was assembled with Flye.")
            print("Extracting Chromosome.")
            logger.info("Extracting Chromosome.")
            chromosome_flag = extract.extract_chromosome(out_dir, args.chromosome, no_plasmids_flag)
            ####################################################################
            # Case 2 - where no chromosome was identified (likely below required depth) - need more long reads or user got chromosome parameter wrong - exit plassembler
            ####################################################################
            if chromosome_flag == False:
                print('No chromosome was idenfitied. Likely, there was insufficient long read depth for Flye to assemble a chromosome. Increasing sequencing depth is recommended. Also please check your -c or --chromosome parameter, it may be too high. ')
                logger.info("No chromosome was idenfitied. Likely, there was insufficient long read depth for Flye to assemble a chromosome. Increasing sequencing depth is recommended. Also please check your -c or --chromosome parameter, it may be too high.")
                cleanup.move_and_copy_files(out_dir, prefix, chromosome_flag)
                cleanup.remove_intermediate_files(out_dir)
            ####################################################################
            # Case 3 - where a chromosome and plasmids were identified in the Flye assembly -> mappeed to plasmids, unmapped to chromosome and assembly
            ####################################################################
            else:
                if args.kmer_mode == False:
                    print('Chromosome Identified. Plassembler will now use long and short reads to assemble plasmids accurately.')
                    logger.info("Chromosome Identified. Plassembler will now use both long and short reads to assemble plasmids accurately.")
                    
                    print('Trimming short reads.')
                    logger.info("Trimming short reads.")
                    qc.trim_short_read(args.short_one, args.short_two, out_dir,  logger)
                    ##### modules
                    # assembly plasmids
                    case_three.case_three(out_dir, args.threads,logger)
                    # get copy number 
                    print('Calculating Plasmid Copy Numbers.')
                    logger.info("Calculating Plasmid Copy Numbers.")
                    depth.get_depth(out_dir, logger,  args.threads, prefix)

                # kmer_mode
                else:
                    print('Chromosome Identified. Plassembler will now use long reads and Unicycler to assemble plasmids accurately.')
                    logger.info("Chromosome Identified. Plassembler will now use long reads and Unicycler to assemble plasmids accurately.")
                    case_three_kmer.case_three_kmer(out_dir, args.threads,logger)
                    print('Calculating Plasmid Copy Numbers.')
                    logger.info("Calculating Plasmid Copy Numbers.")
                    depth.get_depth_kmer(out_dir, logger,  args.threads, prefix)

                # run mash
                print('Calculating mash distances to PLSDB.')
                logger.info('Calculating mash distances to PLSDB.')
                run_mash.mash_sketch(out_dir, logger)
                run_mash.run_mash(out_dir, os.path.join(out_dir,"unicycler_output", "assembly.fasta.msh"), args.database, logger)
                mash_empty = run_mash.process_mash_tsv(out_dir, args.database, prefix)

                # rename contigs and update copy bumber with plsdb
                if args.kmer_mode == False:
                    cleanup.rename_contigs(out_dir, prefix)
                else:
                    cleanup.rename_contigs_kmer(out_dir, prefix)
                cleanup.update_copy_number_summary_plsdb(out_dir, prefix, mash_empty)

                cleanup.move_and_copy_files(out_dir, prefix, chromosome_flag)
                cleanup.remove_intermediate_files(out_dir)


    # Determine elapsed time
    elapsed_time = time.time() - start_time
    elapsed_time = round(elapsed_time, 2)

    # Show elapsed time for the process
    logger.info("plassembler has finished")
    logger.info("Elapsed time: "+str(elapsed_time)+" seconds")

    print("plassembler has finished")
    print("Elapsed time: "+str(elapsed_time)+" seconds")
