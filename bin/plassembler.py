#!/usr/bin/env python3
import os
import logging
import time
import datetime
import input_commands
import qc
import mapping
import extract
import case_one
import case_one_kmer
import case_three
import run_flye
import case_three_kmer
import depth
import extract
import cleanup
import run_mash
import install_database
import assembly
import sys
import bam
from plass_class import Plass
import sam_to_fastq
import concat
import deduplicate
import run_unicycler

from version import __version__

v = __version__

if __name__ == "__main__":

    # get start time
    start_time = time.time()

    # getting time for log file 
    time_for_log = datetime.datetime.now().strftime("%m%d%Y_%H%M%S")

    # get inputs
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
        input_commands.validate_fastas_assembled_mode(args.input_chromosome, args.input_plasmids)
        print("Checking input long fastqs.")
        logger.info("Checking input long fastqs.")
        long_zipped = input_commands.validate_fastq(args.longreads)
        print("Filtering long reads.")
        logger.info("Filtering long reads.")
        qc.chopper(args.longreads, out_dir, args.min_length, args.min_quality, long_zipped)
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
            assembly.get_depth_assembly(args.input_chromosome,args.input_plasmids, out_dir, logger,  args.threads, prefix)
        else:
            assembly.get_depth_assembly_long_only(args.input_chromosome,args.input_plasmids, out_dir, logger,  args.threads, prefix)

        input_fasta = os.path.join(out_dir,"input.fasta")
        # run mash
        print('Calculating mash distances to PLSDB.')
        logger.info('Calculating mash distances to PLSDB.')
        assembly.mash_sketch_assembly(out_dir, input_fasta, logger)
        assembly.run_mash(out_dir, os.path.join(out_dir, "input.fasta.msh"), args.database, logger)
        mash_empty = assembly.process_mash_tsv_assembly(out_dir, args.database, prefix, input_fasta)
        # rename contigs and update copy number with plsdb
        assembly.rename_contigs_assembly(out_dir, input_fasta, prefix, short_reads) 
        cleanup.update_copy_number_summary_plsdb(out_dir, prefix, mash_empty)
        assembly.remove_intermediate_files(out_dir)

#############################
# not in assembled_mode #
#############################

    else:

        print("Checking input fastqs.")
        logger.info("Checking input fastqs")

        # experimental kmer mode - high quality long read only
        if args.kmer_mode == False:
            if args.short_one == 'nothing':
                logger.info("ERROR: You have forgotten to specify short reads fastq files. Please try again and specify these with -1 and -2.")
                sys.exit("ERROR: You have forgotten to specify short reads fastq files. Please try again and specify these with -1 and -2.")
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
        #qc.chopper(args.longreads, out_dir, args.min_length, args.min_quality, long_zipped)

        # running Flye
        print("Running Flye.")
        logger.info("Running Flye")
        #run_flye.run_flye( out_dir, args.threads,args.raw_flag, logger)

        # instanatiate the class with some of the commands
        plass = Plass()
        plass.threads = args.threads
        plass.kmer = args.kmer_mode

        # count contigs and add to the object
        print("Counting Contigs.")
        logger.info("Counting Contigs")
        plass.get_contig_count(out_dir, logger)

        ####################################################################
        # Case 1: where there is only 1 contig -> no plasmids in the long read only assembly
        ####################################################################

        if plass.contig_count == 1:
            logger.info("Only one contig was assembled with Flye.")
            print("Only one contig was assembled with Flye.")

            # no_plasmids_flag = True as no plasmids
            plass.no_plasmids_flag = True

            # identifies chromosome and renames contigs
            plass.identify_chromosome_process_flye(out_dir, args.chromosome)

            # no chromosome identified - cleanup and exit
            if plass.chromosome_flag == False:
                message = 'No chromosome was identified. Likely, there was insufficient long read depth for Flye to assemble a chromosome. \nIncreasing sequencing depth is recommended. \nAlso please check your -c or --chromosome parameter, it may be too high. '
                print(message)
                logger.info(message)
                cleanup.move_and_copy_files(out_dir, prefix, plass.chromosome_flag)
            else: # chromosome identified -> move on 
                if args.kmer_mode == False:
                    message = 'Chromosome Identified. Plassembler will now use long and short reads to assemble plasmids accurately.'
                    print(message)
                    logger.info(message)

                    message = 'Trimming short reads.'
                    print(message)
                    logger.info(message)
                    #qc.trim_short_read(args.short_one, args.short_two, out_dir,  logger)

                    message = 'Mapping Long Reads.'
                    print(message)
                    logger.info(message)
                    #mapping.minimap_long_reads(out_dir, args.threads, logger)

                    #### short reads mapping
                    message = 'Mapping Short Reads.'
                    print(message)
                    logger.info(message)
                    #mapping.minimap_short_reads(out_dir, args.threads, logger)

                    message = 'Processing Sam/Bam Files and extracting Fastqs.'
                    print(message)
                    logger.info(message)
                  

                    # for long, custom function
                    # sam_to_fastq.extract_bin_long_fastqs(out_dir, args.multi_map)

                    # # short
                    # # for short, too slow so use samtools unless multimap is on
                    # if args.multi_map == True:
                    #     sam_to_fastq.extract_bin_short_fastqs(out_dir)
                    # else:
                    #     bam.sam_to_bam_short(out_dir, args.threads, logger)
                    #     bam.split_bams(out_dir, args.threads, logger)
                    #     bam.bam_to_fastq_short(out_dir, args.threads, logger)

                    # print('Concatenating and Deduplicating Fastqs.')
                    # logger.info('Concatenating and Deduplicating Fastqs.')

                    concat.concatenate_all_fastqs(out_dir,logger)

                    # running unicycler
                    message = 'Running Unicycler.'
                    print(message)
                    logger.info(message)

                    # short only flag False
                    short_r1 = os.path.join(out_dir, "short_read_concat_R1.fastq")
                    short_r2 = os.path.join(out_dir, "short_read_concat_R2.fastq")
                    long_reads = os.path.join(out_dir, "plasmid_long.fastq")

                    run_unicycler.run_unicycler(False, args.threads, logger, short_r1, short_r2, long_reads, 
                                                 os.path.join(out_dir, "unicycler_output"))


                else: # kmer mode
                    message = 'Plassembler will now try to use Unicycler to find possible plasmids that Flye may have missed in the long reads.'
                    print(message)
                    logger.info(message)

                    # maybe spades? test it out
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
                plass.process_mash_tsv(out_dir, args.database, prefix)
                # rename contigs and update copy bumber with plsdb
                cleanup.rename_contigs(out_dir, prefix)
                cleanup.update_copy_number_summary_plsdb(out_dir, prefix, mash_empty)

                cleanup.move_and_copy_files(out_dir, prefix, successful_unicycler_recovery)
                #cleanup.remove_intermediate_files(out_dir)
            ####################################################################
            # Case 4: where there are truly no plasmids even after unicycler runs
            ####################################################################
            else: # unicycler did not successfully finish, just touch the files empty for downstream (snakemake)
                print('No plasmids found.')
                logger.info("No plasmids found.")
                cleanup.move_and_copy_files(out_dir, prefix, successful_unicycler_recovery)
                #cleanup.remove_intermediate_files(out_dir)

        # where more than 1 contig was assembled
        else:
            logger.info("More than one contig was assembled with Flye.")
            print("More than one contig was assembled with Flye.")
            print("Extracting Chromosome.")
            logger.info("Extracting Chromosome.")

            # no_plasmids_flag = False as no plasmids
            plass.no_plasmids_flag = False

            # identifies chromosome and renames contigs
            plass.identify_chromosome_process_flye(out_dir, args.chromosome)

            ####################################################################
            # Case 2 - where no chromosome was identified (likely below required depth) - need more long reads or user got chromosome parameter wrong - exit plassembler
            ####################################################################
            if plass.chromosome_flag == False:
                message = 'No chromosome was idenfitied. please check your -c or --chromosome parameter, it may be too high. \nLikely, there was insufficient long read depth for Flye to assemble a chromosome. Increasing sequencing depth is recommended.'
                print(message)
                logger.info(message)
                cleanup.move_and_copy_files(out_dir, prefix, chromosome_flag)
                #cleanup.remove_intermediate_files(out_dir)
            ####################################################################
            # Case 3 - where a chromosome and plasmids were identified in the Flye assembly -> get reads mappeed to plasmids, unmapped to chromosome and assemble
            ####################################################################
            else:
                if args.kmer_mode == False:
                    message = 'Chromosome Identified. Plassembler will now use long and short reads to assemble plasmids accurately.'
                    print(message)
                    logger.info(message)

                    message = 'Trimming short reads.'
                    print(message)
                    logger.info(message)
                    #qc.trim_short_read(args.short_one, args.short_two, out_dir,  logger)

                    message = 'Mapping Long Reads.'
                    print(message)
                    logger.info(message)
                    #mapping.minimap_long_reads( out_dir, args.threads, logger)

                    #### short reads mapping
                    message = 'Mapping Short Reads.'
                    print(message)
                    logger.info(message)
                    #mapping.minimap_short_reads(out_dir, args.threads, logger)

                    message = 'Processing Sam/Bam Files and extracting Fastqs.'
                    print(message)
                    logger.info(message)
                  
                    # for long, custom function
                    # for short, too slow so use samtools
                    # sam_to_fastq.extract_bin_long_fastqs(out_dir, args.multi_map)

                    # # short
                    # # for short, too slow so use samtools unless multimap is on
                    # if args.multi_map == True:
                    #     sam_to_fastq.extract_bin_short_fastqs(out_dir)
                    # else:
                    #     bam.sam_to_bam_short(out_dir, args.threads, logger)
                    #     bam.split_bams(out_dir, args.threads, logger)
                    #     bam.bam_to_fastq_short(out_dir, args.threads, logger)
                    
                    # #print('Concatenating and Deduplicating Fastqs.')
                    # #logger.info('Concatenating and Deduplicating Fastqs.')

                    # concat.concatenate_all_fastqs(out_dir,logger)

                    # # running unicycler
                    # message = 'Running Unicycler.'
                    # print(message)
                    # logger.info(message)

                    # # short only flag False
                    # short_r1 = os.path.join(out_dir, "short_read_concat_R1.fastq")
                    # short_r2 = os.path.join(out_dir, "short_read_concat_R2.fastq")
                    # long_reads = os.path.join(out_dir, "plasmid_long.fastq")

                    # run_unicycler.run_unicycler(False, args.threads, logger, short_r1, short_r2, long_reads, 
                    #                              os.path.join(out_dir, "unicycler_output"))

                # kmer_mode
                else:
                    print('Chromosome Identified. Plassembler will now use long reads and Unicycler to assemble plasmids accurately.')
                    logger.info("Chromosome Identified. Plassembler will now use long reads and Unicycler to assemble plasmids accurately.")
                    case_three_kmer.case_three_kmer(out_dir, args.threads,logger)
                    print('Calculating Plasmid Copy Numbers.')
                    logger.info("Calculating Plasmid Copy Numbers.")
                    depth.get_depth_kmer(out_dir, logger,  args.threads, prefix)
                
                ##################################
                ##### get copy number depths
                ##################################
                message = 'Calculating Plasmid Copy Numbers.'
                print(message)
                logger.info(message)
                
                # as class so saves the depth dataframe nicely
                plass.get_depth(out_dir, logger,  args.threads, prefix)

                # run mash
                print('Calculating mash distances to PLSDB.')
                logger.info('Calculating mash distances to PLSDB.')
                # sketches the plasmids
                run_mash.mash_sketch(out_dir, logger)
                # runs mash 
                run_mash.run_mash(out_dir, os.path.join(out_dir,"unicycler_output", "assembly.fasta.msh"), args.database, logger)
                # processes output
                plass.process_mash_tsv(out_dir, args.database)
                # combine depth and mash tsvs
                plass.combine_depth_mash_tsvs(out_dir, prefix)


                # rename contigs and update copy bumber with plsdb
                if args.kmer_mode == False:
                    cleanup.rename_contigs(out_dir, prefix)
                else:
                    cleanup.rename_contigs_kmer(out_dir, prefix)
                
                cleanup.update_copy_number_summary_plsdb(out_dir, prefix, mash_empty)

                cleanup.move_and_copy_files(out_dir, prefix, chromosome_flag)
                #cleanup.remove_intermediate_files(out_dir)


    # Determine elapsed time
    elapsed_time = time.time() - start_time
    elapsed_time = round(elapsed_time, 2)

    # Show elapsed time for the process
    logger.info("plassembler has finished")
    logger.info("Elapsed time: "+str(elapsed_time)+" seconds")

    print("plassembler has finished")
    print("Elapsed time: "+str(elapsed_time)+" seconds")
