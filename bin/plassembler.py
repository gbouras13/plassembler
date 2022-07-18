#!/usr/bin/env python3
import input_commands
import processes
import os
import subprocess as sp
import logging
import time
import datetime

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
    input_commands.validate_fastq(args.longreads)
    input_commands.validate_fastq(args.short_one)
    input_commands.validate_fastq(args.short_two)

    print("Running Flye.")
    logger.info("Running Flye")
    processes.run_flye(args.longreads, out_dir, args.threads, logger)

    print("Counting Contigs.")
    logger.info("Counting Contigs")
    contig_count = processes.contig_count(args.outdir)

    if contig_count == 1:
        logger.info("Only one contig was assembled. There are no plasmids.")
        print("Only one contig was assembled. There are no plasmids.")
    else:
        print("Extracting Chromosome.")
        chromosome_cirularised_flag = processes.extract_chromosome(args.outdir, args.chromosome)
        if chromosome_cirularised_flag == False:
            print('Insufficient long read depth for chromosome to circularise. Increasing sequencing depth is recommended.')
            logger.info("Insufficient long read depth for chromosome to circularise. Increasing sequencing depth is recommended.")
        else:
            print('Trimming short reads.')
            logger.info("Trimming short reads.")
            processes.trim_short_read(args.short_one, args.short_two, out_dir,  logger)
            print('Indexing Chromosome.')
            logger.info("Indexing Chromosome.")
            processes.index_chromosome( out_dir,  logger)
            print('Mapping Short Reads to Chromosome.')
            logger.info('Mapping Short Reads to Chromosome.')
            processes.bwa_map_chromosome( out_dir,args.threads,  logger)
            print('Converting Sam to Bam.')
            logger.info('Converting Sam to Bam.')
            processes.sam_to_bam( out_dir, args.threads,  logger)
            print('Extracting Unmapped Reads.')
            logger.info('Extracting Unmapped Reads.')
            processes.bam_to_unmap( out_dir, args.threads,  logger)
            processes.extract_unmap_fastq( out_dir, args.threads,  logger)
            print('Running Unicycler.')
            logger.info('Running Unicycler.')
            processes.unicycler( out_dir, args.threads,  logger)
            

            





    # # gene predictor
    # if args.gene_predictor == "phanotate":
    #     logger.info("Starting Phanotate")
    #     processes.run_phanotate(args.infile, out_dir, logger)
    # if gene_predictor == "prodigal":
    #     logger.info("Starting Prodigal")
    #     processes.run_prodigal(args.infile, out_dir, logger)

    # logger.info("Translating gene predicted fastas.")
    # processes.translate_fastas(out_dir,gene_predictor)
    # logger.info("Starting tRNA-scanSE")
    # processes.run_trna_scan(args.infile, out_dir, logger)

    # # set the db dir
    # if args.database == "Default":
    #     DBDIR = os.path.join(os.path.dirname(__file__),'../',"databases/")  
    # else:
    #     DBDIR = args.database

    # processes.remove_delim_fastas(out_dir,gene_predictor)

    # # runnin mmseqs2
    # logger.info("Starting mmseqs2")
    # processes.run_mmseqs(DBDIR, out_dir, args.threads, logger, gene_predictor)
    # logger.info("Starting hhsuite")
    # processes.run_hmmsuite(DBDIR, out_dir, args.threads, logger, args.gene_predictor)

    # # post processing
    # phan_mmseq_merge_df = post_processing.process_results(DBDIR, out_dir, prefix, gene_predictor)
    # logger.info("Post Processing Data")
    # length_df = post_processing.get_contig_name_lengths(args.infile, out_dir, prefix)
    # post_processing.create_gff(phan_mmseq_merge_df, length_df, args.infile, out_dir, prefix, locustag)
    # post_processing.create_tbl(phan_mmseq_merge_df, length_df, out_dir, prefix)
    # post_processing.create_txt(phan_mmseq_merge_df, length_df,out_dir, prefix)
    # logger.info("Converting gff to genbank using seqret")
    # processes.convert_gff_to_gbk(args.infile, out_dir, prefix, logger)
    
    # # delete tmp files
    # sp.run(["rm", "-rf", os.path.join(out_dir, "target_dir") ])
    # sp.run(["rm", "-rf", os.path.join(out_dir, "tmp_dir/") ])
    # sp.run(["rm", "-rf", os.path.join(out_dir, "mmseqs/") ])
    # sp.run(["rm", "-rf", os.path.join(out_dir, "cleaned_" + gene_predictor + ".tsv") ])
    # sp.run(["rm", "-rf", os.path.join(out_dir, "input_fasta_delim.fasta") ])
    # sp.run(["rm", "-rf", os.path.join(out_dir, "mmseqs_results.tsv") ])
    # sp.run(["rm", "-rf", os.path.join(out_dir, "top_hits_hhsuite.tsv") ])
    # sp.run(["rm", "-rf", os.path.join(out_dir, "top_hits_mmseqs.tsv") ])
    # sp.run(["rm", "-rf", os.path.join(out_dir, "hhsuite_target_dir") ])
    # sp.run(["rm", "-rf", os.path.join(out_dir, "phanotate_out.txt") ])
    # sp.run(["rm", "-rf", os.path.join(out_dir, "trnascan_out.gff") ])
    # sp.run(["rm", "-rf", os.path.join(out_dir, gene_predictor + "_aas_tmp.fasta") ])
    # sp.run(["rm", "-rf", os.path.join(out_dir, gene_predictor + "_out_tmp.fasta") ])

    # Determine elapsed time
    elapsed_time = time.time() - start_time

    # Show elapsed time for the process
    logger.info("plassembler has finished")
    logger.info("Elapsed time: "+str(elapsed_time)+" seconds")

    print("plassembler has finished")
    print("Elapsed time: "+str(elapsed_time)+" seconds")

    




