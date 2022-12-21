import plassemblerModules
import os

def case_one(out_dir, threads, logger):

    #### indexing chromosome
    print('Indexing Chromosome.')
    logger.info("Indexing Chromosome.")
    plassemblerModules.index_fasta(os.path.join(out_dir, "chromosome.fasta"),  logger)

    ######################
    ##### Mapping #######
    ######################

    #### long reads mapping to chromosome
    print('Mapping Long Reads to Chromosome.')
    logger.info('Mapping Long Reads to Chromosome.')
    plassemblerModules.minimap_long_reads(True, out_dir, threads, logger)

    #### short reads mapping to chromosome
    print('Mapping Short Reads to Chromosome.')
    logger.info('Mapping Short Reads to Chromosome.')
    plassemblerModules.bwa_map_short_reads( out_dir, True, threads,  logger)

    #### Processing bams ######
    print('Processing Bams.')
    logger.info('Processing Bams.')

    # convert sam to bam
    plassemblerModules.sam_to_bam(out_dir, "chromosome_long", threads,  logger)
    plassemblerModules.sam_to_bam(out_dir, "chromosome_short", threads,  logger)

    ### extractng unmapped bams ###
    plassemblerModules.bam_to_mapped_or_unmapped(out_dir, "chromosome_long", threads, logger)
    plassemblerModules.bam_to_mapped_or_unmapped(out_dir, "chromosome_short", threads, logger)


    ### extracting fastqs
    print('Extracting Fastqs.')
    logger.info('Extracting Fastqs.')

    plassemblerModules.extract_long_fastq(out_dir, "chromosome", logger)
    plassemblerModules.extract_short_fastq( out_dir, "chromosome",  threads,  logger)   

    # running unicycler
    print('Running Unicycler')
    logger.info('Running Unicycler')

    # short only flag
    short_R1 = os.path.join(out_dir, "unmapped_chromosome_R1.fastq.gz")
    short_R2 = os.path.join(out_dir, "unmapped_chromosome_R2.fastq.gz")
    long_reads = os.path.join(out_dir, "long_read_chromosome_unmapped.fastq")
    
    # is not short only
    plassemblerModules.run_unicycler(False, threads, logger, short_R1, short_R2, long_reads, os.path.join(out_dir, "unicycler_output"))
    
    # flag for successful unicycler run
    successful_unicycler = True
    # check if unicycler succeded according to the output (it won't if no plasmids)
    successful_unicycler = os.path.isfile(os.path.join(out_dir, "unicycler_output", "assembly.fasta")) 
    
    if successful_unicycler == False:
        print('Unicycler failed - there are likely no plasmids in your isolate.')
        logger.info('Unicycler failed - there are likely no plasmids in your isolate.')
    
    return successful_unicycler

