import processes
import os

def recover_plasmids(out_dir, threads, logger):
    #### indexing contigs
    print('Indexing Chromosome.')
    logger.info("Indexing Chromosome.")
    processes.index_fasta(os.path.join(out_dir, "chromosome.fasta"),  logger)

    ##### Mapping #######

    #### long reads mapping
    print('Mapping Long Reads to Chromosome.')
    logger.info('Mapping Long Reads to Chromosome.')
    processes.minimap_long_reads(True, out_dir, threads, logger)

    #### short reads mapping
    print('Mapping Short Reads to Chromosome Contig')
    logger.info('Mapping Short Reads to Chromosome Contig')
    processes.bwa_map_short_reads( out_dir, True, threads,  logger)

    #### Processing bams ######
    print('Processing Bams.')
    logger.info('Processing Bams.')

    processes.sam_to_bam( out_dir, "chromosome_long", threads,  logger)
    processes.sam_to_bam( out_dir, "chromosome_short", threads,  logger)

    ### extractng mapped_unmapped bams ###
    processes.bam_to_mapped_or_unmapped(out_dir, "chromosome_long", threads, logger)
    processes.bam_to_mapped_or_unmapped(out_dir, "chromosome_short", threads, logger)

    ### extracting fastqs
    print('Extracting Fastqs.')
    logger.info('Extracting Fastqs.')

    processes.extract_long_fastq(out_dir, "chromosome", logger)
    processes.extract_short_fastq( out_dir, "chromosome",  threads,  logger)   

    # running unicycler
    print('Running Unicycler')
    logger.info('Running Unicycler')
    # short only flag
    short_r1 = os.path.join(out_dir, "unmapped_chromosome_R1.fastq.gz")
    short_r2 = os.path.join(out_dir, "unmapped_chromosome_R2.fastq.gz")
    long_reads = os.path.join(out_dir, "long_read_chromosome_unmapped.fastq")
    successful_unicycler = True
    processes.unicycler(False, threads, logger, short_r1, short_r2, long_reads, os.path.join(out_dir, "unicycler_output"))
    # check if unicycler completely (usually it won't)
    successful_unicycler = os.path.isfile(os.path.join(out_dir, "unicycler_output", "assembly.fasta")) 
    if successful_unicycler == False:
        print('Error with Unicycler - Likely due to insufficient Coverage (likely, there is truly no plasmids)')
        logger.info('Error with Unicycler - Likely due to insufficient Coverage (likely, there is truly no plasmids)')
    return successful_unicycler

