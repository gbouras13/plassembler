import os
import plassemblerModules


####################################################################
# Case 3 - where a chromosome and plasmids were identified in the Flye assembly
####################################################################
def case_three(out_dir, threads, logger):

    #### indexing 
    plassemblerModules.index_fasta( os.path.join(out_dir, "non_chromosome.fasta"),  logger)
    plassemblerModules.index_fasta( os.path.join(out_dir, "chromosome.fasta"),  logger)

    ##### Mapping #######

    #### long reads mapping to plasmids
    print('Mapping Long Reads to Putative Plasmid Contigs.')
    logger.info('Mapping Long Reads to Putative Plasmid Contigs.')
    plassemblerModules.minimap_long_reads(False, out_dir, threads, logger)
    
    #### long reads mapping to chromosome
    print('Mapping Long Reads to Chromosome.')
    logger.info('Mapping Long Reads to Chromosome.')
    plassemblerModules.minimap_long_reads(True, out_dir, threads, logger)

    #### short reads mapping to plasmids
    print('Mapping Short Reads to Putative Plasmid Contigs')
    logger.info('Mapping Short Reads to Putative Plasmid Contigs')
    plassemblerModules.bwa_map_short_reads( out_dir, False, threads,  logger)

    #### short reads mapping to chromosome
    print('Mapping Short Reads to Chromosome Contig')
    logger.info('Mapping Short Reads to Chromosome Contig')
    plassemblerModules.bwa_map_short_reads( out_dir, True, threads,  logger)

    #### Processing bams ######
    print('Processing Bams.')
    logger.info('Processing Bams.')

    plassemblerModules.sam_to_bam( out_dir, "chromosome_long", threads,  logger)
    plassemblerModules.sam_to_bam( out_dir, "non_chromosome_long", threads,  logger)
    plassemblerModules.sam_to_bam( out_dir, "chromosome_short", threads,  logger)
    plassemblerModules.sam_to_bam( out_dir, "non_chromosome_short", threads,  logger)

    ### extractng mapped_unmapped bams ###
    plassemblerModules.bam_to_mapped_or_unmapped(out_dir, "chromosome_long", threads, logger)
    plassemblerModules.bam_to_mapped_or_unmapped(out_dir, "non_chromosome_long", threads, logger)
    plassemblerModules.bam_to_mapped_or_unmapped(out_dir, "chromosome_short", threads, logger)
    plassemblerModules.bam_to_mapped_or_unmapped(out_dir, "non_chromosome_short", threads, logger)

    ### extracting fastqs
    print('Extracting Fastqs.')
    logger.info('Extracting Fastqs.')

    plassemblerModules.extract_long_fastq(out_dir, "chromosome", logger)
    plassemblerModules.extract_long_fastq(out_dir, "non_chromosome", logger)
    plassemblerModules.extract_short_fastq( out_dir, "chromosome",  threads,  logger)   
    plassemblerModules.extract_short_fastq( out_dir, "non_chromosome",  threads,  logger)   

    # concatenating and deduplicating fastqs

    print('Concatenating and Deduplicating Fastqs.')
    logger.info('Concatenating and Deduplicating Fastqs.')

    plassemblerModules.concatenate_fastqs(out_dir,logger)
    plassemblerModules.deduplicate_fastqs(out_dir, threads, logger)

    # running unicycler
    print('Running Unicycler')
    logger.info('Running Unicycler')

    # short only flag False
    short_r1 = os.path.join(out_dir, "short_read_dedup_R1.fastq.gz")
    short_r2 = os.path.join(out_dir, "short_read_dedup_R2.fastq.gz")
    long_reads = os.path.join(out_dir, "long_read_dedup.fastq")

    plassemblerModules.run_unicycler(False, threads, logger, short_r1, short_r2, long_reads, os.path.join(out_dir, "unicycler_output"))


