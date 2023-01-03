
#########################################
##### some code taken from & modified
##### https://github.com/rrwick/Small-plasmid-Nanopore/blob/main/scripts/get_depths.py
#########################################
from Bio import SeqIO
import os
import sys
import subprocess as sp
import statistics
import numpy as np
import pandas as pd
import mapping
import concat

def get_depth(out_dir, logger,  threads, prefix):
    """ wrapper function to get depth of each plasmid
    :param prefix: prefix (default plassembler)
    :param out_dir:  Output Directory
    :param threads: threads
    :param logger: logger
    :return: 
    """
    concatenate_chrom_plasmids(out_dir, logger)
    mapping.index_fasta(os.path.join(out_dir, "combined.fasta"),  logger)
    bwa_map_depth_sort(out_dir, threads)
    minimap_depth_sort(out_dir, threads)
    contig_lengths = get_contig_lengths(out_dir)
    depthsShort = get_depths_from_bam(out_dir, "short", contig_lengths)
    depthsLong = get_depths_from_bam(out_dir, "long", contig_lengths)
    circular_status = get_contig_circularity(out_dir)
    summaryDepthdf_short = collate_depths(depthsShort,"short",contig_lengths)
    summaryDepthdf_long = collate_depths(depthsLong,"long",contig_lengths)
    combine_depth_dfs(out_dir, summaryDepthdf_short, summaryDepthdf_long, prefix, circular_status)


def concatenate_chrom_plasmids(out_dir, logger):
    """ concatenates chromosome and plasmids
    :param out_dir:  Output Directory
    :param logger: logger
    :return: 
    """
    chrom_fasta =os.path.join(out_dir,"chromosome.fasta")
    plas_fasta = os.path.join(out_dir,"unicycler_output", "assembly.fasta")
    concat_fasta = open(os.path.join(out_dir, "combined.fasta"), "w")
    try:
        concat.concatenate_single(chrom_fasta, plas_fasta, concat_fasta, logger)
    except:
        sys.exit("Error with concatenate_fastas\n")  


# get lengths of contigs
def get_contig_lengths(out_dir):
    """ gets contig lengths of combined chrom and plasmids fastas
    :param out_dir:  Output Directory
    :return: contig_lengths: dictionary of headers and lengths
    """
    contig_lengths = {}
    for dna_record in SeqIO.parse(os.path.join(out_dir, "combined.fasta"), 'fasta'):
        plas_len = len(dna_record.seq)
        dna_header = dna_record.id
        contig_lengths[dna_header] = plas_len
    return contig_lengths



# get circular status of contigs
def get_contig_circularity(out_dir):
    """ gets circularity of contigs
    :param out_dir:  Output Directory
    :return: circular_status: dictionary of contig header and circular status
    """
    circular_status = {}
    # add circularity
    for dna_record in SeqIO.parse(os.path.join(out_dir, "combined.fasta"), 'fasta'):
        dna_header = dna_record.id
        # check if circular is in unicycler output description
        if "circular=true" in dna_record.description:
            circular_status[dna_header] = "circular"
        # circular chromsome
        elif "chromosome" in dna_record.id:
            circular_status[dna_header] = "circular"
        else:
            circular_status[dna_header] = "not_circular"
    return circular_status



def bwa_map_depth_sort(out_dir, threads):
    """ maps short reads using bwa to combined fasta and sorts bam
    :param out_dir:  Output Directory
    :return: threads: threads
    """
    trim_one = os.path.join(out_dir, "trimmed_R1.fastq")
    trim_two = os.path.join(out_dir, "trimmed_R2.fastq")
    fasta = os.path.join(out_dir, "combined.fasta")
    bam = os.path.join(out_dir, "combined_sorted.bam")
    try:
        bwa_map = sp.Popen(["bwa", "mem", "-t", threads, fasta, trim_one, trim_two ], stdout=sp.PIPE, stderr=sp.DEVNULL) 
        samtools_sort = sp.Popen(["samtools", "sort", "-@", threads, "-o", bam, "-" ], stdin=bwa_map.stdout, stderr=sp.DEVNULL ) 
        samtools_sort.communicate()[0]
    except:
        sys.exit("Error with bwa mem or samtools sort.\n")  


def minimap_depth_sort(out_dir, threads):
    """ maps long reads using minimap2 to combined fasta and sorts bam
    :param out_dir:  out_dir
    :param: threads: threads
    """
    input_long_reads = os.path.join(out_dir, "filtered_long_reads.fastq.gz")
    fasta = os.path.join(out_dir, "combined.fasta")
    bam = os.path.join(out_dir, "combined_sorted_long.bam")
    try:
        minimap = sp.Popen(["minimap2", "-ax", "map-ont", "-t", threads, fasta, input_long_reads ], stdout=sp.PIPE, stderr=sp.DEVNULL) 
        samtools_sort = sp.Popen(["samtools", "sort", "-@", threads, "-o", bam, "-" ], stdin=minimap.stdout, stderr=sp.DEVNULL ) 
        samtools_sort.communicate()[0]
    except:
        sys.exit("Error with mapping and sorting\n")  


def get_depths_from_bam(out_dir, shortFlag, contig_lengths):
    """ maps runs samtools depth on bam
    :param out_dir:  out_dir
    :param: shortFlag: string either "short" or "long"
    :param: contig_lengths: dictionary of headers and contig lengths
    :return: depths: dictionary of contigs and depths
    """
    depths = {}
    if shortFlag == "short":
        filename = os.path.join(out_dir, "combined_sorted.bam")
    else: # long
        filename = os.path.join(out_dir, "combined_sorted_long.bam")
    for repName, repLength in contig_lengths.items():
        depths[repName] = [0] * repLength
    depthCommand = ['samtools', 'depth', filename]
    with open(os.devnull, 'wb') as devNull:
        depthOutput = sp.check_output(depthCommand, stderr=devNull).decode()
    for line in depthOutput.splitlines(): # parse output
        parts = line.strip().split('\t')
        repName = parts[0]
        depths[repName][int(parts[1])-1] = int(parts[2])
    return depths


def collate_depths(depths, shortFlag, contig_lengths):
    """ calculates summary statistics for all depths
    :param depths:  dictionary of contigs and depths from get_depths_from_bam
    :param: shortFlag: string either "short" or "long"
    :param: contig_lengths: dictionary of headers and contig lengths
    :return: summary_df: pandas df of depth summary statistics
    """
    # define the columns of dataframe
    contig_names = []
    contig_length = []    
    mean_depth_col = []
    sd_depth_col = []
    q25_depth = []
    q75_depth = []
    # iterate over the conitgs
    for replicon_name, base_depths in depths.items():
        replicon_length = contig_lengths[replicon_name]
        unmaskedDepths = []
        for depth in enumerate(base_depths):
            unmaskedDepths.append(depth)
        try:
            mean_depth = round(statistics.mean(base_depths),2)
            depth_stdev = round(statistics.stdev(base_depths),2)
            q25, q75 = np.percentile(base_depths, [25 ,75])
            q25, q75 = int(q25), int(q75)
            # save the chromosome depth 
            if replicon_name == "chromosome":
                chromosome_depth = mean_depth
        except statistics.StatisticsError: # if can't calculate
            mean_depth, depth_stdev, q25, q75 = 'NA', 'NA', 'NA', 'NA'
        # append to list
        contig_names.append(replicon_name)
        contig_length.append(replicon_length)
        mean_depth_col.append(mean_depth)
        sd_depth_col.append(depth_stdev)
        q25_depth.append(q25)
        q75_depth.append(q75)
    # make summary df    
    if shortFlag == "short":
        summary_df = pd.DataFrame(
        {'contig': contig_names,
        'length': contig_length,
        'mean_depth_short': mean_depth_col,
        'sd_depth_short': sd_depth_col, 
        'q25_depth_short': q25_depth,
        'q75_depth_short': q75_depth
        })
        summary_df['plasmid_copy_number_short'] = round(summary_df['mean_depth_short'] / chromosome_depth,2)
    else: # long
        summary_df = pd.DataFrame(
        {'contig': contig_names,
        'mean_depth_long': mean_depth_col,
        'sd_depth_long': sd_depth_col, 
        'q25_depth_long': q25_depth,
        'q75_depth_long': q75_depth
        })
        summary_df['plasmid_copy_number_long'] = round(summary_df['mean_depth_long'] / chromosome_depth,2)
    # return df         
    return(summary_df)

def combine_depth_dfs(out_dir, df_short, df_long, prefix, circular_status):
    """ combines long and short depths
    :param out_dir:  output directory
    :param df_short: short depth summary df
    :param df_long: long depth summary df
    :param: prefix: prefix - default plassembler
    :param: circular_status: dictionary of contig header and circular status
    """
    combined_df = pd.merge(df_short, df_long, on='contig', how='outer')
        # add in circularity info 

    combined_df['circularity'] = combined_df['contig'].map(circular_status)
    out_file = os.path.join(out_dir, prefix + "_copy_number_summary.tsv")
    with open(out_file, 'w') as f:
        combined_df.to_csv(f, sep="\t", index=False, header=True)
    


