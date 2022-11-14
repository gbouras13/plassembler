
#########################################
##### code taken from & modified
##### https://github.com/rrwick/Small-plasmid-Nanopore/blob/main/scripts/get_depths.py
#########################################
from Bio import SeqIO
import os
import processes
import sys
import subprocess as sp
import statistics
import numpy as np
import pandas as pd

# concatenate plasmids and chromosome 

def get_depth(out_dir, logger,  threads, prefix):
    concatenate_chrom_plasmids(out_dir, logger)
    processes.index_fasta(os.path.join(out_dir, "combined.fasta"),  logger)
    bwa_map_depth_sort(out_dir, threads)
    minimap_depth_sort(out_dir, threads)
    contig_lengths = get_contig_lengths(out_dir)
    depths = get_depths_from_bam(out_dir, "short", contig_lengths)
    depths_long = get_depths_from_bam(out_dir, "long", contig_lengths)
    circular_status = get_contig_circularity(out_dir)
    summary_df_short = collate_depths(depths,"short", out_dir)
    summary_df_long = collate_depths(depths_long,"long", out_dir)
    combine_outputs(out_dir, summary_df_short, summary_df_long, prefix, circular_status)

def concatenate_chrom_plasmids(out_dir, logger):
    chrom_fasta =os.path.join(out_dir,"chromosome.fasta")
    plas_fasta = os.path.join(out_dir,"unicycler_output", "assembly.fasta")
    concat_file = open(os.path.join(out_dir, "combined.fasta"), "w")
    try:
        processes.concatenate_single(chrom_fasta, plas_fasta, concat_file, logger)
    except:
        sys.exit("Error with concatenate_fastas\n")  

# get lengths of contigs
def get_contig_lengths(out_dir):
    contig_lengths = {}
    for dna_record in SeqIO.parse(os.path.join(out_dir, "combined.fasta"), 'fasta'):
        plas_len = len(dna_record.seq)
        dna_header = dna_record.id
        contig_lengths[dna_header] = plas_len
    return contig_lengths

# get circular status of contigs
def get_contig_circularity(out_dir):
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
    trim_one = os.path.join(out_dir, "trimmed_R1.fastq")
    trim_two = os.path.join(out_dir, "trimmed_R2.fastq")
    fasta = os.path.join(out_dir, "combined.fasta")
    bam = os.path.join(out_dir, "combined_sorted.bam")
    try:
        bwa_map = sp.Popen(["bwa", "mem", "-t", threads, fasta, trim_one, trim_two ], stdout=sp.PIPE) 
        samtools_sort = sp.Popen(["samtools", "sort", "-@", threads, "-o", bam, "-" ], stdin=bwa_map.stdout ) 
        samtools_sort.communicate()[0]
    except:
        sys.exit("Error with mapping and sorting\n")  

def minimap_depth_sort(out_dir, threads):
    input_long_reads = os.path.join(out_dir, "filtered_long_reads.fastq.gz")
    fasta = os.path.join(out_dir, "combined.fasta")
    bam = os.path.join(out_dir, "combined_sorted_long.bam")
    try:
        minimap = sp.Popen(["minimap2", "-ax", "map-ont", "-t", threads, fasta, input_long_reads ], stdout=sp.PIPE) 
        samtools_sort = sp.Popen(["samtools", "sort", "-@", threads, "-o", bam, "-" ], stdin=minimap.stdout ) 
        samtools_sort.communicate()[0]
    except:
        sys.exit("Error with mapping and sorting\n")  


def get_depths_from_bam(out_dir, flag, contig_lengths):
    depths = {}
    if flag == "short":
        filename = os.path.join(out_dir, "combined_sorted.bam")
    else:
        filename = os.path.join(out_dir, "combined_sorted_long.bam")
    for rep_name, rep_length in contig_lengths.items():
        depths[rep_name] = [0] * rep_length
    depth_command = ['samtools', 'depth', filename]
    with open(os.devnull, 'wb') as dev_null:
        depth_output = sp.check_output(depth_command, stderr=dev_null).decode()
    for line in depth_output.splitlines():
        parts = line.strip().split('\t')
        rep_name = parts[0]
        depths[rep_name][int(parts[1])-1] = int(parts[2])
    return depths


def collate_depths(depths, flag, out_dir):
    replicon_lengths = get_contig_lengths(out_dir)
    # define the columns of dataframe
    contig_names = []
    contig_length = []    
    mean_depth_col = []
    sd_depth_col = []
    q25_depth = []
    q75_depth = []
    # iterate over the conitgs
    for replicon_name, base_depths in depths.items():
        replicon_length = replicon_lengths[replicon_name]
        unmasked_depths = []
        for depth in enumerate(base_depths):
            unmasked_depths.append(depth)
        try:
            mean_depth = round(statistics.mean(base_depths),2)
            depth_stdev = round(statistics.stdev(base_depths),2)
            q25, q75 = np.percentile(base_depths, [25 ,75])
            q25, q75 = int(q25), int(q75)
            # save the chromosome depth 
            if replicon_name == "chromosome":
                chromosome_depth = mean_depth
        except statistics.StatisticsError:
            mean_depth, depth_stdev, q25, q75 = 'NA', 'NA', 'NA', 'NA'
        # print(f'{replicon_name}\t{replicon_length}\t{len(base_depths)}\t{mean_depth}\t'
        #         f'{depth_stdev}\t{q25}\t{q75}\t{iqr}')
        contig_names.append(replicon_name)
        contig_length.append(replicon_length)
        mean_depth_col.append(mean_depth)
        sd_depth_col.append(depth_stdev)
        q25_depth.append(q25)
        q75_depth.append(q75)
    # make summary df    
    if flag == "short":
        summary_df = pd.DataFrame(
        {'contig': contig_names,
        'length': contig_length,
        'mean_depth_short': mean_depth_col,
        'sd_depth_short': sd_depth_col, 
        'q25_depth_short': q25_depth,
        'q75_depth_short': q75_depth
        })
        summary_df['plasmid_copy_number_short'] = round(summary_df['mean_depth_short'] / chromosome_depth,2)
    else:
        summary_df = pd.DataFrame(
        {'contig': contig_names,
        'mean_depth_long': mean_depth_col,
        'sd_depth_long': sd_depth_col, 
        'q25_depth_long': q25_depth,
        'q75_depth_long': q75_depth
        })
        summary_df['plasmid_copy_number_long'] = round(summary_df['mean_depth_long'] / chromosome_depth,2)
    # return df         
    print(summary_df)
    return(summary_df)

def combine_outputs(out_dir, df_short, df_long, prefix, circular_status):
    combined_df = pd.merge(df_short, df_long, on='contig', how='outer')
        # add in circularity info 

    combined_df['circularity'] = combined_df['contig'].map(circular_status)
    out_file = os.path.join(out_dir, prefix + "_copy_number_summary.tsv")
    with open(out_file, 'w') as f:
        combined_df.to_csv(f, sep="\t", index=False, header=True)

