
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

def get_depth(out_dir, logger, chromosome_len, threads):
    #concatenate_chrom_plasmids(out_dir, logger)
    #processes.index_fasta(os.path.join(out_dir, "combined.fasta"),  logger)
    #bwa_map_depth_sort(out_dir, threads, logger)
    contig_lengths = get_contig_lengths(out_dir, chromosome_len)
    depths = get_depths_from_bam(out_dir, contig_lengths)
    collate_depths(depths, out_dir, chromosome_len)

def concatenate_chrom_plasmids(out_dir, logger):
    chrom_fasta =os.path.join(out_dir,"chromosome.fasta")
    plas_fasta = os.path.join(out_dir,"unicycler_output", "assembly.fasta")
    concat_file = open(os.path.join(out_dir, "combined.fasta"), "w")
    try:
        processes.concatenate_single(chrom_fasta, plas_fasta, concat_file, logger)
    except:
        sys.exit("Error with concatenate_fastas\n")  

# get lengths of contigs
def get_contig_lengths(out_dir, chromosome_len):
    contig_lengths = {}
    for dna_record in SeqIO.parse(os.path.join(out_dir, "combined.fasta"), 'fasta'):
        plas_len = len(dna_record.seq)
        dna_header = dna_record.id
        contig_lengths[dna_header] = plas_len
    return contig_lengths

def bwa_map_depth_sort(out_dir, threads, logger):
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


def get_depths_from_bam(out_dir, contig_lengths):
    depths = {}
    filename = os.path.join(out_dir, "combined_sorted.bam")
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


def collate_depths(depths, out_dir, chromosome_len):
    replicon_lengths = get_contig_lengths(out_dir, chromosome_len)
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
    summary_df = pd.DataFrame(
    {'contig': contig_names,
     'length': contig_length,
     'mean_depth': mean_depth_col,
     'sd_depth': sd_depth_col, 
     'q25_depth': q25_depth,
     'q75_depth': q75_depth
    })
    # write output
    summary_df['plasmid_copy_number'] = round(summary_df['mean_depth'] / chromosome_depth,2)
    print(summary_df)
    with open(os.path.join(out_dir, "copy_number_summary.tsv"), 'w') as f:
        summary_df.to_csv(f, sep="\t", index=False, header=True)
