import os
import subprocess as sp
from Bio import SeqIO
import pandas as pd
import sys
from Bio.SeqRecord import SeqRecord

#######################################################
# rename output files to add new depth
#####################################################



def rename_contigs(out_dir, prefix):
    """
    Renames the contigs of unicycler with the new plasmid copy numbers
    :param out_dir: output directory
    :return: 
    """

    depth_df = pd.read_csv(os.path.join(out_dir, prefix + "_copy_number_summary.tsv"), delimiter= '\t', index_col=False, header=0 ) 
    depth_df = depth_df.loc[depth_df['contig'] != 'chromosome'].reset_index(drop=True)
    # get contigs only
    plasmid_fasta = os.path.join(out_dir,"unicycler_output", "assembly.fasta")
    i = 0
    with open(os.path.join(out_dir, prefix + "_plasmids.fasta"), 'w') as dna_fa:
        for dna_record in SeqIO.parse(plasmid_fasta, 'fasta'): 
            if "circular" in dna_record.description:
                id_updated = dna_record.description.split(' ')[0] + " " + dna_record.description.split(' ')[1] + " plasmid_copy_number_short=" + str(depth_df.plasmid_copy_number_short[i]) + "x plasmid_copy_number_long=" + str(depth_df.plasmid_copy_number_long[i]) + "x " + dna_record.description.split(' ')[3]
            else:
                id_updated = dna_record.description.split(' ')[0] + " " + dna_record.description.split(' ')[1] + " plasmid_copy_number_short=" + str(depth_df.plasmid_copy_number_short[i]) + "x plasmid_copy_number_long=" + str(depth_df.plasmid_copy_number_long[i]) + "x " 
            i += 1
            record = SeqRecord(dna_record.seq, id=id_updated, description = "" )
            SeqIO.write(record, dna_fa, 'fasta')

def rename_contigs_kmer(out_dir, prefix):
    """
    Renames the contigs of unicycler with the new plasmid copy numbers kmer mode
    :param out_dir: output directory
    :return: 
    """

    depth_df = pd.read_csv(os.path.join(out_dir, prefix + "_copy_number_summary.tsv"), delimiter= '\t', index_col=False, header=0 ) 
    depth_df = depth_df.loc[depth_df['contig'] != 'chromosome'].reset_index(drop=True)
    # get contigs only
    plasmid_fasta = os.path.join(out_dir,"unicycler_output", "assembly.fasta")
    i = 0
    with open(os.path.join(out_dir, prefix + "_plasmids.fasta"), 'w') as dna_fa:
        for dna_record in SeqIO.parse(plasmid_fasta, 'fasta'): 
            if "circular" in dna_record.description:
                id_updated = dna_record.description.split(' ')[0] + " " + dna_record.description.split(' ')[1] + " plasmid_copy_number_long=" + str(depth_df.plasmid_copy_number_long[i]) + "x " + dna_record.description.split(' ')[3]
            else:
                id_updated = dna_record.description.split(' ')[0] + " " + dna_record.description.split(' ')[1] + " plasmid_copy_number_long=" + str(depth_df.plasmid_copy_number_long[i]) + "x " 
            i += 1
            record = SeqRecord(dna_record.seq, id=id_updated, description = "" )
            SeqIO.write(record, dna_fa, 'fasta')




#######################################################
# add PLSDB hit to copy_number_summary
#####################################################

def update_copy_number_summary_plsdb(out_dir, prefix, mash_empty):
    """
    Updates copy number summary
    :param out_dir: output directory
    :return: 
    """
    depth_df = pd.read_csv(os.path.join(out_dir, prefix + "_copy_number_summary.tsv"), delimiter= '\t', index_col=False, header=0 ) 
    
    if mash_empty == False:

        mash_df = pd.read_csv(os.path.join(out_dir, prefix + "_top_hits_mash_plsdb.tsv"), delimiter= '\t', index_col=False, header=0 )  

        mash_df_reduced =  mash_df[['contig']].copy()
        mash_df_reduced['plsdb_hit'] = 'Yes'

        mash_df_reduced['contig']=mash_df_reduced['contig'].astype(str)
        depth_df['contig']=depth_df['contig'].astype(str)

        combined_df = depth_df.merge(mash_df_reduced, on='contig', how='left')
        combined_df['plsdb_hit'] = combined_df['plsdb_hit'].fillna("No")

        # overwrite the file
        out_file = os.path.join(out_dir, prefix + "_copy_number_summary.tsv")
        with open(out_file, 'w') as f:
            combined_df.to_csv(f, sep="\t", index=False, header=True)
    # empty mash - update with 
    else:
        depth_df['plsdb_hit'] = 'No'
            # overwrite the file
        out_file = os.path.join(out_dir, prefix + "_copy_number_summary.tsv")
        with open(out_file, 'w') as f:
            depth_df.to_csv(f, sep="\t", index=False, header=True)








####################################################
# cleanup
##########################################################

def remove_intermediate_files(out_dir):
    """ removes intermediate files
    :param out_dir:  Output Directory
    :return: 
    """
    sp.run(["rm -rf "+ os.path.join(out_dir,"*.fastq") ], shell=True)
    sp.run(["rm -rf "+ os.path.join(out_dir,"*.fastq.gz") ], shell=True)
    sp.run(["rm -rf "+ os.path.join(out_dir,"*.bam") ], shell=True)
    sp.run(["rm -rf "+ os.path.join(out_dir,"*.sa") ], shell=True)
    sp.run(["rm -rf "+ os.path.join(out_dir,"*.sam") ], shell=True)
    sp.run(["rm -rf "+ os.path.join(out_dir,"*.amb") ], shell=True)
    sp.run(["rm -rf "+ os.path.join(out_dir,"*.ann") ], shell=True)
    sp.run(["rm -rf "+ os.path.join(out_dir,"*.pac") ], shell=True)
    sp.run(["rm -rf "+ os.path.join(out_dir,"*.bwt") ], shell=True)
    sp.run(["rm", "-rf", os.path.join(out_dir,"00-assembly") ])
    sp.run(["rm", "-rf", os.path.join(out_dir,"10-consensus") ])
    sp.run(["rm", "-rf", os.path.join(out_dir,"20-repeat") ])
    sp.run(["rm", "-rf", os.path.join(out_dir,"30-contigger") ])
    sp.run(["rm", "-rf", os.path.join(out_dir,"40-polishing") ])
    sp.run(["rm", "-rf", os.path.join(out_dir,"params.json") ])
    # delete flye assemble files
    sp.run(["rm", "-rf", os.path.join(out_dir,"chromosome.fasta") ])
    sp.run(["rm", "-rf", os.path.join(out_dir,"non_chromosome.fasta") ])
    sp.run(["rm", "-rf", os.path.join(out_dir,"combined.fasta") ])
    # delete mash 
    sp.run(["rm", "-rf", os.path.join(out_dir,"mash.tsv") ])



def move_and_copy_files(out_dir, prefix, unicycler_success_flag):
    """ moves and copies files
    :param out_dir:  Output Directory
    :param prefix: prefix
    :param unicycler_success_flag: whether or not unicycler worked
    :return: 
    """
    # move flye output into dir
    sp.run(["mkdir", "-p", os.path.join(out_dir,"flye_output") ])
    sp.run(["mv",  os.path.join(out_dir,"assembly.fasta"), os.path.join(out_dir,"flye_output") ])
    sp.run(["mv",  os.path.join(out_dir,"assembly_info.txt"), os.path.join(out_dir,"flye_output") ])
    sp.run(["mv",  os.path.join(out_dir,"flye.log"), os.path.join(out_dir,"flye_output") ])
    sp.run(["mv", os.path.join(out_dir,"assembly_graph.gfa"), os.path.join(out_dir,"flye_output") ])
    sp.run(["mv", os.path.join(out_dir,"assembly_graph.gv"), os.path.join(out_dir,"flye_output") ])
    if unicycler_success_flag == True:
         # move unicycler output to main directory
        sp.run(["cp", os.path.join(out_dir,"unicycler_output", "assembly.gfa"), os.path.join(out_dir, prefix + "_plasmids.gfa") ])
    else:
        # to touch empty versions of the output files if no plasmids 
        touch_output_fail_files(out_dir, prefix)


# function to touch create a file 
# https://stackoverflow.com/questions/12654772/create-empty-file-using-python
def touch_file(path):
    with open(path, 'a'):
        os.utime(path, None)

# to create empty plasmids fasta and gfa files
def touch_output_fail_files(out_dir, prefix):
    touch_file(os.path.join(out_dir, prefix + "_plasmids.fasta"))
    touch_file(os.path.join(out_dir, prefix + "_plasmids.gfa"))
    touch_file(os.path.join(out_dir, prefix + "_copy_number_summary.tsv"))
    touch_file(os.path.join(out_dir, prefix + "_top_hits_mash_plsdb.tsv"))

