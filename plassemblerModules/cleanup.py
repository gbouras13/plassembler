import os
import subprocess as sp

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
        sp.run(["cp", os.path.join(out_dir,"unicycler_output", "assembly.fasta"), os.path.join(out_dir, prefix + "_plasmids.fasta") ])
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

