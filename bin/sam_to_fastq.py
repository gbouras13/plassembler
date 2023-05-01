import pysam
import os
from collections import defaultdict


def extract_bin_long_fastqs(out_dir, multi_map):

    #################################################
    #### Define file paths
    #################################################

    sam_name = os.path.join(out_dir, "long_read.sam")

    # reads mapping to plasmids, or not mapping to any contigs
    plasmidfile = open(os.path.join(out_dir, "plasmid_long.fastq"), 'w')

    # Open a FASTQ file for writing reads mapping to multiple contigs
    if multi_map == True:
        multimap_plasmid_chromosome_fastqfile = open(os.path.join(out_dir, "multimap_plasmid_chromosome_long.fastq"), 'w')

    # chromosome fastqs
    chrom_fastqfile = open(os.path.join(out_dir, "chromosome_mapped_long.fastq"), 'w')

    #################################################
    #### Get the single and multiple map reads as lists
    #################################################

    # get names, single and multi as lists
    read_names = []
    single_read_names = []
    multi_read_names = []

    # open samfile
    samfile = pysam.AlignmentFile(sam_name, 'r')

    # get list of all read names
    for read in samfile.fetch():
      read_names.append(read.query_name)

    # Create a defaultdict with int as the default factory
    count_dict = defaultdict(int)

    # Loop through the list and count occurrences
    for item in read_names:
        count_dict[item] += 1

    # Get the counts
    for key, value in count_dict.items():
      if value == 1:
        single_read_names.append(key)
      else:
        multi_read_names.append(key)

    samfile.close()

    #################################################
    #### process all single reads and then get counts of plasmid vs chromosome ####
    #################################################

    samfile = pysam.AlignmentFile(sam_name, 'r')

    # Create a defaultdict with int as the default factory for the multimap reads
    plasmid_mm_dict = defaultdict(int)
    chromosome_mm_dict = defaultdict(int)

    for read in samfile.fetch():
      # Access the read's name, sequence, quality scores, etc.
      read_name = read.query_name
      sequence = read.query_sequence
      quality = read.query_qualities
      flag = read.flag
      # get contig name for the read
      contig_name = samfile.get_reference_name(read.reference_id)

      # print(read_name)
      # print(read.reference_id)
      # print(read.next_reference_id)
      # print(contig_name)
      # print(flag)

      # single reads - easy :)
      if read_name in single_read_names: 
          # gets all reads that plasmid mapped reads and all unmapped reads
          if (contig_name and 'plasmid' in contig_name) or read.is_unmapped: 
            # Write the read to the FASTQ file
            plasmidfile.write(f'@{read_name}\n')  # Write read name
            plasmidfile.write(f'{sequence}\n')    # Write sequence
            plasmidfile.write('+{0}\n'.format(read_name))  # Write quality header
            plasmidfile.write(''.join(chr(q + 33) for q in quality) + '\n')  # Write quality scores in ASCII format
          elif contig_name and 'chromosome' in contig_name:
            # Write the read to the unmapped reads FASTQ file
            chrom_fastqfile.write(f'@{read_name}\n')  # Write read name
            chrom_fastqfile.write(f'{sequence}\n')    # Write sequence
            chrom_fastqfile.write('+{0}\n'.format(read_name))  # Write quality header
            chrom_fastqfile.write(''.join(chr(q + 33) for q in quality) + '\n')  # Write quality scores in ASCII format  
    # create count dictionaries for multimap reads next step
      else: 
        if contig_name and 'plasmid' in contig_name:
          plasmid_mm_dict[read_name] += 1
        elif contig_name and 'chromosome' in contig_name:
          chromosome_mm_dict[read_name] += 1

    samfile.close()

    #################################################
    #### process all multimap reads
    #################################################

    samfile = pysam.AlignmentFile(sam_name, 'r')

    for read in samfile.fetch():

        read_name = read.query_name
        sequence = read.query_sequence
        quality = read.query_qualities
        flag = read.flag

      # multireads
        if read_name in multi_read_names: 
            if plasmid_mm_dict[read_name] > 0 and chromosome_mm_dict[read_name] > 0: # multimap both plasmid and chromosome
                if quality is not None and (flag == 0 or flag == 16):  # get only the primary
                    if multi_map == True: # only write to file if specified
                        multimap_plasmid_chromosome_fastqfile.write(f'@{read_name}\n')  # Write read name
                        multimap_plasmid_chromosome_fastqfile.write(f'{sequence}\n')    # Write sequence
                        multimap_plasmid_chromosome_fastqfile.write('+{0}\n'.format(read_name))  # Write quality header
                        multimap_plasmid_chromosome_fastqfile.write(''.join(chr(q + 33) for q in quality) + '\n')  # Write quality scores in ASCII format
            # write all that map to plasmid to the plasmid file 
            elif plasmid_mm_dict[read_name] > 0: # multimap plasmid 
                if quality is not None and (flag == 0 or flag == 16): # get only the primary
                    plasmidfile.write(f'@{read_name}\n')  # Write read name
                    plasmidfile.write(f'{sequence}\n')    # Write sequence
                    plasmidfile.write('+{0}\n'.format(read_name))  # Write quality header
                    plasmidfile.write(''.join(chr(q + 33) for q in quality) + '\n')  # Write quality scores in ASCII format
                    # write all that map to chromosome to the plasmid file 
            elif chromosome_mm_dict[read_name] > 0: # multimap chromosome 
                if quality is not None and (flag == 0 or flag == 16): # get only the primary
                    chrom_fastqfile.write(f'@{read_name}\n')  # Write read name
                    chrom_fastqfile.write(f'{sequence}\n')    # Write sequence
                    chrom_fastqfile.write('+{0}\n'.format(read_name))  # Write quality header
                    chrom_fastqfile.write(''.join(chr(q + 33) for q in quality) + '\n')  # Write quality scores in ASCII format

    # Close the FASTQ file
    plasmidfile.close()
    chrom_fastqfile.close()
    if multi_map == True:
        multimap_plasmid_chromosome_fastqfile.close()

    # Close the SAM file
    samfile.close()


#### thorough but too slow - no multithreading - at least a good example of what I want
#### maybe rewrite in rust when I get some time
#### only activate if multi_map is true

def extract_bin_short_fastqs(out_dir):

     #################################################
     #### Define file paths
     #################################################

     sam_name = os.path.join(out_dir, "short_read.sam")

     # reads mapping to plasmids, or not mapping to any contigs
     plasmidfile_r1 = open(os.path.join(out_dir, "short_read_concat_R1.fastq"), 'w')
     plasmidfile_r2 = open(os.path.join(out_dir, "short_read_concat_R2.fastq"), 'w')

     # Open a FASTQ file for writing reads mapping to multiple contigs
     multimap_plasmid_chromosome_fastqfile_r1 = open(os.path.join(out_dir, "multimap_plasmid_chromosome_R1.fastq"), 'w')
     multimap_plasmid_chromosome_fastqfile_r2 = open(os.path.join(out_dir, "multimap_plasmid_chromosome_R2.fastq"), 'w')

     # chromosome fastqs
     chrom_fastqfile_r1 = open(os.path.join(out_dir, "chromosome_mapped_R1.fastq"), 'w')
     chrom_fastqfile_r2 = open(os.path.join(out_dir, "chromosome_mapped_R2.fastq"), 'w')

     # singltones 
     singletons_fastqfile = open("singletons_short.fastq", "w")


     #################################################
     #### Get the single and multiple map reads as lists
     #################################################

     # get names, single and multi as lists
     read_names = []
     singletons_read_names =[]
     single_read_names = []
     multi_read_names = []
     unmapped_read_names =[]

     # open samfile
     samfile = pysam.AlignmentFile(sam_name, 'r')

     # get list of all read names
     for read in samfile.fetch():
         if not read.is_unmapped:
             if read.is_paired:
                 read_names.append(read.query_name)
             else:
                 singletons_read_names.append(read.query_name)
         else:
             unmapped_read_names.append(read.query_name)


     # Create a defaultdict with int as the default factory
     count_dict = defaultdict(int)

     # Loop through the list and count occurrences
     for item in read_names:
         count_dict[item] += 1

     # Get the counts
     for key, value in count_dict.items():
         if value == 2:
             single_read_names.append(key)
         else:
             multi_read_names.append(key)

     print(multi_read_names)
     samfile.close()

     #################################################
     #### process all single reads and then get counts of plasmid vs chromosome ####
     #################################################

     samfile = pysam.AlignmentFile(sam_name, 'r')

     # Create a defaultdict with int as the default factory for the multimap reads
     plasmid_mm_dict = defaultdict(int)
     chromosome_mm_dict = defaultdict(int)

     for read in samfile.fetch():
       # Access the read's name, sequence, quality scores, etc.
       read_name = read.query_name
       sequence = read.query_sequence
       quality = read.query_qualities
       flag = read.flag
       # get contig name for the read
       contig_name = samfile.get_reference_name(read.reference_id)

       # print(read_name)
       # print(read.reference_id)
       # print(read.next_reference_id)
       # print(contig_name)
       # print(flag)

       # single reads - easy :)
      # keep all unmapped and single mapped reads
       if (read_name in single_read_names) or (read_name in unmapped_read_names): 
           # gets all paired plasmid mapped reads and all paired unmapped reads
           if (contig_name and 'plasmid' in contig_name) or read.is_unmapped: 
             if read.is_read1:
                 plasmidfile_r1.write("@{}\n{}\n+\n{}\n".format(
                     read_name, read_name, quality))
             # Extract read 2 and write to tp_R2.fastq
             elif read.is_read2:
                 plasmidfile_r2.write("@{}\n{}\n+\n{}\n".format(
                     read_name, read_name, quality))
           # gets all paired chromosome mapped reads 
           elif contig_name and 'chromosome' in contig_name:
             # Write the read to the unmapped reads FASTQ file
             if read.is_read1:
                 chrom_fastqfile_r1.write("@{}\n{}\n+\n{}\n".format(
                     read_name, read_name, quality))
             # Extract read 2 and write to tp_R2.fastq
             elif read.is_read2:
                 chrom_fastqfile_r2.write("@{}\n{}\n+\n{}\n".format(
                     read_name, read_name, quality))
     # create count dictionaries for multimap reads next step
       elif read_name in multi_read_names: 
         if contig_name and 'plasmid' in contig_name:
           plasmid_mm_dict[read_name] += 1
         elif contig_name and 'chromosome' in contig_name:
           chromosome_mm_dict[read_name] += 1
     # singletons
       elif read_name in singletons_read_names: 
         singletons_fastqfile.write("@{}\n{}\n+\n{}\n".format(
                 read_name, sequence, quality))

     samfile.close()

     #################################################
     #### process all multimap reads
     #################################################

     samfile = pysam.AlignmentFile(sam_name, 'r')

     for read in samfile.fetch():

         read_name = read.query_name
         sequence = read.query_sequence
         quality = read.query_qualities
         flag = read.flag

         if read_name in multi_read_names: 
             if plasmid_mm_dict[read_name] > 0 and chromosome_mm_dict[read_name] > 0: # multimap both plasmid and chromosome
                 if read.is_read1 and (flag == 83 or flag == 99): # primary R1s mapped in proper pair
                     multimap_plasmid_chromosome_fastqfile_r1.write("@{}\n{}\n+\n{}\n".format(
                         read_name, read_name, quality))
                 elif read.is_read1 and (flag == 147 or flag == 163): # primary R1s mapped in proper pair
                     multimap_plasmid_chromosome_fastqfile_r2.write("@{}\n{}\n+\n{}\n".format(
                         read_name, read_name, quality))
                     # write all that map to plasmid to the plasmid file 
             elif plasmid_mm_dict[read_name] > 0 : # multimap plasmid 
                 if read.is_read1 and (flag == 83 or flag == 99): # primary R1s mapped in proper pair
                     plasmidfile_r1.write("@{}\n{}\n+\n{}\n".format(
                         read_name, read_name, quality))
                 elif read.is_read1 and (flag == 147 or flag == 163): # primary R1s mapped in proper pair
                     plasmidfile_r2.write("@{}\n{}\n+\n{}\n".format(
                         read_name, read_name, quality))
                 # write all that map to chromosome to the chromosome file 
             elif chromosome_mm_dict[read_name] > 0: # multimap chromosome 
                 if read.is_read1 and (flag == 83 or flag == 99): # primary R1s mapped in proper pair
                     chrom_fastqfile_r1.write("@{}\n{}\n+\n{}\n".format(
                         read_name, read_name, quality))
                 elif read.is_read2 and (flag == 147 or flag == 163): # primary R1s mapped in proper pair
                     chrom_fastqfile_r2.write("@{}\n{}\n+\n{}\n".format(
                         read_name, read_name, quality))

     # Close the FASTQ file
     plasmidfile_r1.close()
     plasmidfile_r2.close()
     chrom_fastqfile_r1.close()
     chrom_fastqfile_r2.close()
     plasmidfile_r1.close()
     plasmidfile_r2.close()
     multimap_plasmid_chromosome_fastqfile_r1.close()
     multimap_plasmid_chromosome_fastqfile_r2.close()
     singletons_fastqfile.close()

     # Close the SAM file
     samfile.close()