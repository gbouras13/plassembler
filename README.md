plassembler
===============

Automated Bacterial Extra-chromosomal Plasmid Assembly Program
------------

plassembler is designed for automated assembly of extra-chromosomal plasmids in haploid bacterial genomes that have been sequenced with hybrid ONT (R9.4.1 or earlier) + paired-end short read sequencing.

If you are assembling a small number of bacterial genomes manually, I would recommend using Trycycler (https://github.com/rrwick/Trycycler/wiki/Generating-assemblies).

Why Does plassembler exist?
----

In long read assembled bacterial genomes, small extra-chromosomal plasmids are often difficult to assemble correctly with long read assemblers such as Flye. They often have circularisation issues and can be duplicated (see https://github.com/rrwick/Trycycler/wiki/Clustering-contigs).

plassembler was created as an automated way to ensure plasmids assemble correctly without duplicated regions.

Method
-------

1. Long reads are filtered using filtlong (https://github.com/rrwick/Filtlong) and porechop (https://github.com/rrwick/Filtlong).
2. Long-read assembly is conducted with Flye (https://github.com/fenderglass/Flye).
3. If the resulting assembly has more than 1 contig, the largest contig is checked. If it is over 90% of the length of the provided chromosome size, or is circular, then it is identified as the chromosome and extracted. All other contigs are also extracted as putative plasmid extra-chromosomal contigs.
4. Short reads are filtered using fastp (https://github.com/OpenGene/fastp).
5. Long reads are mapped to the extra-chromosomal contigs using minimap2 (https://github.com/lh3/minimap2#uguide), and short reads are mapped using bwa (https://github.com/lh3/bwa).
6. All mapped reads are kept assembled using the hybrid assembler Unicycler to generate final plasmid contigs.


Installation
------

The easiest way to install plassembler is via conda.

`conda install -c gbouras13 plassembler`

This will install all the dependencies along with plassembler.

Alternatively, the development version of plassembler can be installed manually via github.

`git clone https://github.com/gbouras13/plassembler.git`

The dependencies found in environment.yml will then need to be installed manually.

For example using conda:

```
git clone https://github.com/gbouras13/plassembler.git
cd plassembler
conda env create -f environment.yml
conda activate plassembler_env
plassembler.py -h
```

**Note for Mac Users**

plassembler should run on Linux and MacOSX machines. For MacOSX machines, the current conda installation method will only install the latest available bioconda Unicycler version of v0.4.8. plassembler should still run without any issue and provide a satisfactory assembly.
However, Ryan Wick (author of Unicycler) suggests that v0.5.0 should be used, as v0.4.8 is not compatible with the latest versions of spades (https://github.com/rrwick/Unicycler/releases/tag/v0.5.0). This will require manual installation.

To install Unicycler v0.5.0, please see the Installation section of the Unicycler github https://github.com/rrwick/Unicycler. In particular, it is recommended that you install plassembler from github:

```
git clone https://github.com/rrwick/Unicycler.git
cd Unicycler
pip3 install git+https://github.com/rrwick/Unicycler.git
```

Mac M1 users may need to change some compiler settings e.g.
```
git clone https://github.com/rrwick/Unicycler.git
cd Unicycler
python3 setup.py install --makeargs "CXX=g++"
```


Running plassembler
--------

To run plassembler

` plassembler.py -l <long read fastq> -o <output dir> -s1 < short read R1 fastq> -s2 < short read R2 fastq>  -c <estimated chromosome length>`

-c will default to 2500000 (a lower bound for the estimated genome length of staphylococcus aureus).

To specify threads:

` plassembler.py -l <long read fastq> -o <output dir> -s1 < short read R1 fastq> -s2 < short read R2 fastq>  -c <estimated chromosome length> -t <threads>`

To specify a prefix for the output files:

` plassembler.py -l <long read fastq> -o <output dir> -s1 < short read R1 fastq> -s2 < short read R2 fastq>  -c <estimated chromosome length> -t <threads> -p <prefix>`

To specify a minimum length and minimum read quality Q-score for filtlong :

` plassembler.py -l <long read fastq> -o <output dir> -s1 < short read R1 fastq> -s2 < short read R2 fastq>  -c <estimated chromosome length> -t <threads> -p <prefix> -m <min length> -q <min quality>`

-m will default to 1000 and -q will default to 9. Note that for some tiny plasmids, -m may need to be reduced (see https://www.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000631#tab2).

To overwrite an existing output directory, use -f

` plassembler.py -l <long read fastq> -o <output dir> -s1 < short read R1 fastq> -s2 < short read R2 fastq>  -c <estimated chromosome length> -t <threads> -p <prefix> -m <min length> -q <min quality> -f`

plassembler defaults to 8 threads.


Outputs
-------
plassembler will output a `_plasmids.fasta` file, which will contain the assembled plasmid sequence(s) in fasta format, and a `_plasmids.gfa` file, which will contain the assembly  graph that can be visualised in Bandage (https://github.com/rrwick/Bandage) - they should all hopefully be circular plasmids.

Acknowledgements
-------

Infinite thanks are owed to Ryan Wick (https://github.com/rrwick), who not only wrote most of the programs used in plassembler, but also gave me numerous ideas about how to approach the plasmid assembly issue.

Version Log
--------
A brief description of what is new in each update of plassembler can be found in the HISTORY.md file.

Bugs and Suggestions
--------
If you come across bugs with plassembler, or would like to make any suggestions to improve the program, please open an issue or email george.bouras@adelaide.edu.au.

Other Possible Directions
------
At the moment, plassembler is designed for users with hybrid ONT long read (R9.4.1 and earlier) and matching short read data. However, with the new Kit 14 chemistry, ONT long reads may be accurate enough that short read sequencing is not required to polish bacterial assemblies. Other approaches may be more appropriate for Kit 14 long read only assemblies - see https://twitter.com/rrwick/status/1548926644085108738?cxt=HHwWhMClvfCk8v4qAAAA - this is a work in progress.
