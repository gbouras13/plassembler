plassembler
===============

Automated Bacterial Extra-chromosomal Plasmid Assembly Program
------------

plassembler is designed for automated assembly of extra-chromosomal plasmids in haploid bacterial genomes that have been sequenced with hybrid ONT (R9.4.1 or earlier) + paired-end short read sequencing.

If you are assembling a small number of bacterial genomes manually, I would recommend using Trycycler (https://github.com/rrwick/Trycycler/wiki/Generating-assemblies).

Generally, I would highly recommend reading the following guides to bacterial genome assembly regardless of whether you want to use plassembler (https://github.com/rrwick/Trycycler/wiki https://github.com/rrwick/Trycycler/wiki/Guide-to-bacterial-genome-assembly).

Why Does plassembler exist?
----

In long read assembled bacterial genomes, small extra-chromosomal plasmids are often difficult to assemble correctly with long read assemblers such as Flye. They often have circularisation issues and can be duplicated or missed (see https://f1000research.com/articles/8-2138 https://github.com/rrwick/Trycycler/wiki/Clustering-contigs).

plassembler was created as an automated way to ensure plasmids assemble correctly without duplicated regions for high-throughput uses.

Method
-------

1. Long reads are filtered using nanofilt (https://github.com/wdecoster/nanofilt) .
2. Long-read assembly is conducted with Flye (https://github.com/fenderglass/Flye).
3. If the resulting assembly has more than 1 contig, the largest contig is checked. If it is over 90% of the length of the provided chromosome size, or is circular, then it is identified as the chromosome and extracted. All other contigs are extracted as putative plasmid extra-chromosomal contigs.
4. Short reads are filtered using fastp (https://github.com/OpenGene/fastp).
5. Long reads are mapped to the extra-chromosomal contigs using minimap2 (https://github.com/lh3/minimap2#uguide), and short reads are mapped using bwa (https://github.com/lh3/bwa).
6. Long reads are mapped to the chromosome using minimap2 and short reads are mapped using bwa. This is done to identify reads that do not map to the chromosome (for any plasmids that Flye may have missed assembling).
7. All reads that map to the extra-chromosomal contigs and all reads that do not map the chromosome are extracted, combined and de-duplicated.
8. These are assembled using the hybrid assembler Unicycler to generate final plasmid contigs.
9. Average read coverage depth for each plasmid is calculating using a modified version of code found in (https://github.com/rrwick/Small-plasmid-Nanopore, https://www.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000631#tab2).
10. Plasmid copy number is calculated by dividing the plasmid read depth by the chromosome read depth.


Other Features (Work in Progress)

1. All reads that map to both the chromosome and plasmid are extracted and assembled (short read only assembly).


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

plassembler should run on Linux and MacOSX machines. For Linux environments, Unicycler v0.5.0 should be installed with the conda installation. For MacOSX environments, the current conda installation method will only install the latest available bioconda Unicycler version of v0.4.8. plassembler should still run without any issue and provide a satisfactory assembly.

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

* -c will default to 2500000 (a lower bound for the estimated genome length of staphylococcus aureus).

To specify threads:

` plassembler.py -l <long read fastq> -o <output dir> -s1 < short read R1 fastq> -s2 < short read R2 fastq>  -c <estimated chromosome length> -t <threads>`

To specify a prefix for the output files:

` plassembler.py -l <long read fastq> -o <output dir> -s1 < short read R1 fastq> -s2 < short read R2 fastq>  -c <estimated chromosome length> -t <threads> -p <prefix>`

To specify a minimum length and minimum read quality Q-score for nanofilt :

` plassembler.py -l <long read fastq> -o <output dir> -s1 < short read R1 fastq> -s2 < short read R2 fastq>  -c <estimated chromosome length> -t <threads> -p <prefix> -m <min length> -q <min quality>`

* -m will default to 1000 and -q will default to 9. Note that for some tiny plasmids, -m may need to be reduced (see https://www.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000631#tab2).

To overwrite an existing output directory, use -f

` plassembler.py -l <long read fastq> -o <output dir> -s1 < short read R1 fastq> -s2 < short read R2 fastq>  -c <estimated chromosome length> -t <threads> -p <prefix> -m <min length> -q <min quality> -f`

plassembler defaults to 8 threads.


Outputs
-------
plassembler will output a `_plasmids.fasta` file, which will contain the assembled plasmid sequence(s) in fasta format, and a `_plasmids.gfa` file, which will contain the assembly  graph that can be visualised in Bandage (https://github.com/rrwick/Bandage) - they should all ideally be circular plasmids.

plassembler also outputs a `copy_number_summary.tsv` files, which gives the estimated copy number for each plasmid, for both short reads and long reads (see https://www.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000631#tab2 for more details about plasmid copy numbers).

Acknowledgements
-------

Infinite thanks are owed to Ryan Wick (https://github.com/rrwick), who not only wrote Unicycler and some other code used in plassembler, but also gave me numerous ideas about how to approach the plasmid assembly issue.

Version Log
--------
A brief description of what is new in each update of plassembler can be found in the HISTORY.md file.

Bugs and Suggestions
--------
If you come across bugs with plassembler, or would like to make any suggestions to improve the program, please open an issue or email george.bouras@adelaide.edu.au.

Other Future Directions
------
At the moment, plassembler is designed for users with hybrid ONT long read (R9.4.1 and earlier) and matching short read data. However, with the new Kit 14 chemistry, ONT long reads may be accurate enough that short read sequencing is not required to polish bacterial assemblies. Other approaches may be more appropriate for Kit 14 long read only assemblies - see https://twitter.com/rrwick/status/1548926644085108738?cxt=HHwWhMClvfCk8v4qAAAA - this is not supported at the moment but may be in the future.

Citations
-----

If you use plassembler, please cite:

* Kolmogorov, M., Yuan, J., Lin, Y. et al. Assembly of long, error-prone reads using repeat graphs. Nat Biotechnol 37, 540–546 (2019). https://doi.org/10.1038/s41587-019-0072-8
* Li H., Minimap2: pairwise alignment for nucleotide sequences, Bioinformatics, Volume 34, Issue 18 Pages 3094–3100 (2018), https://doi.org/10.1093/bioinformatics/bty191
* Li H., Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM. arXiv preprint arXiv:1303.3997 (2013).
* Wick RR, Judd LM, Gorrie CL, Holt KE Unicycler: Resolving bacterial genome assemblies from short and long sequencing reads. PLoS Comput Biol 13(6): e1005595 (2017). https://doi.org/10.1371/journal.pcbi.1005595
* Heng Li, Bob Handsaker, Alec Wysoker, Tim Fennell, Jue Ruan, Nils Homer, Gabor Marth, Goncalo Abecasis, Richard Durbin, 1000 Genome Project Data Processing Subgroup, The Sequence Alignment/Map format and SAMtools, Bioinformatics, Volume 25, Issue 16, 15 August 2009, Pages 2078–2079, https://doi.org/10.1093/bioinformatics/btp352
* Wick RR, Judd LM, Wyres KL, Holt KE. Recovery of small plasmid sequences via Oxford Nanopore sequencing. Microb Genom. 2021 Aug;7(8):000631. doi: 10.1099/mgen.0.000631. PMID: 34431763; PMCID: PMC8549360.
