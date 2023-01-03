plassembler
===============

Automated Bacterial Plasmid Assembly Program
------------

plassembler is designed for automated assembly of plasmids in haploid bacterial genomes that have been hybrid sequenced with Oxford Nanopore Technologies long read (ideally R9.4.1, also should work with R 10.4.1 or earlier chemistries) & paired-end short read sequencing.

If you are assembling a small number of bacterial genomes manually, I would highly recommend starting by using [Trycycler](https://github.com/rrwick/Trycycler). If you have more genomes or want to assemble your genomes in a more automated way, try [dragonflye](https://github.com/rpetit3/dragonflye), especially if you are used to Shovill, or my own pipeline [hybracter](https://github.com/gbouras13/hybracter) that is more targeted at large datasets.  

Additionally, I would recommend reading the following guides to bacterial genome assembly regardless of whether you want to use plassembler:
*  [Trycycler](https://github.com/rrwick/Trycycler/wiki/Guide-to-bacterial-genome-assembly)
*  [Perfect Bacterial Assembly Tutorial](https://github.com/rrwick/Perfect-bacterial-genome-tutorial)
*  [Perfect bacterial assembly preprint](https://preprints.scielo.org/index.php/scielo/preprint/view/5053)

Table of Contents
-----------
- [plassembler](#plassembler)
  - [Automated Bacterial Plasmid Assembly Program](#automated-bacterial-plasmid-assembly-program)
  - [Table of Contents](#table-of-contents)
- [Quick Start](#quick-start)
- [Why Does plassembler exist?](#why-does-plassembler-exist)
- [Why Not Use Unicycler?](#why-not-use-unicycler)
- [Documentation](#documentation)
- [Method](#method)
- [Installation](#installation)
  - [Unicycler v0.5.0 Installation Issues](#unicycler-v050-installation-issues)
- [Running plassembler](#running-plassembler)
- [Outputs](#outputs)
- [Acknowledgements](#acknowledgements)
- [Version Log](#version-log)
- [Bugs and Suggestions](#bugs-and-suggestions)
- [Other Future Directions](#other-future-directions)
- [Citations](#citations)

# Quick Start

The easiest way to install plassembler is via conda:

`conda install -c bioconda plassembler`

Followed by database download and installation:

`install_database.py -o <path/to/databse_dir>`

And finally assembly:

` plassembler.py -d <database directory> -l <long read fastq> -o <output dir> -1 < short read R1 fastq> -2 < short read R2 fastq>  -c <estimated chromosome length>`

Please read below for more details, especially if you are an inexperienced command line user.

# Why Does plassembler exist?

In long-read first assembled bacterial genomes, small plasmids are often difficult to assemble correctly with long read assemblers such as Flye. They often have circularisation issues and can be duplicated or missed (see [this](https://f1000research.com/articles/8-2138) and [this](https://github.com/rrwick/Trycycler/wiki/Clustering-contigs) ).

plassembler was created as an automated tool to ensure plasmids assemble correctly without duplicated regions for high-throughput uses - and to provide some useful statistics as well (like copy number for long and short read sets). Plassembler will likely also recover small plasmids that long read assemblers like Flye simply miss.

As of v 0.1.4, plassembler also uses [mash](https://github.com/marbl/Mash) as a quick way to determine whether each assembled contig has any similar hits in [PLSDB](https://doi.org/10.1093/nar/gkab1111). 

# Why Not Use Unicycler?


Unicycler is awesome and still probably the best way to assemble plasmids from hybrid sequencing - plassembler uses it! But there are a few reasons to use plassembler instead:

1. Time. plassember throws away all the chromosomal reads (i.e. most of them) before running Unicycler, so it is much faster. 
2. Plassembler will output only the likely plasmids, which may be integrated into pipelines. You shouldn't be assembling the chromosome using Unicycler [anymore](https://preprints.scielo.org/index.php/scielo/preprint/view/5053) so plassembler can get you only what is necessary from Unicycler.
3. Plassembler will give you summary coverage stats for both long and short reads.
4. Plassembler can be used as fast-ish quality control to check if your short and long reads come from the same sample - if plassembler results in many non-circular contigs (particularly those that have no hits in PLSDB), it is likely because your read sets do not come from the same isolate! 
5. As of v 0.1.4, you will get information whether each assembled contig has a similar entry in [PLSDB](https://doi.org/10.1093/nar/gkab1111). Especially for common pathogen species that are well represented in databases, this will likely tell you specifically what plasmid you have in your sample. Additionally, if there are many contigs with no PLSDB hits, this can help you determine that your long and short reads do not come from the same isolate, or that there may be some other biological phenomenon in your "isolate" (for example, perhaps there is some other non-plasmid mobile genetic element present in low abundance). For less commonly sequenced species, I would not suggest that that absence of a PLSDB hit is necessary meaningful, especially for circular contigs - those would likely be novel plasmids uncaptured by PLSDB.

# Documentation

Documentation can be found at http://plassembler.readthedocs.io/.

# Method

1. Long reads are filtered using [nanofilt](https://github.com/wdecoster/nanofilt) .
2. Long-read only assembly is conducted with [Flye](https://github.com/fenderglass/Flye).
3. If the resulting assembly is checked. If the largest contig is over 90% of the length of the provided chromosome size -c, then it is identified as the chromosome and extracted. Any other contigs are extracted as putative plasmid contigs, if Flye assembled any. Otherwise, plassembler will exit - you probably need to get some more long reads to complete your assembly (or check -c).
4. Short reads are filtered using [fastp](https://github.com/OpenGene/fastp).
5. Long reads are mapped to any putative plasmid contigs using [minimap2](https://github.com/lh3/minimap2#uguide), and short reads are mapped using [bwa](https://github.com/lh3/bwa).
6. Long reads are mapped to the chromosome using minimap2 and short reads are mapped using bwa. This is done to identify reads that do not map to the chromosome (for any plasmids that Flye may have missed assembling).
7. All reads that map to the putative plasmid contigs and all reads that do not map the chromosome are extracted, combined and de-duplicated using [seqkit](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0163962).
8. The de-deduplicated read sets are assembled using the hybrid assembler [Unicycler](https://github.com/rrwick/Unicycler) to generate final plasmid contigs.
9. Average read coverage depth for each plasmid is calculated using a modified version of code found [here](https://github.com/rrwick/Small-plasmid-Nanopore). See also this [paper](https://www.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000631#tab2).
10. Plasmid copy number is calculated by dividing the plasmid read depth by the chromosome read depth.
11. All plasmid contigs are compared against [PLSDB](https://doi.org/10.1093/nar/gkab1111) using [mash](https://github.com/marbl/Mash) with a cutoff maximum mash distance of 0.1.

Other Features (Work in Progress)

1. All reads that map to both the chromosome and plasmid are extracted and assembled (short read only assembly). This may be useful to identify possible insertion sequences and transposases that are shared between plasmid and chromosome.


# Installation

Plassembler is on bioconda. 

plassembler should run and has been tested on Linux and MacOSX machines. 

The easiest way to install plassembler is via conda.

`conda install -c bioconda plassembler`

or mamba for quicker solving:

`mamba install -c bioconda plassembler`

This will install all the dependencies along with plassembler.

Alternatively, the development version of plassembler can be installed manually via github - it may contain untested changes.

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

Unicycler v0.5.0 Installation Issues
------


For Linux environments, Unicycler v0.5.0 should be installed with the plassembler bioconda installation.

You can force it as follows:

`conda install -c bioconda plassembler unicycler==0.5.0`

or manually install Unicycler v0.5.0 after installing plassembler:

```
conda install -c bioconda plassembler
pip3 install git+https://github.com/rrwick/Unicycler.git
```

**MacOS**

For MacOS environments, the current conda installation method will only install the latest available bioconda Unicycler version of v0.4.8. plassembler should still run without any issue and provide a satisfactory assembly, but you will be warned of this when you run plassembler.

Ryan Wick (the author of Unicycler) suggests that v0.5.0 should be used, as v0.4.8 is not compatible with the latest versions of spades (see [here](https://github.com/rrwick/Unicycler/releases/tag/v0.5.0) ). This will require another installation step on MacOS.

To install Unicycler v0.5.0, it is recommended that you install Unicycler from github after installing Plassembler follows:

```
conda create -n plassemblerENV
conda activate plassemblerENV
conda install -c bioconda plassembler
pip3 install git+https://github.com/rrwick/Unicycler.git
```

Mac M1 users may need to change some compiler settings and install from the Unicycler github repo e.g.

```
conda create -n plassemblerENV
conda activate plassemblerENV
conda install -c bioconda plassembler
git clone https://github.com/rrwick/Unicycler.git
cd Unicycler
python3 setup.py install --makeargs "CXX=g++"
```

# Running plassembler


To run plassembler, first you need to install the database in a directory of your chosing:

`install_database.py -d <database directory>`

Once this is finished, you can run plassembler as followed

` plassembler.py -d <database directory> -l <long read fastq> -o <output dir> -1 < short read R1 fastq> -2 < short read R2 fastq>  -c <estimated chromosome length>`

* -c will default to 2500000 (a lower bound for the estimated genome length of Staphylococcus aureus).

To specify more threads:

` plassembler.py -d <database directory> -l <long read fastq> -o <output dir> -1 < short read R1 fastq> -2 < short read R2 fastq>  -c <estimated chromosome length> -t <threads>`

To specify a prefix for the output files:

` plassembler.py -d <database directory> -l <long read fastq> -o <output dir> -1 < short read R1 fastq> -2 < short read R2 fastq>  -c <estimated chromosome length> -t <threads> -p <prefix>`

To specify a minimum length and minimum read quality Q-score for nanofilt :

` plassembler.py -d <database directory> -l <long read fastq> -o <output dir> -1 < short read R1 fastq> -2 < short read R2 fastq>  -c <estimated chromosome length> -t <threads> -p <prefix> -m <min length> -q <min quality>`

* -m will default to 500 and -q will default to 8. Note that for some tiny plasmids, -m may need to be reduced (see this [paper](https://www.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000631#tab2) ).

To overwrite an existing output directory, use -f

` plassembler.py -d <database directory> -l <long read fastq> -o <output dir> -1 < short read R1 fastq> -2 < short read R2 fastq>  -c <estimated chromosome length> -t <threads> -p <prefix> -m <min length> -q <min quality> -f`

plassembler defaults to 1 thread.

```
usage: plassembler.py [-h] -d DATABASE -l LONGREADS -1 SHORT_ONE -2 SHORT_TWO [-c CHROMOSOME] [-o OUTDIR] [-m MIN_LENGTH]
                      [-t THREADS] [-f] [-r] [-p PREFIX] [-q MIN_QUALITY] [-V]

plassembler: accurate extra-chromosomal plasmid assembler pipeline for haploid bacterial genomes.

optional arguments:
  -h, --help            show this help message and exit
  -d DATABASE, --database DATABASE
                        Directory of PLSDB database downloaded using install_database.py.
  -l LONGREADS, --longreads LONGREADS
                        Fastq File of ONT Long Reads. Required
  -1 SHORT_ONE, --short_one SHORT_ONE
                        R1 short read fastq file. Required.
  -2 SHORT_TWO, --short_two SHORT_TWO
                        R2 short read fastq file. Required.
  -c CHROMOSOME, --chromosome CHROMOSOME
                        Approximate chromosome length of bacteria. Defaults to 2500000.
  -o OUTDIR, --outdir OUTDIR
                        Directory to write the output to. Defaults to output/
  -m MIN_LENGTH, --min_length MIN_LENGTH
                        minimum length for long reads for nanofilt. Defaults to 500.
  -t THREADS, --threads THREADS
                        Number of threads for flye and unicycler. Defaults to 1.
  -f, --force           Overwrites the output directory.
  -r, --raw_flag        Use --nano-raw for Flye Guppy FAST reads. 
                        By default, Flye will assume SUP or HAC reads and use --nano-hq
  -p PREFIX, --prefix PREFIX
                        Prefix for output files. This is not required
  -q MIN_QUALITY, --min_quality MIN_QUALITY
                        minimum quality of long reads for nanofilt. Defaults to 9.
  -V, --version         show plassembler version and exit.

```


# Outputs

plassembler will output a `_plasmids.fasta` file, which will contain the assembled plasmid sequence(s) in FASTA format (including long and short read copy numbers in the header), and a `_plasmids.gfa` file, which will contain the assembly graph from Unicycler that can be visualised in [Bandage](https://github.com/rrwick/Bandage). 

plassembler also outputs a `_copy_number_summary.tsv` file, which gives the estimated copy number for each plasmid, for both short reads and long reads (see this [paper](https://www.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000631#tab2) for more details about plasmid copy numbers) and also a `_top_hits_mash_plsdb.tsv` file, which gives each contig's top hit by mash distance in the PLSDB (if there is a hit), along with all its supporting information. 

If there is no PLSDB hit for a contig, there will be no entry in the `_top_hits_mash_plsdb.tsv` file.

If plassembler fails to assemble any plasmids at all in `_plasmids.fasta`, all these files will still exist, but will be empty (to ensure plassembler can be easily integrated into workflow managers like Snakemake).

plassembler will also output a log file, a `flye_output` directory, which contains the output from Flye (it may be useful to decide whether you need more sequencing reads, or some strange assembly artifact occured) and a `unicycler_output` directory containing the output from Unicycler.


# Acknowledgements


Many thanks are owed to Ryan Wick (https://github.com/rrwick), who not only wrote Unicycler and some other code used in plassembler, but also gave me  ideas about how to approach the plasmid assembly issue. If you are doing any bacterial genome assembly, you should read all of his work.

# Version Log

A brief description of what is new in each update of plassembler can be found in the HISTORY.md file.

# Bugs and Suggestions

If you come across bugs with plassembler, or would like to make any suggestions to improve the program, please open an issue or email george.bouras@adelaide.edu.au.

# Other Future Directions

At the moment, plassembler is designed for users with hybrid ONT long read and matching short read data. However, with the new Kit 14 chemistry, ONT long reads may be accurate enough that short read sequencing is not required to polish bacterial assemblies. However, I am not aware of any studies regarding the recovery of small plasmids - it is possible that Kit 14 chemistries may miss these, much like R9.4.1 chemistries, therefore necessitating short reads for plasmid recovery.

Further, other approaches may be more appropriate for Kit 14 long read only assemblies - see this [tweet](https://twitter.com/rrwick/status/1548926644085108738?cxt=HHwWhMClvfCk8v4qAAAA). 

In theory, plassembler can be applied to Pacbio reads not just ONT reads - it would only take a small tweak to change Flye parameters. If you would like this functionality, please let me know, I just haven't had the usecase for it.

# Citations

If you use plassembler, please cite:

* Kolmogorov, M., Yuan, J., Lin, Y. et al. Assembly of long, error-prone reads using repeat graphs. Nat Biotechnol 37, 540–546 (2019). https://doi.org/10.1038/s41587-019-0072-8
* Li H., Minimap2: pairwise alignment for nucleotide sequences, Bioinformatics, Volume 34, Issue 18 Pages 3094–3100 (2018), https://doi.org/10.1093/bioinformatics/bty191
* Li H., Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM. arXiv preprint arXiv:1303.3997 (2013).
* Wick RR, Judd LM, Gorrie CL, Holt KE Unicycler: Resolving bacterial genome assemblies from short and long sequencing reads. PLoS Comput Biol 13(6): e1005595 (2017). https://doi.org/10.1371/journal.pcbi.1005595
* Heng Li, Bob Handsaker, Alec Wysoker, Tim Fennell, Jue Ruan, Nils Homer, Gabor Marth, Goncalo Abecasis, Richard Durbin, 1000 Genome Project Data Processing Subgroup, The Sequence Alignment/Map format and SAMtools, Bioinformatics, Volume 25, Issue 16, 15 August 2009, Pages 2078–2079, https://doi.org/10.1093/bioinformatics/btp352
* Wick RR, Judd LM, Wyres KL, Holt KE. Recovery of small plasmid sequences via Oxford Nanopore sequencing. Microb Genom. 2021 Aug;7(8):000631. doi: 10.1099/mgen.0.000631. PMID: 34431763; PMCID: PMC8549360.
* Shen W, Le S, Li Y, Hu F (2016) SeqKit: A Cross-Platform and Ultrafast Toolkit for FASTA/Q File Manipulation. PLoS ONE 11(10): e0163962. https://doi.org/10.1371/journal.pone.0163962.
* Schmartz GP, Hartung A, Hirsch P, Kern F, Fehlmann T, Müller R, Keller A, PLSDB: advancing a comprehensive database of bacterial plasmids, Nucleic Acids Research, Volume 50, Issue D1, 7 January 2022, Pages D273–D278, https://doi.org/10.1093/nar/gkab1111.
* Ondov, B.D., Treangen, T.J., Melsted, P. et al. Mash: fast genome and metagenome distance estimation using MinHash. Genome Biol 17, 132 (2016). https://doi.org/10.1186/s13059-016-0997-x.
