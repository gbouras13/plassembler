plassembler
===============

Automated Bacterial Plasmid Assembly Program
------------

Plassembler is a program that is designed for automated & fast assembly of plasmids in  bacterial genomes that have been hybrid sequenced with long read & paired-end short read sequencing. It was originally designed for Oxford Nanopore Technologies long reads, but will also work with Pacbio reads. 

If you are assembling a small number of bacterial genomes manually, I would highly recommend starting by using [Trycycler](https://github.com/rrwick/Trycycler) to recover the chromosome before using Plassembler. If you have more genomes or want to assemble your genomes in a more automated way, try [dragonflye](https://github.com/rpetit3/dragonflye), especially if you are used to Shovill, or my own pipeline [hybracter](https://github.com/gbouras13/hybracter) that is more appropriate for large datasets.  

Additionally, I would recommend reading the following guides to bacterial genome assembly regardless of whether you want to use Plassembler:
*  [Trycycler](https://github.com/rrwick/Trycycler/wiki/Guide-to-bacterial-genome-assembly)
*  [Perfect Bacterial Assembly Tutorial](https://github.com/rrwick/Perfect-bacterial-genome-tutorial)
*  [Perfect bacterial assembly Paper](https://doi.org/10.1371/journal.pcbi.1010905)

Table of Contents
-----------
- [plassembler](#plassembler)
  - [Automated Bacterial Plasmid Assembly Program](#automated-bacterial-plasmid-assembly-program)
  - [Table of Contents](#table-of-contents)
- [Quick Start](#quick-start)
- [Why Does Plassembler Exist?](#why-does-plassembler-exist)
- [Why Not Just Use Unicycler?](#why-not-just-use-unicycler)
- [Documentation](#documentation)
- [Method](#method)
- [Other Features](#other-features)
- [Quality Control](#quality-control)
- [Installation](#installation)
  - [Unicycler v0.5.0 Installation Issues](#unicycler-v050-installation-issues)
- [Running plassembler](#running-plassembler)
- [Outputs](#outputs)
- [Benchmarking](#benchmarking)
- [Acknowledgements](#acknowledgements)
- [Version Log](#version-log)
- [Bugs and Suggestions](#bugs-and-suggestions)
- [Other Future Directions](#other-future-directions)
- [Citations](#citations)

# Quick Start

The easiest way to install plassembler is via conda:

`conda install -c bioconda plassembler`

Followed by database download and installation:

`install_database.py -d <path/to/databse_dir>`

And finally assembly:

`plassembler.py -d <database directory> -l <long read fastq> -o <output dir> -1 < short read R1 fastq> -2 < short read R2 fastq>  -c <estimated chromosome length>`

Please read the [Installation](#installation) section for more details, especially if you are an inexperienced command line user.

# Why Does Plassembler Exist?

In long-read first assembled bacterial genomes, small plasmids are often difficult to assemble correctly with long read assemblers such as Flye. They often have circularisation issues and can be duplicated or missed (see [this](https://doi.org/10.1371/journal.pcbi.1010905) and [this](https://github.com/rrwick/Trycycler/wiki/Clustering-contigs) ).

Plassembler was therefore created as an automated tool to ensure plasmids assemble correctly without duplicated regions for high-throughput uses - and to provide some useful statistics as well (such as estimate plasmid copy numbers for both long and short read sets). Plassembler will likely also recover small plasmids that long read assemblers like Flye simply miss.

As it turns out (though this wasn't a motivation for making it!), Plassembler also assembles plasmids more accurately than Unicycler. I think this is because subsampling short reads sets can improve assembly quality (see [this](https://doi.org/10.1093/bioinformatics/btv311) and [this](https://doi.org/10.1371/journal.pone.0060204)), so throwing away chromosomal reads probably has a similar effect.

Plassembler also uses [mash](https://github.com/marbl/Mash) as a quick way to determine whether each assembled contig has any similar hits in [PLSDB](https://doi.org/10.1093/nar/gkab1111). 

Additionally, due to its mapping approach, Plassembler can also be used as a quality control tool for checking whether your long and short read sets come from the same isolate. This may be particularly useful if your read sets come from different extractions, or you have multiplexed many samples (& want to avoid mislabelling).  

# Why Not Just Use Unicycler?

Unicycler is awesome and still a good way to assemble plasmids from hybrid sequencing - plassembler uses it! But there are a few reasons to use plassembler instead:

1. Time. Plassember throws away all the chromosomal reads (i.e. most of them) before running Unicycler, so it is much faster (4-20x, and will be higher if you have lots of long reads). 
2. Accuracy. Benchmarking has shown Plassembler is more accurate than Unicycler in terms of indels and mismatches. I honestly wasn't even aiming for this when I wrote Plassembler, just a big speed-up, but it is a nice result of course!
3. Plassembler will output only the likely plasmids, and can more easily be integrated into pipelines. You shouldn't be assembling the chromosome using Unicycler [anymore](https://doi.org/10.1371/journal.pcbi.1010905) so plassembler can get you only what is necessary from Unicycler.
4. Plassembler will give you summary depth and copy number stats for both long and short reads.
5. Plassembler can be used as a quality control to check if your short and long reads come from the same sample - if plassembler results in many non-circular contigs (particularly those that have no hits in PLSDB), it is likely because your read sets do not come from the same isolate! 
6. You will get information whether each assembled contig has a similar entry in [PLSDB](https://doi.org/10.1093/nar/gkab1111). Especially for common pathogen species that are well represented in databases, this will likely tell you specifically what plasmid you have in your sample. 
* Note: Especially for less commonly sequenced species, I would not suggest that that absence of a PLSDB hit is necessary meaningful, especially for circular contigs - those would likely be novel plasmids uncaptured by PLSDB.


# Documentation

Documentation can be found at http://plassembler.readthedocs.io/.

# Method

1. Long reads are filtered using [chopper](https://github.com/wdecoster/chopper) .
2. Long reads are subsampled to 30x to reduce runtime (subsampling depth be changed with `-s` or skipped with `--no_subsample`).
3. Long-read only assembly is conducted with [Flye](https://github.com/fenderglass/Flye).
4. If the resulting assembly is checked. Contigs bigger than the provided chromosome size `-c`, are identified as chromosomal and extracted. Any other contigs are extracted as putative plasmid contigs, if Flye assembled any. If no contigs were larger than `-c`, plassembler will exit - you probably need to get some more long reads to complete your assembly (or check `-c` wasn't too big).
5. Short reads are filtered using [fastp](https://github.com/OpenGene/fastp).
6. Long and short reads are mapped to a reference containing the chromosomal contigs plus putative plasmid contigs using [minimap2](https://github.com/lh3/minimap2#uguide).
7. All reads that map to the putative plasmid contigs and all reads that are unmapped are extracted and combined.
8. These reads are assembled using the hybrid assembler [Unicycler](https://github.com/rrwick/Unicycler) to generate final plasmid contigs.
9. Average read coverage depth for each plasmid is calculated using a modified version of code found [here](https://github.com/rrwick/Small-plasmid-Nanopore). See also this [paper](https://www.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000631#tab2).
10. Plasmid copy number is calculated by dividing the plasmid read depth by the chromosome read depth.
11. All plasmid contigs are compared against [PLSDB](https://doi.org/10.1093/nar/gkab1111) using [mash](https://github.com/marbl/Mash) with a cutoff maximum mash distance of 0.1.


# Other Features 

1. Assembled mode.

* Thanks to a suggestion from @gaworj, assembled mode has been added to Plassembler from v1.0.0. This allows you to calculate the copy numbers of already assembled plasmids you may have, skipping assembly. You can specify this with the `-a` flag, along with your chromosome FASTA file using `--input-chromosome` and your plasmids FASTA file `--input_plasmids`.

2. Multi-mapped reads.

* All long reads that map to multiple contigs (mostly, reads that map to both the chromosome and plasmids, but also to multiple putative plasmids) will be extracted when using the `--keep-fastqs` options. These may be of interest if you are looking at shared mobile genetic elements.

3. Multiple chromosome bacteria/megaplasmids/chromids

* Plassembler should work with bacteria with multiple chromosomes, megaplasmids or chromids. In this case, I would treat the megaplasmids etc like chromosomes and assemble them using a long-read first approach with Trycycler or Dragonflye, as they are of approximately chromosome size. 
* However, I'd still use Plassembler to recover small plasmids - for example, for  it managed to recover a 5386bp plasmid in the _Vibrio campbellii DS40M4_ genome (see this [paper](https://doi.org/10.1128/MRA.01187-18) and this [bioproject](https://www.ncbi.nlm.nih.gov/bioproject/479421) ) that was missed with Unicycler v.0.4.4 (or at least not uploaded!).
* I would recommend tweaking the parameters a bit in this use case. -c needs to be smaller than the size of the largest chromosome-like element, and I would increase the subsampling depth `-s` from 30 to something higher (e.g. `-s 100`), because this is based off the `-c` value. 
* For example, for the vibrio example, which had approximately 1.8Mbp and 3.3Mbp chromosomes , I used `-c 1500000 -s 100`.

Please see [here](docs/multiple_chromosomes.md) for more details and an example. 

# Quality Control

* Plassembler can also be used for quality control to test whether your long and short read sets come from the same isolate, even within the same species.

Please see [here](docs/quality_control.md) for more details and some examples.

# Installation

Plassembler should run and has been tested on Linux and MacOS machines. 

The easiest way to install plassembler is via conda - Plassembler is on bioconda. 

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
plassembler.py --help
```

Unicycler v0.5.0 Installation Issues
------

Plassembler works best with Unicycler v0.5.0. With Unicycler v0.4.8, Plassembler should still run without any issue and provide a satisfactory assembly, but you will be warned of this when you run plassembler. Plassembler will not work with any older version.

**Linux**

For Linux environments, Unicycler v0.5.0 should be installed automaticall with the plassembler bioconda installation.

You can force it as follows:

`conda install -c bioconda plassembler unicycler==0.5.0`

or manually install Unicycler v0.5.0 after installing plassembler:

```
conda install -c bioconda plassembler
pip3 install git+https://github.com/rrwick/Unicycler.git
```

**MacOS**

For MacOS environments, the current conda installation method will only install the latest available bioconda Unicycler version of v0.4.8. 

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

Once this is finished, you can run plassembler as follows:

` plassembler.py -d <database directory> -l <long read fastq> -o <output dir> -1 < short read R1 fastq> -2 < short read R2 fastq>  -c <estimated chromosome length>`

* -c will default to 1000000 if it is absent.

To specify more threads:

` plassembler.py -d <database directory> -l <long read fastq> -o <output dir> -1 < short read R1 fastq> -2 < short read R2 fastq>  -c <estimated chromosome length> -t <threads>`

Plassembler defaults to 1 thread.

To specify a prefix for the output files:

` plassembler.py -d <database directory> -l <long read fastq> -o <output dir> -1 < short read R1 fastq> -2 < short read R2 fastq>  -c <estimated chromosome length> -t <threads> -p <prefix>`

To specify a minimum long read length and minimum read quality Q-score for filtering with [chopper](https://github.com/wdecoster/chopper):

` plassembler.py -d <database directory> -l <long read fastq> -o <output dir> -1 < short read R1 fastq> -2 < short read R2 fastq>  -c <estimated chromosome length> -t <threads> -p <prefix> -m <min length> -q <min quality>`

* -m will default to 500 and -q will default to 9. Note that for some tiny plasmids, -m should be reduced or perhaps even set to 1 (see this [paper](https://www.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000631#tab2) ).

To overwrite an existing output directory, use -f

` plassembler.py -d <database directory> -l <long read fastq> -o <output dir> -1 < short read R1 fastq> -2 < short read R2 fastq>  -c <estimated chromosome length> -t <threads> -p <prefix> -m <min length> -q <min quality> -f`

To overwrite an existing output directory, use -f

` plassembler.py -d <database directory> -l <long read fastq> -o <output dir> -1 < short read R1 fastq> -2 < short read R2 fastq>  -c <estimated chromosome length> -t <threads> -p <prefix> -m <min length> -q <min quality> -f`


```
usage: plassembler.py [-h] -d DATABASE [-l LONGREADS] [-1 SHORT_ONE] [-2 SHORT_TWO] [-c CHROMOSOME] [-o OUTDIR] [-m MIN_LENGTH]
                      [-q MIN_QUALITY] [-t THREADS] [-f] [-r] [-p PREFIX] [-s SUBSAMPLE_DEPTH] [--long_only]
                      [--pacbio_model PACBIO_MODEL] [--no_subsample] [--keep_fastqs] [--keep_chromosome] [-a]
                      [--input_chromosome INPUT_CHROMOSOME] [--input_plasmids INPUT_PLASMIDS] [-V]

plassembler: automated bacterial plasmid assembly tool.

optional arguments:
  -h, --help            show this help message and exit
  -d DATABASE, --database DATABASE
                        Directory of PLSDB database downloaded using install_database.py.
  -l LONGREADS, --longreads LONGREADS
                        FASTQ file of long reads.
  -1 SHORT_ONE, --short_one SHORT_ONE
                        R1 short read FASTQ file.
  -2 SHORT_TWO, --short_two SHORT_TWO
                        R2 short read FASTQ file.
  -c CHROMOSOME, --chromosome CHROMOSOME
                        Approximate lower-bound chromosome length of bacteria. 
                        Defaults to 1000000.
  -o OUTDIR, --outdir OUTDIR
                        Directory to write the output to. Defaults to output/
  -m MIN_LENGTH, --min_length MIN_LENGTH
                        minimum length for filtering long reads with chopper. Defaults to 500.
  -q MIN_QUALITY, --min_quality MIN_QUALITY
                        minimum quality for filtering long reads with chopper. Defaults to 9.
  -t THREADS, --threads THREADS
                        Number of threads. Defaults to 1.
  -f, --force           Overwrites the output directory.
  -r, --raw_flag        Use --nano-raw for Flye. 
                         Designed for Guppy fast configuration reads. 
                        By default, Flye will assume SUP or HAC reads and use --nano-hq
  -p PREFIX, --prefix PREFIX
                        Prefix for output files. This is not required. Defaults to plassembler.
  -s SUBSAMPLE_DEPTH, --subsample_depth SUBSAMPLE_DEPTH
                        Subsample long-read depth as an integer. 
                        Used combined with the coverage of the chromosome length provided with -c. 
                        Defaults to 30.
  --long_only           Very high quality Nanopore R10.4 and above reads. 
                        Takes Flye output, extracts contigs under size -c and runs the depth arguments. 
                        No short reads required. 
                        Experimental for now.
  --pacbio_model PACBIO_MODEL
                        Pacbio Flye model. Must be pacbio-raw, pacbio-corr or pacbio-hifi. 
                        Use pacbio-raw for PacBio regular CLR reads (<20 percent error), 
                        pacbio-corr for PacBio reads that were corrected with other methods (<3 percent error) 
                        or pacbio-hifi for PacBio HiFi reads (<1 percent error).
  --no_subsample        Turns off long-read sub-sampling. 
                        Recommended if long-read sets have low N50s/N90s, 
                        or are of a difficult-to-assemble species with lots of repeats.
  --keep_fastqs         Whether you want to keep FASTQ files containing putative plasmid reads 
                        and long reads that map to multiple contigs (plasmid and chromosome).
  --keep_chromosome     Whether you want to keep the polished Flye chromosome assembly.
  -a, --assembled_mode  Activates assembled mode.
  --input_chromosome INPUT_CHROMOSOME
                        Input FASTA file consisting of already assembled chromosome with assembled mode. 
                        Must be 1 complete contig.
  --input_plasmids INPUT_PLASMIDS
                        Input FASTA file consisting of already assembled plasmids with assembled mode. 
                        Requires FASTQ file input (short only, long only or long + short).
  -V, --version         show plassembler version and exit.
```

# Outputs

Plassembler will output a `_plasmids.fasta` file, which will contain the assembled plasmid sequence(s) in FASTA format (including long and short read copy numbers in the header), and a `_plasmids.gfa` file, which will contain the assembly graph from Unicycler that can be visualised in [Bandage](https://github.com/rrwick/Bandage). 

Plassembler also outputs a `_summary.tsv` file, which gives the estimated copy number for each plasmid, for both short reads and long reads (see this [paper](https://www.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000631#tab2) for more details about plasmid copy numbers) and also gives each contig's top hit by mash distance in the PLSDB (if there is a hit), along with all its supporting information. 

If plassembler fails to assemble any plasmids at all in `_plasmids.fasta`, all these files will still exist, but will be empty (to ensure plassembler can be easily integrated into workflow managers like Snakemake).

plassembler will also output a log file, a `flye_output` directory, which contains the output from Flye (it may be useful to decide whether you need more sequencing reads, or some strange assembly artifact occured) and a `unicycler_output` directory containing the output from Unicycler.


# Benchmarking

Please see [benchmarking](docs/benchmarking_results_simulated.md) (and [here](docs/benchmarking_results_real.md))for a full benchmarking analysis. 

Tldr: Plassembler is much faster than Unicycler (4-20x usually), is more accurate (fewer indels and mismatches) and is more likely to recover low coverage plasmids that Unicycler might miss.

# Acknowledgements

Many thanks are owed to Ryan Wick (https://github.com/rrwick), who not only wrote Unicycler and some other code used in Plassembler, but also gave me original ideas about how to approach the plasmid assembly problem originally. If you are doing any bacterial genome assembly, you should read all of his work, but if you have read this far you probably already have.

# Version Log

A brief description of what is new in each update of Plassembler can be found in the HISTORY.md file.

# Bugs and Suggestions

If you come across bugs with Plassembler, or would like to make any suggestions to improve the program, please open an issue or email george.bouras@adelaide.edu.au.

# Other Future Directions

At the moment, plassembler is designed for users with hybrid long read and matching short read data. However, with the new Kit 14 chemistry, ONT long reads may be accurate enough that short read sequencing is not required to polish bacterial assemblies - it may already be there for Pacbio! However, I am not aware of any studies regarding the recovery of small plasmids - it is possible that Kit 14 chemistries may miss these, much like R9.4.1 chemistries, therefore necessitating short reads for plasmid recovery. If you want to try it out, use `--long_only` - Plassembler will just take the Flye output contigs below `-c` and treat them as contigs. This is experimental only and untested.

Further, other approaches may be more appropriate for Kit 14 long read only assemblies - see this [tweet](https://twitter.com/rrwick/status/1548926644085108738?cxt=HHwWhMClvfCk8v4qAAAA). 

# Citations

Plassembler manuscript is under review :).

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
* Hall, M. B., (2022). Rasusa: Randomly subsample sequencing reads to a specified coverage. Journal of Open Source Software, 7(69), 3941, https://doi.org/10.21105/joss.03941.
* Wouter De Coster, Svenn D’Hert, Darrin T Schultz, Marc Cruts, Christine Van Broeckhoven, NanoPack: visualizing and processing long-read sequencing data, Bioinformatics, Volume 34, Issue 15, August 2018, Pages 2666–2669, https://doi.org/10.1093/bioinformatics/bty149.