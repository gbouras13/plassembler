[![Paper](https://img.shields.io/badge/paper-Bioinformatics-teal.svg?style=flat-square&maxAge=3600)](https://doi.org/10.1093/bioinformatics/btad409)
[![CI](https://github.com/gbouras13/plassembler/actions/workflows/ci.yaml/badge.svg)](https://github.com/gbouras13/plassembler/actions/workflows/ci.yaml)
[![BioConda Install](https://img.shields.io/conda/dn/bioconda/plassembler.svg?style=flag&label=BioConda%20install)](https://anaconda.org/bioconda/plassembler)
[![codecov](https://codecov.io/gh/gbouras13/plassembler/branch/main/graph/badge.svg?token=4B1T2PGM9V)](https://codecov.io/gh/gbouras13/plassembler)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

[![Anaconda-Server Badge](https://anaconda.org/bioconda/plassembler/badges/version.svg)](https://anaconda.org/bioconda/plassembler)
[![Bioconda Downloads](https://img.shields.io/conda/dn/bioconda/plassembler)](https://img.shields.io/conda/dn/bioconda/plassembler)
[![PyPI version](https://badge.fury.io/py/plassembler.svg)](https://badge.fury.io/py/plassembler)
[![Downloads](https://static.pepy.tech/badge/plassembler)](https://pepy.tech/project/plassembler)
[![DOI](https://zenodo.org/badge/514596389.svg)](https://zenodo.org/doi/10.5281/zenodo.10035954)


# plassembler

## Automated Bacterial Plasmid Assembly Program

`plassembler` is a program that is designed for automated & fast assembly of plasmids in  bacterial genomes that have been hybrid sequenced with long read & paired-end short read sequencing. It was originally designed for Oxford Nanopore Technologies long reads, but will also work with Pacbio reads. As of v1.3.0, it should also work well for long-read only assembled genomes (although we would still recommend getting short reads too if you can).

If you are assembling a small number of bacterial genomes manually, I would recommend starting by using [Trycycler](https://github.com/rrwick/Trycycler) to recover the chromosome before using Plassembler to recover plasmids, especially the small ones. If you have more genomes or want to assemble your genomes in a more automated way, try [dragonflye](https://github.com/rpetit3/dragonflye), especially if you are used to Shovill, or even better my own pipeline [hybracter](https://github.com/gbouras13/hybracter) that is more appropriate for large datasets and implemented Plassembler in it.  

Additionally, I would recommend reading the following guides to bacterial genome assembly regardless of whether you want to use Plassembler:
*  [Trycycler](https://github.com/rrwick/Trycycler/wiki/Guide-to-bacterial-genome-assembly)
*  [Perfect Bacterial Assembly Tutorial](https://github.com/rrwick/Perfect-bacterial-genome-tutorial)
*  [Perfect bacterial assembly Paper](https://doi.org/10.1371/journal.pcbi.1010905)

## Quick Start

The easiest way to install `plassembler` is via conda:

`conda install -c bioconda plassembler`

Followed by database download and installation:

`plassembler download -d <databse directory>`

And finally run `plassembler`:

`plassembler run -d <database directory> -l <long read fastq> -o <output dir> -1 < short read R1 fastq> -2 < short read R2 fastq>  -c <estimated chromosome length>`

Please read the [Installation](#installation) section for more details, especially if you are an inexperienced command line user.

## Manuscript

`plassembler` has been recently published in *Bioinformatics*:

George Bouras, Anna E. Sheppard, Vijini Mallawaarachchi, Sarah Vreugde, Plassembler: an automated bacterial plasmid assembly tool, Bioinformatics, Volume 39, Issue 7, July 2023, btad409, https://doi.org/10.1093/bioinformatics/btad409.

If you use `plassembler`, please see the full [Citations](#citations) section for a list of all programs `plassembler` uses under the hood, in order to fully recognise the creators of these tools for their work.

## Documentation

The full documentation for Plassembler can be found [here](https://plassembler.readthedocs.io/en/latest).

## Table of Contents 

- [plassembler](#plassembler)
  - [Automated Bacterial Plasmid Assembly Program](#automated-bacterial-plasmid-assembly-program)
  - [Quick Start](#quick-start)
  - [Manuscript](#manuscript)
  - [Documentation](#documentation)
  - [Table of Contents](#table-of-contents)
  - [`plassembler` v1.5.0 Update New Database (21 November 2023)](#plassembler-v150-update-new-database-21-november-2023)
  - [`plassembler` v1.3.0 Updates (24 October 2023)](#plassembler-v130-updates-24-october-2023)
  - [Why Does Plassembler Exist?](#why-does-plassembler-exist)
  - [Why Not Just Use Unicycler?](#why-not-just-use-unicycler)
  - [Other Features](#other-features)
  - [Quality Control](#quality-control)
  - [Metagenomes](#metagenomes)
  - [Installation](#installation)
    - [Conda](#conda)
    - [Pip](#pip)
    - [Source](#source)
  - [Unicycler v0.5.0 Installation Issues](#unicycler-v050-installation-issues)
  - [Running plassembler](#running-plassembler)
  - [Outputs](#outputs)
  - [Benchmarking](#benchmarking)
  - [Acknowledgements](#acknowledgements)
  - [Version Log](#version-log)
  - [Bugs and Suggestions](#bugs-and-suggestions)
  - [Citations](#citations)

## `plassembler` v1.5.0 Update New Database (21 November 2023)

* **If you upgrade to v1.5.0, you will need to update the database using `plassembler download`** 
* Plassembler v1.5.0 incorporates a new expanded database thanks to the recent PLSDB release [2023_11_03_v2](https://ccb-microbe.cs.uni-saarland.de/plsdb/). Thanks @[biobrad](https://github.com/biobrad) for the heads up.


## `plassembler` v1.3.0 Updates (24 October 2023)

* `plassembler long` should yield improved results. It achieves this by treating long reads as both short reads (in the sense of creating a de Brujin graph based short read assembly to begin) and long reads (for scaffolding) in Unicycler.
* While I'd still recommend short reads if you can get them, I am now confident that if your isolate has small plasmids in the long read set, `plassembler long` is very likely to find and recover them.
* For more information, see the [documentation](https://plassembler.readthedocs.io/en/latest/long/).
* The ability to specify a `--flye_assembly` and `--flye_info` if you already have a Flye assembly for your long reads instead of `--flye_directory` has been added. Thanks to @[incoherentian](https://github.com/incoherentian)'s [issue](https://github.com/gbouras13/plassembler/issues/37)
* The ability to specify a `--no_copy_numbers` with `plassembler assembled` if you just want to run some plasmids against the PLSDB has been added. Thanks to @[gaworj](https://github.com/gaworj)'s [issue](https://github.com/gbouras13/plassembler/issues/36).

## Why Does Plassembler Exist?

In long-read assembled bacterial genomes, small plasmids are difficult to assemble correctly with long read assemblers. They commonly have circularisation issues and can be duplicated or missed (see [this](https://doi.org/10.1371/journal.pcbi.1010905), [this](https://f1000research.com/articles/8-2138) and [this](https://github.com/rrwick/Trycycler/wiki/Clustering-contigs)). This recent [paper](https://www.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.001024) in _Microbial Genomics_ by Johnson et al also suggests that long read assemblers particularly miss small plasmids.

`plassembler` was therefore created as a fast automated tool to ensure plasmids are assembled correctly without duplicated regions for high-throughput uses - like Unicycler but a lot laster - and to provide some useful statistics as well (such as estimate plasmid copy numbers for both long and short read sets).  

As it turns out (though this wasn't a motivation for making it), `plassembler` also recovers more small plasmids than the existing gold standard tool Unicycler. I think this is because it throws away chromosomal reads, similar to subsampling short reads sets which can improve recovery. As there are more plasmid reads a proportion of the overall read set, there seems to be a higher chance of recovering smaller plasmids.

You can see this increase in accuracy and speed in the benchmarking results for [simulated](docs/benchmarking_results_sim.md) and [real](docs/benchmarking_results_real.md) datasets.

Plassembler also uses [mash](https://github.com/marbl/Mash) as a quick way to determine whether each assembled contig has any similar hits in [PLSDB](https://doi.org/10.1093/nar/gkab1111). 

Additionally, due to its mapping approach, Plassembler can also be used as a quality control tool for checking whether your long and short read sets come from the same isolate. This may be particularly useful if your read sets come from different extractions, or you have multiplexed many samples (& want to avoid mislabelling).  

## Why Not Just Use Unicycler?

Unicycler is awesome and still a good way to assemble plasmids from hybrid sequencing - `plassembler` uses it! But there are a few reasons to use plassembler instead:

1. Time. Plassember throws away all the chromosomal reads (i.e. most of them) before running Unicycler, so it is much faster (wall clock 3-10x faster generally). 
2. Accuracy. Benchmarking has shown `plassembler` is better than Unicycler in terms of recovering small plasmids.
3. `plassembler` will output only the likely plasmids, and can more easily be integrated into pipelines. You shouldn't be assembling the chromosome using Unicycler [anymore](https://doi.org/10.1371/journal.pcbi.1010905) so `plassembler` can get you only what is necessary from Unicycler.
4. `plassembler` will give you summary depth and copy number stats for both long and short reads.
5. `plassembler` can be used as a quality control to check if your short and long reads come from the same sample - if `plassembler` results in many non-circular contigs (particularly those that have no hits in PLSDB), it is likely because your read sets do not come from the same isolate! See [Quality Control](#quality-control).
6. You will get information whether each assembled contig has a similar entry in [PLSDB](https://doi.org/10.1093/nar/gkab1111). Especially for common pathogen species that are well represented in databases, this will likely tell you specifically what plasmid you have in your sample. 
* Note: Especially for less commonly sequenced species, I would not suggest that that absence of a PLSDB hit is necessary meaningful, especially for circular contigs - those would likely be novel plasmids uncaptured by PLSDB.

## Other Features 

1. Assembled mode.

* Thanks to a suggestion from [gaworj](https://github.com/gaworj), assembled mode has been added to Plassembler. This allows you to calculate the copy numbers of already assembled plasmids you may have, skipping assembly. 

You can use this feature with `plassembler assembled`.

2. Multi-mapped reads.

* All long reads that map to multiple contigs (mostly, reads that map to both the chromosome and plasmids, but also to multiple putative plasmids) will be extracted when using the `--keep-fastqs` options. These may be of interest if you are looking at shared mobile genetic elements.

3. Multiple chromosome bacteria/megaplasmids/chromids

* Plassembler should work with bacteria with multiple chromosomes, megaplasmids or chromids. In this case, I would treat the megaplasmids etc like chromosomes and assemble them using a long-read first approach with Trycycler or Dragonflye, as they are of approximately chromosome size. 
* I'd still use Plassembler to recover small plasmids - for example, for  Plassembler v1.1.0 recovered the 77.5 kbp plasmiod along with a 5386bp contig (coresponding to phage phiX174, a common sequencing spike-in) in the _Vibrio campbellii DS40M4_  (see this [paper](https://doi.org/10.1128/MRA.01187-18) and this [bioproject](https://www.ncbi.nlm.nih.gov/bioproject/479421) ).
* `-c` needs to be smaller than the size of the largest chromosome-like element.
* For example, for the vibrio example, which had approximately 1.8Mbp and 3.3Mbp chromosomes , I used `-c 1500000`.

Please see [here](docs/multiple_chromosomes.md) for more details and an example. 

4. Phages, Phage-Plasmids and Other Extrachromosomal Replicons

* If you have sufficient hybrid sequencing data, Plassembler will theoretically recover assemblies of all non-chromosomal replicons, including phages and phage-plasmids
* A good example of this is the _Vibrio campbellii DS40M4_  example, where Plassembler recovered the assembly of phage phiX174, albeit it was from sequencing spike-in contamination in that case.

## Quality Control

* `plassembler` can also be used for quality control to test whether your long and short read sets come from the same isolate, even within the same species.

Please see [here](docs/quality_control.md) for more details and some examples.

## Metagenomes

* `plassembler` is not currently recommended for metagenomic datasets, because of their high diversity, leading to difficulties in recovering chromosome-length contigs for bacteria. Additionally, Unicycler is not recommended for metagenomes. However,  `plassembler` was tested on a high depth very simple mock community dataset from this [paper](https://www.nature.com/articles/s41592-022-01539-7). It worked quite nicely, recovering the 5 known plasmids, but we don't anticipate it will work as well on your data! If you try it and it works please let us know.

Please see [here](docs/metagenomics.md) for more details.

## Installation

Plassembler has been tested on Linux and MacOS machines. 

### Conda

The easiest way to install `plassembler` is via conda - Plassembler is on bioconda. 

```
conda install -c bioconda plassembler
```

or mamba for quicker solving:

```
mamba install -c bioconda plassembler
```

This will install all the dependencies along with `plassembler`.

### Pip

You can install the Python components of `plassembler` using pip.

```
pip install plassembler
```

You will then need to install the external dependencies separately, which can be found in `build/environment.yaml`

* [Flye](https://github.com/fenderglass/Flye) >=2.9
* [Unicycler](https://github.com/rrwick/Unicycler) >=0.4.8
* [Minimap2](https://github.com/lh3/minimap2) >=2.11
* [fastp](https://github.com/OpenGene/fastp) >=0.18.0
* [chopper](https://github.com/wdecoster/chopper) >=0.5.0
* [mash](https://github.com/marbl/Mash) >=2.2
* [Raven](https://github.com/lbcb-sci/raven) >=1.8
* [Samtools](https://github.com/samtools/samtools) >=0.15.0

### Source

Alternatively, the development version of `plassembler` can be installed manually via github.

```
git clone https://github.com/gbouras13/plassembler.git
cd plassembler
pip install -e .
```

## Unicycler v0.5.0 Installation Issues

`plassembler` works best with Unicycler v0.5.0. With Unicycler v0.4.8, `plassembler` should still run without any issue and provide a satisfactory assembly, but you will be warned of this when you run `plassembler`. `plassembler` will not work with any older version of Unicycler.

**Linux**

For Linux environments, Unicycler v0.5.0 should be installed automaticall with the `plassembler` bioconda installation.

You can force it as follows:

`conda install -c bioconda plassembler unicycler==0.5.0`

or manually install Unicycler v0.5.0 after installing `plassembler`:

```
conda install -c bioconda plassembler
pip3 install git+https://github.com/rrwick/Unicycler.git
```

**MacOS**

For MacOS environments, the current conda installation method will only install the latest available bioconda Unicycler version of v0.4.8. 

Ryan Wick (the author of Unicycler) suggests that v0.5.0 should be used, as v0.4.8 is not compatible with the latest versions of spades (see [here](https://github.com/rrwick/Unicycler/releases/tag/v0.5.0) ). This will require another installation step on MacOS.

To install Unicycler v0.5.0, it is recommended that you install Unicycler from github after installing Plassembler follows:

```
# installs plassembler into an environment called 'plassemblerENV' and activates it
conda create -n plassemblerENV plassembler
conda activate plassemblerENV
# installs Unicycler v0.5.0
pip3 install git+https://github.com/rrwick/Unicycler.git
```

Mac M1 users may need to change some compiler settings and install from the Unicycler github repo e.g.

```
# installs plassembler into an environment called 'plassemblerENV' and activates it
conda create -n plassemblerENV plassembler
conda activate plassemblerENV
# installs Unicycler v0.5.0
git clone https://github.com/rrwick/Unicycler.git
cd Unicycler
python3 setup.py install --makeargs "CXX=g++"
```

## Running plassembler

To run `plassembler`, first you need to install the database in a directory of your chosing:

`plassembler download -d <database directory>`

Once this is finished, you can run plassembler as follows:

`plassembler run -d <database directory> -l <long read fastq> -o <output dir> -1 < short read R1 fastq> -2 < short read R2 fastq>  -c <estimated chromosome length>`

* -c will default to 1000000 if it is absent.

To specify more threads:

`plassembler run -d <database directory> -l <long read fastq> -o <output dir> -1 < short read R1 fastq> -2 < short read R2 fastq>  -c <estimated chromosome length> -t <threads>`

Plassembler defaults to 1 thread.

To specify a prefix for the output files:

`plassembler run -d <database directory> -l <long read fastq> -o <output dir> -1 < short read R1 fastq> -2 < short read R2 fastq>  -c <estimated chromosome length> -t <threads> -p <prefix>`

To specify a minimum long read length and minimum read quality Q-score for filtering with [chopper](https://github.com/wdecoster/chopper):

`plassembler run -d <database directory> -l <long read fastq> -o <output dir> -1 < short read R1 fastq> -2 < short read R2 fastq>  -c <estimated chromosome length> -t <threads> -p <prefix> -m <min length> -q <min quality>`

* -m will default to 500 and -q will default to 9. Note that for some tiny plasmids, -m should be reduced or perhaps even set to 1 (see this [paper](https://www.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000631#tab2) ).

To overwrite an existing output directory, use -f

` plassembler run -d <database directory> -l <long read fastq> -o <output dir> -1 < short read R1 fastq> -2 < short read R2 fastq>  -c <estimated chromosome length> -t <threads> -p <prefix> -m <min length> -q <min quality> -f`

To use Raven instead of Flye as a long read assembler, use `--use_raven`.

`plassembler run -d <database directory> -l <long read fastq> -o <output dir> -1 < short read R1 fastq> -2 < short read R2 fastq>  -c <estimated chromosome length> -t <threads> --use_raven `

Please see the [documentation](docs/run.md) for more options.

```
Usage: plassembler run [OPTIONS]

  Runs Plassembler

Options:
  -h, --help                Show this message and exit.
  -V, --version             Show the version and exit.
  -d, --database PATH       Directory of PLSDB database.  [required]
  -l, --longreads PATH      FASTQ file of long reads.  [required]
  -1, --short_one PATH      R1 short read FASTQ file.  [required]
  -2, --short_two PATH      R2 short read FASTQ file.  [required]
  -c, --chromosome INTEGER  Approximate lower-bound chromosome length of
                            bacteria (in base pairs).  [default: 1000000]
  -o, --outdir PATH         Directory to write the output to.  [default:
                            plassembler.output/]
  -m, --min_length TEXT     minimum length for filtering long reads with
                            chopper.  [default: 500]
  -q, --min_quality TEXT    minimum quality q-score for filtering long reads
                            with chopper.  [default: 9]
  -t, --threads TEXT        Number of threads.  [default: 1]
  -f, --force               Force overwrites the output directory.
  -p, --prefix TEXT         Prefix for output files. This is not required.
                            [default: plassembler]
  --skip_qc                 Skips qc (chopper and fastp).
  --pacbio_model TEXT       Pacbio model for Flye.  Must be one of pacbio-raw,
                            pacbio-corr or pacbio-hifi.  Use pacbio-raw for
                            PacBio regular CLR reads (<20 percent error),
                            pacbio-corr for PacBio reads that were corrected
                            with other methods (<3 percent error) or pacbio-
                            hifi for PacBio HiFi reads (<1 percent error).
  -r, --raw_flag            Use --nano-raw for Flye.  Designed for Guppy fast
                            configuration reads.  By default, Flye will assume
                            SUP or HAC reads and use --nano-hq.
  --keep_fastqs             Whether you want to keep FASTQ files containing
                            putative plasmid reads  and long reads that map to
                            multiple contigs (plasmid and chromosome).
  --keep_chromosome         If you want to keep the chromosome assembly.
  --use_raven               Uses Raven instead of Flye for long read assembly.
                            May be useful if you want to reduce runtime.
  --flye_directory PATH     Directory containing Flye long read assembly.
                            Needs to contain assembly_info.txt and
                            assembly_info.fasta. Allows Plassembler to Skip
                            Flye assembly step.
  --flye_assembly PATH      Path to file containing Flye long read assembly
                            FASTA. Allows Plassembler to Skip Flye assembly
                            step in conjunction with  --flye_info.
  --flye_info PATH          Path to file containing Flye long read assembly
                            info text file. Allows Plassembler to Skip Flye
                            assembly step in conjunction with
                            --flye_assembly.
  --no_chromosome           Run Plassembler assuming no chromosome can be
                            assembled. Use this if your reads only contain
                            plasmids that you would like to assemble.
```

## Outputs

Plassembler will output a `_plasmids.fasta` file, which will contain the assembled plasmid sequence(s) in FASTA format (including long and short read copy numbers in the header), and a `_plasmids.gfa` file, which will contain the assembly graph from Unicycler that can be visualised in [Bandage](https://github.com/rrwick/Bandage). 

Plassembler also outputs a `_summary.tsv` file, which gives the estimated copy number for each plasmid, for both short reads and long reads (see this [paper](https://www.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000631#tab2) for more details about plasmid copy numbers) and also gives each contig's top hit by mash distance in the PLSDB (if there is a hit), along with all its supporting information. 

If `plassembler` fails to assemble any plasmids at all in `_plasmids.fasta`, all these files will still exist, but will be empty (to ensure `plassembler` can be easily integrated into workflow managers like Snakemake).

`plassembler` will also output a log file, a `flye_output` directory, which contains the output from Flye (it may be useful to decide whether you need more sequencing reads, or some strange assembly artifact occured) and a `unicycler_output` directory containing the output from Unicycler. If `--use_raven` is specified, a `raven_output` directory will be present instead.

## Benchmarking

The benchmarking results for [simulated](docs/benchmarking_results_sim.md) and [real](docs/benchmarking_results_real.md) datasets are available. The full benchmarking output can be found [here](https://doi.org/10.5281/zenodo.7996690).

All benchmarking was conducted on a Intel® Core™ i7-10700K CPU @ 3.80GHz on a machine running Ubuntu 20.04.6 LTS. 

Tldr: Plassembler is much faster than Unicycler (3-10x usually) and is more accurate because it is more likely to recover low coverage plasmids that Unicycler might miss.

## Acknowledgements

Many thanks are owed to [Ryan Wick](https://github.com/rrwick), who not only wrote Unicycler and some other code used in Plassembler, but also gave me some initial ideas about how to approach the plasmid assembly problem originally. If you are doing any bacterial genome assembly, you should read all of his work, but if you have read this far you probably already have.

Also thanks to [Vijini Mallawaarachchi](https://github.com/Vini2) who helped refactor the code - if you are interested in recovering phages (especially in the metagenome context) please give [phables](https://github.com/Vini2/phables) a go.

## Version Log

A brief description of what is new in each update of `plassembler` can be found in the HISTORY.md file.

## Bugs and Suggestions

If you come across bugs with `plassembler`, or would like to make any suggestions to improve the program, please open an issue or email george.bouras@adelaide.edu.au.

## Citations

`plassembler` has been recently published in *Bioinformatics*:

George Bouras, Anna E. Sheppard, Vijini Mallawaarachchi, Sarah Vreugde, Plassembler: an automated bacterial plasmid assembly tool, Bioinformatics, Volume 39, Issue 7, July 2023, btad409, https://doi.org/10.1093/bioinformatics/btad409.

If you use `plassembler`, please also consider citing where relevant:

* Kolmogorov, M., Yuan, J., Lin, Y. et al. Assembly of long, error-prone reads using repeat graphs. Nat Biotechnol 37, 540–546 (2019). https://doi.org/10.1038/s41587-019-0072-8
* Li H., Minimap2: pairwise alignment for nucleotide sequences, Bioinformatics, Volume 34, Issue 18 Pages 3094–3100 (2018), https://doi.org/10.1093/bioinformatics/bty191
* Wick RR, Judd LM, Gorrie CL, Holt KE Unicycler: Resolving bacterial genome assemblies from short and long sequencing reads. PLoS Comput Biol 13(6): e1005595 (2017). https://doi.org/10.1371/journal.pcbi.1005595
* Heng Li, Bob Handsaker, Alec Wysoker, Tim Fennell, Jue Ruan, Nils Homer, Gabor Marth, Goncalo Abecasis, Richard Durbin, 1000 Genome Project Data Processing Subgroup, The Sequence Alignment/Map format and SAMtools, Bioinformatics, Volume 25, Issue 16, 15 August 2009, Pages 2078–2079, https://doi.org/10.1093/bioinformatics/btp352
* Wick RR, Judd LM, Wyres KL, Holt KE. Recovery of small plasmid sequences via Oxford Nanopore sequencing. Microb Genom. 2021 Aug;7(8):000631. doi: 10.1099/mgen.0.000631. PMID: 34431763; PMCID: PMC8549360.
* Schmartz GP, Hartung A, Hirsch P, Kern F, Fehlmann T, Müller R, Keller A, PLSDB: advancing a comprehensive database of bacterial plasmids, Nucleic Acids Research, Volume 50, Issue D1, 7 January 2022, Pages D273–D278, https://doi.org/10.1093/nar/gkab1111.
* Ondov, B.D., Treangen, T.J., Melsted, P. et al. Mash: fast genome and metagenome distance estimation using MinHash. Genome Biol 17, 132 (2016). https://doi.org/10.1186/s13059-016-0997-x.
* De Coster,W. and Rademakers,R. (2023) NanoPack2: population-scale evaluation of long-read sequencing data. Bioinformatics, 39, btad311. https://doi.org/10.1093/bioinformatics/btad311.
* Vaser,R. and Šikić,M. (2021) Time-and memory-efficient genome assembly with Raven. Nat. Comput. Sci., 1, 332–336. https://doi.org/10.1038/s43588-021-00073-4.
* Koren S, Walenz BP, Berlin K, Miller JR, Bergman NH, Phillippy AM. (2017) Canu: scalable and accurate long-read assembly via adaptive k-mer weighting and repeat separation. Genome Res. 2017 May;27(5):722-736. doi: https://doi.org/10.1101/gr.215087.116. 
* Bouras, G., Roach, M. J., Mallawaarachchi V., Grigson., S., Papudeshi., B. (2023) Dnaapler: A tool to reorient circular microbial genomes https://github.com/gbouras13/dnaapler
