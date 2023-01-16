
To run Plassembler, first you need to install the database in a directory of your chosing:

`install_database.py -d <database directory>`

Once this is finished, you can run Plassembler as followed

` plassembler.py -d <database directory> -l <long read fastq> -o <output dir> -1 < short read R1 fastq> -2 < short read R2 fastq>  -c <estimated chromosome length>`

* -c will default to 2500000 if not specified (a lower bound for the estimated genome length of Staphylococcus aureus).

To specify more threads to speed up Plassembler:

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
usage: plassembler.py [-h] -d DATABASE -l LONGREADS [-1 SHORT_ONE] [-2 SHORT_TWO] [-c CHROMOSOME] [-o OUTDIR] [-m MIN_LENGTH]
                      [-t THREADS] [-f] [-r] [-p PREFIX] [-q MIN_QUALITY] [-k] [-V]

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
  -k, --kmer_mode       Very high quality Nanopore R10.4 and above reads. No short reads required. Experimental for now.
  -V, --version         show plassembler version and exit.
```
