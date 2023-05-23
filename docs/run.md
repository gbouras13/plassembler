
To run Plassembler, first you need to install the database in a directory of your chosing:

`install_database.py -d <database directory>`

Once this is finished, you can run Plassembler as follows:

` plassembler.py -d <database directory> -l <long read fastq> -o <output dir> -1 < short read R1 fastq> -2 < short read R2 fastq>  -c <estimated chromosome length>`

* -c will default to 1000000 if not specified.

To specify more threads to speed up Plassembler:

` plassembler.py -d <database directory> -l <long read fastq> -o <output dir> -1 < short read R1 fastq> -2 < short read R2 fastq>  -c <estimated chromosome length> -t <threads>`

plassembler defaults to 1 thread.

To specify a prefix for the output files:

` plassembler.py -d <database directory> -l <long read fastq> -o <output dir> -1 < short read R1 fastq> -2 < short read R2 fastq>  -c <estimated chromosome length> -t <threads> -p <prefix>`

To specify a minimum length and minimum read quality Q-score for chopper :

` plassembler.py -d <database directory> -l <long read fastq> -o <output dir> -1 < short read R1 fastq> -2 < short read R2 fastq>  -c <estimated chromosome length> -t <threads> -p <prefix> -m <min length> -q <min quality>`

* -m will default to 500 and -q will default to 9. Note that for some tiny plasmids, -m may need to be reduced (see this [paper](https://www.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000631#tab2) ).

To overwrite an existing output directory, use `-f`

` plassembler.py -d <database directory> -l <long read fastq> -o <output dir> -1 < short read R1 fastq> -2 < short read R2 fastq>  -c <estimated chromosome length> -t <threads> `

To use Raven instead of Flye as a long read assembler, use `--use_raven`.

` plassembler.py -d <database directory> -l <long read fastq> -o <output dir> -1 < short read R1 fastq> -2 < short read R2 fastq>  -c <estimated chromosome length> -t <threads> --use_raven `

To keep the Flye assembled chromosome(s) (as `chromosome.fasta`)

` plassembler.py -d <database directory> -l <long read fastq> -o <output dir> -1 < short read R1 fastq> -2 < short read R2 fastq>  -c <estimated chromosome length> -t <threads>   --keep_chromosome`

To use pacbio reads (e.g. with regular CLR reads so with `pacbio-raw` model specified in Flye):

` plassembler.py -d <database directory> -l <long read fastq> -o <output dir> -1 < short read R1 fastq> -2 < short read R2 fastq>  -c <estimated chromosome length> -t <threads>  --pacbio_model pacbio-raw`

To use assembled mode to calculate plasmid copy numbers, you need `-a`, along with an already assembled chromosome with `--input_chromosome` and plasmids with `--input_plasmids`.

` plassembler.py -d <database directory> -l <long read fastq> -o <output dir> -1 < short read R1 fastq> -2 < short read R2 fastq>  -c <estimated chromosome length> -t <threads>  -a --input_chromosome <path to chromosome FASTA> --input_plasmids <path to plasmids FASTA> `

You can also use `--long_only` mode, which will simply run Flye and keep all contigs below `-c` and denote them as 'plasmids', but this is experimental only for now and I do not vouch for its performance.

```
usage: plassembler.py [-h] -d DATABASE [-l LONGREADS] [-1 SHORT_ONE]
                      [-2 SHORT_TWO] [-c CHROMOSOME] [-o OUTDIR]
                      [-m MIN_LENGTH] [-q MIN_QUALITY] [-t THREADS] [-f]
                      [-p PREFIX] [--pacbio_model PACBIO_MODEL] [-r]
                      [--keep_fastqs] [--keep_chromosome] [-a]
                      [--input_chromosome INPUT_CHROMOSOME]
                      [--input_plasmids INPUT_PLASMIDS] [--long_only]
                      [--use_raven] [-V]

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
  -p PREFIX, --prefix PREFIX
                        Prefix for output files. This is not required. Defaults to plassembler.
  --pacbio_model PACBIO_MODEL
                        Pacbio Flye model. 
                        Must be one of pacbio-raw, pacbio-corr or pacbio-hifi. 
                        Use pacbio-raw for PacBio regular CLR reads (<20 percent error), pacbio-corr for PacBio reads that were corrected with other methods (<3 percent error) or pacbio-hifi for PacBio HiFi reads (<1 percent error).
  -r, --raw_flag        Use --nano-raw for Flye. 
                        Designed for Guppy fast configuration reads. 
                        By default, Flye will assume SUP or HAC reads and use --nano-hq
  --keep_fastqs         Whether you want to keep FASTQ files containing putative plasmid reads 
                        and long reads that map to multiple contigs (plasmid and chromosome).
  --keep_chromosome     If you want to keep the chromosome assembly.
  -a, --assembled_mode  Activates assembled mode.
  --input_chromosome INPUT_CHROMOSOME
                        Input FASTA file consisting of already assembled chromosome with assembled mode. 
                        Must be 1 complete contig.
  --input_plasmids INPUT_PLASMIDS
                        Input FASTA file consisting of already assembled plasmids with assembled mode. 
                        Requires FASTQ file input (short only, long only or long + short).
  --long_only           Experimental for now. 
                        Very high quality Nanopore R10.4 and above reads. 
                        Assembly using Flye, extracts contigs under size -c and runs the depth arguments. 
                        No short reads required.
  --use_raven           Uses Raven instead of Flye for long read assembly. 
                        May be useful if you want to reduce runtime.
  -V, --version         show plassembler version and exit.
```
