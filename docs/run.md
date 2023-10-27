
To run `plassembler`, first you need to install the database in a directory of your chosing:

`plassembler download -d <database directory>`

Once this is finished, you can run Plassembler as follows:

`plassembler run -d <database directory> -l <long read fastq> -o <output dir> -1 < short read R1 fastq> -2 < short read R2 fastq>  -c <estimated lower bound of chromosome length>`

* `-c` or `--chromosome` will default to 1000000 if not specified.

To specify more threads to speed up Plassembler, use `-t` or `--threads`:

`plassembler run -d <database directory> -l <long read fastq> -o <output dir> -1 < short read R1 fastq> -2 < short read R2 fastq>  -c <estimated chromosome length> -t <threads>`

`plassembler` defaults to 1 thread.

To specify a prefix for the output files, use `-p` or `--prefix`:

`plassembler run -d <database directory> -l <long read fastq> -o <output dir> -1 < short read R1 fastq> -2 < short read R2 fastq>  -c <estimated chromosome length> -t <threads> -p <prefix>`

To specify a minimum length and minimum read quality Q-score for chopper, use `-m` and `-q` :

`plassembler run -d <database directory> -l <long read fastq> -o <output dir> -1 < short read R1 fastq> -2 < short read R2 fastq>  -c <estimated chromosome length> -t <threads> -p <prefix> -m <min length> -q <min quality>`

* `-m` will default to 500 and `-q` will default to 9. 

To overwrite an existing output directory, use `-f` or `--force`

`plassembler run -d <database directory> -l <long read fastq> -o <output dir> -1 < short read R1 fastq> -2 < short read R2 fastq>  -c <estimated chromosome length> -t <threads> `

To use Raven instead of Flye as a long read assembler, use `--use_raven`

`plassembler run -d <database directory> -l <long read fastq> -o <output dir> -1 < short read R1 fastq> -2 < short read R2 fastq>  -c <estimated chromosome length> -t <threads> --use_raven`

To keep the Flye assembled chromosome(s) (as `chromosome.fasta`), use `--keep-chromosome`

`plassembler run -d <database directory> -l <long read fastq> -o <output dir> -1 < short read R1 fastq> -2 < short read R2 fastq>  -c <estimated chromosome length> -t <threads>   --keep_chromosome`

To use pacbio reads use `--pacbio_model` (e.g. with regular CLR reads so with `pacbio-raw` model specified in Flye):

`plassembler run -d <database directory> -l <long read fastq> -o <output dir> -1 < short read R1 fastq> -2 < short read R2 fastq>  -c <estimated chromosome length> -t <threads>  --pacbio_model pacbio-raw`

To skip quality control (chopper and fastp), use `--skip_qc` 

`plassembler run -d <database directory> -l <long read fastq> -o <output dir> -1 < short read R1 fastq> -2 < short read R2 fastq>  -c <estimated chromosome length> -t <threads> --skip_qc`

To specify a directory containing an existing Flye assembly for your long reads  use `--flye_directory` 

`plassembler run -d <database directory> -l <long read fastq> -o <output dir> -1 < short read R1 fastq> -2 < short read R2 fastq>  -c <estimated chromosome length> -t <threads> --flye_directory <flye directory>`

To use assembled mode to calculate plasmid copy numbers, you need to use `plassembler assembled`, along with an already assembled chromosome with `--input_chromosome` and plasmids with `--input_plasmids`.

`plassembler assembled -d <database directory> -l <long read fastq> -o <output dir> -1 < short read R1 fastq> -2 < short read R2 fastq>  -c <estimated chromosome length> -t <threads>  -a --input_chromosome <path to chromosome FASTA> --input_plasmids <path to plasmids FASTA> `


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
  --flye_directory PATH     Directory containing Flye long read assembly.
                            Needs to contain assembly_info.txt and
                            assembly_info.fasta. Allows Plassembler to Skip
                            Flye assembly step.
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

All options 

```
Usage: plassembler [OPTIONS] COMMAND [ARGS]...

Options:
  -h, --help     Show this message and exit.
  -V, --version  Show the version and exit.

Commands:
  assembled  Runs assembled mode
  citation   Print the citation(s) for this tool
  download   Downloads Plassembler DB
  long       Plassembler with long reads only - experimental and untested
  run        Runs Plassembler
  ```