##### Plassembler Quality Control

The following example shows how Plassembler can be used for assembling small plasmids in a bacterial isolate with multiple chromosomes (or chromosome sized replicons).

We will assemble _Vibrio ampbellii DS40M4_ from Bioproject PRJNA479421, which is a _Vibrio_ that has 2 chromosomes. You can read more about this bacterium in the [paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6346181/). The authors hybrid sequenced this bacterium with ONT long and illumina short reads. 


#### vibrio


Firstly, I downloaded the fastqs from the SRA using the fantastic [fastq-ql](https://github.com/rpetit3/fastq-dl) program (after installation with mamba).

```
# installation
mamba create -n fastq-dl fastq-dl
conda activate fastq-dl

# downloads all the read sets
fastq-dl PRJNA479421	

conda deactivate
```

#### Running Plassembler 

First install plassembler. If you need Mac see the [instructions](https://github.com/gbouras13/plassembler#installation).
And then download the database.

```
mamba create -n plassembler plassembler 
conda activate plassembler
install_database.py plassembler_db
```

The below was run on my Mac Mini M1 (2021).

Plassembler was run using the default the `-r` flag, as these Nanopore reads were sequenced in 2019. They were assembled with 16 threads.

From the paper, I knew that the smaller chromosome was 1.9Mbp, with the larger being 3.3 Mbp. So I decided on using a `-s` value of 100 along with a `-c` value of 1500000 - this means approximately a 30x coverage of the chromosome (Plassembler will keep 1500000x100=150Mbp of long reads, with a combined chromosome size of approximately 5.2 Mbp). 


```
plassembler.py -d plassembler_db  -l SRR8335319_1.fastq.gz  -1 SRR8335320_1.fastq.gz  -2 SRR8335320_2.fastq.gz  -o vibrio -t 16 -s 100 -c 1500000 -f -r

```

The terminal output looks like this.

```
Starting plassembler v1.0.0
Checking dependencies.
Flye version found is v2.9.2-b1786.
Flye version is ok.
Unicycler version found is v0.5.0.
Unicycler version is ok.
SPAdes v3.15.2 found.
Samtools v1.17 found.
minimap2 v2.24-r1122 found.
fastp v0.23.2 found.
chopper v0.5.0 found.
seqkit v2.4.0 found.
mash v2.3 found.
rasusa v0.7.1 found.
All dependencies found.
Checking database installation.
Database successfully checked.
Checking input fastqs.
FASTQ SRR8335319_1.fastq.gz checked
FASTQ SRR8335320_1.fastq.gz checked
FASTQ SRR8335320_2.fastq.gz checked
Filtering long reads with chopper.
Kept 60090 reads out of 88896 reads
Subsampling long reads with rasusa.
Running Flye.
Counting Contigs.
Flye assembled 4 contigs.
More than one contig was assembled with Flye.
Extracting Chromosome.
Chromosome Identified. Plassembler will now use long and short reads to assemble plasmids accurately.
Mapping Long Reads.
Trimming short reads.
Mapping short reads.
Processing Sam/Bam Files and extracting Fastqs.
Running Unicycler.
Calculating Plasmid Copy Numbers.
Calculating mash distances to PLSDB.
Plassembler has finished.
Elapsed time: 1107.47 seconds.
```

Intrestingly, I found a small 5386bp plasmid that isn't in the PRJNA479421 assembly! 