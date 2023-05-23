##### Plassembler Quality Control

The following example shows how Plassembler can be used for assembling small plasmids in a bacterial isolate with multiple chromosomes (or chromosome sized replicons).

We will assemble _Vibrio ampbellii DS40M4_ from Bioproject PRJNA479421, which is a _Vibrio_ that has 2 chromosomes. You can read more about this bacterium in the [paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6346181/). The authors hybrid sequenced this with ONT long and illumina short reads. 

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

The below was run on my Macbook M1 (2022).

From the paper, I knew that the smaller chromosome was 1.9Mbp, with the larger being 3.3 Mbp. So I decided on using  a `-c` value of 1500000. There were lots of long reads in this sample set, so this was have been a good usecase for `--use_raven` to speed up assembly, as the long read set was quite deep.


```
plassembler.py -d plassembler_db  -l SRR8335319_1.fastq.gz  -1 SRR8335320_1.fastq.gz  -2 SRR8335320_2.fastq.gz  -o vibrio -t 8 -f -r --use_raven

```

The terminal output looks like this.

```
Starting plassembler v1.1.0
Checking dependencies.
Flye version found is v2.9.2-b1786.
Flye version is ok.
Raven v1.8.1 found.
Raven version is ok.
Unicycler version found is v0.5.0.
Unicycler version is ok.
SPAdes v3.15.2 found.
Samtools v1.9 found.
minimap2 v2.26-r1175 found.
fastp v0.22.0 found.
chopper v0.5.0 found.
mash v2.2.2 found.
All dependencies found.
Checking database installation.
Database successfully checked.
Checking input fastqs.
FASTQ SRR8335319_1.fastq.gz checked
FASTQ SRR8335320_1.fastq.gz checked
FASTQ SRR8335320_2.fastq.gz checked
Filtering long reads with chopper.
Kept 60090 reads out of 88896 reads
You have specified --use_raven. Using Raven for long read assembly.
Running Raven.
Counting Contigs.
Raven assembled 3 contigs.
More than one contig was assembled with Raven.
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
Elapsed time: 1277.44 seconds.
```

Intrestingly, I found a small 5386bp contig that isn't in the PRJNA479421 assembly! It has a short read copy number of around 4, but was missed in the long read set. 

It turns out this is phage [phiX174](https://en.wikipedia.org/wiki/Phi_X_174), which is commonly used as a spike in positive control in Illumina sequencing. Therefore, it likely reflects contamination in this sample.

However, this is also a good example showing that Plassembler can also be used to recover phages and phage-plasmids in hybrid sequencing data, so long as they have not integrated into the bacterial chromosome(s) (like prophages).