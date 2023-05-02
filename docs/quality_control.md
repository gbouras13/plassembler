##### Plassembler Quality Control

The following example shows how Plassembler can be used for quality control (checking your long and short read sets come from the same isolate).

We will assemble SAMN32360844, a _S. aureus_ ST22 isolate long read set with a SAMN32360859 short read set, a distinct _S. aureus_ ST22 isolate. These samples are extremely closely related, being part of the same sequence type. We will also assemble the SAMN32360844 long read set with the SAMN32360837 short read set, a more distantly related ST30 isolate. These samples are taken from Houtak et al (2023) _The Intra-Host Evolutionary Landscape And Pathoadaptation Of Persistent Staphylococcus aureus In Chronic Rhinosinusitis_ available on bioRxiv [here](https://doi.org/10.1101/2023.03.28.534496) if you want more details.

Firstly, I downloaded the fastqs from the SRA using the fantastic [fastq-ql](https://github.com/rpetit3/fastq-dl) program (after installation with mamba).

```
# installation
mamba create -n fastq-dl fastq-dl
conda activate fastq-dl

# SAMN32360844 long reads
fastq-dl SRR22859710	

# SAMN32360859 short reads
fastq-dl SRR22859826		


# SAMN32360837 short reads
fastq-dl SRR22859924
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

Plassembler was run using the default --nano-hq mode, as these reads were seuqenced SUP basecalled in 2022 with Guppy v6. They were assembled with 16 threads using -t 16 and a chromosome size -c of 2500000.

# Same Sequence Type

As you can see below, plassembler recovered 10 non-circular contigs with no hits to PLSDB, indictating that the long and short read sets likely don't match.


```
# Wick et al
plassembler.py -d plassembler_db  -l SRR22859710_1.fastq.gz  -1  SRR22859826_1.fastq.gz -2 SRR22859826_2.fastq.gz   -o S_Aureus_QC_Illustration_same_sequence_type -t 16 -f -c 2500000
```

```
Starting plassembler v1.0.0
Checking dependencies.
Flye version found is v2.9.1-b1780.
Flye version is ok.
Unicycler version found is v0.5.0.
Unicycler version is ok.
Samtools v1.16.1 found.
minimap2 v2.24-r1122 found.
fastp v0.23.2 found.
chopper v0.5.0 found.
seqkit v2.3.0 found.
mash v2.3 found.
rasusa v0.7.1 found.
All dependencies found.
Checking database installation.
Database successfully checked.
Checking input fastqs.
FASTQ SRR22859710_1.fastq.gz checked
FASTQ SRR22859826_1.fastq.gz checked
FASTQ SRR22859826_2.fastq.gz checked
Filtering long reads with chopper.
Kept 60568 reads out of 80687 reads
Subsampling long reads with rasusa.
Running Flye.
Counting Contigs.
Flye assembled 2 contigs.
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
WARNING: 10 non-circular contigs with no PLSDB mash hits were detected. 
This indicates your long and short read sets may come from different bacterial isolates. 
Please check this!
plassembler has finished
Elapsed time: 803.24 seconds

```

# Different Sequence Type

For the different sequence type isolates, plassembler took double as long as with the same sequence type isolates, which is not surprising as there would have been many more unmapped reads to assemble with Unicycler.

As you can see below, plassembler recovered 59 non-circular contigs with no hits to PLSDB, indictating that the long and short read sets don't match.

```
plassembler.py -d plassembler_db -l SRR22859710_1.fastq.gz  -1  SRR22859826_1.fastq.gz -2 SRR22859826_2.fastq.gz   -o S_Aureus_QC_Illustration_same_sequence_type -t 16 -f -c 2500000
```

```
Starting plassembler v1.0.0
Checking dependencies.
Flye version found is v2.9.1-b1780.
Flye version is ok.
Unicycler version found is v0.5.0.
Unicycler version is ok.
Samtools v1.16.1 found.
minimap2 v2.24-r1122 found.
fastp v0.23.2 found.
chopper v0.5.0 found.
seqkit v2.3.0 found.
mash v2.3 found.
rasusa v0.7.1 found.
All dependencies found.
Checking database installation.
Database successfully checked.
Checking input fastqs.
FASTQ SRR22859710_1.fastq.gz checked
FASTQ SRR22859924_1.fastq.gz checked
FASTQ SRR22859924_2.fastq.gz checked
Filtering long reads with chopper.
Kept 60568 reads out of 80687 reads
Subsampling long reads with rasusa.
Running Flye.
Counting Contigs.
Flye assembled 2 contigs.
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
WARNING: 59 non-circular contigs with no PLSDB mash hits were detected. 
This indicates your long and short read sets may come from different bacterial isolates. 
Please check this!
plassembler has finished
Elapsed time: 1678.74 seconds
```





