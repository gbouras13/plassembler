##### Plassembler Quality Control

The following example shows how `plassembler` can be used for assembling small plasmids in a bacterial isolate with multiple chromosomes (or chromosome sized replicons).

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

The below was run on a Macbook M1 (2022).

From the paper, it was known that the smaller chromosome was 1.9Mbp, with the larger being 3.3 Mbp. So I decided on using  a `-c` value of 1500000. There were lots of long reads in this sample set, so this was have been a good usecase for `--use_raven` to speed up assembly, as the long read set was quite deep.

```
plassembler run -d Plassembler_DB  -l SRR8335319_1.fastq.gz  -1 SRR8335320_1.fastq.gz  -2 SRR8335320_2.fastq.gz  -o vibrio -t 8 -f -r --use_raven
```




```
2023-06-01 23:53:07.723 | INFO     | plassembler:begin_plassembler:73 - You are using Plassembler version 1.1.0
2023-06-01 23:53:07.723 | INFO     | plassembler:begin_plassembler:74 - Repository homepage is https://github.com/gbouras13/plassembler
2023-06-01 23:53:07.724 | INFO     | plassembler:begin_plassembler:75 - Written by George Bouras: george.bouras@adelaide.edu.au
2023-06-01 23:53:07.724 | INFO     | plassembler:run:312 - Database directory is ../Plassembler_DB
2023-06-01 23:53:07.724 | INFO     | plassembler:run:313 - Longreads file is SRR8335319_1.fastq.gz
2023-06-01 23:53:07.724 | INFO     | plassembler:run:314 - R1 fasta file is SRR8335320_1.fastq.gz
2023-06-01 23:53:07.724 | INFO     | plassembler:run:315 - R2 fasta file is SRR8335320_2.fastq.gz
2023-06-01 23:53:07.724 | INFO     | plassembler:run:316 - Chromosome length threshold is 1500000
2023-06-01 23:53:07.724 | INFO     | plassembler:run:317 - Output directory is vibrio
2023-06-01 23:53:07.724 | INFO     | plassembler:run:318 - Min long read length is 500
2023-06-01 23:53:07.724 | INFO     | plassembler:run:319 - Min long read quality is 9
2023-06-01 23:53:07.724 | INFO     | plassembler:run:320 - Thread count is 16
2023-06-01 23:53:07.724 | INFO     | plassembler:run:321 - --force is True
2023-06-01 23:53:07.724 | INFO     | plassembler:run:322 - --skip_qc is False
2023-06-01 23:53:07.725 | INFO     | plassembler:run:323 - --raw_flag is False
2023-06-01 23:53:07.725 | INFO     | plassembler:run:324 - --pacbio_model is nothing
2023-06-01 23:53:07.725 | INFO     | plassembler:run:325 - --keep_fastqs is False
2023-06-01 23:53:07.725 | INFO     | plassembler:run:326 - --keep_chromosome is False
2023-06-01 23:53:07.725 | INFO     | plassembler:run:330 - Checking dependencies
2023-06-01 23:53:07.969 | INFO     | plassembler.utils.input_commands:check_dependencies:140 - Flye version found is v2.9.2-b1786.
2023-06-01 23:53:07.970 | INFO     | plassembler.utils.input_commands:check_dependencies:150 - Flye version is ok.
2023-06-01 23:53:08.451 | INFO     | plassembler.utils.input_commands:check_dependencies:159 - Raven v1.8.1 found.
2023-06-01 23:53:08.452 | INFO     | plassembler.utils.input_commands:check_dependencies:161 - Raven version is ok.
2023-06-01 23:53:08.838 | INFO     | plassembler.utils.input_commands:check_dependencies:190 - Unicycler version found is v0.5.0.
2023-06-01 23:53:08.839 | INFO     | plassembler.utils.input_commands:check_dependencies:203 - Unicycler version is ok.
2023-06-01 23:53:09.620 | INFO     | plassembler.utils.input_commands:check_dependencies:213 - SPAdes v3.15.2 found.
2023-06-01 23:53:10.420 | INFO     | plassembler.utils.input_commands:check_dependencies:226 - Samtools v1.17 found.
2023-06-01 23:53:10.495 | INFO     | plassembler.utils.input_commands:check_dependencies:237 - minimap2 v2.26-r1175 found.
2023-06-01 23:53:10.652 | INFO     | plassembler.utils.input_commands:check_dependencies:248 - fastp v0.23.4 found.
2023-06-01 23:53:11.114 | INFO     | plassembler.utils.input_commands:check_dependencies:259 - chopper v0.5.0 found.
2023-06-01 23:53:11.372 | INFO     | plassembler.utils.input_commands:check_dependencies:274 - mash v2.3 found.
2023-06-01 23:53:11.373 | INFO     | plassembler.utils.input_commands:check_dependencies:279 - All dependencies found.
2023-06-01 23:53:11.373 | INFO     | plassembler:run:335 - Checking database installation.
2023-06-01 23:53:11.374 | INFO     | plassembler.utils.db:check_db_installation:24 - PLSDB Database at ../Plassembler_DB has already been downloaded
2023-06-01 23:53:11.374 | INFO     | plassembler:run:338 - Database successfully checked.
2023-06-01 23:53:11.374 | INFO     | plassembler:run:341 - Checking input fastqs.
2023-06-01 23:53:11.386 | INFO     | plassembler.utils.input_commands:validate_fastq:24 - FASTQ SRR8335319_1.fastq.gz checked
2023-06-01 23:53:11.387 | INFO     | plassembler.utils.input_commands:validate_fastq:24 - FASTQ SRR8335320_1.fastq.gz checked
2023-06-01 23:53:11.388 | INFO     | plassembler.utils.input_commands:validate_fastq:24 - FASTQ SRR8335320_2.fastq.gz checked
2023-06-01 23:53:11.388 | INFO     | plassembler:run:349 - FASTQ file SRR8335319_1.fastq.gz compression is True
2023-06-01 23:53:11.388 | INFO     | plassembler:run:350 - FASTQ file SRR8335320_1.fastq.gz compression is True
2023-06-01 23:53:11.388 | INFO     | plassembler:run:351 - FASTQ file SRR8335320_2.fastq.gz compression is True
2023-06-01 23:53:11.388 | INFO     | plassembler:run:359 - Filtering long reads with chopper
2023-06-01 23:53:11.388 | INFO     | plassembler.utils.qc:chopper:25 - Started running chopper
2023-06-01 23:55:41.167 | INFO     | plassembler.utils.qc:chopper:82 - Finished running chopper
2023-06-01 23:55:41.169 | INFO     | plassembler:run:372 - You have specified --use_raven. Using Raven for long read assembly.
2023-06-01 23:55:41.169 | INFO     | plassembler:run:375 - Running Raven.
2023-06-01 23:55:41.173 | INFO     | plassembler.utils.external_tools:run_to_stdout:64 - Started running raven -t 16 vibrio/chopper_long_reads.fastq.gz --graphical-fragment-assembly vibrio/assembly_graph.gfa ...
2023-06-02 00:06:08.238 | INFO     | plassembler.utils.external_tools:run_to_stdout:66 - Done running raven -t 16 vibrio/chopper_long_reads.fastq.gz --graphical-fragment-assembly vibrio/assembly_graph.gfa
2023-06-02 00:06:08.240 | INFO     | plassembler:run:388 - Counting Contigs.
2023-06-02 00:06:08.253 | INFO     | plassembler.utils.plass_class:get_contig_count:83 - Assembled 3 contigs.
2023-06-02 00:06:08.253 | INFO     | plassembler:run:566 - More than one contig was assembled with Raven.
2023-06-02 00:06:08.253 | INFO     | plassembler:run:567 - Extracting Chromosome.
2023-06-02 00:06:08.313 | INFO     | plassembler:run:605 - Chromosome Identified. Plassembler will now use long and short reads to assemble plasmids accurately.
2023-06-02 00:06:08.314 | INFO     | plassembler:run:607 - Mapping long reads.
2023-06-02 00:06:08.317 | INFO     | plassembler.utils.external_tools:run_to_stdout:64 - Started running minimap2 -ax map-ont -t 16 vibrio/flye_renamed.fasta vibrio/chopper_long_reads.fastq.gz ...
2023-06-02 00:06:47.157 | INFO     | plassembler.utils.external_tools:run_to_stdout:66 - Done running minimap2 -ax map-ont -t 16 vibrio/flye_renamed.fasta vibrio/chopper_long_reads.fastq.gz
2023-06-02 00:06:47.160 | INFO     | plassembler:run:622 - Trimming short reads.
2023-06-02 00:06:47.161 | INFO     | plassembler.utils.external_tools:run:55 - Started running fastp --out1 vibrio/trimmed_R1.fastq --out2 vibrio/trimmed_R2.fastq --in1 SRR8335320_1.fastq.gz --in2 SRR8335320_2.fastq.gz ...
2023-06-02 00:06:54.904 | INFO     | plassembler.utils.external_tools:run:57 - Done running fastp --out1 vibrio/trimmed_R1.fastq --out2 vibrio/trimmed_R2.fastq --in1 SRR8335320_1.fastq.gz --in2 SRR8335320_2.fastq.gz
2023-06-02 00:06:54.906 | INFO     | plassembler:run:626 - Mapping short reads.
2023-06-02 00:06:54.907 | INFO     | plassembler.utils.external_tools:run_to_stdout:64 - Started running minimap2 -ax sr -t 16 vibrio/flye_renamed.fasta vibrio/trimmed_R1.fastq vibrio/trimmed_R2.fastq ...
2023-06-02 00:07:04.573 | INFO     | plassembler.utils.external_tools:run_to_stdout:66 - Done running minimap2 -ax sr -t 16 vibrio/flye_renamed.fasta vibrio/trimmed_R1.fastq vibrio/trimmed_R2.fastq
2023-06-02 00:07:04.574 | INFO     | plassembler:run:633 - Processing Sam/Bam Files and extracting Fastqs.
2023-06-02 00:08:34.377 | INFO     | plassembler.utils.external_tools:run_to_stdout:64 - Started running samtools view -h -@ 16 -b vibrio/short_read.sam ...
2023-06-02 00:08:36.638 | INFO     | plassembler.utils.external_tools:run_to_stdout:66 - Done running samtools view -h -@ 16 -b vibrio/short_read.sam
2023-06-02 00:08:36.640 | INFO     | plassembler.utils.external_tools:run_to_stdout:64 - Started running samtools view -b -h -@ 16 -L vibrio/non_chromosome.bed vibrio/short_read.bam ...
2023-06-02 00:08:36.986 | INFO     | plassembler.utils.external_tools:run_to_stdout:66 - Done running samtools view -b -h -@ 16 -L vibrio/non_chromosome.bed vibrio/short_read.bam
2023-06-02 00:08:36.988 | INFO     | plassembler.utils.external_tools:run_to_stdout:64 - Started running samtools view -b -h -f 4 -@ 16 vibrio/short_read.bam ...
2023-06-02 00:08:37.298 | INFO     | plassembler.utils.external_tools:run_to_stdout:66 - Done running samtools view -b -h -f 4 -@ 16 vibrio/short_read.bam
2023-06-02 00:08:37.299 | INFO     | plassembler.utils.external_tools:run_to_stdout:64 - Started running samtools view -b -h -@ 16 -L vibrio/chromosome.bed vibrio/short_read.bam ...
2023-06-02 00:08:39.683 | INFO     | plassembler.utils.external_tools:run_to_stdout:66 - Done running samtools view -b -h -@ 16 -L vibrio/chromosome.bed vibrio/short_read.bam
2023-06-02 00:08:39.685 | INFO     | plassembler.utils.external_tools:run:55 - Started running samtools fastq -@ 16 vibrio/unmapped_bam_file.bam -1 vibrio/unmapped_R1.fastq -2 vibrio/unmapped_R2.fastq -0 /dev/null -s /dev/null -n ...
2023-06-02 00:08:39.754 | INFO     | plassembler.utils.external_tools:run:57 - Done running samtools fastq -@ 16 vibrio/unmapped_bam_file.bam -1 vibrio/unmapped_R1.fastq -2 vibrio/unmapped_R2.fastq -0 /dev/null -s /dev/null -n
2023-06-02 00:08:39.755 | INFO     | plassembler.utils.external_tools:run:55 - Started running samtools fastq -@ 16 vibrio/non_chromosome.bam -1 vibrio/mapped_non_chromosome_R1.fastq -2 vibrio/mapped_non_chromosome_R2.fastq -0 /dev/null -s /dev/null -n ...
2023-06-02 00:08:39.831 | INFO     | plassembler.utils.external_tools:run:57 - Done running samtools fastq -@ 16 vibrio/non_chromosome.bam -1 vibrio/mapped_non_chromosome_R1.fastq -2 vibrio/mapped_non_chromosome_R2.fastq -0 /dev/null -s /dev/null -n
2023-06-02 00:08:40.300 | INFO     | plassembler:run:645 - Running Unicycler.
2023-06-02 00:08:40.301 | INFO     | plassembler.utils.external_tools:run:55 - Started running unicycler -1 vibrio/short_read_concat_R1.fastq -2 vibrio/short_read_concat_R2.fastq -l vibrio/plasmid_long.fastq -t 16 -o vibrio/unicycler_output ...
2023-06-02 00:13:25.791 | INFO     | plassembler:run:662 - Unicycler identified plasmids. Calculating Plasmid Copy Numbers.
2023-06-02 00:13:25.846 | INFO     | plassembler.utils.external_tools:run_to_stdout:64 - Started running minimap2 -ax map-ont -t 16 vibrio/combined.fasta vibrio/chopper_long_reads.fastq.gz ...
2023-06-02 00:14:10.006 | INFO     | plassembler.utils.external_tools:run_to_stdout:66 - Done running minimap2 -ax map-ont -t 16 vibrio/combined.fasta vibrio/chopper_long_reads.fastq.gz
2023-06-02 00:14:10.011 | INFO     | plassembler.utils.external_tools:run:55 - Started running samtools sort -@ 16 vibrio/combined_long.sam -o vibrio/combined_sorted_long.bam ...
2023-06-02 00:14:14.531 | INFO     | plassembler.utils.external_tools:run:57 - Done running samtools sort -@ 16 vibrio/combined_long.sam -o vibrio/combined_sorted_long.bam
2023-06-02 00:14:14.533 | INFO     | plassembler.utils.external_tools:run_to_stdout:64 - Started running minimap2 -ax sr -t 16 vibrio/combined.fasta vibrio/trimmed_R1.fastq vibrio/trimmed_R2.fastq ...
2023-06-02 00:14:26.051 | INFO     | plassembler.utils.external_tools:run_to_stdout:66 - Done running minimap2 -ax sr -t 16 vibrio/combined.fasta vibrio/trimmed_R1.fastq vibrio/trimmed_R2.fastq
2023-06-02 00:14:26.053 | INFO     | plassembler.utils.external_tools:run:55 - Started running samtools sort -@ 16 vibrio/combined_short.sam -o vibrio/combined_sorted_short.bam ...
2023-06-02 00:14:28.060 | INFO     | plassembler.utils.external_tools:run:57 - Done running samtools sort -@ 16 vibrio/combined_short.sam -o vibrio/combined_sorted_short.bam
2023-06-02 00:15:02.818 | INFO     | plassembler:run:670 - Calculating mash distances to PLSDB.
2023-06-02 00:15:02.822 | INFO     | plassembler.utils.external_tools:run:55 - Started running mash sketch vibrio/plasmids.fasta -i ...
2023-06-02 00:15:03.149 | INFO     | plassembler.utils.external_tools:run:57 - Done running mash sketch vibrio/plasmids.fasta -i
2023-06-02 00:15:03.151 | INFO     | plassembler.utils.external_tools:run_to_stdout:64 - Started running mash dist vibrio/plasmids.fasta.msh ../Plassembler_DB/plsdb.msh -v 0.1 -d 0.1 -i ...
2023-06-02 00:15:05.353 | INFO     | plassembler.utils.external_tools:run_to_stdout:66 - Done running mash dist vibrio/plasmids.fasta.msh ../Plassembler_DB/plsdb.msh -v 0.1 -d 0.1 -i
2023-06-02 00:15:08.006 | INFO     | plassembler:end_plassembler:90 - Plassembler has finished
2023-06-02 00:15:08.006 | INFO     | plassembler:end_plassembler:91 - Elapsed time: 1335.09 seconds
```

Intrestingly, I found a small 5386bp contig that isn't in the PRJNA479421 assembly! It has a short read copy number of around 4, but was missed in the long read set. 

It turns out this is phage [phiX174](https://en.wikipedia.org/wiki/Phi_X_174), which is commonly used as a spike in positive control in Illumina sequencing. Therefore, it likely reflects contamination in this sample.

However, this is also a good example showing that Plassembler can also be used to recover phages and phage-plasmids in hybrid sequencing data, so long as they have not integrated into the bacterial chromosome(s) (like prophages).