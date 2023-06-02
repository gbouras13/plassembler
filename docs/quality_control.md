##### Plassembler Quality Control

The following example shows how Plassembler can be used for quality control (checking your long and short read sets come from the same isolate).

We will assemble SAMN32360844, a _S. aureus_ ST22 isolate long read set with a SAMN32360859 short read set, a distinct _S. aureus_ ST22 isolate. These samples are extremely closely related, being part of the same sequence type (see this [paper](https://doi.org/10.1128/JCM.38.3.1008-1015.2000) for more information on sequence typing). We will also assemble the SAMN32360844 long read set with the SAMN32360837 short read set, a more distantly related ST30 isolate. These samples are taken from Houtak et al (2023) _The Intra-Host Evolutionary Landscape And Pathoadaptation Of Persistent Staphylococcus aureus In Chronic Rhinosinusitis_ available on bioRxiv [here](https://doi.org/10.1101/2023.03.28.534496) if you want more details.

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

Plassembler was run using the default --nano-hq mode, as these reads were seuqenced SUP basecalled in 2022 with Guppy v6. They were assembled with 16 threads using -t 8 and a chromosome size -c of 2500000.

# Same Sequence Type

As you can see below, plassembler recovered 8 non-circular contigs with no hits to PLSDB, indictating that the long and short read sets likely don't match.


```
plassembler run -d plassembler_db  -l SRR22859710_1.fastq.gz  -1  SRR22859826_1.fastq.gz -2 SRR22859826_2.fastq.gz   -o S_Aureus_QC_Illustration_same_sequence_type -t 8 -f -c 2500000
```

```
(plassembler_env) a1667917@LY0TWV6HTW2 plassembler_qc_vibrio % plassembler run -d ../Plassembler_DB  -l SRR22859710_1.fastq.gz  -1  SRR22859826_1.fastq.gz -2 SRR22859826_2.fastq.gz   -o S_Aureus_QC_Illustration_same_sequence_type -t 8 -f -c 2500000
2023-06-02 00:16:38.738 | INFO     | plassembler:begin_plassembler:73 - You are using Plassembler version 1.1.0
2023-06-02 00:16:38.738 | INFO     | plassembler:begin_plassembler:74 - Repository homepage is https://github.com/gbouras13/plassembler
2023-06-02 00:16:38.738 | INFO     | plassembler:begin_plassembler:75 - Written by George Bouras: george.bouras@adelaide.edu.au
2023-06-02 00:16:38.739 | INFO     | plassembler:run:312 - Database directory is ../Plassembler_DB
2023-06-02 00:16:38.739 | INFO     | plassembler:run:313 - Longreads file is SRR22859710_1.fastq.gz
2023-06-02 00:16:38.739 | INFO     | plassembler:run:314 - R1 fasta file is SRR22859826_1.fastq.gz
2023-06-02 00:16:38.739 | INFO     | plassembler:run:315 - R2 fasta file is SRR22859826_2.fastq.gz
2023-06-02 00:16:38.739 | INFO     | plassembler:run:316 - Chromosome length threshold is 2500000
2023-06-02 00:16:38.739 | INFO     | plassembler:run:317 - Output directory is S_Aureus_QC_Illustration_same_sequence_type
2023-06-02 00:16:38.739 | INFO     | plassembler:run:318 - Min long read length is 500
2023-06-02 00:16:38.739 | INFO     | plassembler:run:319 - Min long read quality is 9
2023-06-02 00:16:38.739 | INFO     | plassembler:run:320 - Thread count is 8
2023-06-02 00:16:38.739 | INFO     | plassembler:run:321 - --force is True
2023-06-02 00:16:38.739 | INFO     | plassembler:run:322 - --skip_qc is False
2023-06-02 00:16:38.739 | INFO     | plassembler:run:323 - --raw_flag is False
2023-06-02 00:16:38.740 | INFO     | plassembler:run:324 - --pacbio_model is nothing
2023-06-02 00:16:38.740 | INFO     | plassembler:run:325 - --keep_fastqs is False
2023-06-02 00:16:38.740 | INFO     | plassembler:run:326 - --keep_chromosome is False
2023-06-02 00:16:38.740 | INFO     | plassembler:run:330 - Checking dependencies
2023-06-02 00:16:38.936 | INFO     | plassembler.utils.input_commands:check_dependencies:140 - Flye version found is v2.9.2-b1786.
2023-06-02 00:16:38.936 | INFO     | plassembler.utils.input_commands:check_dependencies:150 - Flye version is ok.
2023-06-02 00:16:39.010 | INFO     | plassembler.utils.input_commands:check_dependencies:159 - Raven v1.8.1 found.
2023-06-02 00:16:39.010 | INFO     | plassembler.utils.input_commands:check_dependencies:161 - Raven version is ok.
2023-06-02 00:16:39.291 | INFO     | plassembler.utils.input_commands:check_dependencies:190 - Unicycler version found is v0.5.0.
2023-06-02 00:16:39.292 | INFO     | plassembler.utils.input_commands:check_dependencies:203 - Unicycler version is ok.
2023-06-02 00:16:39.619 | INFO     | plassembler.utils.input_commands:check_dependencies:213 - SPAdes v3.15.2 found.
2023-06-02 00:16:39.960 | INFO     | plassembler.utils.input_commands:check_dependencies:226 - Samtools v1.17 found.
2023-06-02 00:16:39.986 | INFO     | plassembler.utils.input_commands:check_dependencies:237 - minimap2 v2.26-r1175 found.
2023-06-02 00:16:40.058 | INFO     | plassembler.utils.input_commands:check_dependencies:248 - fastp v0.23.4 found.
2023-06-02 00:16:40.177 | INFO     | plassembler.utils.input_commands:check_dependencies:259 - chopper v0.5.0 found.
2023-06-02 00:16:40.214 | INFO     | plassembler.utils.input_commands:check_dependencies:274 - mash v2.3 found.
2023-06-02 00:16:40.214 | INFO     | plassembler.utils.input_commands:check_dependencies:279 - All dependencies found.
2023-06-02 00:16:40.215 | INFO     | plassembler:run:335 - Checking database installation.
2023-06-02 00:16:40.215 | INFO     | plassembler.utils.db:check_db_installation:24 - PLSDB Database at ../Plassembler_DB has already been downloaded
2023-06-02 00:16:40.215 | INFO     | plassembler:run:338 - Database successfully checked.
2023-06-02 00:16:40.215 | INFO     | plassembler:run:341 - Checking input fastqs.
2023-06-02 00:16:40.230 | INFO     | plassembler.utils.input_commands:validate_fastq:24 - FASTQ SRR22859710_1.fastq.gz checked
2023-06-02 00:16:40.230 | INFO     | plassembler.utils.input_commands:validate_fastq:24 - FASTQ SRR22859826_1.fastq.gz checked
2023-06-02 00:16:40.231 | INFO     | plassembler.utils.input_commands:validate_fastq:24 - FASTQ SRR22859826_2.fastq.gz checked
2023-06-02 00:16:40.231 | INFO     | plassembler:run:349 - FASTQ file SRR22859710_1.fastq.gz compression is True
2023-06-02 00:16:40.231 | INFO     | plassembler:run:350 - FASTQ file SRR22859826_1.fastq.gz compression is True
2023-06-02 00:16:40.231 | INFO     | plassembler:run:351 - FASTQ file SRR22859826_2.fastq.gz compression is True
2023-06-02 00:16:40.231 | INFO     | plassembler:run:359 - Filtering long reads with chopper
2023-06-02 00:16:40.231 | INFO     | plassembler.utils.qc:chopper:25 - Started running chopper
plassembler.py -d ../Plassembler_DB -l SRR22859710_1.fastq.gz  -1  SRR22859924_1.fastq.gz -2 SRR22859924_2.fastq.gz   -o S_Aureus_QC_Illustration_different_sequence_type -t 8 -f -c 2500000
2023-06-02 00:17:26.134 | INFO     | plassembler.utils.qc:chopper:82 - Finished running chopper
2023-06-02 00:17:26.135 | INFO     | plassembler:run:378 - Running Flye.
2023-06-02 00:17:26.137 | INFO     | plassembler.utils.external_tools:run:55 - Started running flye --nano-hq S_Aureus_QC_Illustration_same_sequence_type/chopper_long_reads.fastq.gz --out-dir S_Aureus_QC_Illustration_same_sequence_type --threads 8 ...
2023-06-02 00:23:44.650 | INFO     | plassembler.utils.external_tools:run:57 - Done running flye --nano-hq S_Aureus_QC_Illustration_same_sequence_type/chopper_long_reads.fastq.gz --out-dir S_Aureus_QC_Illustration_same_sequence_type --threads 8
2023-06-02 00:23:44.652 | INFO     | plassembler:run:388 - Counting Contigs.
2023-06-02 00:23:44.668 | INFO     | plassembler.utils.plass_class:get_contig_count:83 - Assembled 2 contigs.
2023-06-02 00:23:44.669 | INFO     | plassembler:run:570 - More than one contig was assembled with Flye.
2023-06-02 00:23:44.669 | INFO     | plassembler:run:571 - Extracting Chromosome.
2023-06-02 00:23:44.718 | INFO     | plassembler:run:605 - Chromosome Identified. Plassembler will now use long and short reads to assemble plasmids accurately.
2023-06-02 00:23:44.718 | INFO     | plassembler:run:607 - Mapping long reads.
2023-06-02 00:23:44.719 | INFO     | plassembler.utils.external_tools:run_to_stdout:64 - Started running minimap2 -ax map-ont -t 8 S_Aureus_QC_Illustration_same_sequence_type/flye_renamed.fasta S_Aureus_QC_Illustration_same_sequence_type/chopper_long_reads.fastq.gz ...
2023-06-02 00:23:57.652 | INFO     | plassembler.utils.external_tools:run_to_stdout:66 - Done running minimap2 -ax map-ont -t 8 S_Aureus_QC_Illustration_same_sequence_type/flye_renamed.fasta S_Aureus_QC_Illustration_same_sequence_type/chopper_long_reads.fastq.gz
2023-06-02 00:23:57.653 | INFO     | plassembler:run:622 - Trimming short reads.
2023-06-02 00:23:57.654 | INFO     | plassembler.utils.external_tools:run:55 - Started running fastp --out1 S_Aureus_QC_Illustration_same_sequence_type/trimmed_R1.fastq --out2 S_Aureus_QC_Illustration_same_sequence_type/trimmed_R2.fastq --in1 SRR22859826_1.fastq.gz --in2 SRR22859826_2.fastq.gz ...
2023-06-02 00:24:02.913 | INFO     | plassembler.utils.external_tools:run:57 - Done running fastp --out1 S_Aureus_QC_Illustration_same_sequence_type/trimmed_R1.fastq --out2 S_Aureus_QC_Illustration_same_sequence_type/trimmed_R2.fastq --in1 SRR22859826_1.fastq.gz --in2 SRR22859826_2.fastq.gz
2023-06-02 00:24:02.914 | INFO     | plassembler:run:626 - Mapping short reads.
2023-06-02 00:24:02.915 | INFO     | plassembler.utils.external_tools:run_to_stdout:64 - Started running minimap2 -ax sr -t 8 S_Aureus_QC_Illustration_same_sequence_type/flye_renamed.fasta S_Aureus_QC_Illustration_same_sequence_type/trimmed_R1.fastq S_Aureus_QC_Illustration_same_sequence_type/trimmed_R2.fastq ...
2023-06-02 00:24:09.225 | INFO     | plassembler.utils.external_tools:run_to_stdout:66 - Done running minimap2 -ax sr -t 8 S_Aureus_QC_Illustration_same_sequence_type/flye_renamed.fasta S_Aureus_QC_Illustration_same_sequence_type/trimmed_R1.fastq S_Aureus_QC_Illustration_same_sequence_type/trimmed_R2.fastq
2023-06-02 00:24:09.226 | INFO     | plassembler:run:633 - Processing Sam/Bam Files and extracting Fastqs.
2023-06-02 00:24:56.524 | INFO     | plassembler.utils.external_tools:run_to_stdout:64 - Started running samtools view -h -@ 8 -b S_Aureus_QC_Illustration_same_sequence_type/short_read.sam ...
2023-06-02 00:24:58.671 | INFO     | plassembler.utils.external_tools:run_to_stdout:66 - Done running samtools view -h -@ 8 -b S_Aureus_QC_Illustration_same_sequence_type/short_read.sam
2023-06-02 00:24:58.673 | INFO     | plassembler.utils.external_tools:run_to_stdout:64 - Started running samtools view -b -h -@ 8 -L S_Aureus_QC_Illustration_same_sequence_type/non_chromosome.bed S_Aureus_QC_Illustration_same_sequence_type/short_read.bam ...
2023-06-02 00:24:58.976 | INFO     | plassembler.utils.external_tools:run_to_stdout:66 - Done running samtools view -b -h -@ 8 -L S_Aureus_QC_Illustration_same_sequence_type/non_chromosome.bed S_Aureus_QC_Illustration_same_sequence_type/short_read.bam
2023-06-02 00:24:58.978 | INFO     | plassembler.utils.external_tools:run_to_stdout:64 - Started running samtools view -b -h -f 4 -@ 8 S_Aureus_QC_Illustration_same_sequence_type/short_read.bam ...
2023-06-02 00:24:59.268 | INFO     | plassembler.utils.external_tools:run_to_stdout:66 - Done running samtools view -b -h -f 4 -@ 8 S_Aureus_QC_Illustration_same_sequence_type/short_read.bam
2023-06-02 00:24:59.270 | INFO     | plassembler.utils.external_tools:run_to_stdout:64 - Started running samtools view -b -h -@ 8 -L S_Aureus_QC_Illustration_same_sequence_type/chromosome.bed S_Aureus_QC_Illustration_same_sequence_type/short_read.bam ...
2023-06-02 00:25:01.445 | INFO     | plassembler.utils.external_tools:run_to_stdout:66 - Done running samtools view -b -h -@ 8 -L S_Aureus_QC_Illustration_same_sequence_type/chromosome.bed S_Aureus_QC_Illustration_same_sequence_type/short_read.bam
2023-06-02 00:25:01.446 | INFO     | plassembler.utils.external_tools:run:55 - Started running samtools fastq -@ 8 S_Aureus_QC_Illustration_same_sequence_type/unmapped_bam_file.bam -1 S_Aureus_QC_Illustration_same_sequence_type/unmapped_R1.fastq -2 S_Aureus_QC_Illustration_same_sequence_type/unmapped_R2.fastq -0 /dev/null -s /dev/null -n ...
2023-06-02 00:25:01.549 | INFO     | plassembler.utils.external_tools:run:57 - Done running samtools fastq -@ 8 S_Aureus_QC_Illustration_same_sequence_type/unmapped_bam_file.bam -1 S_Aureus_QC_Illustration_same_sequence_type/unmapped_R1.fastq -2 S_Aureus_QC_Illustration_same_sequence_type/unmapped_R2.fastq -0 /dev/null -s /dev/null -n
2023-06-02 00:25:01.551 | INFO     | plassembler.utils.external_tools:run:55 - Started running samtools fastq -@ 8 S_Aureus_QC_Illustration_same_sequence_type/non_chromosome.bam -1 S_Aureus_QC_Illustration_same_sequence_type/mapped_non_chromosome_R1.fastq -2 S_Aureus_QC_Illustration_same_sequence_type/mapped_non_chromosome_R2.fastq -0 /dev/null -s /dev/null -n ...
2023-06-02 00:25:01.637 | INFO     | plassembler.utils.external_tools:run:57 - Done running samtools fastq -@ 8 S_Aureus_QC_Illustration_same_sequence_type/non_chromosome.bam -1 S_Aureus_QC_Illustration_same_sequence_type/mapped_non_chromosome_R1.fastq -2 S_Aureus_QC_Illustration_same_sequence_type/mapped_non_chromosome_R2.fastq -0 /dev/null -s /dev/null -n
2023-06-02 00:25:03.041 | INFO     | plassembler:run:645 - Running Unicycler.
2023-06-02 00:25:03.042 | INFO     | plassembler.utils.external_tools:run:55 - Started running unicycler -1 S_Aureus_QC_Illustration_same_sequence_type/short_read_concat_R1.fastq -2 S_Aureus_QC_Illustration_same_sequence_type/short_read_concat_R2.fastq -l S_Aureus_QC_Illustration_same_sequence_type/plasmid_long.fastq -t 8 -o S_Aureus_QC_Illustration_same_sequence_type/unicycler_output ...
2023-06-02 00:32:30.663 | INFO     | plassembler.utils.external_tools:run:57 - Done running unicycler -1 S_Aureus_QC_Illustration_same_sequence_type/short_read_concat_R1.fastq -2 S_Aureus_QC_Illustration_same_sequence_type/short_read_concat_R2.fastq -l S_Aureus_QC_Illustration_same_sequence_type/plasmid_long.fastq -t 8 -o S_Aureus_QC_Illustration_same_sequence_type/unicycler_output
2023-06-02 00:32:30.664 | INFO     | plassembler:run:662 - Unicycler identified plasmids. Calculating Plasmid Copy Numbers.
2023-06-02 00:32:30.697 | INFO     | plassembler.utils.external_tools:run_to_stdout:64 - Started running minimap2 -ax map-ont -t 8 S_Aureus_QC_Illustration_same_sequence_type/combined.fasta S_Aureus_QC_Illustration_same_sequence_type/chopper_long_reads.fastq.gz ...
2023-06-02 00:32:43.768 | INFO     | plassembler.utils.external_tools:run_to_stdout:66 - Done running minimap2 -ax map-ont -t 8 S_Aureus_QC_Illustration_same_sequence_type/combined.fasta S_Aureus_QC_Illustration_same_sequence_type/chopper_long_reads.fastq.gz
2023-06-02 00:32:43.770 | INFO     | plassembler.utils.external_tools:run:55 - Started running samtools sort -@ 8 S_Aureus_QC_Illustration_same_sequence_type/combined_long.sam -o S_Aureus_QC_Illustration_same_sequence_type/combined_sorted_long.bam ...
2023-06-02 00:32:44.998 | INFO     | plassembler.utils.external_tools:run:57 - Done running samtools sort -@ 8 S_Aureus_QC_Illustration_same_sequence_type/combined_long.sam -o S_Aureus_QC_Illustration_same_sequence_type/combined_sorted_long.bam
2023-06-02 00:32:45.000 | INFO     | plassembler.utils.external_tools:run_to_stdout:64 - Started running minimap2 -ax sr -t 8 S_Aureus_QC_Illustration_same_sequence_type/combined.fasta S_Aureus_QC_Illustration_same_sequence_type/trimmed_R1.fastq S_Aureus_QC_Illustration_same_sequence_type/trimmed_R2.fastq ...
2023-06-02 00:32:51.731 | INFO     | plassembler.utils.external_tools:run_to_stdout:66 - Done running minimap2 -ax sr -t 8 S_Aureus_QC_Illustration_same_sequence_type/combined.fasta S_Aureus_QC_Illustration_same_sequence_type/trimmed_R1.fastq S_Aureus_QC_Illustration_same_sequence_type/trimmed_R2.fastq
2023-06-02 00:32:51.733 | INFO     | plassembler.utils.external_tools:run:55 - Started running samtools sort -@ 8 S_Aureus_QC_Illustration_same_sequence_type/combined_short.sam -o S_Aureus_QC_Illustration_same_sequence_type/combined_sorted_short.bam ...
2023-06-02 00:32:54.312 | INFO     | plassembler.utils.external_tools:run:57 - Done running samtools sort -@ 8 S_Aureus_QC_Illustration_same_sequence_type/combined_short.sam -o S_Aureus_QC_Illustration_same_sequence_type/combined_sorted_short.bam
2023-06-02 00:33:11.514 | INFO     | plassembler:run:670 - Calculating mash distances to PLSDB.
2023-06-02 00:33:11.518 | INFO     | plassembler.utils.external_tools:run:55 - Started running mash sketch S_Aureus_QC_Illustration_same_sequence_type/plasmids.fasta -i ...
2023-06-02 00:33:11.876 | INFO     | plassembler.utils.external_tools:run:57 - Done running mash sketch S_Aureus_QC_Illustration_same_sequence_type/plasmids.fasta -i
2023-06-02 00:33:11.878 | INFO     | plassembler.utils.external_tools:run_to_stdout:64 - Started running mash dist S_Aureus_QC_Illustration_same_sequence_type/plasmids.fasta.msh ../Plassembler_DB/plsdb.msh -v 0.1 -d 0.1 -i ...
2023-06-02 00:33:17.776 | INFO     | plassembler.utils.external_tools:run_to_stdout:66 - Done running mash dist S_Aureus_QC_Illustration_same_sequence_type/plasmids.fasta.msh ../Plassembler_DB/plsdb.msh -v 0.1 -d 0.1 -i
2023-06-02 00:33:18.174 | INFO     | plassembler.utils.test_incompatibility:incompatbility:24 - WARNING: 8 non-circular contigs with no PLSDB mash hits were detected. 
This indicates your long and short read sets may come from different bacterial isolates. 
Please check this!
2023-06-02 00:33:18.316 | INFO     | plassembler:end_plassembler:90 - Plassembler has finished
2023-06-02 00:33:18.316 | INFO     | plassembler:end_plassembler:91 - Elapsed time: 999.6 seconds
```

# Different Sequence Type

For the different sequence type isolates, plassembler took double as long as with the same sequence type isolates, which is not surprising as there would have been many more unmapped reads to assemble with Unicycler.

As you can see below, Plassembler recovered 43 non-circular contigs with no hits to PLSDB, indictating that the long and short read sets don't match.

```
plassembler run -d plassembler_db -l SRR22859710_1.fastq.gz  -1  SRR22859924_1.fastq.gz -2 SRR22859924_2.fastq.gz   -o S_Aureus_QC_Illustration_different_sequence_type -t 8 -f -c 2500000
```

```
2023-06-02 09:53:04.687 | INFO     | plassembler:begin_plassembler:73 - You are using Plassembler version 1.1.0
2023-06-02 09:53:04.687 | INFO     | plassembler:begin_plassembler:74 - Repository homepage is https://github.com/gbouras13/plassembler
2023-06-02 09:53:04.687 | INFO     | plassembler:begin_plassembler:75 - Written by George Bouras: george.bouras@adelaide.edu.au
2023-06-02 09:53:04.687 | INFO     | plassembler:run:312 - Database directory is ../Plassembler_DB
2023-06-02 09:53:04.687 | INFO     | plassembler:run:313 - Longreads file is SRR22859710_1.fastq.gz
2023-06-02 09:53:04.688 | INFO     | plassembler:run:314 - R1 fasta file is SRR22859924_1.fastq.gz
2023-06-02 09:53:04.688 | INFO     | plassembler:run:315 - R2 fasta file is SRR22859924_2.fastq.gz
2023-06-02 09:53:04.688 | INFO     | plassembler:run:316 - Chromosome length threshold is 2500000
2023-06-02 09:53:04.688 | INFO     | plassembler:run:317 - Output directory is S_Aureus_QC_Illustration_different_sequence_type
2023-06-02 09:53:04.688 | INFO     | plassembler:run:318 - Min long read length is 500
2023-06-02 09:53:04.688 | INFO     | plassembler:run:319 - Min long read quality is 9
2023-06-02 09:53:04.688 | INFO     | plassembler:run:320 - Thread count is 8
2023-06-02 09:53:04.688 | INFO     | plassembler:run:321 - --force is True
2023-06-02 09:53:04.688 | INFO     | plassembler:run:322 - --skip_qc is False
2023-06-02 09:53:04.688 | INFO     | plassembler:run:323 - --raw_flag is False
2023-06-02 09:53:04.688 | INFO     | plassembler:run:324 - --pacbio_model is nothing
2023-06-02 09:53:04.688 | INFO     | plassembler:run:325 - --keep_fastqs is False
2023-06-02 09:53:04.688 | INFO     | plassembler:run:326 - --keep_chromosome is False
2023-06-02 09:53:04.689 | INFO     | plassembler:run:330 - Checking dependencies
2023-06-02 09:53:04.942 | INFO     | plassembler.utils.input_commands:check_dependencies:140 - Flye version found is v2.9.2-b1786.
2023-06-02 09:53:04.943 | INFO     | plassembler.utils.input_commands:check_dependencies:150 - Flye version is ok.
2023-06-02 09:53:05.150 | INFO     | plassembler.utils.input_commands:check_dependencies:159 - Raven v1.8.1 found.
2023-06-02 09:53:05.151 | INFO     | plassembler.utils.input_commands:check_dependencies:161 - Raven version is ok.
2023-06-02 09:53:05.625 | INFO     | plassembler.utils.input_commands:check_dependencies:190 - Unicycler version found is v0.5.0.
2023-06-02 09:53:05.625 | INFO     | plassembler.utils.input_commands:check_dependencies:203 - Unicycler version is ok.
2023-06-02 09:53:06.069 | INFO     | plassembler.utils.input_commands:check_dependencies:213 - SPAdes v3.15.2 found.
2023-06-02 09:53:06.814 | INFO     | plassembler.utils.input_commands:check_dependencies:226 - Samtools v1.17 found.
2023-06-02 09:53:06.902 | INFO     | plassembler.utils.input_commands:check_dependencies:237 - minimap2 v2.26-r1175 found.
2023-06-02 09:53:07.039 | INFO     | plassembler.utils.input_commands:check_dependencies:248 - fastp v0.23.4 found.
2023-06-02 09:53:07.166 | INFO     | plassembler.utils.input_commands:check_dependencies:259 - chopper v0.5.0 found.
2023-06-02 09:53:07.562 | INFO     | plassembler.utils.input_commands:check_dependencies:274 - mash v2.3 found.
2023-06-02 09:53:07.563 | INFO     | plassembler.utils.input_commands:check_dependencies:279 - All dependencies found.
2023-06-02 09:53:07.563 | INFO     | plassembler:run:335 - Checking database installation.
2023-06-02 09:53:07.563 | INFO     | plassembler.utils.db:check_db_installation:24 - PLSDB Database at ../Plassembler_DB has already been downloaded
2023-06-02 09:53:07.563 | INFO     | plassembler:run:338 - Database successfully checked.
2023-06-02 09:53:07.563 | INFO     | plassembler:run:341 - Checking input fastqs.
2023-06-02 09:53:07.582 | INFO     | plassembler.utils.input_commands:validate_fastq:24 - FASTQ SRR22859710_1.fastq.gz checked
2023-06-02 09:53:07.583 | INFO     | plassembler.utils.input_commands:validate_fastq:24 - FASTQ SRR22859924_1.fastq.gz checked
2023-06-02 09:53:07.584 | INFO     | plassembler.utils.input_commands:validate_fastq:24 - FASTQ SRR22859924_2.fastq.gz checked
2023-06-02 09:53:07.584 | INFO     | plassembler:run:349 - FASTQ file SRR22859710_1.fastq.gz compression is True
2023-06-02 09:53:07.584 | INFO     | plassembler:run:350 - FASTQ file SRR22859924_1.fastq.gz compression is True
2023-06-02 09:53:07.585 | INFO     | plassembler:run:351 - FASTQ file SRR22859924_2.fastq.gz compression is True
2023-06-02 09:53:07.585 | INFO     | plassembler:run:359 - Filtering long reads with chopper
2023-06-02 09:53:07.585 | INFO     | plassembler.utils.qc:chopper:25 - Started running chopper
2023-06-02 09:53:53.573 | INFO     | plassembler.utils.qc:chopper:82 - Finished running chopper
2023-06-02 09:53:53.575 | INFO     | plassembler:run:378 - Running Flye.
2023-06-02 09:53:53.577 | INFO     | plassembler.utils.external_tools:run:55 - Started running flye --nano-hq S_Aureus_QC_Illustration_different_sequence_type/chopper_long_reads.fastq.gz --out-dir S_Aureus_QC_Illustration_different_sequence_type --threads 8 ...
2023-06-02 10:00:51.040 | INFO     | plassembler.utils.external_tools:run:57 - Done running flye --nano-hq S_Aureus_QC_Illustration_different_sequence_type/chopper_long_reads.fastq.gz --out-dir S_Aureus_QC_Illustration_different_sequence_type --threads 8
2023-06-02 10:00:51.042 | INFO     | plassembler:run:388 - Counting Contigs.
2023-06-02 10:00:51.061 | INFO     | plassembler.utils.plass_class:get_contig_count:83 - Assembled 2 contigs.
2023-06-02 10:00:51.065 | INFO     | plassembler:run:570 - More than one contig was assembled with Flye.
2023-06-02 10:00:51.065 | INFO     | plassembler:run:571 - Extracting Chromosome.
2023-06-02 10:00:51.149 | INFO     | plassembler:run:605 - Chromosome Identified. Plassembler will now use long and short reads to assemble plasmids accurately.
2023-06-02 10:00:51.150 | INFO     | plassembler:run:607 - Mapping long reads.
2023-06-02 10:00:51.154 | INFO     | plassembler.utils.external_tools:run_to_stdout:64 - Started running minimap2 -ax map-ont -t 8 S_Aureus_QC_Illustration_different_sequence_type/flye_renamed.fasta S_Aureus_QC_Illustration_different_sequence_type/chopper_long_reads.fastq.gz ...
2023-06-02 10:01:09.111 | INFO     | plassembler.utils.external_tools:run_to_stdout:66 - Done running minimap2 -ax map-ont -t 8 S_Aureus_QC_Illustration_different_sequence_type/flye_renamed.fasta S_Aureus_QC_Illustration_different_sequence_type/chopper_long_reads.fastq.gz
2023-06-02 10:01:09.113 | INFO     | plassembler:run:622 - Trimming short reads.
2023-06-02 10:01:09.115 | INFO     | plassembler.utils.external_tools:run:55 - Started running fastp --out1 S_Aureus_QC_Illustration_different_sequence_type/trimmed_R1.fastq --out2 S_Aureus_QC_Illustration_different_sequence_type/trimmed_R2.fastq --in1 SRR22859924_1.fastq.gz --in2 SRR22859924_2.fastq.gz ...
2023-06-02 10:01:23.234 | INFO     | plassembler.utils.external_tools:run:57 - Done running fastp --out1 S_Aureus_QC_Illustration_different_sequence_type/trimmed_R1.fastq --out2 S_Aureus_QC_Illustration_different_sequence_type/trimmed_R2.fastq --in1 SRR22859924_1.fastq.gz --in2 SRR22859924_2.fastq.gz
2023-06-02 10:01:23.236 | INFO     | plassembler:run:626 - Mapping short reads.
2023-06-02 10:01:23.237 | INFO     | plassembler.utils.external_tools:run_to_stdout:64 - Started running minimap2 -ax sr -t 8 S_Aureus_QC_Illustration_different_sequence_type/flye_renamed.fasta S_Aureus_QC_Illustration_different_sequence_type/trimmed_R1.fastq S_Aureus_QC_Illustration_different_sequence_type/trimmed_R2.fastq ...
2023-06-02 10:01:36.239 | INFO     | plassembler.utils.external_tools:run_to_stdout:66 - Done running minimap2 -ax sr -t 8 S_Aureus_QC_Illustration_different_sequence_type/flye_renamed.fasta S_Aureus_QC_Illustration_different_sequence_type/trimmed_R1.fastq S_Aureus_QC_Illustration_different_sequence_type/trimmed_R2.fastq
2023-06-02 10:01:36.240 | INFO     | plassembler:run:633 - Processing Sam/Bam Files and extracting Fastqs.
2023-06-02 10:02:26.065 | INFO     | plassembler.utils.external_tools:run_to_stdout:64 - Started running samtools view -h -@ 8 -b S_Aureus_QC_Illustration_different_sequence_type/short_read.sam ...
2023-06-02 10:02:30.618 | INFO     | plassembler.utils.external_tools:run_to_stdout:66 - Done running samtools view -h -@ 8 -b S_Aureus_QC_Illustration_different_sequence_type/short_read.sam
2023-06-02 10:02:30.620 | INFO     | plassembler.utils.external_tools:run_to_stdout:64 - Started running samtools view -b -h -@ 8 -L S_Aureus_QC_Illustration_different_sequence_type/non_chromosome.bed S_Aureus_QC_Illustration_different_sequence_type/short_read.bam ...
2023-06-02 10:02:31.107 | INFO     | plassembler.utils.external_tools:run_to_stdout:66 - Done running samtools view -b -h -@ 8 -L S_Aureus_QC_Illustration_different_sequence_type/non_chromosome.bed S_Aureus_QC_Illustration_different_sequence_type/short_read.bam
2023-06-02 10:02:31.109 | INFO     | plassembler.utils.external_tools:run_to_stdout:64 - Started running samtools view -b -h -f 4 -@ 8 S_Aureus_QC_Illustration_different_sequence_type/short_read.bam ...
2023-06-02 10:02:32.175 | INFO     | plassembler.utils.external_tools:run_to_stdout:66 - Done running samtools view -b -h -f 4 -@ 8 S_Aureus_QC_Illustration_different_sequence_type/short_read.bam
2023-06-02 10:02:32.177 | INFO     | plassembler.utils.external_tools:run_to_stdout:64 - Started running samtools view -b -h -@ 8 -L S_Aureus_QC_Illustration_different_sequence_type/chromosome.bed S_Aureus_QC_Illustration_different_sequence_type/short_read.bam ...
2023-06-02 10:02:36.400 | INFO     | plassembler.utils.external_tools:run_to_stdout:66 - Done running samtools view -b -h -@ 8 -L S_Aureus_QC_Illustration_different_sequence_type/chromosome.bed S_Aureus_QC_Illustration_different_sequence_type/short_read.bam
2023-06-02 10:02:36.402 | INFO     | plassembler.utils.external_tools:run:55 - Started running samtools fastq -@ 8 S_Aureus_QC_Illustration_different_sequence_type/unmapped_bam_file.bam -1 S_Aureus_QC_Illustration_different_sequence_type/unmapped_R1.fastq -2 S_Aureus_QC_Illustration_different_sequence_type/unmapped_R2.fastq -0 /dev/null -s /dev/null -n ...
2023-06-02 10:02:37.640 | INFO     | plassembler.utils.external_tools:run:57 - Done running samtools fastq -@ 8 S_Aureus_QC_Illustration_different_sequence_type/unmapped_bam_file.bam -1 S_Aureus_QC_Illustration_different_sequence_type/unmapped_R1.fastq -2 S_Aureus_QC_Illustration_different_sequence_type/unmapped_R2.fastq -0 /dev/null -s /dev/null -n
2023-06-02 10:02:37.642 | INFO     | plassembler.utils.external_tools:run:55 - Started running samtools fastq -@ 8 S_Aureus_QC_Illustration_different_sequence_type/non_chromosome.bam -1 S_Aureus_QC_Illustration_different_sequence_type/mapped_non_chromosome_R1.fastq -2 S_Aureus_QC_Illustration_different_sequence_type/mapped_non_chromosome_R2.fastq -0 /dev/null -s /dev/null -n ...
2023-06-02 10:02:37.697 | INFO     | plassembler.utils.external_tools:run:57 - Done running samtools fastq -@ 8 S_Aureus_QC_Illustration_different_sequence_type/non_chromosome.bam -1 S_Aureus_QC_Illustration_different_sequence_type/mapped_non_chromosome_R1.fastq -2 S_Aureus_QC_Illustration_different_sequence_type/mapped_non_chromosome_R2.fastq -0 /dev/null -s /dev/null -n
2023-06-02 10:02:53.611 | INFO     | plassembler:run:645 - Running Unicycler.
2023-06-02 10:02:53.613 | INFO     | plassembler.utils.external_tools:run:55 - Started running unicycler -1 S_Aureus_QC_Illustration_different_sequence_type/short_read_concat_R1.fastq -2 S_Aureus_QC_Illustration_different_sequence_type/short_read_concat_R2.fastq -l S_Aureus_QC_Illustration_different_sequence_type/plasmid_long.fastq -t 8 -o S_Aureus_QC_Illustration_different_sequence_type/unicycler_output ...
2023-06-02 10:29:27.975 | INFO     | plassembler.utils.external_tools:run:57 - Done running unicycler -1 S_Aureus_QC_Illustration_different_sequence_type/short_read_concat_R1.fastq -2 S_Aureus_QC_Illustration_different_sequence_type/short_read_concat_R2.fastq -l S_Aureus_QC_Illustration_different_sequence_type/plasmid_long.fastq -t 8 -o S_Aureus_QC_Illustration_different_sequence_type/unicycler_output
2023-06-02 10:29:27.978 | INFO     | plassembler:run:662 - Unicycler identified plasmids. Calculating Plasmid Copy Numbers.
2023-06-02 10:29:28.022 | INFO     | plassembler.utils.external_tools:run_to_stdout:64 - Started running minimap2 -ax map-ont -t 8 S_Aureus_QC_Illustration_different_sequence_type/combined.fasta S_Aureus_QC_Illustration_different_sequence_type/chopper_long_reads.fastq.gz ...
2023-06-02 10:29:48.775 | INFO     | plassembler.utils.external_tools:run_to_stdout:66 - Done running minimap2 -ax map-ont -t 8 S_Aureus_QC_Illustration_different_sequence_type/combined.fasta S_Aureus_QC_Illustration_different_sequence_type/chopper_long_reads.fastq.gz
2023-06-02 10:29:48.784 | INFO     | plassembler.utils.external_tools:run:55 - Started running samtools sort -@ 8 S_Aureus_QC_Illustration_different_sequence_type/combined_long.sam -o S_Aureus_QC_Illustration_different_sequence_type/combined_sorted_long.bam ...
2023-06-02 10:29:50.769 | INFO     | plassembler.utils.external_tools:run:57 - Done running samtools sort -@ 8 S_Aureus_QC_Illustration_different_sequence_type/combined_long.sam -o S_Aureus_QC_Illustration_different_sequence_type/combined_sorted_long.bam
2023-06-02 10:29:50.771 | INFO     | plassembler.utils.external_tools:run_to_stdout:64 - Started running minimap2 -ax sr -t 8 S_Aureus_QC_Illustration_different_sequence_type/combined.fasta S_Aureus_QC_Illustration_different_sequence_type/trimmed_R1.fastq S_Aureus_QC_Illustration_different_sequence_type/trimmed_R2.fastq ...
2023-06-02 10:30:10.966 | INFO     | plassembler.utils.external_tools:run_to_stdout:66 - Done running minimap2 -ax sr -t 8 S_Aureus_QC_Illustration_different_sequence_type/combined.fasta S_Aureus_QC_Illustration_different_sequence_type/trimmed_R1.fastq S_Aureus_QC_Illustration_different_sequence_type/trimmed_R2.fastq
2023-06-02 10:30:10.970 | INFO     | plassembler.utils.external_tools:run:55 - Started running samtools sort -@ 8 S_Aureus_QC_Illustration_different_sequence_type/combined_short.sam -o S_Aureus_QC_Illustration_different_sequence_type/combined_sorted_short.bam ...
2023-06-02 10:30:18.225 | INFO     | plassembler.utils.external_tools:run:57 - Done running samtools sort -@ 8 S_Aureus_QC_Illustration_different_sequence_type/combined_short.sam -o S_Aureus_QC_Illustration_different_sequence_type/combined_sorted_short.bam
2023-06-02 10:30:37.595 | INFO     | plassembler:run:670 - Calculating mash distances to PLSDB.
2023-06-02 10:30:37.603 | INFO     | plassembler.utils.external_tools:run:55 - Started running mash sketch S_Aureus_QC_Illustration_different_sequence_type/plasmids.fasta -i ...
2023-06-02 10:30:37.765 | INFO     | plassembler.utils.external_tools:run:57 - Done running mash sketch S_Aureus_QC_Illustration_different_sequence_type/plasmids.fasta -i
2023-06-02 10:30:37.767 | INFO     | plassembler.utils.external_tools:run_to_stdout:64 - Started running mash dist S_Aureus_QC_Illustration_different_sequence_type/plasmids.fasta.msh ../Plassembler_DB/plsdb.msh -v 0.1 -d 0.1 -i ...
2023-06-02 10:31:04.301 | INFO     | plassembler.utils.external_tools:run_to_stdout:66 - Done running mash dist S_Aureus_QC_Illustration_different_sequence_type/plasmids.fasta.msh ../Plassembler_DB/plsdb.msh -v 0.1 -d 0.1 -i
2023-06-02 10:31:04.857 | INFO     | plassembler.utils.test_incompatibility:incompatbility:24 - WARNING: 43 non-circular contigs with no PLSDB mash hits were detected. 
This indicates your long and short read sets may come from different bacterial isolates. 
Please check this!
2023-06-02 10:31:05.971 | INFO     | plassembler:end_plassembler:90 - Plassembler has finished
2023-06-02 10:31:05.971 | INFO     | plassembler:end_plassembler:91 - Elapsed time: 2281.32 seconds
```





