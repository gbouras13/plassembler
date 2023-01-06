##### Plassembler Benchmarking v0.1.4

* Firstly, I downloaded the 6 tarballs from this [link](https://bridges.monash.edu/articles/dataset/Small_plasmid_Nanopore_data/13543754). 

* These contained tech_rep_1_rapid_reads.fastq.gz, tech_rep_2_rapid_reads.fastq.gz, tech_rep_1_ligation_reads.fastq.gz, tech_rep_2_ligation_reads.fastq.gz, tech_rep_1_illumina_reads  and tech_rep_2_illumina_reads tarballs.

These are the barcodes linked to the species

* barcode01 == Acinetobacter_baumannii_J9
* barcode02 == Citrobacter_koseri_MINF_9D
* barcode03 ==  Enterobacter_kobei_MSB1_1B
* barcode04 == Haemophilus_unknown_M1C132_1
* barcode05 == Klebsiella_oxytoca_MSB1_2C
* barcode07 == Klebsiella_variicola_INF345
 
**Serratia Marscecens**

* barcode08 Serratia Marcescens was excluded from benchmarking due to how devilishly difficult it was for Ryan Wick to assemble it, even with Trycycler - see his [comments](https://github.com/rrwick/Small-plasmid-Nanopore/blob/main/method.md#Long-read). 
* I tried assembling it in different ways with different read subsets, but kept getting different results every time depending on what read subset I chose or even how I randomly sampled the reads, indicating there was enormous heterogeneity in this sample. Even when I used all reads (over 400x chromosome coverage), the Flye step in Plassembler could not resolve the 184447bp plasmid which was concerning! So I excluded it from the benchmarking.


#### Recovering Barcodes 

* These files contained the reads for all barcodes pooled together. The barcode was indicated in the FASTQ header. 
* Accordingly, I split the files using grep as follows - run from within the directory of each unzipped tarball:

```
### tech rep 1
# rapid
gunzip -c tech_rep_1_rapid_reads.fastq.gz | grep -A3 "barcode=barcode01"  | gzip >  barcode01.fastq.gz
gunzip -c tech_rep_1_rapid_reads.fastq.gz | grep -A3 "barcode=barcode02"   | gzip >  barcode02.fastq.gz
gunzip -c tech_rep_1_rapid_reads.fastq.gz | grep -A3 "barcode=barcode03"  | gzip >  barcode03.fastq.gz
gunzip -c tech_rep_1_rapid_reads.fastq.gz | grep -A3 "barcode=barcode04"  | gzip >  barcode04.fastq.gz
gunzip -c tech_rep_1_rapid_reads.fastq.gz | grep -A3 "barcode=barcode05" | gzip >  barcode05.fastq.gz
gunzip -c tech_rep_1_rapid_reads.fastq.gz | grep -A3 "barcode=barcode07" | gzip >  barcode07.fastq.gz


# ligation
gunzip -c tech_rep_1_ligation_reads.fastq.gz | grep -A3 "barcode=barcode01"  | gzip >  barcode01.fastq.gz
gunzip -c tech_rep_1_ligation_reads.fastq.gz | grep -A3 "barcode=barcode02"   | gzip >  barcode02.fastq.gz
gunzip -c tech_rep_1_ligation_reads.fastq.gz | grep -A3 "barcode=barcode03"  | gzip >  barcode03.fastq.gz
gunzip -c tech_rep_1_ligation_reads.fastq.gz | grep -A3 "barcode=barcode04"  | gzip >  barcode04.fastq.gz
gunzip -c tech_rep_1_ligation_reads.fastq.gz | grep -A3 "barcode=barcode05" | gzip >  barcode05.fastq.gz
gunzip -c tech_rep_1_ligation_reads.fastq.gz | grep -A3 "barcode=barcode07" | gzip >  barcode07.fastq.gz


### tech rep 2
# rapid
gunzip -c tech_rep_2_rapid_reads.fastq.gz | grep -A3 "barcode=barcode01" | gzip >  barcode01.fastq.gz
gunzip -c tech_rep_2_rapid_reads.fastq.gz | grep -A3 "barcode=barcode02" | gzip >  barcode02.fastq.gz
gunzip -c tech_rep_2_rapid_reads.fastq.gz | grep -A3 "barcode=barcode03" | gzip >  barcode03.fastq.gz
gunzip -c tech_rep_2_rapid_reads.fastq.gz | grep -A3 "barcode=barcode04" | gzip >  barcode04.fastq.gz
gunzip -c tech_rep_2_rapid_reads.fastq.gz | grep -A3 "barcode=barcode05" | gzip >  barcode05.fastq.gz
gunzip -c tech_rep_2_rapid_reads.fastq.gz | grep -A3 "barcode=barcode07" | gzip >  barcode07.fastq.gz


# ligation
gunzip -c tech_rep_2_ligation_reads.fastq.gz | grep -A3 "barcode=barcode01" | gzip >  barcode01.fastq.gz
gunzip -c tech_rep_2_ligation_reads.fastq.gz | grep -A3 "barcode=barcode02" | gzip >  barcode02.fastq.gz
gunzip -c tech_rep_2_ligation_reads.fastq.gz | grep -A3 "barcode=barcode03" | gzip >  barcode03.fastq.gz
gunzip -c tech_rep_2_ligation_reads.fastq.gz | grep -A3 "barcode=barcode04" | gzip >  barcode04.fastq.gz
gunzip -c tech_rep_2_ligation_reads.fastq.gz | grep -A3 "barcode=barcode05" | gzip >  barcode05.fastq.gz
gunzip -c tech_rep_2_ligation_reads.fastq.gz | grep -A3 "barcode=barcode07" | gzip >  barcode07.fastq.gz
```

  
#### Repairing fastqs

* There were some issues with the resulting long read FASTQ files as lengths of quality and sequence did not match for some reads
* There seemed to be some stray "--" characters. Accordingly I wrote a little script called validate_fastq.py to get all the validated fastq reads

```
# repeat  for all 4 read sets

validate_fastq.py -i barcode01.fastq.gz -o barcode01_fixed.fastq.gz
validate_fastq.py -i barcode02.fastq.gz -o barcode02_fixed.fastq.gz
validate_fastq.py -i barcode03.fastq.gz -o barcode03_fixed.fastq.gz
validate_fastq.py -i barcode04.fastq.gz -o barcode04_fixed.fastq.gz
validate_fastq.py -i barcode05.fastq.gz -o barcode05_fixed.fastq.gz
validate_fastq.py -i barcode07.fastq.gz -o barcode07_fixed.fastq.gz
```

* This recovered around 99% of the original reads (eyeballing by file size) indicating only a few reads were corrupt (sp not enough to impact downstream given how deeply these were sequenced).


#### Subsampling Fastqs

* I ran [rasusa](https://github.com/mbhall88/rasusa) to get 50x coverage on each fixed fastq file to reduce the read sizes and to make sure when I pooled the reads, some subsets with higher depth did not predominate.

```
conda create -n rasusa rasusa
conda install rasusa

# repeat for all 4 long read sets within the directory

rasusa --coverage 50 --genome-size 4mb --input barcode01_fixed.fastq.gz | gzip > barcode01_subsampled.fastq.gz
rasusa --coverage 50 --genome-size 5mb --input barcode02_fixed.fastq.gz | gzip > barcode02_subsampled.fastq.gz
rasusa --coverage 50 --genome-size 5mb --input barcode03_fixed.fastq.gz | gzip > barcode03_subsampled.fastq.gz
rasusa --coverage 50 --genome-size 2.1mb --input barcode04_fixed.fastq.gz | gzip > barcode04_subsampled.fastq.gz
rasusa --coverage 50 --genome-size 6mb --input barcode05_fixed.fastq.gz | gzip > barcode05_subsampled.fastq.gz
rasusa --coverage 50 --genome-size 6mb --input barcode07_fixed.fastq.gz | gzip > barcode07_subsampled.fastq.gz
conda deactivate
```


#### Creating pooled reads 

```
# pool the different library preps

mkdir -p tech_rep_1_pooled

cat tech_rep_1_rapid_reads/barcode01_subsampled.fastq.gz tech_rep_1_ligation_reads/barcode01_subsampled.fastq.gz > tech_rep_1_pooled/barcode01.fastq.gz
cat tech_rep_1_rapid_reads/barcode02_subsampled.fastq.gz tech_rep_1_ligation_reads/barcode02_subsampled.fastq.gz > tech_rep_1_pooled/barcode02.fastq.gz
cat tech_rep_1_rapid_reads/barcode03_subsampled.fastq.gz tech_rep_1_ligation_reads/barcode03_subsampled.fastq.gz > tech_rep_1_pooled/barcode03.fastq.gz
cat tech_rep_1_rapid_reads/barcode04_subsampled.fastq.gz tech_rep_1_ligation_reads/barcode04_subsampled.fastq.gz > tech_rep_1_pooled/barcode04.fastq.gz
cat tech_rep_1_rapid_reads/barcode05_subsampled.fastq.gz tech_rep_1_ligation_reads/barcode05_subsampled.fastq.gz > tech_rep_1_pooled/barcode05.fastq.gz
cat tech_rep_1_rapid_reads/barcode07_subsampled.fastq.gz tech_rep_1_ligation_reads/barcode07_subsampled.fastq.gz > tech_rep_1_pooled/barcode07.fastq.gz


mkdir -p tech_rep_2_pooled

cat tech_rep_2_rapid_reads/barcode01_subsampled.fastq.gz tech_rep_2_ligation_reads/barcode01_subsampled.fastq.gz > tech_rep_2_pooled/barcode01.fastq.gz
cat tech_rep_2_rapid_reads/barcode02_subsampled.fastq.gz tech_rep_2_ligation_reads/barcode02_subsampled.fastq.gz > tech_rep_2_pooled/barcode02.fastq.gz
cat tech_rep_2_rapid_reads/barcode03_subsampled.fastq.gz tech_rep_2_ligation_reads/barcode03_subsampled.fastq.gz > tech_rep_2_pooled/barcode03.fastq.gz
cat tech_rep_2_rapid_reads/barcode04_subsampled.fastq.gz tech_rep_2_ligation_reads/barcode04_subsampled.fastq.gz > tech_rep_2_pooled/barcode04.fastq.gz
cat tech_rep_2_rapid_reads/barcode05_subsampled.fastq.gz tech_rep_2_ligation_reads/barcode05_subsampled.fastq.gz > tech_rep_2_pooled/barcode05.fastq.gz
cat tech_rep_2_rapid_reads/barcode07_subsampled.fastq.gz tech_rep_2_ligation_reads/barcode07_subsampled.fastq.gz > tech_rep_2_pooled/barcode07.fastq.gz


# pool the tech reps
# then I combined the technical reps

mkdir -p pooled_long_reads

cat tech_rep_1_pooled/barcode01.fastq.gz tech_rep_2_pooled/barcode01.fastq.gz > pooled_long_reads/Acinetobacter_baumannii.fastq.gz 
cat tech_rep_1_pooled/barcode02.fastq.gz tech_rep_2_pooled/barcode02.fastq.gz > pooled_long_reads/Citrobacter_koseri.fastq.gz 
cat tech_rep_1_pooled/barcode03.fastq.gz tech_rep_2_pooled/barcode03.fastq.gz > pooled_long_reads/Enterobacter_kobei.fastq.gz 
cat tech_rep_1_pooled/barcode04.fastq.gz tech_rep_2_pooled/barcode04.fastq.gz > pooled_long_reads/Haemophilus_unknown.fastq.gz 
cat tech_rep_1_pooled/barcode05.fastq.gz tech_rep_2_pooled/barcode05.fastq.gz > pooled_long_reads/Klebsiella_oxytoca.fastq.gz 
cat tech_rep_1_pooled/barcode07.fastq.gz tech_rep_2_pooled/barcode07.fastq.gz > pooled_long_reads/Klebsiella_variicola.fastq.gz 
```

* Finally, I subsampled the reads to 60x coverage with rasusa

```
# from inside the pooled_long_reads directory

cd pooled_long_reads

rasusa --coverage 60 --genome-size 4mb --input Acinetobacter_baumannii.fastq.gz | gzip > Acinetobacter_baumannii_subsampled60x.fastq.gz 
rasusa --coverage 60 --genome-size 5mb --input Citrobacter_koseri.fastq.gz | gzip > Citrobacter_koseri_subsampled60x.fastq.gz
rasusa --coverage 60 --genome-size 5mb --input Enterobacter_kobei.fastq.gz | gzip > Enterobacter_kobei_subsampled60x.fastq.gz
rasusa --coverage 60 --genome-size 2.1mb --input Haemophilus_unknown.fastq.gz | gzip > Haemophilus_unknown_subsampled60x.fastq.gz
rasusa --coverage 60 --genome-size 6mb --input Klebsiella_oxytoca.fastq.gz | gzip > Klebsiella_oxytoca_subsampled60x.fastq.gz
rasusa --coverage 60 --genome-size 6mb --input Klebsiella_variicola.fastq.gz | gzip > Klebsiella_variicola_subsampled60x.fastq.gz
```



##### Creating Illumina pooled reads 

```
mkdir -p pooled_illumina

cat tech_rep_2_illumina_reads/Acinetobacter_baumannii_J9/J9_S142_L001_R1_001.fastq.gz tech_rep_1_illumina_reads/Acinetobacter_baumannii_J9/plasmids-J9_1.fastq.gz > pooled_illumina/Acinetobacter_baumannii_1.fastq.gz
cat tech_rep_2_illumina_reads/Acinetobacter_baumannii_J9/J9_S142_L001_R2_001.fastq.gz tech_rep_1_illumina_reads/Acinetobacter_baumannii_J9/plasmids-J9_2.fastq.gz > pooled_illumina/Acinetobacter_baumannii_2.fastq.gz

cat tech_rep_2_illumina_reads/Citrobacter_koseri_MINF_9D/MINF-9D_S143_L001_R1_001.fastq.gz tech_rep_1_illumina_reads/Citrobacter_koseri_MINF_9D/plasmids-MINF-9D_1.fastq.gz > pooled_illumina/Citrobacter_koseri_1.fastq.gz
cat tech_rep_2_illumina_reads/Citrobacter_koseri_MINF_9D/MINF-9D_S143_L001_R2_001.fastq.gz tech_rep_1_illumina_reads/Citrobacter_koseri_MINF_9D/plasmids-MINF-9D_2.fastq.gz > pooled_illumina/Citrobacter_koseri_2.fastq.gz

cat tech_rep_2_illumina_reads/Enterobacter_kobei_MSB1_1B/MSB1-1B_S144_L001_R1_001.fastq.gz tech_rep_1_illumina_reads/Enterobacter_kobei_MSB1_1B/plasmids-MSB1-1B_1.fastq.gz > pooled_illumina/Enterobacter_kobei_1.fastq.gz
cat tech_rep_2_illumina_reads/Enterobacter_kobei_MSB1_1B/MSB1-1B_S144_L001_R2_001.fastq.gz tech_rep_1_illumina_reads/Enterobacter_kobei_MSB1_1B/plasmids-MSB1-1B_2.fastq.gz > pooled_illumina/Enterobacter_kobei_2.fastq.gz

cat tech_rep_2_illumina_reads/Haemophilus_unknown_M1C132_1/M1C132-1_S145_L001_R1_001.fastq.gz tech_rep_1_illumina_reads/Haemophilus_unknown_M1C132_1/plasmids-MIC132-1_1.fastq.gz > pooled_illumina/Haemophilus_unknown_1.fastq.gz
cat tech_rep_2_illumina_reads/Haemophilus_unknown_M1C132_1/M1C132-1_S145_L001_R2_001.fastq.gz tech_rep_1_illumina_reads/Haemophilus_unknown_M1C132_1/plasmids-MIC132-1_2.fastq.gz > pooled_illumina/Haemophilus_unknown_2.fastq.gz

cat tech_rep_2_illumina_reads/Klebsiella_oxytoca_MSB1_2C/MSB1-2C_S146_L001_R1_001.fastq.gz  tech_rep_1_illumina_reads/Klebsiella_oxytoca_MSB1_2C/plasmids-MSB1-2C_1.fastq.gz  > pooled_illumina/Klebsiella_oxytoca_1.fastq.gz
cat tech_rep_2_illumina_reads/Klebsiella_oxytoca_MSB1_2C/MSB1-2C_S146_L001_R2_001.fastq.gz  tech_rep_1_illumina_reads/Klebsiella_oxytoca_MSB1_2C/plasmids-MSB1-2C_2.fastq.gz  > pooled_illumina/Klebsiella_oxytoca_2.fastq.gz

cat tech_rep_2_illumina_reads/Klebsiella_variicola_INF345/INF345_S147_L001_R1_001.fastq.gz  tech_rep_1_illumina_reads/Klebsiella_variicola_INF345/plasmids-INF345_1.fastq.gz   > pooled_illumina/Klebsiella_variicola_1.fastq.gz
cat tech_rep_2_illumina_reads/Klebsiella_variicola_INF345/INF345_S147_L001_R2_001.fastq.gz  tech_rep_1_illumina_reads/Klebsiella_variicola_INF345/plasmids-INF345_2.fastq.gz   > pooled_illumina/Klebsiella_variicola_2.fastq.gz
```

  
#### C222 


* The original files of these reads can also be found on the SRA SAMN32360844 in BioProject PRJNA914892
* These were my local reads hence not called SRR_...
* I subsampled the long reads to 60x to match the other pathogens, and also 30x to illustrate how Plassembler's speedup will be higher at lower read depths (with this read set, a chromosome will be assembled down to around 15x - this may be even lower for reads with higher N50s).

```
conda install rasusa
rasusa --coverage 30 --genome-size 2.5mb --input C222.fastq.gz | gzip > C222_subsampled30x.fastq.gz
rasusa --coverage 60 --genome-size 2.5mb --input C222.fastq.gz | gzip > C222_subsampled60x.fastq.gz
conda deactivate
```

* all reads were placed in the C222 directory


#### Running Plassembler 


*  Install plassembler - Linux. If you need Mac see the [instructions](https://github.com/gbouras13/plassembler#installation)

* First download the database

```
conda create -n plassembler plassembler unicycler==0.5.0
conda activate plassembler
install_database.py plassembler_db
```

* Then run Plassembler
* Plassembler was run using raw mode (-r) as these reads were a few years old and originally assembled with flye --nano-raw by Ryan. 
* They were assembled with 16 threads using -t 16 and with a minimum read depth of 500 (-m 500) and with the default minimum q-score of 9. 

```
mkdir -p Output_Final
mkdir -p Output_Final/Plassembler
mkdir -p Output_Final/Unicycler

P_DIR="Output_Final/Plassembler"
U_DIR="Output_Final/Unicycler"

# Wick et al
plassembler.py -d plassembler_db  -l pooled_long_reads/Acinetobacter_baumannii_subsampled60x.fastq.gz  -1 pooled_illumina/Acinetobacter_baumannii_1.fastq.gz  -2 pooled_illumina/Acinetobacter_baumannii_2.fastq.gz   -o $P_DIR/Acinetobacter_baumannii_pooled_plassembler -t 16 -p Acinetobacter_baumannii_pooled_plassembler -r -f -c 3500000 -m 500 
plassembler.py -d plassembler_db  -l pooled_long_reads/Citrobacter_koseri_subsampled60x.fastq.gz  -1 pooled_illumina/Citrobacter_koseri_1.fastq.gz  -2 pooled_illumina/Citrobacter_koseri_2.fastq.gz   -o $P_DIR/Citrobacter_koseri_pooled_plassembler -t 16 -p Citrobacter_koseri_pooled_plassembler -r -f -c 4500000 -m 500 
plassembler.py -d plassembler_db  -l pooled_long_reads/Enterobacter_kobei_subsampled60x.fastq.gz  -1 pooled_illumina/Enterobacter_kobei_1.fastq.gz  -2 pooled_illumina/Enterobacter_kobei_2.fastq.gz   -o $P_DIR/Enterobacter_kobei_pooled_plassembler -t 16 -p Enterobacter_kobei_pooled_plassembler -r -f -c 4500000 -m 500 
plassembler.py -d plassembler_db  -l pooled_long_reads/Haemophilus_unknown_subsampled60x.fastq.gz  -1 pooled_illumina/Haemophilus_unknown_1.fastq.gz  -2 pooled_illumina/Haemophilus_unknown_2.fastq.gz   -o $P_DIR/Haemophilus_unknown_pooled_plassembler -t 16 -p Haemophilus_unknown_pooled_plassembler -r -f -c 2000000 -m 500 
plassembler.py -d plassembler_db  -l pooled_long_reads/Klebsiella_oxytoca_subsampled60x.fastq.gz  -1 pooled_illumina/Klebsiella_oxytoca_1.fastq.gz  -2 pooled_illumina/Klebsiella_oxytoca_2.fastq.gz   -o $P_DIR/Klebsiella_oxytoca_pooled_plassembler -t 16 -p Klebsiella_oxytoca_pooled_plassembler -r -f -c 5500000 -m 500 
plassembler.py -d plassembler_db  -l pooled_long_reads/Klebsiella_variicola_subsampled60x.fastq.gz  -1 pooled_illumina/Klebsiella_variicola_1.fastq.gz  -2 pooled_illumina/Klebsiella_variicola_2.fastq.gz   -o $P_DIR/Klebsiella_variicola_pooled_plassembler -t 16 -p Klebsiella_variicola_pooled_plassembler -r -f -c 5000000 -m 500 
```

* C222
* These were not assembled with -r, as they were basecalled with Guppy 6.2.11 SUP mode and hence --nano-hq is more appropriate for Flye.

```
plassembler.py -d plassembler_db  -l C222/C222_subsampled30x.fastq.gz -1 C222/C222_S17_R1_001.fastq.gz -2 C222/C222_S17_R2_001.fastq.gz  -o $P_DIR/C222_30x -t 16 -p C222 -f -c 2500000 -m 500
plassembler.py -d plassembler_db  -l C222/C222_subsampled60x.fastq.gz -1 C222/C222_S17_R1_001.fastq.gz -2 C222/C222_S17_R2_001.fastq.gz  -o $P_DIR/C222_60x -t 16 -p C222 -f -c 2500000 -m 500
conda deactivate
```


#### Running Unicycler 

```
conda create -n unicycler unicycler 
conda activate unicycler 
U_DIR="Output_Final/Unicycler"
unicycler -l pooled_long_reads/Acinetobacter_baumannii_subsampled60x.fastq.gz  -1 pooled_illumina/Acinetobacter_baumannii_1.fastq.gz  -2 pooled_illumina/Acinetobacter_baumannii_2.fastq.gz   -o $U_DIR/Acinetobacter_baumannii_pooled_plassembler -t 16 
unicycler  -l pooled_long_reads/Citrobacter_koseri_subsampled60x.fastq.gz  -1 pooled_illumina/Citrobacter_koseri_1.fastq.gz  -2 pooled_illumina/Citrobacter_koseri_2.fastq.gz   -o $U_DIR/Citrobacter_koseri_pooled_plassembler -t 16
unicycler -l pooled_long_reads/Enterobacter_kobei_subsampled60x.fastq.gz  -1 pooled_illumina/Enterobacter_kobei_1.fastq.gz  -2 pooled_illumina/Enterobacter_kobei_2.fastq.gz   -o $U_DIR/Enterobacter_kobei_pooled_plassembler -t 16
unicycler -l pooled_long_reads/Haemophilus_unknown_subsampled60x.fastq.gz  -1 pooled_illumina/Haemophilus_unknown_1.fastq.gz  -2 pooled_illumina/Haemophilus_unknown_2.fastq.gz   -o $U_DIR/Haemophilus_unknown_pooled_plassembler -t 16 
unicycler -l pooled_long_reads/Klebsiella_oxytoca_subsampled60x.fastq.gz  -1 pooled_illumina/Klebsiella_oxytoca_1.fastq.gz  -2 pooled_illumina/Klebsiella_oxytoca_2.fastq.gz   -o $U_DIR/Klebsiella_oxytoca_pooled_plassembler -t 16 
unicycler -l pooled_long_reads/Klebsiella_variicola_subsampled60x.fastq.gz  -1 pooled_illumina/Klebsiella_variicola_1.fastq.gz  -2 pooled_illumina/Klebsiella_variicola_2.fastq.gz   -o $U_DIR/Klebsiella_variicola_pooled_plassembler -t 16 

# C222
unicycler -l C222/C222_subsampled30x.fastq.gz -1 C222/C222_S17_R1_001.fastq.gz -2 C222/C222_S17_R2_001.fastq.gz  -o $U_DIR/C222_30x -t 16
unicycler -l C222/C222_subsampled60x.fastq.gz -1 C222/C222_S17_R1_001.fastq.gz -2 C222/C222_S17_R2_001.fastq.gz  -o $U_DIR/C222_60x -t 16
conda deactivate
```


* All benchmarking output can be found at the Zenodo repository ____.


# Output 

Benchmarking Output is as follows:

Time & Accuracy
------

|                               | **Plassembler**    | **Unicycler**     | **Ground Truth**   |
|-------------------------------|--------------------|-------------------|--------------------|
| **_Acinetobacter baumannii_** |                    |                   |                    |
| Time (sec)                    | 1330               | 3938              |                    |
| Plasmids (bp)                 | 145059, 6078       | 145059, 6078      | 145059, 6078       |
| **_Citrobacter koseri_**      |                    |                   |                    |
| Time (sec)                    | 1321               | 4106              |                    |
| Plasmids (bp)                 | 64962, 9294        | 64962, 9294       | 64962, 9294        |
| **_Enterobacter kobei_**      |                         |                   |                    |
| Time (sec)                    | 2097               | 2097              |                    |
| Plasmids (bp)                 | 136482, 108411, 4665, 3715, 2370      | 136482, 108411, 4665, 3715, 2370  | 136482, 108411, 4665, 3715, 2370      |
| **_Haemophilus sp002998595_**      |                         |                   |                    |
| Time (sec)                    | 1325               | 3221              |                    |
| Plasmids (bp)                 | 39345, 10719, 9975     | 39345, 10719, 9975  | 39398, 10719, 9975, 7392, 5675     |
| **_Klebsiella oxytoca_**      |                         |                   |                    |
| Time (sec)                    | 1467               | 5552              |                    |
| Plasmids (bp)                 | 118161, 58472, 4574    | 118161, 58472, 4574 | 118161, 58472, 4574   |
| **_Klebsiella variicola_**      |                         |                   |                    |
| Time (sec)                    | 1816               | 4527              |                    |
| Plasmids (bp)                 | 250884, 243620, 31078 (linear), 5783, 3514  | 250902, 243534, 31078 (linear), 5783, 3514 | 250980, 243620, 31780 (linear), 5783, 3514  |
| **_Staphylococcus aureus_ 30x**     |                         |                   |                    |
| Time (sec)                    | 548               | 2600              |                    |
| Plasmids (bp)                 | 2473 | 2473 | 2473 |
| **_Staphylococcus aureus_ 60x**     |                         |                   |                    |
| Time (sec)                    | 897               | 3158              |                    |
| Plasmids (bp)                 | 2473 | 2473 | 2473 |



Small Plasmid Duplication
------


| **Small Plasmid Duplication**  | **Plassembler**   | **Flye (Output from Plassembler)**  |
|-------------------------------|--------------------|-------------------|
| **_Acinetobacter baumannii_** |                    |                   |                    
| Plasmids (bp)                 | 6078       | 12147    |
| **_Citrobacter koseri_**      |                    |                   |                    
| Plasmids (bp)                 | 9294        | 27773      |
| **_Enterobacter kobei_**      |                         |                   |                    
| Plasmids (bp)               | 4665, 3715, 2370                | 9652, (3715 plasmid missing), 4676              |                    
| **_Haemophilus sp002998595_**      |                         |                   |                    
| Plasmids (bp)                 | 10719, 9975     | 21402, 9962  | 
| **_Klebsiella oxytoca_**      |                         |                   |                               
| Plasmids (bp)                 | 4574    | 4566 | 
| **_Klebsiella variicola_**      |                         |                   |                                 
| Plasmids (bp)                 |  5783, 3514  |  11573, (3514 plasmid missing) |
| **_Staphylococcus aureus_ 30x**     |                         |                   |                    
| Plasmids (bp)                 | 2473 | 2471 | 
| **_Staphylococcus aureus_ 60x**     |                         |                   |                    
| Plasmids (bp)                 | 2473 | 1611 |
