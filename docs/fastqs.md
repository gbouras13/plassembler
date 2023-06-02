##### Plassembler Benchmarking 

# FASTQ Set Up for Wick et al 

* Firstly, I downloaded the tarballs from this [link](https://bridges.monash.edu/articles/dataset/Small_plasmid_Nanopore_data/13543754). 

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
* There seemed to be some stray "--" characters. Accordingly I wrote a little script called validate_fastq.py to get all the validated fastq reads. You can find this in the Zenodo output directory here.

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

Other Read Sets
===========================

I didn't end up including these in the benchmarking (due to the lack of independent ground truth), but here is how you would go about getting these if you want them.

#### De Maio

* The fastqs were really easy to get - just use fastq-dl

```
fastq-dl --cpus 8  PRJNA422511  
```

* I got the assemblies from [here](https://figshare.com/articles/dataset/Hybrid_Enterobacteriaceae_assemblies_using_PacBio_Illumina_or_ONT_Illumina_sequencing/7649051), using the 'subsampled Nanopore' assemblies as the ground truth.

#### CAV1217 Assemblies

For this isolate (see the [paper](https://doi.org/10.1128/AAC.01823-16) for more information), I only have the assemblies available, not the reads. I downloaded them as follows:

```
conda activate ncbi-acc-download

# https://github.com/kblin/ncbi-acc-download

ncbi-acc-download --format fasta CP018676.1
ncbi-acc-download --format fasta CP018674.1
ncbi-acc-download --format fasta CP018672.1
ncbi-acc-download --format fasta CP018675.1
ncbi-acc-download --format fasta CP018673.1
```

Then I concatenated the files and manually edited the FASTA headers to add `circular=true`.




