##### Metagenomics

As it currently stands, we do not recommend `plassembler` for metagenomic sequences. This is because of their high diversity, leading to difficulties in recovering chromosome-length contigs for bacteria. Additionally, Unicycler (a core dependency of Plassembler) is not recommended for metagenomes.

However, we anticipate that as sequencing becomes more accurate and cheaper, it will be increasingly possible to assemble plasmids using a `plassembler` like approach from metagenomes - it's a work in progress.

So as a test, we tried assembling the ZYMO HMW DNA Standard dataset from this [paper](https://www.nature.com/articles/s41592-022-01539-7), under ENA accession PRJEB48692. This mock community contains 7 bacteria and 1 fungus isolate. Notably, this dataset had extremely had deep (all bacterial chromosomes >100x coverage) and long (N50 > 20kbp) reads, so is unlikely to reflect your real-world metagenomic data as of 2023.

## Get Data

```
# installation
mamba create -n fastq-dl fastq-dl
conda activate fastq-dl

# downloads all the read sets
fastq-dl PRJEB48692	

conda deactivate
```

## Run Plassembler

We decided to use `-m 10000`, because we figures that smalll plasmids would be missed by Flye anyway, and wanted complete chromosome assemblies, and a `-c 500000`. We used 32 threads on 16 cores and allocated 80 GB of RAM.

```
plassembler run -d Plassembler_DB -l ERR7287988.fastq.gz -1 ERR7255689_1.fastq.gz -2 ERR7255689_2.fastq.gz \
-f -t 32 -q 10 -o zymo_R10.4_flye -m 10000 -c 500000
```

`plassembler` took around 8 hours (wall clock) to finish and excitingly we assembled all 7 bacterial chromosomes using Flye (unsurprising!) along with the 5 plasmids indicated in the ground truth (1 _E. coli_ 100kbp, 1 _S. enterica_ 49kbp and 3 small _S. aureus_ plasmids (6, 2 and 2 kbp)) with genome fraction 100% from QUAST. 

So in theory `plassembler`  might work on metagenomes, but I would caution against using it, for now.
