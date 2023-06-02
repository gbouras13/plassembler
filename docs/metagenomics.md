# Metagenomics

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


Contigs 34, 61, 87, 101 and 109 match what was found in the ground truth.

| contig | length | mean_depth_short | circularity  | PLSDB_hit | ACC_NUCCORE   | Description_NUCCORE                                                                                                 | plasmid_copy_number_short | plasmid_copy_number_long |
| ------ | ------ | ---------------- | ------------ | --------- | ------------- | ------------------------------------------------------------------------------------------------------------------- | ------------------------- | ------------------------ |
| 34     | 110007 | 336.5            | circular     | Yes       | NZ_CP061531.1 | Escherichia coli strain WEM25 plasmid p1, complete sequence                                                         | 1.67                      | 1.37                     |
| 61     | 49661  | 357.16           | not_circular | Yes       | NZ_CP012345.2 | Salmonella enterica subsp. enterica serovar Choleraesuis str. ATCC 10708 plasmid pCFSAN000679_01, complete sequence | 1.78                      | 3.17                     |
| 83     | 9628   | 39.11            | not_circular | Yes       | NZ_CP069918.1 | Klebsiella oxytoca strain FDAARGOS_1334 plasmid unnamed7                                                            | 0.19                      | 0                        |
| 87     | 6367   | 9554.23          | circular     | Yes       | NZ_CP013628.1 | Staphylococcus aureus strain RIVM4293 plasmid pRIVM4293, complete sequence.                                         | 47.54                     | 12.67                    |
| 91     | 5355   | 14.61            | not_circular | Yes       | NZ_CP068597.1 | Paenibacillus sonchi strain LMG 24727 plasmid unnamed2, complete sequence                                           | 0.07                      | 0                        |
| 93     | 5010   | 13.32            | not_circular | Yes       | NZ_CP068597.1 | Paenibacillus sonchi strain LMG 24727 plasmid unnamed2, complete sequence                                           | 0.07                      | 0                        |
| 101    | 2993   | 2561.92          | circular     | Yes       | NZ_MH785226.1 | Staphylococcus aureus strain ph1 plasmid pRIVM1295-2, complete sequence                                             | 12.75                     | 1.65                     |
| 103    | 2789   | 1018.1           | not_circular | Yes       | CP048737.1    | Enterobacter sp. T2 plasmid unnamed1, complete sequence                                                             | 5.07                      | 4.75                     |
| 106    | 2667   | 1045.66          | not_circular | Yes       | NZ_CP066061.1 | Actinomyces oris strain FDAARGOS_1051 plasmid unnamed                                                               | 5.2                       | 4.79                     |
| 108    | 2337   | 16.17            | not_circular | Yes       | NZ_CP069918.1 | Klebsiella oxytoca strain FDAARGOS_1334 plasmid unnamed7                                                            | 0.08                      | 0                        |
| 109    | 2216   | 2188.82          | circular     | Yes       | NZ_CP013624.1 | Staphylococcus aureus strain RIVM1076 plasmid pRIVM1076, complete sequence.                                         | 10.89                     | 1.08                     |
| 141    | 1049   | 1015.41          | not_circular | Yes       | NZ_CP066061.1 | Actinomyces oris strain FDAARGOS_1051 plasmid unnamed                                                               | 5.05                      | 4.77                     |
| 145    | 830    | 930.78           | not_circular | Yes       | NZ_CP066061.1 | Actinomyces oris strain FDAARGOS_1051 plasmid unnamed                                                               |