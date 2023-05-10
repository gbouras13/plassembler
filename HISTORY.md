# History

1.0.0 (2023-05-10)
------------------

* Large overhaul.
* Adds `--pacbio_model` for pacbio data
* Replaces nanofilt with chopper
* Adds `-a` for assembled mode
* Adds subsampling with rasusa, and `-s` to change subsampling depth, and `--no_subsample` to turn off subsampling
* Adds `--keep_fastqs`
* Adds `--keep_chromosome`
* Refactors mapping code 
* Adds custom function to identify multimapped reads
* Changes output formats - consolidates all output into `summary.tsv`.


0.1.5 (2023-01-16)
------------------

* Adds checks for dependencies.
* Adds samtools to bioconda recipe (thanks Jan/gaworj).
* Adds long-only kmer_mode for very high quality Nanopore reads (R10.4 and above) - experimental until I get more data especially with small plasmids. Does exactly the same (Flye -> Unicycler). Seems to work pretty well. 


0.1.4 (2023-01-03)
------------------

* Adds install_database.py, the Plassembler database and functionality for mapping plasmid contigs to PLDSB.
* Update the API to -1 and -2 for short reads, matching Unicycler.
* Adds mash as dependency.
* Adds plassembler_top_hits_mash_plsdb.tsv output.
* Adds plasmid_copy_number_short and plasmid_copy_number_long to fasta header for each plasmid.
* Checks dependencies upon initialisation, with message to install Unicycler manually if required.

0.1.3 (2022-12-27)
------------------

* Fix bugs in bioconda release - unicycler.py and flye.py conflicts

0.1.2 (2022-12-22)
------------------

* Code refactored.
* Tests added.
* Bioconda release.

0.1.1 (2022-08-09)
------------------

* Adds module find plasmids where Flye assembles only 1 complete chromosome.


0.1.0 (2022-08-08)
------------------

* First release
