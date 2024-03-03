# History

1.6.1 (2024-03-02)
------------------

* Bug fixes and added tests for `--depth_filter` which would crash in some scenarios.

1.6.0 (2024-03-01)
------------------

* Adds `--depth_filter`. This will filter out all contigs that have long- (and short-read for `plassembler run`) read copy numbers that are less than the specified depth filter. Defaults to 0.25x.
* Adds `--unicycler_options`  and `--spades_options` which allows passing extra Unicycler options (#46)

1.5.1 (2024-02-01)
------------------

* Fixes #44, where `--use_raven` was not working (Flye would be used instead). Thanks @[wanyuac](https://github.com/wanyuac).

1.5.0 (2023-11-21)
------------------

* **If you upgrade to v1.5.0, you will need to update the database using `plassembler download`** 
* Plassembler v1.5.0 incorporates a new database thanks to the recent PLSDB release [2023_11_03_v2](https://ccb-microbe.cs.uni-saarland.de/plsdb/). Thanks @[biobrad](https://github.com/biobrad) for the heads up.

1.4.1 (2023-10-30)
------------------

* Fixes bug with `plassembler run`, which would exit ungracefully if the isolate had more than 1 chromosome, but no plasmids were recovered by Unicycler (e.g. ATCC [17802](https://www.atcc.org/products/17802)).

1.4.0 (2023-10-27)
------------------

* Adds `--no_chromosome` option to `plassembler long` and `plassembler run` after a request to allow for the assembly of read sets that have only plasmids.
* Using this will skip Flye and create a dummy 3MB chromosome full of A's.
* Fixes another bug here [issue](https://github.com/gbouras13/plassembler/issues/37).


1.3.0 (2023-10-24)
------------------

* `plassembler long` should yield improved results. It achieves this by treating long reads as both short reads (in the sense of creating a de Brujin graph based assembly) and long reads (for scaffolding) in Unicycler.
* While I'd still recommend short reads if you can get them, I am now confident that if your isolate has small plasmids in the long read set, `plassembler long` should find them.
* For more information, see the [documentation](https://plassembler.readthedocs.io/en/latest/long/).
* The ability to specify a `--flye_assembly` and `--flye_info` if you already have a Flye assembly for your long reads instead of `--flye_directory` has been added. Thanks to @[incoherentian](https://github.com/incoherentian)'s [issue](https://github.com/gbouras13/plassembler/issues/37)
* The ability to specify a `--no_copy_numbers` with `plassembler assembled` if you just want to run some plasmids against the PLSDB has been added. Thanks to @[gaworj](https://github.com/gaworj)'s [issue](https://github.com/gbouras13/plassembler/issues/36).


1.2.0 (2023-09-12)
------------------

`plassembler` v1.2.0 implements the following features:

* `plassembler long` officially released and implemented using [Canu](https://github.com/marbl/canu) and [dnaapler](https://github.com/gbouras13/dnaapler) to reassemble unmapped reads in place of Unicycler for `plassembler run`. While we'd still recommend getting short reads if you really want to recover plasmids, as long as your long reads are short enough (i.e. not size selected), `plassembler long` should hopefully recover most small plasmids.
* For more information, see the [documentation](https://plassembler.readthedocs.io/en/latest/long/).
* Faster mapping thanks to @[fanvanf](https://github.com/fanvanf)'s [issue](https://github.com/gbouras13/plassembler/issues/29).
* The ability to specify a `--flye_directory` if you already have a Flye assembly for your long reads, which will tell `plassembler` to skip the long read assembly step.

1.1.0 (2023-06-02)
------------------

* Refactored codebase and release on pypi
* Adds unit tests and CI
* Replace `argparse` with `click`
* Adds option to skip chopper and fastq to `plassembler run` using `--skip_qc`
* Breaking CLI changes to be compatible with click
* `plassembler.py` changed to `plassembler run`
* Adds Raven long read assembly option to `plassembler run` using `--use_raven`
* `install_database.py` changed to `plassembler download`
* Assembled mode now `plassembler assembled`
* Untested/experimental long read only mode using `plassembler long`
* Removes rasusa, `-s` and `--no_subsample`. If users want faster runtimes, we recommend `--use_raven` or conduct subsampling prior.


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
