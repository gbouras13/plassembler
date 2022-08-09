To run plassembler

` plassembler.py -l <long read fastq> -o <output dir> -s1 < short read R1 fastq> -s2 < short read R2 fastq>  -c <estimated chromosome length>`

* -c will default to 2500000 (a lower bound for the estimated genome length of staphylococcus aureus).

To specify threads:

` plassembler.py -l <long read fastq> -o <output dir> -s1 < short read R1 fastq> -s2 < short read R2 fastq>  -c <estimated chromosome length> -t <threads>`

To specify a prefix for the output files:

` plassembler.py -l <long read fastq> -o <output dir> -s1 < short read R1 fastq> -s2 < short read R2 fastq>  -c <estimated chromosome length> -t <threads> -p <prefix>`

To specify a minimum length and minimum read quality Q-score for nanofilt :

` plassembler.py -l <long read fastq> -o <output dir> -s1 < short read R1 fastq> -s2 < short read R2 fastq>  -c <estimated chromosome length> -t <threads> -p <prefix> -m <min length> -q <min quality>`

* -m will default to 1000 and -q will default to 9. Note that for some tiny plasmids, -m may need to be reduced (see [https://www.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000631#tab2](https://www.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000631#tab2)).

To overwrite an existing output directory, use -f

` plassembler.py -l <long read fastq> -o <output dir> -s1 < short read R1 fastq> -s2 < short read R2 fastq>  -c <estimated chromosome length> -t <threads> -p <prefix> -m <min length> -q <min quality> -f`

plassembler defaults to 8 threads.
