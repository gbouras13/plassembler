The easiest way to install plassembler is via conda.

`conda install -c gbouras13 plassembler`

This will install all the dependencies along with plassembler.

Alternatively, plassembler can be installed manually via github

`git clone https://github.com/gbouras13/plassembler.git`

The dependencies found in environment.yml will then need to be installed manually.

For example using conda:

```
git clone https://github.com/gbouras13/plassembler.git
cd plassembler
conda env create -f environment.yml
conda activate plassembler_env
plassembler.py -h
```
