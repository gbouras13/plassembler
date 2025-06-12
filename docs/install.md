## Installation

Plassembler has been tested on Linux and MacOS machines. 

### Conda

The easiest way to install plassembler is via conda - Plassembler is on bioconda. 

```
conda install -c bioconda plassembler
```

or mamba for quicker solving:

```
mamba install -c bioconda plassembler
```

This will install all the dependencies along with plassembler.

### Pip

You can install the Python components of `plassembler` using pip.

```
pip install plassembler
```

You will then need to install the external dependencies separately, which can be found in `build/environment.yml`

* [Flye](https://github.com/fenderglass/Flye) >=2.9
* [Unicycler](https://github.com/rrwick/Unicycler) >=0.4.8
* [Minimap2](https://github.com/lh3/minimap2) >=2.11
* [fastp](https://github.com/OpenGene/fastp) >=0.24.2
* [chopper](https://github.com/wdecoster/chopper) >=0.5.0
* [mash](https://github.com/marbl/Mash) >=2.2
* [Raven](https://github.com/lbcb-sci/raven) >=1.8
* [Samtools](https://github.com/samtools/samtools) >=0.15.0

### Source

Alternatively, the development version of plassembler can be installed manually via github.

```
git clone https://github.com/gbouras13/plassembler.git
cd plassembler
pip install -e .
```

## Unicycler v0.5.0 Installation Issues

`plassembler` works best with Unicycler v0.5.0. With Unicycler v0.4.8, `plassembler` should still run without any issue and provide a satisfactory assembly, but you will be warned of this when you run `plassembler`. `plassembler` will not work with any older version of Unicycler.

**Linux**

For Linux environments, Unicycler v0.5.0 should be installed automaticall with the plassembler bioconda installation.

You can force it as follows:

`conda install -c bioconda plassembler unicycler==0.5.0`

or manually install Unicycler v0.5.0 after installing plassembler:

```
conda install -c bioconda plassembler
pip3 install git+https://github.com/rrwick/Unicycler.git
```

**MacOS**

For MacOS environments, the current conda installation method will only install the latest available bioconda Unicycler version of v0.4.8. 

Ryan Wick (the author of Unicycler) suggests that v0.5.0 should be used, as v0.4.8 is not compatible with the latest versions of spades (see [here](https://github.com/rrwick/Unicycler/releases/tag/v0.5.0) ). This will require another installation step on MacOS.

To install Unicycler v0.5.0, it is recommended that you install Unicycler from github after installing Plassembler follows:

```
# installs plassembler into an environment called 'plassemblerENV' and activates it
conda create -n plassemblerENV plassembler
conda activate plassemblerENV
# installs Unicycler v0.5.0
pip3 install git+https://github.com/rrwick/Unicycler.git
```

Mac M1 users may need to change some compiler settings and install from the Unicycler github repo e.g.

```
# installs plassembler into an environment called 'plassemblerENV' and activates it
conda create -n plassemblerENV plassembler
conda activate plassemblerENV
# installs Unicycler v0.5.0
git clone https://github.com/rrwick/Unicycler.git
cd Unicycler
python3 setup.py install --makeargs "CXX=g++"
```

# Beginner Conda Installation

If you are new to using the command-line, please install conda using the following instructions.

1. Install [Anaconda](https://www.anaconda.com/products/distribution). I would recommend [miniconda](https://docs.conda.io/en/latest/miniconda.html).
2. Assuming you are using a Linux x86_64 machine (for other architectures, please replace the URL with the appropriate one on the [miniconda](https://docs.conda.io/en/latest/miniconda.html) website).

`curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh`

For Mac (Intel, will also work with M1):

`curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh`

3. Install miniconda and follow the prompts.

`sh Miniconda3-latest-Linux-x86_64.sh`

4. After installation is complete, you should add the following channels to your conda configuration:

```
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

5. After this, conda should be installed (you may need to restart your terminal). It is recommended that mamba is also installed, as it will solve the enviroment quicker than conda:

`conda install mamba`

 6. Finally, I would recommend installing plassembler into a fresh environment. For example to create an environment called plassemblerENV with plassembler installed:

```
mamba create -n plassemblerENV plassembler
conda activate plassemblerENV
plassembler --help

```

