"""
conda mambabuild . --variants "{python: [3.8,3.9,3.10]}"
"""

from setuptools import setup, find_packages
import os

with open("README.md", "r") as fh:
    long_description = fh.read()

def get_version():
    with open(
        os.path.join(
            os.path.dirname(os.path.realpath(__file__)),
            "src",
            "VERSION",
        )
    ) as f:
        return f.readline().strip()
    

packages = find_packages()
# the directories where the code is
# so knows it should be included
package_data = {"src": ["src/*"]}

# for pypi
data_files = [(".", ["LICENCE", "README.md"])]

CLASSIFIERS = [
    "Environment :: Console",
    "Environment :: MacOS X",
    "Intended Audience :: Science/Research",
    "License :: OSI Approved :: MIT License",
    "Natural Language :: English",
    "Operating System :: POSIX :: Linux",
    "Operating System :: MacOS :: MacOS X",
    "Topic :: Scientific/Engineering :: Bio-Informatics",
]

# set version


setup(
    name="plassembler",
    version=get_version(),
    zip_safe=True,
    author="George Bouras",
    description="Plassembler: Automated Bacterial Plasmid Assembly Program",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author_email="george.bouras@adelaide.edu.au",
    packages=find_packages(),
    package_data=package_data,
    data_files=data_files,
    include_package_data=True,
    scripts=["plassembler"],
    url="https://github.com/gbouras13/plassembler",
    python_requires=">=3.8,<3.11",
    classifiers=CLASSIFIERS,
    install_requires=[
        "pyyaml>=6.0",
        "pytest-runner >= 5.0.0",
        "biopython >=1.76",
        "pytest>=6.2.5",
        "pandas",
        "loguru>=0.5.3",
        "Click>=8.0.0",
        "pytest-cov>=3.0.0",
        "pysam>=0.16.0",
    ],
)
