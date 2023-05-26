"""
conda mambabuild . --variants "{python: [3.6,3.7,3.8,3.9,3.10,3.11]}"
"""

from setuptools import setup, find_packages


with open("README.md", "r") as fh:
    long_description = fh.read()


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
    "License :: OSI Approved :: MIT license",
    "Natural Language :: English",
    "Operating System :: POSIX :: Linux",
    "Operating System :: MacOS :: MacOS X",
    "Programming Language :: Python :: 3.9",
    "Topic :: Scientific/Engineering :: Bio-Informatics",
]

# set version

v = "1.1.0"

setup(
    name="plassembler",
    version=v,
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
    entry_points={"console_scripts": ["install_database = src.install_database:main"] },
    url="https://github.com/gbouras13/plassembler",
    python_requires=">=3.6",
    classifiers=CLASSIFIERS,
    install_requires=[
        "pyyaml>=6.0",
        "pytest-runner >= 5.0.0",
        "biopython >=1.76",
        "pytest>=6.2.5",
        "pandas",
    ],
)
