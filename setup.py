from setuptools import setup, find_packages


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

v='0.1.1'


def main():
    setup(
        name="plassembler",
        packages=find_packages(),
        url="https://github.com/gbouras13/plassembler",
        python_requires=">=3.7",
        description="Automated Bacterial Plasmid Assembly Program",
        version=v,
        author="George Bouras",
        author_email="george.bouras@adelaide.edu.au",
        py_modules=["plassembler"],
        classifiers = CLASSIFIERS,
        install_requires=[
            "snakemake>=7.14.0",
            "pyyaml>=6.0",
            "Click>=8.1.3",
            'pytest-runner >= 5.0.0'
        ],
        entry_points={
            "console_scripts": ["plassembler = plassemblerModules.main:run", "plassembler.py = plassemblerModules.main:run"]
        },
        include_package_data=True,
    )

if __name__ == "__main__":
    main()