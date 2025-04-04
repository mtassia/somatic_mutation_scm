
#!/usr/bin/env python3

"""
TITLE:
    Stochastic character mapping of somatic mutations 

AUTHORS:
    M. Tassia

DESCRIPTION:
    Application of stochastic character mapping to infer the evolutionary history of somatic mutations given:
    1. a phylogenetic tree (.nwk)
    2. a set of somatic mutations (.vcf)
    3. genotype substitution model (.model)
    This method leverages both the phylgoenetic tree and best-fit 10-state substitution model parameters output
    from CellPhy. 

DIRECTORY STRUCTURE:
    .
    ├── logs
    ├── profile
    │   └── config.yaml
    ├── scm.smk
    └── workflow
        |── bin
        |   └── ...
        └── envs
            └── ...

USAGE: 
    snakemake --configfile config/config.yaml --profile /path/to/profile_dir/ -p # slurm
    snakemake --configfile config/config.yaml --profile /path/to/profile_dir/ --use-conda --cores 16 -p # local
"""

#--- IMPORTS ---#
import pathlib

#--- CONFIGURATION ---#
# configure shell behavior for all rules
shell.executable("/bin/bash")
shell.prefix("set -euo pipefail;")

# # Check if config file is provided
# if not config:
#     raise ValueError("ERROR: No config file provided. Please provide a config file using --configfile.")

# Create a log directory
pathlib.Path("logs/").mkdir(parents=True, exist_ok=True)

# Include rules:
# include: "<path/to/snakefile>"

#--- VARIABLES ---#
chrom = ["chr" + str(c) for c in list(range(1, 23))]

#--- RULES ---#
localrules:
    all,

rule all:
    input:
        config['output_dir'] + "/stats/sample_read_counts.jpg",

