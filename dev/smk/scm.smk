
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
    snakemake --configfile config.yaml --profile /path/to/profile_dir/ --snakefile scm.smk -p # slurm
    snakemake --configfile config.yaml --snakefile scm.smk --use-conda --cores 16 -p # local
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
        expand("output/{sample}.{chrom}.h5",
                sample=config["input"]["vcfs"].keys(),
                chrom=chrom)

rule scm_chromosome:
    input:
        tree = lambda wildcards: config["input"]["trees"][wildcards.sample],
        vcf = lambda wildcards: config["input"]["vcfs"][wildcards.sample],
        model = lambda wildcards: config["input"]["models"][wildcards.sample],
        script = "workflow/bin/chromosome_scm.R",
        functions = "workflow/src/SCM.smk.R"
    output:
        h5 = "output/{sample}.{chrom}.h5"
    params:
        outgroup = lambda wildcards: config["input"]["outgroup"][wildcards.sample],
        prefix = "output/{sample}.{chrom}",
        reps = config["run_params"]["n_replicates"]
    resources:
        time = "24:00:00",
        mem_mb = config["run_params"]["mem_mb"]
    threads:
        config["run_params"]["n_threads"]
    conda:
        "workflow/envs/scm.yaml"
    shell:
        """
        Rscript {input.script} \
            -f {input.functions} \
            -v {input.vcf} \
            -t {input.tree} \
            -m {input.model} \
            -p {params.prefix} \
            -o {params.outgroup} \
            -n {threads} \
            -c {wildcards.chrom}
        """

#rule merge_h5:
#    input:
#        expand("output/{sample}.{chrom}.h5",
#                sample=config["input"]["vcfs"].keys(),
#                chrom=chrom)
#    output:
#        "output/{sample}.h5"
#    params:
#        prefix = "output/{sample}"
#    resources:
#        time = "24:00:00",
#        mem_mb = config["run_params"]["mem_mb"]
#    threads:
#        config["run_params"]["n_threads"]
#    conda:
#        "workflow/envs/scm.yaml"
#    shell:
#        """
#        """
