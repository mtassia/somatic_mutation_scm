# Snakemake workflow

This directory contains the Snakemake pipeline used to run stochastic character mapping (SCM) on somatic SNVs.

## Overview

The workflow runs per-sample SCM in three stages:

1. Run SCM for each autosome (`chr1`-`chr22`) and write chromosome-level HDF5 files.
2. Merge chromosome-level HDF5 files into one sample-level HDF5 file.
3. Write a mutation-burden-scaled Newick tree for each sample.

## Directory layout

```text
smk/
|-- config/
|   `-- config_template.yaml        # Example multi-sample configuration
|-- example_data/                   # Example input files for two samples (test1 and test2)
|   |-- test1.vcf.gz
|   |-- test1.nwk
|   |-- test1.cellphy.bestModel
|   |-- test2.vcf.gz
|   |-- test2.nwk
|   `-- test2.cellphy.bestModel
|-- workflow/
|   |-- Snakefile                   # Main Snakemake workflow
|   |-- bin/                        # Scripts and helper functions called by rules
|   `-- envs/
|       `-- scm.yaml                # Rule-level conda environment for SCM scripts
`-- smk7_scm.yaml                   # Environment for Snakemake itself
```

Example data files have been included in the `example_data/` directory for testing purposes. Each VCF contains 10 loci per chromosome (chr1-22), randomly selected from a larger dataset.

## Input requirements

Each sample must provide:

- VCF of somatic SNVs (`.vcf.gz`)
- Newick-formatted phylogeny with branch lengths in substitutions per site (`.supportFBP` from `CellPhy`)
- Best-fit substitution model parameters (`.bestModel` from `CellPhy`)
- Outgroup label present in the tree

## Configuration

Edit `config/config_template.yaml` with your samples and runtime settings.

Expected structure:

```yaml
input:
  vcfs:
    sampleA: /path/to/sampleA.vcf.gz
  trees:
    sampleA: /path/to/sampleA.supportFBP
  outgroup:
    sampleA: outgroup_label
  models:
    sampleA: /path/to/sampleA.bestModel
run_params:
  n_threads: 16
  n_replicates: 1000
  mem_mb: 64000
```

Sample keys (e.g., `sampleA`) must be consistent across all input sections and can be named according to user preferences. 
Adjust `run_params` as needed for your computational environment.

## Outputs

Per sample, the workflow generates:

- `results/merged/{sample}.h5`
- `results/trees/{sample}.scm.nwk`

## Usage

Run from inside this directory (`smk/`).

### 1) Create Snakemake environment

```bash
mamba env create -f smk7_scm.yaml
conda activate smk7_scm
```

### 2) Dry-run (recommended)

```bash
snakemake --configfile config/config.yaml --profile profile/ -pn
```

### 3) Execute with profile (cluster)

```bash
snakemake --configfile config/config.yaml --profile profile/ -p
```

### 4) Execute locally

```bash
snakemake --configfile config/config.yaml --use-conda --cores 16 -p
```

## Rulegraph

![Rulegraph](rg.png)
