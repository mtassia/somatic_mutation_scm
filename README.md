# Stochastic character mapping of somatic mutations

## Table of contents
- [Stochastic character mapping of somatic mutations](#stochastic-character-mapping-of-somatic-mutations)
  - [Table of contents](#table-of-contents)
  - [1. About](#1-about)
  - [2. Dependencies](#2-dependencies)
  - [3. Data prerequisites](#3-data-prerequisites)
  - [4. Installation](#4-installation)
  - [5. Tutorial](#5-tutorial)
    - [5.1 CLI interface:](#51-cli-interface)
    - [5.2 `R` interface:](#52-r-interface)
  - [6. To-do list](#6-to-do-list)
  - [7. Directory structure](#7-directory-structure)
 
## 1. About
This repository contains the tools for inferring the evolution of somatic mutations given an underlying phylogenetic hypothesis. 

*WIP*

## 2. Dependencies

- `phytools`

*WIP*

## 3. Data prerequisites
 - **Multisample VCF** of biallelic somatic variants with genotype likelihoods (`PL`; [see GATK specification for guidance](https://gatk.broadinstitute.org/hc/en-us/articles/360035890451-Calculation-of-PL-and-GQ-by-HaplotypeCaller-and-GenotypeGVCFs)) included in the `FORMAT` fields.
 - **Phylogenetic tree** of sample relationships.
   - Currently, this tree must include an outgroup for rooting the tree.
 - **Best-fit 10-state genotype substitution model ($Q$)** inferred by `cellphy`. 

*WIP*

## 4. Installation

*WIP*

## 5. Tutorial
### 5.1 CLI interface:

### 5.2 `R` interface:

*WIP*

## 6. To-do list
 1. [X] Add HPD intervals to SCM summary output.
 2. [X] Update software to handle multiple loci in a single run.
 3. [ ] Write a CLI wrapper for SCM.
 4. [X] Write function that will consolidate all output into an HDF5 file.
    - Summary output should have the following data:
       `c("locus","root","derived","n_state_changes","class","logL_Q","edge_profile")`
 5. [ ] Consolidate SCM SNP and indel functions into a single function that can handle both SNP and indel data (data ingress, data processing, and data summary).

## 7. Directory structure

```
|-- dev/
|   |-- bin/ # (TEMPORARY) R scripts for developing SCM tools
|   |-- smk/ # Snakemake workflow 
|       |-- config/
|           |-- config_template.yaml # Template config file for running SCM workflow
|       |-- workflow/ 
|           |-- bin/ # Scripts used in SCM workflow
|           |-- envs/ # Conda environments for SCM workflow
|           |-- src/ # Custom R functions used in SCM workflow
|       |-- scm.smk # Snakefile
|       |-- snakemake_env.yaml # Conda environment recipe for running Snakemake
|-- .gitignore
|-- LICENSE
|-- r_env.yaml # (TEMPORARY) environment for developing SCM tools
|-- README.md
```