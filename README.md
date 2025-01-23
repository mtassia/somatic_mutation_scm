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