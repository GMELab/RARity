
![logo](https://github.com/TheNaziaPathan/rarity_demo1/assets/113630451/868de257-61a0-4260-9bc1-173097a331fb)



Rare variant heritability (RARity) estimator is a framework to assess RV heritability (h<sup>2</sup><sub>RV</sub>) without assuming a particular genetic architecture, in a fast, accurate, and unbiased manner. It enables computation of both gene-level and exome-wide heritability estimates of continuous traits.

# Table of Contents
- [Method Overview](#1)
- [Demo requirements](#2)
    - [Hardware requirements](#2a)
    - [Software requirements](#2b)
- [Demo (with instructions for use)](#3)
    - [Input files](#3a)
    - [Steps](#3b)
- [Sample output files](#4)
- [License](#5)
- [Citation](#6)
- [Contact Information](#7)
  
# Method Overview <a name="1"></a>

The RARity method entails parallel computing of the adjusted R<sup>2</sup> based on an ordinary least square (OLS) multiple linear regression as an unbiased estimator of block-wise heritability for each consecutive genetic block. Adjusted R<sup>2</sup> estimates are then summed over all blocks as the overall heritability estimate of a trait.

To reduce the influence of long-range LD between variants that would otherwise inflate the overall heritability estimate, highly correlated RVs should be removed using PLINK 1.9, by LD pruning with a Pearson’s r<sup>2</sup> threshold \> 0.1 and, within a window of 50Mb that are shifted by 500 bases at the end of each step.

The following figure provides a summary of the RARity pipeline:

![rv flow v2](https://github.com/TheNaziaPathan/rarity_demo1/assets/113630451/6052132e-6ce9-4c16-a50b-91733a2697af)

# Demo Requirements  <a name="2"></a>

## Hardware Requirements (Full-scale version)  <a name="2a"></a>

Any standard computer (macOS, Linux, Windows).

Requires a unix-like virtual environment supporting a minimum of 250GB RAM space for in-memory operations.

## Software Requirements  <a name="2b"></a>

#### Essential Dependencies: programs
| Program | Description | Download |
| --- | --- | --- |
| BASH (≥ 5.0) | a unix shell and command language | [https://ubuntu.com/download/desktop](https://www.gnu.org/software/bash/) |
| R (≥ 3.6.3) or newer | R programming language | https://cran.r-project.org/ |

#### Essential Dependencies: R packages
| R package | Install | Reference |
| --- | --- | --- |
| dplyr | install.packages("dplyr") | https://www.r-project.org/nosvn/pandoc/dplyr.html |
| data.table | install.packages("data.table") | https://cran.r-project.org/package=data.table |

# Demo (with instructions for use)  <a name="3"></a>

The following is a demonstration of the RARity algorithm to estimate rare coding variant heritability.

The algorithm uses the **rarity_demo.r** script from the RARity repository. The R script that outputs the RV heritability for each genotype block, as well as the overall heritability estimates, based on all genotype blocks.

## Input files  <a name="3a"></a>

For demonstrative and efficiency purposes, we will use the following simulated data as input.

1.  5 Genotype matrices: geno_block1.RData, geno_block2.RData, geno_block3.RData, geno_block4.RData and geno_block5.RData, each, with a dimension of 10,000 individuals x 1,000 RVs. Assume that the genotype matrices are pre-processed, as per figure 1 and standardized to mean=0, sd=1.

2.  1 Phenotype matrix: pheno.Rdata, with a dimension of 10,000 individuals x 3 phenotypes. Assume that the phenotypes are pre-processed as per figure1 above and standardized to mean=0, sd =1.

## Steps  <a name="3b"></a>
**Step 1:** Create the following working directory in bash:

```
mkdir rarity_practice
```

**Step 2:** Download the “rarity_demo.r” as well as the input data (5 genotype block matrices and 1 phenotype matrix) to rarity_practice directory.

**Step 3:** In bash, tell the system that it has permission to execute the scripts:

```
chmod +x /your_directory/rarity_practice/\*
```

**Step 4:** Rscript requires the following arguments:

File pattern to identify the genotype block data= “geno_block”

pheno_file=“phenos.RData”

Run the script in bash as such:

```
cd your_directory/rarity_practice

Rscript ./rarity_demo.r geno_block phenos.RData
```


**Run time:** The run time for this script on a standard computer should be around 10 seconds.

# Sample Output  <a name="4"></a>
**A successful output will produce:**

1) A single file, BLOCK_HERITABILITY.txt computing RV heritability for each block.


Here is an example of the fields produced and their meaning:

| phenotype | block       | N     | N_RV | r2       | adj_r2   | block_r2_variance | block_adj_r2_variance |
|-----------|-------------|-------|------|----------|----------|-------------------|-----------------------|
| pheno_1   | geno_block1 | 10000 | 1000 | 0.106309 | 0.007    | 2.75E-05          | 3.39E-05              |
| pheno_2   | geno_block1 | 10000 | 1000 | 0.104668 | 0.005175 | 2.72E-05          | 3.35E-05              |
| pheno_3   | geno_block1 | 10000 | 1000 | 0.099288 | -0.0008  | 2.61E-05          | 3.22E-05              |

Acronym descriptions in output: N=number of individuals; N_RV= number of rare variants; adj_r2= adjusted R<sup>2</sup>; block_r2_variance= variance of R<sup>2</sup> in each block; Block_adj_r2_variance= Variance of adjusted R<sup>2</sup> in each block

2) A single file, TOTAL_HERITABILITY.txt showing the total heritability for each trait. The expected output of this file

| phenotype | N     | N_RV | Heritability | LCL      | UCL      |
|-----------|-------|------|--------------|----------|----------|
| pheno_1   | 10000 | 5000 | 0.005413     | -0.01962 | 0.030442 |
| pheno_2   | 10000 | 5000 | 0.009938     | -0.01517 | 0.035046 |
| pheno_3   | 10000 | 5000 | 0.008781     | -0.01631 | 0.033873 |

Acronym descriptions in output: N=number of individuals; N_RV= number of rare variants; LCL= Lower confidence level; UCL=Upper confidence level

# License  <a name="5"></a>

GNU General Public License v3.0

# Citation <a name="6"></a>

### *Nature Communications*

Under review.

# Contact Information  <a name="7"></a>

Any queries pertaining to RARity scripts or methodological framework can be addressed to either: Nazia Pathan ([pathann@mcmaster.ca](mailto:pathann@mcmaster.ca)) or Guillaume Pare (pareg@mcmaster.ca)

