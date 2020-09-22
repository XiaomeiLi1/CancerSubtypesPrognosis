# CancerSubtypesPrognosis
CancerSubtypesPrognosis R package contains functions to predict breast cancer subtypes and prognosis.

## Introduction

Breast cancer is an extremely complex disease. Accurate prognosis and identification of subtypes of breast cancer at the molecular level are important steps towards effective and personalised treatments of the cancer. To this end, many computational methods have been developed to use gene expression data for breast cancer subtype discovery and prognosis. Meanwhile, microRNAs (miRNAs) have been extensively studied in the last two decades and their association with breast cancer subtypes and prognosis has been evidenced. However, it is not clear whether using miRNA data helps improve the performance of the gene expression based subtyping and prognosis methods, and this raises challenges as to how and when to use these methods in practice. In this paper, we conduct a comparative study of 34 methods, including 11 breast cancer subtyping methods and 23 breast cancer prognostic methods, on a collection of 13 independent breast cancer datasets. We aim to uncover the roles of miRNAs in breast cancer subtyping and prognosis and provide a set of recommendations for practical use of the computational methods. In addition, we create an R package, CancerSubtypesPrognosis, to include all the 34 methods to facilitate the reproducibility of the methods and streamline the evaluation. We expect this work can provide biomedical researchers with a practical guide to select the appropriate methods,  a software tool to apply the methods to their data, and facilitate the development of new computational methods for breast cancer subtyping and prognosis.

## Installation
CancerSubtypesPrognosis runs in the R statistical computing environment. R version 3.6.1 or higher and Bioconductor version 3.11 or higher are required.
1. Please install Bioconductor, you can use the following code in R

```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.11")
```

2. Please install Bioconductor dependencies required by CancerSubtypesPrognosis

```
BiocManager::install(c('limma', 'ConsensusClusterPlus', 'survcomp'))
```

3. Please install the devtools package and other dependencies required by CancerSubtypesPrognosis

```
install.packages(c("devtools","SNFtool","iCluster", "PINSPlus", "iC10", "doParallel", "foreach"), dependencies = TRUE)
```

4. Please install CancerSubtypesPrognosis and CIMILR packages from github repository

```
library("devtools")
install_github("danro9685/CIMLR", ref = 'R')

install_github('XiaomeiLi1/CancerSubtypesPrognosis')
```

5. Load CancerSubtypesPrognosis library

```
library("CancerSubtypesPrognosis")
```

## Documentation
Detailed information about the functions implemented in CancerSubtypesPrognosis can be found in the user [manual](https://github.com/XiaomeiLi1/CancerSubtypesPrognosis/blob/master/CancerSubtypesPrognosis_1.0.0.pdf)

Experiments implemented in our paper can be found in [inst/main_cancersubtypes.R](https://github.com/XiaomeiLi1/CancerSubtypesPrognosis/blob/master/inst/main_cancersubtypes.R) and [inst/main_cancerprognosis.R](https://github.com/XiaomeiLi1/CancerSubtypesPrognosis/blob/master/inst/main_cancerprognosis.R)

## Data
The METABRIC data need to be download from the [EMBL-EBI repository](https://www.ebi.ac.uk/ega/) (accession number EGAS00000000083 and EGAS00000000122, require individual access agreement)

Please install packages for datasets TRANSBIG, UNT, UPP, MAINZ, NKI, and load the packages into R
```
BiocManager::install(c('breastCancerTRANSBIG', 'breastCancerUNT', 'breastCancerUPP', 'breastCancerMAINZ', 'breastCancerNKI'))
library(breastCancerTRANSBIG)
library(breastCancerUNT)
library(breastCancerUPP)
library(breastCancerMAINZ)
library(breastCancerNKI)
```

Please download the GSE6532 dataset [here](https://github.com/XiaomeiLi1/CancerSubtypesPrognosis/releases/tag/V1.0)

Please find the other datasets employed in our paper in the folder [data](https://github.com/XiaomeiLi1/CancerSubtypesPrognosis/tree/master/data)
