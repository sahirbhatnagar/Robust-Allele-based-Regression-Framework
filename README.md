# Robust-Allele-based-Regression
Robust Allele-based (RA) Regression Framework for Genetic Association Studies


## Table of Contents
- [Import RA source functions](#import_functions)
- [Independent samples with multiple Ys and Zs](#indep)
- [Mixture of Indepedent and Sibling Pairs](#indep_sib_mix)

## :heavy_plus_sign: Examples
- <a name="import_functions"></a> Import [RA functions](https://github.com/lzhangdc/Robust-Allele-based-Regression-Framework/blob/main/RA_functions.R)
```R
source('RA_functions.R')
```

- <a name="indep"></a> [Independent samples] Independent samples, Y can be continuous or binary
```R
########################################################################
########################################################################
#### Function inputs:
#### - g: genotype vector (G=0, 1, 2), n x 1 vector
#### - y: phenotype traits, n x J matrix 
#### - z: covariates, n x K matrix
#### - HWE: logical, TRUE or FALSE
########################################################################
#### Function output:
#### - p-value of the RA association test 
########################################################################
########################################################################

nIndep <- 1000; geno_prob <- c(0.64, 0.32, 0.04)
g_test <- sample(0:2, size=nIndep, replace=T, prob=geno_prob)
y_test <- cbind(replicate(3, rnorm(nIndep)), sample(0:1, size=nIndep, replace=T, prob=c(0.8, 0.2)))
z_test <- cbind(replicate(2, rnorm(nIndep, mean=3, sd=3)), replicate(2, rnorm(nIndep, mean=4, sd=4)))

## Independent samples, 3 continuous and 1 binary phenotype traits, no covariates, no assumption of HWE
RA_assoc_indep(g=g_test, y=y_test, z=NULL, HWE=F)

## Independent samples, 3 continuous and 1 binary phenotype traits, no covariates, assuming HWE
RA_assoc_indep(g=g_test, y=y_test, z=NULL, HWE=T)

## Independent samples, 3 continuous and 1 binary phenotype traits, 4 covariates, no assumption of HWE
RA_assoc_indep(g=g_test, y=y_test, z=z_test, HWE=F)

## Independent samples, 3 continuous and 1 binary phenotype traits, 4 covariates, assuming HWE
RA_assoc_indep(g=g_test, y=y_test, z=z_test, HWE=T)

```
- <a name="indep_sib_mix"></a> [Mixture of Indepedent and Sibling Pairs] 

```R
########################################################################
########################################################################
#### Function inputs:
#### - gIndep: genotype vector (G=0, 1, 2) of the independent samples
#### - yMat_indep: phenotype traits (nIndep x J matrix) of the independent samples
#### - gSib: genotype vector (G=0, 1, 2) of the sibling pairs
#### - yMat_sib: phenotype traits (nSib x J matrix) of the sibling pairs
########################################################################
#### Function output:
#### - p-value of the RA association test 
########################################################################
########################################################################

nIndep <- 1000; nSib <- 500
gInd <- sample(0:2, size=nIndep, replace=T, prob=c(0.64, 0.32, 0.04))
gSib <- sample(0:2, size=nSib, replace=T, prob=c(0.64, 0.32, 0.04))

y1Indep <- rnorm(nIndep); y2Indep <- rnorm(nIndep)
yMatIndep <- matrix(c(y1Indep, y2Indep), ncol=2)

y1Sib <- rnorm(nSib); y2Sib <- rnorm(nSib)
yMatSib <- matrix(c(y1Sib, y2Sib), ncol=2)

RA_pedi_assoc(gIndep=gInd, gSib=gSib, yMat_indep=yMatIndep, yMat_sib=yMatSib)

```
