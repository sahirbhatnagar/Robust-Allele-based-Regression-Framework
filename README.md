# Robust Allele-based Regression Framework
Robust Allele-based (RA) Regression Framework for Genetic Association Studies


## Table of Contents
- [Import RA source functions](#import_functions)
- [Independent samples with multiple Ys and Zs](#indep)
- [Mixture of indepedent and sibling pairs](#indep_sib_mix)
- [Test of HWE](#indep_HWE)
- [Samples from multiple population association testing](#multi_pop)

## :heavy_plus_sign: Examples
- <a name="import_functions"></a> Save [RA functions](https://github.com/lzhangdc/Robust-Allele-based-Regression-Framework/blob/main/RA_functions.R) in the working directory
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
########################################################################################
########################################################################################
#### Function inputs:
#### - gIndep: genotype vector (G=0, 1, 2) of the independent samples
#### - yMat_indep: phenotype traits (nIndep x J matrix) of the independent samples
#### - gSib: genotype vector (G=0, 1, 2) of the sibling pairs
#### - yMat_sib: phenotype traits (nSib x J matrix) of the sibling pairs
########################################################################################
#### Function output:
#### - p-value of the RA association test 
########################################################################################
########################################################################################

nIndep <- 1000; nSib <- 500

gInd <- sample(0:2, size=nIndep, replace=T, prob=c(0.64, 0.32, 0.04))
y1Indep <- rnorm(nIndep); y2Indep <- rnorm(nIndep)
yMatIndep <- matrix(c(y1Indep, y2Indep), ncol=2)

gSib_test <- simSib(nRep=1, f=nSib/2, p=0.2, p2=0)
gSib1_causal <- simSib(nRep=1, f=nSib/2, p=0.2, p2=0)
gSib2_causal <- simSib(nRep=1, f=nSib/2, p=0.3, p2=0)
y1Sib <- 0.03 + 0.02*gSib1_causal + rnorm(nSib)
y2Sib <- 0.01 + 0.03*gSib2_causal + rnorm(nSib)
yMatSib <- matrix(c(y1Sib, y2Sib), ncol=2)

RA_pedi_assoc(gIndep=gInd, gSib=gSib_test, yMat_indep=yMatIndep, yMat_sib=yMatSib)
```
- <a name="indep_HWE"></a> [Test of HWE] with independent samples

```R
########################################################################
########################################################################
#### Function inputs:
#### - g: genotypes of independent samples
########################################################################
#### Function output:
#### - p-value of the test of HWE
########################################################################
########################################################################

g1 <- sample(0:2, size=1000, replace=T, prob=c(0.64, 0.32, 0.04))
RA_HWE(g=g1)

g2 <- sample(0:2, size=1000, replace=T, prob=c(0.6, 0.3, 0.1))
RA_HWE(g=g2)

```
- <a name="multi_pop"></a> [Association testing] with independent samples from multiple population
```R
########################################################################
########################################################################
#### Function inputs:
#### - g: genotypes of independent samples
#### - y: phenotype traits (nIndep x J matrix) of the independent samples
#### - z: population factor
########################################################################
#### Function output:
#### - p-value of the RA association test
########################################################################
########################################################################

nI <- 500; nII <- 1500

yI <- sample(0:1, size=nI, replace=T, prob=c(0.8, 0.2))
yII <- sample(0:1, size=nII, replace=T, prob=c(0.5, 0.5))
y1 <- matrix(c(yI, yII), ncol=1)
y2 <- matrix(rnorm((nI+nII)), ncol=1)

g1 <- sample(0:2, size=nI, replace=T, prob=c(0.64, 0.32, 0.04))
g2 <- sample(0:2, size=nII, replace=T, prob=c(0.6, 0.3, 0.1))

z <- rep(c('I', 'II'), c(nI, nII))
g <- c(g1, g2)
y <- cbind(y1, y2)

RA_multi_pop(g, y, z)

```



