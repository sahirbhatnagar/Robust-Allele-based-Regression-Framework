# Robust-Allele-based-Regression
Robust Allele-based (RA) Regression Framework for Genetic Association Studies

## :heavy_plus_sign: Quick Start
- [Import RA functions](https://github.com/lzhangdc/Robust-Allele-based-Regression-Framework/blob/main/RA_functions.R)
```R
source('RA_functions.R')
```

- [Independent samples] Independent samples, Y can be continuous or binary
```R
########################################################################
########################################################################
#### Function inputs:
#### - G: genotype vector (G=0, 1, 2), n x 1 vector
#### - Y: phenotype traits, n x J matrix 
#### - Z: covariates, n x K matrix
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
- [Mixture of Indepedent and Sibling Pairs] 

## :heavy_plus_sign: Examples
* [Case-control study with indepedent samples](https://github.com/lzhangdc/Robust-Allele-based-Regression-Framework/blob/main/vignette/case_control_study.md)
* [Association test with multiple phenotypes and multiple covariates with indepedent samples]()
* [Association test with a mixture of independent samples and sibling pairs]()
