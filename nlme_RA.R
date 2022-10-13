source("RA_functions.R")


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

nIndep <- 50000; geno_prob <- c(0.64, 0.32, 0.04)
g_test <- sample(0:2, size=nIndep, replace=T, prob=geno_prob)
y_test <- cbind(replicate(3, rnorm(nIndep)), sample(0:1, size=nIndep, replace=T, prob=c(0.8, 0.2)))
z_test <- cbind(replicate(5, rnorm(nIndep, mean=3, sd=3)), replicate(10, rnorm(nIndep, mean=4, sd=4)))

## Independent samples, 3 continuous and 1 binary phenotype traits, no covariates, no assumption of HWE
RA_assoc_indep(g=g_test, y=y_test, z=NULL, HWE=F, nlme = F)
RA_assoc_indep(g=g_test, y=y_test, z=NULL, HWE=F, nlme = T)

## Independent samples, 3 continuous and 1 binary phenotype traits, 4 covariates, no assumption of HWE
RA_assoc_indep(g=g_test, y=y_test, z=z_test, HWE=F, nlme = F)
RA_assoc_indep(g=g_test, y=y_test, z=z_test, HWE=F, nlme = T)


install.packages("PhenotypeSimulator")
library(PhenotypeSimulator)
# simulate phenotype with the same five phenotype components
# and settings as above; display progress via verbose=TRUE
# total SNP effect on phenotype: 0.01
# read genotypes from external file
# use one of the sample genotype file provided in the 
# extdata/genotypes/subfolders (e.g.extdata/genotypes/hapgen )

PhenotypeSimulator::simulateGenotypes()

genotypefile <- system.file("extdata/genotypes/hapgen", 
                            "genotypes_hapgen.controls.gen", 
                            package = "PhenotypeSimulator")
# remove the .gen ending (oxgen specific endings .gen and .sample are added 
# automatically )
genotypefile <- gsub("\\.gen","", genotypefile)
genVar <- 0.6
noiseVar <- 1 - genVar
totalSNPeffect <- 0.25
h2s <- totalSNPeffect/genVar
phi <- 0.6 
rho <- 0.1
delta <- 0.3
shared <- 0.8
independent <- 1 - shared
phenotype <- runSimulation(N = 10000, P = 15, 
                           # genotypefile = genotypefile, 
                           format = "oxgen", cNrSNP = 30, genVar = genVar, h2s = h2s, 
                           phi = 0.6, delta = 0.3, distBetaGenetic = "unif", mBetaGenetic = 0.5, 
                           sdBetaGenetic = 1, NrFixedEffects = 4, NrConfounders = c(1, 2, 1, 2), 
                           pIndependentConfounders = c(0, 1, 1, 0.5), 
                           distConfounders = c("bin", "cat_norm", "cat_unif", "norm"), 
                           probConfounders = 0.2, catConfounders = c(3, 4), pcorr = 0.8, 
                           verbose = TRUE)

y_test <- phenotype$phenoComponentsFinal$Y
g_test <- phenotype$phenoComponentsIntermediate$genFixed$cov
phenotype$phenoComponentsIntermediate$cov_ge
# str(phenotype)

RA_assoc_indep(g=g_test[,2], y=y_test, z=NULL, HWE=F, nlme = F)
RA_assoc_indep(g=g_test[,2], y=y_test, z=NULL, HWE=F, nlme = T)

#> Set seed: 219453
#> The total noise variance (noiseVar) is: 0.4
#> The noise model is: noiseFixedAndBgAndCorrelated
#> Proportion of non-genetic covariate variance (delta): 0.3
#> Proportion of variance of shared non-genetic covariate effects (gamma): 0.8
#> Proportion of non-genetic covariates to have a trait-independent effect (pIndependentConfounders ): 0 1 1 0.5
#> Proportion of traits influenced by independent non-genetic covariate effects (pTraitIndependentConfounders): 0.2
#> Proportion of variance of correlated noise effects (rho): 0.1
#> Proportion of observational noise variance (phi): 0.6
#> Variance of shared observational noise effect (alpha): 0.8
#> 
#> The total genetic variance (genVar) is: 0.6
#> The genetic model is: geneticFixedAndBg
#> Proportion of variance of genetic variant effects (h2s): 0.0166666666666667
#> Proportion of variance of shared genetic variant effects (theta): 0.8
#> Proportion of genetic variant effects to have a trait-independent fixed effect (pIndependentGenetic): 0.4
#> Proportion of traits influenced by independent genetic variant effects (pTraitIndependentGenetic): 0.2
#> Proportion of variance of infinitesimal genetic effects (h2bg): 0.983333333333333
#> Proportion of variance of shared infinitesimal genetic effects (eta): 0.8
#> Proportion of non-linear phenotype transformation (proportionNonlinear): 0
#> 
#> Simulate genetic effects (genetic model: geneticFixedAndBg)
#> Simulate genetic variant effects
#> Out of 15 total phenotypes, 15 traits will be affected by genetic variant effects
#> Out of these affected traits (15), 3 trait(s) will have independent genetic variant effects
#> Standardising the 1000 SNPs provided
#> Estimating kinship from 1000 SNPs provided
#> Normalising kinship
#> Simulate infinitesimal genetic effects
#> Simulate noise terms (noise model: noiseFixedAndBgAndCorrelated)
#> Simulate correlated background effects
#> Simulate observational noise effects
#> Simulate confounder effects
#> Out of 15 total phenotypes, 15 trait(s) will be affected by the  1  covariate effect
#> Out of 15 total phenotypes, 15 trait(s) will be affected by the  2  covariate effect
#> Out of these affected traits (15), 3 trait(s) will have independent covariate effects
#> Out of 15 total phenotypes, 15 trait(s) will be affected by the  3  covariate effect
#> Out of these affected traits (15), 3 trait(s) will have independent covariate effects
#> Out of 15 total phenotypes, 15 trait(s) will be affected by the  4  covariate effect
#> Out of these affected traits (15), 3 trait(s) will have independent covariate effects
#> Construct final simulated phenotype










library(bench)
bench::mark(
  RA = RA_assoc_indep(g=g_test, y=y_test, z=NULL, HWE=F, nlme = F),
  nlme = RA_assoc_indep(g=g_test, y=y_test, z=NULL, HWE=F, nlme = T),
  iterations = 10,
  check = FALSE
)

simdf <- function(n){
  nIndep <- n; geno_prob <- c(0.64, 0.32, 0.04)
  g_test <- sample(0:2, size=nIndep, replace=T, prob=geno_prob)
  y_test <- cbind(replicate(3, rnorm(nIndep)), sample(0:1, size=nIndep, replace=T, prob=c(0.8, 0.2)))
  z_test <- cbind(replicate(2, rnorm(nIndep, mean=3, sd=3)), replicate(2, rnorm(nIndep, mean=4, sd=4)))
  list(g_test = g_test, y_test = y_test, z_test = z_test)
}

simdf(10)$z_test

result <- bench::press(
  n = c(1e2,1e3,1e4,1e5,4e5), 
  {
    dat <- simdf(n)
    g_test <- dat$g_test
    y_test <- dat$y_test
    z_test <- dat$z_test
    bench::mark(
      RA = RA_assoc_indep(g=g_test, y=y_test, z=NULL, HWE=F, nlme = F),
      nlme = RA_assoc_indep(g=g_test, y=y_test, z=NULL, HWE=F, nlme = T),
      iterations = 10,
      check = FALSE
      )
  }
)


