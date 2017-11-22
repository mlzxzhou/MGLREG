# MGLREG  
**MGLREG** develops a general regression framework for the analysis of secondary phenotypes collected in multiple group association studies. Our regression framework is built on a conditional model for the secondary outcome given the multiple group status and covariates and its relationship with the population regression of interest of the secondary outcome given covariates. Then, we develop generalized estimating equations to estimate parameters of interest. We use simulations and a large-scale imaging genetic data analysis of the ADNI data to evaluate the effects of the multiple group sampling scheme on genome-wide association study (GWAS) results based on some standard statistical methods, such as linear regression methods, while comparing it with our statistical methods that appropriately adjust for the multiple group sampling scheme.

## Features
* It considers the multiple-group case instead of the traditinal case-control status
* We also include the regression model to deal with the two-phase sample other than the simple single-phase study.


## Installation
To install the stable version from CRAN, simply run the following from an R console: **We have not uploaded the package to CRAN and we expect it will be avaiable in the near future.**

```r
install.packages("MGLREG")
```

To install the latest development builds directly from GitHub, run this instead:

```r
if (!require("devtools"))
  install.packages("devtools")
devtools::install_github("./MGLREG")
```

## Using MGLREG
Function 'MGLReg' works for the single-phase three-group sample. By running the example we provide in the document 'MGLReg.Rd', you could exactly copy the simulation study of the reference shown below. On the other hand, we also offer a function called 'MGLReg2', which could handle the two-phase situation, which may be used when you want to do real data analysis on some two-phase studies such as ADNI.  

* Zhou, F., Zhou, H., Li, T. & Zhu, HT. Analysis of Secondary Phenotypes in Multiple-group Association Studies


## Simulation
Below is the example of how to replicate the simulation result given in the paper when the simulated SNP has MAF equal to 0.3. By running the code, you are able to generate the whole population data with three groups and sample 800 observations with proportions (0.25, 0.5, 0.25), and then achieve the unbiased estimator og genetic factor and corresponding p-value.

```r
library("MGLREG")
# suppose the size of the whole population is 100000 with three subgruops and the proportions of the three groups are (10%, 15%, 75%)
# Generate genetic covariate G with MAF=0.2
MAF <- 0.3
G <- rbinom(100000,2,MAF)
s <- 800
beta <- c(5, 2, 1)
gamma1 <- c(-2,-1,-1)
gamma2 <- c(1,1,1)
# Generate X, which is mixture of normal variables
x_p <- c(0.2,0.8)
x_mean <- c(0,2)
x_1 <- rnorm(20000, mean = 2, sd = 2)
x_2 <- rnorm(80000, mean=0, sd=2)
X <- c(x_1,x_2)
data <- generate_data3(beta, gamma1, gamma2, X, G, s)
reg <- MGLReg(data[,2],data[,3],data[,4],data[,1],c(0.1,0.15,0.75))
bias <- abs(reg[[1]]-beta[2]) ## the bias by using our approach
bias_lm <- abs(summary(lm(data[,2]~data[,3:4]))$coefficients[2]-beta[2]) ## the bias by directly running linear model as comparison
```

## Real data analysis
We also include a small demo to play, which replicates the GWAS study on the volume of right hippocampus described in the paper. For simplicity, we only include the SNPs from Chromosome 22, and the required genetic and phenotype data could be achieved by loading the Rdata contained in the package. The output will be a mtrix where the first column corresponds to the SNP names from Chromosome 22 and the second comlumn is the p-value of the SNP effect on the volume of right hipopcampus. It may take some time for the whole running, so you may consider trying it on the cluster. 

```r
p_til=c(0.8,0.1,0.1) # the pre-set global prevelance of three groups
data('real_data')
total_snp=dim(Geno)[2]
p_value=matrix(0,total_snp,2)
p_value[,1]=colnames(Geno)
for (i in 1:total_snp){
 p_value[i,2]=MGLReg2(y,Geno[,i],C,D,p,p_til)[[2]]
}
```
