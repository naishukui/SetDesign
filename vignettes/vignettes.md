

# What is SetDesign?

SetDesign was designed to study the impact of potential data processing choices at study design stages by performing power and bias calculation of Sequece Kernal Association Test (SKAT)  when multi-allelic positions are present. The package supports both linear and binary outcomes and accounts for correlated SNPs by specifying a correlation coefficient for more precise results.


 It compares the simulated and analytical power and bias between:

- **True Model**: account for two alternative allelesfor  tri-allelic loci.
- **Misspecified Model**: Treating tri-allelic loci as a bi-allelic loci.

The package also developed a more efficient method for analystical power calculation based on the original power calculation for SKAT by Lee et al, making power calculatin more efficient and versatile. 

Users can customize parameters such as:

- Number of subjects and SNPs.\
- Correlation between SNPs.\
- Significance level for the SKAT test.\
- Effect sizes of the two alternative alleles.\




# Why use SetDesign?

Traditional methods for SKAT power calculations often oversimplify multi-allelic positions by treating them as single alleles—either considering only the first alternative allele or combining the first and second alternative alleles in a simplified manner. SetDesign highlights the differences in power and effect size estimation when these positions are treated individually versus when they are combined.

Additionally, the package builds on and improves the power calculation method from Lee et al.'s paper, significantly reducing computational complexity by directly specifying the correlation structure of the genotype matrix.


# Example

Suppose we want to calculate the power of the SKAT test for a binary outcome, assuming correlated SNPs and one multi-allelic position among `kk` SNPs.

``` r
library(SetDesign)

# Set number of  SNPs
kk = 50

# Set number of subjects
n = 2000

# Set significance level
alpha = 0.05

# Set correlation coefficient between SNPs
rho = 0.15

# Set minor allele frequency
p = 0.01

# Set effect sizes of the first alternative allele
list1=c(-0.1,-0.25,-0.5,-0.75,-0.9)

# Set effect sizes of the second alternative allele
list2=c(rep(0.5,5))

# Pairwise correlation
rho = 0.1

# get the power for SKAT test
powerD_derive(kk,n,p,alpha,list1,list2)
#> [1] 0.08653423 0.09248395 0.11143999 0.13961267 0.15997840
```



Next we show how to get the parameters in misspecifed model  given the parameters in true models, showing the bias resulted from model misspecificatoin where tri-allelic variants are treated as bi-allelic. 

Suppose we have the effect sizes and allele frequencies for true model, we can estimate the bias in effect size estimation in misspecifed model. For binary outcomes: 

``` r

#Set the true coefficients (intercept, minor allele 1, minor allele 2)
betaVec = c(-1, -0.5, 1)

#Set the probabilities of each of the three alleles (major allele, minor 1, minor 2)
probVec=c(0.98, 0.01, 0.01)

#get a List containing solved coefficients (alpha0, alphaG) for misspecified model
bias_logistic_nleqslv(betaVec = c(-1, -0.5, 1), probVec=c(0.98, 0.01, 0.01))$x
```
For linear outcomes: 

```r
#set beta_0 (intercept), beta_X (covariate effects), beta_G (genetic effects)
beta_list <- list(beta_0 = 1, beta_X = c(0.5, -0.3), beta_G = c(0.2, 0.4))

#set  mu_X (mean of X)
mu_list <- list(mu_X = c(2, 1.5))

#set cov_mat_list : A list of covariance matrices
 cov_mat_list <- list(
   mu_XX = matrix(c(5, 3, 3, 4), nrow = 2),
   mu_X_tildeG = matrix(c(0.1, 0.2, 0.3, 0.4), nrow = 2),
   mu_XG = matrix(c(0.5, 0.6, 0.7, 0.8), nrow = 2),
   mu_tildeG_tildeG = matrix(c(1, 0.5, 0.5, 1), nrow = 2),
   mu_tildeG_G = matrix(c(0.9, 0.8, 0.7, 0.6), nrow = 2))

#get A list with solutions for α₀, α_x, and α_G
bias_linear_nleqslv(beta_list, cov_mat_list, mu_list)
```
