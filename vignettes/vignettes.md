

# What is SetDesign?

SetDesign was developed to assess the impact of data processing choices at the study design stage of genome-wide association studies (GWAS). It focuses on misspecified models, which can arise due to differences in data pipelines. The package performs power calculations for the sequence kernel association test (SKAT) or bias calculations for estimated regression coefficients when models are misspecified. The package supports both continuous and binary outcomes and can account for correlated SNPs.

 SetDesign can perform both simulated and analytical calculations for power and bias between:

- **True Model**: (e.g. modeling a tri-allelic variant with two terms)
- **Misspecified Model**: (e.g. modeling a tri-allelic variant with only one term)

The package also introduces a more computationally efficient method for analytical power calculation that does not require actual genotype matrices or bootstraping a full model.

Users can customize parameters such as:

- Number of subjects and SNPs.
- Correlation between SNPs.
- Significance level for the SKAT test.
- Effect sizes of alleles.


# Why use SetDesign?

Researchers make many important - but often unreported - data processing choices in set-based genetic association analysis. These choices can greatly affect the stability and reproducibility of set-based results. SetDesign helps researchers easily understand the (potential) large impacts of their choices at the start of their study.

# Example

Suppose we want to calculate the power of the SKAT test for a binary outcome in a misspecified model. Assume we are interested in the situation where there are `kk=50` total SNPs and `n=2000` subjects. Assume further that we want to consider the setting where the one true causal SNP is a tri-allelic SNP where both non-reference alleles have a non-zero effect. However, we misspecify the model by modeling the causal SNP as a single covariate and simply summing the number of non-reference alleles. 

``` r
library(SetDesign)

# Set number of  SNPs
kk = 50

# Set number of subjects
n = 2000

# Set significance level
alpha = 0.05

# Set minor allele frequency
p = 0.01

# Set true effect sizes of the first alternative allele
list1=c(-0.1,-0.25,-0.5,-0.75,-0.9)

# Set true effect sizes of the second alternative allele
list2=c(rep(0.5,5))

# get the power for SKAT test under the misspecified model
powerD_derive(kk,n,p,alpha,list1,list2)
#> [1] 0.08653423 0.09248395 0.11143999 0.13961267 0.15997840
```


Next we consider estimating regression coefficients under model misspecification. We show how to get the regression coefficient estimates in the misspecified model when tri-allelic variants are treated as bi-allelic. 
 For binary outcomes: 

``` r

# Set the true coefficients (intercept, coefficient 1, coefficient 2)
betaVec = c(-1, -0.5, 1)

# Set the probabilities of each of the three alleles (reference allele, effect alelle 1, effect allele 2)
probVec=c(0.98, 0.01, 0.01)

# return a list containing solved coefficients (alpha0, alphaG) for misspecified model
bias_logistic_nleqslv(betaVec = c(-1, -0.5, 1), probVec=c(0.98, 0.01, 0.01))$x
```
For linear outcomes: 

```r
# set beta_0 (intercept), beta_X (covariate effects), beta_G (genetic effects)
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
