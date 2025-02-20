
# What is SetDesign?

SKATmultiallele was designed for power calculation of Sequece Kernal
Association Test (SKAT) when there are multi-allele positions exist.
Users can specify the number of subjects, SNPs and the minor allele
frequency, significance level and effect size of both alternative
alleles to calculate the power of SKAT test. The package can handle both
linear and binary outcomes. For correlated SNPs, users can specify the
correlation coefficient between SNPs for a more accurate result.

# Why use SetDesign?

In many traditional ways of SKAT power calculation, multi-allelic
position are often treated as single allele where only the first
alternative allele are considered, or the first and second alternative
allele are combined in some way. We want to show the difference of power
when doing so, as opposed to treating them as separated alleles.
Further, we improved the power calculation method based on Lee et al’s
paper to decrease the computational burden and made it possible to
calculate the power by simply specifying the correaltion structure of
genotype matrix.

# Example

We show an example with binary outcomes and correlated SNPs

``` r
library(SKATmultiallele)

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
