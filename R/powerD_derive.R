#' powerD_derive.R
#' calculate the  theoretical power for SKAT test for binary outcome, uncorrelated genotypes
#' @importFrom bindata bincorr2commonprob commonprob2sigma rmvbin
#' @importFrom stats dbeta rnorm qchisq pchisq rbinom
#' @importFrom SKAT  SKAT_Null_Model
#' @param kk Number of SNPs.
#' @param n Number of subjects.
#' @param p Minor allele frequency.
#' @param alpha Significance level.
#' @param list1 Vector of effect sizes for the first alternative allele in  multi-allelic position.
#' @param list2 Vector of effect sizes for the second alternative allele in  multi-allelic position.
#'
#' @return \code{power}:Vector of derived power for SKAT test using our own method without heavy simulation.
#'
#' @export
#' @examples
#' beta1<-c(-0.9, -0.75, -0.5, -0.25, -0.1)
#' beta2<-c(rep(0.5, 5))
#' powerD_derive(k=10,n=2000,alpha=0.05,p=0.01,list1=beta1,list2=beta2)
#' # Should return 0.9991246 0.9912256 0.8957736 0.6698220 0.5704163

powerD_derive<-function(kk,n,p,alpha,list1,list2){
power<-c()
#position of first alternative allele in multi-allelic position
k1=kk-1
#position of second alternative allele in multi-allelic position
k2=kk
alpha0=-1
pi_0=exp(alpha0)/(1+exp(alpha0))

  for (ii in 1:length(list1)){
beta<-c(rep(0,kk-2),list1[ii], list2[ii])
G <- c(0, 1, 2)
#marginal probability of P(G)
prob_G <- c((1 - p)^2, 2 * p * (1 - p), p^2)
Gc <- G - 2 * p

genotypes <- expand.grid(Gk1 = G, Gk2 = G)
genotypes$prob_Gk1 <- prob_G[genotypes$Gk1 + 1]
genotypes$prob_Gk2 <- prob_G[genotypes$Gk2 + 1]
genotypes$prob <- genotypes$prob_Gk1 * genotypes$prob_Gk2
genotypes$Gck1 <- genotypes$Gk1 - 2 * p
genotypes$Gck2 <- genotypes$Gk2 - 2 * p

sigma_prime <- function(x) {
  exp(x) / (1 + exp(x))^2
}

genotypes$eta <- beta[k1] * genotypes$Gck1 + beta[k2] * genotypes$Gck2 + alpha0
genotypes$v1 <- sigma_prime(genotypes$eta)
genotypes$pi_i <- exp(genotypes$eta) / (1 + exp(genotypes$eta))
#pi_0 <- sum(genotypes$pi_i * genotypes$prob)
genotypes$mu_beta <- genotypes$pi_i - pi_0
E_v1 <- sum(genotypes$v1 * genotypes$prob)
E_mu_beta2 <- sum(genotypes$mu_beta^2 * genotypes$prob)
Var_G <- 2 * p * (1 - p)

## -----------------get A----------------##

compute_element <- function(k, l) {

  result <- 0

  for (idx in 1:nrow(genotypes)) {
    if (k == k1) {
      Gk <- genotypes$Gck1[idx]
    } else if (k == k2) {
      Gk <- genotypes$Gck2[idx]
    }

    if (l == k1) {
      Gl <- genotypes$Gck1[idx]
    } else if (l == k2) {
      Gl <- genotypes$Gck2[idx]
    }

    v1 <- genotypes$v1[idx]
    prob <- genotypes$prob[idx]

    result <- result + Gk * v1 * Gl * prob
  }
  return(result)

}

A <- matrix(0, nrow = k2, ncol = k2)

for (k in 1:k2) {
  for (l in 1:k2) {
    if ((k < k1 || l < k1) && k != l) {
      A[k, l] <- 0
    } else if (k < k1 && l < k1 && k == l) {
      A[k, l] <- Var_G * E_v1
    } else {
      A[k, l] <- compute_element(k, l)
    }
  }
}
## ----------------get B----------------##

B_matrix <- matrix(0, nrow = k2, ncol = k2)
#first k2-2 diagonal elements of B
for (k in 1:(k2-2)) {
  B_matrix[k, k] <- Var_G * E_mu_beta2 / n
}
#E(G_1_k1*v1*G_1_k1)
E_k1k1 <-(n*sum(genotypes$Gck1^2 * genotypes$mu_beta^2 * genotypes$prob)
+n*(n-1)*(sum(genotypes$Gck1 * genotypes$mu_beta * genotypes$prob))^2)
B_matrix[k1, k1] <- E_k1k1 / n^2
#E(G_1_k2*v1*G_1_k2)
E_k2k2 <- (n*sum(genotypes$Gck2^2 * genotypes$mu_beta^2 * genotypes$prob)
+n*(n-1)*(sum(genotypes$Gck2 * genotypes$mu_beta * genotypes$prob))^2)
B_matrix[k2, k2] <- E_k2k2 / n^2
#E(G_1_k1*v1*G_1_k2)
E_k1k2 <- (n*sum(genotypes$Gck1 * genotypes$mu_beta^2 * genotypes$Gck2 * genotypes$prob)
           +n*(n-1)*sum(genotypes$Gck1 * genotypes$mu_beta * genotypes$prob)
                   *sum(genotypes$Gck2 * genotypes$mu_beta * genotypes$prob))
B_matrix[k1, k2] <- E_k1k2/n^2
B_matrix[k2, k1] <- B_matrix[k1, k2]
B<-B_matrix

  c1<-rep(0,4)
  c2<-rep(0,4)
  c_a<-rep(0,4)

  A2<-A %*% A
  A3<-A2 %*% A

  c1[1] = sum(diag(A)) * n
  c1[2] = sum(A *t(A)) * n^2
  c1[3] = sum(A2 * t(A) ) * n^3
  c1[4] = sum(A2 * t(A2)) * n^4

  c2[1] = sum(diag(B)) * n^2
  c2[2] = 2 * sum(A * t(B)) * n^3
  c2[3] = 3 * sum(A2 * t(B)) * n^4
  c2[4] = 4 * sum(A3* t(B)) * n^5

    for(mm in 1:4){
      c_a[mm]<-c1[mm] + c2[mm]
    }

    #under null
    null<-lee(c1)
    quantile<-qchisq(1-alpha, df = null$l, ncp = null$delta)
    qc=(quantile-null$mux)*null$sigmaq/null$sigmax+null$muq
    #under alternative
    alt<-lee(c_a)
    temp<-(qc-alt$muq)*alt$sigmax/alt$sigmaq+alt$mux
    power[ii]<-pchisq(temp,df=alt$l,ncp=alt$delta,lower.tail = FALSE)
  }
  return(power)
}
