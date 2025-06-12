#' powerC_deriveCor.R
#' calculate theoretical power for SKAT test for continuous outcome, correlated genotypes
#' @importFrom bindata bincorr2commonprob commonprob2sigma rmvbin
#' @importFrom stats dbeta rnorm qchisq pchisq rbinom qnorm
#' @importFrom mvtnorm pmvnorm
#' @importFrom SKAT SKAT SKAT_Null_Model
#' @param k Number of SNPs.
#' @param n Number of subjects.
#' @param alpha Significance level.
#' @param p Minor allele frequency.
#' @param rho Correlation coefficient between SNPs.
#' @param list1 Vector of effect sizes for the first alternative allele in the multi-allelic position.
#' @param list2 Vector of effect sizes for the second alternative allele in the multi-allelic position.
#'
#' @return \code{power}:Vector of derived power for SKAT test using our own method without heavy simulation.
#'
#' @export
#' @examples
#' beta1<-c(-0.9, -0.75, -0.5, -0.25, -0.1)
#' beta2<-c(rep(0.5, 5))
#' powerC_deriveCor(k=10,n=2000,alpha=0.05,p=0.01,rho=0.1,list1=beta1,list2=beta2)
#' # Should return 0.9940736 0.9611341 0.7585572 0.5327774 0.5094248

powerC_deriveCor<-function(k,n,alpha,p,rho,list1,list2) {
  power<-c()
  for (ii in 1:length(list1)){
    # set effect sizes for each SNP
    beta<-c(rep(0,k-2),list1[ii], list2[ii])

   # calculate Φ₃ and Φ₃ and Φ₄
    margprob <- rep(p, 4)
    bincorr <- matrix(rho, nrow = 4, ncol = 4)
    diag(bincorr) <- 1
    commonprob <- bincorr2commonprob(margprob, bincorr)
    sigma <- commonprob2sigma(commonprob)
    mu <- rep(qnorm(p), 4)
    Phi4 <- pmvnorm(lower = rep(0, 4), upper = rep(Inf, 4), mean = mu, sigma = sigma)[1]
    Phi3 <- pmvnorm( lower = rep(0, 3),upper = rep(Inf, 3), mean = mu[1:3],sigma = sigma[1:3, 1:3])[1]

    # calculate elements a to h

    # a=E(G11*G12*G21*G22)
    a=(rho*2*p*(1-p))^2

    # b=E(G11^2*G12*G13)
    b<-2 * (1 - 2*p) * Phi3 + 2 * rho * p^3 * (1 - p) - 4 * (1 - 2*p) * rho * p^2 * (1 - p)- 2 * (1 - 2*p) * p^3+4 * (rho * p * (1 - p))^2

    # c=E(G11^2*G21*G22)
    c=rho*(2*p*(1-p))^2

    # d=E(G11^2*G12^2)
    #d=4 * p^2 * (1 - p)^2 + 2 * rho * p * (1 - p) * (1 - 2*p)^2
    d= 2 * (1 - 2*p)^2 * rho * p * (1 - p) + (4 + 4*rho^2) * p^2 * (1 - p)^2

    # e=E(G11^2*G21^2)
    e=4*p^2*(1-p)^2

    # f=E(G11^3*G12)
    f=2 * rho * p * (1 - p) * (3*p^2 - 3*p + 1)

    # g=E(G11^4)
    g=2*p*(1-p)

    # h=E(G11*G12*G13*G14)
    h <- 2 * Phi4 -8 * p * Phi3 + 6 * p^4 + 12 * rho * p^3 * (1 - p) + 6 * (rho * p * (1 - p))^2


    # construct matrix A and B using elements a to h
    A <- matrix(rho*2*p*(1-p), nrow = k, ncol = k)
    diag(A) <-2*p*(1-p)
    B_diag <- function(x) {
      term1 <- n * (g + (n - 1) * e) * beta[x]^2
      term2 <- n * (d + (n - 1) * a) * sum(beta[-x]^2)
      term3 <- 2 * n * (f + (n - 1) * c) * sum(beta[-x] * beta[x])
      term4 <-  n * (b + (n - 1) * a) * (sum(outer(beta[-x], beta[-x], "*"))-sum(beta[-x]^2))

      result <- (term1 + term2 + term3 + term4 )/n^2
      return(result)
    }

    B_offdiag <- function(x,y) {
      term1 <- (n * f + n * (n - 1) * c) * (beta[x]^2 + beta[y]^2)
      term2 <- (n * b + n * (n - 1) * a) * sum(beta[-c(x, y)]^2)
      term3 <- (2*n*d+n * (n - 1) * ( a + e)) * beta[x] * beta[y]
      term4 <- (2*n*b+n * (n - 1) * (a + c)) * sum(outer(beta[-c(x, y)], beta[x], "*"),outer(beta[-c(x, y)], beta[y], "*"))
      term5<-n*h+n*(n-1)*a*(sum(outer(beta[-c(x,y)], beta[-c(x,y)], "*"))-sum(beta[-c(x,y)]^2))

      result <- (term1 + term2 + term3 + term4 + term5)/n^2
      return(result)
    }


    B <- matrix(0, k,k)

    for (x in 1:k) {
      for (y in 1:k) {
        if (x == y) {
          B[x, y] <- B_diag(x)
        } else {
          B[x, y] <- B_offdiag(x,y)
        }
      }
    }
    suml=sumld=cnull=calt<- numeric(4)



    A2<-A %*% A
    A3<-A2 %*% A

    c1<-rep(0,4)
    c2<-rep(0,4)
    c_a<-rep(0,4)
    c1[1] = sum(diag(A)) * n
    c1[2] = sum(A *t(A)) * n^2
    c1[3] = sum(A2 * t(A) ) * n^3
    c1[4] = sum(A2 * t(A2)) * n^4

    c2[1] = sum(diag(B)) * n^2
    c2[2] = 2 * sum(A * t(B)) * n^3
    c2[3] = 3 * sum(A2 * t(B)) * n^4
    c2[4] = 4 * sum(A3* t(B)) * n^5

    for(i in 1:4){
      c_a[i]<-c1[i] + c2[i]
    }
    #under null
    null<-lee(c1)
    q0<-qchisq(1-alpha, df = null$l, ncp = null$delta)
    qc=(q0-null$mux)*null$sigmaq/null$sigmax+null$muq
    #under alternative
    alt<-lee(c_a)
    temp<-(qc-alt$muq)*alt$sigmax/alt$sigmaq+alt$mux
    power[ii]<-pchisq(temp,df=alt$l,ncp=alt$delta,lower.tail = FALSE)
  }
  return(power)
}
