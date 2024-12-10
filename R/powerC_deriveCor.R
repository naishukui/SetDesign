#' powerC_deriveCor.R
#' calculate theoretical power for SKAT test for continuous outcome, correlated genotypes
#' @importFrom bindata bincorr2commonprob commonprob2sigma rmvbin
#' @importFrom stats dbeta rnorm qchisq pchisq rbinom
#' @importFrom SKAT SKAT SKAT_Null_Model
#' @param k Number of SNPs.
#' @param n Number of subjects.
#' @param alpha Significance level.
#' @param p Minor allele frequency.
#' @param rho Correlation coefficient between SNPs.
#' @param list1 Vector of effect sizes for the first alternative allele in the multi-allelic position.
#' @param list2 Vector of effect sizes for the second alternative allele in the multi-allelic position.
#'
#' @return \code{power}:Vector of derived power for SKAT test.
#'
#' @export
#' @examples
#' beta1<-c(-0.9, -0.75, -0.5, -0.25, -0.1)
#' beta2<-c(rep(0.5, 5))
#' powerC_deriveCor(k=10,n=2000,alpha=0.05,p=0.01,rho=0.1,list1=beta1,list2=beta2)
#' # Should return 0.9962518 0.9717550 0.7914034 0.5562763 0.5183587

powerC_deriveCor<-function(k,n,alpha,p,rho,list1,list2) {
  power<-c()
   for (ii in 1:length(list1)){
# set effect sizes for each SNP
beta<-c(rep(0,k-2),list1[ii], list2[ii])
# calculate elements a to h
a=(rho*2*p*(1-p))^2
c=rho*(2*p*(1-p))^2
e=4*p^2*(1-p)^2
g=2*p*(1-p)

# conditional probability for correlated SNPs P(G11|G12=0)
p0_0=(1-p+rho*p)^2
p1_0=2*p*(1-rho)*(1-p+rho*p)
p2_0=(p-rho*p)^2
# p2_0+p1_0+p0_0==1

# conditional probability for correlated SNPs P(G11|G12=1)
p0_1=(1-p-rho+rho*p)*(1-p+rho*p)
p1_1=(p+rho-rho*p)*(1-p+rho*p)+(p-p*rho)*(1-p-rho+p*rho)
p2_1=(p+rho-rho*p)*(p-p*rho)
# p0_1+p1_1+p2_1==1

# conditional probability for correlated SNPs P(G11|G12=1)
p0_2=(1-p-rho+rho*p)^2
p1_2=2*(p+rho-rho*p)*(1-p-rho+p*rho)
p2_2=(p+rho-rho*p)^2
# p0_2+p1_2+p2_2==1

# marginal probability P(G12)
p0=(1-p)^2
p1=2*p*(1-p)
p2=p^2

# joint probability P(G11,G12)
p00=p0_0*p0
p11=p1_1*p1
p22=p2_2*p2
p20=p2_0*p0
p02=p0_2*p2
p12=p1_2*p2
p21=p2_1*p1
p10=p1_0*p0
p01=p0_1*p1

# centered G
g0=0-2*p
g1=1-2*p
g2=2-2*p

# d=E(G11^2 * G12^2)
d=e_g2g2=(g0^2*p0*(g0^2*p0_0+g1^2*p1_0+g2^2*p2_0)
          +g1^2*p1*(g0^2*p0_1+g1^2*p1_1+g2^2*p2_1)
          +g2^2*p2*(g0^2*p0_2+g1^2*p1_2+g2^2*p2_2))
# f=E(G11 * G12^3)
f=e_g1g3=(g0^3*p0*(g0*p0_0+g1*p1_0+g2*p2_0)
          +g1^3*p1*(g0*p0_1+g1*p1_1+g2*p2_1)
          +g2^3*p2*(g0*p0_2+g1*p1_2+g2*p2_2))
# b=E(G11^2 * G12 * G13)
 b=e_g2g1g1=(g0^2*g0*p00+g0^2*g1*p01+g0^2*g2*p02)*(g0*p0_0+g1*p1_0+g2*p2_0)+ (g1^2*g0*p10+g1^2*g1*p11+g1^2*g2*p12)*(g0*p0_1+g1*p1_1+g2*p2_1)+(g2^2*g0*p20+g2^2*g1*p21+g2^2*g2*p22)*(g0*p0_2+g1*p1_2+g2*p2_2)

# simulate h=E(G11 * G12 * G13 * G14)
set.seed(100)
margprob <- c(rep(p, 4))
bincorr <- (1 - rho) * diag(4) + rho
commonprob<-bincorr2commonprob(margprob, bincorr)
sigma<-commonprob2sigma(commonprob)
C <- rmvbin(n,margprob=diag(commonprob),sigma=sigma)
C1 <- C[, 1]
C2 <- C[, 2]
C3 <- C[, 3]
C4 <- C[, 4]
D <- rmvbin(n,margprob=diag(commonprob),sigma=sigma)
D1 <- D[, 1]
D2 <- D[, 2]
D3 <- D[, 3]
D4 <- D[, 4]

G1=G11 <- C1 + D1
G2=G12 <- C2 + D2
G3=G13 <- C3 + D3
G4=G14 <- C4 + D4
GG1=G1-2*p
GG2=G2-2*p
GG3=G3-2*p
GG4=G4-2*p

h=mean(GG1*GG2*GG3*GG4)

# construct A and B using elements a to h
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
