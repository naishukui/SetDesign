#' powerC_derive.R
#' calculate theoretical power for SKAT test for continuous outcome, uncorrelated genotypes.
#' @importFrom bindata bincorr2commonprob commonprob2sigma rmvbin
#' @importFrom stats dbeta rnorm qchisq pchisq rbinom
#' @importFrom SKAT SKAT SKAT_Null_Model
#' @param k Number of SNPs.
#' @param n Number of subjects.
#' @param p Minor allele frequency.
#' @param alpha Significance level.
#' @param list1 Vector of effect sizes for the first alternative allele in the multi-allelic position.
#' @param list2 Vector of effect sizes for the second alternative allele in the multi-allelic position.
#'
#' @return \code{power}:Vector of derived power for SKAT test.
#'
#' @export
#' @examples
#' beta1<-c(-0.9, -0.75, -0.5, -0.25, -0.1)
#' beta2<-c(rep(0.5, 5))
#' powerC_derive(k=10,n=2000,alpha=0.05,p=0.01,list1=beta1,list2=beta2)
#' # Should return 0.9991246 0.9912256 0.8957736 0.6698220 0.5704163

powerC_derive<-function(k,n,alpha,p,list1,list2){
  power<-c()
  for (ii in 1:length(list1)){
    # effect size of all snp
    beta<-c(rep(0,k-2),list1[ii], list2[ii])
    # theoretical value of A matrix based on derivation
    A<-diag(2*p*(1-p),nrow=k,ncol=k)
    # theoretical value of B matrix based on derivation
    B<-diag((2*n*p*(1-3*p+4*p^2-2*p^3)+4*n^2*(p^2-2*p^3+p^4))*beta^2/n^2,nrow=k,ncol=k)
    A2<-A %*% A
    A3<-A2 %*% A
    #  under the null distribution
    c1<-rep(0,4)
    c2<-rep(0,4)
    #  under the alternative distribution
    c_a<-rep(0,4)

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

