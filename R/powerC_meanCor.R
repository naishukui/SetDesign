#' powerC_meanCor.R
#' For continuous outcome and correalted genotypes, calculate simulated power for combined and seperated modles using SKAT, and also calculate the theoretical power using  Lee's method
#' @importFrom bindata bincorr2commonprob commonprob2sigma rmvbin
#' @importFrom stats dbeta rnorm qchisq pchisq rbinom
#' @importFrom SKAT SKAT SKAT_Null_Model
#' @param k Number of SNPs.
#' @param n Number of subjects.
#' @param alpha Significance level.
#' @param rho Correlatoin coefficient between SNPs.
#' @param p Minor allele frequency.
#' @param runs Number of simulations.
#' @param list1 Vector of effect sizes for the first alternative allele in the multi-allelic position.
#' @param list2 Vector of effect sizes for the second alternative allele in the multi-allelic position.
#'
#' @return \code{powerE,powerSeparated,powerCombined}:Vector of three powers for SKAT test.
#' @export
#' @examples
#' beta1=-0.5
#' beta2=0.75
#' powerC_meanCor(k=10,n=200,alpha=0.05,p=0.1,rho=0.1,runs=10,list1=beta1,list2=beta2)
#' # Should return powerE=0.9650208,powerSeparated=0.7,powerCombined=0.1

powerC_meanCor<-function(n,k,p,alpha,rho, runs,list1,list2){
  skat1=skat2=sig1=sig2<-c()
  sig <- matrix(0, nrow = length(list1), ncol = 2)
  powerE<-c()
  tempA=tempB<- vector("list", runs)
  result<-vector("list", length(list1))
  for (ii in seq_along(list1)){
    A=B<- matrix(0, nrow = k, ncol = k)
    for (i in 1:runs) {
      set.seed(i)
      tempA[[i]]=tempB[[i]] <- matrix(0, nrow = k, ncol = k)
      # construct correlated G matrix in separated model
      margprob<-c(rep(p,k))
      bincorr=(1-rho)*diag(k)+rho
      commonprob<-bincorr2commonprob(margprob, bincorr)
      sigma<-commonprob2sigma(commonprob)
      Ga<-rmvbin(n,margprob=diag(commonprob),sigma=sigma)
      Gb<-rmvbin(n,margprob=diag(commonprob),sigma=sigma)
      G<-Ga+Gb
      # combine the last two alleles and assume they are from the same multi-allelic position
      last <- G[, k-1] + G[, k]
      last[last > 2] <- 2
      # construct correalted G2 matrix in combined model
      G2 <- cbind(G[,1:(k-2)], last)
      #effect sizes for each SNP
      beta<-c(rep(0,k-2),list1[ii], list2[ii])
      #center G
      Gc<- t(t(G) - colMeans(G))
      mubeta<-Gc%*%beta
      # simulate linear outcome y
      y<-rnorm(n,mubeta,1)
      # calculate A and B
      tempA[[i]]<-t(Gc)%*%Gc/n
      tempB[[i]]<-t(Gc)%*%mubeta%*%t(mubeta)%*%Gc/n^2

      # separated model
      # get minor allele frequency for seperated G matrix
      MAF1<-colSums(G)/n/2
      # set weight
      w1 <- dbeta(MAF1,1,1)
      all<-as.data.frame(cbind(y,G))
      # fit null model under SKAT
      null_mod<-SKAT_Null_Model(y ~ 1, out_type="C", data=all)
      # get p value for SKAT
      skat1[i]<-SKAT(G, obj=null_mod, weights=w1)$p.value
      # get significance result
      sig1[i]<-skat1[i]< alpha

      # combined model
      MAF2<-colSums(G2)/n/2
      w2 <- dbeta(MAF2,1,1)
      all2<-as.data.frame(cbind(y,G2))
      null_mod2<-SKAT_Null_Model(y ~ 1, out_type="C", data=all2)
      skat2[i]<-SKAT(G2, obj=null_mod2, weights=w2)$p.value
      sig2[i]<-skat2[i]< alpha

    }

    # calculate theoretical power using Lee's method
    # calculate A and B using the mean of all simulation results
    A <- Reduce(`+`, tempA) / length(tempA)
    B <- Reduce(`+`, tempB) / length(tempB)
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
    q1<-(qc-alt$muq)*alt$sigmax/alt$sigmaq+alt$mux
    powerE[ii]<-pchisq(q1,df=alt$l,ncp=alt$delta,lower.tail = FALSE)
    result[[ii]]<-list( powerE=powerE[ii],powerSeparated=mean(sig1),powerCombined=mean(sig2))
  }
  return(result)
}
