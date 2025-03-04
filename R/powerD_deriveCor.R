#' powerD_deriveCor.R
#' calculate theoretical power for SKAT test for binary outcome, correlated genotypes
#' @importFrom bindata bincorr2commonprob commonprob2sigma rmvbin
#' @importFrom stats dbeta rnorm qchisq pchisq rbinom
#' @importFrom SKAT SKAT SKAT_Null_Model
#' @param kk Number of SNPs.
#' @param n Number of subjects.
#' @param p Minor allele frequency.
#' @param rho Correlation coefficient between SNPs.
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
#' powerD_deriveCor(k=10,n=2000,alpha=0.05,p=0.01,rho=0.15,list1=beta1,list2=beta2)
#' # Should return 0.25564418 0.15141926 0.08757907 0.13321503 0.21503152

powerD_deriveCor<-function(kk,n,alpha,p,rho,list1,list2) {
  power<-c()
  k1=kk-1
  k2=kk
  alpha0=-1
  pi_0=exp(alpha0)/(1+exp(alpha0))　

  for (ii in 1:length(list1)){
    #effect sizes for each SNPs
    beta<-c(rep(0,kk-2),list1[ii], list2[ii])
    G_value <- c(0, 1, 2)
    prob_G <- c((1 - p)^2, 2 * p * (1 - p), p^2)
    #centered G
    Gc <- G_value - 2 * p

    genotypes <- expand.grid(Gk1 = G_value, Gk2 = G_value)
    genotypes$prob_Gk1 <- prob_G[genotypes$Gk1 + 1]
    genotypes$prob_Gk2 <- prob_G[genotypes$Gk2 + 1]
    genotypes$prob <- genotypes$prob_Gk1 * genotypes$prob_Gk2

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

    genotypes$probrho<-c(p00,p10,p20,
                         p01,p11,p21,
                         p02,p12,p22)

    genotypes$Gck1 <- genotypes$Gk1 - 2 * p
    genotypes$Gck2 <- genotypes$Gk2 - 2 * p

    sigma_prime <- function(x) {
      exp(x) / (1 + exp(x))^2
    }

    genotypes$eta <- beta[k1] * genotypes$Gck1 + beta[k2] * genotypes$Gck2 + alpha0
    genotypes$v1 <- sigma_prime(genotypes$eta)
    genotypes$pi_i <- exp(genotypes$eta) / (1 + exp(genotypes$eta))
    # pi_0 <- sum(genotypes$pi_i * genotypes$prob)
    genotypes$mu1<- genotypes$pi_i - pi_0
    genotypes$mu_beta <- genotypes$pi_i - pi_0
    E_v1 <- sum(genotypes$v1 * genotypes$prob)
    E_mu_beta2 <- sum(genotypes$mu_beta^2 * genotypes$prob)
    Var_G <- 2 * p * (1 - p)

    E_k1k1<-sum(genotypes$Gck1^2 * genotypes$v1 * genotypes$probrho)
    E_k2k2<-sum(genotypes$Gck2^2 * genotypes$v1 * genotypes$probrho)
    E_k1k2<-sum(genotypes$Gck1 * genotypes$Gck2^2 * genotypes$v1 * genotypes$probrho)
    # simulate some probabilities
    set.seed(100)
    margprob <- c(rep(p, 4))
    bincorr <- (1 - rho) * diag(4) + rho
    commonprob<-bincorr2commonprob(margprob, bincorr)
    sigma<-commonprob2sigma(commonprob)
    C <- rmvbin(n*1000,margprob=diag(commonprob),sigma=sigma)
    D <- rmvbin(n*1000,margprob=diag(commonprob),sigma=sigma)
    G<-C+D
    G11 <- G[, 1]
    G12 <- D[, 2]
    G13 <- D[, 3]
    G14 <- D[, 4]

  #---------------------------------------------------------------------------------------------------------------get A
    # P(G13, G14 | G11)
    joint_con_prob1 <- function(G11_value) {
      subset_G <- G[G[,1] == G11_value, ]
      joint_con_prob <- matrix(0, nrow = 3, ncol = 3)
      for (G13_value in 0:2) {
        for (G14_value in 0:2) {
          count_G13_G14 <- sum(subset_G[,3] == G13_value & subset_G[,4] == G14_value)
          joint_con_prob[G13_value+1, G14_value+1] <- count_G13_G14 / nrow(subset_G)
        }
      }

      return(joint_con_prob)
    }

    # P(G13, G14 | G11, G12)
    joint_con_prob2 <- function(G11_value, G12_value) {
      subset_G <- G[G[,1] == G11_value & G[,2] == G12_value, ]
      joint_con_prob <- matrix(0, nrow = 3, ncol = 3)
      for (G13_value in 0:2) {
        for (G14_value in 0:2) {
          count_G13_G14 <- sum(subset_G[,3] == G13_value & subset_G[,4] == G14_value)
          joint_con_prob[G13_value+1, G14_value+1] <- count_G13_G14 / nrow(subset_G)
        }
      }

      return(joint_con_prob)
    }

    # E(v1|G11,G12)/E(v1|G11)
    inner_e0 <- function(joint_con_prob) {
      inner_e <- 0
      for (G13_value in 0:2) {
        for (G14_value in 0:2) {
          # P(G13, G14 | G11, G12) or P(G13, G14 | G11)
          prob_G13_G14 <- joint_con_prob[G13_value+1, G14_value+1]
          v1 <-   exp(alpha0 + G13_value * beta[k1] + G14_value * beta[k2]) / (1 + exp(alpha0 + G13_value * beta[k1] + G14_value * beta[k2]))^2
          inner_e <- inner_e + v1 * prob_G13_G14

        }
      }
      return(inner_e)
    }

    # E(G13*v1*G11|G11)
    inner_e13 <- function(joint_con_prob) {
      inner_e <- 0
      for (G13_value in 0:2) {
        for (G14_value in 0:2) {
          # P(G13, G14 | G11, G12)
          prob_G13_G14 <- joint_con_prob[G13_value+1, G14_value+1]
          v1 <-   exp(alpha0 + G13_value * beta[k1] + G14_value * beta[k2]) / (1 + exp(alpha0 + G13_value * beta[k1] + G14_value * beta[k2]))^2
          inner_e <- inner_e + (G13_value-2*p)*v1 * prob_G13_G14

        }
      }
      return(inner_e)
    }

    # E(G14*v1*G11|G11)
    inner_e14 <- function(joint_con_prob) {
      inner_e <- 0
      for (G13_value in 0:2) {
        for (G14_value in 0:2) {
          # P(G13, G14 | G11, G12)
          prob_G13_G14 <- joint_con_prob[G13_value+1, G14_value+1]
          v1 <-   exp(alpha0 + G13_value * beta[k1] + G14_value * beta[k2]) / (1 + exp(alpha0 + G13_value * beta[k1] + G14_value * beta[k2]))^2
          inner_e <- inner_e + (G14_value-2*p)*v1 * prob_G13_G14

        }
      }
      return(inner_e)
    }

    # E(G11*G11*v1)
    outer_e11<- function() {
      outer_e <- 0
      for (Gi_val in 0:2) {
        # P(G11)
        prob_Gi <-prob_G[Gi_val+1]
        # p(G13,G14|G11)
        joint_con_prob <- joint_con_prob1(Gi_val)
        #  E(v1 | G11)
        inner_e <- inner_e0(joint_con_prob)
        # E(G11*G11*v1)
        outer_e <- outer_e + (Gi_val-2*p)^2 * prob_Gi * inner_e
      }
      return(outer_e)
    }

    # E(G11*G12*v1)
    outer_e12<- function() {
      outer_e <- 0
      for (Gi_val in 0:2) {
        for (Gj_val in 0:2) {
          # P(G11, G12)
          prob_Gi_Gj <- genotypes$probrho[genotypes$Gk1==Gi_val & genotypes$Gk2==Gj_val]
          # p(G13,G14|G11,G12)
          joint_con_prob <- joint_con_prob2(Gi_val, Gj_val)
          #  E(v1 | G11, G12)
          inner_e <- inner_e0(joint_con_prob)
          # E(G11*G12*v1)
          outer_e <- outer_e + (Gi_val-2*p) * (Gj_val-2*p)  * prob_Gi_Gj * inner_e
        }
      }
      return(outer_e)
    }

    # E(G11*G13*v1)
    outer_e13<- function() {
      outer_e <- 0
      for (Gi_val in 0:2) {
        # P(G11)
        prob_Gi <-prob_G[Gi_val+1]
        # p(G13,G14|G11)
        joint_con_prob <- joint_con_prob1(Gi_val)
        #  E(G13*v1 | G11)
        inner_e <- inner_e13(joint_con_prob)
        # E(G11*G13*v1)
        outer_e <- outer_e + (Gi_val-2*p)  * prob_Gi * inner_e
      }
      return(outer_e)
    }

    # E(G11*G14*v1)
    outer_e14<- function() {
      outer_e <- 0
      for (Gi_val in 0:2) {
        # P(G11)
        prob_Gi <-prob_G[Gi_val+1]
        # p(G13,G14|G11)
        joint_con_prob <- joint_con_prob1(Gi_val)
        # E(G14*v1 | G11)
        inner_e <- inner_e14(joint_con_prob)
        # E(G11*G14*v1)
        outer_e <- outer_e + (Gi_val-2*p)  * prob_Gi * inner_e
      }
      return(outer_e)
    }

    generate_matrix_X <- function(a, b, c, d, e, f, g) {
      # Initialize a kk x kk matrix with zeros
      X <- matrix(0, nrow = kk, ncol = kk)

      # Fill in the diagonal elements based on the conditions
      for (k in 1:k2) {
        if (k < k1) {
          X[k, k] <- a  # For k = 1, 2: X[k, k] = a
        } else if (k == k1) {
          X[k, k] <- b  # For k = 3: X[k, k] = b
        } else if (k == k2) {
          X[k, k] <- c  # For k = 4: X[k, k] = c
        }
      }

      # Fill in the off-diagonal elements
      for (k in 1:k2) {
        for (l in 1:k2) {
          if (k != l) {
            if (k < k1 && l < k1) {
              X[k, l] <- d
            } else if ( (k == k1 & l<k2) || (l == k1 && k<k2)) {
              X[k, l] <- e
            } else if (k == k2 || l == k2) {
              if ((k == k2 && l == k1) || (l == k2 && k == k1)) {
                X[k, l] <- g
              } else {
                X[k, l] <- f  # For k = k2 or l = k2
              }
            }
          }
        }
      }

      return(X)
    }


    a <- outer_e11()  # Diagonal element when k < k1
    b <- E_k1k1  # Diagonal element when k = k1
    c <- E_k2k2  # Diagonal element when k = k2
    d <- outer_e12()  # Off-diagonal element when k < k1 and l < k1
    e <- outer_e13()   # Off-diagonal element when k = k1 or l = k1
    f <- outer_e14()   # Off-diagonal element when k = k2 or l = k2
    g <- E_k1k2   # Special case for (k = k2, l = k1) or (k = k1, l = k2)

    A <- generate_matrix_X(a, b, c, d, e, f, g)

    #---------------------------------------------------------------------------------------------------------------get B
    # E(mu1|G11,G12)/E(v1|G11)
    inner_e0<- function(joint_con_prob) {
      inner_e <- 0
      for (G13_value in 0:2) {
        for (G14_value in 0:2) {
          # P(G13, G14 | G11, G12) or P(G13, G14 | G11)
          prob_G13_G14 <- joint_con_prob[G13_value+1, G14_value+1]
          pi_1 <-   exp(alpha0 + G13_value * beta[k1] + G14_value * beta[k2]) / (1 + exp(alpha0 + G13_value * beta[k1] + G14_value * beta[k2]))
          mu1=pi_1-pi_0
          inner_e <- inner_e + mu1 * prob_G13_G14
        }
      }
      return(inner_e)
    }

    # E(mu1^2|G11,G12)/E(v1|G11)
    inner_e0sq <- function(joint_con_prob) {
      inner_e <- 0
      for (G13_value in 0:2) {
        for (G14_value in 0:2) {
          # P(G13, G14 | G11, G12) or P(G13, G14 | G11)
          prob_G13_G14 <- joint_con_prob[G13_value+1, G14_value+1]
          pi_1 <-   exp(alpha0 + G13_value * beta[k1] + G14_value * beta[k2]) / (1 + exp(alpha0 + G13_value * beta[k1] + G14_value * beta[k2]))
          mu1=pi_1-pi_0
          inner_e <- inner_e + mu1^2 * prob_G13_G14
        }
      }
      return(inner_e)
    }


    # E(G13*mu1^2|G11)
    inner_e13 <- function(joint_con_prob) {
      inner_e <- 0
      for (G13_value in 0:2) {
        for (G14_value in 0:2) {
          # P(G13, G14 | G11, G12)
          prob_G13_G14 <- joint_con_prob[G13_value+1, G14_value+1]
          pi_1 <-   exp(alpha0 + G13_value * beta[k1] + G14_value * beta[k2]) / (1 + exp(alpha0 + G13_value * beta[k1] + G14_value * beta[k2]))
          mu1=pi_1-pi_0
          inner_e <- inner_e + (G13_value-2*p)*mu1^2 * prob_G13_G14
        }
      }
      return(inner_e)
    }


    # E(G14*mu1^2|G11)
    inner_e14<- function(joint_con_prob) {
      inner_e <- 0
      for (G13_value in 0:2) {
        for (G14_value in 0:2) {
          # P(G13, G14 | G11, G12)
          prob_G13_G14 <- joint_con_prob[G13_value+1, G14_value+1]
          pi_1 <-   exp(alpha0 + G13_value * beta[k1] + G14_value * beta[k2]) / (1 + exp(alpha0 + G13_value * beta[k1] + G14_value * beta[k2]))
          mu1=pi_1-pi_0
          inner_e <- inner_e + (G14_value-2*p)*mu1^2 * prob_G13_G14
        }
      }
      return(inner_e)
    }

    # E(G11^2*mu1^2)
    outer_e11sq<- function() {
      outer_e <- 0
      for (Gi_val in 0:2) {
        # P(G11)
        prob_Gi <-prob_G[Gi_val+1]
        # p(G13,G14|G11)
        joint_con_prob <- joint_con_prob1(Gi_val)
        # E(mu1^2 | G11)
        inner_e <- inner_e0sq(joint_con_prob)
        # E(G11*G11*v1)
        outer_e <- outer_e + (Gi_val-2*p)^2 * prob_Gi * inner_e
      }
      return(outer_e)
    }

    # E(G11*mu1)
    outer_e11<- function() {
      outer_e <- 0
      for (Gi_val in 0:2) {
        # P(G11)
        prob_Gi <-prob_G[Gi_val+1]
        # p(G13,G14|G11)
        joint_con_prob <- joint_con_prob1(Gi_val)
        # E(mu | G11)
        inner_e <- inner_e0(joint_con_prob)
        # E(G11*mu1)
        outer_e <- outer_e + (Gi_val-2*p) * prob_Gi * inner_e
      }
      return(outer_e)
    }

    # E(G11*G12*mu1^2)
    outer_e12<- function() {
      outer_e <- 0
      for (Gi_val in 0:2) {
        for (Gj_val in 0:2) {
          # P(G11, G12)
          prob_Gi_Gj <- genotypes$probrho[genotypes$Gk1==Gi_val & genotypes$Gk2==Gj_val]
          # p(G13,G14|G11,G12)
          joint_con_prob <- joint_con_prob2(Gi_val, Gj_val)
          # E(mu1^2| G11, G12)
          inner_e <- inner_e0sq(joint_con_prob)
          # E(G11*G12*mu1^2)
          outer_e <- outer_e + (Gi_val-2*p) * (Gj_val-2*p)  * prob_Gi_Gj * inner_e
        }
      }

      return(outer_e)
    }

    # E(G11*G13*mu1^2)
    outer_e13<- function() {
      outer_e <- 0
      for (Gi_val in 0:2) {
        # P(G11)
        prob_Gi <-prob_G[Gi_val+1]
        # p(G13,G14|G11)
        joint_con_prob <- joint_con_prob1(Gi_val)
        # E(G13*mu1^2|G11)
        inner_e <- inner_e13(joint_con_prob)
        # E(G11*G13*mu1^2)
        outer_e <- outer_e + (Gi_val-2*p)  * prob_Gi * inner_e
      }
      return(outer_e)
    }

    # E(G11*G14*mu1^2)
    outer_e14<- function() {
      outer_e <- 0
      for (Gi_val in 0:2) {
        # P(G11)
        prob_Gi <-prob_G[Gi_val+1]
        # p(G13,G14|G11)
        joint_con_prob <- joint_con_prob1(Gi_val)
        #  E(G14*mu1^2 | G11)
        inner_e <- inner_e14(joint_con_prob)
        # E(G11*G14*mu1^2)
        outer_e <- outer_e + (Gi_val-2*p)  * prob_Gi * inner_e
      }
      return(outer_e)
    }

   # calculate elements a to j
    a<-outer_e11sq()
    b<-outer_e11()
    c<-sum(genotypes$Gck1^2 * genotypes$mu1^2 * genotypes$probrho)
    d<-sum(genotypes$Gck1 * genotypes$mu1 * genotypes$probrho)
    e<-sum(genotypes$Gck2^2 * genotypes$mu1^2 * genotypes$probrho)
    f<-sum(genotypes$Gck2 * genotypes$mu1 * genotypes$probrho)
    g<-outer_e12()
    h<-outer_e13()
    i<-outer_e14()
    j<-sum(genotypes$Gck1 * genotypes$Gck2^2 * genotypes$mu1^2 * genotypes$probrho)

    # calculate sub-elements aa to gg

    aa<-n*a+n*(n-1)*b^2
    bb<-n*c+n*(n-1)*d^2
    cc<-n*e+n*(n-1)*f^2
    dd<-n*g+n*(n-1)*b^2
    ee<-n*h+n*(n-1)*d*b
    ff<-n*i+n*(n-1)*f*b
    gg<-n*j+n*(n-1)*f*d

    generate_matrix_Y <- function() {
      # Initialize a 4x4 matrix with zeros
      X <- matrix(0, nrow = kk, ncol = kk)

      # Fill in the diagonal elements based on the conditions
      for (k in 1:k2) {
        if (k < k1) {
          X[k, k] <- aa  # For k = 1, 2: X[k, k] = a
        } else if (k == k1) {
          X[k, k] <- bb # For k = 3: X[k, k] = b
        } else if (k == k2) {
          X[k, k] <- cc  # For k = 4: X[k, k] = c
        }
      }

      # Fill in the off-diagonal elements
      for (k in 1:k2) {
        for (l in 1:k2) {
          if (k != l) {
            if (k < k1 && l < k1) {
              X[k, l] <- dd  # For k < 3 and l < 3: X[k, l] = d
            } else if ( (k == k1 & l<k2) || (l == k1 && k<k2)) {
              X[k, l] <- ee  # For k = 3 or l = 3: X[k, l] = e
            } else if (k == k2 || l == k2) {
              if ((k == k2 && l == k1) || (l == k2 && k == k1)) {
                X[k, l] <- gg  # Special case: (k = 4, l = 3) or (k = 3, l = 4): X[k, l] = g
              } else {
                X[k, l] <- ff  # For k = 4 or l = 4: X[k, l] = f
              }
            }
          }
        }
      }
      return(X)
    }

    B <- generate_matrix_Y()/n^2

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
    # under null
    null<-lee(c1)
    q0<-qchisq(1-alpha, df = null$l, ncp = null$delta)
    qc=(q0-null$mux)*null$sigmaq/null$sigmax+null$muq
    # under alternative
    alt<-lee(c_a)
    temp<-(qc-alt$muq)*alt$sigmax/alt$sigmaq+alt$mux
    power[ii]<-pchisq(temp,df=alt$l,ncp=alt$delta,lower.tail = FALSE)
  }
  return(power)
}
