#' get_param_mis.R
#' For binary outcome, calculate parameter for misspecified model given parameter of the true model
#' @importFrom rootSolve
#' @param cor correlated genoyypes or not.
#' @param rho Correlation coefficient between SNPs.
#' @param p Minor allele frequency.
#' @param rho Correlation coefficient between SNPs.
#' @param beta1list Vector of effect sizes for the first alternative allele in the multi-allelic position.
#' @param beta2list Vector of effect sizes for the second alternative allele in the multi-allelic position.
#'
#' @return \code{alpha0,alpha1}:Vector of parameters for misspecified model.
#'
#' @export
#' @examples
#' beta1list=c(-0.9, -0.5)
#' beta2list=c(0.9,0.9)
#' get_param_mis(cor=T,p=0.05,rho=0.15,beta1list,beta2list)
#' Should return alpha0=c(-0.9735817, -0.9858034),alpha1=c(0.1056052,0.2539339)


get_param_mis<-function(cor,p,rho,beta1list,beta2list){

  alpha0<-c()
  alpha1<-c()
  beta0 <- -1.0

# Uncentered X in {0,1,2} with pmf:
#   P(X=0)=(1-p)^2,  P(X=1)=2*p*(1-p),  P(X=2)=p^2
X_vals  <- c(0,1,2)
X_probs <- c((1-p)^2, 2*p*(1-p), p^2)

# Centered G = X - 2p
G_vals <- X_vals - 2*p   # numeric vector of length 3
names(G_vals) <- paste0("g=", G_vals)

# Probability for each G-value:
G_probs <- setNames(X_probs, as.character(G_vals))

#when G1, G2 are independent
pairs <- expand.grid(g1 = G_vals, g2 = G_vals, KEEP.OUT.ATTRS=FALSE)
pairs$prob <- mapply(function(x,y) {
  G_probs[as.character(x)] * G_probs[as.character(y)]
}, x = pairs$g1, y = pairs$g2)
if (cor==T){
#when G1, G2 are correlated
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

pairs$probrho<-c(p00,p10,p20,
                 p01,p11,p21,
                 p02,p12,p22)
}
# True two-score logistic model parameters


for (ii in seq_along(beta1list)){
beta1 <-beta1list[ii]
beta2 <-beta2list[ii]

logistic <- function(u) 1 / (1 + exp(-u))
true_mu <- function(g1,g2) {
  logistic(beta0 + beta1*g1 + beta2*g2)
}


residual_fun <- function(par, extra=NULL) {
  alpha0 <- par[1]  # alpha0
  alpha1  <- par[2]  # alpha1

  #  compute E[Y] and E[G * Y], using the TRUE model
  #         E[Y] = sum_{g1,g2} mu(...) * prob(g1,g2)
  #         E[G Y] = sum_{g1,g2} (g1+g2)*mu(...) * prob(g1,g2)
  EY_true    <- 0
  EGY_true   <- 0

  EY_model   <- 0
  EGY_model  <- 0

  for(i in seq_len(nrow(pairs))) {
    g1i  <- pairs$g1[i]
    g2i  <- pairs$g2[i]
    if (cor==F){
    prob <- pairs$prob[i]}
    else {
    prob <- pairs$probrho[i]
    }
    mu_val <- true_mu(g1i, g2i)  # logistic( beta0 + beta1*g1 + beta2*g2 )

    G_val  <- g1i + g2i
    pi_val <- logistic(alpha0 + alpha1*G_val)

    EY_true   <- EY_true    + mu_val * prob
    EGY_true  <- EGY_true   + (g1i+g2i)*mu_val * prob

    EY_model  <- EY_model   + pi_val * prob
    EGY_model <- EGY_model  + (g1i+g2i)*pi_val * prob
  }

  # eqnA = E[Y] - E[ pi(...) ]
  eqnA <- EY_true - EY_model

  # eqnB = E[G * Y] - E[G * pi(...)]
  eqnB <- EGY_true - EGY_model

  return(c(eqnA, eqnB))
}

# Use rootSolve::multiroot to solve eqnA=0, eqnB=0
# initial guess for (alpha0, alpha1)
start_par <- c(0, 0)

sol <- multiroot(f = residual_fun, start = start_par)

alpha0[ii]=sol$root[1]
alpha1[ii]=sol$root[2]

}
return(list(alpha0=alpha0,alpha1=alpha1))
}

