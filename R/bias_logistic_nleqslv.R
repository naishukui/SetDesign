#' bias_logistic_nleqslv.R
#'
#' A function to solve score equations for bias in misspecified logistic regression models
#' This code considers the case where the true model has one tri-allelic variant
#' (and thus there are two true covariates in the model), however we fit a misspecified
#' model where the two minior alleles are summed together and thus there is only one covariate in
#' the fitted model.
#' Solves the score equations for α₀ and α_G
#' @param betaVec is a 3*1 vector of the true coefficients (intercept, minor allele 1, minor allele 2)
#' @param probVec is a 3*1 vector of the probabilities of each of the three alleles (major allele, minor 1, minor 2)
#'
#' @return List containing solved coefficients (alpha0, alphaG)
#' @export


bias_logistic_nleqslv <- function(betaVec, probVec) {

  # here we define the score equation to be solved. x is a 2*1 vector of initial guesses
  # for the misspecified model coefficients
  scoreEq <- function(x) {
    # the two misspecified coefficients
    alpha0 <- x[1]
    alphaG <- x[2]
    # the value of the score equations
    scoreVec <- c(0, 0)
    # here we expand the expectations of the score equations. there are 3*3=9 possible combinations
    # of the genotypes of the tri-alleleic variant.
    scoreVec[1] <- probVec[1]^2 * (rje::expit(betaVec[1]) - rje::expit(alpha0)) +
      probVec[1] * probVec[2] * (rje::expit(betaVec[1] + betaVec[2]) - rje::expit(alpha0 + alphaG)) +
      probVec[1] * probVec[3] * (rje::expit(betaVec[1] + betaVec[3]) - rje::expit(alpha0 + alphaG)) +
      probVec[2] * probVec[1] * (rje::expit(betaVec[1] + betaVec[2]) - rje::expit(alpha0 + alphaG)) +
      probVec[2] * probVec[2] * (rje::expit(betaVec[1] + 2*betaVec[2]) - rje::expit(alpha0 + 2*alphaG)) +
      probVec[2] * probVec[3] * (rje::expit(betaVec[1] + betaVec[2] + betaVec[3]) - rje::expit(alpha0 + 2*alphaG)) +
      probVec[3] * probVec[1] * (rje::expit(betaVec[1] + betaVec[3]) - rje::expit(alpha0 + alphaG)) +
      probVec[3] * probVec[2] * (rje::expit(betaVec[1] + betaVec[2] + betaVec[3]) - rje::expit(alpha0 + 2*alphaG)) +
      probVec[3] * probVec[3] * (rje::expit(betaVec[1] + 2*betaVec[3]) - rje::expit(alpha0 + 2*alphaG))

    scoreVec[2] <- probVec[1] * probVec[2] * (rje::expit(betaVec[1] + betaVec[2]) - rje::expit(alpha0 + alphaG)) +
      probVec[1] * probVec[3] * (rje::expit(betaVec[1] + betaVec[3]) - rje::expit(alpha0 + alphaG)) +
      probVec[2] * probVec[1] * (rje::expit(betaVec[1] + betaVec[2]) - rje::expit(alpha0 + alphaG)) +
      2 * probVec[2] * probVec[2] * (rje::expit(betaVec[1] + 2*betaVec[2]) - rje::expit(alpha0 + 2*alphaG)) +
      2 * probVec[2] * probVec[3] * (rje::expit(betaVec[1] + betaVec[2] + betaVec[3]) - rje::expit(alpha0 + 2*alphaG)) +
      probVec[3] * probVec[1] * (rje::expit(betaVec[1] + betaVec[3]) - rje::expit(alpha0 + alphaG)) +
      2 * probVec[3] * probVec[2] * (rje::expit(betaVec[1] + betaVec[2] + betaVec[3]) - rje::expit(alpha0 + 2*alphaG)) +
      2 * probVec[3] * probVec[3] * (rje::expit(betaVec[1] + 2*betaVec[3]) - rje::expit(alpha0 + 2*alphaG))

    scoreVec

  }

  # this part of the code just runs the nleqslv() function
  solveEqs <- nleqslv::nleqslv(x=c(0, 0), fn=scoreEq)
  solveEqs
}

