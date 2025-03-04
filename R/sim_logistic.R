#' sim_logistic.R
#' this function simulates the above setting so that we can check to make the
#' solver code solve_logistic() is working correctly.
#' @param betaVec is a 3*1 vector of the true coefficients (intercept, minor allele 1, minor allele 2)
#' @param probVec is a 3*1 vector of the probabilities of each of the three alleles (major allele, minor 1, minor 2)
#' @param nSubs is the number of subjects to use in simulation
#' @param nSims is the number of simulations to perform
#' @param checkpoint is just a progress update
#' @return List containing solved coefficients (alpha0, alphaG)
#' @export

sim_logistic <- function(betaVec, nSubs=5000, probVec=c(0.98, 0.01, 0.01),
                         nSims=5000, checkpoint=TRUE) {

  # hold the fitted coefficients in this matrix
  fittedCoef <- matrix(data=NA, nrow=nSims, ncol=2)
  for (sim_it in 1:nSims) {
    # simulate the tri-allelic genotype
    gVec <- rmultinom(n=nSubs, size=2, prob=probVec)
    # add the two minor alleles
    gMinor <- as.numeric(gVec[2, ] + gVec[3, ])
    # true model
    trueEta <- cbind(1, t(gVec[2:3, ])) %*% betaVec
    trueMu <- rje::expit(trueEta)
    yVec <- rbinom(n=nSubs, size=1, prob=trueMu)
    # fit the misspecified model
    fittedMod <- glm(yVec ~ gMinor, family=binomial)
    # record
    fittedCoef[sim_it, ] <- summary(fittedMod)$coefficients[, 1]

    if(checkpoint) {
      if (sim_it%%100 == 0) {cat(sim_it)}
    }
  }

  return(fittedCoef)
}

# # test the simulation
# set.seed(0)
# test = sim_logistic(betaVec = c(-1, -0.5, 1), probVec=c(0.9, 0.08, 0.02))
# apply(test, 2, mean)
# # test the solver
# test2 <- solve_logistic(betaVec = c(-1, -0.5, 1), probVec=c(0.9, 0.08, 0.02))
# test2$x
