#' bias_linear.R
#'
#' A function to calculate the bias in misspecified SKAT model for continuous outcomes
#' Solves the score equations for α₀, α_x, and α_G
#' @param beta_list List containing beta_0 (intercept), beta_G (genetic effects), beta_X (covariate effects)
#' @param cov_mat_list List with matrices: mu_GG (cov(G)), mu_GX (cov(G,X)), mu_XX (E[XX^T])
#' @param mu_list List with mu_X (E[X])
#'
#' @return List containing solved coefficients (alpha_0, alpha_G, alpha_X)
#' @export
#'
#' @examples
#' beta_list <- list(beta_0 = 1, beta_X = c(0.7, -0.3), beta_G = c(0.2, 0.4))
#' mu_list <- list(mu_X = c(4, 2))
#' cov_mat_list <- list(
#' mu_XX = matrix(c(5, 3, 3, 4), nrow = 2),
#' mu_X_tildeG = matrix(c(0.1, 0.2, 0.3, 0.6), nrow = 2),
#' mu_XG = matrix(c(0.5, 0.6, 0.7, 0.8), nrow = 2),
#' mu_tildeG_tildeG = matrix(c(1, 0.5, 0.5, 1), nrow = 2),
#' mu_tildeG_G = matrix(c(0.9, 0.8, 0.7, 0.6), nrow = 2)
#' solution <- bias_linear(beta_list, cov_mat_list, mu_list)
#' print(solution)


bias_linear <- function(beta_list, cov_mat_list, mu_list) {
  beta_0 <- beta_list$beta_0
  beta_X <- as.matrix(beta_list$beta_X)
  beta_G <- as.matrix(beta_list$beta_G)

  mu_X <- as.matrix(mu_list$mu_X)
  mu_XX <- cov_mat_list$mu_XX
  mu_X_tildeG <- cov_mat_list$mu_X_tildeG
  mu_XG <- cov_mat_list$mu_XG
  mu_tildeG_tildeG <- cov_mat_list$mu_tildeG_tildeG
  mu_tildeG_G <- cov_mat_list$mu_tildeG_G

  mu_tildeG_X <- t(mu_X_tildeG)

  # Compute O = μ_XX^T - μ_X μ_X^T and its inverse
  O <- mu_XX - (mu_X %*% t(mu_X))
  solve_O <- solve(O)

  # Compute Q and R for α_G
  Q <- mu_tildeG_tildeG - mu_tildeG_X %*% solve_O %*% mu_X_tildeG
  R <- (mu_tildeG_G -  mu_tildeG_X %*% solve_O %*% mu_XG) %*% beta_G
  alpha_G <- solve(Q, R)  # Solves Q α_G = R

  # Compute P and α_x
  P <- (mu_X_tildeG %*% alpha_G) - (mu_XG %*% beta_G)
  alpha_X <- beta_X - solve_O %*% P


  # Compute α_0
  alpha_0 <- beta_0 + (t(mu_X) %*% beta_X) - (t(mu_X) %*% alpha_X)

  return(list(alpha_0 = alpha_0, alpha_X = alpha_X, alpha_G = alpha_G))
}

