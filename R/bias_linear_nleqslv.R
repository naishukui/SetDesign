#' bias_linear_nleqslv.R
#'
#' Solves the score equations for α₀, α_x, and α_G using nleqslv.
#' this could be used to check our direct solution from bias_linear()
#' @param beta_list A list containing beta_0 (intercept), beta_X (covariate effects), beta_G (genetic effects).
#' @param cov_mat_list A list of covariance matrices:
#'                    - mu_XX: E[XXᵀ]
#'                    - mu_X_tildeG: E[XG̃ᵀ]
#'                    - mu_XG: E[XGᵀ]
#'                    - mu_tildeG_tildeG: E[G̃G̃ᵀ]
#'                    - mu_tildeG_G: E[G̃Gᵀ]
#' @param mu_list A list containing mu_X (mean of X).
#'
#' @return A list with solutions for α₀, α_x, and α_G.
#'
#' @export
#' @examples
#' # Example inputs (replace with actual data)
#' beta_list <- list(beta_0 = 1, beta_X = c(0.5, -0.3), beta_G = c(0.2, 0.4))
#' mu_list <- list(mu_X = c(2, 1.5))
#' cov_mat_list <- list(
#'   mu_XX = matrix(c(5, 3, 3, 4), nrow = 2),
#'   mu_X_tildeG = matrix(c(0.1, 0.2, 0.3, 0.4), nrow = 2),
#'   mu_XG = matrix(c(0.5, 0.6, 0.7, 0.8), nrow = 2),
#'   mu_tildeG_tildeG = matrix(c(1, 0.5, 0.5, 1), nrow = 2),
#'   mu_tildeG_G = matrix(c(0.9, 0.8, 0.7, 0.6), nrow = 2)
#' )
#' result <- bias_linear_nleqslv(beta_list, cov_mat_list, mu_list)
#' print(result)

bias_linear_nleqslv <- function(beta_list, cov_mat_list, mu_list) {
  # Load required package
  if (!requireNamespace("nleqslv", quietly = TRUE)) {
    stop("Package 'nleqslv' required but not installed.")
  }

  # Extract parameters and matrices
  beta_0 <- beta_list$beta_0
  beta_X <- as.matrix(beta_list$beta_X)
  beta_G <- as.matrix(beta_list$beta_G)

  mu_X <- as.matrix(mu_list$mu_X)
  mu_XX <- cov_mat_list$mu_XX
  mu_X_tildeG <- cov_mat_list$mu_X_tildeG
  mu_XG <- cov_mat_list$mu_XG
  mu_tildeG_tildeG <- cov_mat_list$mu_tildeG_tildeG
  mu_tildeG_G <- cov_mat_list$mu_tildeG_G

  # Define score equations
  score_eqs <- function(x) {
    alpha_0 <- x[1]
    alpha_X <- as.matrix(x[2:(length(mu_X) + 1)])
    alpha_G <- as.matrix(x[(length(mu_X) + 2):length(x)])

    # Score equations
    y <- numeric(length(x))

    # Equation for α₀: E[Y - α₀ - Xᵀα_x - G̃ᵀα_G] = 0
    y[1] <- alpha_0 + t(mu_X) %*% alpha_X - beta_0 - t(mu_X) %*% beta_X

    # Equation for α_X: E[X(Y - α₀ - Xᵀα_x - G̃ᵀα_G)] = 0
    O <- mu_XX - (mu_X %*% t(mu_X))  # Compute O here
    term_x <- O %*% (alpha_X - beta_X) + (mu_X_tildeG %*% alpha_G - mu_XG %*% beta_G)
    y[2:(length(mu_X) + 1)] <- term_x  # Remove extraneous mu_X term

    # Equation for α_G: E[G̃(Y - α₀ - Xᵀα_x - G̃ᵀα_G)] = 0

    term_G <- mu_tildeG_tildeG %*% alpha_G + t(mu_X_tildeG) %*% alpha_X -
      (mu_tildeG_G %*% beta_G + t(mu_X_tildeG) %*% beta_X)
    y[(length(mu_X) + 2):length(x)] <- term_G

    return(y)
  }

  # Initial guess (α₀, α_x, α_G)
  initial_guess <- c(0, rep(0, length(mu_X)), rep(0, length(beta_G)))

  # Solve numerically
  solution <- nleqslv::nleqslv(
    x = initial_guess,
    fn = score_eqs,
    method = "Newton",
    control = list(ftol = 1e-8)
  )

  # Extract results
  result <- list(
    alpha_0 = solution$x[1],
    alpha_X = solution$x[2:(length(mu_X) + 1)],
    alpha_G = solution$x[(length(mu_X) + 2):length(solution$x)]
  )

  return(result)
}
