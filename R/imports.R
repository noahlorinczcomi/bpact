#' @importFrom dplyr filter arrange mutate group_by summarise inner_join select any_of rename
#' @importFrom dplyr "%>%"
#' @importFrom rlang sym
#' @importFrom stats coef dbeta dgamma lm median na.omit optim qgamma quantile rbeta runif toeplitz var
#' @importFrom graphics abline curve hist legend lines par points
#' @importFrom methods as
#' @importClassesFrom Matrix sparseMatrix
#' @importFrom utils tail
NULL

utils::globalVariables(c(
  "chr", "gene", "mid", "shape", "rate",
  "x", "iteration", "pp1", "pp2", "posterior_probability",
  "tau_a", "tau_x1", "tau_x2", "tau_b"
))
