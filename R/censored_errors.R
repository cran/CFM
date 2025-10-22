#' Censored Factor Models Data Generation
#'
#' Generate multivariate data that follow a latent factor structure
#' with censoring errors drawn from Normal, Student-t or Logistic
#' distributions.  Convenience wrapper around \code{\link[crch]{rcnorm}},
#' \code{\link[crch]{rct}}, and \code{\link[crch]{rclogis}}.
#'
#' @param n        sample size (\eqn{n \times 1} observations).
#' @param p        number of manifest variables.
#' @param m        number of latent factors.
#' @param cens.dist censoring error distribution:
#'   \code{"normal"}, \code{"t"}, or \code{"logistic"}.
#' @param df       degrees of freedom when \code{cens.dist = "t"}.
#' @param seed     optional random seed for reproducibility.
#'
#' @return A named list with components:
#'   \item{data}{numeric \eqn{n \times p} matrix of observations.}
#'   \item{F}{factor scores matrix (\eqn{n \times m}).}
#'   \item{A}{factor loadings matrix (\eqn{p \times m}).}
#'   \item{D}{unique variances diagonal matrix (\eqn{p \times p}).}
#'
#' @examples
#' \donttest{
#' set.seed(2025)
#' # Normal censoring
#' obj <- CFM(n = 200, p = 10, m = 3, cens.dist = "normal")
#' head(obj$data)
#'
#' # t-censoring with 6 d.f.
#' obj <- CFM(n = 300, p = 12, m = 4, cens.dist = "t", df = 6)
#' psych::KMO(obj$data)
#' }
#'
#' @export
CFM <- function(n, p, m,
                cens.dist = c("normal", "t", "logistic"),
                df        = 5,
                seed      = NULL) {

  cens.dist <- match.arg(cens.dist)
  if (!is.null(seed)) set.seed(seed)

  ## latent factors
  mu0 <- runif(m, 0, 1)
  Sig <- diag(runif(m, 0.5, 1.5))
  F   <- MASS::mvrnorm(n, mu0, Sig)          # (n x m)

  ## loadings and uniquenesses
  A <- matrix(runif(p * m, -1, 1), p, m)     # (p x m)
  D <- diag(runif(p, 0.2, 1))                # (p x p)

  ## censoring errors
  eps <- switch(cens.dist,
                normal   = crch::rcnorm(n * p, mean = 0, sd = 1),
                t        = crch::rct(n * p, df = df, location = 0, scale = 1),
                logistic = crch::rclogis(n * p, location = 0, scale = 1))
  eps <- matrix(eps, n, p)

  ## observed data
  data <- tcrossprod(F, A) + eps             # (n x p)

  list(data = data, F = F, A = A, D = D)
}
