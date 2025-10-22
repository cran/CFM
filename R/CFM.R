#' Basic censored-factor data simulator
#'
#' Generates multivariate data that follow a latent factor structure
#' with censored errors (Normal, Student-t or Logistic).
#'
#' @param n       Sample size (> 0).
#' @param p       Number of observed variables (> 0).
#' @param m       Number of latent factors (< p).
#' @param distribution  Error distribution: "normal" (default), "t", "logistic".
#' @param df      Degrees of freedom when distribution = "t".
#' @param seed    Optional random seed.
#'
#' @return A list with components:
#'   \item{data}{numeric n × p matrix of observations}
#'   \item{loadings}{p × m factor loadings matrix}
#'   \item{uniqueness}{p × p diagonal uniqueness matrix}
#'   \item{KMO}{KMO measure of sampling adequacy}
#'   \item{Bartlett_p}{p-value of Bartlett's test}
#'   \item{distribution}{error distribution used}
#'   \item{seed}{random seed}
#'
#' @examples
#' \donttest{
#' set.seed(2025)
#' obj <- censored_factor_models(200, 6, 2)
#' psych::KMO(obj$data)
#'}
#' @importFrom MASS mvrnorm
#' @importFrom psych KMO cortest.bartlett
#' @importFrom crch rcnorm rct rclogis
#' @export
#'  
censored_factor_models <- function(n, p, m,
                                   distribution = c("normal", "t", "logistic"),
                                   df        = NULL,
                                   seed      = NULL) {
  distribution <- match.arg(distribution)
  if (!is.null(seed)) set.seed(seed)
  
  stopifnot(
    length(n) == 1L, n == as.integer(n), n > 0,
    length(p) == 1L, p == as.integer(p), p > 0,
    length(m) == 1L, m == as.integer(m), 0 < m, m < p
  )
  if (distribution == "t") {
    stopifnot(length(df) == 1L, df == as.integer(df), df > 2)
  }
  
  F <- MASS::mvrnorm(n, mu = rep(0, m), Sigma = diag(m))
  Lambda <- matrix(runif(p * m, -1, 1), p, m)
  Psi    <- diag(runif(p, .2, 1))
  
  eps <- switch(
    distribution,
    normal   = crch::rcnorm(n * p),
    t        = crch::rct(n * p, df = df),
    logistic = crch::rclogis(n * p)
  )
  eps <- matrix(eps, n, p)
  
  X <- F %*% t(Lambda) + eps
  
  kmo  <- psych::KMO(X)
  bart <- psych::cortest.bartlett(R = cor(X), n = n)
  
  list(
    data         = X,
    loadings     = Lambda,
    uniqueness   = Psi,
    KMO          = kmo$KMOS,
    Bartlett_p   = bart$p.value,
    distribution = distribution,
    seed         = seed
  )
}