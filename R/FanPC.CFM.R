#' Censored Factor Analysis via Principal Component (FanPC, pure R)
#'
#' @param data  Numeric matrix or data frame of dimension \eqn{n \times p}.
#' @param m  Number of factors (< p).
#' @param A     Optional true loading matrix, used only for error calculation.
#' @param D     Optional true unique-variance diagonal matrix, used only for error calculation.
#' @param p     Number of variables (deprecated; detected automatically).
#' @param cens.dist  Error distribution, reserved for future use.
#' @param df    Degrees of freedom, reserved for future use.
#' @param cens.method  Censoring handling method; currently only \code{"winsorise"} is implemented. 
#'   Defaults to \code{"winsorise"}.
#' @param cens_prop    Winsorisation proportion, default 0.01.
#' @param surv.obj     Reserved for future use.
#' @param ctrl         Reserved for future use.
#' @param verbose      Reserved for future use.
#'
#' @return 
#' \describe{
#'   \item{AF}{Estimated loading matrix, p × m.}
#'   \item{DF}{Estimated unique-variance diagonal matrix, p × p.}
#'   \item{MSESigmaA}{Mean squared error of loadings (if A is provided).}
#'   \item{MSESigmaD}{Mean squared error of unique variances (if D is provided).}
#'   \item{LSigmaA}{Relative error of loadings (if A is provided).}
#'   \item{LSigmaD}{Relative error of unique variances (if D is provided).}
#' }
#'
#' @examples
#' \donttest{
#' library(CFM)
#' obj <- CFM(n = 500, p = 10, m = 2, cens.dist = "normal")
#' res <- FanPC.CFM(obj$data, m = 2, A = obj$A, D = obj$D, cens.method = "winsorise")
#' print(res$MSESigmaA)
#' }
#'
#' @importFrom stats cov quantile
#' @importFrom matrixcalc frobenius.norm
#' @export
FanPC.CFM <- function(data,
                      m,
                      A       = NULL,
                      D       = NULL,
                      p       = NULL,
                      cens.dist = c("normal", "t", "logistic"),
                      df        = NULL,
                      cens.method = c("winsorise", "em"),
                      cens_prop   = 0.01,
                      surv.obj  = NULL,
                      ctrl      = NULL,
                      verbose   = NULL) {
  
  if (!is.matrix(data) && !is.data.frame(data))
    stop("'data' must be a matrix or data frame.")
  X <- as.matrix(data)
  n <- nrow(X);  pp <- ncol(X)
  if (!missing(p) && p != pp)
    warning("'p' is deprecated and will be ignored.")
  m <- as.integer(m)
  if (m <= 0 || m >= pp)
    stop("'m' must be positive and < number of variables.")
  cens.dist   <- match.arg(cens.dist)
  cens.method <- match.arg(cens.method)
  
  if (cens.method != "winsorise")
    stop("Only cens.method = 'winsorise' is implemented in pure-R version.")
  
  
  for (j in 1:pp) {
    y <- X[, j]
    if (all(is.na(y))) next
    qlo <- quantile(y, probs = cens_prop,   na.rm = TRUE)
    qhi <- quantile(y, probs = 1 - cens_prop, na.rm = TRUE)
    y <- pmin(pmax(y, qlo), qhi)
    X[, j] <- y
  }
  
  X <- scale(X, center = TRUE, scale = TRUE)
  Sigma <- cov(X, use = "pairwise.complete.obs")
  eig   <- eigen(Sigma)
  idx   <- order(eig$values, decreasing = TRUE)[1:m]
  AF    <- eig$vectors[, idx] %*% diag(sqrt(pmax(eig$values[idx], 0)), m, m)
  
  h2 <- rowSums(AF^2)
  DF <- diag(Sigma - diag(h2))
  
  if (!is.null(A) && !is.null(D)) {
    MSESigmaA <- frobenius.norm(AF - A)^2 / (pp^2)
    MSESigmaD <- frobenius.norm(DF - D)^2 / (pp^2)
    LSigmaA   <- frobenius.norm(AF - A)^2 / frobenius.norm(A)^2
    LSigmaD   <- frobenius.norm(DF - D)^2 / frobenius.norm(D)^2
  } else {
    MSESigmaA <- MSESigmaD <- LSigmaA <- LSigmaD <- NULL
  }
  
  list(AF        = AF,
       DF        = DF,
       MSESigmaA = MSESigmaA,
       MSESigmaD = MSESigmaD,
       LSigmaA   = LSigmaA,
       LSigmaD   = LSigmaD)
}
