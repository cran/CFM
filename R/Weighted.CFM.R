#' Weighted Censored Factor Model
#'
#' Implementation of weighted censored factor model with weighted PCA initialization
#' for handling left-censored data with observation weights.
#'
#' @param X Data matrix (n x p)
#' @param m Number of factors
#' @param weights Observation weights vector of length n (optional)
#' @param max_iter Maximum number of ECM iterations (default: 100)
#' @param tol Convergence tolerance (default: 1e-4)
#' @param nugget Numerical stability term (default: 1e-6)
#' 
#' @return A list containing model results
#' @export
#' @importFrom mvtnorm dmvnorm
#' @importFrom truncnorm rtruncnorm
#' @importFrom stats cov2cor prcomp
weighted_censored_factor_model <- function(X, m, weights = NULL, 
                                           max_iter = 100, tol = 1e-4, 
                                           nugget = 1e-6) {
  
  # Input validation
  if (!is.matrix(X)) stop("X must be a matrix")
  if (m <= 0 || m >= ncol(X)) stop("m must be between 1 and p-1")
  if (!is.null(weights) && length(weights) != nrow(X)) {
    stop("weights length must match number of observations")
  }
  
  # Weighted PCA initialization
  weighted_pca_init <- function(data, m, w = NULL) {
    X_scaled <- scale(data)
    n <- nrow(X_scaled)
    p <- ncol(X_scaled)
    
    if (is.null(w)) w <- rep(1, n)
    w <- w / sum(w)  # Normalize weights
    
    # Weighted covariance matrix
    X_centered <- scale(X_scaled, center = TRUE, scale = FALSE)
    Sigma_w <- crossprod(sqrt(w) * X_centered) / (1 - sum(w^2))
    Sigma_w <- stats::cov2cor(Sigma_w)  # Convert to correlation matrix
    
    # Eigen decomposition
    eig <- eigen(Sigma_w, symmetric = TRUE)
    V <- eig$vectors[, 1:m, drop = FALSE]
    h <- rowSums(V^2)
    D <- pmax(1 - h, nugget)
    
    list(loadings = V, uniqueness = D)
  }
  
  # Initialize parameters
  init <- weighted_pca_init(X, m, weights)
  L <- init$loadings
  Psi <- init$uniqueness
  n <- nrow(X)
  p <- ncol(X)
  
  # Initial factor scores
  pca_result <- stats::prcomp(X, center = TRUE, scale. = TRUE)
  F <- scale(pca_result$x[, 1:m, drop = FALSE])
  
  # Censoring detection
  cens <- matrix(FALSE, n, p)
  for (j in 1:p) {
    min_val <- min(X[, j])
    cens[X[, j] <= min_val + sqrt(.Machine$double.eps), j] <- TRUE
  }
  
  loglik <- numeric(max_iter)
  
  # ECM algorithm
  for (it in 1:max_iter) {
    # E-step: Imputation
    Ximp <- X
    for (j in which(colSums(cens) > 0)) {
      idx <- which(cens[, j])
      if (length(idx) > 0) {
        mu <- as.numeric(F[idx, , drop = FALSE] %*% L[j, , drop = TRUE])
        sd_val <- sqrt(Psi[j] + nugget)
        Ximp[idx, j] <- truncnorm::rtruncnorm(
          length(idx), a = X[idx, j], b = Inf, mean = mu, sd = sd_val
        )
      }
    }
    
    # CM-step 1: Update factor scores
    Dinv <- diag(1 / (Psi + nugget))
    M <- crossprod(L, Dinv) %*% L + diag(m)
    M_inv <- solve(M)
    
    for (i in 1:n) {
      rhs <- crossprod(L, Dinv) %*% Ximp[i, ]
      F[i, ] <- as.numeric(M_inv %*% rhs)
    }
    
    # CM-step 2: Update parameters
    FtF <- crossprod(F)
    FtF_inv <- solve(FtF)
    
    for (j in 1:p) {
      L[j, ] <- as.numeric(FtF_inv %*% crossprod(F, Ximp[, j]))
      resid <- Ximp[, j] - F %*% L[j, ]
      Psi[j] <- pmax(mean(resid^2), nugget)
    }
    
    # Log-likelihood
    Sigma <- tcrossprod(L) + diag(Psi)
    tryCatch({
      loglik[it] <- sum(mvtnorm::dmvnorm(Ximp, sigma = Sigma, log = TRUE))
    }, error = function(e) {
      loglik[it] <- -Inf
    })
    
    # Convergence check
    if (it > 1 && abs(loglik[it] - loglik[it-1]) < tol) break
  }
  
  list(
    loadings = L,
    uniqueness = Psi,
    scores = F,
    Sigma = Sigma,
    iterations = it,
    loglik = loglik[it],
    censored = cens
  )
}
