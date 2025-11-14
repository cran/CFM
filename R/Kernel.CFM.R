#' Kernel Censored Factor Model
#'
#' Implementation of kernel-based censored factor model using kernel PCA
#' initialization for nonlinear factor analysis with censored data.
#'
#' @param X Data matrix (n x p)
#' @param m Number of factors
#' @param kernel_type Kernel type: "rbf" or "linear" (default: "rbf")
#' @param gamma Gamma parameter for RBF kernel (default: 1/p)
#' @param max_iter Maximum number of ECM iterations (default: 100)
#' @param tol Convergence tolerance (default: 1e-4)
#' @param nugget Numerical stability term (default: 1e-6)
#' 
#' @return A list containing model results
#' @export
#' @importFrom mvtnorm dmvnorm
#' @importFrom truncnorm rtruncnorm
#' @importFrom stats dist prcomp
censored_kernel_factor_model <- function(X, m, kernel_type = "rbf", 
                                         gamma = NULL, max_iter = 100, 
                                         tol = 1e-4, nugget = 1e-6) {
  
  # Input validation
  if (!is.matrix(X)) stop("X must be a matrix")
  if (m <= 0 || m >= ncol(X)) stop("m must be between 1 and p-1")
  if (is.null(gamma)) gamma <- 1/ncol(X)
  
  # Kernel PCA initialization
  kernel_pca_init <- function(data, m, k_type = "rbf", gam = 1/ncol(data)) {
    X_scaled <- scale(data)
    n <- nrow(X_scaled)
    
    # Compute kernel matrix
    if (k_type == "rbf") {
      # RBF kernel
      dist_sq <- as.matrix(stats::dist(X_scaled))^2
      K <- exp(-gam * dist_sq)
    } else {
      # Linear kernel
      K <- tcrossprod(X_scaled)
    }
    
    # Center kernel matrix
    H <- diag(n) - 1/n
    K_centered <- H %*% K %*% H
    
    # Eigen decomposition
    eig <- eigen(K_centered, symmetric = TRUE)
    lambda <- pmax(eig$values, 0)
    alpha <- eig$vectors[, 1:m, drop = FALSE]
    
    # Recover loadings in feature space
    loadings_approx <- crossprod(X_scaled, alpha) / sqrt(n)
    
    # Ensure reasonable uniqueness values
    h <- rowSums(loadings_approx^2)
    D <- pmax(1 - h, nugget)
    
    list(loadings = loadings_approx, uniqueness = D)
  }
  
  # Initialize parameters
  init <- kernel_pca_init(X, m, kernel_type, gamma)
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

