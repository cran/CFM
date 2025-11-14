#' Incremental Censored Factor Model
#'
#' Implementation of incremental censored factor model for streaming data
#' with left-censored observations. Processes data in batches.
#'
#' @param batch_data List of data matrices (batches)
#' @param m Number of factors
#' @param max_iter Maximum number of ECM iterations (default: 100)
#' @param tol Convergence tolerance (default: 1e-4)
#' @param nugget Numerical stability term (default: 1e-6)
#' 
#' @return A list containing model results
#' @export
#' @importFrom mvtnorm dmvnorm
#' @importFrom truncnorm rtruncnorm
#' @importFrom stats cor prcomp
incremental_censored_factor_model <- function(batch_data, m, 
                                              max_iter = 100, tol = 1e-4, 
                                              nugget = 1e-6) {
  
  # Input validation
  if (!is.list(batch_data)) stop("batch_data must be a list of matrices")
  if (length(batch_data) == 0) stop("batch_data cannot be empty")
  if (m <= 0 || m >= ncol(batch_data[[1]])) {
    stop("m must be between 1 and p-1")
  }
  
  # Incremental PCA initialization
  incremental_pca_init <- function(batches, m) {
    # Process first batch
    X1 <- scale(batches[[1]])
    p <- ncol(X1)
    Sigma <- stats::cor(X1)
    
    # Initialize with first batch
    eig <- eigen(Sigma, symmetric = TRUE)
    V <- eig$vectors[, 1:m, drop = FALSE]
    n_total <- nrow(X1)
    
    # Process subsequent batches incrementally
    if (length(batches) > 1) {
      for (i in 2:length(batches)) {
        X_batch <- scale(batches[[i]])
        n_batch <- nrow(X_batch)
        
        # Update covariance matrix incrementally
        Sigma_batch <- stats::cor(X_batch)
        Sigma <- (n_total * Sigma + n_batch * Sigma_batch) / (n_total + n_batch)
        
        # Update eigenvectors
        eig <- eigen(Sigma, symmetric = TRUE)
        V <- eig$vectors[, 1:m, drop = FALSE]
        n_total <- n_total + n_batch
      }
    }
    
    # Compute uniqueness
    h <- rowSums(V^2)
    D <- pmax(1 - h, nugget)
    
    list(loadings = V, uniqueness = D)
  }
  
  # Initialize parameters
  init <- incremental_pca_init(batch_data, m)
  L <- init$loadings
  Psi <- init$uniqueness
  
  # Combine all batches
  X_all <- do.call(rbind, batch_data)
  X_all <- scale(X_all)
  n <- nrow(X_all)
  p <- ncol(X_all)
  
  # Initial factor scores
  pca_result <- stats::prcomp(X_all, center = TRUE, scale. = TRUE)
  F <- scale(pca_result$x[, 1:m, drop = FALSE])
  
  # Censoring detection
  cens <- matrix(FALSE, n, p)
  for (j in 1:p) {
    min_val <- min(X_all[, j])
    cens[X_all[, j] <= min_val + sqrt(.Machine$double.eps), j] <- TRUE
  }
  
  loglik <- numeric(max_iter)
  
  # ECM algorithm
  for (it in 1:max_iter) {
    # E-step: Imputation
    Ximp <- X_all
    for (j in which(colSums(cens) > 0)) {
      idx <- which(cens[, j])
      if (length(idx) > 0) {
        mu <- as.numeric(F[idx, , drop = FALSE] %*% L[j, , drop = TRUE])
        sd_val <- sqrt(Psi[j] + nugget)
        Ximp[idx, j] <- truncnorm::rtruncnorm(
          length(idx), a = X_all[idx, j], b = Inf, mean = mu, sd = sd_val
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
    censored = cens,
    n_batches = length(batch_data)
  )
}
