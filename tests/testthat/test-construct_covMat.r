# library(testthat)

## Helper functions --------------------------------------------------------

is_psd <- function(M, tol = 1e-8) {
  ev <- eigen(M, symmetric = TRUE, only.values = TRUE)$values
  min(ev) > -tol
}

is_toeplitz <- function(M, tol = 1e-8) {
  all(abs(M[row(M) - col(M)] - M[1, 1]) < tol) # simplified check
}

is_banded <- function(M, bw, tol = 1e-8) {
  idx <- abs(row(M) - col(M)) > bw
  all(abs(M[idx]) < tol)
}

is_block <- function(M, groups, tol = 1e-12) {
  d <- sum(groups)
  if (!isSymmetric(M, tol = tol)) return(FALSE)
  if (nrow(M) != d) return(FALSE)
  block_id <- rep(seq_along(groups), groups)

  # Check all pairs (i,j)
  for (i in 1:d) {
    for (j in 1:d) {
      same_block <- (block_id[i] == block_id[j])

      if (same_block) {
        if (i != j) {
          b <- block_id[i]
          idx_block <- which(block_id == b)
          ref <- M[idx_block[1], idx_block[min(2, length(idx_block))]]
          if (abs(M[i, j] - ref) > tol) return(FALSE)
        }
      } else {
        if (abs(M[i, j]) > tol) return(FALSE)
      }
    }
  }
  TRUE
}

is_sparse <- function(M, sparsity) {
  offdiag <- M[row(M) != col(M)]
  prop_zero <- mean(abs(offdiag) < 1e-10)
  prop_zero >= sparsity
}

get_diag <- function(M, k = 0) {
  n <- nrow(M)
  if (k >= n) return(numeric(0))
  i <- 1:(n - k)
  j <- i + k
  M[cbind(i, j)]
}

check_acf <- function(M, phi = numeric(0), theta = numeric(0), sigma2 = 1, corr = TRUE) {
  n <- nrow(M)
  acf_true <- stats::ARMAacf(ar = phi, ma = theta, lag.max = n-1)
  if (!corr) acf_true <- acf_true * sigma2
  expect_equal(M[1, ], unname(acf_true), tolerance = 1e-6)
}

## Tests for Basic structures, including block --------------------------------------------------------


test_that("diag stays diagonal even with random variances", {
  v <- runif(10, 0.1, 2)
  M <- construct_covMat(10, corr=FALSE, type="diag", params=list(variances=v))

  expect_equal(M, diag(v))
  expect_true(isSymmetric(M,tol=0))
})


test_that("CS eigenvalues match analytic form", {
  d  <- 8
  rho <- 0.3
  M <- construct_covMat(d, corr=T, type="CS", params=list(rho=rho, variances=1:d))

  expect_true(isSymmetric(M,tol=0))
  expect_equal(diag(M), rep(1,d))
  expect_true(all(M[lower.tri(M)] == rho))

  ev <- sort(eigen(M)$values)
  expect_equal(ev[1], 1 - rho)
  expect_equal(ev[d], 1 + (d-1)*rho)
})



test_that("banded zeros outside bandwidth", {
  bw <- 2
  d<- 12
  M  <- construct_covMat(12, corr=TRUE, type="CS",
                            params=list(bandwidth=bw, rho=0.3))

  idx <- abs(row(M)-col(M)) > bw
  expect_equal(M[idx], rep(0, sum(idx)))

  expect_equal(diag(M), rep(1, d))
  expect_true(is_banded(M, bw=bw))
})



test_that("block off-diagonal blocks are zero", {
  groups <- c(3,4)
  M <- construct_covMat(7, corr=TRUE, type="block",
                           params=list(groups=groups, rho=c(.6,.2)))
  G1 <- 1:3
  G2 <- 4:7

  expect_true(is_block(M, groups=groups))
  expect_equal(M[G1,G2], matrix(0,3,4))
})



## Tests for Time series cov --------------------------------------------------------


test_that("toeplitz exact subdiagonals", {

  d<-4
  rho <- c(1, 0.8, 0.5, 0.1)  # first-row rho, including diagonal
  M <- construct_covMat(d, corr = TRUE, type = "toeplitz", params = list(rho = rho))

  expect_equal(diag(M), rep(1, d))
  expect_true(isSymmetric(M, tol=0))

  # Test Toeplitz structure
  expect_true(all(M[row(M) - col(M) == 1] == rho[2]))


  get_diag <- function(M, k = 0) {
    n <- nrow(M)
    if (k >= n) return(numeric(0))
    i <- 1:(n - k)
    j <- i + k
    M[cbind(i, j)]
  }

  for (k in 0:3) {
    expect_true(all(abs(get_diag(M, k) - rho[k + 1]) ==0))
  }
})

test_that("AR1 satisfies the Markov identity", {
  rho <- 0.8
  M <- construct_covMat(6, corr=TRUE, type="AR", params=list(rho=rho))
  expect_equal(M[1,2], rho)
  expect_equal(M[1,3], rho^2)
  expect_true(is_psd(M, tol=0))
  expect_equal(M[1,3], M[1,2] * M[2,3])
})

test_that("ARMA matches theoretical ACF", {
  phi <- c(.5, -.3)
  theta <- 0
  d <- 10
  M <- construct_covMat(d, corr=TRUE, type="ARMA",
                           params=list(arma_phi=phi, arma_theta=theta))

  check_acf(M, phi, theta, corr=TRUE)

  expect_true(isSymmetric(M, tol=0))
  expect_true(is_psd(M,tol=0))
  expect_true(all(eigen(M, only.values=TRUE)$values >= 0))
  expect_equal(diag(M), rep(1, d))
})

test_that("Pure AR covariance", {
  phi <- c(0.7)
  M <- construct_covMat(6, params=list(arma_phi=phi, sigma2=1), type="AR")

  check_acf(M, phi, theta=numeric(0), sigma2=1, corr=TRUE)

  expect_true(isSymmetric(M, tol=0))
  expect_true(all(eigen(M, only.values=TRUE)$values >= 0))
})

test_that("Pure MA covariance", {
  theta <- c(0.5, 0.2)
  M <- construct_covMat(5, params=list(arma_theta=theta), type="MA")

  check_acf(M, phi=numeric(0), theta=theta,  corr=TRUE)

  expect_true(isSymmetric(M, tol=0))
  expect_true(all(eigen(M, only.values=TRUE)$values >= 0))
})

test_that("Correlation matrix output", {
  phi <- c(0.5)
  theta <- c(0.3)
  sigma2 <- 2
  M <- construct_covMat(5, params=list(arma_phi=phi, arma_theta=theta, sigma2=sigma2), type="ARMA", corr=TRUE)

  check_acf(M, phi, theta, corr=TRUE)
  expect_true(isSymmetric(M, tol=0))
  expect_equal(diag(M), rep(1,5))
})

test_that("Zero AR/MA orders", {
  M <- construct_covMat(4, params=list(arma_phi=numeric(0), arma_theta=numeric(0), sigma2=1), type="ARMA")

  expect_true(isSymmetric(M))
  expect_equal(diag(M), rep(1, 4))
  expect_true(all(abs(M[upper.tri(M)]) == 0))
})

test_that("ARIMA covariance", {
  phi <- c(0.3)
  theta <- c(0.2)
  d <- 1
  dim <- 6

  # Construct covariance matrix
  M <- construct_covMat(dim, params=list(
    arma_phi=phi, arma_theta=theta, d=d, sigma2=1
  ), type="ARIMA")

  # --- Build AR/MA polynomials ---
  ar_poly <- c(1, -phi)
  ma_poly <- c(1, theta)
  ar_combined <- -ar_poly[-1]
  ma_combined <- ma_poly[-1]

  # --- Expand to MA(∞) for differencing ---
  ma_length <- dim + 50
  psi <- stats::ARMAtoMA(ar = ar_combined, ma = ma_combined, lag.max = ma_length)

  # Non-seasonal differencing (1-B)^d
  base_diff <- c(1, -1)
  difff <- base_diff
  if (d > 1) {
    for (k in 2:d) difff <- convolve(difff, base_diff, type="open")
  }
  psi <- convolve(psi, difff, type="open")

  # Compute autocovariances
  gamma <- function(h) {
    if (h < length(psi))
      sum(psi[1:(length(psi)-h)] * psi[(1+h):length(psi)])
    else 0
  }
  acvs <- sapply(0:(dim-1), gamma)

  expect_equal(M[1, 1:length(acvs)], unname(acvs))
  expect_true(isSymmetric(M, tol=0))
  expect_true(all(eigen(M, only.values=TRUE)$values >= 0))
})




## Tests for Factor --------------------------------------------------------


test_that("factor model produces valid correlation when corr=TRUE", {
  L <- matrix(rnorm(15), 5, 3)
  psi <- runif(5, .2, .5)

  M <- construct_covMat(5, type="factor",
                           params=list(loadings=L, psi=psi))

  # Covariance = L L' + diag(psi)
  Cov <- L %*% t(L) + diag(psi)
  expect_equal(norm(M - Cov, "F"), 0)

  ## For correlation
  M <- construct_covMat(5, corr=TRUE, type="factor",
                           params=list(loadings=L, psi=psi))

  # recompute covariance then convert to correlation
  Corr <- Cov / sqrt(outer(diag(Cov), diag(Cov)))
  expect_equal(M, Corr)
})


## Tests for Graph  --------------------------------------------------------


test_that("graph", {
  A <- matrix(0,6,6)
  A[cbind(1:5, 2:6)] <- 1
  A[cbind(2:6, 1:5)] <- 1

  M <- construct_covMat(6, corr=FALSE, type="graph",
                           params=list(adjacency=A))
  expect_true(is_psd(M, tol=0))
  expect_true(isSymmetric(M,tol=0))

  M <- construct_covMat(6, corr=TRUE, type="graph",
                           params=list(adjacency=A))
  expect_true(is_psd(M, tol=0))
  expect_true(isSymmetric(M,tol=0))
})





## Tests for Spatial cov --------------------------------------------------------


test_that("spatial exponential", {
  dim<-6
  D <- as.matrix(dist(seq(0,1,length.out=6)))
  vars<-1:dim # rep(2,dim)

  M <- construct_covMat(6, corr=TRUE, type="spatial",
                           params=list(distances=D, kernel="exp", phi=3))
  expect_equal(diag(M), rep(1,6))
  expect_true(is_psd(M,0))
  expect_true(isSymmetric(M,0))

  M <- construct_covMat(6, corr=FALSE, type="spatial",
                           params=list(distances=D, kernel="exp", phi=3, variances=vars))
  expect_equal(diag(M), vars)
  expect_true(is_psd(M,0))
  expect_true(isSymmetric(M,0))
})

test_that("spatial matern", {
  D <- as.matrix(dist(seq(0,1,length.out=6)))
  M <- construct_covMat(6, corr=TRUE, type="spatial",
                           params=list(distances=D, kernel="Matern", phi=1, nu=0.5))
  expect_equal(diag(M), rep(1,6))
  expect_true(is_psd(M,0))
})

test_that("spatial Gaussian kernel decreases with distance", {
  D <- as.matrix(dist(seq(0,1,length=10)))
  M <- construct_covMat(10, corr=TRUE, type="spatial",
                           params=list(distances=D, kernel="norm", phi=.3))

  expect_true(M[1,2] > M[1,5])   # farther distances → lower corr
})


test_that("spatial custom kernel", {
  D <- as.matrix(dist(seq(0,1,length.out=6)))
  M <- construct_covMat(6, corr=TRUE, type="spatial",
                           params=list(distances=D, kernel=function(d) exp(-d^2)))
  expect_equal(diag(M), rep(1,6))
  expect_true(is_psd(M,0))
})


test_that("custom monotone kernel decreases", {
  ker <- function(d) exp(-d^2)
  D <- as.matrix(dist(seq(0,1,len=10)))
  M <- construct_covMat(10, corr=TRUE, type="spatial",
                           params=list(distances=D, kernel=ker))

  expect_true(M[1,2] > M[1,5])
})



## Tests for Constrained cov --------------------------------------------------------

test_that("rank projection yields unit diagonal", {
  rank <- 6
  d<- 10
  M <- suppressWarnings(construct_covMat(d, corr=TRUE, type="specConstr", params=list(rank=rank)))

  ev <- eigen(M, symmetric=TRUE)$values
  expect_equal(sum(ev > 1e-12), rank)
  expect_true(all(ev >= -1e-12))
  expect_true(isSymmetric(M,tol=0))
  expect_equal(diag(M), rep(1, d))

  M <- suppressWarnings(construct_covMat(d, corr=FALSE, type="specConstr", params=list(rank=rank, variances=1:d)))
  ev <- eigen(M, symmetric=TRUE)$values
  expect_equal(sum(ev > 1e-12), rank)
  expect_true(all(ev >= -1e-12))
  expect_true(isSymmetric(M,tol=0))
  expect_equal(diag(M), 1:d)
})

test_that("CN eigenvalues preserve ordering", {
  CN<- 10
  dim<- 8
  vars<-1:dim # rep(2,dim)

  M <- construct_covMat(dim, corr=TRUE, type="specConstr", params=list(cond_num=CN))
  ev_cn <- eigen(M, symmetric=TRUE)$values
  expect_true(isSymmetric(M,tol=0))
  expect_equal(diag(M), rep(1, dim))
  expect_true(all(ev_cn >= 0))

  nonzero_ev_cn <- sort(ev_cn[ev_cn > 0], decreasing = TRUE)
  actual_cond_cn <- nonzero_ev_cn[1] / nonzero_ev_cn[length(nonzero_ev_cn)]
  expect_equal(actual_cond_cn, CN)

  M <- construct_covMat(dim, corr=FALSE, type="specConstr", params=list(cond_num=CN, variances=vars))
  ev_cn <- eigen(M, symmetric=TRUE)$values
  expect_true(isSymmetric(M,tol=0))
  expect_equal(diag(M), vars)
  expect_true(all(ev_cn >= 0))

  nonzero_ev_cn <- sort(ev_cn[ev_cn > 0], decreasing = TRUE)
  actual_cond_cn <- nonzero_ev_cn[1] / nonzero_ev_cn[length(nonzero_ev_cn)]
  if (length(unique(vars))==1) expect_equal(actual_cond_cn, CN) else expect_lte(actual_cond_cn, CN)
})

test_that("both CN and rank works properly", {
  # set.seed(42)

  dim <- 10
  rank_target <- 5
  cond_num <- 1000
  vars<-1:dim # rep(2,dim)
  eps<-1e-12

  # Call for correlation
  mat <- suppressWarnings(construct_covMat(dim = dim, corr = TRUE, type = "specConstr",
                             params = list(rank = rank_target, cond_num = cond_num, variances=vars)))

  # 1. Symmetric
  expect_true(isSymmetric(mat, tol = 0))
  # 2. Diagonal = 1 (correlation matrix)
  expect_equal(diag(mat), rep(1, dim))
  # 3. PSD
  ev <- eigen(mat, symmetric = TRUE)$values
  expect_true(all(ev >= -eps))
  # 4. Rank
  effective_rank <- sum(ev > eps)
  expect_equal(effective_rank, rank_target)
  # 5. Condition number of non-zero eigenvalues
  nonzero_ev <- sort(ev[ev > eps], decreasing = TRUE)
  actual_cond <- nonzero_ev[1] / nonzero_ev[length(nonzero_ev)]
  expect_equal(actual_cond, cond_num)  # small numerical tolerance

  # Call for covarince
  mat <- suppressWarnings(construct_covMat(dim = dim, corr = FALSE, type = "specConstr",
                             params = list(rank = rank_target, cond_num = cond_num, variances=vars)))

  # 1. Symmetric
  expect_true(isSymmetric(mat, tol = 0))
  # 2. Diagonal = 1 (correlation matrix)
  expect_equal(diag(mat), vars)
  # 3. PSD
  ev <- eigen(mat, symmetric = TRUE)$values
  expect_true(all(ev >= -eps))
  # 4. Rank
  effective_rank <- sum(ev > eps)
  expect_equal(effective_rank, rank_target)
  # 5. Condition number of non-zero eigenvalues
  nonzero_ev <- sort(ev[ev > eps], decreasing = TRUE)
  actual_cond <- nonzero_ev[1] / nonzero_ev[length(nonzero_ev)]
  if (length(unique(vars))==1) expect_equal(actual_cond, cond_num) else expect_lte(actual_cond, cond_num)
})

test_that("eigen works correctly", {

  dim <- 6
  eigen_target <- c(2, 1.5, 1, 0.5, 0.7, 0.3)  # will be padded to sum to dim=6
  vars<-1:dim # rep(2,dim)
  eps<- 1e-8

  R <- construct_covMat(dim, corr=TRUE, type="specConstr",
                           params=list(eigenvalues=eigen_target))
  # 1. Symmetry
  expect_true(isSymmetric(R, tol = 0))
  # 2. Diagonal = 1
  expect_true(max(abs(diag(R) - 1)) < eps)
  # 3. Eigenvalues match target (up to sorting + tolerance)
  eig <- eigen(R, symmetric = TRUE)$values
  target_full <- c(eigen_target, rep(0, dim - length(eigen_target)))
  expect_true(max(abs(sort(eig, decreasing = TRUE) -
                        sort(target_full, decreasing = TRUE))) < eps)
  # 4. PSD
  expect_true(all(eig > -eps))
  # 5. Correlation matrix sanity check
  expect_true(all(R <= 1 + eps))
  expect_true(all(R >= -1 - eps))


  eigen_target <- c(2, 2.5, 1, 0.5)*sum(vars)/dim  # will be padded to sum to dim=6
  R <- suppressWarnings(construct_covMat(dim, corr=FALSE, type="specConstr",
                           params=list(eigenvalues=eigen_target, variances=vars)))
  # 1. Symmetry
  expect_true(isSymmetric(R, tol = 0))
  # 2. Diagonal = 1
  expect_true(max(abs(diag(R) - vars)) < eps)
  # 3. Eigenvalues match target (up to sorting + tolerance)
  eig <- eigen(R, symmetric = TRUE)$values
  target_full <- c(eigen_target, rep(0, dim - length(eigen_target)))
  expect_true(max(abs(sort(eig, decreasing = TRUE) -
                        sort(target_full, decreasing = TRUE))) < eps)
  # 4. PSD
  expect_true(all(eig > -eps))
})


