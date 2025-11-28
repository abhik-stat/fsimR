#' @title Construct Structured Covariance or Correlation Matrices (non-random)
#'
#' @description
#' Constructs a \eqn{dim \times dim} covariance or correlation matrix according to a specified structural type.
#' Supported types include diagonal, compound symmetry, block, Toeplitz, autoregressive/moving-average processes
#' (AR/MA/ARMA/ARIMA), factor models, graph-based precision structures, spatial kernels, and spectrum-constrained matrices.
#'
#' @param dim Integer. Dimension of the resulting matrix.
#' @param corr Logical. If `TRUE`, returns a correlation matrix.
#'   If `FALSE` (default), returns a covariance matrix scaled by `variances`.
#'   For types `factor` and `graph`, `variances` are ignored.
#'
#' @param type Character. Structure type of the covariance/correlation matrix. Supported values:
#'   \describe{
#'     \item{\code{"diag"}}{Diagonal matrix (identity correlation).}
#'     \item{\code{"CS"}}{Compound symmetry (exchangeable correlation).}
#'     \item{\code{"block"}}{Block-diagonal matrix.}
#'     \item{\code{"toeplitz"}}{Toeplitz matrix with lag correlations (AR(1) structure).}
#'     \item{\code{"AR"}}{Autoregressive AR(p) process.}
#'     \item{\code{"MA"}}{Moving average MA(q) process.}
#'     \item{\code{"ARMA"}}{ARMA(p, q) process.}
#'     \item{\code{"ARIMA"}}{ARIMA(p, d, q) process.}
#'     \item{\code{"SAR"}}{Seasonal AR process.}
#'     \item{\code{"SMA"}}{Seasonal MA process.}
#'     \item{\code{"SARMA"}}{Seasonal ARMA process.}
#'     \item{\code{"SARIMA"}}{Seasonal ARIMA process.}
#'     \item{\code{"factor"}}{Factor model: loadings + uniquenesses.}
#'     \item{\code{"graph"}}{Graph-structured precision matrix (GMRF-style).
#'                Uses precision = adjacency + diagonal dominance; covariance = inverse of precision.}
#'     \item{\code{"spatial"}}{Spatial correlation via a distance matrix and kernel function.}
#'     \item{\code{"specConstr"}}{Spectrum-constrained matrix with specified eigenvalues, rank or condition number.}
#'   }
#' @param params Named list of parameters required for the selected type. Common parameters:
#'   \describe{
#'     \item{`variances`}{Numeric scalar or vector; used if `corr = FALSE` and `type` is neither `factor` nor `graph`.
#'     Defaults to 1.}
#'     \item{`rho`}{Numeric scalar, vector, or matrix/list. Correlation parameter(s)}
#'     \item{`groups`}{Vector of group sizes for block-diagonal matrices.}
#'     \item{`arma_phi`, `arma_theta`}{AR and MA coefficients for time series types.}
#'     \item{`d`}{Differencing order for ARIMA/SARIMA processes.}
#'     \item{`D`, `seasonal_period`}{Seasonal differencing order and period for SARIMA/SARMA/SMA/SAR processes.}
#'     \item{`seasonal_phi`, `seasonal_theta`}{Seasonal AR/MA coefficients.}
#'     \item{`loadings`, `psi`}{Factor loadings and uniquenesses for factor models.}
#'     \item{`adjacency`}{Symmetric adjacency matrix for graph-based covariance.}
#'     \item{`distances`}{Distance matrix for spatial correlation. Defaults to absolute differences.}
#'     \item{`kernel`}{Kernel function for spatial correlation.
#'     Can be a function or one of: `"exponential"`, `"gaussian"`, `"matern"`, `"rational_quadratic"`,
#'     `"powered_exponential"`, `"cauchy"`.}
#'     \item{`phi`, `nu`, `alpha`, `beta`, `power`}{Spatial kernel parameters.}
#'     \item{`eigenvalues`, `rank`, `cond_num`}{Parameters for spectrum-constrained matrices.
#'            Specified eigenvalues, rank and condition numbers respectively.}
#'     \item{bandwidth}{Integer; optional for CS or Toeplitz matrices to zero out entries beyond a certain lag.}
#'     \item{max_iter, tol}{Maximum iterations and tolerance for spectrum-constrained iterative projections. See Details.}
#'   }
#' @param eigen_check Logical. If `TRUE` (default), checks if the resulting matrix is PSD/PD.
#'
#' @return A symmetric `dim` \eqn{\times} `dim` matrix of class `covMat`. Either correlation or covariance depending on `corr`.
#'
#' @details
#' The `params` list must contain type-specific elements as follows:
#' \describe{
#'   \item{"diag"}{Optional: `variances` (numeric vector or scalar) if `corr = FALSE`. Defaults to 1.}
#'   \item{"CS"}{`rho` (scalar correlation), optional `variances` (if `corr = FALSE`), optional `bandwidth`.}
#'   \item{"block"}{`groups` (vector of integers specifying block sizes), `rho` (scalar or vector for each block), optional `variances`.}
#'   \item{"toeplitz"}{`rho` (vector or scalar), optional `variances`, optional `bandwidth`.}
#'   \item{"AR"}{`arma_phi` (AR coefficients), optional `variances`.}
#'   \item{"MA"}{`arma_theta` (MA coefficients), optional `variances`.}
#'   \item{"ARMA"}{`arma_phi` (AR coefficients), `arma_theta` (MA coefficients), optional `variances`.}
#'   \item{"ARIMA"}{`arma_phi` (AR), `arma_theta` (MA), `d` (differencing order), optional `variances`.}
#'   \item{"SAR"}{`seasonal_phi` (seasonal AR), `D` (seasonal differencing), `seasonal_period`.}
#'   \item{"SMA"}{`seasonal_theta` (seasonal MA), `D`, `seasonal_period`.}
#'   \item{"SARMA"}{`seasonal_phi`, `seasonal_theta`, `D`, `seasonal_period`.}
#'   \item{"SARIMA"}{`arma_phi`, `arma_theta`, `d`, `seasonal_phi`, `seasonal_theta`, `D`, `seasonal_period`.}
#'   \item{"factor"}{`loadings` (matrix `dim` \eqn{\times} `k`), `psi` (vector of uniquenesses length dim).}
#'   \item{"graph"}{`adjacency` (symmetric matrix `dim` \eqn{\times} `dim`).}
#'   \item{"spatial"}{`distances` (`dim` \eqn{\times} `dim` matrix), `kernel` (function or string: `"exponential"`, `"gaussian"`, `"matern"`, `"rational_quadratic"`, `"powered_exponential"`, `"cauchy"`), kernel-specific parameters (`phi`, `nu`, `alpha`, `beta`, `power`).}
#'   \item{"specConstr"}{`eigenvalues` (numeric vector), `rank` (integer <= `dim`), `cond_num` (numeric), optional `max_iter` and `tol`.}
#' }
#'
#' **Time series types:**
#' Non-stationary models use MA(\eqn{\infty}) expansion with differencing.
#' Roots of AR/MA polynomials are stabilized for invertibility.
#'
#' **Spectrum-Constrained constructions:**
#'            Alternating projections enforce target eigenvalues/rank/condition number while preserving diagonal.
#'            The matrix supplied in `loadings` is used to generate an initial positive semi-definite matrix
#'            for iterative projection onto the target spectrum.
#'            Default base: a `dim` \eqn{\times} `rank_target` Vandermonde-like matrix
#'            of sin(0 : \eqn{\pi}/2) values, scaled by 10
#'
#' **PSD/PD Guarantee:**
#' All matrices constructed by this function are guaranteed to be symmetric
#' and either positive semi-definite (PSD) or positive definite (PD), depending on the type and parameterization:
#' Supplied parameters are validated to ensure PSD/PD properties for all structures,  expect for `type = specConstr`.
#' In spectrum-constrained cases, particularly when `rank < dim` or user-specified `eigenvalues` are provided,
#' the matrix may be near PSD.
#' Before returning the matrix, the function performs a *final validation step* (if `eigen_check = TRUE`)
#' to verify that the output is truly PSD/PD and issues errors, warnings, or messages accordingly.
#'
#'
#' @examples
#' # Diagonal covariance
#' construct_covMat(5, type = "diag", params = list(variances = 1:5))
#'
#' # AR(1) correlation
#' construct_covMat(6, corr = TRUE, type = "AR", params = list(rho = 0.8))
#'
#' # Block matrix
#' construct_covMat(5, type = "block", params = list(groups = c(2,3), rho = c(0.5, 0.3)))
#'
#' # Factor model
#' L <- matrix(rnorm(10), 5, 2); psi <- rep(0.2, 5)
#' construct_covMat(5, type = "factor", params = list(loadings = L, psi = psi))
#'
#' # Spatial correlation with Matern kernel
#' D <- as.matrix(dist(matrix(runif(25), ncol = 5)))
#' construct_covMat(5, type = "spatial",
#'                     params = list(distances = D, kernel = "matern", phi = 1, nu = 0.5))
#'
#' # Spectrum-constrained correlation matrix
#' construct_covMat(4, type = "specConstr", params = list(eigenvalues = c(2,1,0.5,0.5)))
#'
#' @export
construct_covMat <- function(dim, corr = FALSE, type = "diag", params = list(), eigen_check = TRUE) {

  ## Validate dimension
  stopifnot(dim > 0)

  ## Validate argument 'type'
  valid_types <- c("diag", "CS", "block", "toeplitz",
                   "AR", "MA", "ARMA","ARIMA", "SAR", "SMA","SARMA","SARIMA",
                   "factor", "graph", "spatial", "specConstr")

  if(length(type) != 1 || !type %in% valid_types){
    stop(sprintf("Invalid 'type': %s. Must be one of: %s",
                 type, paste(valid_types, collapse=", ")))
  }

  ## Extract main parameters - variances and correlation
  variances <- NULL
  rho <- NULL
  if(!is.list(params)) {
    if(corr) rho <- params else variances <- params
  } else {
    variances    <- params$variances
    rho          <- params$rho
  }

  ## Adjust variances if needed (correlation = FALSE)
  if(!corr){
    if(!is.null(variances)){
      if(length(variances) == 1) variances <- rep(variances, dim)
      else if(length(variances) != dim) stop("Length of variances must match `dim` or be a scalar.")
    } else variances <- rep(1, dim)

    if(all(variances == 1) && !(type %in% c("factor", "graph"))) {
      corr = TRUE
    } else {
        D <- diag(sqrt(variances))
    }
  }


  ## Start matrix construction
  mat <- NULL

  ## Diagonal ----
  if(type == "diag") {
    mat <- diag(1, dim)
    if(!corr)  {
      mat <- D %*% mat %*% D  # Scale to covariance if requested
      mat <- (mat + t(mat)) / 2
    }
  }

  ## Compound Symmetry ----
  if(type == "CS"){
    rho <- rho[1]
    cut_off <- -1/(dim-1)
    if(rho < cut_off || rho > 1)
      stop(paste0("rho must be in [", cut_off, ", 1] for dim = ", dim))

    mat <- matrix(rho, dim, dim)
    diag(mat) <- 1

    bandwidth <- params$bandwidth %||% dim
    if(bandwidth < dim){
      ind <- 1:dim
      dist <- abs(outer(ind,ind,"-"))
      mat[dist > bandwidth] <- 0

    if(!corr)  {
      mat <- D %*% mat %*% D  # Scale to covariance if requested
      mat <- (mat + t(mat)) / 2
      }
    }
  }

  ## Block ----
  if(type == "block") {

    groups <- params$groups
    if(is.null(groups)) stop("Provide 'groups' for Block structure.")
    if(sum(groups) != dim) stop("Sum of groups must equal dim.")
    if (!requireNamespace("Matrix", quietly = TRUE)) {
      stop("Package 'Matrix' is required for constructiong block matrices. Please install it first.")
    }

    len <- length(groups)
    if (is.matrix(rho)) rho <- list(rho)

    if(is.numeric(rho)) {

      if(length(rho) == 1) rho <- rep(rho, len)
      if(length(rho) != len) stop("Length of rho must match number of groups (if not 1).")

      blocks <- lapply(seq_along(groups), function(i) {
        k <- groups[i]
        r <- rho[i]
        cut_off <- -1 / (k-1)
        if(r < cut_off || r > 1)
          stop(paste0("Block rho must be in [", cut_off, ", 1] for k =", k))

        b <- matrix(r, k, k)
        diag(b) <- 1
        b
      })
    } else if(is.list(rho)) {

      if(length(rho) == 1) rho <- rep(rho, len)
      if(length(rho) != len) stop("Length of rho list must match number of groups (if not 1).")

      blocks <- lapply(seq_along(groups), function(i) {
        rmat <- rho[[i]]
        k <- groups[i]
        if(!all(dim(rmat) == c(k,k))) stop(paste0("Block matrix size mismatch for group ", k))
        diag(rmat) <- 1
        eig <- eigen(rmat, symmetric = TRUE, only.values = TRUE)$values
        if(any(eig < 0)) stop("A block produces a non-PD matrix")
        rmat
      })
    } else stop("rho must be numeric, matrix or list of matrices.")

    # Assemble block-diagonal matrix
    mat <- as.matrix(Matrix::bdiag(blocks))

    if(!corr)  {
      mat <- D %*% mat %*% D  # Scale to covariance if requested
      mat <- (mat + t(mat)) / 2
    }
  }

  ## Toeplitz ----
  if(type=="toeplitz"){

    if(length(rho)==1) rho <- rho^(0:(dim-1))
    if(length(rho)!=dim) stop("rho length must match dim")
    mat <- toeplitz(rho)

    bandwidth <- params$bandwidth %||% dim
    if(bandwidth < dim){
      ind <- 1:dim
      dist <- abs(outer(ind,ind,"-"))
    mat[dist > bandwidth] <- 0
    }

   if(!corr)  {
      mat <- D %*% mat %*% D  # Scale to covariance if requested
      mat <- (mat + t(mat)) / 2
    }
  }


  ## ARIMA, SARIMA and its special cases AR/SAR/MA/SMA/ARMA/SARMA ----
  if (type %in% c("AR", "MA", "SAR", "SMA", "ARMA", "SARMA", "ARIMA", "SARIMA")) {

    ## Extract parameters
    phi        <- params$arma_phi        %||% as.numeric(rho[1]) %||% numeric(0)
    theta      <- params$arma_theta      %||% numeric(0)
    d          <- params$d               %||% 0
    seasonal_d <- params$D               %||% 0   # renamed
    s          <- params$seasonal_period %||% 1    # REQUIRED for SARIMA

    ## ------------ Root stabilization (corrected) ------------
    stabilize <- function(coef, invert = FALSE) {
      if (length(coef) == 0) return(coef)

      # AR polynomial: 1 - phi_1 z - ...
      # MA polynomial: 1 + theta_1 z + ...
      poly <- if (invert) c(1, coef) else c(1, -coef)
      rts <- polyroot(poly)

      # push roots outside unit circle
      if (any(Mod(rts) <= 1)) {
        rts <- rts / (Mod(rts) + 1e-8)
        # reconstruct coefficients
        newpoly <- Re(poly(rts))
        coef <- if (invert) newpoly[-1] else -newpoly[-1]
      }
      coef
    }

    phi   <- stabilize(phi, invert = FALSE)
    theta <- stabilize(theta, invert = TRUE)

    ## Seasonal components
    if (startsWith(toupper(type), "S")) {
      seasonal_phi   <- params$seasonal_phi   %||% numeric(0)
      seasonal_theta <- params$seasonal_theta %||% numeric(0)

      seasonal_phi   <- stabilize(seasonal_phi, invert = FALSE)
      seasonal_theta <- stabilize(seasonal_theta, invert = TRUE)
    } else {
      seasonal_phi <- numeric(0)
      seasonal_theta <- numeric(0)
      seasonal_d <- 0
    }

    ## Override for AR/SAR/MA/SMA
    if (type %in% c("AR", "SAR")) {
      theta <- numeric(0)
      seasonal_theta <- numeric(0)
    }
    if (type %in% c("MA", "SMA")) {
      phi <- numeric(0)
      seasonal_phi <- numeric(0)
    }

    ## ------------ Construct AR and MA polynomials ------------
    # (Φ(B^s) φ(B))  and  (Θ(B^s) θ(B))

    # Non-seasonal polynomials
    ar_poly <- c(1, -phi)
    ma_poly <- c(1,  theta)

    # Seasonal polynomials
    if (length(seasonal_phi) > 0) {
      Phi <- c(1, rep(0, s * length(seasonal_phi)))
      Phi[c(1 + s*(1:length(seasonal_phi)))] <- -seasonal_phi
      ar_poly <- stats::convolve(ar_poly, Phi, type="open")
    }
    if (length(seasonal_theta) > 0) {
      Theta <- c(1, rep(0, s * length(seasonal_theta)))
      Theta[c(1 + s*(1:length(seasonal_theta)))] <- seasonal_theta
      ma_poly <- stats::convolve(ma_poly, Theta, type="open")
    }

    ar_combined <- -ar_poly[-1]      # strip leading 1, convert back to AR coefficients
    ma_combined <-  ma_poly[-1]      # MA coefficients

    p <- length(ar_combined)
    q <- length(ma_combined)


    ## ------------ Special case: ARMA(0,0) = white noise ------------
    if (p == 0 && q == 0 && d == 0 && seasonal_d == 0) {
      mat <- diag(dim)       # correlation matrix (variance = 1)
    } else if (d == 0 && seasonal_d == 0) {
      ## ------------ Stationary ARMA covariance ------------
      rho <- stats::ARMAacf(ar = ar_combined, ma = ma_combined, lag.max = dim - 1)
      mat <- toeplitz(rho)
    } else {
      # ------------ Differenced / nonstationary: MA(∞) expansion ------------
      ma_length <- dim + 50
      psi <- stats::ARMAtoMA(ar = ar_combined, ma = ma_combined, lag.max = ma_length)

      ## Non-seasonal differencing (1 - B)^d
      if (d > 0) {
        # Base differencing polynomial: (1 - B)
        base_diff <- c(1, -1)

        # Convolve it d times to get (1 - B)^d
        difff <- base_diff
        if (d > 1) {
          for (k in 2:d) {
            difff <- stats::convolve(difff, base_diff, type = "open")
          }
        }

        # Apply non-seasonal differencing to psi
        psi <- stats::convolve(psi, difff, type = "open")
      }

      ## Seasonal differencing (1 - B^s)^D
      if (seasonal_d > 0 && s > 1) {
        # Base seasonal difference polynomial: (1 - B^s)
        base_diff <- c(1, rep(0, s - 1), -1)

        # Convolve it D times to get (1 - B^s)^D
        diffs <- base_diff
        if (seasonal_d > 1) {
          for (k in 2:seasonal_d) {
            diffs <- stats::convolve(diffs, base_diff, type = "open")
          }
        }

        # Apply seasonal differencing to psi
        psi <- stats::convolve(psi, diffs, type = "open")
      }

      psi <- psi[!is.na(psi)]

      acf_fun <- function(h) {
        if (h < length(psi))
          sum(psi[1:(length(psi)-h)] * psi[(1+h):length(psi)])
        else 0
      }
      acvs <- sapply(0:(dim-1), acf_fun)
      mat <- toeplitz(acvs)
    }

    if(!corr)  {
      mat <- D %*% mat %*% D  # Scale to covariance if requested
      mat <- (mat + t(mat)) / 2
    }
  }


  ## Factor Structure ----
  if(type == "factor"){

    loadings <- params$loadings
    psi      <- params$psi

    if(is.null(loadings) || is.null(psi))
      stop("Provide loadings and psi within 'params'.")
    if(length(psi) == 1) psi <- rep(psi, dim)
    if(length(psi) != dim) stop("psi length mismatch.")
    if(any(psi < 0)) stop("All uniquenesses 'psi' must be non-negative.")

    mat <- tcrossprod(loadings) + diag(psi)

    if (corr) {
        D_sqrt <- sqrt(diag(mat))
        if (all(D_sqrt>0)) {
          mat <- mat / (D_sqrt %o% D_sqrt)
          } else {
            warning("At least one variance is zero; returning covariance matrix.")
          }
    }
    mat <- (mat + t(mat)) / 2
  }


  ## Graph Structure ----
  if (type == "graph") {

    adj <- params$adjacency
    if (is.null(adj) || !is.matrix(adj) )
      stop("adjacency must be a dim X dim 0/1 matrix.")

    if (nrow(adj) != ncol(adj)) stop("adjacency matrix must be square.")
    if (!all(adj == t(adj))) stop("adjacency must be symmetric.")
    if (!all(adj %in% c(0,1))) stop("adjacency matrix must contain only 0/1.")

    # Precision = adjacency structure + diagonal dominance
    Q <- adj
    epsilon <- 1e-3 * max(rowSums(abs(adj)))
    diag(Q) <- rowSums(abs(Q)) + epsilon

    mat <-  chol2inv(chol(Q))   # covariance

    if (corr) {
      D <- sqrt(diag(mat))
      mat <- mat / (D %o% D)
    }
    mat <- (mat + t(mat)) / 2
  }


  ## Spatial correlation (kernel-based) ----
  if(type == "spatial") {

    distances <- params$distances %||% abs(outer(1:dim, 1:dim, "-"))

    # Validate distances matrix
    if (!is.matrix(distances) || nrow(distances) != dim || ncol(distances) != dim) {
      stop("`distances` must be a numeric square matrix of size dim X dim")
    }
    if (any(distances < 0)) stop("`distances` must be non-negative")
    if (!isSymmetric(distances, tol = sqrt(.Machine$double.eps), trans = "T")) {
      warning("`distances` is not exactly symmetric — symmetrizing it.")
      distances <- (distances + t(distances)) / 2
    }
    diag(distances) <- 0

    kernel    <- params$kernel  # either a string name or a function
    phi       <- params$phi     %||% 1
    if(phi <= 0) stop("'phi' must be positive for spatial correlation")

    # Decide kernel type
    if(is.function(kernel)) {
      # If it's a function, use it directly

      mat <- outer(1:dim, 1:dim, Vectorize(function(i, j) {kernel(distances[i, j])}))

      # PSD check for custom kernel
      eig <- eigen(mat, symmetric = TRUE, only.values = TRUE)$values
      if(any(eig < 0)) stop("Custom kernel produced a non-PSD matrix")
      if(any(eig == 0)) message("Custom kernel produces zero eigenvalues.")

    } else if (is.character(kernel) && length(kernel) == 1) {
      # Predefined kernels by name
      kernel <- tolower(kernel)

      if(kernel == "exponential" || kernel == "exp" ) mat <- exp(-distances / phi)

      if(kernel == "gaussian" || kernel == "norm") mat <- exp(-(distances^2) / (2 * phi^2))

      if(kernel == "matern") {
        nu        <- params$nu      %||% 0.5
        if(nu <= 0) stop("'nu' must be > 0 for Matern kernel in Spatial Correlation.")

        arg <- sqrt(2 * nu) * distances / phi
        mat <- (2^(1-nu) / gamma(nu)) * (arg^nu) * besselK(arg, nu)
        mat[distances == 0] <- 1
      }

      if(kernel == "rational_quadratic" || kernel == "rq") {
        alpha     <- params$alpha   %||% 1

        if(alpha <= 0)
          stop("'alpha' must be > 0 for rational-quadratic kernel in Spatial Correlation.")
        mat <- (1 + (distances^2) / (2 * alpha * phi^2))^(-alpha)
      }

      if(kernel == "powered_exponential" || kernel == "pexp") {
        power     <- params$power   %||% 1.5
        if(power <= 0 || power > 2)
          stop("Powered exponential requires 0 < power <= 2 in Spatial Correlation." )
        mat <- exp(-(distances / phi)^power)
      }

      if(kernel == "cauchy") {
        alpha     <- params$alpha   %||% 1
        beta      <- params$beta    %||% 1
        if(alpha <= 0) stop("'alpha' must be > 0 for Cauchy kernel in Spatial Correlation.")
        if(beta <= 0 || beta > 2)
          stop("'beta' must be in (0,2] for PSD Cauchy kernel in Spatial Correlation.")
        mat <- (1 + (distances / phi)^beta)^(-alpha)
      }
      # diag(mat) <- 1

    } else {
      stop("`kernel` must be either a string name of a kernel, or a function.")
    }

    if (!is.null(mat)) mat <- unname(mat) else stop("Invalid kernel name.")

    if(!corr)  {
      mat <- D %*% mat %*% D  # Scale to covariance if requested
      mat <- (mat + t(mat)) / 2
    }
  }


  ## Spectrum constrained Corr/Cov --- ----------------
  if (type == "specConstr") {

    max_iter <- params$max_iter %||% 1000
    tol <- params$tol %||% 1e-10
    eigenvalues <- sort(params$eigenvalues, decreasing = TRUE)

    rank_target <- if (!is.null(eigenvalues)) sum(abs(eigenvalues) > tol) else {params$rank %||% dim}
    if (rank_target < 1 || rank_target > dim)
      stop("rank must be in 1:dim")

    base <- params$loadings %||% {10 * sapply(0:(rank_target-1), function(k) sin(seq(0, pi/2, length.out = dim))^k)}
    if ((nrow(base) != dim) || (ncol(base)!=rank_target))
      stop("Dimension mismatch in supplied 'loadings' matrix.")

    if (!is.null(eigenvalues)){

      # Pad or truncate eigenvalues to match dimension
      p <- length(eigenvalues)
      if (p < dim) eigenvalues <- c(eigenvalues, rep(0, dim - p))
      if (p > dim) eigenvalues <- eigenvalues[1:dim]

      target_diag <- if (corr) rep(1, dim) else variances

      # Correlation matrices require trace = sum(diag)
      if (abs(sum(eigenvalues) - sum(target_diag)) > 1e-12)
        stop("Eigenvalues must sum to the sum of diagonals of a correlation/covariance matrix.")

      Lambda <- diag(eigenvalues)

      # Random PSD initialization
      A <- tcrossprod(base)
      R <- A / mean(diag(A))   # scale diagonal roughly to 1

      for (iter in seq_len(max_iter)) {

        R[!is.finite(R)] <- 0

        # Project onto correlation/covariance diagonal
        d <- sqrt(diag(R))
        # d[d == 0] <- 1e-12   # prevent division by zero
        R <- R / (d %o% d)
        R <-  diag(sqrt(target_diag)) %*% R %*% diag(sqrt(target_diag))

        # Project onto target eigenvalues
        eigA <- eigen(R, symmetric = TRUE)
        U <- eigA$vectors
        R <- U %*% Lambda %*% t(U)

        # Convergence check
        if (max(abs(diag(R) - target_diag)) < tol) break
      }
    } else {

      cond_num   <- params$cond_num
      if (!is.null(cond_num)){
        if (!corr) cond_num <- cond_num * min(variances)/max(variances)
        if (cond_num < 1) stop("Condition number must be >= 1")
      }

      # Initial PSD matrix of desired rank
      A <- tcrossprod(base)
      R <- A / mean(diag(A))   # scale diagonal roughly to 1

      for (iter in seq_len(max_iter)) {
        # PSD + rank projection
        eig1 <- eigen(R, symmetric = TRUE)
        vals1 <- eig1$values
        vecs1 <- eig1$vectors
        # vals1[vals1 < 0 & vals1 > -1e-12] <- 0

        # Keep exactly rank_target nonzero eigenvalues (>= eps)
        vals1_proj <- c(pmax(vals1[1:rank_target], 1e-8), rep(0, length(vals1) - rank_target))
        R1 <- vecs1 %*% diag(vals1_proj) %*% t(vecs1)

        # Scale to unit diagonal (correlation)
        d <- sqrt(diag(R1))
        R <- R1 / (d %o% d)

        ## Step 3: Enforce condition number on the non-zero eigenvalue subspace
        if (!is.null(cond_num)) {
          eig2 <- eigen(R, symmetric = TRUE)
          vecs2 <- eig2$vectors
          vals2 <- eig2$values
          idx <- order(vals2, decreasing = TRUE)
          vals2_sorted <- vals2[idx]

          # Only consider the top `rank_target` eigenvalues for condition number
          nonzero_vals <- vals2_sorted[1:rank_target]
          mini <- min(nonzero_vals)
          maxi <- max(nonzero_vals)
          target_max <- cond_num * mini
          # Rescale the non-zero block linearly
          if (maxi > mini) {
            nonzero_vals_new <- mini + (nonzero_vals - mini) * ((target_max - mini) / (maxi - mini))
          } else {
            nonzero_vals_new <- nonzero_vals  # all equal, nothing to scale
          }

          # Put back into full eigenvalue vector (including zeros)
          vals2_new_sorted <- c(nonzero_vals_new, rep(0, length(vals2) - rank_target))
          # Reorder to match original ordering
          vals2_new <- vals2_new_sorted[order(idx)]

          R <- vecs2 %*% diag(vals2_new) %*% t(vecs2)
        }

        ## Convergence check: how close is the diagonal to the target diagonals
        if (max(abs(diag(R) - 1)) < tol) break
      }
    }

    if (iter == max_iter) warning("Did not converge in max_iter for the given spectrum constraints.")
    mat <- R
    if (is.null(eigenvalues) && (!corr)) mat <-  diag(sqrt(variances)) %*% mat %*% diag(sqrt(variances))

    mat <- (mat +t(mat))/2
  }

  # PSD check
  eig <- eigen(mat, symmetric=TRUE, only.values=TRUE)$values
  if(any(eig < -1e-8)) stop("The constructed matrix not PSD.")
  if(any(eig < 0)) {
    warning("The constructed matrix is near PSD.")
  } else {
    if(any(eig == 0)) warning("The constructed matrix is PSD.") else message("The constructed matrix is PD.")
    }

  class(mat) <- c("covMat", class(mat))
  return(mat)
}




