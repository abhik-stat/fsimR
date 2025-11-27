#' @title Simulate Data from a Linear Mixed Model (LMM)
#' @md
#'
#' @description
#' Generates datasets from a Linear Mixed Model (LMM)
#' with flexible control over fixed effects, random effects,
#' correlation structures, design matrices, and error distributions.
#'
#' All stochastic components (random effects, error terms,
#' and design matrices when unspecified)  are generated via \code{\link{simulate_IIDdata}},
#' a **highly general and flexible sample generator**
#' supporting almost all univariate and multivariate distributions.
#'
#' The function supports general, possibly correlated and/or heteroscedastic,
#' error distributions. When the error distribution *is*
#' independent and identically distributed (**IID**, meaning all error terms share
#' the same distribution and are mutually independent),
#' users may specify a target signal-to-noise ratio (SNR).
#' In this case, the error variance is automatically scaled to match the desired SNR.
#' SNR control is not supported for non-IID errors.
#'
#'
#'
#' @param n_rep Integer scalar or length-2 vector.
#' If a vector of length 2, the first value gives the number of
#' independent simulated datasets (`n_sim`)
#' and the second gives the number of within-dataset iterations (`n_iter`).
#' Within a dataset, the design matrices `X` and `Z` remain fixed
#' but random effects (`RE`) and error terms vary across iterations,
#' producing different response vectors (`y`).
#' If a single value is given, it is interpreted as `n_sim`, with `n_iter = 1` (default).
#'
#' @param sim_settings A named list of simulation parameters,
#' usually constructed via \code{\link{sim_spec}} (default).
#' Common fields (all optional) for the mixed models include:
#' \describe{
#'   \item{n_subj}{Number of subjects or clusters.}
#'
#'   \item{n_obs}{Integer scalar or vector (>0).
#'   Specifies the number of observations per subject/cluster.
#'   - If a **scalar** is supplied, all `n_subj` clusters have the same number of observations.
#'   - If a **vector**, its length must be `n_subj`, and each entry gives the size of the respective cluster.
#'   Default: `5` observations per cluster/subject.}
#'
#'   \item{beta_coeff}{Fixed-effect coefficient vector.
#'   Default: a zero vector of appropriate length.}
#'
#'   \item{SNR}{Target signal-to-noise ratio for IID errors.
#'       Overrides the error variance in \code{distr_settings$error_distr}.
#'       Set to \code{NULL} (default) to disable SNR control.}
#'
#'   \item{X, Z}{Optional pre-specified design matrices.
#'               May also be lists of length `n_sim` to allow dataset-specific designs.
#'               If supplied, they override stochastic generation.}
#'
#'   \item{p, q}{Dimensions of fixed- and random-effect design matrices
#'              (including intercept).}
#'
#'   \item{is.CorrelatedXZ}{Logical; If `TRUE`, `X` and `Z` are generated *jointly*
#'                           from `distr_settings$XZ_distr`,
#'                           and any user-supplied \code{X} or \code{Z} is ignored.}
#'
#'   \item{is.ZInX}{Logical. If `TRUE`, and `Z` is not supplied, the first `q`
#'       columns of `X` (after intercept handling) are used to construct `Z`.}
#'
#'  \item{include.Xintercept, include.Zintercept}{Logical;
#'      whether to add intercept columns in `X`/`Z`.}
#'
#'   \item{orthogonalize.X, orthogonalize.Z}{Logical;
#'       If `TRUE`, the non-intercept columns of `X` and/or `Z`
#'       are orthogonalized using QR decomposition.}
#' }
#'
#' @param distr_settings A named list of distribution specifications for stochastic components.
#' Valid entries are: \code{error_distr}, \code{RE_distr}, \code{X_distr},
#' \code{Z_distr}, and \code{XZ_distr}, specifying the distributions used to
#' generate errors, random effects, and design matrices, respectively.
#' Each entry is itself a list with elements (see \code{\link{simulate_IIDdata}}):
#' \describe{
#'   \item{distr_name}{Distribution name (e.g., \code{"norm"}, \code{"mvnorm"}, \code{"copula"}).}
#'   \item{distr_params}{List of parameters for the specified distribution.}
#'   \item{generator}{Optional user-supplied function for random generation,
#'   enabling simulation from arbitrary or fully custom distributions.}
#' }
#' Missing or incomplete entries are automatically completed via
#'   \code{\link{distr_spec}} with suitable dimension defaults.
#'   The default distribution is multivariate normal with zero mean and identity covariance.
#'   See @details for more on specification of the error distribution.
#'
#' @param seed Optional integer. If provided, the random number generator (RNG)
#'    state is temporarily replaced by this seed and restored upon exit.
#'
#'
#' @return
#' A named list (of class \code{"LMMdata"}) containing:
#' \describe{
#'    \item{n}{Total number of observations in each simulated dataset.}
#'   \item{y}{List of simulated response matrices of length `n_sim`.
#'            Each matrix has dimension `n × n_iter`.}
#'   \item{X, Z}{List of design matrices of length `n_sim`,
#'            one per simulated dataset.}
#'   \item{RE}{List of random effect realizations of length `n_sim`.
#'            Each entry is either a matrix (when `n_iter = 1`) or
#'            a list of matrices (one per iteration).
#'            These matrices are of order `n_subj × q`.}
#'   \item{sigma_e}{List of error standard deviations (length `n_sim`).
#'            Each entry is a vector of length `n_iter`.}
#'   \item{SNR}{List of achieved signal-to-noise ratios (length `n_sim`).
#'            Each entry is a vector of length `n_iter`.}
#'   \item{ID}{Subject identifiers.
#'        A factor of length equal to the total number of observations,
#'        with \code{n_subj} levels, that aligns with the ordering of rows
#'        in \code{y}, \code{X}, and \code{Z}.
#'        Each level corresponds to one subject/cluster,
#'        and repeated occurrences of a level indicate multiple observations for that subject/cluster.
#'        }
#'   \item{sim_settings, distr_settings}{Lists of (possibly updated) simulation and
#'     distribution settings used to generate the data.}
#' }
#' If a list has only one single element (i.e., `n_sim = 1` or `n_iter  = 1`),
#' it is automatically flattened.
#'
#'
#' @details
#' **Construction of design matrices.**
#' The design matrices `X` and `Z` are constructed in the following order:
#' \itemize{
#'    \item If \code{is.CorrelatedXZ = TRUE},
#'        `X` and `Z` are generated *jointly* from  the distribution
#'        specified in \code{distr_settings$XZ_distr}.
#'        User-supplied `X` or `Z` are **ignored** in this case.
#'
#'    \item If `is.CorrelatedXZ = FALSE`, user-supplied matrices take precedence:
#'        - If \code{sim_settings$X} is supplied: `X` is used exactly as provided
#'        and its column dimension updates `p`.
#'        - If \code{sim_settings$Z} is supplied: `Z` is used exactly as provided
#'        and its column dimension updates `q`.
#'
#'    \item If `X` is not supplied: `X` is generated using `distr_settings$X_distr`.
#'
#'    \item If `Z` is not supplied:
#'        - If `is.ZInX = FALSE`, `Z` is generated from \code{distr_settings$Z_distr};
#'        - otherwise `Z` is inherited from the *first `q` columns of `X`*
#'        (with correct intercept handling).
#'
#'    \item Intercepts are added **after** stochastic generation, if requested.
#'    \item Orthogonalization of non-intercept columns is performed last
#'        using the \code{\link[base]{qr}} function.
#'    }
#'
#' **Random effects**
#' Random effects are drawn per subject from `distr_settings$RE_distr`,
#' with dimensions matched to the columns of `Z`.
#'
#' **Error distribution specification**
#' The error distribution (\code{error_distr}) may be:
#' \itemize{
#'   \item \strong{Univariate} with scalar \code{sigma},
#'          which corresponds to IID errors.
#'
#'   \item \strong{Multivariate}, with a covariance matrix `sigma`
#'      of dimension equal to the total number of observations,
#'      enabling arbitrary correlation and heteroscedasticity.
#' }
#' IID errors can also be represented through a multivariate specification
#' with diagonal covariance matrix in `sigma`.
#' If \code{sigma} is not specified in `error_distr$distr_params`,
#' a default IID specification is constructed.
#'
#'
#' **SNR control**
#' When `SNR` is supplied and the error distribution is IID,
#' the error variance is automatically scaled to achieve
#' the target signal-to-noise ratio
#' (defined as the ratio of variances of signals and errors).
#' SNR control is not supported for non-IID errors.
#'
#'
#' **Dependencies**
#' Certain distribution choices (e.g., copulas) may require additional packages.
#' See \code{\link{simulate_IIDdata}} for details on required packages
#' and supported distributions.
#'
#'
#' @examples
#'
#' # Basic simulation
#' set.seed(123)
#' sim1 <- simulate_LMMdata()
#' str(sim1, max.level = 1)
#'
#'
#' # Multiple replications and iterations
#' sim <- simulate_LMMdata(n_rep = c(2, 3)) # 2 reps × 3 iterations
#' length(sim$y)         # 2 replications
#' dim(sim$y[[1]])       # 50 × 3 iterations (default sample size 50)
#'
#'
#' # Non-Gaussian error and multivariate random effects
#' distr_settings <- list(
#'   error_distr = list(distr_name = "t", distr_params = distr_spec(df = 4)),
#'   RE_distr    = list(distr_name = "mvnorm", distr_params = list(dim = 2, sigma = diag(2)))
#' )
#' sim <- simulate_LMMdata(
#'   n_rep = 1,
#'   sim_settings = sim_spec(n_subj = 15, n_obs = 5),
#'   distr_settings = distr_settings
#' )
#'
#'
#' # Correlated X and Z design matrices
#' Sigma_X <- matrix(c(1, 0.5, 0.5, 1), ncol = 2)
#' Sigma_Z <- diag(2)
#' Sigma_XZ <- matrix(c(0.3, 0.2, 0, 0.4), ncol = 2)
#' Sigma <- rbind(cbind(Sigma_X, Sigma_XZ), cbind(t(Sigma_XZ), Sigma_Z))
#'
#' sim_settings <- sim_spec(
#'   n_subj = 10,
#'   n_obs  = 4,
#'   is.CorrelatedXZ = TRUE,
#'   include.Xintercept = TRUE,
#'   include.Zintercept = TRUE
#' )
#' distr_settings <- list(
#'   X_distr  = list(distr_name = "mvnorm", distr_params = list(sigma = Sigma_X)),
#'   Z_distr  = list(distr_name = "mvnorm", distr_params = list(sigma = Sigma_Z)),
#'   XZ_distr = list(distr_name = "mvnorm", distr_params = list(sigma = Sigma))
#' )
#' sim <- simulate_LMMdata(sim_settings = sim_settings, distr_settings = distr_settings)
#' pairs(cbind(sim$X, sim$Z), main = "Correlated X predictors")
#'
#'
#' # User-supplied design matrices
#' n_subj <- 5; n_obs <- 6; n <- n_subj * n_obs
#' X_user <- matrix(rnorm(n*3), ncol=3)
#' Z_user <- matrix(rnorm(n*2), ncol=2)
#' sim_settings <- sim_spec(n_subj = n_subj, n_obs = n_obs, X = X_user, Z = Z_user, beta_coeff = c(1,0.5,-1))
#' sim <- simulate_LMMdata(sim_settings = sim_settings)
#'
#'
#' # Copula-based generation
#' library(copula)
#' normal_cop <- normalCopula(param = 0.6, dim = 2)
#' distr_settings <- list(
#'   X_distr = list(
#'     distr_name = "copula",
#'     distr_params = distr_spec(
#'       copula = normal_cop,
#'       margins = list(
#'         list(dist="norm", params=list(mean=0, sd=1)),
#'         list(dist="norm", params=list(mean=0, sd=2))
#'       )
#'     )
#'   )
#' )
#' sim <- simulate_LMMdata(sim_settings = sim_spec(n_subj = 10, n_obs = 3),
#'                          distr_settings = distr_settings)
#' plot(sim$X[,2:3], main = "X from Copula-based distribution")
#'
#'
#' @seealso
#' \code{\link{simulate_IIDdata}}, \code{\link{distr_spec}}, \code{\link{sim_spec}}.
#' Also \code{\link{summary.LMMdata}}, \code{\link{print.LMMdata}}, \code{\link{plot.LMMdata}}
#' for S3 methods applicable to the returned object
#'
#' @export
simulate_LMMdata <- function(n_rep = 1, sim_settings = list(),
                             distr_settings = list(), seed = NULL) {

  ## Save RNG state when seed is provided
  if (!is.null(seed)) {
    old_seed <- .GlobalEnv$.Random.seed %||% NULL
    set.seed(seed)
    on.exit({
      if (!is.null(old_seed)) .GlobalEnv$.Random.seed <- old_seed
      else rm(".Random.seed", envir = .GlobalEnv)
    }, add = TRUE)
  }


  ## Update default arguments with user-supplied arguments
  sim_settings <- do.call(sim_spec, sim_settings)

  p_from_distr <- ncol(distr_settings$X_distr$distr_params$sigma) + sim_settings$include.Xintercept

  q_from_distr <-
    ncol(distr_settings$RE_distr$distr_params$sigma) %||%
    (ncol(distr_settings$Z_distr$distr_params$sigma) + sim_settings$include.Zintercept)

  p <- max(sim_settings$p, p_from_distr, ncol(sim_settings$X))
  q <- max(sim_settings$q, q_from_distr, ncol(sim_settings$Z))
  p0 <- max(0, p - sim_settings$include.Xintercept)
  q0 <- max(0, q - sim_settings$include.Zintercept)

  X_distr <- distr_settings$X_distr
  Z_distr <- distr_settings$Z_distr
  XZ_distr <- distr_settings$XZ_distr
  RE_distr <- distr_settings$RE_distr
  error_distr <- distr_settings$error_distr

  X_distr$distr_name <- X_distr$distr_name %||% "mvnorm"
  Z_distr$distr_name <- Z_distr$distr_name %||% "mvnorm"
  XZ_distr$distr_name <- XZ_distr$distr_name %||% "mvnorm"
  RE_distr$distr_name <- RE_distr$distr_name %||% "mvnorm"

  if (X_distr$distr_name!="copula")
    X_distr$distr_params <- do.call(distr_spec, X_distr$distr_params %||% list(dim = p0))
  if (Z_distr$distr_name!="copula")
    Z_distr$distr_params <- do.call(distr_spec, Z_distr$distr_params %||% list(dim = q0))
  if (XZ_distr$distr_name!="copula")
    XZ_distr$distr_params <- do.call(distr_spec, XZ_distr$distr_params %||% list(dim = p0 + q0))
  if (RE_distr$distr_name!="copula")
    RE_distr$distr_params <- do.call(distr_spec, RE_distr$distr_params %||% list(dim = q))

  error_distr$distr_name <- error_distr$distr_name %||% "norm"
  error_distr$distr_params <- do.call(distr_spec, error_distr$distr_params %||% list(dim = 1))
  error_params_WOsd <- error_distr$distr_params[names(error_distr$distr_params)!="sd"]
  error_sigma <- error_distr$distr_params$sd

  ## Extract Other Settings and parameters
  n_rep <- c(n_rep, 1)
  sim_num <- n_rep[1]
  it_num <- n_rep[2]

  nrsubj <- sim_settings$n_subj
  nrobs <- sim_settings$n_obs
  if (length(nrobs) == 1) nrobs <- rep(nrobs, nrsubj)
  n <- sum(nrobs)
  if (n <= 0) stop("Error! n must be a positive integer")

  subj_idx <- rep(seq_len(nrsubj), times = nrobs)
  ID <- factor(subj_idx)

  beta_coeff <- sim_settings$beta_coeff
  SNR <- sim_settings$SNR

  SNR_control <- FALSE
  IID_error <- (identical(as.matrix(error_sigma), t(error_sigma)))  &&
              (all(error_sigma[upper.tri(error_sigma)] == 0)) &&
              (length(unique(diag(error_sigma))) == 1)
  if (!is.null(SNR)){
    if (!IID_error) stop("SNR control is allowed only for IID error distribution!")
    SNR_control <- TRUE
  }
  if (IID_error) error_sigma <- as.numeric(error_sigma[1,1])


  ## Simulate data
  Xs <- vector("list", sim_num)
  Zs <- vector("list", sim_num)
  yss <- vector("list", sim_num)
  bss <- vector("list", sim_num)
  sigma_list <- vector("list", sim_num)
  SNRs <- vector("list", sim_num)


  for (si in seq_len(sim_num)) {

    if(sim_settings$is.CorrelatedXZ) {
      # Generate X and Z together
      xz_data <- simulate_IIDdata(n, XZ_distr$distr_name, XZ_distr$distr_params,
                                  generator = XZ_distr$generator)

      if (p0 > 0) X <- xz_data[, seq_len(p0), drop=FALSE] else X <- NULL
      if (q0 > 0) Z <- xz_data[, p0 + seq_len(q0), drop=FALSE] else Z <- NULL

      if (sim_settings$include.Xintercept) X <- cbind(1, X)
      if (sim_settings$include.Zintercept) Z <- cbind(1, Z)

    } else {

      # Generate X if required
      X_list <- sim_settings$X

      if (!is.null(X_list)) X <- if(is.list(X_list)) X_list[[si]] else X_list
      else {
        X <- if (sim_settings$include.Xintercept)  matrix(1, n, 1) else NULL
        if (p0 > 0) X <- cbind(X,  simulate_IIDdata(n, distr_name = X_distr$distr_name,
                                                    distr_params  = X_distr$distr_params,
                                                    generator     = X_distr$generator))
      }

      # Generate Z as required
      Z_list <- sim_settings$Z
      if (!is.null(Z_list)) Z <- if(is.list(Z_list)) Z_list[[si]] else Z_list
      else if (!sim_settings$is.ZInX) {
        Z <- if (sim_settings$include.Zintercept)  matrix(1, n, 1) else NULL
        if (q0 > 0) Z <- cbind(Z,  simulate_IIDdata(n, distr_name = Z_distr$distr_name,
                                                    distr_params  = Z_distr$distr_params,
                                                    generator     = Z_distr$generator))
      } else {
        if (sim_settings$include.Xintercept) {
          Z <- if (sim_settings$include.Zintercept) X[,seq_len(q)] else X[,seq_len(q)+1]
        } else
          Z <- if (!sim_settings$include.Zintercept) X[,seq_len(q0)] else cbind(1, X[,seq_len(q0)])
      }
    }


    # Dimension checks
    beta_coeff <- beta_coeff %||% rep(0, ncol(X))
    if (length(beta_coeff) != ncol(X)) {
      stop(sprintf("Length of beta_coeff (%d) does not match number of X columns (%d)",
                   length(beta_coeff), ncol(X)))
    }
    if (ncol(Z) != ncol(RE_distr$distr_params$sigma)) {
      stop(sprintf("Number of columns in Z (%d) does not match random effects dim (%d)",
                   ncol(Z), ncol(RE_distr$distr_params$sigma)))
    }


    # Orthogonalize if required
    if (sim_settings$orthogonalize.X && p0 > 0) {
      if (sim_settings$include.Xintercept) {
        X1 <- X[,1,drop=FALSE]
        Xrest <- X[,-1,drop=FALSE]
        X <- cbind(X1, qr.Q(qr(Xrest)))
      } else {
        X <- qr.Q(qr(X))
      }
    }

    if (sim_settings$orthogonalize.Z && q0 > 0) {
      if (sim_settings$include.Zintercept) {
        Z1 <- Z[,1,drop=FALSE]
        Zrest <- Z[,-1,drop=FALSE]
        Z <- cbind(Z1, qr.Q(qr(Zrest)))
      } else {
        Z <- qr.Q(qr(Z))
      }
    }

    # Name the columns of X and Z if required
    if (is.null(colnames(Z))) colnames(Z) <- if (sim_settings$include.Zintercept) paste0("Z", 0:(q - 1)) else paste0("Z", 1:q)
    if (is.null(colnames(X))) colnames(X) <- if (sim_settings$include.Xintercept) paste0("X", 0:(ncol(X) - 1)) else paste0("X", 1:ncol(X))


    # Generate responses
    ys_mat <- matrix(NA_real_, nrow = n, ncol = it_num)
    bs_list <- vector("list", it_num)
    local_sigmas <- numeric(it_num)
    etas_var <- numeric(it_num)

    # Prepare error parameters
    if (!SNR_control) {
      if (IID_error){
        error_params_updated <- do.call(distr_spec, c(list(sd = error_sigma), error_params_WOsd))
        eps <- matrix(simulate_IIDdata(n*it_num, distr_name = error_distr$distr_name,
                                       distr_params  = error_params_updated,
                                       generator     = error_distr$generator),
                      nrow = n, ncol = it_num)
      } else {
        eps <- t(simulate_IIDdata(it_num, distr_name = error_distr$distr_name,
                                  distr_params  = error_distr$distr_params,
                                  generator     = error_distr$generator))
      }
    }

    # Loop over n_iter
    for (it in seq_len(it_num)) {

      # Simulate Random Effects
      b_mat <- simulate_IIDdata(nrsubj,
                                distr_name   = RE_distr$distr_name,
                                distr_params = RE_distr$distr_params,
                                generator    = RE_distr$generator)

      if (!all(dim(b_mat) == c(nrsubj, ncol(Z)))) {
        stop("Random effects matrix dimension does not match Z columns and subjects")
      }
      bs_list[[it]] <- b_mat

      # Linear Predictor
      eta <- drop(X %*% beta_coeff) + rowSums(Z * b_mat[subj_idx, ,drop = FALSE])
      etas_var[it] <- var(eta)

      # Simulate error and response with SNR control for IID errors
      if(SNR_control){
        local_sigma <- sqrt(var(eta) / SNR)
        error_params_updated <- do.call(distr_spec, c(list(sd=local_sigma), error_params_WOsd))
        eps <- simulate_IIDdata(n, distr_name = error_distr$distr_name,
                                distr_params  = error_params_updated,
                                generator     = error_distr$generator)
        ys_mat[, it] <- eta + eps
        local_sigmas[it] <- local_sigma
      } else {
        ys_mat[, it]<- eta + eps[,it]
        local_sigmas[it] <- if (IID_error) error_sigma else sqrt(mean(diag(error_sigma)))
      }
    }


    # Store Results
    class(X) <- c("IIDdata", class(X))
    class(Z) <- c("IIDdata", class(Z))
    Xs[[si]] <- X
    Zs[[si]] <- Z
    yss[[si]] <- ys_mat
    bss[[si]] <- if (it_num == 1) bs_list[[1]] else bs_list
    sigma_list[[si]] <- local_sigmas
    SNRs[[si]] <- etas_var / (local_sigmas^2)
  }

  ## Output
  out <- list(
    y = yss, X = Xs, Z = Zs, RE = bss,
    sigma_e = sigma_list, SNR = SNRs
  )

  if (sim_num == 1) out <- lapply(out, `[[`, 1)
  out$ID=ID
  out$n <- n


  sim_settings$p <- p
  sim_settings$q <- q
  sim_settings$beta_coeff <- beta_coeff
  out$sim_settings <- sim_settings

  distr_settings$X_distr <- X_distr
  distr_settings$Z_distr <- Z_distr
  distr_settings$XZ_distr <- XZ_distr
  distr_settings$RE_distr <- RE_distr
  distr_settings$error_distr <- error_distr
  out$distr_settings = distr_settings

  class(out) <- c("LMMdata", class(out))
  return(out)
}

