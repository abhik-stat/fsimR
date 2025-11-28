#' @title Simulate Data from a Linear Regression Model
#' @md
#'
#' @description
#' Generates datasets from a Linear Model (LM)
#' with flexible control over design matrix,
#' correlation structures, and error distributions.
#'
#' All stochastic components (error terms, and also design matrices if unspecified)
#' are generated via [simulate_IIDdata],
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
#' Within a dataset, the design matrix `X` remains fixed
#' but the error terms vary across iterations,
#' producing different response vectors (`y`).
#' If a single value is given, it is interpreted as `n_sim`, with `n_iter = 1` (default).
#'
#' @param sim_settings A named list of simulation parameters,
#' usually constructed via \code{\link{sim_spec}} (Default).
#' Common fields (all optional) for the regression models include:
#' \describe{
#'   \item{n_subj}{Number of subjects (sample size here).}
#'
#'   \item{beta_coeff}{Regression coefficient vector.
#'   Default: a zero vector of appropriate length.}
#'
#'   \item{SNR}{Target signal-to-noise ratio for IID errors.
#'       Overrides the error variance in \code{distr_settings$error_distr}.
#'       Set to \code{NULL} (default) to disable SNR control.}
#'
#'   \item{X}{Optional pre-specified design matrix.
#'            May also be lists of length `n_sim` to allow dataset-specific designs.
#'            If supplied, they override stochastic generation.}
#'
#'   \item{p}{Dimensions of the design matrix (including intercept).}
#'
#'  \item{include.Xintercept}{Logical; whether to add intercept columns in `X`.}
#'
#'   \item{orthogonalize.X}{Logical; If `TRUE`, the non-intercept columns of `X`
#'       are orthogonalized using QR decomposition.}
#' }
#' All other components of `sim_settings` are ignored for regression models.
#'
#'
#' @param distr_settings A named list of distribution specifications for stochastic components.
#' Valid entries are: \code{error_distr} and \code{X_distr},
#' specifying the distributions used to generate random errors,
#' and the design matrix, respectively.
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
#' @param seed Optional integer.If provided, the random number generator (RNG)
#'    state is temporarily replaced by this seed and restored upon exit.
#'
#'
#' @return
#' A named list (of class \code{"LMdata"}) containing:
#' \describe{
#'    \item{n}{Total number of observations in each simulated dataset.}
#'   \item{y}{List of simulated response matrices of length `n_sim`.
#'            Each matrix has dimension `n × n_iter`.}
#'   \item{X}{List of design matrices of length `n_sim`,
#'            one per simulated dataset.}
#'   \item{sigma_e}{List of error standard deviations (length `n_sim`).
#'            Each entry is a vector of length `n_iter`.}
#'   \item{SNR}{List of achieved signal-to-noise ratios (length `n_sim`).
#'            Each entry is a vector of length `n_iter`.}
#'   \item{sim_settings, distr_settings}{Lists of (possibly updated) simulation and
#'     distribution settings used to generate the data.}
#' }
#' If a list has only one single element (i.e., `n_sim = 1` or `n_iter  = 1`),
#' it is automatically flattened.
#'
#'
#' @details
#' **Construction of the design matrix.**
#'  - If \code{sim_settings$X} is supplied: `X` is used exactly as provided
#'        and its column dimension updates `p`.
#'  - If `X` is not supplied: `X` is generated using `distr_settings$X_distr`.
#'  - Intercepts are added **after** stochastic generation, if requested.
#'  - Orthogonalization of non-intercept columns is performed last
#'        using the \code{\link[base]{qr}} function.
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
#' sim1 <- simulate_LMdata()
#' str(sim1, max.level = 1)
#'
#'
#' # Multiple replications and iterations
#' sim <- simulate_LMdata(n_rep = c(2, 3)) # 2 reps × 3 iterations
#' length(sim$y)         # 2 replications
#' dim(sim$y[[1]])       # 50 × 3 iterations (default sample size 50)
#'
#'
#' # Non-Gaussian error and multivariate predictors
#' distr_settings <- list(
#'   error_distr = list(distr_name = "t", distr_params = distr_spec(df = 4)),
#'   X_distr    = list(distr_name = "mvtnorm::rmvnorm", distr_params = list(dim = 2, sigma = diag(2)))
#' )
#' sim <- simulate_LMdata(
#'   n_rep = 1,
#'   sim_settings = sim_spec(n_subj = 15),
#'   distr_settings = distr_settings
#' )
#'
#' # User-supplied design matrices
#' n_subj <- 10
#' n <- n_subj
#' X_user <- matrix(rnorm(n*3), ncol=3)
#' Z_user <- matrix(rnorm(n*2), ncol=2)
#' sim_settings <- sim_spec(n_subj = n_subj, X = X_user, Z = Z_user,
#'                           beta_coeff = c(1,0.5,-1))
#' sim <- simulate_LMdata(sim_settings = sim_settings)
#'
#'
#' # Copula-based generation
#' if (requireNamespace("copula", quietly = TRUE)) {
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
#' sim <- simulate_LMdata(sim_settings = sim_spec(n_subj = 10),
#'                          distr_settings = distr_settings)
#' plot(sim$X[,2:3], main = "X from Copula-based distribution")
#' }
#'
#' @seealso
#' \code{\link{simulate_IIDdata}}, \code{\link{distr_spec}}, \code{\link{sim_spec}}.
#' Also \code{\link{summary.LMdata}}, \code{\link{print.LMdata}}, \code{\link{plot.LMdata}}
#' for S3 methods applicable to the returned object
#'
#' @export
simulate_LMdata <- function(n_rep = 1, sim_settings = list(),
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

  p <- max(sim_settings$p, p_from_distr, ncol(sim_settings$X))
  p0 <- max(0, p - sim_settings$include.Xintercept)

  X_distr <- distr_settings$X_distr
  error_distr <- distr_settings$error_distr

  X_distr$distr_name <- X_distr$distr_name %||% "mvnorm"
  if (X_distr$distr_name!="copula")
    X_distr$distr_params <- do.call(distr_spec, X_distr$distr_params %||% list(dim = p0))

  error_distr$distr_name <- error_distr$distr_name %||% "norm"
  error_distr$distr_params <- do.call(distr_spec, error_distr$distr_params %||% list(dim = 1))
  error_params_WOsd <- error_distr$distr_params[names(error_distr$distr_params)!="sd"]
  error_sigma <- error_distr$distr_params$sd

  ## Extract Other Settings and parameters
  n_rep <- c(n_rep, 1)
  sim_num <- n_rep[1]
  it_num <- n_rep[2]

  n <- sim_settings$n_subj
  if (n <= 0) stop("Error! n must be a positive integer")


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


  ## Simulate data
  Xs <- vector("list", sim_num)
  yss <- vector("list", sim_num)
  sigma_list <- vector("list", sim_num)
  SNRs <- vector("list", sim_num)


  for (si in seq_len(sim_num)) {

    # Generate X if required
    X_list <- sim_settings$X

    if (!is.null(X_list)) {
      X <- if(is.list(X_list)) X_list[[si]] else X_list
      n <- nrow(X)
    } else {
      X <- if (sim_settings$include.Xintercept)  matrix(1, n, 1) else NULL
      if (p0 > 0) X <- cbind(X,  simulate_IIDdata(n, distr_name = X_distr$distr_name,
                                                  distr_params  = X_distr$distr_params,
                                                  generator     = X_distr$generator))
    }

    # Dimension checks
    beta_coeff <- beta_coeff %||% rep(0, ncol(X))
    if (length(beta_coeff) != ncol(X)) {
      stop(sprintf("Length of beta_coeff (%d) does not match number of X columns (%d)",
                   length(beta_coeff), ncol(X)))
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

    # Name the columns of X if required
    if (is.null(colnames(X))) colnames(X) <- if (sim_settings$include.Xintercept) paste0("X", 0:(ncol(X) - 1)) else paste0("X", 1:ncol(X))

    # linear predictor
    eta <- drop(X %*% beta_coeff)

    # Simulate errors
    if (IID_error){

      local_sigma <- if (SNR_control) sqrt(var(eta) / SNR) else as.numeric(error_sigma[1,1])
      error_params_updated <- do.call(distr_spec, c(list(sd=local_sigma), error_params_WOsd))

      eps <- matrix(simulate_IIDdata(n*it_num, distr_name = error_distr$distr_name,
                                     distr_params  = error_params_updated,
                                     generator     = error_distr$generator),
                  nrow = n, ncol = it_num)
      } else {
        local_sigma <- sqrt(mean(diag(error_sigma)))
        eps <- t(simulate_IIDdata(it_num, distr_name = error_distr$distr_name,
                                  distr_params  = error_distr$distr_params,
                                  generator     = error_distr$generator))
      }


    # Generate responses
    ys_mat <- matrix(NA_real_, nrow = n, ncol = it_num)

    # Loop over n_iter
    for (it in seq_len(it_num)) {
      ys_mat[, it]<- eta + eps[,it]
    }

    # Store Results
    class(X) <- c("IIDdata", class(X))
    Xs[[si]] <- X
    yss[[si]] <- ys_mat
    sigma_list[[si]] <- local_sigma
    SNRs[[si]] <- var(eta)/ (local_sigma^2)
  }

  ## Output
  out <- list(
    y = yss, X = Xs,
    sigma_e = sigma_list, SNR = SNRs
  )

  if (sim_num == 1) out <- lapply(out, `[[`, 1)
  out$n <- n


  sim_settings$p <- p
  sim_settings$beta_coeff <- beta_coeff
  out$sim_settings <- sim_settings

  distr_settings$X_distr <- X_distr
  distr_settings$error_distr <- error_distr
  out$distr_settings = distr_settings

  class(out) <- c("LMdata", class(out))
  return(out)
}

