##===========================================================
#' @title Distribution Specification Constructor
#'
#' @description
#' Construct a unified distribution specification object (list) for use in
#' all simulator functions, e.g., \code{\link{simulate_IIDdata}},
#' \code{\link{simulate_LMdata}}, etc.
#'
#' The function resolves parameter *alias groups* (location, scale, shape),
#' infers the dimension of the distribution, validates
#' compatibility across parameter groups, and returns a fully specified
#' list of distribution parameters. See @details.
#'
#'
#' @param dim Integer (>0). Initial dimension of the random variable.
#'   Default: `1`.
#' @param mean,mu,location,xi Numeric. Aliases for the *location* parameter.
#'   May be scalars or vectors. Default: \code{rep(0, dim)}.
#' @param sd,scale,sigma,omega,Omega,V Numeric or matrix.
#'   Aliases for the *scale* or *covariance* parameter.
#'   Scalars are expanded to \code{scalar \eqn{\times} I(dim)};
#'   matrices must be square.Default: \code{I(dim)}.
#' @param alpha Numeric. Slant / skewness parameter.
#'   May be scalar or vector. Default: \code{rep(0, dim)}.
#' @param lambda Numeric (>=0). Mean parameter for Poisson-type distributions.
#'   Default: \code{1}.
#' @param rate Numeric (>0). Rate parameter for distributions with positive supports
#'   (e.g., exponential, gamma). Default: \code{1}.
#' @param df,nu Numeric (>=0). Degree of freedom for distributions with a
#'   single df parameter (e.g., \eqn{t}, skew-\eqn{t}). Default: \code{5}.
#' @param df1,df2 Numeric (>=0). Degrees of freedom for distributions with two
#'   df parameters (e.g., \eqn{F}). Default: \code{5} for both.
#' @param ncp Numeric (>=0). Non-centrality parameter. Default: \code{0}.
#' @param shape Numeric (>=0). Shape parameter for distributions with one shape parameter
#'   (e.g., gamma, Weibull). Default: \code{1}.
#' @param shape1,shape2 Numeric (>=0). Shape parameters for distributions with
#'   two shape parameters (e.g., beta). Default: \code{1} for both.
#' @param size Integer (>0). Size parameter for binomial-type distributions.
#'   Default: \code{1}.
#' @param prob Numeric in [0,1]. Probability parameter. Default: \code{0.5}.
#' @param meanlog,sdlog Numeric. Mean and standard deviation on the log-scale.
#'   Defaults: \code{0}, \code{1}.
#' @param min,max Numeric. Lower and upper bounds of the support for bounded distributions,
#'   (e.g., uniform). Defaults: \code{0}, \code{1}.
#'
#' @param copula Optional copula object (from the \pkg{copula} package),
#'    necessary for copula-based multivariate simulations.
#' @param margins Optional list of distributional specifications for each margins in copula models
#'    (see \code{\link{simulate_IIDdata}}).
#'
#' @param generator Optional user-defined function for custom data generation.
#'
#' @param ... Additional named arguments. These are appended verbatim to
#' the output list. If names conflict with existing parameters, the
#' values supplied in \code{...} will overwrite the defaults.
#'
#' @return
#' A named list containing:
#' \itemize{
#'   \item Resolved dimension (\code{dim})
#'   \item Resolved aliases (\code{mean}, \code{mu}, \code{location}, \code{xi}, \code{sd}, etc.)
#'   \item All remaining formal parameters
#'   \item Additional user-supplied entries from \code{...}
#' }
#'
#'
#' @details
#' **Alias groups:**
#' The following alias groups are supported (to date) through the internal function \code{.resolve_alias_group()}:
#'
#' \itemize{
#'   \item \strong{Location:} \code{mean}, \code{mu}, \code{location}, \code{xi}
#'   \item \strong{Scale / covariance:} \code{sd}, \code{scale}, \code{sigma},
#'         \code{omega}, \code{Omega}
#'   \item \strong{Shape / skewness:} \code{alpha}
#' }
#'
#' Within each group, the first non-\code{NULL} argument in the order of
#' formal arguments is used. Scalars are automatically expanded.
#' If no alias in a group is supplied, the default is set to
#' `rep(0, dim)` for vectors, and to `diag(dim)` for matrices.
#'
#' **Dimension inference and interaction:**
#' The working dimension is initialized from \code{dim} and may subsequently
#' be updated by:
#' \enumerate{
#'   \item location aliases (vector length),
#'   \item scale aliases (matrix order),
#'   \item alpha aliases (vector length).
#' }
#' If a later group implies a different dimension than earlier groups,
#' strict compatibility checks are enforced using \code{.check_dim_strict()}.
#'
#' When the dimension changes due to a later alias group, dependent
#' parameters are reset to their default values at the new dimension.
#'
#'
#' **Error conditions:** Errors are raised when:
#' \itemize{
#'   \item scale aliases are non-scalar non-square matrices,
#'   \item vector aliases have incompatible lengths,
#'   \item multiple alias groups imply incompatible dimensions,
#'   \item user-supplied dimensions conflict with resolved dimension.
#' }
#'
#'
#' @examples
#' # Default univariate distribution specifications
#' distr_spec()
#'
#' # Multivariate distribution setup (dimension inferred from inputs)
#' distr_spec(mean = c(1, 2, 3), sigma = diag(3))
#'
#' # Specification for Gamma distribution
#' distr_spec(shape = 2, rate = 1)
#'
#' # Specification for Beta distribution
#' distr_spec(shape1 = 2, shape2 = 5)
#'
#' # Custom additional arguments
#' distr_spec(mean = 10, extra1 = "hello", extra2 = 42)
#'
#' @seealso \code{\link{simulate_IIDdata}}
#' @export
distr_spec <- function(
    dim = 1,
    mean = NULL, mu = NULL, location = NULL, xi = NULL,
    sd = NULL, scale = NULL, sigma = NULL, omega = NULL, Omega = NULL, V = NULL,
    alpha = NULL,
    lambda = 1, rate = 1,
    df = 5, nu = 5, df1 = 5, df2 = 5,
    ncp = 0,
    shape = 1, shape1 = 1, shape2 = 1,
    size = 1, prob = 0.5,
    meanlog = 0, sdlog = 1,
    min = 0, max = 1,
    copula = NULL,
    margins = NULL,
    generator = NULL,
    ...
) {

  # 1. LOCATION GROUP (mean / mu / location / xi)
  loc_out <- .resolve_alias_group(
    values  = list(mean, mu, location, xi),
    dim     = dim,
    default = 0,
    type    = "vector"
  )

  loc <- loc_out$value
  dim <- loc_out$dim
  loc_user_dim <- if (loc_out$has_alias) length(loc) else NULL


  # 2. SCALE GROUP (sd / scale / sigma / omega / Omega)
  scale_out <- .resolve_alias_group(
    values  = list(sd, scale, sigma, omega, Omega, V),
    dim     = dim,
    default = 1,
    type    = "matrix"
  )

  scale_mat <- scale_out$value
  scale_user_dim <- if (scale_out$has_alias) scale_out$dim else NULL


  # If scale implies new dimension → STRICT check
  if (scale_out$has_alias && scale_out$dim != dim) {
    new_dim <- scale_out$dim

    .check_dim_strict("Location", loc_user_dim, new_dim)

    dim <- new_dim
    loc <- rep(0, dim)
  }


  # 3. ALPHA GROUP
  alpha_out <- .resolve_alias_group(
    values  = list(alpha),
    dim     = dim,
    default = 0,
    type    = "vector"
  )

  alpha_val <- alpha_out$value
  alpha_user_dim <- if (alpha_out$has_alias) alpha_out$dim else NULL


  # If alpha determines dimension → STRICT check
  if (alpha_out$has_alias && alpha_out$dim != dim) {
    new_dim <- alpha_out$dim

    .check_dim_strict("Location", loc_user_dim, new_dim)
    .check_dim_strict("Scale",    scale_user_dim, new_dim)

    dim <- new_dim
    loc <- rep(0, dim)
    scale_mat <- diag(dim)
  }

  # scale_mat <- if(dim == 1) drop(scale_mat) else scale_mat

  # Final Output
  out <- list(
    dim = dim,
    mean = loc, mu = loc, location = loc, xi = loc,
    sd = scale_mat, scale = scale_mat, sigma = scale_mat,
    omega = scale_mat, Omega = scale_mat, V = scale_mat,
    alpha = alpha_val,
    lambda = lambda, rate = rate,
    df = df, nu = nu, df1 = df1, df2 = df2,
    ncp = ncp,
    shape = shape, shape1 = shape1, shape2 = shape2,
    size = size, prob = prob,
    meanlog = meanlog, sdlog = sdlog,
    min = min, max = max,
    copula = copula,
    margins = margins,
    generator = generator
  )

  c(out, list(...))
}

#' Alias for distr_spec
#'
#' @inherit distr_spec
#' @rdname distr_spec
#' @export
distrSpec <- distr_spec



##===========================================================
#' @title Simulation Specification Constructor
#'
#' @description
#' Constructs a unified simulation specification object (list) for generating data
#' from several structured statistical models, e.g., linear models (LM),
#' generalized linear models (GLM), linear mixed models (LMM),
#' generalized LMMs (GLMM), or other regression-type models.
#' The output may directly be passed to the corresponding simulator functions,
#' such as \code{\link{simulate_LMdata}}, \code{\link{simulate_LMMdata}}, etc.
#'
#'
#' @param name Character. Label for the simulation setup. Default: \code{"Setup0"}.
#' @param n_subj Integer (>0). Number of subjects or clusters. Default: `10`.
#' @param n_obs Integer scalar or integer vector (>0).
#'   Specifies the number of observations per subject/cluster.
#'   - If a **scalar** is supplied, all `n_subj` clusters are assigned the same number of observations.
#'   - If a **vector**, it must be of length `n_subj`, with each entry giving the size of the respective cluster.
#'   Useful for longitudinal or clustered designs.  Default: `5` for all clusters.
#'
#' @param p Integer (>0). Number of fixed-effect predictors (including intercept).
#'          If \code{NULL}, it is inferred as \code{length(beta_coeff)} whenever possible.
#'          Default fallback value: `1`
#'
#' @param q Integer (>0). Number of random-effect predictors. Default: `1`.
#'
#' @param beta_coeff Numeric vector. Fixed-effect coefficients of length `p`.
#' If supplied, its length is used to infer \code{p} when \code{p = NULL}.
#'
#' @param SNR Numeric (>0). Signal-to-noise ratio.
#'
#' @param X Numeric matrix or list. User-supplied design matrix for fixed effects.
#'          If provided, it is stored verbatim.
#'
#' @param Z Numeric matrix or list. User-supplied design matrix for random effects.
#'          If provided, it is stored verbatim.
#'
#' @param is.CorrelatedXZ Logical. Whether `X` and `Z` are correlated
#'          (and hence they need to be generated jointly). Default: `FALSE`.
#' @param is.ZInX Logical. Whether the columns of `Z` are the same as
#'          the initial columns of `X`. Default: `FALSE`.
#'
#' @param orthogonalize.X Logical. Whether orthogonalize (non-intercept) columns
#'          of `X` after construction. Default: `FALSE`.
#' @param orthogonalize.Z Logical. Whether orthogonalize (non-intercept) columns
#'          of `Z` after construction. Default: `FALSE`.
#'
#' @param include.Xintercept Logical. Whether an intercept column should be included in `X`. Default: `TRUE`.
#' @param include.Zintercept Logical. Whether an intercept column should be included in `Z`. Default: `TRUE`.
#'
#' @param ... Additional named arguments to append to the returned specification list.
#'
#'
#' #' @return
#' A named list containing:
#' \itemize{
#'   \item \code{name} — Setup label.
#'   \item \code{n_subj}, \code{n_obs} — Sample size arguments.
#'   \item \code{p}, \code{q} — Predictor dimensions (as supplied or inferred).
#'   \item \code{beta_coeff} — Fixed-effect coefficients.
#'   \item \code{SNR} - Desired signal-to-noise ratio.
#'   \item \code{X}, \code{Z} — Optional user-supplied design matrices.
#'   \item \code{is.CorrelatedXZ}, \code{is.ZInX} — Logical flags,
#'          specifying desired relationships between the design matrices `X` and `Z` .
#'   \item \code{orthogonalize.X}, \code{orthogonalize.Z} — Orthogonalization flags.
#'   \item \code{include.Xintercept}, \code{include.Zintercept} — Intercept controls.
#'   \item Any additional entries from \code{...}.
#' }
#'
#' @examples
#' # Basic specification
#' spec1 <- sim_spec()
#'
#' # Custom fixed effects
#' spec2 <- sim_spec(beta_coeff = c(1, -1, 0.5))
#'
#' # Mixed-effects style specification
#' spec3 <- sim_spec(
#'   n_subj = 15, n_obs = 6,
#'   beta_coeff = c(1, 0.5),
#' )
#'
#' @seealso
#' \code{\link{simulate_LMdata}},
#' \code{\link{simulate_LMMdata}}
#'
#' @export
sim_spec <- function(
    name = "Setup0",
    n_subj = 10, n_obs = 5,
    p = NULL, q = 1,
    beta_coeff = NULL,
    SNR = NULL,
    X = NULL, Z = NULL,
    is.CorrelatedXZ = FALSE,
    is.ZInX = FALSE,
    orthogonalize.X = FALSE,
    orthogonalize.Z = FALSE,
    include.Xintercept = TRUE,
    include.Zintercept = TRUE,
    ...
    ) {

  # Default for `p`
  # ---------------------------
  p <- p %||% length(beta_coeff)
  if (p == 0) p = 1

  # Output
  c(list(
    name = name,
    n_subj = n_subj, n_obs = n_obs,
    p = p, q = q,
    beta_coeff = beta_coeff,
    SNR = SNR,
    X = X, Z = Z,
    is.CorrelatedXZ = is.CorrelatedXZ ,
    is.ZInX = is.ZInX,
    orthogonalize.X = orthogonalize.X,
    orthogonalize.Z = orthogonalize.Z,
    include.Xintercept = include.Xintercept,
    include.Zintercept = include.Zintercept
  ), list(...))
}


#' @rdname sim_spec
#' @export
simSpec <- sim_spec



##=========================================================================================
#' @title Find a Function Across Base R and All Installed R Packages
#'
#' @description
#' Extends the behavior of \code{\link[base]{match.fun}} by searching not only
#' the current environment and base R, but also all installed packages
#' (namespaces). It additionally supports the explicit \code{"pkg::fun"} syntax
#' to retrieve a function from a specific package namespace.
#'
#' @param FUN A function,  character string, or symbol naming a function.
#'   May also be of the form \code{"pkg::fun"} to explicitly reference a
#'   function within a particular package.
#'
#' @param descend Logical; if \code{TRUE} (default), the lookup in the current
#'   environment behaves like \code{\link[base]{match.fun}} and descends
#'   the inheritance tree when resolving S3/S4 methods. If \code{FALSE},
#'   only non-generic function objects are accepted.
#'
#' @param list_all Logical; if \code{TRUE}, returns a list of all matching
#'   function objects found across installed R packages. If \code{FALSE}
#'   (default), returns only the first match found. See @details
#'
#' @return A function (if a single match is found or `list_all = FALSE`)
#'   or a named list of functions (if multiple matches are found and `list_all = TRUE`).
#'   If a list is returned, the names correspond to the package
#'   or environment where each function is found.
#'   An informative error is raised if the function cannot be found.
#'
#'
#' @details
#' It searches for \code{FUN} in the following order:
#' \enumerate{
#'   \item If \code{FUN} is already a function, it is returned unchanged.
#'   \item If \code{FUN} is of the form \code{"pkg::fun"},
#'         then \code{pkg} must be installed;
#'         the function is retrieved directly from \code{asNamespace(pkg)}.
#'   \item The current calling environment is searched (similar to
#'         \code{\link[base]{match.fun}}). If \code{descend = TRUE},
#'         method dispatch is allowed; otherwise only plain functions are accepted.
#'
#'   \item All installed package namespaces are searched (both exports and internal functions).
#'         The first match found—according to the alphabetical order
#'         of installed packages—is returned if \code{list_all = FALSE}.
#'         If   \code{list_all = TRUE}, all matches are returned.
#'         *Note:* This may yield surprising matches when multiple packages
#'         define functions with identical names and \code{list_all = FALSE};
#'         in such cases, the explicit \code{"pkg::fun"} syntax is recommended.
#' }
#'
#' @examples
#' # Base R function
#' f1 <- match.fun.allR("rnorm")
#' f1(5)
#'
#' # Function from installed package (if sn is installed)
#' if (requireNamespace("sn", quietly = TRUE)) {
#'   f2 <- match.fun.allR("rsn")       # search all installed packages
#'   f3 <- match.fun.allR("sn::rsn")   # explicit package reference
#'   f2(5); f3(5)
#' }
#'
#' f4 <- match.fun.allR("rnorm", list_all = TRUE)
#' class(f4)
#'
#' @seealso
#' \code{\link[base]{match.fun}} for the base R version,
#' \code{getNamespace}, \code{getExportedValue}
#'
#' @export
match.fun.allR <- function(FUN, descend = TRUE, list_all = FALSE) {

  # If already a function, return it
  if (is.function(FUN)) return(FUN)

  # If not character or symbol, try evaluating
  if (!(is.character(FUN) && length(FUN) == 1L || is.symbol(FUN))) {
    FUN <- eval.parent(substitute(substitute(FUN)))
    if (!is.symbol(FUN)) stop(sprintf("'%s' is not a function, character or symbol", deparse(FUN)))
  }

  fname <- as.character(FUN)

  # Check specified package if 'fname' is of the form pkg::fun
  if (grepl("::", fname, fixed = TRUE)) {
    parts <- strsplit(fname, "::", fixed = TRUE)[[1]]
    pkg <- parts[1]; f <- parts[2]
    if (!requireNamespace(pkg, quietly = TRUE)) stop("Package '", pkg, "' is required but not installed.")
    return(get(f, envir = asNamespace(pkg)))
  }

  # Check the current environment (like match.fun)
  envir <- parent.frame(2)
  if (descend) f <- tryCatch(get(fname, mode = "function", envir = envir), error = function(e) NULL)
  else {
    f <- tryCatch(get(fname, mode = "any", envir = envir), error = function(e) NULL)
    if (!is.function(f)) stop(sprintf("found non-function '%s'", fname))
  }
  if (!is.null(f)) return(f)

  # Search across all installed R packages
  pkgs <- installed.packages()[, "Package"]
  matches <- list()

  for (pkg in pkgs) {
    ns <- tryCatch(getNamespace(pkg), error = function(e) NULL)
    if (!is.null(ns) && exists(fname, envir = ns, mode = "function", inherits = FALSE)) {
      f_obj <- get(fname, envir = ns)
      if (list_all) matches[[pkg]] <- f_obj
      else return(f_obj)  # Return first match immediately
    }
  }

  if (length(matches) > 1) return(matches)
  if (length(matches) == 1) return(matches[[1]])

  # Not found
  stop(sprintf("Object '%s' of mode 'function' is not found in Base R or any installed packages.", fname))
}

