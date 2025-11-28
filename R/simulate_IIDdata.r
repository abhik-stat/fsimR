#' @title Simulate IID Data from a Wide Range of Distributions
#'
#' @description
#' Generates independent and identically distributed (IID) random samples
#' with flexibility and ease of control.
#'
#' The function supports standard R random generators (e.g., `rnorm`, `rpois`),
#' external generator functions from *all installed R packages*
#' (e.g., `\pkg{sn}`, `\pkg{mvtnorm}`, `\pkg{VGAM}`,`\pkg{extraDistr}`, `\pkg{gamlss.dist}`, `\pkg{gamlss.dist}`),
#' copula-based multivariate distributions (via `\pkg{copula}`)
#' as well as user-defined custom generators (supplied via the `generator` argument).
#' User-supplied parameters override defaults produced by [distr_spec()].
#'
#'
#' @param n Integer (>0). Number of IID observations to simulate.
#'          Fractional values are rounded down to the nearest integer.
#'
#' @param distr_name Character string specifying the target distribution.
#'   See \link[=simulate_IIDdata]{Details} for full information on aliases,
#'   parameter handling, and examples.  Supported options (case-insensitive) include:
#'   \itemize{
#'     \item **Base R distributions:** `"normal"`, `"gaussian"`, `"binomial"`, `"poisson"`,
#'           `"uniform"`, `"gamma"`, `"beta"`, `"chisq"`, etc.
#'           See the full list of base R distributions at:
#'           \url{https://stat.ethz.ch/R-manual/R-devel/library/stats/html/Distributions.html}
#'     \item **Multivariate distributions:** `"mvnorm"`, `"mvt"`, `"mvst"`, `"mvsnorm"`.
#'     \item **Copula-based generation:** \code{"copula"} → see @details for usage.
#'     \item **Custom generator:**  \code{"custom"} → requires providing \code{generator}.
#'     \item **Explicit package-qualified functions:**  any \code{"pkg::fun"} form
#'            (e.g., \code{"MASS::mvrnorm"})
#'     \item **External functions without package prefix:**
#'           the simulator attempts to locate a matching function across
#'           all installed packages via \code{\link{match.fun.allR}}.
#'   }
#'
#' @param distr_params Named list of distribution parameters,
#'   usually created with [distr_spec()].
#'   Can include an optional argument `dim` representing the dimension of the random variables.
#'   For multivariate distributions, `dim` must either be supplied explicitly
#'   or be inferable from the user-supplied location, scale, or shape parameters
#'   (or from the 'copula' object if `distr_name = "copula"`).
#'   If `distr_params` is not a list, it is interpreted as:
#'   \itemize{
#'     \item The **copula object** if `distr_name = "copula"`.
#'     \item The **dimension** (`dim`) if `distr_name` is neither "copula" nor "custom".
#'   }
#'
#' @param generator A user-defined random generator function, used only when
#'   \code{distr_name = "custom"}. Must accept at least the
#'   argument \code{n}, and return either a vector of length \code{n}
#'   or a matrix of dimension \code{n × dim}).
#'
#' @param seed Optional integer. If provided, the RNG state is temporarily
#'   replaced with this seed, and restored upon exit.
#'
#'
#' @return
#' A numeric matrix (with class `"IIDdata"`) with `n` rows and `dim` columns.
#' Each row is an IID observation of dimension `dim`.
#'
#' @details
#' \itemize{
#'   \item **Copula-based simulations:** `distr_params` must contain a `copula` object (from \pkg{copula})
#'         and an optional list `margins` of length `dim` (>1).
#'         Each element of `margins`, if provided, must be a list containing
#'           - `dist`: univariate distribution name (e.g. `"norm"`),
#'           - `params`: named list of distribution parameters.
#'         The simulator first draws \eqn{U}  from the specified copula \eqn{C}
#'         and applies the appropriate quantile functions to obtain the final multivariate sample.
#'         *Defaults:* If `dim` is missing, it is inferred from the copula object;
#'         if `margins` is missing, all margins are taken to be standard normal.
#'
#'         An informative error is raised if there is a dimension mismatch between the copula
#'         and the user-supplied arguments.
#'
#'   \item **Custom generators:** Only arguments matching \code{generator}'s formal arguments are passed.
#'   Sample size is forwarded as \code{n}.
#'
#'   \item **External generators with specified r package:**
#'   If `distr_name` is provided in the form `"pkg::function"`,
#'   the simulator dynamically loads the specified function from the indicated package.
#'   An error is raised if the package is not installed.
#'
#'   \item **External generators without package prefix:**
#'   If `distr_name` is **not** of the form `"pkg::function"`,
#'   the function first looks for a base R random generator (e.g., `rnorm`, `rpois`).
#'   If none is found, it searches all installed packages for a matching exported function,
#'   prepending `r` if needed (e.g., `norm` → `rnorm`).
#'   An informative error is raised if no suitable function is found.
#'   **Warning:** naming conflicts across packages may lead to unexpected behavior.
#'
#'   \item **Sample-size detection:**
#'         The function automatically detects `n` or `nn` within formal arguments.
#'         If neither exists, `n` is passed as the first positional argument
#'         (specifying the sample size).
#' }
#'
#' **Distribution name aliases:**
#' The argument `distr_name` supports the following (case-insensitive) aliases for common distributions.
#'    \tabular{llll}{
#' *Distribution*             \tab Aliases                    → Internal function         \tab `distr_params` (with defaults)\cr
#' Binomial                   \tab binomial, binom            → [rbinom]             \tab size = 1, prob = 0.5\cr
#' Poisson                    \tab poisson, pois              → [rpois]              \tab lambda = 1\cr
#' Normal                     \tab normal, gaussian, norm     → [rnorm]              \tab mean = 0, sd = 1\cr
#' Multivariate Normal        \tab mvnorm, mvn                → [mvtnorm::rmvnorm]   \tab mean = rep(0, dim), sigma = diag(dim)\cr
#' Multivariate Normal (Singular)   \tab smvnorm, smvn        → [mgcv::rmvn]         \tab mu = rep(0, dim), V = diag(dim)\cr
#' Multivariate t             \tab mvt                        → [mvtnorm::rmvt]      \tab sigma = diag(dim), df = 5, delta = rep(0, dim)\cr
#' Multivariate t (Singular)  \tab smvt                       → \link[mgcv]{r.mvnt}  \tab mu = rep(0, dim), V = diag(dim), df = 5\cr
#' Skew-normal (univariate)   \tab skewnorm, sn               → [sn::rsn]            \tab xi = 0, omega = 1, alpha = 0, tau = 0\cr
#' Skew-normal (multivariate) \tab mvsnorm, mvsn              → [sn::rmsn]           \tab xi = alpha = rep(0, dim), Omega = diag(dim), tau = 0 \cr
#' Skew-t (univariate)        \tab skewt, st                  → [sn::rst]            \tab xi = 0, omega = 1, alpha = 0, nu = 5\cr
#' Skew-t (multivariate)      \tab mvskewt, mvst              → [sn::rmst]           \tab xi = alpha = rep(0, dim), Omega = diag(dim), nu = 5\cr
#' Pareto (type I)            \tab pareto, pareto_I           → [VGAM::rpareto]      \tab scale = 1, shape = 1 (k parameter)  \cr
#' Laplace                    \tab laplace, double_exp        → [VGAM::rlaplace]     \tab location = 0, scale = 1  \cr
#'          }
#' Functions from external packages (e.g., `sn::rmsn`) may also be specified
#' without the package prefix  if no naming conflict exists,
#' though this can lead to ambiguity or unexpected behavior.
#' Omitting the leading `r` is allowed for internal R function names
#' when the package name is not explicitly specified using `::`.
#'
#' See the \href{https://CRAN.R-project.org/view=Distributions}{CRAN Task View on Distributions}
#' for a comprehensive list of R packages providing additional distributions,
#' all of which can be used within this simulator with suitable `distr_name` and `distr_params`.
#'
#' @section Dependencies:
#' Depending on the distribution specified in `distr_name`,
#' one or more of the following packages may be optionally required
#' (depending on the case):
#' \itemize{
#'   \item \pkg{copula} – For copula-based generation (\code{distr_name = "copula"}).
#'         CRAN: \url{https://CRAN.R-project.org/package=copula}
#'   \item \pkg{mvtnorm} – For multivariate normal (\code{"mvnorm", "mvt"}) and t distributions.
#'         CRAN: \url{https://CRAN.R-project.org/package=mvtnorm}
#'   \item \pkg{sn} – For skew-normal and skew-t distributions (both univariate and multivariate:
#'         \code{"skewnorm", "skewt", "mvsnorm", "mvskewt"}).
#'         CRAN: \url{https://CRAN.R-project.org/package=sn}
#'   \item \pkg{VGAM} – For non-standard univariate distributions such as Pareto, Laplace, and others:
#'         \code{"pareto", "laplace"}.
#'         CRAN: \url{https://CRAN.R-project.org/package=VGAM}
#'   \item \pkg{extraDistr} – Provides additional discrete and continuous distributions
#'         (e.g., generalized beta, zero-inflated distributions).
#'         CRAN: \url{https://CRAN.R-project.org/package=extraDistr}
#'   \item \pkg{MASS} – For distributions like negative binomial (\code{"MASS::rnegbin"}).
#'         CRAN: \url{https://CRAN.R-project.org/package=MASS}
#'   \item Any other R package \code{"pkg"} that provides a function specified
#'          via \code{"pkg::fun"} in \code{distr_name}.
#'          The function will attempt to dynamically load the specified function from the package.
#' }
#'
#' @references
#' \itemize{
#'   \item Azzalini A, Capitanio A (2013). \emph{The Skew-Normal and Related Families} Cambridge University Press. ISBN 978-1-107-02927-9. \url{https://doi.org/10.1017/CBO9781139248891}
#'   \item Adcock C, Dias A, Salmon M (Eds) (2013). \emph{Copulae & Multivariate Probability Distributions in Finance} Routledge. ISBN 978-0-415-81485-0. \url{https://www.routledge.com/Copulae-and-Multivariate-Probability-Distributions-in-Finance/Adcock-Dias-Salmon/p/book/9780415814850}
#'   \item Durante F, Sempi C (2015). \emph{Principles of Copula Theory} CRC Press. ISBN 978-1-439-88442-3. \url{https://www.crcpress.com/Principles-of-Copula-Theory/Durante-Sempi/p/book/9781439884423}
#'   \item Genton MG (Ed) (2004). \emph{Skew-Elliptical Distributions and Their Applications: A Journey Beyond Normality} Chapman & Hall / CRC Press. ISBN 978-0-415-36559-7. \url{https://www.routledge.com/Skew-Elliptical-Distributions-and-Their-Applications-A-Journey-Beyond-Normality/Genton/p/book/9780367578312}
#'   \item Johnson NL, Kotz S, Balakrishnan N (1994). \emph{Continuous Univariate Distributions, Volume 1} (2nd ed.) Wiley. ISBN 978-0-471-58495-7. \url{https://www.wiley.com/en-us/Continuous+Univariate+Distributions\%2C+Volume+1\%2C+2nd+Edition-p-9780471584957}
#'   \item Johnson NL, Kotz S, Balakrishnan N (1995). \emph{Continuous Univariate Distributions, Volume 2} (2nd ed.) Wiley. ISBN 978-0-471-58494-0. \url{https://www.wiley.com/en-us/Continuous+Univariate+Distributions\%2C+Volume+2\%2C+2nd+Edition-p-9780471584940}
#'   \item Johnson NL, Kotz S, Balakrishnan N (1997). \emph{Discrete Multivariate Distributions} Wiley. ISBN 978-0-471-12844-1. \url{https://www.wiley.com/en-us/Discrete+Multivariate+Distributions-p-9780471128441}
#'   \item Johnson NL, Kotz S, Kemp AW (2005). \emph{Univariate Discrete Distributions} (3rd ed.) Wiley. \url{https://www.wiley.com/en-us/Univariate+Discrete+Distributions\%2C+3rd+Edition-p-9780471712384}
#'   \item Kotz S, Balakrishnan N, Johnson NL (2000). \emph{Continuous Multivariate Distributions, Volume 1: Models & Applications} (2nd ed.) Wiley. ISBN 978-0-471-18387-7. \url{https://www.wiley.com/en-us/Continuous+Multivariate+Distributions\%2C+Volume+1\%3A+Models+and+Applications-p-9780471183877}
#'   \item Lai CD, Balakrishnan N (2002). \emph{Continuous Bivariate Distributions} Springer. \url{https://link.springer.com/book/10.1007/b101765}
#'   \item Nelsen RB (2006). \emph{An Introduction to Copulas} (2nd ed.) Springer. \url{https://link.springer.com/book/10.1007/0-387-28678-0}
#' }
#'
#'
#' @examples
#' # A Basic Example
#' set.seed(123)
#' x <- simulate_IIDdata(500, "laplace")
#' mean(x); var(x); hist(x)
#'
#' # From a recommended package (MASS)
#' x <- simulate_IIDdata(5, "MASS::rnegbin", list(mu = 5, theta = 2))
#' mean(x); var(x); hist(x)
#'
#' # Custom generator
#' my_custom <- function(n, min = 1, max = 100) sample(min:max, n, replace = TRUE)
#' x <- simulate_IIDdata(500, "custom", distr_params = list(min = 1, max = 10), generator = my_custom)
#' mean(x); var(x); plot(x)
#'
#' # Multivariate example
#' params <- distr_spec(dim = 3, mean = 1:3, sigma = diag(3))
#' x <- simulate_IIDdata(5, "mvnorm", params)
#' plot(x)
#'
#' # Skew-normal (requires \pkg{sn})
#' if (requireNamespace("sn", quietly = TRUE)) {
#'   x <- simulate_IIDdata(50, "skewnorm")
#'   plot(x)
#' }
#'
#' # Copula-based Examples
#' if (requireNamespace("copula", quietly = TRUE)) {
#'   library(copula)
#'
#'   # Gaussian copula with Normal margins
#'   cop <- normalCopula(param = 0.5, dim = 3)
#'   margins <- list(
#'     list(dist = "norm", params = list(mean = 0, sd = 1)),
#'     list(dist = "norm", params = list(mean = 1, sd = 2)),
#'     list(dist = "norm", params = list(mean = -1, sd = 0.5))
#'   )
#'   x<-simulate_IIDdata(5, "copula", list(dim = 3, copula = cop, margins = margins))
#'   plot(x)
#'
#'   # Gumbel copula with mixed margins
#'   gumbel_cop <- gumbelCopula(param = 2, dim = 4)
#'   mixed_margins <- list(
#'     list(dist = "norm", params = list(mean = 0, sd = 1)),
#'     list(dist = "gamma", params = list(shape = 3, rate = 1)),
#'     list(dist = "exp", params = list(rate = 2)),
#'     list(dist = "unif", params = list(min = 0, max = 1))
#'   )
#'   x <- simulate_IIDdata(10, "copula", list(dim = 4, copula = gumbel_cop, margins = mixed_margins))
#'   plot(x)
#'
#'   # Large-dimensional Gaussian copula with exponential margins
#'   d <- 20
#'   exp_margins <- replicate(d, list(dist = "exp", params = list(rate = 1)), simplify = FALSE)
#'   rho <- 0.3
#'   Sigma <- matrix(rho, nrow = d, ncol = d); diag(Sigma) <- 1
#'   norm_cop <- normalCopula(P2p(Sigma), dim = d, dispstr = "un")
#'   x_large <- simulate_IIDdata(10, "copula", list(dim = d, copula = norm_cop, margins = exp_margins))
#'   dim(x_large)  # should be 10 x 20
#'   plot(x)
#' }
#'
#' @seealso [distr_spec()], \code{\link{match.fun.allR}}, [summary.IIDdata], [plot.IIDdata]
#'
#' @importFrom MASS rnegbin
#' @importFrom stats setNames rbinom rpois rnorm runif rgamma rbeta rchisq
#' @importFrom copula rCopula
#' @importFrom sn rsn rst rmsn rmst
#' @importFrom mvtnorm rmvnorm rmvt
#' @importFrom VGAM rpareto rlaplace
#'
#' @export
simulate_IIDdata <- function(n, distr_name, distr_params = distr_spec(),
                             generator = NULL, seed = NULL) {

  ## Save RNG state when seed is provided
  if (!is.null(seed)) {
    old_seed <- .GlobalEnv$.Random.seed %||% NULL
    set.seed(seed)
    on.exit({
      if (!is.null(old_seed)) .GlobalEnv$.Random.seed <- old_seed
      else rm(".Random.seed", envir = .GlobalEnv)
      }, add = TRUE)
  }


  ## Input validations
  # Input validation check for n: n will be rounded below if it is a positive fraction
  if(!is.numeric(n) || length(n)!=1 || n <= 0) stop("'n' must be a positive integer.")
  n <- as.integer(n)

  # Input validation check for distribution name: rename 'distr_name' suitably
  if (!grepl("::", distr_name, fixed = TRUE)) {
    distr_name <- tolower(distr_name)
    distr_name <- switch(distr_name,
                       gaussian = "rnorm",  normal = "rnorm",
                       binomial = "rbinom",
                       poisson  = "rpois",
                       skewnorm = "sn::rsn",  sn = "sn::rsn",
                       skewt    = "sn::rst",  st = "sn::rst",
                       mvnorm   = "mvtnorm::rmvnorm", mvn = "mvtnorm::rmvnorm",
                       mvt      = "mvtnorm::rmvt",
                       mvsnorm  = "sn::rmsn",  mvsn = "sn::rmsn",
                       mvskewt  = "sn::rmst",  mvst = "sn::rmst",
                       pareto   = "VGAM::rpareto",  pareto_I = "VGAM::rpareto",
                       laplace  = "VGAM::rlaplace",  double_exp = "VGAM::rlaplace",
                       distr_name)
    }


  ## Custom generator based simulation: requires specified 'generator' arguments in distr_paranms
  if (distr_name == "custom") {
    if (is.null(generator)) stop("Custom generator function must be provided for distr_name = 'custom'")
    valid_args <- names(formals(generator)) %||% character(0)
    params_filtered <- distr_params[names(distr_params) %in% valid_args]
    out <- as.matrix(do.call(generator, c(list(n = n), params_filtered)))
    class(out) <- c("IIDdata", class(out))
    return(out)
  }



  ## Copula based simulation: requires package 'copula' and data dimension d > 1
  if (distr_name == "copula") {

    if(!is.list(distr_params)) distr_params <- list(copula=distr_params)

    if (is.null(distr_params$copula))
      stop("`copula` must be provided in `distr_params` for copula based simulations.")
    if(!requireNamespace("copula", quietly = TRUE)) stop("Package `copula` is not installed.")

    d <- distr_params$dim %||% 1
    if (d == 1) d <- dim(distr_params$copula)
    if (dim(distr_params$copula) != d) stop("Dimension of `copula` must match `dim`.")

    if (is.null(distr_params$margins)) distr_params$margins <- rep(list(list(dist = "norm")), d)
    if (length(distr_params$margins) != d) stop("Length of `margins` must match `dim`.")

    U <- copula::rCopula(n, distr_params$copula)
    margins <- distr_params$margins

    # Precompute quantile functions
    qfun_list <- lapply(margins, function(margin) {
      if (is.null(margin$dist)) stop("Each margin must have a `dist`")
      match.fun.allR(if (grepl("::", margin$dist, fixed = TRUE)) margin$dist else paste0("q", margin$dist))
    })
    params_list <- lapply(margins, function(m) m$params %||% list())
    # simulate data
    out <- vapply(seq_len(d), function(i) {
                                          do.call(qfun_list[[i]], c(list(p = U[, i]), params_list[[i]]))
                                          }, numeric(n))
    out <- as.matrix(out)
    class(out) <- c("IIDdata", class(out))
    return(out)
  }


  ## For simulations based on R functions (Base R or from any installed packages)
  # Merge user-supplied distribution parameters with default specifications
  if(!is.list(distr_params) && !is.null(distr_params))
    distr_params=list(dim = distr_params)
  distr_params <- do.call(distr_spec, distr_params)
  d = distr_params$dim %||% 1

  # identify function name and load that function
  rfun_name <- if (grepl("::", distr_name, fixed = TRUE) || startsWith(distr_name, "r")) distr_name else paste0("r", distr_name)
  rfun <- match.fun.allR(rfun_name)

  # Extract valid arguments
  valid_args <- names(formals(rfun)) %||% character(0)
  args <- distr_params[names(distr_params) %in% valid_args]

  # Detect sample size argument as 'n' or 'nn'
  sample_arg <- intersect(c("nn", "n"), valid_args)
  if (length(sample_arg) > 0) {
    args[[sample_arg[1]]] <- n
  } else if (length(valid_args) > 0) {
    warning("Sample size argument `n` not found in function; passing `n` as first argument")
    args[[valid_args[1]]] <- n
  }

  # Call function and return as matrix
  # return(as.matrix(do.call(rfun, args)))
  out <- do.call(rfun, args)
  class(out) <- c("IIDdata", class(out))
  return(out)

}

