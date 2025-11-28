# Simulate Data from a Linear Mixed Model (LMM)

Generates datasets from a Linear Mixed Model (LMM) with flexible control
over fixed effects, random effects, correlation structures, design
matrices, and error distributions.

All stochastic components (random effects, error terms, and design
matrices when unspecified) are generated via
[`simulate_IIDdata`](https://abhik-stat.github.io/fsimR/reference/simulate_IIDdata.md),
a **highly general and flexible sample generator** supporting almost all
univariate and multivariate distributions.

The function supports general, possibly correlated and/or
heteroscedastic, error distributions. When the error distribution *is*
independent and identically distributed (**IID**, meaning all error
terms share the same distribution and are mutually independent), users
may specify a target signal-to-noise ratio (SNR). In this case, the
error variance is automatically scaled to match the desired SNR. SNR
control is not supported for non-IID errors.

## Usage

``` r
simulate_LMMdata(
  n_rep = 1,
  sim_settings = list(),
  distr_settings = list(),
  seed = NULL
)
```

## Arguments

- n_rep:

  Integer scalar or length-2 vector. If a vector of length 2, the first
  value gives the number of independent simulated datasets (`n_sim`) and
  the second gives the number of within-dataset iterations (`n_iter`).
  Within a dataset, the design matrices `X` and `Z` remain fixed but
  random effects (`RE`) and error terms vary across iterations,
  producing different response vectors (`y`). If a single value is
  given, it is interpreted as `n_sim`, with `n_iter = 1` (default).

- sim_settings:

  A named list of simulation parameters, usually constructed via
  [`sim_spec`](https://abhik-stat.github.io/fsimR/reference/sim_spec.md)
  (default). Common fields (all optional) for the mixed models include:

  n_subj

  :   Number of subjects or clusters.

  n_obs

  :   Integer scalar or vector (\>0). Specifies the number of
      observations per subject/cluster.

      - If a **scalar** is supplied, all `n_subj` clusters have the same
        number of observations.

      - If a **vector**, its length must be `n_subj`, and each entry
        gives the size of the respective cluster. Default: `5`
        observations per cluster/subject.

  beta_coeff

  :   Fixed-effect coefficient vector. Default: a zero vector of
      appropriate length.

  SNR

  :   Target signal-to-noise ratio for IID errors. Overrides the error
      variance in `distr_settings$error_distr`. Set to `NULL` (default)
      to disable SNR control.

  X, Z

  :   Optional pre-specified design matrices. May also be lists of
      length `n_sim` to allow dataset-specific designs. If supplied,
      they override stochastic generation.

  p, q

  :   Dimensions of fixed- and random-effect design matrices (including
      intercept).

  is.CorrelatedXZ

  :   Logical; If `TRUE`, `X` and `Z` are generated *jointly* from
      `distr_settings$XZ_distr`, and any user-supplied `X` or `Z` is
      ignored.

  is.ZInX

  :   Logical. If `TRUE`, and `Z` is not supplied, the first `q` columns
      of `X` (after intercept handling) are used to construct `Z`.

  include.Xintercept, include.Zintercept

  :   Logical; whether to add intercept columns in `X`/`Z`.

  orthogonalize.X, orthogonalize.Z

  :   Logical; If `TRUE`, the non-intercept columns of `X` and/or `Z`
      are orthogonalized using QR decomposition.

- distr_settings:

  A named list of distribution specifications for stochastic components.
  Valid entries are: `error_distr`, `RE_distr`, `X_distr`, `Z_distr`,
  and `XZ_distr`, specifying the distributions used to generate errors,
  random effects, and design matrices, respectively. Each entry is
  itself a list with elements (see
  [`simulate_IIDdata`](https://abhik-stat.github.io/fsimR/reference/simulate_IIDdata.md)):

  distr_name

  :   Distribution name (e.g., `"norm"`, `"mvnorm"`, `"copula"`).

  distr_params

  :   List of parameters for the specified distribution.

  generator

  :   Optional user-supplied function for random generation, enabling
      simulation from arbitrary or fully custom distributions.

  Missing or incomplete entries are automatically completed via
  [`distr_spec`](https://abhik-stat.github.io/fsimR/reference/distr_spec.md)
  with suitable dimension defaults. The default distribution is
  multivariate normal with zero mean and identity covariance. See
  @details for more on specification of the error distribution.

- seed:

  Optional integer. If provided, the random number generator (RNG) state
  is temporarily replaced by this seed and restored upon exit.

## Value

A named list (of class `"LMMdata"`) containing:

- n:

  Total number of observations in each simulated dataset.

- y:

  List of simulated response matrices of length `n_sim`. Each matrix has
  dimension `n × n_iter`.

- X, Z:

  List of design matrices of length `n_sim`, one per simulated dataset.

- RE:

  List of random effect realizations of length `n_sim`. Each entry is
  either a matrix (when `n_iter = 1`) or a list of matrices (one per
  iteration). These matrices are of order `n_subj × q`.

- sigma_e:

  List of error standard deviations (length `n_sim`). Each entry is a
  vector of length `n_iter`.

- SNR:

  List of achieved signal-to-noise ratios (length `n_sim`). Each entry
  is a vector of length `n_iter`.

- ID:

  Subject identifiers. A factor of length equal to the total number of
  observations, with `n_subj` levels, that aligns with the ordering of
  rows in `y`, `X`, and `Z`. Each level corresponds to one
  subject/cluster, and repeated occurrences of a level indicate multiple
  observations for that subject/cluster.

- sim_settings, distr_settings:

  Lists of (possibly updated) simulation and distribution settings used
  to generate the data.

If a list has only one single element (i.e., `n_sim = 1` or
`n_iter = 1`), it is automatically flattened.

## Details

**Construction of design matrices.** The design matrices `X` and `Z` are
constructed in the following order:

- If `is.CorrelatedXZ = TRUE`, `X` and `Z` are generated *jointly* from
  the distribution specified in `distr_settings$XZ_distr`. User-supplied
  `X` or `Z` are **ignored** in this case.

- If `is.CorrelatedXZ = FALSE`, user-supplied matrices take
  precedence: - If `sim_settings$X` is supplied: `X` is used exactly as
  provided and its column dimension updates `p`. - If `sim_settings$Z`
  is supplied: `Z` is used exactly as provided and its column dimension
  updates `q`.

- If `X` is not supplied: `X` is generated using
  `distr_settings$X_distr`.

- If `Z` is not supplied: - If `is.ZInX = FALSE`, `Z` is generated from
  `distr_settings$Z_distr`; - otherwise `Z` is inherited from the *first
  `q` columns of `X`* (with correct intercept handling).

- Intercepts are added **after** stochastic generation, if requested.

- Orthogonalization of non-intercept columns is performed last using the
  [`qr`](https://rdrr.io/r/base/qr.html) function.

**Random effects** Random effects are drawn per subject from
`distr_settings$RE_distr`, with dimensions matched to the columns of
`Z`.

**Error distribution specification** The error distribution
(`error_distr`) may be:

- **Univariate** with scalar `sigma`, which corresponds to IID errors.

- **Multivariate**, with a covariance matrix `sigma` of dimension equal
  to the total number of observations, enabling arbitrary correlation
  and heteroscedasticity.

IID errors can also be represented through a multivariate specification
with diagonal covariance matrix in `sigma`. If `sigma` is not specified
in `error_distr$distr_params`, a default IID specification is
constructed.

**SNR control** When `SNR` is supplied and the error distribution is
IID, the error variance is automatically scaled to achieve the target
signal-to-noise ratio (defined as the ratio of variances of signals and
errors). SNR control is not supported for non-IID errors.

**Dependencies** Certain distribution choices (e.g., copulas) may
require additional packages. See
[`simulate_IIDdata`](https://abhik-stat.github.io/fsimR/reference/simulate_IIDdata.md)
for details on required packages and supported distributions.

## See also

[`simulate_IIDdata`](https://abhik-stat.github.io/fsimR/reference/simulate_IIDdata.md),
[`distr_spec`](https://abhik-stat.github.io/fsimR/reference/distr_spec.md),
[`sim_spec`](https://abhik-stat.github.io/fsimR/reference/sim_spec.md).
Also
[`summary.LMMdata`](https://abhik-stat.github.io/fsimR/reference/summarySimData.md),
[`print.LMMdata`](https://abhik-stat.github.io/fsimR/reference/printSimData.md),
[`plot.LMMdata`](https://abhik-stat.github.io/fsimR/reference/plotSimData.md)
for S3 methods applicable to the returned object

## Examples

``` r
# Basic simulation
set.seed(123)
sim1 <- simulate_LMMdata()
str(sim1, max.level = 1)
#> List of 10
#>  $ y             : num [1:50, 1] -0.3072 0.0231 1.812 0.3238 0.3826 ...
#>  $ X             : 'IIDdata' num [1:50, 1] 1 1 1 1 1 1 1 1 1 1 ...
#>   ..- attr(*, "dimnames")=List of 2
#>  $ Z             : 'IIDdata' num [1:50, 1] 1 1 1 1 1 1 1 1 1 1 ...
#>   ..- attr(*, "dimnames")=List of 2
#>  $ RE            : 'IIDdata' num [1:10, 1] 0.2533 -0.0285 -0.0429 1.3686 -0.2258 ...
#>  $ sigma_e       : num 1
#>  $ SNR           : num 0.674
#>  $ ID            : Factor w/ 10 levels "1","2","3","4",..: 1 1 1 1 1 2 2 2 2 2 ...
#>  $ n             : num 50
#>  $ sim_settings  :List of 15
#>  $ distr_settings:List of 5
#>  - attr(*, "class")= chr [1:2] "LMMdata" "list"


# Multiple replications and iterations
sim <- simulate_LMMdata(n_rep = c(2, 3)) # 2 reps × 3 iterations
length(sim$y)         # 2 replications
#> [1] 2
dim(sim$y[[1]])       # 50 × 3 iterations (default sample size 50)
#> [1] 50  3


# Non-Gaussian error and multivariate random effects
distr_settings <- list(
  error_distr = list(distr_name = "t", distr_params = distr_spec(df = 4)),
  RE_distr    = list(distr_name = "mvnorm", distr_params = list(dim = 2, sigma = diag(2)))
)
sim <- simulate_LMMdata(
  n_rep = 1,
  sim_settings = sim_spec(n_subj = 15, n_obs = 5),
  distr_settings = distr_settings
)


# Correlated X and Z design matrices
Sigma_X <- matrix(c(1, 0.5, 0.5, 1), ncol = 2)
Sigma_Z <- diag(2)
Sigma_XZ <- matrix(c(0.3, 0.2, 0, 0.4), ncol = 2)
Sigma <- rbind(cbind(Sigma_X, Sigma_XZ), cbind(t(Sigma_XZ), Sigma_Z))

sim_settings <- sim_spec(
  n_subj = 10,
  n_obs  = 4,
  is.CorrelatedXZ = TRUE,
  include.Xintercept = TRUE,
  include.Zintercept = TRUE
)
distr_settings <- list(
  X_distr  = list(distr_name = "mvnorm", distr_params = list(sigma = Sigma_X)),
  Z_distr  = list(distr_name = "mvnorm", distr_params = list(sigma = Sigma_Z)),
  XZ_distr = list(distr_name = "mvnorm", distr_params = list(sigma = Sigma))
)
sim <- simulate_LMMdata(sim_settings = sim_settings, distr_settings = distr_settings)
pairs(cbind(sim$X, sim$Z), main = "Correlated X predictors")



# User-supplied design matrices
n_subj <- 5; n_obs <- 6; n <- n_subj * n_obs
X_user <- matrix(rnorm(n*3), ncol=3)
Z_user <- matrix(rnorm(n*2), ncol=2)
sim_settings <- sim_spec(n_subj = n_subj, n_obs = n_obs, X = X_user, Z = Z_user, beta_coeff = c(1,0.5,-1))
sim <- simulate_LMMdata(sim_settings = sim_settings)


# Copula-based generation
library(copula)
normal_cop <- normalCopula(param = 0.6, dim = 2)
distr_settings <- list(
  X_distr = list(
    distr_name = "copula",
    distr_params = distr_spec(
      copula = normal_cop,
      margins = list(
        list(dist="norm", params=list(mean=0, sd=1)),
        list(dist="norm", params=list(mean=0, sd=2))
      )
    )
  )
)
sim <- simulate_LMMdata(sim_settings = sim_spec(n_subj = 10, n_obs = 3),
                         distr_settings = distr_settings)
plot(sim$X[,2:3], main = "X from Copula-based distribution")


```
