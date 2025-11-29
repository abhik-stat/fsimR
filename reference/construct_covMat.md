# Construct Structured Covariance or Correlation Matrices (non-random)

Constructs a \\dim \times dim\\ covariance or correlation matrix
according to a specified structural type. Supported types include
diagonal, compound symmetry, block, Toeplitz,
autoregressive/moving-average processes (AR/MA/ARMA/ARIMA), factor
models, graph-based precision structures, spatial kernels, and
spectrum-constrained matrices.

## Usage

``` r
construct_covMat(
  dim,
  corr = FALSE,
  type = "diag",
  params = list(),
  eigen_check = TRUE
)
```

## Arguments

- dim:

  Integer. Dimension of the resulting matrix.

- corr:

  Logical. If `TRUE`, returns a correlation matrix. If `FALSE`
  (default), returns a covariance matrix scaled by `variances`. For
  types `factor` and `graph`, `variances` are ignored.

- type:

  Character. Structure type of the covariance/correlation matrix.
  Supported values:

  `"diag"`

  :   Diagonal matrix (identity correlation).

  `"CS"`

  :   Compound symmetry (exchangeable correlation).

  `"block"`

  :   Block-diagonal matrix.

  `"toeplitz"`

  :   Toeplitz matrix with lag correlations (AR(1) structure).

  `"AR"`

  :   Autoregressive AR(p) process.

  `"MA"`

  :   Moving average MA(q) process.

  `"ARMA"`

  :   ARMA(p, q) process.

  `"ARIMA"`

  :   ARIMA(p, d, q) process.

  `"SAR"`

  :   Seasonal AR process.

  `"SMA"`

  :   Seasonal MA process.

  `"SARMA"`

  :   Seasonal ARMA process.

  `"SARIMA"`

  :   Seasonal ARIMA process.

  `"factor"`

  :   Factor model: loadings + uniquenesses.

  `"graph"`

  :   Graph-structured precision matrix (GMRF-style). Uses precision =
      adjacency + diagonal dominance; covariance = inverse of precision.

  `"spatial"`

  :   Spatial correlation via a distance matrix and kernel function.

  `"specConstr"`

  :   Spectrum-constrained matrix with specified eigenvalues, rank or
      condition number.

- params:

  Named list of parameters required for the selected type. Common
  parameters:

  `variances`

  :   Numeric scalar or vector; used if `corr = FALSE` and `type` is
      neither `factor` nor `graph`. Defaults to 1.

  `rho`

  :   Numeric scalar, vector, or matrix/list. Correlation parameter(s)

  `groups`

  :   Vector of group sizes for block-diagonal matrices.

  `arma_phi`, `arma_theta`

  :   AR and MA coefficients for time series types.

  `d`

  :   Differencing order for ARIMA/SARIMA processes.

  `D`, `seasonal_period`

  :   Seasonal differencing order and period for SARIMA/SARMA/SMA/SAR
      processes.

  `seasonal_phi`, `seasonal_theta`

  :   Seasonal AR/MA coefficients.

  `loadings`, `psi`

  :   Factor loadings and uniquenesses for factor models.

  `adjacency`

  :   Symmetric adjacency matrix for graph-based covariance.

  `distances`

  :   Distance matrix for spatial correlation. Defaults to absolute
      differences.

  `kernel`

  :   Kernel function for spatial correlation. Can be a function or one
      of: `"exponential"`, `"gaussian"`, `"matern"`,
      `"rational_quadratic"`, `"powered_exponential"`, `"cauchy"`.

  `phi`, `nu`, `alpha`, `beta`, `power`

  :   Spatial kernel parameters.

  `eigenvalues`, `rank`, `cond_num`

  :   Parameters for spectrum-constrained matrices. Specified
      eigenvalues, rank and condition numbers respectively.

  bandwidth

  :   Integer; optional for CS or Toeplitz matrices to zero out entries
      beyond a certain lag.

  max_iter, tol

  :   Maximum iterations and tolerance for spectrum-constrained
      iterative projections. See Details.

- eigen_check:

  Logical. If `TRUE` (default), checks if the resulting matrix is
  PSD/PD.

## Value

A symmetric `dim` \\\times\\ `dim` matrix of class `covMat`. Either
correlation or covariance depending on `corr`.

## Details

The `params` list must contain type-specific elements as follows:

- "diag":

  Optional: `variances` (numeric vector or scalar) if `corr = FALSE`.
  Defaults to 1.

- "CS":

  `rho` (scalar correlation), optional `variances` (if `corr = FALSE`),
  optional `bandwidth`.

- "block":

  `groups` (vector of integers specifying block sizes), `rho` (scalar or
  vector for each block), optional `variances`.

- "toeplitz":

  `rho` (vector or scalar), optional `variances`, optional `bandwidth`.

- "AR":

  `arma_phi` (AR coefficients), optional `variances`.

- "MA":

  `arma_theta` (MA coefficients), optional `variances`.

- "ARMA":

  `arma_phi` (AR coefficients), `arma_theta` (MA coefficients), optional
  `variances`.

- "ARIMA":

  `arma_phi` (AR), `arma_theta` (MA), `d` (differencing order), optional
  `variances`.

- "SAR":

  `seasonal_phi` (seasonal AR), `D` (seasonal differencing),
  `seasonal_period`.

- "SMA":

  `seasonal_theta` (seasonal MA), `D`, `seasonal_period`.

- "SARMA":

  `seasonal_phi`, `seasonal_theta`, `D`, `seasonal_period`.

- "SARIMA":

  `arma_phi`, `arma_theta`, `d`, `seasonal_phi`, `seasonal_theta`, `D`,
  `seasonal_period`.

- "factor":

  `loadings` (matrix `dim` \\\times\\ `k`), `psi` (vector of
  uniquenesses length dim).

- "graph":

  `adjacency` (symmetric matrix `dim` \\\times\\ `dim`).

- "spatial":

  `distances` (`dim` \\\times\\ `dim` matrix), `kernel` (function or
  string: `"exponential"`, `"gaussian"`, `"matern"`,
  `"rational_quadratic"`, `"powered_exponential"`, `"cauchy"`),
  kernel-specific parameters (`phi`, `nu`, `alpha`, `beta`, `power`).

- "specConstr":

  `eigenvalues` (numeric vector), `rank` (integer \<= `dim`), `cond_num`
  (numeric), optional `max_iter` and `tol`.

**Time series types:** Non-stationary models use MA(\\\infty\\)
expansion with differencing. Roots of AR/MA polynomials are stabilized
for invertibility.

**Spectrum-Constrained constructions:** Alternating projections enforce
target eigenvalues/rank/condition number while preserving diagonal. The
matrix supplied in `loadings` is used to generate an initial positive
semi-definite matrix for iterative projection onto the target spectrum.
Default base: a `dim` \\\times\\ `rank_target` Vandermonde-like matrix
of sin(0 : \\\pi\\/2) values, scaled by 10

**PSD/PD Guarantee:** All matrices constructed by this function are
guaranteed to be symmetric and either positive semi-definite (PSD) or
positive definite (PD), depending on the type and parameterization:
Supplied parameters are validated to ensure PSD/PD properties for all
structures, expect for `type = specConstr`. In spectrum-constrained
cases, particularly when `rank < dim` or user-specified `eigenvalues`
are provided, the matrix may be near PSD. Before returning the matrix,
the function performs a *final validation step* (if
`eigen_check = TRUE`) to verify that the output is truly PSD/PD and
issues errors, warnings, or messages accordingly.

## Examples

``` r
# Diagonal covariance
construct_covMat(5, type = "diag", params = list(variances = 1:5))
#> The constructed matrix is PD.
#>      [,1] [,2] [,3] [,4] [,5]
#> [1,]    1    0    0    0    0
#> [2,]    0    2    0    0    0
#> [3,]    0    0    3    0    0
#> [4,]    0    0    0    4    0
#> [5,]    0    0    0    0    5
#> attr(,"class")
#> [1] "covMat" "matrix" "array" 

# AR(1) correlation
construct_covMat(6, corr = TRUE, type = "AR", params = list(rho = 0.8))
#> The constructed matrix is PD.
#>         [,1]   [,2]  [,3]  [,4]   [,5]    [,6]
#> [1,] 1.00000 0.8000 0.640 0.512 0.4096 0.32768
#> [2,] 0.80000 1.0000 0.800 0.640 0.5120 0.40960
#> [3,] 0.64000 0.8000 1.000 0.800 0.6400 0.51200
#> [4,] 0.51200 0.6400 0.800 1.000 0.8000 0.64000
#> [5,] 0.40960 0.5120 0.640 0.800 1.0000 0.80000
#> [6,] 0.32768 0.4096 0.512 0.640 0.8000 1.00000
#> attr(,"class")
#> [1] "covMat" "matrix" "array" 

# Block matrix
construct_covMat(5, type = "block", params = list(groups = c(2,3), rho = c(0.5, 0.3)))
#> The constructed matrix is PD.
#>      [,1] [,2] [,3] [,4] [,5]
#> [1,]  1.0  0.5  0.0  0.0  0.0
#> [2,]  0.5  1.0  0.0  0.0  0.0
#> [3,]  0.0  0.0  1.0  0.3  0.3
#> [4,]  0.0  0.0  0.3  1.0  0.3
#> [5,]  0.0  0.0  0.3  0.3  1.0
#> attr(,"class")
#> [1] "covMat" "matrix" "array" 

# Factor model
L <- matrix(rnorm(10), 5, 2); psi <- rep(0.2, 5)
construct_covMat(5, type = "factor", params = list(loadings = L, psi = psi))
#> The constructed matrix is PD.
#>            [,1]        [,2]       [,3]       [,4]        [,5]
#> [1,]  0.4412270 -0.47116779 -0.7049116  0.4740768 -0.72719326
#> [2,] -0.4711678  1.66880924  1.5084694 -0.6968435  0.04974662
#> [3,] -0.7049116  1.50846940  2.2914734 -1.3303597  1.79609481
#> [4,]  0.4740768 -0.69684353 -1.3303597  1.2274031 -2.00167352
#> [5,] -0.7271933  0.04974662  1.7960948 -2.0016735  5.81701618
#> attr(,"class")
#> [1] "covMat" "matrix" "array" 

# Spatial correlation with Matern kernel
D <- as.matrix(dist(matrix(runif(25), ncol = 5)))
construct_covMat(5, type = "spatial",
                    params = list(distances = D, kernel = "matern", phi = 1, nu = 0.5))
#> The constructed matrix is PD.
#>           [,1]      [,2]      [,3]      [,4]      [,5]
#> [1,] 1.0000000 0.5297695 0.5202922 0.6258200 0.4647797
#> [2,] 0.5297695 1.0000000 0.3675388 0.6379802 0.3467122
#> [3,] 0.5202922 0.3675388 1.0000000 0.5453145 0.4331569
#> [4,] 0.6258200 0.6379802 0.5453145 1.0000000 0.4632974
#> [5,] 0.4647797 0.3467122 0.4331569 0.4632974 1.0000000
#> attr(,"class")
#> [1] "covMat" "matrix" "array" 

# Spectrum-constrained correlation matrix
construct_covMat(4, type = "specConstr", params = list(eigenvalues = c(2,1,0.5,0.5)))
#> The constructed matrix is PD.
#>           [,1]      [,2]      [,3]      [,4]
#> [1,] 1.0000000 0.4314820 0.1874848 0.1097512
#> [2,] 0.4314820 1.0000000 0.3959925 0.3411828
#> [3,] 0.1874848 0.3959925 1.0000000 0.4933676
#> [4,] 0.1097512 0.3411828 0.4933676 1.0000000
#> attr(,"class")
#> [1] "covMat" "matrix" "array" 
```
