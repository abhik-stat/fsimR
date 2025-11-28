# Distribution Specification Constructor

Construct a unified distribution specification object (list) for use in
all simulator functions, e.g.,
[`simulate_IIDdata`](https://abhik-stat.github.io/fsimR/reference/simulate_IIDdata.md),
[`simulate_LMdata`](https://abhik-stat.github.io/fsimR/reference/simulate_LMdata.md),
etc.

The function resolves parameter *alias groups* (location, scale, shape),
infers the dimension of the distribution, validates compatibility across
parameter groups, and returns a fully specified list of distribution
parameters.

## Usage

``` r
distr_spec(
  dim = 1,
  mean = NULL,
  mu = NULL,
  location = NULL,
  xi = NULL,
  sd = NULL,
  scale = NULL,
  sigma = NULL,
  omega = NULL,
  Omega = NULL,
  V = NULL,
  alpha = NULL,
  lambda = 1,
  rate = 1,
  df = 5,
  nu = 5,
  df1 = 5,
  df2 = 5,
  ncp = 0,
  shape = 1,
  shape1 = 1,
  shape2 = 1,
  size = 1,
  prob = 0.5,
  meanlog = 0,
  sdlog = 1,
  min = 0,
  max = 1,
  copula = NULL,
  margins = NULL,
  generator = NULL,
  ...
)
```

## Arguments

- dim:

  Integer (\>0). Initial dimension of the random variable. Default: `1`.

- mean, mu, location, xi:

  Numeric. Aliases for the *location* parameter. May be scalars or
  vectors. Default: `rep(0, dim)`.

- sd, scale, sigma, omega, Omega, V:

  Numeric or matrix. Aliases for the *scale* or *covariance* parameter.
  Scalars are expanded to `scalar \(\times\) I(dim)`; matrices must be
  square.Default: `I(dim)`.

- alpha:

  Numeric. Slant / skewness parameter. May be scalar or vector. Default:
  `rep(0, dim)`.

- lambda:

  Numeric (\>=0). Mean parameter for Poisson-type distributions.
  Default: `1`.

- rate:

  Numeric (\>0). Rate parameter for distributions with positive supports
  (e.g., exponential, gamma). Default: `1`.

- df, nu:

  Numeric (\>=0). Degree of freedom for distributions with a single df
  parameter (e.g., \\t\\, skew-\\t\\). Default: `5`.

- df1, df2:

  Numeric (\>=0). Degrees of freedom for distributions with two df
  parameters (e.g., \\F\\). Default: `5` for both.

- ncp:

  Numeric (\>=0). Non-centrality parameter. Default: `0`.

- shape:

  Numeric (\>=0). Shape parameter for distributions with one shape
  parameter (e.g., gamma, Weibull). Default: `1`.

- shape1, shape2:

  Numeric (\>=0). Shape parameters for distributions with two shape
  parameters (e.g., beta). Default: `1` for both.

- size:

  Integer (\>0). Size parameter for binomial-type distributions.
  Default: `1`.

- prob:

  Numeric in `[0,1]`. Probability parameter. Default: `0.5`.

- meanlog, sdlog:

  Numeric. Mean and standard deviation on the log-scale. Defaults: `0`,
  `1`.

- min, max:

  Numeric. Lower and upper bounds of the support for bounded
  distributions, (e.g., uniform). Defaults: `0`, `1`.

- copula:

  Optional copula object (from the copula package), necessary for
  copula-based multivariate simulations.

- margins:

  Optional list of distributional specifications for each margins in
  copula models (see
  [`simulate_IIDdata`](https://abhik-stat.github.io/fsimR/reference/simulate_IIDdata.md)).

- generator:

  Optional user-defined function for custom data generation.

- ...:

  Additional named arguments. These are appended verbatim to the output
  list. If names conflict with existing parameters, the values supplied
  in `...` will overwrite the defaults.

## Value

A named list containing:

- Resolved dimension (`dim`)

- Resolved aliases (`mean`, `mu`, `location`, `xi`, `sd`, etc.)

- All remaining formal parameters

- Additional user-supplied entries from `...`

## Details

**Alias groups:** The following alias groups are supported (to date)
through the internal function `.resolve_alias_group()`:

- **Location:** `mean`, `mu`, `location`, `xi`

- **Scale / covariance:** `sd`, `scale`, `sigma`, `omega`, `Omega`

- **Shape / skewness:** `alpha`

Within each group, the first non-`NULL` argument in the order of formal
arguments is used. Scalars are automatically expanded. If no alias in a
group is supplied, the default is set to `rep(0, dim)` for vectors, and
to `diag(dim)` for matrices.

**Dimension inference and interaction:** The working dimension is
initialized from `dim` and may subsequently be updated by:

1.  location aliases (vector length),

2.  scale aliases (matrix order),

3.  alpha aliases (vector length).

If a later group implies a different dimension than earlier groups,
strict compatibility checks are enforced using `.check_dim_strict()`.

When the dimension changes due to a later alias group, dependent
parameters are reset to their default values at the new dimension.

**Error conditions:** Errors are raised when:

- scale aliases are non-scalar non-square matrices,

- vector aliases have incompatible lengths,

- multiple alias groups imply incompatible dimensions,

- user-supplied dimensions conflict with resolved dimension.

## See also

[`simulate_IIDdata`](https://abhik-stat.github.io/fsimR/reference/simulate_IIDdata.md)

## Examples

``` r
# Default univariate distribution specifications
distr_spec()
#> $dim
#> [1] 1
#> 
#> $mean
#> [1] 0
#> 
#> $mu
#> [1] 0
#> 
#> $location
#> [1] 0
#> 
#> $xi
#> [1] 0
#> 
#> $sd
#>      [,1]
#> [1,]    1
#> 
#> $scale
#>      [,1]
#> [1,]    1
#> 
#> $sigma
#>      [,1]
#> [1,]    1
#> 
#> $omega
#>      [,1]
#> [1,]    1
#> 
#> $Omega
#>      [,1]
#> [1,]    1
#> 
#> $V
#>      [,1]
#> [1,]    1
#> 
#> $alpha
#> [1] 0
#> 
#> $lambda
#> [1] 1
#> 
#> $rate
#> [1] 1
#> 
#> $df
#> [1] 5
#> 
#> $nu
#> [1] 5
#> 
#> $df1
#> [1] 5
#> 
#> $df2
#> [1] 5
#> 
#> $ncp
#> [1] 0
#> 
#> $shape
#> [1] 1
#> 
#> $shape1
#> [1] 1
#> 
#> $shape2
#> [1] 1
#> 
#> $size
#> [1] 1
#> 
#> $prob
#> [1] 0.5
#> 
#> $meanlog
#> [1] 0
#> 
#> $sdlog
#> [1] 1
#> 
#> $min
#> [1] 0
#> 
#> $max
#> [1] 1
#> 
#> $copula
#> NULL
#> 
#> $margins
#> NULL
#> 
#> $generator
#> NULL
#> 

# Multivariate distribution setup (dimension inferred from inputs)
distr_spec(mean = c(1, 2, 3), sigma = diag(3))
#> $dim
#> [1] 3
#> 
#> $mean
#> [1] 1 2 3
#> 
#> $mu
#> [1] 1 2 3
#> 
#> $location
#> [1] 1 2 3
#> 
#> $xi
#> [1] 1 2 3
#> 
#> $sd
#>      [,1] [,2] [,3]
#> [1,]    1    0    0
#> [2,]    0    1    0
#> [3,]    0    0    1
#> 
#> $scale
#>      [,1] [,2] [,3]
#> [1,]    1    0    0
#> [2,]    0    1    0
#> [3,]    0    0    1
#> 
#> $sigma
#>      [,1] [,2] [,3]
#> [1,]    1    0    0
#> [2,]    0    1    0
#> [3,]    0    0    1
#> 
#> $omega
#>      [,1] [,2] [,3]
#> [1,]    1    0    0
#> [2,]    0    1    0
#> [3,]    0    0    1
#> 
#> $Omega
#>      [,1] [,2] [,3]
#> [1,]    1    0    0
#> [2,]    0    1    0
#> [3,]    0    0    1
#> 
#> $V
#>      [,1] [,2] [,3]
#> [1,]    1    0    0
#> [2,]    0    1    0
#> [3,]    0    0    1
#> 
#> $alpha
#> [1] 0 0 0
#> 
#> $lambda
#> [1] 1
#> 
#> $rate
#> [1] 1
#> 
#> $df
#> [1] 5
#> 
#> $nu
#> [1] 5
#> 
#> $df1
#> [1] 5
#> 
#> $df2
#> [1] 5
#> 
#> $ncp
#> [1] 0
#> 
#> $shape
#> [1] 1
#> 
#> $shape1
#> [1] 1
#> 
#> $shape2
#> [1] 1
#> 
#> $size
#> [1] 1
#> 
#> $prob
#> [1] 0.5
#> 
#> $meanlog
#> [1] 0
#> 
#> $sdlog
#> [1] 1
#> 
#> $min
#> [1] 0
#> 
#> $max
#> [1] 1
#> 
#> $copula
#> NULL
#> 
#> $margins
#> NULL
#> 
#> $generator
#> NULL
#> 

# Specification for Gamma distribution
distr_spec(shape = 2, rate = 1)
#> $dim
#> [1] 1
#> 
#> $mean
#> [1] 0
#> 
#> $mu
#> [1] 0
#> 
#> $location
#> [1] 0
#> 
#> $xi
#> [1] 0
#> 
#> $sd
#>      [,1]
#> [1,]    1
#> 
#> $scale
#>      [,1]
#> [1,]    1
#> 
#> $sigma
#>      [,1]
#> [1,]    1
#> 
#> $omega
#>      [,1]
#> [1,]    1
#> 
#> $Omega
#>      [,1]
#> [1,]    1
#> 
#> $V
#>      [,1]
#> [1,]    1
#> 
#> $alpha
#> [1] 0
#> 
#> $lambda
#> [1] 1
#> 
#> $rate
#> [1] 1
#> 
#> $df
#> [1] 5
#> 
#> $nu
#> [1] 5
#> 
#> $df1
#> [1] 5
#> 
#> $df2
#> [1] 5
#> 
#> $ncp
#> [1] 0
#> 
#> $shape
#> [1] 2
#> 
#> $shape1
#> [1] 1
#> 
#> $shape2
#> [1] 1
#> 
#> $size
#> [1] 1
#> 
#> $prob
#> [1] 0.5
#> 
#> $meanlog
#> [1] 0
#> 
#> $sdlog
#> [1] 1
#> 
#> $min
#> [1] 0
#> 
#> $max
#> [1] 1
#> 
#> $copula
#> NULL
#> 
#> $margins
#> NULL
#> 
#> $generator
#> NULL
#> 

# Specification for Beta distribution
distr_spec(shape1 = 2, shape2 = 5)
#> $dim
#> [1] 1
#> 
#> $mean
#> [1] 0
#> 
#> $mu
#> [1] 0
#> 
#> $location
#> [1] 0
#> 
#> $xi
#> [1] 0
#> 
#> $sd
#>      [,1]
#> [1,]    1
#> 
#> $scale
#>      [,1]
#> [1,]    1
#> 
#> $sigma
#>      [,1]
#> [1,]    1
#> 
#> $omega
#>      [,1]
#> [1,]    1
#> 
#> $Omega
#>      [,1]
#> [1,]    1
#> 
#> $V
#>      [,1]
#> [1,]    1
#> 
#> $alpha
#> [1] 0
#> 
#> $lambda
#> [1] 1
#> 
#> $rate
#> [1] 1
#> 
#> $df
#> [1] 5
#> 
#> $nu
#> [1] 5
#> 
#> $df1
#> [1] 5
#> 
#> $df2
#> [1] 5
#> 
#> $ncp
#> [1] 0
#> 
#> $shape
#> [1] 1
#> 
#> $shape1
#> [1] 2
#> 
#> $shape2
#> [1] 5
#> 
#> $size
#> [1] 1
#> 
#> $prob
#> [1] 0.5
#> 
#> $meanlog
#> [1] 0
#> 
#> $sdlog
#> [1] 1
#> 
#> $min
#> [1] 0
#> 
#> $max
#> [1] 1
#> 
#> $copula
#> NULL
#> 
#> $margins
#> NULL
#> 
#> $generator
#> NULL
#> 

# Custom additional arguments
distr_spec(mean = 10, extra1 = "hello", extra2 = 42)
#> $dim
#> [1] 1
#> 
#> $mean
#> [1] 10
#> 
#> $mu
#> [1] 10
#> 
#> $location
#> [1] 10
#> 
#> $xi
#> [1] 10
#> 
#> $sd
#>      [,1]
#> [1,]    1
#> 
#> $scale
#>      [,1]
#> [1,]    1
#> 
#> $sigma
#>      [,1]
#> [1,]    1
#> 
#> $omega
#>      [,1]
#> [1,]    1
#> 
#> $Omega
#>      [,1]
#> [1,]    1
#> 
#> $V
#>      [,1]
#> [1,]    1
#> 
#> $alpha
#> [1] 0
#> 
#> $lambda
#> [1] 1
#> 
#> $rate
#> [1] 1
#> 
#> $df
#> [1] 5
#> 
#> $nu
#> [1] 5
#> 
#> $df1
#> [1] 5
#> 
#> $df2
#> [1] 5
#> 
#> $ncp
#> [1] 0
#> 
#> $shape
#> [1] 1
#> 
#> $shape1
#> [1] 1
#> 
#> $shape2
#> [1] 1
#> 
#> $size
#> [1] 1
#> 
#> $prob
#> [1] 0.5
#> 
#> $meanlog
#> [1] 0
#> 
#> $sdlog
#> [1] 1
#> 
#> $min
#> [1] 0
#> 
#> $max
#> [1] 1
#> 
#> $copula
#> NULL
#> 
#> $margins
#> NULL
#> 
#> $generator
#> NULL
#> 
#> $extra1
#> [1] "hello"
#> 
#> $extra2
#> [1] 42
#> 
```
