# Simulation Specification Constructor

Constructs a unified simulation specification object (list) for
generating data from several structured statistical models, e.g., linear
models (LM), generalized linear models (GLM), linear mixed models (LMM),
generalized LMMs (GLMM), or other regression-type models. The output may
directly be passed to the corresponding simulator functions, such as
[`simulate_LMdata`](https://abhik-stat.github.io/fsimR/reference/simulate_LMdata.md),
[`simulate_LMMdata`](https://abhik-stat.github.io/fsimR/reference/simulate_LMMdata.md),
etc.

## Usage

``` r
sim_spec(
  name = "Setup0",
  n_subj = 10,
  n_obs = 5,
  p = NULL,
  q = 1,
  beta_coeff = NULL,
  SNR = NULL,
  X = NULL,
  Z = NULL,
  is.CorrelatedXZ = FALSE,
  is.ZInX = FALSE,
  orthogonalize.X = FALSE,
  orthogonalize.Z = FALSE,
  include.Xintercept = TRUE,
  include.Zintercept = TRUE,
  ...
)
```

## Arguments

- name:

  Character. Label for the simulation setup. Default: `"Setup0"`.

- n_subj:

  Integer (\>0). Number of subjects or clusters. Default: `10`.

- n_obs:

  Integer scalar or integer vector (\>0). Specifies the number of
  observations per subject/cluster.

  - If a **scalar** is supplied, all `n_subj` clusters are assigned the
    same number of observations.

  - If a **vector**, it must be of length `n_subj`, with each entry
    giving the size of the respective cluster. Useful for longitudinal
    or clustered designs. Default: `5` for all clusters.

- p:

  Integer (\>0). Number of fixed-effect predictors (including
  intercept). If `NULL`, it is inferred as `length(beta_coeff)` whenever
  possible. Default fallback value: `1`

- q:

  Integer (\>0). Number of random-effect predictors. Default: `1`.

- beta_coeff:

  Numeric vector. Fixed-effect coefficients of length `p`. If supplied,
  its length is used to infer `p` when `p = NULL`.

- SNR:

  Numeric (\>0). Signal-to-noise ratio.

- X:

  Numeric matrix or list. User-supplied design matrix for fixed effects.
  If provided, it is stored verbatim.

- Z:

  Numeric matrix or list. User-supplied design matrix for random
  effects. If provided, it is stored verbatim.

- is.CorrelatedXZ:

  Logical. Whether `X` and `Z` are correlated (and hence they need to be
  generated jointly). Default: `FALSE`.

- is.ZInX:

  Logical. Whether the columns of `Z` are the same as the initial
  columns of `X`. Default: `FALSE`.

- orthogonalize.X:

  Logical. Whether orthogonalize (non-intercept) columns of `X` after
  construction. Default: `FALSE`.

- orthogonalize.Z:

  Logical. Whether orthogonalize (non-intercept) columns of `Z` after
  construction. Default: `FALSE`.

- include.Xintercept:

  Logical. Whether an intercept column should be included in `X`.
  Default: `TRUE`.

- include.Zintercept:

  Logical. Whether an intercept column should be included in `Z`.
  Default: `TRUE`.

- ...:

  Additional named arguments to append to the returned specification
  list.

  \#' @return A named list containing:

  - `name` — Setup label.

  - `n_subj`, `n_obs` — Sample size arguments.

  - `p`, `q` — Predictor dimensions (as supplied or inferred).

  - `beta_coeff` — Fixed-effect coefficients.

  - `SNR` - Desired signal-to-noise ratio.

  - `X`, `Z` — Optional user-supplied design matrices.

  - `is.CorrelatedXZ`, `is.ZInX` — Logical flags, specifying desired
    relationships between the design matrices `X` and `Z` .

  - `orthogonalize.X`, `orthogonalize.Z` — Orthogonalization flags.

  - `include.Xintercept`, `include.Zintercept` — Intercept controls.

  - Any additional entries from `...`.

## See also

[`simulate_LMdata`](https://abhik-stat.github.io/fsimR/reference/simulate_LMdata.md),
[`simulate_LMMdata`](https://abhik-stat.github.io/fsimR/reference/simulate_LMMdata.md)

## Examples

``` r
# Basic specification
spec1 <- sim_spec()

# Custom fixed effects
spec2 <- sim_spec(beta_coeff = c(1, -1, 0.5))

# Mixed-effects style specification
spec3 <- sim_spec(
  n_subj = 15, n_obs = 6,
  beta_coeff = c(1, 0.5),
)
```
