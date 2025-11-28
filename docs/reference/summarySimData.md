# Summarize Simulated Data

Computes and prints a comprehensive summary for simulated datasets,
supporting multiple types:

- **Mixed Models** (of class `LMMdata` or `GLMMdata`)

- **Regression Models** (of class `LMdata` or `GLMdata`)

- **IID samples** (of class `IIDdata`)

The function also works for named lists containing:

- response vectors as `y`

- fixed-effects design matrix as `X`

- optional random-effects design matrix as `Z`

- optional random-effects realizations as `RE`

Numeric vectors or matrices are treated as IID samples (columns =
variables).

For details of summary statistics, see Details.

## Usage

``` r
# S3 method for class 'LMMdata'
summary(object, replication = 1, iteration = 1)

# S3 method for class 'GLMMdata'
summary(object, replication = 1, iteration = 1)

# S3 method for class 'LMdata'
summary(object, replication = 1, iteration = 1)

# S3 method for class 'GLMdata'
summary(object, replication = 1, iteration = 1)

# S3 method for class 'IIDdata'
summary(data)

# S3 method for class 'numeric'
summary(object)

# S3 method for class 'matrix'
summary(object)

# S3 method for class 'list'
summary(object, replication = 1, iteration = 1)
```

## Arguments

- object:

  Numeric vector, matrix, named list, or simulation object (of class
  `LMMdata`, `GLMMdata`, `LMdata`, `GLMdata`, `IIDdata`). If a list, it
  must include the components required for the selected `type`.

- replication:

  Integer; which replication (dataset) to plot (LMM/GLMM, LM/GLM only;
  default `1`).

- iteration:

  Integer; which iteration within the replication to plot (LMM/GLMM,
  LM/GLM only; default `1`).

## Value

Invisibly returns a numeric matrix of computed summary statistics:

- **Rows:** variables present in the object (e.g., `y`, `X`, `Z`, `RE`).

- **Columns:** `Mean`, `SD`, `Min`, `Q1`, `Median`, `Q3`, `Max`, and
  total count `n`.

## Details

The function prints a comprehensive summary for the simulated dataset:

- Simulation metadata such as number of subjects, observation counts,
  dimensions of `X`, `Z`, `RE` (for LMM/GLMM), fixed-effects
  coefficients, signal-to-noise ratio (`SNR`), and residual standard
  deviation (`sigma_e`).

- Optional distribution names if `distr_settings` is provided within
  object.

- Standard summary statistics (Mean, SD, Min, Q1, Median, Q3, Max) and
  total number of observations (`n`) for all variables in the object
  (`y`, `X`, `Z`, `RE` for mixed models; `y`, `X` for regression models;
  all columns for IID samples) for the specified `replication` and
  `iteration`.

- Correlation matrices for `X`, `Z`, and `RE` (if present), excluding
  intercept or any constant columns.

## Examples

``` r
if (FALSE) { # \dontrun{
# Mixed Model
sim_settings <- sim_spec(p = 5)
lmm_data <- simulate_LMMdata(c(2,3), sim_settings)
summary(lmm_data, replication = 1, iteration = 1)

# Regression Model
lm_data <- simulate_LMdata(c(2,3), sim_settings)
summary(lm_data, replication = 1, iteration = 1)

# IID Samples (numeric matrix)
iid_data <- matrix(rnorm(100), ncol = 5)
summary(iid_data)

# Named list with components like LMM
mylist <- list(
  y = lmm_data$y,
  X = lmm_data$X,
  Z = lmm_data$Z,
  RE = lmm_data$RE
)
summary(mylist, replication = 1, iteration = 1)
} # }
```
