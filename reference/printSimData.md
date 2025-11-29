# Print Simulated Dataset

Prints a dataset from a simulation object, showing metadata and actual
data. Supports multiple types of simulated datasets:

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

## Usage

``` r
# S3 method for class 'LMMdata'
print(x, replication = 1, iteration = 1, ...)

# S3 method for class 'GLMMdata'
print(x, replication = 1, iteration = 1, ...)

# S3 method for class 'LMdata'
print(x, replication = 1, iteration = 1, ...)

# S3 method for class 'GLMMdata'
print(x, replication = 1, iteration = 1, ...)
```

## Arguments

- x:

  Numeric vector, matrix, named list, or simulation object (of class
  `LMMdata`, `GLMMdata`, `LMdata`, `GLMdata`, `IIDdata`). If a list, it
  must include the components required for the selected `type`.

- replication:

  Integer; which replication (dataset) to plot (LMM/GLMM, LM/GLM only;
  default `1`).

- iteration:

  Integer; which iteration within the replication to plot (LMM/GLMM,
  LM/GLM only; default `1`).

- ...:

  Any additional arguments.

## Details

The function prints:

- Metadata including number of subjects/observations, dimensions of
  design matrices, regression or fixed-effects coefficients, `SNR`, and
  `sigma_e`.

- Optional distribution names if `distr_settings` is provided within
  object.

- The actual dataset (`y`, `X`, `Z`, `RE`) for the specified replication
  and iteration.

## Examples

``` r
if (FALSE) { # \dontrun{
# Mixed model
sim_data <- simulate_LMMdata(c(2,3), sim_spec(p=5))
print(sim_data, replication = 1, iteration = 2)

# Regression model
sim_data_lm <- simulate_LMdata(c(2,3), sim_spec(p=5))
print(sim_data_lm)

# IID samples
iid_data <- matrix(rnorm(50), ncol=5)
print(iid_data)

# Named list
mylist <- list(y = rnorm(10), X = matrix(rnorm(20), ncol=2))
print(mylist)
} # }
```
