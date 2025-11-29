# Plot Simulated Data

Generates visualizations for simulated datasets, supporting multiple
types:

- **Mixed Models** (of class `LMMdata` or `GLMMdata`)

- **Regression Models** (of class `LMdata` or `GLMdata`)

- **IID samples** (of class `IIDdata`)

The function also works for named lists containing:

- response vectors as `y`

- fixed-effects design matrix as `X`

- optional random-effects design matrix as `Z`

- optional random-effects realizations as `RE`

Numeric vectors or matrices are treated as IID samples where columns
correspond to variables, rows to observations.

For plot types and layout, see Details.

## Usage

``` r
# S3 method for class 'LMMdata'
plot(
  x,
  replication = 1,
  iteration = 1,
  type = c("all", "y", "X", "Z", "RE"),
  main = TRUE,
  plot_args = list(),
  ...
)

# S3 method for class 'GLMMdata'
plot(x, ...)

# S3 method for class 'LMdata'
plot(
  x,
  replication = 1,
  iteration = 1,
  type = c("all", "y", "X"),
  main = TRUE,
  plot_args = list(),
  ...
)

# S3 method for class 'GLMdata'
plot(x, ...)

# S3 method for class 'IIDdata'
plot(x, main = TRUE, plot_args = list(), ...)

# S3 method for class 'list'
plot(x, ...)
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

- type:

  Character vector specifying which components to plot:

  - `"y"`: response vector

  - `"X"`: fixed-effects or design matrix

  - `"Z"`: random-effects design matrix (LMM/GLMM only)

  - `"RE"`: random effects (LMM/GLMM only)

  - `"all"`: all available components (default)

  Ignored for IID samples, where all columns are always plotted.

- main:

  Logical; if `TRUE`, descriptive titles are shown (default `TRUE`).

- plot_args:

  Optional list of graphical parameters to customize appearance of the
  plots. Supported arguments:

  breaks

  :   Number of histogram bins (default `10`)

  hist_col

  :   Histogram fill color (default `"lightsteelblue"`)

  dens_col

  :   Line color for density curve (default `"firebrick"`)

  dens_lwd

  :   Line width for density curve (default `2`)

  smooth_col

  :   Smoothing line color (default `"blue"`)

  smooth_lwd

  :   Smoothing line width (default `2`)

  cex_axis

  :   Character expansion for axis labels (default `1.4`)

  label_line

  :   Line offset for margin labels (default `0`)

  low_col

  :   Background color for `r = 0` in correlation panels (default
      `"#FFFFFF"`)

  high_col

  :   Background color for `|r| = 1` in correlation panels (default
      `"#2166AC"`)

  cex_text

  :   Character expansion for correlation/p-value text (default `1.4`)

  digits

  :   Number of digits for correlation/p-value (default `4`)

- ...:

  Additional graphical arguments passed to the standard plotting
  functions, such as `hist`, `lines`, `pairs` and `text` (via internal
  functions `.plot_hist_density` and `.plot_pairs()`).

## Value

Invisibly returns `NULL`. Generates plots for visual inspection.

## Details

The plotting behavior depends on the class and selected components (in
`type`):

1.  **Histograms with density overlay** for response vectors (`y`) or
    any numeric vector (treated as a univariate IID samples). Skipped
    for constant vectors.

2.  **Enhanced pairs plots** for matrices (`X`, `Z`, `RE`) or
    multivariate IID samples:

    - Diagonal panels: histograms with density scaled to `[0,1]`.

    - Upper panels: scatter plots with smoothing lines
      ([loess](https://rdrr.io/r/stats/loess.html)).

    - Lower panels: correlation coefficients, p-values, and significance
      stars (`*`, `**`, `***`). Background color reflects the absolute
      correlation; negative correlations appear in red.

    - Constant columns are automatically removed.

## Examples

``` r
if (FALSE) { # \dontrun{
# Mixed Model
sim_settings <- sim_spec(p = 3)
sim_data <- simulate_LMMdata(c(2,3), sim_settings)
plot(sim_data, type = c("y", "RE"))

# Regression Model
sim_data <- simulate_LMdata(c(2,3), sim_settings)
plot(sim_data, type = "X")

# IID data
sim_data <- simulate_IIDdata(n = 100)
plot(sim_data)

# Named list wrapper
mylist <- list(y = rnorm(50), X = matrix(rnorm(100), ncol=2))
plot(mylist)

} # }
```
