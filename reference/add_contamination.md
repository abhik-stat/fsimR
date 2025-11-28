# Add Contamination to Simulated Data

Introduces controlled contamination into simulated datasets.
Contamination can be applied to:

- **Mixed Model datasets** (of class `LMMdata` or `GLMMdata`)

- **Regression datasets** (of class `LMdata` or `GLMdata`)

- **IID samples** (of class `IIDdata`)

- **Named lists** containing components such as:

  - `y`: response vector or matrix

  - `X`: fixed-effects design matrix

  - `Z`: random-effects design matrix (optional) Numeric vectors or
    matrices are treated as IID samples (columns = variables).

## Usage

``` r
add_contamination(simData, cont_pos = list(), cont_settings = list())
```

## Arguments

- simData:

  Numeric vector, matrix, named list, or simulation object (of class
  `LMMdata`, `GLMMdata`, `LMdata`, `GLMdata`, `IIDdata`).

  - If a **list**, it must include the component specified in
    `cont_pos$name`.

  - If a **vector/matrix**, it is treated as the target component to
    contaminate.

- cont_pos:

  Named list specifying the target component to contaminate:

  name

  :   Character string, one of `"y"`, `"X"`, or `"Z"` (default: `y`).
      Determines which component of the simulation object is modified.

  col

  :   Vector of column indices to contaminate. If `col = NULL` (default)
      or `name = "y"`, all columns are contaminated. For `y`, multiple
      columns usually represent replications in the simulation object,
      and all such columns are contaminated consistently.

- cont_settings:

  Named list specifying contamination details:

  cont_mode

  :   Mode of contamination: `"casewise"` (default) or `"cellwise"`. -
      `"casewise"`: Same rows contaminated across all specified
      columns. - `"cellwise"`: Different rows may be contaminated
      independently per column.

  cont_prop

  :   Proportion of contamination `[0,1]`. In `casewise` mode, this is
      the proportion of rows. In `cellwise` mode, this is the proportion
      of cells across all selected columns. (default: `c(0,1)`)

  cont_value

  :   Optional numeric scalar or vector specifying contamination values.
      If provided, this overrides `cont_distr` and no random generation
      occurs. *Length must be either 1 or equal to the number of
      contaminated columns*.

  cont_distr

  :   Named list describing the distribution of contamination values
      (only when `cont_value` is not provided). Accepts:

      - `distr_name` Distribution name (e.g., `"norm"`, `"t"`,
        `"unif"`).

      - `distr_params` List of distribution parameters (see
        [simulate_IIDdata](https://abhik-stat.github.io/fsimR/reference/simulate_IIDdata.md)).

      - `generator` Optional custom generator function.

      The default is the standard normal distribution.

  cont_type

  :   Character string specifying the contamination mechanism:

      - "additive": Adds contamination values to selected rows/cells
        (default).

      - "multiplicative": Multiplies values in selected rows/cells by
        contamination values.

      - "replace": Replaces values in selected rows/cells with
        contamination values.

      - "outlier": Sets selected rows/cells to extreme values of the
        form: \$\$\text{mean} \pm (\text{sd} \times
        \text{outlier_factor}).\$\$ This generates classical point
        outliers.

      - "leverage": Generates high-leverage points in `X` or `Z`.
        Corresponding rows in `y` (if present) are shifted in the
        opposite direction. See *Details* for full description.

  leverage_factor

  :   Numeric multiplier controlling the magnitude of leverage
      contamination (default: `5`).

  outlier_factor

  :   Numeric multiplier controlling extreme values for `"outlier"` or
      `y` shifts in `"leverage"` (default: `5`).

  outlier_dir

  :   Probability that an outlier is generated in the right (positive)
      tail versus the left (negative) tail. Default: `1` (all outliers
      on the right).

## Value

Modified simulated datasets of the same class/structure as the input,
with the requested contamination applied.

- If `simData` is a list, the contaminated component is replaced.

- If a vector or matrix is supplied, the returned object is the
  contaminated matrix.

## Details

The function provides a general framework for injecting various forms of
contamination into simulated data, useful for robustness studies and
sensitivity analyses.

**Contamination mode:**

- *Casewise*: A total of `ceil(n × cont_prop)` rows are selected and
  contaminated across all columns.

- *Cellwise*: A proportion of individual cells is contaminated
  independently in each specified column. Rows are re-sampled
  independently per column.

*Leverage contamination:* For `cont_type = "leverage"` applied to `X` or
`Z`, a subset of rows/cells is randomly selected for contamination.
Their values are set to the column mean plus a random shift, drawn from
a normal distribution with the column's standard deviation scaled by
`leverage_factor`. If a `y` component exists, the corresponding rows are
shifted in the opposite direction, scaled by `outlier_factor`, creating
influential points within the regression framework.

*Distribution-based contamination:* When `cont_value = NULL`, random
contamination values are generated using
[simulate_IIDdata](https://abhik-stat.github.io/fsimR/reference/simulate_IIDdata.md),
allowing flexible custom distributions or user-provided generators.

## Examples

``` r
# Additive contamination to an IID sample
set.seed(1)
x <- rnorm(100)
contaminated <- add_contamination(
  x,
  cont_pos = list(name="y", col=1),
  cont_settings = list(
    cont_prop = 0.1,
    cont_type = "additive",
    cont_distr = list(distr_name="norm", distr_params=list(mean=5, sd=1))
  )
)


# Outlier contamination in regression response
set.seed(2)
sim <- list(y = rnorm(200), X = cbind(1, rnorm(200)))

sim_outliers <- add_contamination(
  sim,
  cont_pos = list(name="y"),
  cont_settings = list(
    cont_type="outlier",
    outlier_factor=8,
    outlier_dir=0.8
  )
)


# High-leverage contamination in design matrix X
set.seed(3)
sim <- list(
  y = rnorm(150),
  X = cbind(1, rnorm(150), rnorm(150))
)

sim_lev <- add_contamination(
  sim,
  cont_pos = list(name="X", col=2:3),
  cont_settings = list(
    cont_type="leverage",
    leverage_factor=10,
    outlier_factor=6
  )
)

# Cellwise contamination (different rows per column)
set.seed(4)
mat <- matrix(rnorm(50*3), ncol=3)
cont_cell <- add_contamination(
  mat,
  cont_pos = list(name="X", col=1:3),
  cont_settings = list(
    cont_prop=0.1,
    cont_type="replace",
    cont_value=99,
    cont_mode="cellwise"
  )
)
```
