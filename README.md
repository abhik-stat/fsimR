
<!-- README.md is generated from README.Rmd. Please edit that file -->

# fsimR

`fsimR` is an R package for simulating data from a wide range of
statistical models. It provides a unified flexible framework for
simulating datasets from various statistical models, along with
summarizing, printing and visualizing simulated data. The package
support generation of independent and identically distributed (IID) data
under a variety of distributions, including user-defined generators and
copula-based multivariate sampling. It also includes functions for
generating data under linear models, generalized linear models,
mixed-effects models, and many more, with customizable covariate
distributions, regression coefficients, and error structures. The
package is designed to support reproducible simulation studies,
teaching, benchmarking, and Monte Carlo experiments. Its consistent
interface and tidyverse-friendly design make it easy to extend and
integrate into larger statistical workflows.

------------------------------------------------------------------------

## Installation

You can install the development version from GitHub using `devtools`:

``` r
# Install devtools if you haven't already
install.packages("devtools")

# Install fsimR from GitHub
devtools::install_github("abhik-stat/fsimR")
```

------------------------------------------------------------------------

## Usage Example

Here is a simple example demonstrating how to define a simulation
specification, simulate data, and plot the results:

``` r
library(fsimR)

# Define a linear model simulation specification
spec <- simSpec(
  formula = y ~ x1 + x2,
  distribution = list(x1 = "normal", x2 = "normal", y = "normal"),
  params = list(x1 = c(0,1), x2 = c(0,1), y = c(0,1))
)

# Simulate data for a linear model
sim_data <- simulate_LMdata(spec)

# Print a summary of simulated data
print(sim_data)

# Plot simulated data
plot(sim_data)
```

------------------------------------------------------------------------

## Functions & Documentation

The `fsimR` package provides the following functions. Internal functions
(starting with `.`) are not intended for direct use.

### Helper Functions

| Function Name | Documentation |
|----|----|
| `match.fun.allR` | [match.fun.allR.Rd](https://github.com/abhik-stat/fsimR/blob/master/man/match.fun.allR.Rd) |
| `distr_spec` | [distr_spec.Rd](https://github.com/abhik-stat/fsimR/blob/master/man/distr_spec.Rd) |
| `sim_spec` | [sim_spec.Rd](https://github.com/abhik-stat/fsimR/blob/master/man/sim_spec.Rd) |

### Exported Functions & Aliases

| Function Name | Documentation | Description |
|----|----|----|
| `distrSpec` | [distr_spec.Rd](https://github.com/abhik-stat/fsimR/blob/master/man/distr_spec.Rd) | User-friendly alias for `distr_spec`. |
| `simSpec` | [sim_spec.Rd](https://github.com/abhik-stat/fsimR/blob/master/man/sim_spec.Rd) | Exported simulation specification object. |
| `add_contamination` | [add_contamination.Rd](https://github.com/abhik-stat/fsimR/blob/master/man/add_contamination.Rd) | Adds contamination/noise to simulated data. |
| `simulate_IIDdata` | [simulate_IIDdata.Rd](https://github.com/abhik-stat/fsimR/blob/master/man/simulate_IIDdata.Rd) | Simulates IID data. |
| `simulate_LMdata` | [simulate_LMdata.Rd](https://github.com/abhik-stat/fsimR/blob/master/man/simulate_LMdata.Rd) | Simulates linear model data. |
| `simulate_LMMdata` | [simulate_LMMdata.Rd](https://github.com/abhik-stat/fsimR/blob/master/man/simulate_LMMdata.Rd) | Simulates linear mixed model data. |

### Plot Methods for Simulated Data

| Function Name | Documentation | Description |
|----|----|----|
| `plotSimData` | [plotSimData.Rd](https://github.com/abhik-stat/fsimR/blob/master/man/plotSimData.Rd) | Main plotting function for simulated data objects. |
| `plot.GLMdata` | [plotSimData.Rd](https://github.com/abhik-stat/fsimR/blob/master/man/plotSimData.Rd) | Plot method for GLM simulation output. |
| `plot.GLMMdata` | [plotSimData.Rd](https://github.com/abhik-stat/fsimR/blob/master/man/plotSimData.Rd) | Plot method for GLMM simulation output. |
| `plot.IIDdata` | [plotSimData.Rd](https://github.com/abhik-stat/fsimR/blob/master/man/plotSimData.Rd) | Plot method for IID simulation output. |
| `plot.LMdata` | [plotSimData.Rd](https://github.com/abhik-stat/fsimR/blob/master/man/plotSimData.Rd) | Plot method for LM simulation output. |
| `plot.LMMdata` | [plotSimData.Rd](https://github.com/abhik-stat/fsimR/blob/master/man/plotSimData.Rd) | Plot method for LMM simulation output. |
| `plot.list` | [plotSimData.Rd](https://github.com/abhik-stat/fsimR/blob/master/man/plotSimData.Rd) | Plot method for lists of simulated data. |

### Print Methods for Simulated Data

| Function Name | Documentation | Description |
|----|----|----|
| `printSimData` | [printSimData.Rd](https://github.com/abhik-stat/fsimR/blob/master/man/printSimData.Rd) | Main print method for simulated data objects. |
| `print.GLMMdata` | [printSimData.Rd](https://github.com/abhik-stat/fsimR/blob/master/man/printSimData.Rd) | Print method for GLMM simulation output. |
| `print.LMdata` | [printSimData.Rd](https://github.com/abhik-stat/fsimR/blob/master/man/printSimData.Rd) | Print method for LM simulation output. |
| `print.LMMdata` | [printSimData.Rd](https://github.com/abhik-stat/fsimR/blob/master/man/printSimData.Rd) | Print method for LMM simulation output. |

------------------------------------------------------------------------

## Citation

If you use `fsimR` in your research, please cite it:

``` r
citation("fsimR")
```

------------------------------------------------------------------------

## Contributing

Contributions are welcome! Please open issues or submit pull requests
via the [GitHub repository](https://github.com/abhik-stat/fsimR).

------------------------------------------------------------------------

## License

# `fsimR` is licensed under the Apache License 2.0.

<!-- README.md is generated from README.Rmd. Please edit that file -->

# fsimR

`fsimR` is an R package for simulating data from a wide range of
statistical models. It provides a unified flexible framework for
simulating datasets from various statistical models, along with
summarizing, printing and visualizing simulated data. The package
support generation of independent and identically distributed (IID) data
under a variety of distributions, including user-defined generators and
copula-based multivariate sampling. It also includes functions for
generating data under linear models, generalized linear models,
mixed-effects models, and many more, with customizable covariate
distributions, regression coefficients, and error structures. The
package is designed to support reproducible simulation studies,
teaching, benchmarking, and Monte Carlo experiments. Its consistent
interface and tidyverse-friendly design make it easy to extend and
integrate into larger statistical workflows.

------------------------------------------------------------------------

## Installation

You can install the development version from GitHub using `devtools`:

``` r
# Install devtools if you haven't already
install.packages("devtools")

# Install fsimR from GitHub
devtools::install_github("abhik-stat/fsimR")
```

------------------------------------------------------------------------

## Usage Example

Here is a simple example demonstrating how to define a simulation
specification, simulate data, and plot the results:

``` r
library(fsimR)

# Define a linear model simulation specification
spec <- simSpec(
  formula = y ~ x1 + x2,
  distribution = list(x1 = "normal", x2 = "normal", y = "normal"),
  params = list(x1 = c(0,1), x2 = c(0,1), y = c(0,1))
)

# Simulate data for a linear model
sim_data <- simulate_LMdata(spec)

# Print a summary of simulated data
print(sim_data)

# Plot simulated data
plot(sim_data)
```

------------------------------------------------------------------------

## Functions & Documentation

The `fsimR` package provides the following functions. Internal functions
(starting with `.`) are not intended for direct use.

### Helper Functions

| Function Name | Documentation |
|----|----|
| `match.fun.allR` | [match.fun.allR.Rd](https://github.com/abhik-stat/fsimR/blob/master/man/match.fun.allR.Rd) |
| `distr_spec` | [distr_spec.Rd](https://github.com/abhik-stat/fsimR/blob/master/man/distr_spec.Rd) |
| `sim_spec` | [sim_spec.Rd](https://github.com/abhik-stat/fsimR/blob/master/man/sim_spec.Rd) |

### Exported Functions & Aliases

| Function Name | Documentation | Description |
|----|----|----|
| `distrSpec` | [distr_spec.Rd](https://github.com/abhik-stat/fsimR/blob/master/man/distr_spec.Rd) | User-friendly alias for `distr_spec`. |
| `simSpec` | [sim_spec.Rd](https://github.com/abhik-stat/fsimR/blob/master/man/sim_spec.Rd) | Exported simulation specification object. |
| `add_contamination` | [add_contamination.Rd](https://github.com/abhik-stat/fsimR/blob/master/man/add_contamination.Rd) | Adds contamination/noise to simulated data. |
| `simulate_IIDdata` | [simulate_IIDdata.Rd](https://github.com/abhik-stat/fsimR/blob/master/man/simulate_IIDdata.Rd) | Simulates IID data. |
| `simulate_LMdata` | [simulate_LMdata.Rd](https://github.com/abhik-stat/fsimR/blob/master/man/simulate_LMdata.Rd) | Simulates linear model data. |
| `simulate_LMMdata` | [simulate_LMMdata.Rd](https://github.com/abhik-stat/fsimR/blob/master/man/simulate_LMMdata.Rd) | Simulates linear mixed model data. |

### Plot Methods for Simulated Data

| Function Name | Documentation | Description |
|----|----|----|
| `plotSimData` | [plotSimData.Rd](https://github.com/abhik-stat/fsimR/blob/master/man/plotSimData.Rd) | Main plotting function for simulated data objects. |
| `plot.GLMdata` | [plotSimData.Rd](https://github.com/abhik-stat/fsimR/blob/master/man/plotSimData.Rd) | Plot method for GLM simulation output. |
| `plot.GLMMdata` | [plotSimData.Rd](https://github.com/abhik-stat/fsimR/blob/master/man/plotSimData.Rd) | Plot method for GLMM simulation output. |
| `plot.IIDdata` | [plotSimData.Rd](https://github.com/abhik-stat/fsimR/blob/master/man/plotSimData.Rd) | Plot method for IID simulation output. |
| `plot.LMdata` | [plotSimData.Rd](https://github.com/abhik-stat/fsimR/blob/master/man/plotSimData.Rd) | Plot method for LM simulation output. |
| `plot.LMMdata` | [plotSimData.Rd](https://github.com/abhik-stat/fsimR/blob/master/man/plotSimData.Rd) | Plot method for LMM simulation output. |
| `plot.list` | [plotSimData.Rd](https://github.com/abhik-stat/fsimR/blob/master/man/plotSimData.Rd) | Plot method for lists of simulated data. |

### Print Methods for Simulated Data

| Function Name | Documentation | Description |
|----|----|----|
| `printSimData` | [printSimData.Rd](https://github.com/abhik-stat/fsimR/blob/master/man/printSimData.Rd) | Main print method for simulated data objects. |
| `print.GLMMdata` | [printSimData.Rd](https://github.com/abhik-stat/fsimR/blob/master/man/printSimData.Rd) | Print method for GLMM simulation output. |
| `print.LMdata` | [printSimData.Rd](https://github.com/abhik-stat/fsimR/blob/master/man/printSimData.Rd) | Print method for LM simulation output. |
| `print.LMMdata` | [printSimData.Rd](https://github.com/abhik-stat/fsimR/blob/master/man/printSimData.Rd) | Print method for LMM simulation output. |

------------------------------------------------------------------------

## Citation

If you use `fsimR` in your research, please cite it:

``` r
citation("fsimR")
```

------------------------------------------------------------------------

## Contributing

Contributions are welcome! Please open issues or submit pull requests
via the [GitHub repository](https://github.com/abhik-stat/fsimR).

------------------------------------------------------------------------

## License

`fsimR` is licensed under the Apache License 2.0.
