# fsimR

This R package provides a unified flexible framework for simulating data
from a wide range of structured statistical models with ease of control.
The package support generation of independent and identically
distributed (IID) data under a variety of distributions, including
user-defined generators and copula-based multivariate sampling. It also
includes functions for generating data from different regression and
mixed-effects models, with customizable covariate distributions and
error structures. It also allows summarizing, printing and visualizing
simulated data, as well as adding artificial contamination or censoring
of various types.

The package is designed to support reproducible simulation studies,
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

## Functions & Documentation

The `fsimR` package provides the following functions.

### Helper Functions

| Function Name    | Documentation                                                                           | Description                                                |
|------------------|-----------------------------------------------------------------------------------------|------------------------------------------------------------|
| `distr_spec`     | [`distr_spec.Rd`](https://abhik-stat.github.io/fsimR/reference/distr_spec.html)         | Distribution Specification Constructor                     |
| `sim_spec`       | [`sim_spec.Rd`](https://abhik-stat.github.io/fsimR/reference/sim_spec.html)             | Simulation Specification Constructor                       |
| `match.fun.allR` | [`match.fun.allR.Rd`](https://abhik-stat.github.io/fsimR/reference/match.fun.allR.html) | Find a Function Across Base R and All Installed R Packages |

### Main Simulation Functions

| Function Name       | Documentation                                                                               | Description                                                                       |
|---------------------|---------------------------------------------------------------------------------------------|-----------------------------------------------------------------------------------|
| `add_contamination` | [add_contamination.Rd](https://abhik-stat.github.io/fsimR/reference/add_contamination.html) | Adds contamination/noise to simulated data.                                       |
| `simulate_IIDdata`  | [simulate_IIDdata.Rd](https://abhik-stat.github.io/fsimR/reference/simulate_IIDdata.html)   | Simulate IID Data from a Wide Range of Univariate and Multivariate Distributions. |
| `simulate_LMdata`   | [simulate_LMdata.Rd](https://abhik-stat.github.io/fsimR/reference/simulate_LMdata.html)     | Simulate Data from a Linear Regression Model                                      |
| `simulate_LMMdata`  | [simulate_LMMdata.Rd](https://abhik-stat.github.io/fsimR/reference/simulate_LMMdata.html)   | Simulate Data from a Linear Mixed Model (LMM)                                     |

### S3 Generic Methods for Simulated Data

| Function Name | Documentation                                                                         | Description                             |
|---------------|---------------------------------------------------------------------------------------|-----------------------------------------|
| `as.<class>`  | [asClass.Rd](https://abhik-stat.github.io/fsimR/reference/asClass.html)               | Coercion into package-specific classes. |
| `plot`        | [plotSimData.Rd](https://abhik-stat.github.io/fsimR/reference/plotSimData.html)       | Plot simulated data objects.            |
| `print`       | [printSimData.Rd](https://abhik-stat.github.io/fsimR/reference/printSimData.html)     | Print simulated data objects.           |
| `summary`     | [summarySimData.Rd](https://abhik-stat.github.io/fsimR/reference/summarySimData.html) | Summarize simulated data objects.       |

------------------------------------------------------------------------

## Usage Example

Here is a simple example demonstrating how to define a simulation
specification, simulate data, and plot the results:

``` r
library(fsimR)

# Define a linear model simulation specification with Non-Gaussian error and multivariate predictors
sim_settings = sim_spec(n_subj = 15)

distr_settings <- list(
error_distr = list(distr_name = "t", distr_params = distr_spec(df = 4)),
X_distr    = list(distr_name = "mvnorm", distr_params = list(dim = 4, sigma = diag(4)))
)

# Simulate data for a linear model
sim_data <- simulate_LMdata(n_rep = 1, sim_settings = sim_settings, distr_settings = distr_settings)

# Print a summary of simulated data
print(sim_data)

# Plot simulated data
plot(sim_data)
```

More complex and realistic examples can be found in:

- [Multivariate
  Data](https://abhik-stat.github.io/fsimR/articles/examples_multivar.html)
- [Linear Regression
  Models](https://abhik-stat.github.io/fsimR/articles/examples_LM.html)
- [Linear Mixed
  Models](https://abhik-stat.github.io/fsimR/articles/examples_LMM.html)

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
