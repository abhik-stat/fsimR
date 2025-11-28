# Enhanced Pairs Plot with histogram, smoothers, and correlation panels

Creates a pairs plot with:

- histograms on the diagonal

- smoothing curves in the upper panels

- correlation coefficients and p-values in the lower panels

Constant columns are removed automatically.

## Usage

``` r
.plot_pairs(mat, main_text = NULL, plot_args = list(), ...)
```

## Arguments

- mat:

  Numeric matrix or data frame to plot.

- main_text:

  Character: main title for the plot.

- plot_args:

  List of graphical parameters (colors, line widths, label offsets)

- ...:

  Additional arguments passed to `pairs` and panel functions
