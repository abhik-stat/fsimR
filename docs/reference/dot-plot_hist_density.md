# Plot Histogram With Optional Density Overlay

Internal helper for plotting histograms with an optional kernel density
curve. For diagonal panels in pairs plots, counts and density are scaled
to 0,1.

Skips plotting if data has only one unique value.

## Usage

``` r
.plot_hist_density(
  x_vec,
  main_text = NULL,
  label = NULL,
  is_diag = FALSE,
  main = TRUE,
  plot_args = list(),
  ...
)
```

## Arguments

- x_vec:

  Numeric vector to plot.

- label:

  Character; label for both x- and y-axis (or diagonal panels).

- is_diag:

  Logical: whether plotting in a diagonal panel within pairs plot

- main:

  Logical; if TRUE, display the main title specified in `main_text`.

- plot_args:

  List of graphical parameters (colors, breaks, line widths, etc.)

- ...:

  Additional graphical parameters passed to
  [`hist()`](https://rdrr.io/r/graphics/hist.html) or
  [`lines()`](https://rdrr.io/r/graphics/lines.html)
