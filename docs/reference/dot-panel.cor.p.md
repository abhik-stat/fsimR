# Correlation Panel for Pairs Plots

Displays correlation coefficient `r`, p-value (`p`), and significance
stars in the selected panels of pairs plots, using math expressions.
Negative correlations are marked in red, and background color is scaled
by the absolute correlation.

## Usage

``` r
.panel.cor.p(
  x,
  y,
  digits = 4,
  cex = 1.2,
  low_col = "#FFFFFF",
  high_col = "#2166AC",
  ...
)
```

## Arguments

- x, y:

  Numeric vectors to plot.

- digits:

  Number of digits used for correlation coefficient and p-value rounding
  (default 4).

- cex:

  Character expansion factor for text (default 1.2).

- low_col:

  Background color for `r = 0` (default "#FFFFFF").

- high_col:

  Background color for `|r| = 1` (default "#2166AC").

- ...:

  Additional graphical arguments passed to text().
