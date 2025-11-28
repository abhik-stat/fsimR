# Find a Function Across Base R and All Installed R Packages

Extends the behavior of
[`match.fun`](https://rdrr.io/r/base/match.fun.html) by searching not
only the current environment and base R, but also all installed packages
(namespaces). It additionally supports the explicit `"pkg::fun"` syntax
to retrieve a function from a specific package namespace.

## Usage

``` r
match.fun.allR(FUN, descend = TRUE, list_all = FALSE)
```

## Arguments

- FUN:

  A function, character string, or symbol naming a function. May also be
  of the form `"pkg::fun"` to explicitly reference a function within a
  particular package.

- descend:

  Logical; if `TRUE` (default), the lookup in the current environment
  behaves like [`match.fun`](https://rdrr.io/r/base/match.fun.html) and
  descends the inheritance tree when resolving S3/S4 methods. If
  `FALSE`, only non-generic function objects are accepted.

- list_all:

  Logical; if `TRUE`, returns a list of all matching function objects
  found across installed R packages. If `FALSE` (default), returns only
  the first match found. See @details

## Value

A function (if a single match is found or `list_all = FALSE`) or a named
list of functions (if multiple matches are found and `list_all = TRUE`).
If a list is returned, the names correspond to the package or
environment where each function is found. An informative error is raised
if the function cannot be found.

## Details

It searches for `FUN` in the following order:

1.  If `FUN` is already a function, it is returned unchanged.

2.  If `FUN` is of the form `"pkg::fun"`, then `pkg` must be installed;
    the function is retrieved directly from `asNamespace(pkg)`.

3.  The current calling environment is searched (similar to
    [`match.fun`](https://rdrr.io/r/base/match.fun.html)). If
    `descend = TRUE`, method dispatch is allowed; otherwise only plain
    functions are accepted.

4.  All installed package namespaces are searched (both exports and
    internal functions). The first match found—according to the
    alphabetical order of installed packages—is returned if
    `list_all = FALSE`. If `list_all = TRUE`, all matches are returned.
    *Note:* This may yield surprising matches when multiple packages
    define functions with identical names and `list_all = FALSE`; in
    such cases, the explicit `"pkg::fun"` syntax is recommended.

## See also

[`match.fun`](https://rdrr.io/r/base/match.fun.html) for the base R
version, `getNamespace`, `getExportedValue`

## Examples

``` r
# Base R function
f1 <- match.fun.allR("rnorm")
f1(5)
#> [1] -0.93199903 -1.42789198  0.97576509 -1.54634119  0.01770348

# Function from installed package (if sn is installed)
if (requireNamespace("sn", quietly = TRUE)) {
  f2 <- match.fun.allR("rsn")       # search all installed packages
  f3 <- match.fun.allR("sn::rsn")   # explicit package reference
  f2(5); f3(5)
}
#> Registered S3 methods overwritten by 'RcppEigen':
#>   method               from         
#>   predict.fastLm       RcppArmadillo
#>   print.fastLm         RcppArmadillo
#>   summary.fastLm       RcppArmadillo
#>   print.summary.fastLm RcppArmadillo
#> Registered S3 method overwritten by 'fit.models':
#>   method       from 
#>   vcov.default Hmisc
#> Registered S3 method overwritten by 'resample':
#>   method         from  
#>   print.resample modelr
#> [1] -0.7997961 -1.0250161  0.3817112  1.9005427  0.3804597
#> attr(,"family")
#> [1] "SN"
#> attr(,"parameters")
#> [1] 0 1 0 0

f4 <- match.fun.allR("rnorm", list_all = TRUE)
class(f4)
#> [1] "function"
```
