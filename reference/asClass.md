# Coercion helpers for package-specific S3 classes

These functions provide lightweight coercion into the package’s custom
classes, such as `covMat`, `IIDdata`, `LMdata`, `LMMdata`, etc. Each
coercion function has the form: `as.<ClassName>(object)` and simply
prepends the corresponding class name to the object's existing class
vector. This enables S3 dispatch while preserving the underlying data
structure.

## Usage

``` r
as.covMat(object)

as.IIDdata(object)

as.LMMdata(object)

as.LMdata(object)
```

## Arguments

- object:

  Any R object.

## Value

The same object, with the target class prepended to its class vector (if
the object is compatible with the class).
