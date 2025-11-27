# =========== Internal Helper Functions ================

#' @keywords internal
#' @title Plot Histogram With Optional Density Overlay
#' @description
#' Internal helper for plotting histograms with an optional kernel density curve.
#' For diagonal panels in pairs plots, counts and density are scaled to [0,1].
#'
#' Skips plotting if data has only one unique value.
#'
#' @param x_vec Numeric vector to plot.
#' @param label Character; label for both x- and y-axis (or diagonal panels).
#' @param is_diag Logical: whether plotting in a diagonal panel within pairs plot
#' @param main Logical; if TRUE, display the main title specified in `main_text`.
#' @param plot_args List of graphical  parameters (colors, breaks, line widths, etc.)
#' @param ... Additional graphical parameters passed to `hist()` or `lines()`
#'
.plot_hist_density <- function(x_vec, main_text = NULL, label = NULL,
                               is_diag = FALSE, main = TRUE,
                               plot_args = list(), ...) {

  if (length(unique(x_vec)) <= 1) {
    # Skip plotting constant data
    warning("Column has only one unique value; skipping histogram.")
    return(invisible(NULL))
  }

  if (!is_diag) {

    # Compute objects first
    h <- hist(x_vec, plot = FALSE, breaks = plot_args$breaks)
    d <- density(x_vec, na.rm = TRUE)

    # Determine y-limit based on both
    ymax <- max(c(h$density, d$y), na.rm = TRUE)

    # Now plot histogram using the adjusted ylim
    hist(x_vec, probability = TRUE, col = plot_args$hist_col,
         breaks = plot_args$breaks,
         main = if (main) main_text else "",
         cex.axis = plot_args$cex_axis,
         ylim = c(0, ymax),
         ...)

    # Add density curve
    lines(d$x, d$y, col = plot_args$dens_col, lwd = plot_args$dens_lwd, ...)

    # Add labels
    if (!is.null(label))
      title(xlab = label, ylab = label, line = plot_args$label_line)
  } else {
    # Diagonal histogram for pairs plot
    usr_old <- par("usr")
    on.exit(par(usr = usr_old), add = TRUE)
    par(usr = c(usr_old[1:2], 0, 1.1))

    h <- hist(x_vec, plot = FALSE, breaks = plot_args$breaks)
    yscaled <- h$counts / max(h$counts)
    rect(h$breaks[-length(h$breaks)], 0,
         h$breaks[-1], yscaled, col = plot_args$hist_col)

    if (length(unique(x_vec)) > 1) {
      d <- density(x_vec)
      lines(d$x, d$y / max(d$y), col = plot_args$dens_col, lwd = plot_args$dens_lwd)
    }
  }
}



#' @keywords internal
#'
#' @title Correlation Panel for Pairs Plots
#'
#' @description
#' Displays correlation coefficient `r`, p-value (`p`),
#' and significance stars in the selected panels of pairs plots,
#' using math expressions. Negative correlations are marked in red,
#' and background color is scaled by the absolute correlation.
#'
#' @param x,y Numeric vectors to plot.
#' @param digits Number of digits used for correlation coefficient and p-value rounding (default 4).
#' @param cex Character expansion factor for text (default 1.2).
#' @param low_col Background color for `r = 0` (default "#FFFFFF").
#' @param high_col Background color for `|r| = 1` (default "#2166AC").
#' @param ... Additional graphical arguments passed to text().
#'
.panel.cor.p <- function(x, y,
                         digits = 4,
                         cex = 1.2,
                         low_col = "#FFFFFF",
                         high_col = "#2166AC",
                         ...) {

  old_usr <- par("usr"); on.exit(par(usr = old_usr), add = TRUE)
  par(usr = c(0,1,0,1))

  # compute correlation and p-value
  ct <- suppressWarnings(cor.test(x, y))
  r <- round(ct$estimate, digits)
  p <- round(ct$p.value, digits)

  # determine significance stars
  stars <- if (p < .001) "***" else if (p < .01) "**" else if (p < .05) "*" else ""

  # background color scaled by |r|
  abs_r <- abs(r)
  bg_col <- grDevices::colorRampPalette(c(low_col, high_col))(100)[
    max(1, min(100, round(abs_r * 100)))
  ]
  rect(0, 0, 1, 1, col = bg_col, border = NA)


  # text color: red for negative r
  txt_col <- if (r < 0) "red" else "black"

  # math-mode r
  expr_r <- if (r < 0) {
    substitute(italic(r) == -V, list(V = abs(r)))
  } else {
    substitute(italic(r) == V, list(V = r))
  }

  # math-mode p
  expr_p <- if (p == 0) {
    bquote(italic(p) == 0)
  } else if (p < 1e-4) {
    sci <- format(p, scientific = TRUE)
    parts <- strsplit(sci, "e")[[1]]
    mant <- as.numeric(parts[1])
    expo <- as.integer(parts[2])
    bquote(italic(p) == .(mant) %*% 10^.(expo))
  } else {
    bquote(italic(p) == .(p))
  }

  # draw r
  text(0.5, 0.65, expr_r, cex = cex, col = txt_col, ...)

  # draw stars properly to the right of r
  if (stars != "") {
    # width of r in user coordinates
    str_w <- strwidth(as.expression(expr_r), cex = cex)
    text(0.5 + str_w/2 + 0.01, 0.65, stars, cex = cex*0.8,
         col = txt_col, adj = 0, ...)
  }

  # draw p
  text(0.5, 0.35, expr_p, cex = cex, col = txt_col, ...)

  invisible(NULL)
}



#' @keywords internal
#' @title Enhanced Pairs Plot with histogram, smoothers, and correlation panels
#' @description
#' Creates a pairs plot with:
#' - histograms on the diagonal
#' - smoothing curves in the upper panels
#' - correlation coefficients and p-values in the lower panels
#'
#' Constant columns are removed automatically.
#'
#' @param mat Numeric matrix or data frame to plot.
#' @param main_text Character: main title for the plot.
#' @param plot_args List of graphical parameters (colors, line widths, label offsets)
#' @param ... Additional arguments passed to `pairs` and panel functions
#'
.plot_pairs <- function(mat, main_text = NULL, plot_args = list(), ...) {

  mat <- as.data.frame(mat)

  # Remove constant columns
  constant_cols <- sapply(mat, function(col) length(unique(col)) <= 1)
  if (all(constant_cols)) {
    warning("All columns are constant; skipping pairs plot.")
    return(invisible(NULL))
  }
  if (any(constant_cols)) {
    mat <- mat[, !constant_cols, drop = FALSE]
    # warning("Removed constant columns from pairs plot.")
  }
  p <- ncol(mat)

  if (ncol(mat) == 1) {
    # If only one column, just plot histogram
    .plot_hist_density(mat[[1]], main_text, label = colnames(mat),
                       is_diag = FALSE, plot_args = plot_args, ...)
    return()
  }

  # set margins --
  # base outer margin (in lines)
  base_margin <- plot_args$label_line  # adjust this if needed

  # scale with number of columns and axis text size
  oma_bottom <- base_margin + plot_args$cex_axis
  oma_left   <- base_margin + plot_args$cex_axis
  oma_top    <- base_margin
  oma_right  <- base_margin

  par(oma = c(oma_bottom, oma_left, oma_top, oma_right), mar=c(1,1,1,1))

  upper_panel <- function(x, y, ...) {
    usr_old <- par("usr"); on.exit(par(usr = usr_old), add = TRUE)
    par(usr = c(range(x, na.rm = TRUE), range(y, na.rm = TRUE)))

    # panel.smooth ignores axes parameter, so call without it
    panel.smooth(x, y, col.smooth = plot_args$smooth_col,
                 lwd = plot_args$smooth_lwd)

    # manually remove axes using axis()/box()
    box()             # optional if you want a border
    row_idx <- par("mfg")[1]
    if (row_idx == 1) axis(side = 3, cex.axis = plot_args$cex_axis)

    col_idx <- par("mfg")[2]
    if (col_idx == p) axis(side = 4, cex.axis = plot_args$cex_axis)
  }

  # Lower panel: correlation + p-value, right ticks in last column
  lower_panel <- function(x, y, ...) {
    usr_old <- par("usr"); on.exit(par(usr = usr_old), add = TRUE)
    par(usr = c(0, 1, 0, 1))

    .panel.cor.p(x, y,
                 digits = plot_args$digits,
                 cex = plot_args$cex_text,
                 low_col = plot_args$low_col,
                 high_col = plot_args$high_col)

  }

  # Standard pairs plot, upper panel uses .panel.cor.p
  pairs(mat,
        labels = rep("", ncol(mat)),
        diag.panel = function(x, ...) .plot_hist_density(x, is_diag = TRUE, plot_args = plot_args, ...),
        upper.panel = upper_panel,
        lower.panel = lower_panel,
        xaxt="n", yaxt="n",
        main = if (!is.null(main_text)) main_text else "",
        ...)

  # Add margin labels, shifted outward to avoid overlap with axis ticks
  label_names <- colnames(mat)
  # compute vertical positions
  panel_centers <- seq(from = 1 - 1/(2*p), to = 1/(2*p), length.out = p)

  for (i in 1:p) {
    # horizontal labels
    mtext(label_names[i], side = 1, line = base_margin-1,
          at = (i - 0.5)/p - (0.01*i), outer = TRUE, cex = plot_args$cex_axis)

    # vertical labels
    mtext(label_names[i], side = 2, line = base_margin-1,
          at = panel_centers[i] - (0.1/i), outer = TRUE, cex = plot_args$cex_axis)
  }
}





# =========== Common documentation for all plot methods ================
#' @title Plot Simulated Data
#'
#' @description
#' Generates visualizations for simulated datasets, supporting multiple types:
#'
#' - **Mixed Models** (of class `LMMdata` or `GLMMdata`)
#' - **Regression Models** (of class `LMdata` or `GLMdata`)
#' - **IID samples** (of class `IIDdata`)
#'
#' The function also works for named lists containing:
#' - response vectors as `y`
#' - fixed-effects design matrix as `X`
#' - optional random-effects design matrix as `Z`
#' - optional random-effects realizations as `RE`
#'
#' Numeric vectors or matrices are treated as IID samples
#' where columns correspond to variables, rows to observations.
#'
#' For plot types and layout, see Details.
#'
#' @param object Numeric vector, matrix, named list, or simulation object
#'    (of class `LMMdata`, `GLMMdata`, `LMdata`, `GLMdata`, `IIDdata`).
#'    If a list, it must include the components required for the selected `type`.
#'
#' @param replication Integer; which replication (dataset) to plot
#'    (LMM/GLMM, LM/GLM only; default `1`).
#' @param iteration Integer; which iteration within the replication to plot
#'    (LMM/GLMM, LM/GLM only; default `1`).
#'
#' @param type Character vector specifying which components to plot:
#'   - `"y"`: response vector
#'   - `"X"`: fixed-effects or design matrix
#'   - `"Z"`: random-effects design matrix (LMM/GLMM only)
#'   - `"RE"`: random effects (LMM/GLMM only)
#'   - `"all"`: all available components (default)
#'
#'  Ignored for IID samples, where all columns are always plotted.
#'
#' @param main Logical; if `TRUE`, descriptive titles are shown (default `TRUE`).
#'
#' @param plot_args Optional list of graphical parameters to customize
#' appearance of the plots. Supported arguments:
#'   \describe{
#'     \item{breaks}{Number of histogram bins (default `10`)}
#'     \item{hist_col}{Histogram fill color (default `"lightsteelblue"`)}
#'     \item{dens_col}{Line color for density curve (default `"firebrick"`)}
#'     \item{dens_lwd}{Line width for density curve (default `2`)}
#'     \item{smooth_col}{Smoothing line color (default `"blue"`)}
#'     \item{smooth_lwd}{Smoothing line width (default `2`)}
#'     \item{cex_axis}{Character expansion for axis labels (default `1.4`)}
#'     \item{label_line}{Line offset for margin labels (default `0`)}
#'     \item{low_col}{Background color for `r = 0` in correlation panels (default `"#FFFFFF"`)}
#'     \item{high_col}{Background color for `|r| = 1` in correlation panels (default `"#2166AC"`)}
#'     \item{cex_text}{Character expansion for correlation/p-value text (default `1.4`)}
#'     \item{digits}{Number of digits for correlation/p-value (default `4`)}.
#'   }
#'
#'
#' @param ... Additional graphical arguments passed to
#' the standard plotting functions, such as `hist`, `lines`, `pairs` and `text`
#' (via internal functions `.plot_hist_density` and `.plot_pairs()`).
#'
#' @details
#' The plotting behavior depends on the class and selected components (in `type`):
#'
#' 1. **Histograms with density overlay** for response vectors (`y`)
#'    or any numeric vector (treated as a univariate IID samples).
#'    Skipped for constant vectors.
#'
#' 2. **Enhanced pairs plots** for matrices (`X`, `Z`, `RE`) or multivariate IID samples:
#'    - Diagonal panels: histograms with density scaled to [0, 1].
#'    - Upper panels: scatter plots with smoothing lines ([loess]).
#'    - Lower panels: correlation coefficients, p-values,
#'      and significance stars (`*`, `**`, `***`).
#'      Background color reflects the absolute correlation;
#'      negative correlations appear in red.
#'    - Constant columns are automatically removed.
#'
#' @return
#' Invisibly returns `NULL`. Generates plots for visual inspection.
#'
#' @examples
#' \dontrun{
#' # Mixed Model
#' sim_settings <- sim_spec(p = 3)
#' sim_data <- simulate_LMMdata(c(2,3), sim_settings)
#' plot(sim_data, type = c("y", "RE"))
#'
#' # Regression Model
#' sim_data <- simulate_LMdata(c(2,3), sim_settings)
#' plot(sim_data, type = "X")
#'
#' # IID data
#' sim_data <- simulate_IIDdata(n = 100)
#' plot(sim_data)
#'
#' # Named list wrapper
#' mylist <- list(y = rnorm(50), X = matrix(rnorm(100), ncol=2))
#' plot(mylist)
#'
#' }
#'
#' @importFrom graphics axis box hist lines mtext pairs panel.smooth par plot.new rect strwidth text title
#' @name plotSimData
NULL



# =========== Plot method for LMMdata/GLMMdata objects ================
#' @rdname plotSimData
#' @method plot LMMdata
#' @export
plot.LMMdata <- function(object, replication = 1, iteration = 1,
                         type = c("all", "y", "X", "Z", "RE"),
                         main = TRUE,
                         plot_args = list(), ...) {

  # Default plot arguments
  default_args <- list(
    breaks = 10,
    hist_col = "lightsteelblue",
    dens_col = "firebrick",
    dens_lwd = 2,
    smooth_col = "blue",
    smooth_lwd = 2,
    cex_axis = 1.4,
    label_line = 0,
    low_col = "#FFFFFF",
    high_col = "#2166AC",
    cex_text = 1.4,
    digits = 4
  )

  # Merge user-provided args
  plot_args <- modifyList(default_args, plot_args)

  # Save current graphics parameters and restore at exit
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))

  par(mfrow = c(1,1))
  plot.new()
  par(mai = c(1, 1, 1, 1))

  type <- match.arg(type, several.ok = TRUE)
  plot_types <- if ("all" %in% type) c("y", "X", "Z", "RE") else type

  # Validate required components
  required <- c("y", "X", "Z", "RE")
  missing_comp <- setdiff(required, names(object))
  if (length(missing_comp) > 0) {
    stop("Object is missing required components: ", paste(missing_comp, collapse = ", "))
  }

  # Ensure components are lists
  object[c("y", "X", "Z", "RE")] <- lapply(object[c("y", "X", "Z", "RE")],
                                                  function(x) if (!is.list(x)) list(x) else x)

  # Determine replication and iteration
  n_sim <- length(object$y)
  n_iter <- if (is.matrix(object$y[[1]])) ncol(object$y[[1]]) else length(object$y[[1]])

  if (replication > n_sim || replication < 1) stop("Replication index out of range")
  if (iteration > n_iter || iteration < 1) stop("Iteration index out of range")

  y_sel <- if (n_iter == 1) object$y[[replication]] else object$y[[replication]][, iteration, drop = TRUE]
  X_sel <- object$X[[replication]]
  Z_sel <- object$Z[[replication]]

  RE_raw <- object$RE[[replication]]
  RE_sel <- if (is.matrix(RE_raw)) {
    if (n_iter > 1) RE_raw[, iteration, drop = FALSE] else RE_raw
  } else {
    if (n_iter > 1) RE_raw[[iteration]] else RE_raw
  }

  if (is.vector(RE_sel)) RE_sel <- matrix(RE_sel, ncol = 1)
  colnames(RE_sel) <- paste0("RE", seq_len(ncol(RE_sel)))

  # plot selected components
  for (tp in plot_types) {
    switch(tp,
           y  = .plot_hist_density(y_sel, "Response vector (y)", label = "", is_diag = FALSE, plot_args = plot_args, ...),
           X  = .plot_pairs(X_sel, "Fixed-effects design (X)", plot_args = plot_args, ...),
           Z  = .plot_pairs(Z_sel, "Random-effects design (Z)", plot_args = plot_args, ...),
           RE = .plot_pairs(RE_sel, "Random effects (RE)", plot_args = plot_args, ...)
    )
  }

  invisible(NULL)
}



# GLMMdata uses same function as LMMdata
#' @rdname plotSimData
#' @method plot GLMMdata
#' @export
plot.GLMMdata <- function(object, ...) plot.LMMdata(object, ...)



# =========== Plot method for LMdata/GLMdata objects ================
#' @rdname plotSimData
#' @method plot LMdata
#' @export
plot.LMdata <- function(object, replication = 1, iteration = 1,
                         type = c("all", "y", "X"),
                         main = TRUE,
                         plot_args = list(), ...) {

  # Default plot arguments
  default_args <- list(
    breaks = 10,
    hist_col = "lightsteelblue",
    dens_col = "firebrick",
    dens_lwd = 2,
    smooth_col = "blue",
    smooth_lwd = 2,
    cex_axis = 1.4,
    label_line = 0,
    low_col = "#FFFFFF",
    high_col = "#2166AC",
    cex_text = 1.4,
    digits = 4
  )

  # Merge user-provided args
  plot_args <- modifyList(default_args, plot_args)

  # Save current graphics parameters and restore at exit
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))

  par(mfrow = c(1,1))
  plot.new()
  par(mai = c(1, 1, 1, 1))

  type <- match.arg(type, several.ok = TRUE)
  plot_types <- if ("all" %in% type) c("y", "X") else type

  # Validate required components
  required <- c("y", "X")
  missing_comp <- setdiff(required, names(object))
  if (length(missing_comp) > 0) {
    stop("Object is missing required components: ", paste(missing_comp, collapse = ", "))
  }

  # Ensure components are lists
  object[c("y", "X")] <- lapply(object[c("y", "X")],
                                function(x) if (!is.list(x)) list(x) else x)

  # Determine replication and iteration
  n_sim <- length(object$y)
  n_iter <- if (is.matrix(object$y[[1]])) ncol(object$y[[1]]) else length(object$y[[1]])

  if (replication > n_sim || replication < 1) stop("Replication index out of range")
  if (iteration > n_iter || iteration < 1) stop("Iteration index out of range")

  y_sel <- if (n_iter == 1) object$y[[replication]] else object$y[[replication]][, iteration, drop = TRUE]
  X_sel <- object$X[[replication]]

  # plot selected components
  for (tp in plot_types) {
    switch(tp,
           y  = .plot_hist_density(y_sel, "Response vector (y)", label = "", is_diag = FALSE, plot_args = plot_args, ...),
           X  = .plot_pairs(X_sel, "Design matrix (X)", plot_args = plot_args, ...),
    )
  }

  invisible(NULL)
}


# GLMdata uses same function as LMdata
#' @rdname plotSimData
#' @method plot GLMdata
#' @export
plot.GLMdata <- function(object, ...) plot.LMdata(object, ...)




# =========== Plot method for IIDdata objects ================
#' @rdname plotSimData
#' @method plot IIDdata
#' @param object Numeric matrix or data frame of IID samples.
#' @export
plot.IIDdata <- function(data, main = TRUE, plot_args = list(), ...) {

  # format data suitably
  data <- as.matrix(data)
  n <- nrow(data)
  d <- ncol(data)
  if(is.null(colnames(data))) colnames(data) <- paste0("X", 1:d)

  # Default plot arguments
  default_args <- list(
    breaks = 10,
    hist_col = "lightsteelblue",
    dens_col = "firebrick",
    dens_lwd = 2,
    smooth_col = "blue",
    smooth_lwd = 2,
    cex_axis = 1.4,
    label_line = 0,
    low_col = "#FFFFFF",
    high_col = "#2166AC",
    cex_text = 1.4,
    digits = 4
  )

  # Merge user-provided args
  plot_args <- modifyList(default_args, plot_args)

  # Save current graphics parameters and restore at exit
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))

  par(mfrow = c(1,1))
  plot.new()
  par(mai = c(1, 1, 1, 1))


  # .plot_hist_density(y_sel, "Response vector (y)", label = "", is_diag = FALSE, plot_args = plot_args, ...),
  .plot_pairs(data, "Sample IID data (X)", plot_args = plot_args, ...)

  invisible(NULL)
}



# =========== Plot method for Other classes ================

# List wrapper
#' @rdname plotSimData
#' @method plot list
#' @param object List containing components y, X, Z, RE (LMM/GLMM) or y, X (LM/GLM).
#' @export
plot.list <- function(object, ...) {
  if (all(c("y","X","Z","RE") %in% names(object))) {
    plot.LMMdata(object, ...)
  } else if (all(c("y","X") %in% names(object))) {
    plot.LMdata(object, ...)
  } else {
    stop("Cannot determine plot method for this list. Ensure it contains valid components.")
  }
}
