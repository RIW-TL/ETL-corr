# ============================================================
# Non-Gaussian simulation: covariance, correlation and variance
# Requires f1, f2, f3, f4, f5 in the R environment.
# No estimation results are hard-coded.
# ============================================================

library(ggplot2)
library(ggpubr)
library(latex2exp)

# ------------------------------------------------------------
# 1. Data and x-axis
# ------------------------------------------------------------
fits <- list(f1, f2, f3, f4, f5)

# If known, replace NULL with the five actual tau values,
# in the order corresponding to f1, ..., f5.
# Otherwise, the x-axis displays f1, ..., f5.
tau_values <- seq(0,0.8,0.2)

if (is.null(tau_values)) {
  x_values <- seq_along(fits)
  x_labels <- paste0("f", seq_along(fits))
  x_axis_label <- "Non-Gaussian setting"
} else {
  stopifnot(is.numeric(tau_values), length(tau_values) == length(fits),
            !anyNA(tau_values), !anyDuplicated(tau_values))
  x_values <- tau_values
  x_labels <- as.character(tau_values)
  x_axis_label <- TeX("$\\tau$")
}

cov_methods <- c("ETL-cov", "ETL-Scov", "Target-cov",
                 "Target-Scov", "Pool-cov")
var_methods <- c("ETL-var", "Target-var", "Pool-var")

extract_losses <- function(field, methods) {
  out <- lapply(fits, function(fit) {
    vals <- fit[[field]]
    if (is.null(vals) || !all(methods %in% names(vals))) {
      stop(sprintf("Missing or unnamed methods in '%s'", field))
    }
    as.numeric(unlist(vals[methods], use.names = FALSE))
  })
  mat <- do.call(rbind, out)
  colnames(mat) <- methods
  mat
}

cov_mat <- extract_losses("mean_cov_loss", cov_methods)
cor_mat <- extract_losses("mean_cor_loss", cov_methods)
var_mat <- extract_losses("mean_var_loss", var_methods)

# Use the same five legend entries in all three panels.
# ETL-var corresponds to the diagonal estimates of ETL-cov/ETL-Scov.
# Target-var corresponds to those of Target-cov/Target-Scov.
legend_order <- c("ETL-cov", "ETL-Scov", "Pool-cov",
                  "Target-cov", "Target-Scov")
color_Values <- c(
  "ETL-cov" = "red", "ETL-Scov" = "orange",
  "Pool-cov" = "purple", "Target-cov" = "black",
  "Target-Scov" = "blue"
)

make_long <- function(mat, mapped_methods = colnames(mat)) {
  stopifnot(length(mapped_methods) == ncol(mat))
  data.frame(
    x = rep(x_values, times = ncol(mat)),
    method = factor(rep(mapped_methods, each = nrow(mat)),
                    levels = legend_order),
    loss = as.vector(mat)
  )
}

cov_df <- make_long(cov_mat)
cor_df <- make_long(cor_mat)
var_df <- make_long(var_mat, c("ETL-cov", "Target-cov", "Pool-cov"))

# ------------------------------------------------------------
# 2. Piecewise compressed y scale (first two panels only)
# ------------------------------------------------------------
compress_interval_trans <- function(from = 0.5, to = 1,
                                    compression = 5) {
  stopifnot(from < to, compression > 1)
  transform <- function(x) {
    y <- x
    middle <- which(x > from & x < to)
    upper <- which(x >= to)
    y[middle] <- from + (x[middle] - from) / compression
    y[upper] <- from + (to - from) / compression + (x[upper] - to)
    y
  }
  inverse <- function(x) {
    y <- x
    transformed_to <- from + (to - from) / compression
    middle <- which(x > from & x < transformed_to)
    upper <- which(x >= transformed_to)
    y[middle] <- from + (x[middle] - from) * compression
    y[upper] <- to + (x[upper] - transformed_to)
    y
  }
  scales::trans_new("compressed_interval", transform, inverse,
                    domain = c(-Inf, Inf))
}

y_trans <- compress_interval_trans(from = 0.5, to = 1.9,
                                   compression = 10)

# Five major ticks per panel; third panel remains linear.
y_scale_cov <- scale_y_continuous(
  trans = y_trans,
  breaks = c(0.10, 0.30, 0.50, 1.80),
  limits = c(0.05, 1.90),
  labels = function(x) formatC(x, format = "f", digits = 2),
  expand = expansion(mult = c(0.01, 0.015)),
  minor_breaks = NULL
)

y_scale_cor <- scale_y_continuous(
  trans = y_trans,
  breaks = c(0.10, 0.30, 0.50, 1.60),
  limits = c(0.05, 1.65),
  labels = function(x) formatC(x, format = "f", digits = 2),
  expand = expansion(mult = c(0.01, 0.015)),
  minor_breaks = NULL
)

y_scale_var <- scale_y_continuous(
  breaks = c(1.5, 2.5, 3.5, 4.5, 5.5),
  limits = c(1.5, 5.55),
  labels = function(x) formatC(x, format = "f", digits = 1),
  expand = expansion(mult = c(0.01, 0.015)),
  minor_breaks = NULL
)

# ------------------------------------------------------------
# 3. Three panels
# ------------------------------------------------------------
make_panel <- function(dat, panel_title, y_axis_label, y_scale,
                       show_legend = TRUE) {
  dat$panel <- panel_title
  ggplot(dat, aes(x = x, y = loss, colour = method, group = method)) +
    geom_line(linewidth = 1, show.legend = show_legend) +
    geom_point(size = 2, show.legend = show_legend) +
    scale_colour_manual(
      name = NULL, breaks = legend_order,
      limits = legend_order, values = color_Values,
      drop = FALSE
    ) +
    scale_x_continuous(breaks = x_values, labels = x_labels) +
    y_scale +
    labs(x = x_axis_label, y = y_axis_label) +
    theme_bw() +
    theme(
      plot.margin = margin(5, 6, 5, 5, unit = "pt"),
      panel.grid.minor = element_blank(),
      legend.position = "top"
    ) +
    facet_grid(. ~ panel)
}

k1 <- make_panel(
  cov_df, "Covariance",
  TeX("$\\|\\hat{\\Sigma}-\\Sigma^{(0)}\\|_F^2/p$"),
  y_scale_cov
)

k2 <- make_panel(
  cor_df, "Correlation",
  TeX("$\\|\\hat{R}-R^{(0)}\\|_F^2/p$"),
  y_scale_cor
)

k3 <- make_panel(
  var_df, "Marginal variances",
  TeX("$\\|\\hat{\\Lambda}-\\Lambda^{(0)}\\|_F^2$"),
  y_scale_var,
  show_legend = FALSE
)

p_NG <- ggarrange(
  k1, k2, k3,
  ncol = 3, nrow = 1,
  common.legend = TRUE,
  legend = "top"
)
print(p_NG)

# ggsave("nonGaussian-ETL.pdf", p_NG, width = 14, height = 4.5)
# ggsave("nonGaussian-ETL.png", p_NG, width = 14, height = 4.5, dpi = 300)
