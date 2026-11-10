library(ggplot2)
library(dplyr)
library(tidyr)
library(latex2exp)
library(cowplot)
library(scales)

#============================================================
# 0. Compression transformation
#============================================================

compress_interval_trans <- function(
    from,
    to,
    compression = 5) {
  
  stopifnot(
    to > from,
    compression > 1
  )
  
  trans_new(
    name = paste0(
      "compress_", from, "_", to, "_", compression
    ),
    
    transform = function(x) {
      y <- x
      
      id_mid <- is.finite(x) & x > from & x <= to
      id_high <- is.finite(x) & x > to
      
      y[id_mid] <-
        from + (x[id_mid] - from) / compression
      
      y[id_high] <-
        from + (to - from) / compression +
        (x[id_high] - to)
      
      y
    },
    
    inverse = function(y) {
      x <- y
      
      to_new <- from + (to - from) / compression
      
      id_mid <- is.finite(y) & y > from & y <= to_new
      id_high <- is.finite(y) & y > to_new
      
      x[id_mid] <-
        from + compression * (y[id_mid] - from)
      
      x[id_high] <-
        to + (y[id_high] - to_new)
      
      x
    }
  )
}


#============================================================
# 1. Basic settings
#============================================================

eta_vec <- c(0, 0.2, 0.4, 0.6, 0.8)

res_list <- list(
  aa1, aa2, aa3, aa4, aa5
)

# Methods displayed in the common legend

breaks.new <- c(
  "ETL-cov",
  "ETL-Scov",
  "Target-cov",
  "Target-Scov",
  "Pool-cov"
)

# Colors

color_Values <- c(
  "ETL-cov"     = "red",
  "ETL-Scov"    = "orange",
  "Target-cov"  = "black",
  "Target-Scov" = "blue",
  "Pool-cov"    = "purple",
  "ETL-var"     = "red",
  "Target-var"  = "black",
  "Pool-var"    = "purple"
)

# Y-axis breaks

cov_breaks <- c(0.10, 0.35, 0.60, 0.85, 1.10)
cor_breaks <- c(0.10, 0.30, 0.50, 0.70, 0.90)
var_breaks <- c(1.5, 2.0, 2.5, 6.0, 10.0)

# Compression settings

compress_from <- 2.5
compress_to   <- 12
compression   <- 5


#============================================================
# 2. Data extraction
#============================================================

make_plot_data <- function(
    loss_name,
    methods) {
  
  stopifnot(length(res_list) == length(eta_vec))
  
  out <- do.call(
    rbind,
    Map(
      function(x, eta) {
        
        values <- x[[loss_name]]
        
        if (is.null(values)) {
          stop(paste("Missing loss:", loss_name))
        }
        
        missing_methods <- setdiff(
          methods,
          names(values)
        )
        
        if (length(missing_methods) > 0) {
          stop(
            paste(
              "Missing methods:",
              paste(missing_methods, collapse = ", ")
            )
          )
        }
        
        data.frame(
          eta = eta,
          Method = methods,
          Loss = as.numeric(values[methods])
        )
      },
      res_list,
      eta_vec
    )
  )
  
  out
}


#============================================================
# 3. Data
#============================================================

dat_cov <- make_plot_data(
  loss_name = "mean_cov_loss",
  methods = c(
    "ETL-cov",
    "ETL-Scov",
    "Target-cov",
    "Target-Scov",
    "Pool-cov"
  )
)

dat_cor <- make_plot_data(
  loss_name = "mean_cor_loss",
  methods = c(
    "ETL-cov",
    "ETL-Scov",
    "Target-cov",
    "Target-Scov",
    "Pool-cov"
  )
)

dat_var <- make_plot_data(
  loss_name = "mean_var_loss",
  methods = c(
    "ETL-var",
    "Target-var",
    "Pool-var"
  )
)


#============================================================
# 4. Common plotting function
#============================================================

plot_one_panel <- function(
    dat,
    panel_title,
    ylab_text,
    ylim = NULL,
    y_breaks = waiver(),
    y_trans = "identity",
    show_xlab = FALSE) {
  
  dat$panel_title <- panel_title
  
  p <- ggplot(
    dat,
    aes(
      x = eta,
      y = Loss,
      colour = Method,
      group = Method
    )
  ) +
    
    geom_line(
      linewidth = 1
    ) +
    
    geom_point(
      size = 2
    ) +
    
    scale_colour_manual(
      values = color_Values,
      guide = "none"
    ) +
    
    scale_x_continuous(
      breaks = eta_vec
    ) +
    
    scale_y_continuous(
      breaks = y_breaks,
      trans = y_trans,
      labels = function(x) {
        formatC(
          x,
          format = "f",
          digits = 2
        )
      },
      expand = expansion(
        mult = c(0.01, 0.015)
      ),
      minor_breaks = NULL
    ) +
    
    xlab(
      if (show_xlab) {
        TeX("$\\eta$")
      } else {
        NULL
      }
    ) +
    
    ylab(
      TeX(ylab_text)
    ) +
    
    theme_bw() +
    
    facet_grid(. ~ panel_title) +
    
    theme(
      legend.position = "none",
      
      axis.title = element_text(
        size = 12
      ),
      
      axis.text = element_text(
        size = 11
      ),
      
      panel.grid.major = element_line(
        linewidth = 0.4
      ),
      
      panel.grid.minor = element_line(
        linewidth = 0.2
      ),
      
      plot.margin = grid::unit(
        c(0.2, 0.2, 0.2, 0.2),
        "cm"
      )
    )
  
  if (!is.null(ylim)) {
    p <- p +
      coord_cartesian(
        ylim = ylim
      )
  }
  
  p
}


#============================================================
# 5. Three panels
#============================================================

p1 <- plot_one_panel(
  dat = dat_cov,
  panel_title = "Covariance",
  ylab_text =
    "$\\|\\hat{\\Sigma}-\\Sigma^{(0)}\\|_F^2/p$",
  ylim = c(0.05, 1.15),
  y_breaks = cov_breaks,
  show_xlab = FALSE
)


p2 <- plot_one_panel(
  dat = dat_cor,
  panel_title = "Correlation",
  ylab_text =
    "$\\|\\hat{R}-R^{(0)}\\|_F^2/p$",
  ylim = c(0.02, 1.00),
  y_breaks = cor_breaks,
  show_xlab = TRUE
)


p3 <- plot_one_panel(
  dat = dat_var,
  panel_title = "Marginal variances",
  ylab_text =
    "$\\|\\hat{\\Lambda}-\\Lambda^{(0)}\\|_F^2$",
  ylim = c(1.45, 12),
  y_breaks = var_breaks,
  
  y_trans = compress_interval_trans(
    from = compress_from,
    to = compress_to,
    compression = compression
  ),
  
  show_xlab = FALSE
)


#============================================================
# 6. Construct a standalone common legend
#============================================================

# This legend is drawn explicitly.
# It does not depend on cowplot::get_legend().

legend_dat <- data.frame(
  Method = breaks.new,
  x = seq_along(breaks.new)
)

legend_dat$Method <- factor(
  legend_dat$Method,
  levels = breaks.new
)

legend_plot <- ggplot(
  legend_dat
) +
  
  # Colored line segments
  geom_segment(
    aes(
      x = x - 0.25,
      xend = x - 0.07,
      y = 1,
      yend = 1,
      colour = Method
    ),
    linewidth = 1.2
  ) +
  
  # Colored points
  geom_point(
    aes(
      x = x - 0.16,
      y = 1,
      colour = Method
    ),
    size = 2.5
  ) +
  
  # Method names
  geom_text(
    aes(
      x = x - 0.01,
      y = 1,
      label = Method
    ),
    hjust = 0,
    vjust = 0.5,
    size = 3.6,
    colour = "black"
  ) +
  
  scale_colour_manual(
    values = color_Values,
    guide = "none"
  ) +
  
  coord_cartesian(
    xlim = c(0.5, 5.65),
    ylim = c(0.7, 1.3),
    clip = "off"
  ) +
  
  theme_void() +
  
  theme(
    plot.margin = grid::unit(
      c(0, 0, 0, 0),
      "cm"
    )
  )


#============================================================
# 7. Combine three panels
#============================================================

panels <- cowplot::plot_grid(
  p1,
  p2,
  p3,
  nrow = 1,
  rel_widths = c(1, 1, 1),
  align = "hv",
  axis = "tblr"
)


#============================================================
# 8. Place a compact legend above the panels
#============================================================

# Center the legend and reduce its width
legend_centered <- cowplot::plot_grid(
  NULL,
  legend_plot,
  NULL,
  nrow = 1,
  rel_widths = c(0.20, 0.60, 0.20)
)

# Combine legend and panels
p_all <- cowplot::plot_grid(
  legend_centered,
  panels,
  ncol = 1,
  rel_heights = c(0.12, 1)
)

print(p_all)


#============================================================
# 9. Display
#============================================================

print(p_all)


#============================================================
# 10. Save figure
#============================================================

# ggsave(
#   filename = "ETL_sensitivity_eta.pdf",
#   plot = p_all,
#   width = 13,
#   height = 4.5,
#   units = "in"
# )