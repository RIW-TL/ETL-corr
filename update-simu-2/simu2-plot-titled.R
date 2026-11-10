library(ggplot2)
library(dplyr)
library(tidyr)
library(latex2exp)
library(patchwork)

#============================================================
# 1. Settings
#============================================================

eta_vec <- c(
  0.1, 0.3, 0.5, 0.7, 0.9
)

res_list <- list(
  s1, s3, s5, s7, s9
)

# Union of all methods
breaks.new <- c(
  "ETL-cov",
  "ETL-Scov",
  "Target-band",
  "ADT"
)

color_Values <- c(
  "ETL-cov"     = "red",
  "ETL-Scov"    = "orange",
  "Target-band" = "black",
  "ADT"         = "purple"
)


#============================================================
# 2. Extract plotting data
#============================================================

make_plot_data <- function(
    loss_name,
    methods) {
  
  out <- do.call(
    rbind,
    lapply(
      res_list,
      function(x) {
        x[[loss_name]][methods]
      }
    )
  )
  
  out <- as.data.frame(out)
  out$eta <- eta_vec
  
  out <- out %>%
    pivot_longer(
      cols = -eta,
      names_to = "Method",
      values_to = "Loss"
    )
  
  out$Method <- factor(
    out$Method,
    levels = breaks.new
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
    "Target-band",
    "ADT"
  )
)

dat_cor <- make_plot_data(
  loss_name = "mean_cor_loss",
  methods = c(
    "ETL-cov",
    "ETL-Scov",
    "Target-band",
    "ADT"
  )
)

dat_var <- make_plot_data(
  loss_name = "mean_var_loss",
  methods = c(
    "ETL-cov",
    "Target-band",
    "ADT"
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
      linewidth = 1,
      show.legend = TRUE
    ) +
    
    geom_point(
      size = 2,
      show.legend = TRUE
    ) +
    
    # Common color scale for all panels
    scale_colour_manual(
      name = NULL,
      values = color_Values,
      breaks = breaks.new,
      limits = breaks.new,
      drop = FALSE,
      guide = guide_legend(
        nrow = 1,
        byrow = TRUE,
        override.aes = list(
          linewidth = 1.2,
          size = 2.5
        )
      )
    ) +
    
    scale_x_continuous(
      breaks = eta_vec
    ) +
    
    scale_y_continuous(
      breaks = y_breaks,
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
      # Legend settings
      legend.position = "top",
      legend.direction = "horizontal",
      legend.justification = "center",
      legend.title = element_blank(),
      
      legend.text = element_text(
        size = 11
      ),
      
      legend.key.width = grid::unit(
        1.3, "cm"
      ),
      
      legend.spacing.x = grid::unit(
        0.25, "cm"
      ),
      
      # Axis settings
      axis.title = element_text(
        size = 12
      ),
      
      axis.text = element_text(
        size = 11
      ),
      
      # Grid settings
      panel.grid.major = element_line(
        linewidth = 0.4
      ),
      
      panel.grid.minor = element_line(
        linewidth = 0.2
      ),
      
      # Margins
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
  
  return(p)
}


#============================================================
# 5. Three panels
#============================================================

p1 <- plot_one_panel(
  dat = dat_cov,
  panel_title = "Covariance",
  ylab_text =
    "$\\|\\hat{\\Sigma} - \\Sigma^{(0)}\\|_F^2/p$",
  ylim = c(
    0.06,
    0.20
  ),
  y_breaks = c(
    0.08,
    0.11,
    0.14,
    0.18
  ),
  show_xlab = FALSE
)


p2 <- plot_one_panel(
  dat = dat_cor,
  panel_title = "Correlation",
  ylab_text =
    "$\\|\\hat{R} - R^{(0)}\\|_F^2/p$",
  ylim = c(
    0.05,
    0.18
  ),
  y_breaks = c(
    0.07,
    0.10,
    0.13,
    0.16
  ),
  show_xlab = TRUE
)


p3 <- plot_one_panel(
  dat = dat_var,
  panel_title = "Marginal variances",
  ylab_text =
    "$\\|\\hat{\\Lambda} - \\Lambda^{(0)}\\|_F^2$",
  show_xlab = FALSE
)


#============================================================
# 6. Combine panels with ONE shared legend
#============================================================

p_all <- (
  p1 + p2 + p3
) +
  plot_layout(
    ncol = 3,
    widths = c(1, 1, 1),
    guides = "collect"
  ) &
  theme(
    legend.position = "top",
    legend.direction = "horizontal",
    legend.justification = "center"
  )


#============================================================
# 7. Display
#============================================================

print(p_all)


#============================================================
# 8. Save figure (optional)
#============================================================

# ggsave(
#   filename = "ETL_sensitivity_eta.pdf",
#   plot = p_all,
#   width = 12,
#   height = 4,
#   units = "in"
# )