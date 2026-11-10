library(ggplot2)
library(ggthemes)
library(latex2exp)
library(dplyr)
library(tidyr)

#============================================================
# 0. Parameters to modify
#============================================================

#------------------------------------------------------------
# y-axis compression
#------------------------------------------------------------

compress_from <- 0.5
compress_to   <- 1.2
compression   <- 5

#------------------------------------------------------------
# y-axis breaks and limits
#------------------------------------------------------------

y_breaks <- c(0.1, 0.3, 0.5, 1)

y_limits <- c(0, 1.2)

#------------------------------------------------------------
# Panel labels
#------------------------------------------------------------

d_labels <- c(
  a = "d = 0.2",
  b = "d = 0.4",
  c = "d = 0.6"
)

#============================================================
# 1. Compressed y-axis transformation
#============================================================

compress_interval_trans <- function(from,
                                    to,
                                    compression) {
  
  stopifnot(
    is.numeric(from),
    length(from) == 1L,
    is.numeric(to),
    length(to) == 1L,
    is.numeric(compression),
    length(compression) == 1L,
    from < to,
    compression > 1
  )
  
  transform <- function(x) {
    
    y <- x
    
    middle <- which(x > from & x < to)
    upper  <- which(x >= to)
    
    y[middle] <-
      from + (x[middle] - from) / compression
    
    y[upper] <-
      from +
      (to - from) / compression +
      (x[upper] - to)
    
    y
  }
  
  inverse <- function(x) {
    
    y <- x
    
    transformed_to <-
      from + (to - from) / compression
    
    middle <- which(
      x > from &
        x < transformed_to
    )
    
    upper <- which(
      x >= transformed_to
    )
    
    y[middle] <-
      from +
      (x[middle] - from) * compression
    
    y[upper] <-
      to +
      (x[upper] - transformed_to)
    
    y
  }
  
  scales::trans_new(
    name = "compressed_interval",
    transform = transform,
    inverse = inverse,
    domain = c(-Inf, Inf)
  )
}

y_trans <- compress_interval_trans(
  from = compress_from,
  to = compress_to,
  compression = compression
)

#============================================================
# 2. Simulation results
#============================================================

res_list <- list(
  a = list(a0, a1, a2, a3, a4, a5),
  b = list(b0, b1, b2, b3, b4, b5),
  c = list(c0, c1, c2, c3, c4, c5)
)

#============================================================
# 3. Method appearance
#============================================================

breaks.new <- c(
  "ETL-cov",
  "ETL-Scov",
  "Pool-cov",
  "Target-cov",
  "Target-Scov"
)

color_Values <- c(
  "ETL-cov"     = "red",
  "ETL-Scov"    = "orange",
  "Target-cov"  = "black",
  "Pool-cov"    = "purple",
  "Target-Scov" = "blue"
)

#============================================================
# 4. Extract all three settings
#============================================================

make_plot_data <- function(loss_name) {
  
  bind_rows(
    
    lapply(names(res_list), function(setting) {
      
      res <- res_list[[setting]]
      
      out <- do.call(
        rbind,
        lapply(
          res,
          function(x) x[[loss_name]]
        )
      )
      
      out <- as.data.frame(out)
      out$k0 <- 0:5
      
      out %>%
        pivot_longer(
          cols = -k0,
          names_to = "Method",
          values_to = "Loss"
        ) %>%
        mutate(
          Setting = d_labels[setting]
        )
    })
  )
}

#============================================================
# 5. Plot all three panels in ONE ggplot
#============================================================

plot_three_settings <- function(loss_name,
                                ylab_text) {
  
  dat <- make_plot_data(
    loss_name = loss_name
  )
  
  dat$Setting <- factor(
    dat$Setting,
    levels = c(
      "d = 0.2",
      "d = 0.4",
      "d = 0.6"
    )
  )
  
  ggplot(
    dat,
    aes(
      x = k0,
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
    
    facet_grid(
      . ~ Setting
    ) +
    
    scale_colour_manual(
      "",
      breaks = breaks.new,
      values = color_Values
    ) +
    
    scale_x_continuous(
      breaks = 0:5
    ) +
    
    scale_y_continuous(
      trans = y_trans,
      breaks = y_breaks,
      labels = function(x)
        formatC(
          x,
          format = "f",
          digits = 2
        ),
      limits = y_limits,
      expand = expansion(
        mult = c(0.01, 0.015)
      ),
      minor_breaks = NULL
    ) +
    
    xlab(
      TeX("$k_0$")
    ) +
    
    ylab(
      TeX(ylab_text)
    ) +
    
    theme_bw() +
    
    theme(
      
      # legend
      legend.position = "top",
      legend.direction = "horizontal",
      legend.justification = "center",
      
      # equal spacing between the three panels
      panel.spacing.x =
        grid::unit(0.8, "cm"),
      
      # # facet titles
      # strip.background =
      #   element_rect(
      #     fill = "white"
      #   ),
      
      strip.text =
        element_text(
          size = 12
        ),
      
      # margins
      plot.margin =
        grid::unit(
          c(0.2, 0.2, 0.2, 0.2),
          "cm"
        )
    )
}

#============================================================
# 6. Covariance loss
#============================================================

p_cov <- plot_three_settings(
  loss_name = "mean_cor_loss",
  ylab_text =
    "$\\|\\hat{R} - R^{(0)}\\|_F^2/p$"
)

p_cov