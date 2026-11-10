library(ggplot2)
library(ggthemes)
library(latex2exp)
library(dplyr)
library(tidyr)

#============================================================
# 0. Parameters
#============================================================

# x-axis values
k0_vec <- c(2, 4, 6, 8, 10)

# y-axis breaks and limits
# adjust slightly if needed after viewing the final figure
y_breaks <- c(0.15, 0.25, 0.35, 0.45)
y_limits <- c(0.14, 0.48)

# panel labels
d_labels <- c(
  a = "d = 0.2",
  b = "d = 0.3",
  c = "d = 0.4"
)

#============================================================
# 1. Simulation results
#    Five k0 values: 2, 4, 6, 8, 10
#    Objects required in the environment:
#    a1-a5, b1-b5, c1-c5
#============================================================

res_list <- list(
  a = list(a1, a2, a3, a4, a5),
  b = list(b1, b2, b3, b4, b5),
  c = list(c1, c2, c3, c4, c5)
)

#============================================================
# 2. Methods retained in the figure
#============================================================

breaks.new <- c(
  "ETL-Glasso",
  "Target-Glasso",
  "Pool-Glasso",
  "Trans-CLIME"
)

color_Values <- c(
  "ETL-Glasso"    = "red",
  "Target-Glasso" = "black",
  "Pool-Glasso"   = "purple",
  "Trans-CLIME"   = "blue"
)

#============================================================
# 3. Extract plotting data
#============================================================

make_plot_data <- function(loss_name = "mean_precision_loss") {
  bind_rows(
    lapply(
      names(res_list),
      function(setting) {
        res <- res_list[[setting]]

        out <- do.call(
          rbind,
          lapply(
            res,
            function(x) {
              x[[loss_name]]
            }
          )
        )

        out <- as.data.frame(out)

        # keep the 4 displayed methods
        # if the stored object is already a 4-column matrix, this still works
        if (ncol(out) >= 4) {
          out <- out[, 1:4, drop = FALSE]
        }

        colnames(out) <- c(
          "ETL-Glasso",
          "Target-Glasso",
          "Pool-Glasso",
          "Trans-CLIME"
        )

        out$k0 <- k0_vec

        out %>%
          pivot_longer(
            cols = -k0,
            names_to = "Method",
            values_to = "Loss"
          ) %>%
          mutate(
            Setting = d_labels[setting],
            Method = factor(Method, levels = breaks.new)
          )
      }
    )
  )
}

#============================================================
# 4. Plot three panels
#============================================================

plot_three_settings <- function() {
  dat <- make_plot_data(loss_name = "mean_precision_loss")

  dat$Setting <- factor(
    dat$Setting,
    levels = c("d = 0.2", "d = 0.3", "d = 0.4")
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
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    facet_grid(. ~ Setting) +
    scale_colour_manual(
      "",
      breaks = breaks.new,
      values = color_Values
    ) +
    scale_x_continuous(
      breaks = c(2, 4, 6, 8, 10),
      limits = c(2, 10)
    ) +
    scale_y_continuous(
      breaks = y_breaks,
      labels = function(x) formatC(x, format = "f", digits = 2),
      limits = y_limits,
      expand = expansion(mult = c(0.01, 0.015)),
      minor_breaks = NULL
    ) +
    xlab(TeX("$k_0$")) +
    ylab("Precision matrix estimation error") +
    theme_bw() +
    theme(
      legend.position = "top",
      legend.direction = "horizontal",
      legend.justification = "center",
      legend.text = element_text(size = 11),
      panel.spacing.x = grid::unit(0.8, "cm"),
      strip.background = element_rect(fill = "grey90", colour = "black"),
      strip.text = element_text(size = 12),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      plot.margin = grid::unit(c(0.2, 0.2, 0.2, 0.2), "cm")
    )
}

#============================================================
# 5. Final plot
#============================================================

p_precision <- plot_three_settings()
p_precision

