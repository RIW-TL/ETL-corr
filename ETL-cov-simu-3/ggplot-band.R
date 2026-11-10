library(ggplot2)
library(tidyr)
library(dplyr)

############################################################
## 1. Extract results
############################################################

eta.vec <- c(0.1, 0.3, 0.5, 0.7, 0.9)

res.list <- list(a1, a3, a5, a7, a9)

# Extract out.F.ave from each list
res.mat <- do.call(
  rbind,
  lapply(res.list, function(x) as.numeric(x$out.F.ave))
)

# Method names
method.names <- names(a1$out.F.ave)

colnames(res.mat) <- method.names

# Convert to data frame
res.df <- data.frame(
  eta = eta.vec,
  res.mat,
  check.names = FALSE
)

# res.df

############################################################
## 2. Long format
############################################################

plot.df <- res.df %>%
  pivot_longer(
    cols = -eta,
    names_to = "Method",
    values_to = "Error"
  )

# plot.df
############################################################
## 3. Line plot
############################################################

ggplot(
  plot.df,
  aes(
    x = eta,
    y = Error,
    group = Method,
    color = Method,
    shape = Method
  )
) +
  geom_line(
    linewidth = 1
  ) +
  geom_point(
    size = 3
  ) +
  scale_x_continuous(
    breaks = eta.vec
  ) +
  labs(
    x = expression(eta),
    y = "Estimation error",
    color = "Method",
    shape = "Method"
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_blank()
  )
