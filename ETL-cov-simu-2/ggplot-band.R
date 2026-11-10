library(ggplot2)
library(tidyr)
library(dplyr)
library(ggpubr)
library(ggthemes)
library(scales)
#------------------------------------------------
squash_axis = function(from, to, factor) { 
  # Args:
  #   from: left end of the axis
  #   to: right end of the axis
  #   factor: the compression factor of the range [from, to]
  
  trans <- function(x) {    
    # get indices for the relevant regions
    isq <- x > from & x < to
    ito <- x >= to
    
    # apply transformation
    x[isq] <- from + (x[isq] - from)/factor
    x[ito] <- from + (to - from)/factor + (x[ito] - to)
    
    return(x)
  }
  
  inv <- function(x) {
    # get indices for the relevant regions
    isq <- x > from & x < from + (to - from)/factor
    ito <- x >= from + (to - from)/factor
    
    # apply transformation
    x[isq] <- from + (x[isq] - from) * factor
    x[ito] <- to + (x[ito] - (from + (to - from)/factor))
    
    return(x)
  }
  
  # return the transformation
  return(trans_new("squash_axis", trans, inv))
}

############################################################
## 1. Extract results
############################################################

eta.vec <- c(0.1, 0.3, 0.5, 0.7, 0.9)

res.list <- list(c1, c2, c3, c4, c5)

##-============================================================
# Extract out.F.ave from each list
res.mat <- do.call(
  rbind,
  lapply(res.list, function(x) as.numeric(x$out.cov.F.ave))
)

# Method names
method.names <- names(c1$out.cov.F.ave)

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
  ) %>%
  filter(
    !Method %in% c(
      "Target-cov",
      "Pool-cov",
      "Target-scov"
    )
  )

# plot.df
############################################################
## 3. Line plot
############################################################

p1 = ggplot(
  plot.df,
  aes(
    x = eta,
    y = Error,
    # group = Method,
    color = Method,
    # shape = Method
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
  # coord_trans(y = squash_axis(0.35, 1.2, 10)) + ## 将1.5-3压缩6倍
  # scale_y_continuous(breaks = c(0.1,0.2,0.3,0.9)) + 
  labs(
    x = expression(eta),
    y = "Covariance Estimation error",
    color = "Method",
    shape = "Method"
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_blank()
  )

#===========================================================

# Extract out.F.ave from each list
res.mat <- do.call(
  rbind,
  lapply(res.list, function(x) as.numeric(x$out.corr.F.ave))
)

# Method names
method.names <- names(c1$out.cov.F.ave)

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
  ) %>%
  filter(
    !Method %in% c(
      "Target-cov",
      "Pool-cov",
      "Target-scov"
    )
  )

# plot.df
############################################################
## 3. Line plot
############################################################

p2 = ggplot(
  plot.df,
  aes(
    x = eta,
    y = Error,
    # group = Method,
    color = Method,
    # shape = Method
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
  # coord_trans(y = squash_axis(0.35, 1.2, 10)) + ## 将1.5-3压缩6倍
  # scale_y_continuous(breaks = c(0.1,0.2,0.3,0.9)) + 
  labs(
    x = expression(eta),
    y = "Correlation Estimation error",
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

ggarrange(p1, p2, ncol = 2, nrow = 1,
          common.legend = TRUE,legend = "top") 
