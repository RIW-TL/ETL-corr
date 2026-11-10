library(ggplot2)
library(ggthemes)
library(latex2exp)
library(dplyr)
library(tidyr)

#============================================================
# 0. Parameters to modify
#============================================================

#------------------------------------------------------------
# k0 values
#------------------------------------------------------------

k0_vec <- 0:5


#------------------------------------------------------------
# y-axis compression
#
# Pool-PCA is much larger than the other methods.
# Compress the upper part of the y-axis.
#------------------------------------------------------------

compress_from <- 1.0
compress_to   <- 3.0
compression   <- 8


#------------------------------------------------------------
# y-axis breaks and limits
#
# Keep about 4--5 ticks.
# Modify here if necessary after seeing the final results.
#------------------------------------------------------------

y_breaks <- c(
  0.5,
  0.75,
  1.0,
  # 2.0,
  3.0
)

y_limits <- c(
  0.45,
  3.0
)


#------------------------------------------------------------
# Panel labels
#------------------------------------------------------------

d_labels <- c(
  a = "d = 4",
  b = "d = 6",
  c = "d = 8"
)


#============================================================
# 1. Compressed y-axis transformation
#============================================================

compress_interval_trans <- function(
    from,
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
    
    middle <- which(
      x > from &
        x < to
    )
    
    upper <- which(
      x >= to
    )
    
    y[middle] <-
      from +
      (x[middle] - from) / compression
    
    y[upper] <-
      from +
      (to - from) / compression +
      (x[upper] - to)
    
    y
  }
  
  
  inverse <- function(x) {
    
    y <- x
    
    transformed_to <-
      from +
      (to - from) / compression
    
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
#
# Six k0 values:
# k0 = 0, 1, 2, 3, 4, 5
#============================================================

res_list <- list(
  
  a = list(
    a0,
    a1,
    a2,
    a3,
    a4,
    a5
  ),
  
  b = list(
    b0,
    b1,
    b2,
    b3,
    b4,
    b5
  ),
  
  c = list(
    c0,
    c1,
    c2,
    c3,
    c4,
    c5
  )
)


#============================================================
# 3. Methods retained in the figure
#
# Original columns:
# 1st  = ETL-PCA
# 3rd  = Target-PCA
# 5th  = Pool-PCA
# last = GB-PCA
#============================================================

breaks.new <- c(
  "ETL-PCA",
  "Target-PCA",
  "Pool-PCA",
  "GB-PCA"
)


color_Values <- c(
  "ETL-PCA"    = "red",
  "Target-PCA" = "black",
  "Pool-PCA"   = "purple",
  "GB-PCA"     = "blue"
)


#============================================================
# 4. Extract plotting data
#============================================================

make_plot_data <- function(
    loss_name = "mean_cov_loss") {
  
  bind_rows(
    
    lapply(
      names(res_list),
      function(setting) {
        
        res <- res_list[[setting]]
        
        
        #----------------------------------------------------
        # Combine k0 = 0,...,5
        #----------------------------------------------------
        
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
        
        
        #----------------------------------------------------
        # Keep only:
        # 1st column
        # 3rd column
        # 5th column
        # last column
        #----------------------------------------------------
        
        out <- out[
          ,
          c(
            1,
            3,
            5,
            ncol(out)
          ),
          drop = FALSE
        ]
        
        
        #----------------------------------------------------
        # Rename methods
        #----------------------------------------------------
        
        colnames(out) <- c(
          "ETL-PCA",
          "Target-PCA",
          "Pool-PCA",
          "GB-PCA"
        )
        
        
        #----------------------------------------------------
        # Add k0
        #----------------------------------------------------
        
        out$k0 <- k0_vec
        
        
        #----------------------------------------------------
        # Wide -> long
        #----------------------------------------------------
        
        out %>%
          
          pivot_longer(
            cols = -k0,
            names_to = "Method",
            values_to = "Loss"
          ) %>%
          
          mutate(
            
            Setting =
              d_labels[setting],
            
            Method =
              factor(
                Method,
                levels = breaks.new
              )
          )
      }
    )
  )
}


#============================================================
# 5. Plot three panels
#============================================================

plot_three_settings <- function() {
  
  dat <- make_plot_data(
    loss_name = "mean_cov_loss"
  )
  
  
  #----------------------------------------------------------
  # Panel order
  #----------------------------------------------------------
  
  dat$Setting <- factor(
    dat$Setting,
    levels = c(
      "d = 4",
      "d = 6",
      "d = 8"
    )
  )
  
  
  #----------------------------------------------------------
  # Plot
  #----------------------------------------------------------
  
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
    
    
    #--------------------------------------------------------
  # Three panels
  #--------------------------------------------------------
  
  facet_grid(
    . ~ Setting
  ) +
    
    
    #--------------------------------------------------------
  # Legend
  #--------------------------------------------------------
  
  scale_colour_manual(
    "",
    breaks = breaks.new,
    values = color_Values
  ) +
    
    
    #--------------------------------------------------------
  # x-axis
  #--------------------------------------------------------
  
  scale_x_continuous(
    breaks = 0:5,
    limits = c(0, 5)
  ) +
    
    
    #--------------------------------------------------------
  # Compressed y-axis
  #--------------------------------------------------------
  
  scale_y_continuous(
    trans = y_trans,
    
    breaks = y_breaks,
    
    labels = function(x) {
      formatC(
        x,
        format = "f",
        digits = 2
      )
    },
    
    limits = y_limits,
    
    expand = expansion(
      mult = c(
        0.01,
        0.015
      )
    ),
    
    minor_breaks = NULL
  ) +
    
    
    #--------------------------------------------------------
  # Axis labels
  #--------------------------------------------------------
  
  xlab(
    TeX("$k_0$")
  ) +
    
    ylab("Subspace estimation error"
    ) +
    
    
    #--------------------------------------------------------
  # Theme
  #--------------------------------------------------------
  
  theme_bw() +
    
    theme(
      
      # legend
      legend.position = "top",
      legend.direction = "horizontal",
      legend.justification = "center",
      
      legend.text =
        element_text(
          size = 11
        ),
      
      
      # equal spacing between panels
      panel.spacing.x =
        grid::unit(
          0.8,
          "cm"
        ),
      
      
      # facet title background
      strip.background =
        element_rect(
          fill = "grey90",
          colour = "black"
        ),
      
      strip.text =
        element_text(
          size = 12
        ),
      
      
      # axis
      axis.title =
        element_text(
          size = 12
        ),
      
      axis.text =
        element_text(
          size = 10
        ),
      
      
      # margins
      plot.margin =
        grid::unit(
          c(
            0.2,
            0.2,
            0.2,
            0.2
          ),
          "cm"
        )
    )
}


#============================================================
# 6. Final plot
#============================================================

p_PCA <- plot_three_settings()

p_PCA
