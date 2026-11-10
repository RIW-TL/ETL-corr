library(ggplot2)
library(ggpubr)
library(ggthemes)
library(latex2exp)
library(patchwork)

# -----------------------------
# 2. Shared compressed y scale
# -----------------------------

# Keep the scientifically important range below 0.60 at full size,
# compress only 0.60--1.02, and then return to a 1:1 scale so that the
# two target baselines (about 1.055 and 1.098) remain visibly separated.
compress_interval_trans <- function(from = 0.60, to = 1.02,
                                    compression = 5) {
  stopifnot(
    is.numeric(from), length(from) == 1L,
    is.numeric(to), length(to) == 1L,
    is.numeric(compression), length(compression) == 1L,
    from < to,
    compression > 1
  )
  
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
  
  scales::trans_new(
    name = "compressed_interval",
    transform = transform,
    inverse = inverse,
    domain = c(-Inf, Inf)
  )
}

y_trans <- compress_interval_trans(
  from = 0.5,
  to = 1,
  compression = 5
)

# Exactly the same breaks, limits, and transformation are applied to all facets.
y_breaks <- c(0.10, 0.30, 0.50, 1.4)

#--------------------------------------------------------
out.a = rbind(aa1$out.cov.F.ave,aa2$out.cov.F.ave,
              aa3$out.cov.F.ave,aa4$out.cov.F.ave,aa5$out.cov.F.ave)
out.a = as.data.frame(cbind(out.a, seq(0,0.8,0.2)))
colnames(out.a) = c("ETL-cov","ETL-Scov","Target-cov","Pool-cov","Target-Scov","kp")

breaks.new = c("ETL-cov", "ETL-Scov", "Pool-cov", "Target-cov", "Target-Scov")
color_Values = c("ETL-cov" = "red",
                 "ETL-Scov" = "orange",
                 "Target-cov" = "black",
                 "Pool-cov" = "purple",
                 "Target-Scov" = "blue")

out.a$title = "d = 0.3"
k1 = ggplot()+
  geom_line(data = out.a,aes(x = kp,y = `ETL-cov`,colour = "ETL-cov"),size=1)+
  geom_point(data = out.a,aes(x = kp,y = `ETL-cov`,colour = "ETL-cov"),size=2)+
  geom_line(data = out.a,aes(x = kp,y = `ETL-Scov`,colour = "ETL-Scov"),size=1)+
  geom_point(data = out.a,aes(x = kp,y = `ETL-Scov`,colour = "ETL-Scov"),size=2)+
  geom_line(data = out.a,aes(x = kp,y = `Pool-cov`,colour = "Pool-cov"),size=1)+
  geom_point(data = out.a,aes(x = kp,y = `Pool-cov`,colour = "Pool-cov"),size=2)+
  geom_line(data = out.a,aes(x = kp,y = `Target-cov`,colour = "Target-cov"),size=1)+
  geom_point(data = out.a,aes(x = kp,y = `Target-cov`,colour = "Target-cov"),size=2)+
  geom_line(data = out.a,aes(x = kp,y = `Target-Scov`,colour = "Target-Scov"),size=1)+
  geom_point(data = out.a,aes(x = kp,y = `Target-Scov`,colour = "Target-Scov"),size=2)+
  scale_colour_manual("",breaks = breaks.new, values = color_Values)+
  xlab(TeX("$\\eta$")) +
  ylab(TeX("$\\|\\hat{R} - R^{(0)}\\|_F^2/p$")) +
  scale_x_continuous(breaks = seq(0,0.8,0.2)) + 
  scale_y_continuous(
    trans = y_trans,
    breaks = y_breaks,
    labels = function(x) formatC(x, format = "f", digits = 2),
    limits = c(0.02, 1.4),
    expand = expansion(mult = c(0.01, 0.015)),
    minor_breaks = NULL
  ) +
  theme_bw() +
  theme(plot.margin=unit(c(0.2,0.2,0.2,0),'cm')) + 
  facet_grid(. ~ title)

#--------------------------------------------------------
out.b = rbind(bb1$out.cov.F.ave,bb2$out.cov.F.ave,
              bb3$out.cov.F.ave,bb4$out.cov.F.ave,bb5$out.cov.F.ave)
out.b = as.data.frame(cbind(out.b, seq(0,0.8,0.2)))
colnames(out.b) = c("ETL-cov","ETL-Scov","Target-cov","Pool-cov","Target-Scov","kp")

breaks.new = c("ETL-cov", "ETL-Scov", "Pool-cov", "Target-cov", "Target-Scov")
color_Values = c("ETL-cov" = "red",
                 "ETL-Scov" = "orange",
                 "Target-cov" = "black",
                 "Pool-cov" = "purple",
                 "Target-Scov" = "blue")

out.b$title = "d = 0.5"
k2 = ggplot()+
  geom_line(data = out.b,aes(x = kp,y = `ETL-cov`,colour = "ETL-cov"),size=1)+
  geom_point(data = out.b,aes(x = kp,y = `ETL-cov`,colour = "ETL-cov"),size=2)+
  geom_line(data = out.b,aes(x = kp,y = `ETL-Scov`,colour = "ETL-Scov"),size=1)+
  geom_point(data = out.b,aes(x = kp,y = `ETL-Scov`,colour = "ETL-Scov"),size=2)+
  geom_line(data = out.b,aes(x = kp,y = `Pool-cov`,colour = "Pool-cov"),size=1)+
  geom_point(data = out.b,aes(x = kp,y = `Pool-cov`,colour = "Pool-cov"),size=2)+
  geom_line(data = out.b,aes(x = kp,y = `Target-cov`,colour = "Target-cov"),size=1)+
  geom_point(data = out.b,aes(x = kp,y = `Target-cov`,colour = "Target-cov"),size=2)+
  geom_line(data = out.b,aes(x = kp,y = `Target-Scov`,colour = "Target-Scov"),size=1)+
  geom_point(data = out.b,aes(x = kp,y = `Target-Scov`,colour = "Target-Scov"),size=2)+
  scale_colour_manual("",breaks = breaks.new, values = color_Values)+
  xlab(TeX("$\\eta$")) +
  ylab("") +
  scale_x_continuous(breaks = seq(0,0.8,0.2)) + 
  scale_y_continuous(
    trans = y_trans,
    breaks = y_breaks,
    labels = function(x) formatC(x, format = "f", digits = 2),
    limits = c(0.02, 1.4),
    expand = expansion(mult = c(0.01, 0.015)),
    minor_breaks = NULL
  ) +
  theme_bw() +
  theme(plot.margin=unit(c(0.2,0.2,0.2,0),'cm')) + 
  facet_grid(. ~ title)

#--------------------------------------------------------
out.c = rbind(cc1$out.cov.F.ave,cc2$out.cov.F.ave,
              cc3$out.cov.F.ave,cc4$out.cov.F.ave,cc5$out.cov.F.ave)
out.c = as.data.frame(cbind(out.c, seq(0,0.8,0.2)))
colnames(out.c) = c("ETL-cov","ETL-Scov","Target-cov","Pool-cov","Target-Scov","kp")

breaks.new = c("ETL-cov", "ETL-Scov", "Pool-cov", "Target-cov", "Target-Scov")
color_Values = c("ETL-cov" = "red",
                 "ETL-Scov" = "orange",
                 "Target-cov" = "black",
                 "Pool-cov" = "purple",
                 "Target-Scov" = "blue")

out.c$title = "d = 0.7"
k3 = ggplot()+
  geom_line(data = out.c,aes(x = kp,y = `ETL-cov`,colour = "ETL-cov"),size=1)+
  geom_point(data = out.c,aes(x = kp,y = `ETL-cov`,colour = "ETL-cov"),size=2)+
  geom_line(data = out.c,aes(x = kp,y = `ETL-Scov`,colour = "ETL-Scov"),size=1)+
  geom_point(data = out.c,aes(x = kp,y = `ETL-Scov`,colour = "ETL-Scov"),size=2)+
  geom_line(data = out.c,aes(x = kp,y = `Pool-cov`,colour = "Pool-cov"),size=1)+
  geom_point(data = out.c,aes(x = kp,y = `Pool-cov`,colour = "Pool-cov"),size=2)+
  geom_line(data = out.c,aes(x = kp,y = `Target-cov`,colour = "Target-cov"),size=1)+
  geom_point(data = out.c,aes(x = kp,y = `Target-cov`,colour = "Target-cov"),size=2)+
  geom_line(data = out.c,aes(x = kp,y = `Target-Scov`,colour = "Target-Scov"),size=1)+
  geom_point(data = out.c,aes(x = kp,y = `Target-Scov`,colour = "Target-Scov"),size=2)+
  scale_colour_manual("",breaks = breaks.new, values = color_Values)+
  xlab(TeX("$\\eta$")) +
  ylab("") +
  scale_x_continuous(breaks = seq(0,0.8,0.2)) + 
  scale_y_continuous(
    trans = y_trans,
    breaks = y_breaks,
    labels = function(x) formatC(x, format = "f", digits = 2),
    limits = c(0.02, 1.4),
    expand = expansion(mult = c(0.01, 0.015)),
    minor_breaks = NULL
  ) +
  theme_bw() +
  theme(plot.margin=unit(c(0.2,0.2,0.2,0),'cm')) + 
  facet_grid(. ~ title)

ggarrange(k1, k2, k3 ,ncol = 3, nrow = 1,
          common.legend = TRUE,legend = "top")
