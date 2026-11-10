library(ggplot2)
library(ggpubr)
library(scales)
library(ggthemes)
library(latex2exp)

library(scales)

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

color_Values = c("ETL-cov" = "red",
                 "ETL-Scov" = "orange",
                 "Target-band" = "black",
                 # "Pooled-PCA" = "purple",
                 "ADT" = "blue")
out.a = rbind(c1$out.cov.F.ave,c2$out.cov.F.ave,c3$out.cov.F.ave,
              c4$out.cov.F.ave,c5$out.cov.F.ave)[,c(1,2,3,7)]
colnames(out.a)[2] = c("ETL-Scov")
out.a = cbind(out.a, seq(0.1,0.9,0.2))
colnames(out.a) = c("ETL-cov","ETL-scov","ADT","Target-band","kp")


ggplot()+
  geom_line(data = out.a,aes(x = kp,y = `ETL-cov`,colour = "ETL-cov"),size=1)+
  geom_point(data = out.a,aes(x = kp,y = `ETL-cov`,colour = "ETL-cov"),size=2)+
  geom_line(data = out.a,aes(x = kp,y = `ETL-scov`,colour = "ETL-Scov"),size=1)+
  geom_point(data = out.a,aes(x = kp,y = `ETL-scov`,colour = "ETL-Scov"),size=2)+
  # geom_line(data = out.a,aes(x = kp,y = `Pooled-PCA`,colour = "Pooled-PCA"),size=1)+
  # geom_point(data = out.a,aes(x = kp,y = `Pooled-PCA`,colour = "Pooled-PCA"),size=2)+
  geom_line(data = out.a,aes(x = kp,y = `Target-band`,colour = "Target-band"),size=1)+
  geom_point(data = out.a,aes(x = kp,y = `Target-band`,colour = "Target-band"),size=2)+
  geom_line(data = out.a,aes(x = kp,y = `ADT`,colour = "ADT"),size=1)+
  geom_point(data = out.a,aes(x = kp,y = `ADT`,colour = "ADT"),size=2)+
  scale_colour_manual("",values = color_Values)+
  xlab(TeX("$\\eta$")) +
  ylab(TeX("$\\|\\hat{\\Sigma} - \\Sigma^{(0)}\\|_F^2/p$")) +
  scale_x_continuous(breaks = seq(0.1,0.9,0.2)) + 
  # coord_trans(y = squash_axis(0.6,1.8,10)) +
  # scale_y_continuous(limits = c(0.1,1.8),
  #                    breaks = c(0.2,0.4,0.6,1.5)) +
  theme_bw() + 
  theme(
    legend.position = "top")
  # theme(plot.margin=unit(c(0.2,0.2,0.2,0),'cm'))

##############################################################################
##############################################################################

color_Values = c("ETL-cov" = "red",
                 "ETL-Scov" = "orange",
                 "Target-band" = "black",
                 # "Pooled-PCA" = "purple",
                 "ADT" = "blue")
out.a = rbind(c1$out.loss.var.ave,c2$out.loss.var.ave,c3$out.loss.var.ave,
              c4$out.loss.var.ave,c5$out.loss.var.ave)[,c(1,2,3,7)]
colnames(out.a)[2] = c("ETL-Scov")
out.a = cbind(out.a, seq(0.1,0.9,0.2))
colnames(out.a) = c("ETL-cov","ETL-scov","ADT","Target-band","kp")


p2 = ggplot()+
  geom_line(data = out.a,aes(x = kp,y = `ETL-cov`,colour = "ETL-cov"),size=1)+
  geom_point(data = out.a,aes(x = kp,y = `ETL-cov`,colour = "ETL-cov"),size=2)+
  geom_line(data = out.a,aes(x = kp,y = `ETL-scov`,colour = "ETL-Scov"),size=1)+
  geom_point(data = out.a,aes(x = kp,y = `ETL-scov`,colour = "ETL-Scov"),size=2)+
  # geom_line(data = out.a,aes(x = kp,y = `Pooled-PCA`,colour = "Pooled-PCA"),size=1)+
  # geom_point(data = out.a,aes(x = kp,y = `Pooled-PCA`,colour = "Pooled-PCA"),size=2)+
  geom_line(data = out.a,aes(x = kp,y = `Target-band`,colour = "Target-band"),size=1)+
  geom_point(data = out.a,aes(x = kp,y = `Target-band`,colour = "Target-band"),size=2)+
  geom_line(data = out.a,aes(x = kp,y = `ADT`,colour = "ADT"),size=1)+
  geom_point(data = out.a,aes(x = kp,y = `ADT`,colour = "ADT"),size=2)+
  scale_colour_manual("",values = color_Values)+
  xlab(TeX("$\\eta$")) +
  ylab(TeX("$\\|\\hat{\\Lambda^{(0)}} - \\Lambda^{(0)}\\|_F^2$")) +
  scale_x_continuous(breaks = seq(0.1,0.9,0.2)) +
  coord_trans(y = squash_axis(0.4,2,10)) +
  scale_y_continuous(limits = c(0.32,2),breaks = c(0.35,0.5,1,1.5,2)) +
  theme_bw() + 
  theme(
    legend.position = "top")
# theme(plot.margin=unit(c(0.2,0.2,0.2,0),'cm'))

############################################################################

color_Values = c("ETL-cov" = "red",
                 "ETL-Scov" = "orange",
                 "Target-band" = "black",
                 # "Pooled-PCA" = "purple",
                 "ADT" = "blue")
out.a = rbind(c1$out.corr.F.ave,c2$out.corr.F.ave,c3$out.corr.F.ave,
              c4$out.corr.F.ave,c5$out.corr.F.ave)[,c(1,2,3,7)]
colnames(out.a)[2] = c("ETL-Scov")
out.a = cbind(out.a, seq(0.1,0.9,0.2))
colnames(out.a) = c("ETL-cov","ETL-scov","ADT","Target-band","kp")


p3 = ggplot()+
  geom_line(data = out.a,aes(x = kp,y = `ETL-cov`,colour = "ETL-cov"),size=1)+
  geom_point(data = out.a,aes(x = kp,y = `ETL-cov`,colour = "ETL-cov"),size=2)+
  geom_line(data = out.a,aes(x = kp,y = `ETL-scov`,colour = "ETL-Scov"),size=1)+
  geom_point(data = out.a,aes(x = kp,y = `ETL-scov`,colour = "ETL-Scov"),size=2)+
  # geom_line(data = out.a,aes(x = kp,y = `Pooled-PCA`,colour = "Pooled-PCA"),size=1)+
  # geom_point(data = out.a,aes(x = kp,y = `Pooled-PCA`,colour = "Pooled-PCA"),size=2)+
  geom_line(data = out.a,aes(x = kp,y = `Target-band`,colour = "Target-band"),size=1)+
  geom_point(data = out.a,aes(x = kp,y = `Target-band`,colour = "Target-band"),size=2)+
  geom_line(data = out.a,aes(x = kp,y = `ADT`,colour = "ADT"),size=1)+
  geom_point(data = out.a,aes(x = kp,y = `ADT`,colour = "ADT"),size=2)+
  scale_colour_manual("",values = color_Values)+
  xlab(TeX("$\\eta$")) +
  ylab(TeX("$\\|\\hat{\\R^{(0)}} - \\R^{(0)}\\|_F^2/p$")) +
  scale_x_continuous(breaks = seq(0.1,0.9,0.2)) +
  # coord_trans(y = squash_axis(0.4,2,10)) +
  # scale_y_continuous(limits = c(0.32,2),breaks = c(0.35,0.5,1,1.5,2)) +
  theme_bw() + 
  theme(
    legend.position = "top")
# theme(plot.margin=unit(c(0.2,0.2,0.2,0),'cm'))

ggarrange(p3, p2 ,ncol = 2, nrow = 1,
          common.legend = TRUE,legend = "top")
