# One entry point: open this file in RStudio and click Source.
# Recreates both figures from saved analytical results; no ETL refitting.
local({
  candidates <- unlist(lapply(sys.frames(), function(frame) {
    keys <- intersect(c("ofile","file","filename","fileName","path",
                        "script","scriptPath","filePath"), ls(frame, all.names=TRUE))
    values <- mget(keys, envir=frame, inherits=FALSE)
    unlist(Filter(is.character, values), use.names=FALSE)
  }), use.names=FALSE)
  candidates <- candidates[basename(candidates)=="MAIN.R" & file.exists(candidates)]
  entry <- if (length(candidates)) tail(candidates,1) else NULL
  if (is.null(entry)) {
    arg <- grep("^--file=", commandArgs(), value=TRUE)
    entry <- if (length(arg)) sub("^--file=", "", arg[1]) else "MAIN.R"
  }
  root <- dirname(normalizePath(entry, winslash="/", mustWork=TRUE))
  read_data <- function(name) read.csv(file.path(root, "data", name),
                                     stringsAsFactors=FALSE, check.names=FALSE)
  out <- file.path(root, "results")
  dir.create(out, showWarnings=FALSE)
  pairs <- read_data("correlation_pairs.csv")
  coordinates <- read_data("variance_coordinates.csv")
  scores <- read_data("split_errors.csv")
  p <- nrow(coordinates)
  sources <- c("Ages 20-39", "Ages 40-59")
  source_colors <- c("#0072B2", "#D55E00")
  nearest <- max.col(-as.matrix(pairs[c("difference_20_39", "difference_40_59")]),
                     ties.method="first")
  ratios <- as.matrix(coordinates[c("ratio_20_39", "ratio_40_59")])
  methods <- c("ETL-scov", "ETL-cov", "Target-cov", "Pool-cov")
  metrics <- c("covariance", "correlation", "variance")
  stopifnot(p==54L, nrow(pairs)==choose(p,2), nrow(scores)==400L,
    !anyDuplicated(paste(pairs$j,pairs$l)), all(pairs$j>pairs$l),
    all(pairs$l>=1 & pairs$j<=p), all(is.finite(ratios)), all(ratios>0),
    !anyDuplicated(paste(scores$repetition,scores$method)),
    setequal(scores$method, methods), setequal(scores$repetition,1:100),
    all(is.finite(as.matrix(scores[metrics]))))
  bars <- do.call(rbind, lapply(metrics, function(metric)
    do.call(rbind, lapply(methods, function(method) {
      x <- scores[scores$method==method, metric]
      stopifnot(length(x)==100L)
      data.frame(method=method, metric=metric, mean=mean(x), sd=sd(x),
                 se=sd(x)/sqrt(length(x)), nrepeat=length(x))
    }))))
  bars$low <- pmax(0, bars$mean-1.96*bars$se)
  bars$high <- bars$mean+1.96*bars$se
  stopifnot(identical(scores$variance[scores$method=="ETL-cov"],
                      scores$variance[scores$method=="ETL-scov"]))
  write.csv(bars, file.path(out,"summary_metrics.csv"), row.names=FALSE)

  figure6 <- function() {
    ticks <- unique(c(1L, seq.int(10L,p,10L), p))
    left_ticks <- unique(c(1L, seq.int(10L,p-10L,10L), p))
    layout(matrix(1:2,1L), widths=c(1,1.12))
    par(mar=c(5.1,4.6,3.2,.6), mgp=c(2.9,.7,0), cex=.95, pty="s")
    plot(NA, xlim=c(.5,p+.5), ylim=c(.5,p+.5), asp=1, xaxs="i", yaxs="i",
      xaxt="n", yaxt="n", xlab="Clinical feature (l)", ylab="Clinical feature (j)", bty="n")
    rect(pairs$l-.5,pairs$j-.5,pairs$l+.5,pairs$j+.5,
         col=source_colors[nearest],border=NA)
    axis(1,left_ticks,sprintf("F%02d",left_ticks),cex.axis=.82)
    axis(2,left_ticks,sprintf("F%02d",left_ticks),las=1,cex.axis=.82)
    lines(c(.5,.5,p+.5,.5),c(.5,p+.5,p+.5,.5),col="#444444",lwd=.7)
    title("(a) Closest source for each correlation",line=1.2,cex.main=.97)
    legend_args <- list(x=29*p/54,y=19*p/54,
      legend=sprintf("%s (%.1f%%)",sources,100*tabulate(nearest,2L)/length(nearest)),
      fill=source_colors,bty="n",cex=.83,xpd=FALSE,x.intersp=.7,y.intersp=1.2)
    box <- do.call(graphics::legend,c(legend_args,list(plot=FALSE)))$rect
    legend_args$x <- max(.5,min(legend_args$x,p+.5-box$w-.6))
    legend_args$y <- max(.5+box$h+.6,min(legend_args$y,legend_args$x-.8,p-.1))
    do.call(graphics::legend,legend_args)

    log_ratio <- log2(ratios)
    y_range <- range(c(-1.85,.9,log_ratio))
    y_ticks <- pretty(y_range,n=6L)
    y_ticks <- y_ticks[y_ticks>=y_range[1] & y_ticks<=y_range[2]]
    par(mar=c(5.1,4.8,3.2,1),mgp=c(3.1,.7,0),pty="m")
    plot(NA,xlim=c(.1,p+.9),ylim=y_range,xaxs="i",yaxs="i",xaxt="n",yaxt="n",
      xlab="Clinical feature (j)",ylab=expression(log[2](s[kj]^2/s["0j"]^2)),bty="n")
    abline(h=y_ticks[y_ticks!=0],col="#E8E8E8",lwd=.7)
    abline(h=0,col="#555555",lty=2,lwd=1)
    segments(seq_len(p),log_ratio[,1],seq_len(p),log_ratio[,2],col="#BDBDBD",lwd=.8)
    points(seq_len(p),log_ratio[,1],col=source_colors[1],pch=16,cex=.7)
    points(seq_len(p),log_ratio[,2],col=source_colors[2],pch=17,cex=.75)
    axis(1,ticks,sprintf("F%02d",ticks),cex.axis=.82)
    axis(2,y_ticks,las=1,cex.axis=.83)
    graphics::box(col="#888888",lwd=.6)
    legend("bottomleft",legend=sources,col=source_colors,pch=c(16,17),
      bty="n",cex=.83,horiz=TRUE,inset=c(.01,.025))
    title("(b) Coordinate-wise variance heterogeneity",line=1.2,cex.main=.97)
  }

  figure7 <- function() {
    colors <- setNames(c("#0072B2","#56B4E9","#D55E00","#999999"),methods)
    labels <- setNames(c("ETL-Scov","ETL-cov","Target-cov","Pool-cov"),methods)
    strips <- c("Covariance estimation error","Correlation estimation error",
                "Diagonal variance estimation error")
    par(mfrow=c(1,3),mar=c(4.8,3.7,2.2,.7),mgp=c(3,.65,0),las=1,cex=.96,
        oma=rep(0,4))
    for (i in seq_along(metrics)) {
      order <- if (metrics[i]=="variance") methods[-1L] else methods
      a <- bars[bars$metric==metrics[i],]
      a <- a[match(order,a$method),]
      upper <- max(a$high,a$mean)*1.3
      mid <- barplot(a$mean,plot=FALSE,width=1,space=.45)
      xmax <- max(mid)+.7
      plot(NA,xlim=c(0,xmax),ylim=c(0,upper),xaxs="i",yaxs="i",xaxt="n",yaxt="n",
           xlab="",ylab="",bty="n")
      ticks <- pretty(c(0,upper),n=5); ticks <- ticks[ticks>=0 & ticks<=upper]
      abline(h=ticks,col="#E3E3E3",lwd=.6)
      rect(mid-.5,0,mid+.5,a$mean,col=colors[a$method],border=NA)
      axis(2,ticks,las=1,cex.axis=.85); graphics::box(col="#A0A0A0",lwd=.6)
      nz <- a$high>a$low
      arrows(mid[nz],a$low[nz],mid[nz],a$high[nz],angle=90,code=3,length=.04,lwd=1)
      text(mid,a$high+.035*upper,formatC(a$mean,digits=3,format="f"),cex=.82)
      text(mid,-.04*upper,labels[a$method],srt=28,adj=1,xpd=NA,cex=.8)
      rect(0,upper,xmax,upper*1.075,col="#F2F2F2",border="black",xpd=NA,lwd=.6)
      text(xmax/2,upper*1.036,strips[i],xpd=NA,cex=.84)
    }
  }
  save_pdf <- function(name,width,height,draw) {
    pdf(file.path(out,name),width=width,height=height,useDingbats=FALSE)
    tryCatch(draw(),finally=dev.off())
  }
  save_pdf("Figure6_fullsource.pdf",12,4.8,figure6)
  save_pdf("Figure7_components.pdf",13.2,4.2,figure7)
  cat("Done: two PDFs and summary_metrics.csv in",out,"\n")
})
