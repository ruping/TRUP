.smoothScatterCalcDensity1 <- function(x, nbin, bandwidth, range.x) {
  
  if (length(nbin) == 1)
    nbin <- c(nbin, nbin)
  if (!is.numeric(nbin) || (length(nbin)!=2))
    stop("'nbin' must be numeric of length 1 or 2")

  if (missing(bandwidth)) {
    bandwidth <- diff(apply(x, 2, quantile, probs=c(0.05, 0.95), na.rm=TRUE)) / 25
  } else {
    if(!is.numeric(bandwidth))
      stop("'bandwidth' must be numeric")
  }
  ## create density map
  if(missing(range.x))
     rv <- bkde2D(x, gridsize=nbin, bandwidth=bandwidth)
  else
     rv <- bkde2D(x, gridsize=nbin, bandwidth=bandwidth, range.x=range.x) 
  rv$bandwidth <- bandwidth
  return(rv)
}

image.plot2 = function (..., add = FALSE, nlevel = 64, legend.shrink = 0.9, 
    legend.width = 1.2, legend.mar = NULL, graphics.reset = FALSE, 
    horizontal = FALSE, bigplot = NULL, smallplot = NULL, legend.only = FALSE, 
    col = tim.colors(nlevel)) 
{
    old.par <- par(no.readonly = TRUE)
    info <- image.plot.info(...)
    if (add) {
        big.plot <- old.par$plt
    }
    if (legend.only) {
        graphics.reset <- TRUE
    }
    if (is.null(legend.mar)) {
        legend.mar <- ifelse(horizontal, 3.1, 5.1)
    }
    temp <- image.plot.plt(add = add, legend.shrink = legend.shrink, 
        legend.width = legend.width, legend.mar = legend.mar, 
        horizontal = horizontal, bigplot = bigplot, smallplot = smallplot)
    smallplot <- temp$smallplot
    bigplot <- temp$bigplot
    if (!legend.only) {
        if (!add) {
            par(plt = bigplot)
        }
	##par(bty = 'n')
        image(..., add = add, col = col)
        box()
	
        big.par <- par(no.readonly = TRUE)
	
    }
    if ((smallplot[2] < smallplot[1]) | (smallplot[4] < smallplot[3])) {
        par(old.par)
        stop("plot region too small to add legend\n")
    }
    ix <- 1
    minz <- info$zlim[1]
    maxz <- info$zlim[2]
    binwidth <- (maxz - minz)/nlevel
    midpoints <- seq(minz + binwidth/2, maxz - binwidth/2, by = binwidth)
    iy <- midpoints
    iz <- matrix(iy, nrow = 1, ncol = length(iy))
    breaks <- list(...)$breaks
    par(new = TRUE, pty = "m", plt = smallplot, err = -1)
    if (!horizontal) {
        if (is.null(breaks)) {
            image(ix, iy, iz, xaxt = "n", yaxt = "n", xlab = "", 
                ylab = "", col = col)
            axis(4, mgp = c(3, 1, 0), las = 2)
        }
        else {
            image(ix, iy, iz, xaxt = "n", yaxt = "n", xlab = "", 
                ylab = "", col = col, breaks = breaks)
            axis(4, at = breaks, labels = format(breaks), mgp = c(3, 
                1, 0), las = 2)
        }
    }
    else {
        if (is.null(breaks)) {
            image(iy, ix, t(iz), xaxt = "n", yaxt = "n", xlab = "", 
                ylab = "", col = col)
            axis(1, mgp = c(3, 1, 0))
        }
        else {
            image(iy, ix, t(iz), xaxt = "n", yaxt = "n", xlab = "", 
                ylab = "", col = col, breaks = breaks)
            axis(1, at = breaks, labels = format(breaks), mgp = c(3, 
                1, 0))
        }
    }
    box()
    mfg.save <- par()$mfg
    if (graphics.reset | add) {
        par(old.par)
        par(mfg = mfg.save, new = FALSE)
        invisible()
    }
    else {
        par(big.par)
        par(plt = big.par$plt, xpd = TRUE)
        par(mfg = mfg.save, new = FALSE)
        invisible()
    }
}


smkey <- function(x, y=NULL, 
                          nbin=128,
                          bandwidth,
                          colramp=colorRampPalette(c("white", brewer.pal(9, "Blues"))),
                          nrpoints=100,
                          transformation=function(x) x^.25,
                          xlab=NULL, ylab=NULL, postPlotHook=box,
                          pch=".", cex=1,
                          xlim, ylim, col="black",
                          xaxs=par("xaxs"), yaxs=par("yaxs"), ...) {
  
  if (!is.numeric(nrpoints) | (nrpoints<0) | (length(nrpoints)!=1) )
    stop("'nrpoints' should be numeric scalar with value >= 0.")

  ## similar as in plot.default
  xlabel <- if (!missing(x)) 
    deparse(substitute(x))
  ylabel <- if (!missing(y)) 
    deparse(substitute(y))
  xy <- xy.coords(x, y, xlabel, ylabel)
  xlab <- if (is.null(xlab)) 
    xy$xlab
  else xlab
  ylab <- if (is.null(ylab)) 
    xy$ylab
  else ylab


  ## eliminate NA
  x <- cbind(xy$x, xy$y)[!(is.na(xy$x)|is.na(xy$y)), ]

  ## xlim and ylim
  if(!missing(xlim)) {
    stopifnot(is.numeric(xlim), length(xlim)==2, !any(is.na(xlim)))
    x <- x[ (x[,1]>=xlim[1]) & (x[,1]<=xlim[2]), ]
  } else {
    xlim <- range(x[,1], na.rm=TRUE)
  }
  if(!missing(ylim)) {
    stopifnot(is.numeric(ylim), length(ylim)==2, !any(is.na(ylim)))
    x <- x[ (x[,2]>=ylim[1]) & (x[,2]<=ylim[2]), ]
  } else {
    ylim <- range(x[,2], na.rm=TRUE)
  }
  
  ## create density map
  map  <- .smoothScatterCalcDensity1(x, nbin, bandwidth)
  xm   <- map$x1
  ym   <- map$x2
  dens <- map$fhat
  dens <- array(transformation(dens), dim=dim(dens))	
  

	
  	  
  ## plot color image
  image.plot2(xm, ym, z=dens, legend.shrink = 1.0,
        xlab = xlab, ylab = ylab, nlevel = 256,...)
  if(!is.null(postPlotHook)) postPlotHook()
  
  ## plot selection of dots
  if (nrpoints!=0){
    ## we assume that map$x1 and map$x2 go linearly from
    ## their first to their last value in nbin steps
    stopifnot(length(xm)==nrow(dens), length(ym)==ncol(dens))
    ixm <- round((x[,1]-xm[1])/(xm[length(xm)]-xm[1])*(length(xm)-1))
    iym <- round((x[,2]-ym[1])/(ym[length(ym)]-ym[1])*(length(ym)-1))
    idens <- dens[1 + iym*length(xm) + ixm]
    nrpoints <- min(nrow(x), ceiling(nrpoints))
    sel <- order(idens, decreasing=FALSE)[1:nrpoints]
    points(x[sel,1:2], pch=pch, cex=cex, col=col)
  }
}
