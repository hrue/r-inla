## Nothing to Export.

`inla.display.matrix` = function(x, wrap=TRUE, xaxt=FALSE, yaxt=FALSE, col=gray(seq(0, 1, len=256)), ...)
{
    ## display a matrix as an image with correct layout and autoscaling

    x = as.matrix(x)
    n = dim(x)
    y = x
    if (wrap) {
        ii = 1L:n[1L]
        jj = n[1L]:1L
        for(j in 1L:n[2L])
            y[ii, j] = x[jj, j]
    }

    ## use the image.plot-function in package fields; its much better...
    inla.squishplot(c(0, 1), c(0, 1), n[1]/n[2])
    if (inla.require("fields")) {
        fields::image.plot(t(y), col=col, bty="n", xaxt="n", yaxt="n", ...)
    } else {
        warning("Please install package `fields'")
        image(t(y), col=col, bty="n", xaxt="n", yaxt="n", ...)
    }

    box()
    if (xaxt) {
        title(xlab="")
        nn = (n[2]+1)%/%2
        axis(1, at=c(0,(nn-1)/(n[2]-1), 1), labels=as.character(c(1, nn, n[2])))
    }
    if (yaxt) {
        title(ylab="")
        nn = (n[1]+1)%/%2
        axis(2, at=c(0,(nn-1)/(n[1]-1), 1), labels=as.character(c(1, nn, n[1])))
    }
}

`inla.squishplot` = function (xlim, ylim, asp = 1, newplot = TRUE)
{
    ## This function is a copy from package TeachingDemos

    if (length(xlim) < 2)
        stop("xlim must be a vector of length 2")
    if (length(ylim) < 2)
        stop("ylim must be a vector of length 2")
    if (newplot)
        plot.new()
    tmp <- par(c("plt", "pin", "xaxs", "yaxs"))
    if (tmp$xaxs == "i") {
        xlim <- range(xlim)
    } else {
        tmp.r <- diff(range(xlim))
        xlim <- range(xlim) + c(-1, 1) * 0.04 * tmp.r
    }
    if (tmp$yaxs == "i") {
        ylim <- range(ylim)
    } else {
        tmp.r <- diff(range(ylim))
        ylim <- range(ylim) + c(-1, 1) * 0.04 * tmp.r
    }
    tmp2 <- (ylim[2] - ylim[1])/(xlim[2] - xlim[1])
    tmp.y <- tmp$pin[1] * tmp2 * asp
    if (tmp.y < tmp$pin[2]) {
        par(pin = c(tmp$pin[1], tmp.y))
        par(plt = c(tmp$plt[1:2], par("plt")[3:4]))
    } else {
        tmp.x <- tmp$pin[2]/tmp2/asp
        par(pin = c(tmp.x, tmp$pin[2]))
        par(plt = c(par("plt")[1:2], tmp$plt[3:4]))
    }
    return(invisible(tmp["plt"]))
}
