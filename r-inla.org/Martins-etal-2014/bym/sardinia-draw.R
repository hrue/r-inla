##source("sardinia.map")
load("sardinia.map.RData")

sardinia.map = function(data,
        cutpoints=seq(min(data),max(data),length=256),
        autoscale=FALSE,
        legend=TRUE,
        append=FALSE)
{
    if (autoscale)
        data = (data-min(data))/(max(data)-min(data)+1e-8)
    
    farben <- gray(as.numeric(cut(data,cutpoints,include.lowest=T))/
                   length(cutpoints))

    xmin <- 1:length(the.sardinia.map)
    xmax <- 1:length(the.sardinia.map)
    ymin <- 1:length(the.sardinia.map)
    ymax <- 1:length(the.sardinia.map)

    for(i in 1:length(the.sardinia.map)) {
        xmin[i] <- min(the.sardinia.map[[i]][,2],na.rm=T)
        xmax[i] <- max(the.sardinia.map[[i]][,2],na.rm=T)
        ymin[i] <- min(the.sardinia.map[[i]][,3],na.rm=T)
        ymax[i] <- max(the.sardinia.map[[i]][,3],na.rm=T)
    }

    breite <- c(min(xmin),max(xmax))
    hoehe <- c(min(ymin),max(ymax))

    if (TRUE) {
        ## correct the aspect-ratio,  this seems quite ok
        print(paste("aspect.ratio", diff(hoehe)/diff(breite)))
        fac = 0.58
        for(k in 1:length(the.sardinia.map))
            the.sardinia.map[[k]][, 2] = breite[1] + fac*(the.sardinia.map[[k]][, 2]-breite[1])
    }

    breite[1] = breite[1] - 0.15 * diff(breite)
    if (!append) plot(breite,hoehe,type="n",axes=F, xlab=" ", ylab=" ")

    for(k in length(the.sardinia.map):1) {
        polygon(the.sardinia.map[[k]][,2],the.sardinia.map[[k]][,3],col=farben[k])
    }
    
    if (legend) {
        nc = 20
        x1 = breite[1] - 0.025*diff(breite)
        x2 = x1 + 0.05*diff(breite)
        y1 = hoehe[1] + 0.0*diff(hoehe)
        dy = diff(hoehe)*0.05
        for(i in 1:nc)
            polygon(c(x1, x1, x2, x2),
                    c(y1+dy*(i-1),y1+dy*i,y1+dy*i,y1+dy*(i-1)),col=gray(i/(nc+1)))
        for(i in 2:nc)
            text(x2 + 0.04*diff(breite), y1  + dy*(i-1),
                 as.character(round(cutpoints[i/nc*length(cutpoints)], 2)),
                 cex=.7,col=rgb(0,0,0))
    }
}
