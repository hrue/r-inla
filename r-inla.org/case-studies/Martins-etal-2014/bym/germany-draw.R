germany.make.graph <- function(filename)
{
    graph <- list(n=1, nbs = list(), nnbs = c(), map=germany)
    xx <- scan(filename)
    graph$n = xx[1]
    jj = 2
    for(i in 1:graph$n)
    {
        idx = xx[jj]+1
        nnbs = xx[jj+1]
        if (nnbs)
        {
            nbs = xx[(jj+2):(jj+1+nnbs)]+1
        }
        else
        {
            nbs = NULL
        }
        graph$nbs[[idx]] = nbs
        graph$nnbs[idx] = nnbs
        jj = jj + nnbs + 2
    }

    graph
}
germany.map <- function(data,
                        cutpoints=seq(min(data),max(data),length=256),
                        autoscale=FALSE,
                        legend=TRUE,
                        append=FALSE)
{
    if (autoscale)
    {
        data = (data-min(data))/(max(data)-min(data)+1e-8)
    }
    #cutpoints = c(-1e9,cutpoints, 1e9)
    
    farben <- gray(as.numeric(cut(data,cutpoints,include.lowest=T))/length(cutpoints))

    xmin <- 1:length(germany)
    xmax <- 1:length(germany)
    ymin <- 1:length(germany)
    ymax <- 1:length(germany)

    for(i in 1:length(germany))
    {
        xmin[i] <- min(germany[[i]][,2],na.rm=T)
        xmax[i] <- max(germany[[i]][,2],na.rm=T)
        ymin[i] <- min(germany[[i]][,3],na.rm=T)
        ymax[i] <- max(germany[[i]][,3],na.rm=T)
    }

    breite <- c(min(xmin),max(xmax))
    hoehe <- c(min(ymin),max(ymax))

    if (TRUE) {
        ## correct the aspect-ratio,  this seems quite ok
        print(paste("aspect.ratio", diff(hoehe)/diff(breite)))
        fac = 0.75
        for(k in 1:length(germany))
            germany[[k]][, 2] = breite[1] + fac*(germany[[k]][, 2]-breite[1])
        ##breite[2] = breite[1] + fac*diff(breite)
    }

    breite[2] = breite[2] + 0.1 * diff(breite)

    if (!append) plot(breite,hoehe,type="n",axes=F, xlab=" ", ylab=" ")

    for(k in length(germany):1)
    {
        polygon(germany[[k]][,2],germany[[k]][,3],col=farben[k])
    }
    
    if (legend)
    {
        nc = 20
        x1 = 5100
        x2 = x1 + 300
        y1 = 500
        dy = 350

        for(i in 1:nc)
        {
            polygon(c(x1, x1,x2, x2),c(y1+dy*(i-1),y1+dy*i,y1+dy*i,y1+dy*(i-1)),col=gray(i/(nc+1)))
        }
        for(i in 2:nc)
        {
            text(x2+dy/2, y1+dy*(i-1), as.character(round(cutpoints[(i/nc)*length(cutpoints)], 2)),
                 cex=.7,col=rgb(0,0,0))
        }
    }
}
germany.graph.point <- function(x, y)
{
    ##
    require("splancs",quietly=TRUE)

    if (F)
    {
        list(x=mean(x[!is.na(x)]), y=mean(y[!is.na(y)]))
    }
    else
    {
        n.seg = sum(is.na(x))+1
        if (n.seg > 1) 
        {
            na.pos = c(which(is.na(x)),length(x))
        }
        else
        {
            na.pos = length(x)
        }

        pos = list(x=0,y=0,area=0)
        jj = 1
        weig = 0
        for(j in 1:n.seg)
        {
            xx = x[jj:(na.pos[j]-1)]
            yy = y[jj:(na.pos[j]-1)]
            jj = na.pos[j]+1

            area = as.double(areapl(matrix(data=c(xx,yy),length(xx),2)))
            pos$x = pos$x + area*mean(xx)
            pos$y = pos$y + area*mean(yy)
            weig = weig + area
        }
        pos$x = pos$x/weig
        pos$y = pos$y/weig
        pos$area = weig
        pos
    }
}
germany.map.graph <- function(data, cutpoints=seq(0,1,length=256), autoscale=T,
                          legend=F, append=F, pch=20, cex=1.2, ...)
{
    if (autoscale)
    {
        data = (data-min(data))/(max(data)-min(data)+1e-8)
    }
    cutpoints = c(-1e9,cutpoints, 1e9)
    
    farben <- gray(as.numeric(cut(data,cutpoints,include.lowest=T))/length(cutpoints))


    xmin <- 1:length(germany.graph$map)
    xmax <- 1:length(germany.graph$map)
    ymin <- 1:length(germany.graph$map)
    ymax <- 1:length(germany.graph$map)

    for(i in 1:length(germany.graph$map))
    {
        xmin[i] <- min(germany.graph$map[[i]][,2],na.rm=T)
        xmax[i] <- max(germany.graph$map[[i]][,2],na.rm=T)
        ymin[i] <- min(germany.graph$map[[i]][,3],na.rm=T)
        ymax[i] <- max(germany.graph$map[[i]][,3],na.rm=T)
    }

    breite <- c(min(xmin),max(xmax))
    hoehe <- c(min(ymin),max(ymax))

    if (!append) plot(breite,hoehe,type="n",axes=F, xlab=" ", ylab=" ")


    for(k in length(germany.graph$map):1)
    {
        pos = germany.graph.point(germany.graph$map[[k]][,2], germany.graph$map[[k]][,3])
        points(pos, col=farben[k], pch=pch, cex=cex, ...)
    }
    
    for(k in 1:germany.graph$n)
    {
        pos = germany.graph.point(germany.graph$map[[k]][,2],germany.graph$map[[k]][,3])

        if (length(germany.graph$nbs[[k]]) > 0)
        {
            for(kk in germany.graph$nbs[[k]])
            {
                poss = germany.graph.point(germany.graph$map[[kk]][,2], germany.graph$map[[kk]][,3])
                lines(x=c(pos$x,poss$x), y=c(pos$y, poss$y), ...)
            }
        }
    }
    
    if (legend)
    {
        print("AAAAAAAAAAAAAAAAAAAA")
        nc = 20
        for(i in 1:20)
        {
            polygon(c(5800,5800,6100,6100),c(500+350*(i-1),500+350*i,500+350*i,500+350*(i-1)),col=gray(i/(nc+1)))
        }
        for(i in 2:nc)
        {
            text(6375,500+350*(i-1),as.character(round(cutpoints[i/8*length(cutpoints)], 2)),cex=.7,col=rgb(0,0,0))
        }
    }
}
germany.map.add.text <- function(data)
{
    for(k in length(germany.graph$map):1)
    {
        pos = germany.graph.point(germany.graph$map[[k]][,2], germany.graph$map[[k]][,3])
        text(pos, as.character(data[k]))
    }
}

## load the map
##source("~/p/gmrf/draw-maps/germany.txt")
load("germany.map.RData")
## load and build the graph
germany.graph <- germany.make.graph("germany.graph")

