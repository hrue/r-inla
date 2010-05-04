rr = scan(system.file("demodata/Leuk.map", package="INLA"))
Leuk.newengland=vector("list",24)
idx=1
for(i in 1:51526) {
    if(is.na(rr[i]))
    {
        Leuk.newengland[[idx]]=matrix(rr[(i+2):((i+2)+2*rr[i+1]-1)],rr[i+1],2,byrow=T)
        idx=idx+1
    }
}

Leuk.names.region=scan(system.file("demodata/Leuk.names.regions", package="INLA"))


Leuk.map <- function(data, cutpoints=seq(min(data),max(data),length=256),
                           autoscale=FALSE, legend=TRUE, append=FALSE,
                           map=Leuk.newengland, id=Leuk.names.region)
{
    ## This is not clever done; I know....
    
    if (autoscale)
    {
        data = (data-min(data))/(max(data)-min(data)+1e-8)
    }
    if(length(unique(data))>1) 
        farben <- gray(1-as.numeric(cut(data,cutpoints,include.lowest=T))/length(cutpoints))
    else  farben <- rep(0,24)	

    xmin <- 1:length(map)
    xmax <- 1:length(map)
    ymin <- 1:length(map)
    ymax <- 1:length(map)

    for(i in 1:length(map))
    {
        xmin[i] <- min(map[[i]][,1],na.rm=T)
        xmax[i] <- max(map[[i]][,1],na.rm=T)
        ymin[i] <- min(map[[i]][,2],na.rm=T)
        ymax[i] <- max(map[[i]][,2],na.rm=T)
    }

    breite <- c(min(xmin),max(xmax))
    hoehe <- c(min(ymin),max(ymax))

    if (!append) plot(breite,hoehe,type="n",axes=F, xlab=" ", ylab=" ")


    for(k in length(map):1)
    {
        polygon(map[[k]][,1],map[[k]][,2],col=farben[id[k]])
    }
    
    if (legend)
    {
        for(i in 1:8)
        {
            polygon(c(0.7,0.7,0.74,0.74),c(0.73+0.03*(i-1),0.73+0.03*i,0.73+0.03*i,0.73+0.03*(i-1)),col=gray((9-i)/9))
        }
        for(i in 2:8)
        {
            text(0.79,0.73+0.03*(i-1),as.character(round(cutpoints[i/8*length(cutpoints)], 2)),cex=.7,col=rgb(0,0,0))
        }
    }
}
