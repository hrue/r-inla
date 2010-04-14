
rr=scan("zambia.map")

zambia=vector("list",57)

idx=1
for(i in 1:40224)
if(is.na(rr[i]))
{
    zambia[[idx]]=matrix(rr[(i+2):((i+2)+2*rr[i+1]-1)],rr[i+1],2,byrow=T)
    idx=idx+1
}


zambia.map <- function(data, cutpoints=seq(min(data),max(data),length=256), autoscale=FALSE, legend=TRUE, append=FALSE)
{
    if (autoscale)
    {
        data = (data-min(data))/(max(data)-min(data)+1e-8)
    }
    #cutpoints = c(-1e9,cutpoints, 1e9)
    
    farben <- gray(as.numeric(cut(data,cutpoints,include.lowest=T))/length(cutpoints))


    xmin <- 1:length(zambia)
    xmax <- 1:length(zambia)
    ymin <- 1:length(zambia)
    ymax <- 1:length(zambia)

    for(i in 1:length(zambia))
    {
        xmin[i] <- min(zambia[[i]][,1],na.rm=T)
        xmax[i] <- max(zambia[[i]][,1],na.rm=T)
        ymin[i] <- min(zambia[[i]][,2],na.rm=T)
        ymax[i] <- max(zambia[[i]][,2],na.rm=T)
    }

    breite <- c(min(xmin),max(xmax))
    hoehe <- c(min(ymin),max(ymax))

    if (!append) plot(breite,hoehe,type="n",axes=F, xlab=" ", ylab=" ")


    for(k in length(zambia):1)
    {
        polygon(zambia[[k]][,1],zambia[[k]][,2],col=farben[k])
    }
    
    if (legend)
    {
        for(i in 1:8)
        {
             polygon(c(31.5,31.5,32,32),c(-18.2+0.39*(i-1),-18.2+0.39*i,-18.2+0.39*i,-18.2+0.39*(i-1)),col=gray((9-i)/9))
         }
        for(i in 2:8)
        {
            text(32.5,-18.2+0.39*(i-1),as.character(round(cutpoints[i/8*length(cutpoints)], 2)),cex=.7,col=rgb(0,0,0))
        }
    }
}
    
