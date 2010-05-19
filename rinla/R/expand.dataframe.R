`inla.expand.dataframe.1` = function(response, dataframe, control.hazard = inla.set.control.hazard.default())
{
    n.intervals = control.hazard$n.intervals
    cutpoints = control.hazard$cutpoints

    if(class(response) != "inla.surv")
        stop("Response has to be an object of class `inla.surv'")
    class(response) = NULL
    event = response$event
    nn = length(event)
    
    time = numeric(nn)
    time[response$event==1] = response$time[response$event==1]
    time[response$event==0] = response$lower[response$event==0]

    truncation = response$truncation
    
    ##create cutpoints if not provided
    if(is.null(cutpoints)) 
        cutpoints = seq(0,max(time), len = n.intervals + 1)

    new.data = inla.get.poisson.data.1(time=time, truncation=truncation,
            event=event, cutpoints=cutpoints)
   
    expand = table(new.data$indicator)
    if(!missing(dataframe)) {
        new.dataframe = matrix(0,length(new.data$y),dim(dataframe)[2])
        for(i in 1:dim(dataframe)[2])
            new.dataframe[,i] = rep(dataframe[,i],expand)
        new.dataframe = as.data.frame(new.dataframe)
        ## just give some name that we are able to recognise afterwards
        names(new.dataframe) = paste("fake.dataframe.names",1:dim(new.dataframe)[2])
    }
    else
        new.dataframe=NULL

    res = data.frame(y.surv=new.data$y, E=new.data$E, baseline.hazard=new.data$baseline.hazard,
            dataframe=new.dataframe)
    names(res)[grep("fake.dataframe.names",names(res))] = names(dataframe)

    return (list(data = res, cutpoints = cutpoints))
}

`inla.get.poisson.data.1` = function(time, truncation, event, cutpoints)
{
    data.new = numeric(0)
    nn = length(event)
    start = as.numeric(cut(truncation,cutpoints,include.lowest=FALSE))
    end = as.numeric(cut(time,cutpoints,include.lowest=TRUE))
    ds = diff(cutpoints)

    for(i in 1:length(time)) {
        if(is.na(start[i])) {
            if(end[i]>1)
                dc = cbind(ds[1:(end[i]-1)],rep(0,(end[i]-1)),rep(i,(end[i]-1)),c(1:(end[i]-1)))
            else dc = numeric(0)
            dc = rbind(dc, cbind(time[i]-(cutpoints[end[i]]), event[i], i, end[i]))
            data.new = rbind(data.new,dc)
        }
        else {
            if(start[i]<end[i]) {
                dc = cbind((truncation[i]-cutpoints[start[i]]),0,i,start[i])
                dc = rbind(dc, cbind(ds[(start[i]+1):(end[i]-1)],rep(0,(end[i]-start[i])),
                        rep(i,(end[i]-start[i])),c((start[i]+1):(end[i])-1)))
                dc = rbind(dc,cbind(time[i]-(cutpoints[end[i]]),event[i],i,end[i]))
            } else if(start[i]==end[i]) {
                dc = cbind(time[i]-(cutpoints[end[i]]),event[i],i,end[i])
            }
            else
                stop("Truncation cannot be greater than time")
        }
    }
    data.new = data.frame(E=data.new[,1], y=data.new[,2], indicator=data.new[,3],
            baseline.hazard=data.new[,4])

    return(data.new)
}

`inla.expand.dataframe.2` =function(response, dataframe, control.hazard = inla.set.control.hazard.default())
{
    n.intervals = control.hazard$n.intervals
    cutpoints = control.hazard$cutpoints

    if(class(response) != "inla.surv")
        stop("Response has to be an object of class `inla.surv'")
    class(response) = NULL
    event = response$event
    nn = length(event)
    
    ## nhpp models do not work for interval censoring
    if(sum(event==3)>0 || sum(event==2)>0)
        stop("family = 'coxph' and subject presented in the observations, cannot have event = 2 or 3.")
    
    ## time, event and subjects related to response
    time = numeric(nn)
    time[response$event==1] = response$time[response$event==1]
    time[response$event==0] = response$time[response$event==0]
    subject = response$subject
    covariate = response$covariate 

    ##create cutpoints if not provided
    if(is.null(cutpoints))
        cutpoints = seq(0, max(time), len = n.intervals + 1)

    ## check if covariate like treatment is provided
    if (missing(covariate))
        event = rep(1,max(subject) * (length(cutpoints)-1))

    new.data = inla.get.poisson.data.2(time=time, subject=subject, event=event, cutpoints=cutpoints)
    res = data.frame(y.surv=new.data$event, E=new.data$E, baseline.hazard=new.data$baseline.haz, 
            covariate=new.data$covariate)
   
   
    return (list(data = res, cutpoints = cutpoints))
}

`inla.get.poisson.data.2` = function( subject,time, event, cutpoints)
{
    data.new = numeric(0)
    nn = max(subject)
    ds = diff(cutpoints)
    ris = matrix(0,max(subject),(length(cutpoints)-1))
    dataF = cbind(subject,time,event)
    end=0
    length(end)=0
    for(i in 1:nn)
    {
        da = matrix(dataF[dataF[,1]==i,],ncol=3)
                                        # to find the interval for each recurrent time 
        c = cut(da[,2],cutpoints,labels=1:(length(cutpoints)-1))
        ris[i,] = tapply(da[,3],c,sum)
        ris[i,][is.na(ris[i,])]=0 
        end =c(end, as.numeric(cut( max(time[subject==i]) ,cutpoints,include.lowest=TRUE)))
    }
    
                                        # counting number of events in each interval for every subject
    totalevent = 0
    length(totalevent) = 0
    for(i in 1:nn)
        totalevent = c(totalevent,ris[i,])
    
    ## checking for interval lengths 
    E=numeric(0)
    length(E)=0
    for(i in 1: nn)
    {
        if(end[i]==1)
            E = c(E, c(max(time[subject==i])-cutpoints[1]), rep(0, length(cutpoints)-1-end[i]))
        else{
            E = c(E, ds[1:end[i]-1], max(time[subject==i])-cutpoints[end[i]], rep(0, 
                                                      length(cutpoints)-1-end[i]))
        }
    }
    
    index = rep(1:nn,each=length(cutpoints)-1)
    covariate = rep(covariate, each= (length(cutpoints)-1))
    baseline.haz=rep(1:(length(cutpoints)-1), nn)
    dc=cbind(index,E,totalevent,baseline.haz,covariate)
    data.new=data.frame(indicator=dc[,1], E=dc[,2],event=dc[,3],baseline.haz=dc[,4], covariate=dc[,5])   
    data.new=data.new[data.new[,2]!=0,]

    return(data.new)
}
