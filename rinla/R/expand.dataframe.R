## Nothing to export from here

`inla.expand.dataframe.1` = function(response, dataframe, control.hazard = inla.set.control.hazard.default(), suffix="")
{
    n.intervals = control.hazard$n.intervals
    cutpoints = control.hazard$cutpoints

    if(class(response) != "inla.surv")
        stop("Response has to be an object of class `inla.surv'")
    class(response) = NULL
    event = response$event
    nn = length(event)
    
    time = numeric(nn)
    time[response$event == 1L] = response$time[response$event == 1L]
    time[response$event == 0L] = response$lower[response$event == 0L]

    truncation = response$truncation
    
    ##create cutpoints if not provided
    if(is.null(cutpoints)) 
        cutpoints = seq(0L, max(time), len = n.intervals +1L) 

    new.data = inla.get.poisson.data.1(time=time, truncation=truncation, event=event, cutpoints=cutpoints)
    expand.df = table(new.data$indicator)

    if(!missing(dataframe) && prod(dim(dataframe)) > 0L) {
        new.dataframe = as.data.frame(matrix(0.0, length(new.data$y), dim(dataframe)[2L]))
        for(i in 1L:dim(dataframe)[2L])
            new.dataframe[, i] = rep(dataframe[, i], expand.df)

        ## alternative code,  not faster...
        ## new.dataframe = apply(dataframe, 2, rep, times = expand.df)  

        new.dataframe = as.data.frame(new.dataframe)
        ## just give some name that we are able to recognise afterwards
        names(new.dataframe) = paste("fake.dataframe.names", 1L:dim(new.dataframe)[2L])
    } else {
        new.dataframe=NULL
    }

    res = data.frame(
        y..coxph = new.data$y,
        E..coxph = new.data$E,
        expand..coxph = rep(1:nrow(dataframe), expand.df), 
        baseline.hazard = cutpoints[new.data$baseline.hazard],
        baseline.hazard.idx = new.data$baseline.hazard,
        baseline.hazard.time = cutpoints[new.data$baseline.hazard],
        baseline.hazard.length = diff(cutpoints)[new.data$baseline.hazard],
        new.dataframe)
    
    names(res)[grep("fake.dataframe.names", names(res))] = names(dataframe)

    return (list(data = res, data.list = list(baseline.hazard.values = cutpoints)))
}

`inla.get.poisson.data.1` = function(time, truncation, event, cutpoints)
{
    data.new = numeric(0L)
    nn = length(event)
    start = as.numeric(cut(truncation, cutpoints, include.lowest=FALSE))
    end = as.numeric(cut(time, cutpoints, include.lowest=TRUE))
    ds = diff(cutpoints)

    for(i in 1L:length(time)) {
        if(is.na(start[i])) {
            if(end[i]>1.0) {
                dc = cbind(ds[1L:(end[i]-1L)], rep(0L,(end[i]-1L)), rep(i,(end[i]-1L)),
                           c(1L:(end[i]-1L)))
            } else {
                dc = numeric(0L)
            }
            dc = rbind(dc, cbind(time[i]-(cutpoints[end[i]]), event[i], i, end[i]))
            data.new = rbind(data.new, dc)
        }
        else {
            if(start[i]<end[i]) {
                dc = cbind((cutpoints[start[i]+1L]-truncation[i]), 0L, i, start[i])
                if(end[i]>(start[i]+1L)) {
                    dc = rbind(dc, cbind(ds[(start[i]+1L):(end[i]-1L)], rep(0L,(end[i]-start[i]-1L)),
                                         rep(i,(end[i]-start[i]-1L)),
                                         c((start[i]+2L):(end[i])-1L)))
                }
                dc = rbind(dc, cbind(time[i]-(cutpoints[end[i]]), event[i], i, end[i]))
                data.new = rbind(data.new, dc)

            } else if(start[i]==end[i]) {
                dc = cbind(time[i]-(cutpoints[end[i]]), event[i], i, end[i])
                data.new = rbind(data.new, dc)
            } else {
                stop("Truncation cannot be greater than time")
            }
        }
    }
    data.new = data.frame(
        E=data.new[, 1L],
        y=data.new[, 2L],
        indicator=data.new[, 3L],
        baseline.hazard=data.new[, 4L])

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
    
    ##nhpp models do not work for interval censoring
    if(sum(event==3L)>0L || sum(event==2L)>0L)
        stop("coxph model does not work for event type 2 and 3.")

    ## time, event and subjects related to response
    time = numeric(nn)
    time[response$event==1L] = response$time[response$event==1L]
    time[response$event==0L] = response$time[response$event==0L]
    subject = response$subject
   
    ## checking for fixed covariates:: nhpp model works for fixed covariates only
    aa1=which(names(dataframe)=="subject")
    aa2=which(names(dataframe)=="time")
    aa3=which(names(dataframe)=="event")
    jj=unique(dataframe$subject)
    subject.first.line = numeric(length(jj))
    for(i in 1L: length(jj))
    {
        rows = which(dataframe[, aa1]==i)
        sem = dataframe[rows,-c(aa1, aa2, aa3), drop=FALSE]
        if(mode(apply(sem, 2L, unique))=="list")
            stop("coxph with subject only works for fixed covariates")
        subject.first.line[i] = rows[1L]
    }
    dataframe.copy = dataframe[subject.first.line, -c(aa1, aa2, aa3), drop=FALSE]
    if(is.null(cutpoints)) 
        cutpoints = seq(0.0, max(time), len = n.intervals +1L) 

    new.data = inla.get.poisson.data.2(time=time, subject=subject, event=event, cutpoints=cutpoints)
   
    ## we want to expand  only covariates 
    aa = table(new.data$indicator)
    ind = unique(new.data$indicator)
    stopifnot(dim(dataframe)[2L] > 3L)
    new.dataframe = as.data.frame(matrix(0.0, length(new.data$y), dim(dataframe)[2L]-3L))
    
    col.data = grep("(subject)|(time)|(event)", names(dataframe))
    if (length(col.data) != 3L)
        stop("data.frame does not contains columns with names `subject', `time' and `event'")

    ##  rewriting the covariates as per new data
    for(i in 1L: dim(dataframe.copy)[2L])
        new.dataframe[, i] = rep(dataframe.copy[, i], aa)
    names(new.dataframe) = names(dataframe)[-col.data]
   
    res = data.frame(
            y..coxph = new.data$y,
            E..coxph = new.data$E,
            baseline.hazard = cutpoints[new.data$baseline.haz], 
            baseline.hazard.idx = new.data$baseline.haz,
            baseline.hazard.time = cutpoints[new.data$baseline.haz],
            baseline.hazard.length = diff(cutpoints)[new.data$baseline.haz],
            subject = new.data$indicator,
            new.dataframe)
    
    return (list(data = res, data.list = list(baseline.hazard.values = cutpoints)))
}

`inla.get.poisson.data.2` = function( subject, time, event, cutpoints)
{
    data.new = numeric(0L)
    nn = max(subject)
    ds = diff(cutpoints)
    ris = matrix(0.0, max(subject),(length(cutpoints)-1L))
    dataF = cbind(subject, time, event)
    end=integer(0L)
    for(i in 1L:nn)
    {
        da = matrix(dataF[dataF[, 1L]==i,], ncol=3L)
        ## to find the interval for each recurrent time
        rec = cut(da[, 2L], cutpoints, labels=1L:(length(cutpoints)-1L))
        ris[i,] = tapply(da[, 3L], rec, sum)
        ris[i,][is.na(ris[i,])]=0 
        end =c(end, as.numeric(cut( max(time[subject==i]) , cutpoints, include.lowest=TRUE)))
    }
    
    ## counting number of events in each interval for every subject
    totalevent = 0
    length(totalevent) = 0
    for(i in 1L:nn)
        totalevent = c(totalevent, ris[i,])
    
    ## checking for interval lengths 
    E=numeric(0)
    length(E)=0
    for(i in 1L: nn)
    {
        if(end[i]==1L)
            E = c(E, c(max(time[subject==i])-cutpoints[1L]), rep(0.0, length(cutpoints)-1L-end[i]))
        else{
            E = c(E, ds[1L:end[i]-1L], max(time[subject==i])-cutpoints[end[i]],
                    rep(0.0, length(cutpoints)-1L-end[i]))
        }
    }
    
    ## combining the number of events, interval lengths, baseline.hazard and subject.
    index = rep(1L:nn, each=length(cutpoints)-1L)
    baseline.haz=rep(1L:(length(cutpoints)-1L), nn)
    dc=cbind(index, E, totalevent, baseline.haz)
    data.new=data.frame(indicator=dc[, 1L], E=dc[, 2L], y=dc[, 3L], baseline.haz=dc[, 4L]) #y= no.of events
    data.new=data.new[data.new[, 2L] != 0L, ]

    return(data.new)
}
