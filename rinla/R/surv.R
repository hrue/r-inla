`inla.surv` = function(time, event, time2, truncation, subject)
{
    ret = NULL
    if(missing(subject)) {
        ret = inla.surv.1(time, event, time2, truncation)
    } else {
        if (!missing(time2))
            stop("Argument 'time2' is not allowed when 'subject' is used")
        if (!missing(truncation))
            stop("Argument 'truncation' is not allowed when 'subject' is used")
        ret = inla.surv.2(time, event, subject)
    }
    return (ret)
}

`plot.inla.surv` = function(object,...)
{
    if(!is.null(object$lower)) {
        plot.inla.surv.1(object,...)
    } else {
        plot.inla.surv.2(object,...)
    }
}

`print.inla.surv` = function(object,...)
{
    if(!is.null(object$lower)) {
        print.inla.surv.1(object,...)
    } else {
        print.inla.surv.2(object,...)
    }
}

`as.inla.surv` = function(object,...)
{
    if(is.null(object$subject)) {
        as.inla.surv.1(object,...)
    } else {
        as.inla.surv.2(object,...)
    }
}

`inla.surv.1` = function(time, event, time2, truncation)
{
    ## check that time is present
    if (missing(time))
        stop ("Must have a 'time' argument")
    if (!is.numeric(time))
        stop ("'Time' variable is not numeric")
    nn = length(time)
    if (any(is.na(time))) {
        ## fixup the rest of the argument if some of 'time' is NA; in
        ## this case just fill in some dummy values so the rest of the
        ## routine still goes through
        idx = is.na(time)
        if (!missing(event))
            event[idx] = 1
        if (!missing(time2))
            time2[idx] = 0
        if (!missing(truncation))
            truncation[idx] = 0
    }

    ## check that no time varible is negative
    if (any(time[!is.na(time)] < 0))
        stop("Negative times are not allowed")
    if (!missing(time2)) {
        if (any(time2[!is.na(time2)] < 0))
            stop("Negative times2 are not allowed")
    }
    if (!missing(truncation)) {
        if (any(is.na(truncation)))
            stop("Non valid values for 'truncation")
        if (any(truncation<0))
            stop("Negative truncation times are not allowed")
    }
    ## if event is totally missing assume that all data are failures
    if (missing(event)) {
        event = rep(1, nn)
        warning("'event' is missing: assuming all are observed failures")
    } else if (any(is.na(event))) {
        ## if some of the element in event are missing assume that
        ## they're observed failures (give a warning)
        event[is.na(event)] = 1
        warning("Some elements in `event' are NA: set them to observed failures.")
    }
      
    ## check that event is 0, 1, 2, 3
    if (!all(is.element(event, 0:3)))
        stop("Invalid value for event")
  
    ## check that for event=3 time2 is present
    interval = (event==3)
    if (sum(interval)>0 && missing(time2))
        stop("'time2' has to be present for interval censored data")
    ## and warn that time2 is ignored for even!=3
    if (sum(interval)==0 && !missing(time2))
        warning("'time2' is ignored for data that are not interval censored")
    if (missing(time2))
        time2 = rep(0, nn)

    ## if truncation is missing set it to 0
    if (missing(truncation))
        truncation = rep(0, nn)
    ## check that it has the correct length
    if (length(truncation) != nn)
        stop("'truncation' is of the wrong dimension")
  
    surv.time = numeric(nn)
    surv.upper = numeric(nn)
    surv.lower = numeric(nn)

    observed = (event==1)
    right = (event==0)
    left = (event==2)
    interval = (event==3)

    surv.time[observed] = time[observed]
    surv.lower[right] = time[right]
    surv.upper[left] = time[left]
    surv.lower[interval] = time[interval]
    surv.upper[interval] = time2[interval]

    ss = list(time=surv.time, lower=surv.lower, upper=surv.upper,
            event=event, truncation=truncation)
    class(ss) = "inla.surv"

    return (ss)
}

`plot.inla.surv.1` = function(object, legend=TRUE,...)
{
    time = object$time
    upper = object$upper
    lower = object$lower
    event = object$event
    truncation = object$truncation
    
    nn = length(time)

    xmax = max(time, upper, lower, event)
    xax = c(0, xmax+xmax/8)
    yax = c(0, nn)
    plot(xax, yax, type="n", xlab="time", ylab="", axes=FALSE)
    axis(1)
    axis(2, labels=FALSE, tick=FALSE)
    for(i in 1:nn) {
        if (event[i]==1) {
            lines(c(truncation[i], time[i]), c( i, i), type="l")
            text(time[i], i,"*")
        } else if (event[i]==0) {
            lines(c(truncation[i], lower[i]), c( i, i), type="l")
            text(lower[i], i,">")
        } else if (event[i]==2) {
            lines(c(truncation[i], upper[i]), c( i, i), type="l")
            text(upper[i], i,"<")
        } else if (event[i]==3) {
            lines(c(truncation[i], upper[i]), c( i, i), type="l")
            text(upper[i], i,"|")
            text(lower[i], i,"|")
        }
    }

    if (legend)
    {
        leg.text = "failure"
        leg.symb = "*"
        if (any(event==0)) {
            leg.text = c(leg.text,"right cens")
            leg.symb = paste(leg.symb,">", sep="")
        }
        if (any(event==2)) {
            leg.text = c(leg.text,"left cens")
            leg.symb = paste(leg.symb,"<", sep="")
        }
        if (any(event==3)) {
            leg.text = c(leg.text,"interval cens")
            leg.symb = paste(leg.symb,"|", sep="")
        }
        legend("topright", leg.text, pch=leg.symb, inset=0.02)
    }
}

`print.inla.surv.1` = function(object, quote=FALSE, ...)
{
    invisible(print(as.character.inla.surv.1(object), quote=quote, ...))
}

`as.character.inla.surv.1` = function(object, ...)
{
    class(object) = NULL
    nn = length(object$event)

    interval = (object$event==3)
    out = character(nn)
    out[interval] = paste("[", object$lower[interval],",", object$upper[interval],"]", sep="")
    right = (object$event==0)
    out[right] = paste(object$lower[right],"+", sep="")
    left = (object$event==2)
    out[left] = paste(object$upper[left],"-", sep="")
    failure = (object$event==1)
    out[failure] = paste(object$time[failure], sep="")

    return (out)
}

`is.inla.surv` = function(object) inherits(object, "inla.surv")

`as.inla.surv.1` = function(object)
{
    if (is.list(object)) {
        for (nm in names(object))
            if (!is.element(nm, c("event", "time", "lower", "upper", "truncation")))
                stop(paste("Wrong name:", nm))

        class(object) = "inla.surv"
        return(object)
    }
    if (is.data.frame(object))
        return (as.inla.surv.1(as.list(object)))

    stop("Argument must be a list or a data.frame")
}

`inla.surv.2` = function(time, event, subject)
{
    ## check that time is present
    if (missing(time))
        stop ("Must have a 'time' argument")
    if (!is.numeric(time))
        stop ("'Time' variable is not numeric")
    nn = length(time)

    if (any(is.na(time)))
    {
        ## fixup the rest of the argument if some of 'time' is NA; in
        ## this case just fill in some dummy values so the rest of the
        ## routine still goes through
        idx = is.na(time)
        if (!missing(event))
            event[idx] = 0
    }

    ## check that no time varible is negative
    if (any(time[!is.na(time)] < 0))
        stop("Negative times are not allowed")
    
    ## here event is multiple events if event is totally missing
    ## assume that no tumor/event is found
    if (missing(event)) {
        event = rep(0, nn)
        warning("'event' is missing: assuming no events")
    } else if (any(is.na(event))) {
        ## if some of the element in event are missing assume that
        ## they're no event observed (give a warning)
        event[is.na(event)] = 0
        warning("Some elements in `event' are NA: assume no event detected for these cases.")
    }
      
    
    ## check that event is 0, 1
    if (!all(is.element(event, 0:1)))
        stop("Invalid value for event")
  
    if (missing(subject))
        stop("'subject' is missing")
    
    surv.time = numeric(nn)
    detect = (event ==1)
    notdetect= (event== 0)

    surv.time = time
    surv.time[detect] = time[detect]
    surv.time[notdetect] = time[notdetect]
    
    ss = list(time=surv.time, event=event, subject=subject)
    class(ss) = "inla.surv"

    return (ss)
}


## plotting the time according to event detect or not detect
 
`plot.inla.surv.2` = function(object, legend=TRUE,...)
{
    time = object$time
    event = object$event
    subject = object$subject
    
    nn = max(subject)

    xmax = max(time)
    xmin = min(time)
    xax = c(xmin, xmax+xmax/8)
    yax = c(0, nn+0.5)
    plot(xax, yax, type="n", xlab="time", ylab="", axes=FALSE)
    axis(1)
    axis(2, labels=FALSE, tick=FALSE)
     
    for(i in 1 :nn) {
        lines(time[subject==i], rep(i, length(time[subject==i])), type="l")
        points(time[subject==i & event==1],
               rep(i, length(time[subject==i & event==1])),         
               pch=8, col="red")
        points(time[subject==i & event==0],
               rep(i, length(time[subject==i & event==0])), pch=10, 
               col="green")
    }

    if (legend) {
        leg.text = "detect"
        leg.symb = "*"
        leg.col="red"
        if (any(event==0)) {
            leg.text = c(leg.text,"not detect")
            leg.symb = paste(leg.symb,"o", sep="")
            leg.col=c(leg.col,"green")
        }
        legend("topright", leg.text, pch=leg.symb, col=leg.col, inset=0.02)
    }
}

`print.inla.surv.2` = function(object, quote=FALSE, ...)
{
    invisible(print(as.character.inla.surv.2(object), quote=quote, ...))
}

`as.character.inla.surv.2` = function(object, ...)
{
    class(object) = NULL
    nn = length(object$event)
    out = character(nn)
    detect = (object$event==1)
    out[detect] = paste(object$time[detect], sep="")
    notdetect  = (object$event==0)
    out[notdetect] = paste(object$time[notdetect],"+", sep="")

    return (out)
}

`as.inla.surv.2` = function(object)
{
    if (is.list(object)) {
        for (nm in names(object)) {
            if (!is.element(nm, c("event", "time", "subject")))
                stop(paste("Wrong name:", nm))
        }
        class(object) = "inla.surv"

        return(object)
    }
    if (is.data.frame(object))
        return (as.inla.surv.2(as.list(object)))

    stop("Argument must be a list or a data.frame")
}

`inla.strata` = function(object)
{
    ## similar to survival::strata but with levels = 1, 2, ... 

    xf = as.factor(object)
    return (list(strata=as.numeric(xf), coding = levels(xf)))
}
