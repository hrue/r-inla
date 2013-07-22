## Export: inla.surv is.inla.surv as.inla.surv
## Export: plot!inla.surv print!inla.surv

##!\name{inla.surv}
##!\alias{inla.surv}
##!\alias{is.inla.surv}
##!\alias{as.inla.surv}
##!\alias{plot.inla.surv}
##!\alias{print.inla.surv}
##!
##!\title{
##!Create a Survival Object for INLA
##!}
##!
##!\description{
##!Create a survival object, to be used as a response variable in a
##!model  formula for the \code{\link{inla}} function for survival models.
##!}
##!\usage{
##!inla.surv(time, event, time2, truncation,subject)
##!\method{plot}{inla.surv}(x, y, ...)
##!\method{print}{inla.surv}(x, ...)
##!is.inla.surv(object)
##!as.inla.surv(object, ...)
##!}
##!
##!\arguments{
##!  \item{time}{For right censored data, this is the follow up time.  For
##!    interval data, this is the starting time for the interval.  }
##!  \item{event}{The status indicator, 1=observed event, 0=right censored
##!    event, 2=left censored event, 3=interval censored event.  Although
##!    unusual, the event indicator can be omitted, in which case all
##!    subjects are assumed to have an event.}
##!  \item{time2}{Ending time for the interval censured data.}
##!  \item{truncation}{Left truncation. If missing it is assumed to be 0.}
##!  \item{subject}{Patient number in multiple event data, not needed otherwise. }
##!  \item{object}{Any \code{R}-object}
##!  \item{x}{Object to plot or print}
##!  \item{y}{Object to plot (not in use)}
##!  \item{...}{Additional argument}
##!  }
##!
##!\value{An object of class \code{inla.surv}.  There are methods for \code{print}, 
##!  \code{plot} for \code{inla.surv} objects.
##!
##!     \code{is.inla.surv} returns \code{TRUE} if \code{object}
##!     inherits from class \code{inla.surv}, otherwise \code{FALSE}.
##!
##!     \code{as.inla.surv} returns an object of class \code{inla.surv}
##!}
##!
##!\author{
##!Sara Martino and Rupali Akerkar
##!}
##!
##!\seealso{
##!\code{\link{inla}}
##!}
##!
##!\examples{
##!  ## First example
##!  trt = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
##!          0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
##!          1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
##!  time = c(17,42,44,48,60,72,74,95,103, 108, 122, 144, 167, 170, 183, 185,
##!           193, 195, 197, 208, 234, 235, 254, 307, 315, 401, 445, 464, 484,  528, 542, 567,
##!           577, 580, 795, 855, 1174, 1214, 1232, 1366, 1455, 1585, 1622, 1626, 1736, 1,63, 
##!           105, 125, 182, 216, 250, 262, 301, 301, 342, 354, 356, 358, 380, 383, 383, 388, 
##!           394, 408, 460, 489, 499, 523, 524, 535, 562, 569, 675, 676, 748, 778, 786, 797,
##!           955, 968, 977, 1245, 1271, 1420, 1460, 1516, 1551, 1690, 1694)
##!  event = c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
##!            1,1,1,1,0,1,0,1,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
##!            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,1,0,1)
##!  y = inla.surv(time, event)
##!
##!  ## Second example
##!  time = c(182,182,63,68,182,152,182,130,134,145,152,182,98,152,182,88,95,105,130,137,167,182,
##!           152,182,81,182,71,84,126,134,152,182)
##!  event = c(1,0,1,1,0,1,0,1,1,1,1,0,1,1,1,1,1,1,1,1,1,0,1,0,1,0,1,1,1,1,1,0)
##!  subject = c(1,2,3,3,3,4,4,5,5,5,5,5,6,6,6,7,7,7,7,7,7,7,8,8,9,9,10,10,10,10,10,10)
##!  y = inla.surv(time, event, subject=subject)
##!}
##!
##!\keyword{Survival models}

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

`plot.inla.surv` = function(x, y,...)
{
    ## argument y is not used
    if(!is.null(x$lower)) {
        inla.plot.inla.surv.1(x,...)
    } else {
        inla.plot.inla.surv.2(x,...)
    }
}

`print.inla.surv` = function(x,...)
{
    if(!is.null(x$lower)) {
        inla.print.inla.surv.1(x,...)
    } else {
        inla.print.inla.surv.2(x,...)
    }
}

`as.inla.surv` = function(object, ...)
{
    ## the '...' are not used.
    if(is.null(object$subject)) {
        inla.as.inla.surv.1(object)
    } else {
        inla.as.inla.surv.2(object)
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

`inla.plot.inla.surv.1` = function(object, legend=TRUE,...)
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

`inla.print.inla.surv.1` = function(object, quote=FALSE, ...)
{
    invisible(print(inla.as.character.inla.surv.1(object), quote=quote, ...))
}

`inla.as.character.inla.surv.1` = function(x, ...)
{
    class(x) = NULL
    nn = length(x$event)

    interval = (x$event==3)
    out = character(nn)
    out[interval] = paste("[", x$lower[interval],",", x$upper[interval],"]", sep="")
    right = (x$event==0)
    out[right] = paste(x$lower[right],"+", sep="")
    left = (x$event==2)
    out[left] = paste(x$upper[left],"-", sep="")
    failure = (x$event==1)
    out[failure] = paste(x$time[failure], sep="")

    return (out)
}

`is.inla.surv` = function(object) inherits(object, "inla.surv")

`inla.as.inla.surv.1` = function(object)
{
    if (is.list(object)) {
        for (nm in names(object))
            if (!is.element(nm, c("event", "time", "lower", "upper", "truncation")))
                stop(paste("Wrong name:", nm))

        class(object) = "inla.surv"
        return(object)
    }
    if (is.data.frame(object))
        return (inla.as.inla.surv.1(as.list(object)))

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
`inla.plot.inla.surv.2` = function(object, legend=TRUE,...)
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

`inla.print.inla.surv.2` = function(object, quote=FALSE, ...)
{
    invisible(print(inla.as.character.inla.surv.2(object), quote=quote, ...))
}

`inla.as.character.inla.surv.2` = function(x, ...)
{
    class(x) = NULL
    nn = length(x$event)
    out = character(nn)
    detect = (x$event==1)
    out[detect] = paste(x$time[detect], sep="")
    notdetect  = (x$event==0)
    out[notdetect] = paste(x$time[notdetect],"+", sep="")

    return (out)
}

`inla.as.inla.surv.2` = function(object)
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
        return (inla.as.inla.surv.2(as.list(object)))

    stop("Argument must be a list or a data.frame")
}

`inla.strata` = function(object)
{
    ## similar to survival::strata but with levels = 1, 2, ... 

    xf = as.factor(object)
    return (list(strata=as.numeric(xf), coding = levels(xf)))
}
