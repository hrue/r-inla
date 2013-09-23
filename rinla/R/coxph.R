## Export: inla.coxph
## Export: inla.cbind.data.frames

##! \name{inla.coxph}
##! \name{inla.cbind.data.frames}
##! \alias{inla.coxph}
##! \alias{inla.cbind.data.frames}
##! \alias{coxph}
##! \alias{cbind.data.frames}
##! 
##! \title{Convert a Cox proportional hazard model into Poisson regression}
##! 
##! \description{Tools to convert a Cox proportional hazard model into Poisson regression}
##! \usage{
##!     inla.coxph(formula, data, control.hazard = list())
##!     inla.cbind.data.frames(...)
##! }
##! 
##! \arguments{
##!   \item{formula}{The formula for the coxph model where the reponse must be a \code{inla.surv}-object.}
##!   \item{data}{All the data used in the formula,  as a list.}
##!   \item{control.hazard}{Control the model for the baseline-hazard; see \code{?control.hazard}.}
##!   \item{...}{Data.frames to be cbind-ed,  padding with NA's.}
##!}
##!\value{
##!      \code{inla.coxph} returns a list of new expanded variables to be used in the \code{inla}-call.
##!      \code{inla.cbind.data.frames} returns a new data.frame.
##!}
##!\author{Havard Rue \email{hrue@math.ntnu.no}}
##! 
##!\examples{
##!## Standard example of how to convert a coxph into a Poisson regression
##!n = 100
##!x = runif(n) 
##!lambda = exp(1+x)
##!y = rexp(n, rate=lambda)
##!event = rep(1,n)
##!data = list(y=y, event=event, x=x)
##!y.surv = inla.surv(y, event)
##!p = inla.coxph(yy.surv ~ x, list(y.surv = y.surv,  x=x))
##!model = inla(p$formula, 
##!        family = p$family, 
##!        data=p$data,
##!        E = p$E,
##!        .internal = p$.internal)
##!summary(model)
##!
##!## Doing it manually makes it possible to do a coxph model
##!## jointly with other models
##!interc2 = rep(1, n)
##!y = 1 + x + rnorm(n, sd=0.1)
##!df = data.frame(interc2, x, y)
##!
##!df.joint = inla.cbind.data.frames(p$data, df)
##!## Make the new respose
##!Y = cbind(df.joint$y..coxph, df.joint$y)
##!## Add it to the data.frame which must now become a list
##!df.joint = as.list(df.joint)
##!df.joint$Y = Y
##!## Add the second model into the formula
##!formula = update(p$formula, Y ~ interc2 -1 + .)
##!## And we're done...
##!rr = inla(formula,
##!        family = c("poisson", "gaussian"),
##!        data = df.joint,
##!        E = df.joint$E,
##!        .internal = p$.internal)
##!}

`inla.coxph` = function(formula, data, control.hazard = list())
{
    ## convert a coxph-model into poisson-regression and return a new
    ## data-list with the expand variables and new variables to use in
    ## the poisson regression

    if (is.data.frame(data)) {
        stop("'data' must be a 'list' not 'data.frame'.")
    }

    name.y = inla.formula2character(formula[2])
    tmp = (name.y %in% names(data))
    if (tmp) {
        y.surv = data[[which(tmp)]]
        data[[which(tmp)]] = NULL
    } else {
        stop(inla.paste(c("The reponse '", name.y, "' is not in 'data'."), sep=""))
    }
    data = as.data.frame(data)
    
    if (class(y.surv) != "inla.surv") {
        stop(paste("For survival models, then the reponse has to be of class `inla.surv'; you have `",
                   class(y.surv), "'", sep=""))
    }
    control.hazard = inla.check.control(control.hazard, data)
    cont.hazard = inla.set.control.hazard.default()
    cont.hazard[names(control.hazard)] = control.hazard
    cont.hazard$hyper = inla.set.hyper(cont.hazard$model, "hazard", cont.hazard$hyper, 
            cont.hazard$initial, cont.hazard$fixed, cont.hazard$prior, cont.hazard$param)

    if (is.null(y.surv$subject)) {
        res = inla.expand.dataframe.1(y.surv, data, control.hazard = cont.hazard)
    } else {
        res = inla.expand.dataframe.2(y.surv, data, control.hazard = cont.hazard)
    }

    idx = min(which(names(res$data) %in% ".y.surv"))
    names(res$data)[idx] = "y..coxph"
    idx = min(which(names(res$data) %in% ".E"))
    names(res$data)[idx] = "E..coxph"

    if (!is.null(cont.hazard$cutpoints)) {
        baseline.hazard.values = seq(1, length(cont.hazard$cutpoints)-1)
    } else if (!is.null(cont.hazard$n.intervals)) {
        baseline.hazard.values = seq(1, cont.hazard$n.intervals)
    } else {
        baseline.hazard.values = NULL
    }

    strata.var = NULL
    if (!is.null(cont.hazard$strata.name)) {
        if (is.character(cont.hazard$strata.name) && length(cont.hazard$strata.name)==1) {
            ## strata = "x"
            strata.var = cont.hazard$strata
        } else {
            stop("Argument to `strata.name' must be the name of a variable in the data.frame.")
        }
    }
    
    ##if -1 the intercept is not included
    intercept = inla.ifelse(attr(terms(formula), "intercept") == 0, FALSE, TRUE)

    f.hazard = paste(
            "~",
            inla.ifelse(intercept, "1 +", "-1 +"), 
            ". + f(baseline.hazard, model=\"", cont.hazard$model,"\"",
            inla.ifelse(!is.null(baseline.hazard.values),
                               inla.paste(c(", values = ", inla.2list(baseline.hazard.values))), ""),
            ", hyper = ", enquote(cont.hazard$hyper),
            ", constr = ", cont.hazard$constr,
            ", si = ", inla.ifelse(cont.hazard$si, "TRUE", "FALSE"),
            inla.ifelse(is.null(strata.var), "", paste(", replicate=", strata.var)),
            ")", sep="")[2L]
    
    new.formula = update(update(formula, as.formula(f.hazard)), y..coxph ~ .)
    
    return (list(formula = new.formula, 
                 data = res$data,
                 family = "poisson", 
                 E = res$data$E..coxph, 
                 .internal = list(baseline.hazard.cutpoints = res$cutpoints)))
}

`inla.cbind.data.frames` = function(...)
{
    ## cbind data.frames padding with NA

    cbind.two.data.frames = function(df1, df2)
    {
        df = NULL
        nr1 = nrow(df1)
        nr2 = nrow(df2)
        nams = sort(unique(c(names(df1), names(df2))), decreasing=TRUE)
        
        for(nam in nams) {
            if (nam %in% names(df1)) {
                idx = which(names(df1) %in% nam)
                if (nam %in% names(df)) {
                    idxx = which(names(df) %in% nam)
                    df[1:nr1, idxx] = df1[1:nr1, idx]
                } else {
                    if (!is.null(df)) {
                        df = data.frame(c(df1[1:nr1, idx], rep(NA, nr2)), df)
                    } else {
                        df = data.frame(c(df1[1:nr1, idx], rep(NA, nr2)))
                    }
                    names(df)[1] = nam
                }
            }
            if (nam %in% names(df2)) {
                idx = which(names(df2) %in% nam)
                if (nam %in% names(df)) {
                    idxx = which(names(df) %in% nam)
                    df[nr1 + 1:nr2, idxx] = df2[1:nr2, idx]
                } else {
                    if (!is.null(df)) {
                        df = data.frame(c(rep(NA, nr1), df2[1:nr2, idx]), df)
                    } else {
                        df = data.frame(c(rep(NA, nr1), df2[1:nr2, idx]))
                    }
                    names(df)[1] = nam
                }
            }
        }

        return (df)
    }

    args = list(...)
    stopifnot(length(args) >= 2L)

    for(i in 2L:length(args)) {
        if (i == 2L) {
            df = cbind.two.data.frames(args[[1]], args[[2]])
        } else {
            df = cbind.two.data.frames(df, args[[i]])
        }
    }

    return (df)
}
