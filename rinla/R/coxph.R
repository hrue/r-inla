## Export: inla.coxph
## Export: inla.rbind.data.frames

##!\name{inla.coxph}
##!\alias{inla.coxph}
##!\alias{inla.rbind.data.frames}
##!\alias{coxph}
##!\alias{cbind.data.frames}
##!\title{Convert a Cox proportional hazard model into Poisson regression}
##!\description{Tools to convert a Cox proportional hazard model into Poisson regression}
##!\usage{
##!    inla.coxph(formula, data, control.hazard = list(), debug=FALSE)
##!    inla.rbind.data.frames(...)
##!}
##!\arguments{
##!  \item{formula}{The formula for the coxph model where the reponse must be a \code{inla.surv}-object.}
##!  \item{data}{All the data used in the formula,  as a list.}
##!  \item{control.hazard}{Control the model for the baseline-hazard; see \code{?control.hazard}.}
##!  \item{debug}{Print debug-information}
##!  \item{...}{Data.frames to be \code{cbind}-ed,  padding with \code{NA}.}
##!}
##!\value{
##!      \code{inla.coxph} returns a list of new expanded variables to be used in the \code{inla}-call.
##!      Note that element \code{data} and \code{data.list} needs to be merged into
##!      a \code{list} to be passed as the \code{data} argument. See the example for details.
##!      \code{inla.rbind.data.frames} returns the new data.frame.
##!}
##!\author{Havard Rue \email{hrue@r-inla.org}}
##!\examples{
##!## How the cbind.data.frames works:
##!df1 = data.frame(x=1:2, y=2:3, z=3:4)
##!df2 = data.frame(x=3:4, yy=4:5, zz=5:6)
##!inla.rbind.data.frames(df1, df2)
##!
##!## Standard example of how to convert a coxph into a Poisson regression
##!n = 1000
##!x = runif(n)
##!lambda = exp(1+x)
##!y = rexp(n, rate=lambda)
##!event = rep(1,n)
##!data = list(y=y, event=event, x=x)
##!y.surv = inla.surv(y, event)
##!intercept1 = rep(1, n)
##!p = inla.coxph(y.surv ~ -1 + intercept1 + x,
##!               list(y.surv = y.surv,  x=x, intercept1 = intercept1))
##!
##!r = inla(p$formula, 
##!        family = p$family, 
##!        data=c(as.list(p$data), p$data.list), 
##!        E = p$E)
##!summary(r)
##!
##!## How to use this in a joint model
##!intercept2 = rep(1, n)
##!y = 1 + x + rnorm(n, sd=0.1)
##!df = data.frame(intercept2, x, y)
##!
##!## new need to cbind the data.frames, and then add the list-part of
##!## the data
##!df.joint = c(as.list(inla.rbind.data.frames(p$data, df)), p$data.list)
##!df.joint$Y = cbind(df.joint$y..coxph, df.joint$y)
##!
##!## merge the formulas, recall to add '-1' and to use the new joint
##!## reponse 'Y'
##!formula = update(p$formula, Y ~ intercept2 -1 + .)
##!
##!rr = inla(formula,
##!        family = c(p$family, "gaussian"),
##!        data = df.joint,
##!        E = df.joint$E)
##!}

`inla.coxph` = function(formula, data, control.hazard = list(), debug=FALSE)
{
    ## convert a coxph-model into poisson-regression and return a new
    ## data-list with the expand variables and new variables to use in
    ## the poisson regression

    if (is.data.frame(data)) {
        data = as.list(data)
    }

    name.y = inla.formula2character(formula[2])
    tmp = (names(data) %in% name.y)
    y.surv = NULL
    if (any(tmp)) {
        if (sum(tmp) > 1) {
            stop(inla.paste(c("Several entries in 'data' match the name of the response:",
                              "response=", name.y, ", matches=", sum(tmp),".")))
        }
        y.surv = data[[which(tmp)]]
        data[[which(tmp)]] = NULL
    } else {
        try.res = try(inla.eval(inla.paste(c("y.surv = with(data,", name.y,")"))), silent=TRUE)
        if (inherits(try.res, "try-error")) {
            stop(inla.paste(c("The reponse '", name.y, "' is not in 'data' and trying to expand it, failed."),
                            sep=""))
        }
    }
    len.y.surv = max(sapply(y.surv, length))
    data.f = inla.fix.data(data, len.y.surv, revert=FALSE)
    data.l = inla.fix.data(data, len.y.surv, revert=TRUE)
    data.f = try(as.data.frame(data.f),  silent=TRUE)
    if (inherits(data.f, "try-error")) {
        stop("Fail to convert 'data' into a 'data.frame' even after 'inla.fix.data'.")
    }

    if (class(y.surv) != "inla.surv") {
        stop(paste("For survival models, then the reponse has to be of class `inla.surv'; you have `",
                   class(y.surv), "'", sep=""))
    }
    control.hazard = inla.check.control(control.hazard, data.f)
    cont.hazard = inla.set.control.hazard.default()
    cont.hazard[names(control.hazard)] = control.hazard
    cont.hazard$hyper = inla.set.hyper(cont.hazard$model, "hazard", cont.hazard$hyper, 
                                       cont.hazard$initial, cont.hazard$fixed, cont.hazard$prior, cont.hazard$param)

    if (is.null(y.surv$subject)) {
        res = inla.expand.dataframe.1(y.surv, data.f, control.hazard = cont.hazard)
    } else {
        res = inla.expand.dataframe.2(y.surv, data.f, control.hazard = cont.hazard)
    }

    strata.var = NULL
    strata.tmp = NULL
    if (!is.null(cont.hazard$strata.name)) {
        if (is.character(cont.hazard$strata.name) && length(cont.hazard$strata.name)==1) {
            ## strata = "x"
            strata.var = cont.hazard$strata
        } else {
            stop("Argument to `strata.name' must be the name of a variable in the data.frame.")
        }
    }
    if (!is.null(strata.var)) {
        if (!is.element(strata.var, names(data.f))) {
            stop(inla.paste(c("Variable `", strata.var,
                              "' in control.hazard=list(strata=...) needs to be in the data.frame: names(data) = ",
                              names(data.f))))
        }
        if (debug) print("apply inla.strata() on strata.var")
        ## cleaner code...
        strata.tmp = inla.strata(inla.get.var(strata.var, data.f))
        data.f[[strata.var]] = strata.tmp$strata
        data.l$baseline.hazard.strata.coding = strata.tmp$coding
        ## old code:
        ##inla.eval(paste("strata.tmp = inla.strata(data.f$", strata.var, ")", sep=""))
        ##inla.eval(paste("data.f$", strata.var, " = strata.tmp$strata"))
        ##inla.eval(paste("data.l$baseline.hazard.strata.coding = strata.tmp$coding"))
    }
    
    ## if -1 the intercept is not included
    intercept = inla.ifelse(attr(terms(formula), "intercept") == 0, FALSE, TRUE)

    f.hazard = paste(
        "~",
        inla.ifelse(intercept, "1 +", "-1 +"), 
        ". + f(baseline.hazard, model=\"", cont.hazard$model,"\"",
        ", values = baseline.hazard.values", 
        ", hyper = ", enquote(cont.hazard$hyper),
        ", constr = ", cont.hazard$constr,
        ", diagonal = ", inla.ifelse((is.null(cont.hazard$diagonal) && cont.hazard$constr),
                                     inla.set.f.default()$diagonal, cont.hazard$diagonal), 
        ", scale.model = ", inla.ifelse(is.null(cont.hazard$scale.model),
                                        inla.getOption("scale.model.default"),
                                        cont.hazard$scale.model), 
        inla.ifelse(is.null(strata.var), "", paste(", replicate=", strata.var)),
        ")", sep="")[2L]
    
    new.formula = update(update(formula, as.formula(f.hazard)), y..coxph ~ .)
    
    data.list = c(res$data.list, data.l)
    if (!is.null(strata.tmp)) {
        data.list = c(data.list, list(baseline.hazard.strata.coding = strata.tmp$coding))
    }

    return (list(formula = new.formula, 
                 data = res$data,
                 data.list = data.list, 
                 family = "poisson", 
                 E = res$data$E..coxph,
                 expand.df = res$data$expand..coxph, 
                 control.hazard = cont.hazard))
}

`inla.rbind.data.frames` = function(...)
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
