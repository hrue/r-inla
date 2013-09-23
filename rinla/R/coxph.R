##! Export: inla.coxph

##! \name{inla.coxph}
##! \alias{inla.coxph}
##! \alias{coxph}
##! 
##! \title{Convert a Cox proportional hazard model into Poisson regression}
##! 
##! \description{This function is used to convert manually a Cox proportional hazard model into Poisson regression}
##! \usage{
##!     inla.coxph(formula, data, control.hazard = list())
##! }
##! 
##! \arguments{
##!   \item{formula}{The formula for the coxph model where the reponse must be a \code{inla.surv}-object.}
##!   \item{data}{All the data used in the formula,  as a list.}
##!   \item{control.hazard}{Control the model for the baseline-hazard.}
##!}
##!\value{
##!      A list of new expanded variables to be used in the \code{inla}-call.
##!}
##!\author{Havard Rue \email{hrue@math.ntnu.no}}
##! 
##!\examples{
##!n = 1000
##!x = runif(n)
##!lambda = exp(1+x)
##!y = rexp(n, rate=lambda)
##!event = rep(1,n)
##!data = list(y=y, event=event, x=x)
##!y.surv = inla.surv(y, event)
##!p = inla.coxph(yy.surv ~ x,
##!        list(y.surv = y.surv,  x=x))
##!model = inla(p$formula, 
##!        family = p$family, 
##!        data=p$data,
##!        E = p$E,
##!        .internal = p$.internal)
##!summary(model)
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

    if (is.null(y.surv$subject)) {
        res = inla.expand.dataframe.1(y.surv, data, control.hazard = cont.hazard)
    } else {
        res = inla.expand.dataframe.2(y.surv, data, control.hazard = cont.hazard)
    }

    idx = min(which(names(res$data) %in% ".y.surv"))
    names(res$data)[idx] = "y.poisson"
    idx = min(which(names(res$data) %in% ".E"))
    names(res$data)[idx] = "E.poisson"

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
    
    f.hazard = paste(
            "~ . + f(baseline.hazard, model=\"", cont.hazard$model,"\"",
            inla.ifelse(!is.null(baseline.hazard.values),
                               inla.paste(c(", values = ", inla.2list(baseline.hazard.values))), ""),
            ", hyper = ", enquote(cont.hazard$hyper),
            ", constr = ", cont.hazard$constr,
            ", si = ", inla.ifelse(cont.hazard$si, "TRUE", "FALSE"),
            inla.ifelse(is.null(strata.var), "", paste(", replicate=", strata.var)),
            ")", sep="")[2L]

    new.formula = update(update(formula, as.formula(f.hazard)), y.poisson ~ . ),     
    return (list(
        formula = new.formula, 
        data = res$data,
        family = "poisson", 
        E = res$data$E.poisson, 
        .internal = list(baseline.hazard.cutpoints = res$cutpoints)))
}
