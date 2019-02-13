## Export: summary!inla print!summary.inla

##!\name{summary.inla}
##!\alias{summary.inla}
##!\alias{summary.surv.inla}
##!\alias{print.summary.inla}
##!\title{ Summary for a INLA fit }
##!\description{
##!Takes a fitted \code{inla} or \code{surv.inla} object produced by
##!\code{inla} or  \code{surv.inla} and produces
##!a summary from it.
##!
##!}
##!\usage{
##!\method{summary}{inla}(object, digits = 3L, include.lincomb = TRUE, ...)
##!\method{print}{summary.inla}(x, digits = 3L, ...)
##!}
##!%- maybe also 'usage' for other objects documented here.
##!\arguments{
##!  \item{object}{  a fitted \code{inla} object as produced by
##!    \code{inla}.}
##!  \item{x}{a \code{summary.inla} object produced by \code{summary.inla}}
##!  \item{digits}{Integer Number of digits}
##!  \item{include.lincomb}{Logcial Include the summary for the the linear combinations or not}
##!  \item{...}{ other arguments.}
##!
##!}
##!\details{
##!Posterior mean and standard deviation (together with quantiles or
##!cdf) are printed for the fixed effects in the model.
##!
##!For the random effects the function \code{summary()} prints the
##!posterior mean and standard deviations for the hyperparameters
##!}
##!\value{
##!  \code{summary.inla} returns an object of call \code{summaryinla}, a
##!  list with components: 
##!  \item{call}{the component from \code{object}.}
##!  \item{fixed}{the component from \code{object}.}
##!  \item{random}{the component from \code{object}.}
##!  \item{neffp}{the component from \code{object}.}          
##!  \item{linear.predictor}{the component from \code{object}.}
##!  \item{lincomb}{the component from \code{object}.}
##!  \item{lincomb.derived}{the component from \code{object}.}
##!  \item{family}{the component from \code{object}.}   
##!}
##!\author{Sara Martino and Havard Rue}
##!\seealso{ \code{\link{inla}} }

`summary.inla` = function(object, digits = 3L, include.lincomb = TRUE, ...)
{
    inla.eval.dots(...)
    
    ## provides a summary for a inla object
    ret = list()

    maxlen = 2048L
    if (sum(nchar(object$call)) > maxlen) {
        ret = c(ret,
                list(call=paste(substr(inla.paste(object$call), 1L, maxlen),
                             "<<<the rest is not shown>>>")))
    } else {
        ret = c(ret, list(call=object$call))
    }

    ## might not be if using collect directly
    if (inla.is.element("cpu.used", object)) {
        ret = c(ret, list(cpu.used =
                              paste0(names(object$cpu.used)[1], " = ",
                                     format(object$cpu.used[1], digits = digits), ", ", 
                                     names(object$cpu.used)[2], " = ",
                                     format(object$cpu.used[2], digits = digits), ", ", 
                                     names(object$cpu.used)[3], " = ",
                                     format(object$cpu.used[3], digits = digits), ", ",
                                     names(object$cpu.used)[4], " = ",
                                     format(object$cpu.used[4], digits = digits))))
    } else {
        ret = c(ret,  list(cpu.used = NA))
    }

    if(!is.null(object$summary.fixed) && length(object$summary.fixed) > 0) 
        ret = c(ret, list(fixed=round(as.matrix(object$summary.fixed), digits = digits)))

    if (include.lincomb) {
        if(!is.null(object$summary.lincomb) && any(names(object) == "summary.lincomb")
           && (length(object$summary.lincomb) > 0)) 
            ret = c(ret, list(lincomb=round(as.matrix(object$summary.lincomb), digits = digits)))
        if(!is.null(object$summary.lincomb.derived) && length(object$summary.lincomb.derived) > 0)
            ret = c(ret, list(lincomb.derived=round(as.matrix(object$summary.lincomb.derived),
                                  digits = digits)))
    } else {
        if(!is.null(object$summary.lincomb) && any(names(object) == "summary.lincomb")
           && (length(object$summary.lincomb) > 0)) {
            m = nrow(as.matrix(object$summary.lincomb))
            ret = c(ret, list(lincomb=paste("<", m, " lincomb not included>", sep="")))
        }
        if(!is.null(object$summary.lincomb.derived) && length(object$summary.lincomb.derived) > 0) {
            m = nrow(as.matrix(object$summary.lincomb.derived))
            ret = c(ret, list(lincomb.derived=paste("<", m, " lincomb.derived not included>", sep="")))
        }
    }

    if(!is.null(object$summary.hyperpar) && length(object$summary.hyperpar) > 0)
        ret = c(ret, list(hyperpar=round(object$summary.hyperpar, digits = digits)))
    
    if(!is.null(object$summary.random) && length(object$summary.random) > 0) 
        ret = c(ret, list(random.names=names(object$summary.random), random.model=object$model.random))
    
    neffp = object$neffp
    if (!is.null(neffp)) {
        ret = c(ret, list(neffp = round(neffp, digits = digits)))
    }

    if (!is.null(object$dic)) {
        ret = c(ret, list(dic = lapply(object$dic, round, digits = digits)))
    }

    if (!is.null(object$waic)) {
        ret = c(ret, list(waic = lapply(object$waic, round, digits = digits)))
    }

    if(!is.null(object$mlik))
        ret = c(ret, list(mlik = round(object$mlik, digits = digits)))
    
    if(!is.null(object$cpo$cpo) && length(object$cpo$cpo) > 0L)
        ret = c(ret, list(cpo = lapply(object$cpo, round, digits = digits)))

    if(!is.null(object$summary.linear.predictor))
        ret = c(ret, list(linear.predictor= round(object$summary.linear.predictor, digits = digits)))
    
    ret = c(ret, list(family=object$family))
    class(ret) = "summary.inla"

    return (ret)
}

`print.summary.inla` = function(x, digits = 3L, ...)
{
    form = strwrap(inla.formula2character(x$call))
    cat("\nCall:\n")
    for (i in seq_along(form)) {
        cat("  ", form[i], "\n")
    }

    if (inla.is.element("cpu.used",  x)) {
        cat("Time used:\n", "  ", x$cpu.used, "\n")
    }
    
    if (inla.is.element("fixed", x)) {
        cat("Fixed effects:\n")
        print.default(x$fixed)
        cat("\n")
    } else 
        cat("The model has no fixed effects\n\n")

    if (inla.is.element("lincomb.derived", x)) {
        cat("Linear combinations (derived):\n")
        if (is.character(x$lincomb.derived)) {
            cat("\t", x$lincomb.derived)
        } else {
            print.default(x$lincomb.derived)
        }
        cat("\n")
    }
    
    if (inla.is.element("lincomb", x)) {
        cat("Linear combinations:\n")
        if (is.character(x$lincomb)) {
            cat("\t", x$lincomb)
        } else {
            print.default(x$lincomb)
        }
        cat("\n")
    }

    if (inla.is.element("random.names", x)) {
        cat("Random effects:\n")
        cat("  Name\t ", "Model\n ")
        for(i in 1:length(x$random.names))
            cat("  ", paste0(inla.nameunfix(x$random.names[i])," ", x$random.model[i], "\n"))
        cat("\n")
    } else 
        cat("The model has no random effects\n\n")

    if (inla.is.element("hyperpar", x)) {
        cat("Model hyperparameters:\n")
        print(format(x$hyperpar, digits=digits, nsmall=2), quote=FALSE)
        cat("\n")
    } else 
        cat("The model has no hyperparameters\n\n")
    
    if (inla.is.element("neffp", x)) {
        cat("Expected number of effective parameters(stdev): ", format(x$neffp[1], digits=digits, nsmall=2),"(",
            format(x$neffp[2], digits=digits, nsmall=2),")\n", sep="")
        cat("Number of equivalent replicates :", format(x$neffp[3], digits=digits, nsmall=2),"\n\n")
    } else {
        cat("Expected number of effective parameters and\n",
            "number of equivalent replicates are not computed\n\n")
    }

    if (inla.is.element("dic", x)) 
        cat(paste("Deviance Information Criterion (DIC) ...............: ",
                  format(x$dic$dic, digits=digits, nsmall=2), "\n", 
                  "Deviance Information Criterion (DIC, saturated) ....: ",
                  format(x$dic$dic.sat, digits=digits, nsmall=2), "\n", 
                  "Effective number of parameters .....................: ",
                  format(x$dic$p.eff, digits=digits, nsmall=2), "\n\n", sep=""))

    if (inla.is.element("waic", x)) 
        cat(paste("Watanabe-Akaike information criterion (WAIC) ...: ",
                  format(x$waic$waic, digits=digits, nsmall=2), "\n", 
                  "Effective number of parameters .................: ",
                  format(x$waic$p.eff, digits=digits, nsmall=2), "\n\n", sep=""))
     
    if (inla.is.element("mlik", x)) 
        cat(paste("Marginal log-Likelihood: ", format(x$mlik[2], digits=digits, nsmall=2),"\n"))

    if (inla.is.element("cpo", x)) 
        cat("CPO and PIT are computed\n\n")

    if (inla.is.element("po", x)) 
        cat("PO is computed\n\n")

    if (inla.is.element("linear.predictor", x)) 
        cat("Posterior marginals for the linear predictor and\n",
            "the fitted values are computed\n\n")
}
