## Export: print!inla

##!\name{print.inla}
##!\alias{print.inla}
##!\title{Print a INLA fit }
##!\description{
##!Print a INLA fit
##!}
##!\usage{
##!\method{print}{inla}(x, digits = 3L, ...)
##!}
##!\arguments{
##!  \item{x}{An inla-object (output from an \code{\link{inla}}-call).}
##!  \item{digits}{Number of digits to print}
##!  \item{...}{ other arguments.}
##!}
##!\details{
##! None
##!}
##!\value{
##! None
##!}
##!\author{Havard Rue}
##!
##!\seealso{ \code{\link{inla}} }
##!\examples{
##!## None
##!}

`print.inla` = function(x, digits = 3L, ...)
{
    nsmall = 2L
    form = strwrap(inla.formula2character(x$call))
    cat("\nCall:\n")
    for (i in seq_along(form)) {
        cat("  ", form[i], "\n")
    }

    cat("Time used:\n  ")
    if (inla.is.element("cpu.used", x)) {
        cat(sep = "", 
            names(x$cpu.used)[1], " = ",
            format(x$cpu.used[1], digits = digits), ", ", 
            names(x$cpu.used)[2], " = ",
            format(x$cpu.used[2], digits = digits), ", ", 
            names(x$cpu.used)[3], " = ",
            format(x$cpu.used[3], digits = digits), ", ",
            names(x$cpu.used)[4], " = ",
            format(x$cpu.used[4], digits = digits), "\n")
    }

    cat("\nIntegration Strategy: ")
    
    if(is.null(x$.args$control.inla$int.strategy))
        cat("default\n\n ")
    else if(x$.args$control.inla$int.strategy=="eb")
        cat("Empirical Bayes\n\n")
    else if(x$.args$control.inla$int.strategy=="ccd")
        cat("Central Composit Design\n\n")
    else if(x$.args$control.inla$int.strategy=="grid") {
        cat("Integration on a regular grid\n")
        cat("  ",
            "with parameters dz=",
            format(x$.args$control.inla$dz, digits=digits, nsmall = nsmall), 
            " and diff.logdens=",
            format(x$.args$control.inla$diff.logdens,digits = digits, nsmall = nsmall),
            "\n\n")
    }
    
    cat(paste("Model contains ", x$nhyper," hyperparameters\n", sep=""))

    if(!is.null(x$summary.fixed) && !(is.list(x$summary.fixed) && length(x$summary.fixed) == 0)) {
        cat(paste("The model contains ", dim(x$summary.fixed)[1],
                  " fixed effect (including a possible intercept)\n\n", sep=""))
    } else {
        cat("The model has no fixed effects and no intercept\n\n")
    }
    
    ## must check x$control.family that its a list of lists, so we can
    ## treat n.family = 1 the same was as for n.family > 1.
    if (length(x$.args$family) == 1) {
        contfamily = list(x$.args$control.family)
    } else {
        contfamily = x$.args$control.family
    }
    if (length(x$.args$family) > 0L) {
        for(ii in 1:length(x$.args$family)) {
            cat(paste("Likelihood model",
                      inla.ifelse(length(x$.args$family)>1, paste("[", ii, "]", sep=""), ""),
                      ": ", x$.args$family[ii],"\n", sep=""))
        }
        cat("\n")
    }
    if(!is.null(x$summary.random) &&
       !(is.list(x$summary.random) && length(x$summary.random) == 0)) {
        labels <- names(x$summary.random)
        cat(paste("The model has ", length(labels)," random effects:\n", sep=""))
        for(i in 1:length(labels)) {
            cat(paste(i,".'", inla.nameunfix(labels[i]),"' is a " , x$model.random[i],"\n", sep=""))
            cat("\n")
        }
    } else {
        cat("The model has no random effects\n\n")
    }

    return(invisible())
}
