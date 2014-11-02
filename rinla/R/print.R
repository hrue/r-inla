## Export: print!inla

##!\name{print.inla}
##!\alias{print.inla}
##!\title{Print a INLA fit }
##!\description{
##!Print a INLA fit
##!}
##!\usage{
##!\method{print}{inla}(x,...)
##!}
##!\arguments{
##!  \item{x}{An inla-object (output from an \code{\link{inla}}-call).}
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

`print.inla` = function(x, ...)
{
    cat("\nCall:\n", inla.formula2character(x$call), "\n\n", sep = "")
    cat("Time used:\n")
    print(x$cpu.used)
    cat("\nIntegration Strategy: ")
   
    if(is.null(x$.args$control.inla$int.strategy))
        cat("default\n\n ")
    else if(x$.args$control.inla$int.strategy=="eb")
        cat("Empirical Bayes\n\n")
    else if(x$.args$control.inla$int.strategy=="ccd")
        cat("Central Composit Design\n\n")
    else if(x$.args$control.inla$int.strategy=="grid") {
        cat("Integration on a regular grid\n")
        cat("with parameters dz=", x$.args$control.inla$dz,
            " and diff.logdens=", x$.args$control.inla$diff.logdens,"\n\n")
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
            prop = inla.model.properties(x$.args$family[ii], "likelihood")
            ntheta = length(prop$hyper)
            if (ntheta > 0) {
                for (i in 1:ntheta) {
                    ## need a fix here for the 'numeric(0)' problem.
                    if (!is.null(contfamily[[ii]]$hyper[[i]]$fixed)) {
                        if (contfamily[[ii]]$hyper[[i]]$fixed) {
                            cat(paste("\t", prop$hyper[[i]]$name, "in the", x$family[ii], "likelihood is fixed\n"))
                        } else {
                            cat(paste("\t", prop$hyper[[i]]$name, "in the", x$family[ii], "likelihood is random\n"))
                        }
                    }
                }
            }
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
