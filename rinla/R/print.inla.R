
`print.inla` = function(x, ...)
{
    cat("\nCall:\n", inla.formula2character(x$call), "\n\n", sep = "")
    cat("Time used:\n")
    print(x$cpu.used)
    cat("\nIntegration Strategy: ")
   
    if(is.null(x$control.inla$int.strategy))
        cat("DEFAULT\n\n ")
    else if(x$control.inla$int.strategy=="EB")
        cat("Empirical Bayes\n\n")
    else if(x$control.inla$int.strategy=="CCD")
        cat("Central Composit Design\n\n")
    else if(x$control.inla$int.strategy=="GRID") {
        cat("Integration on a regular grid\n")
        cat("with parameters dz=", x$control.inla$dz," and diff.logdens=", x$control.inla$diff.logdens,"\n\n")
    }
    
    cat(paste("Model contains ", x$nhyper," hyperparameters\n", sep=""))

    if(!is.null(x$summary.fixed))
        cat(paste("The model contains ", dim(x$summary.fixed)[1]," fixed effect (including a possible intercept)\n\n", sep=""))
    else
        cat("The model has no fixed effects and no intercept\n\n")
    
    ## must check x$control.data that its a list of lists, so we can
    ## treat n.family = 1 the same was as for n.family > 1.
    if (length(x$family) == 1) {
        contdata = list(x$control.data)
    } else {
        contdata = x$control.data
    }
    for(ii in 1:length(x$family)) {
        cat(paste("Likelihood model",
                  inla.ifelse(length(x$family)>1, paste("[", ii, "]", sep=""), ""),
                  ": ", x$family[ii],"\n", sep=""))
        prop = inla.model.properties(x$family[ii], "likelihood")
        ntheta = length(prop$hyper)
        if (ntheta > 0) {
            for (i in 1:ntheta) {
                ## need a fix here for the 'numeric(0)' problem.
                if (!is.null(contdata[[ii]]$hyper[[i]]$fixed)) {
                    if (contdata[[ii]]$hyper[[i]]$fixed) {
                        cat(paste("\t", prop$hyper[[i]]$name, "in the", x$family[ii], "likelihood is fixed\n"))
                    } else {
                        cat(paste("\t", prop$hyper[[i]]$name, "in the", x$family[ii], "likelihood is random\n"))
                    }
                }
            }
        }
    }
    cat("\n")
  
    if(!is.null(x$summary.random)) {
        labels <- names(x$summary.random)
        cat(paste("The model has ", length(labels)," random effects:\n", sep=""))
        for(i in 1:length(labels)) {
            cat(paste(i,".'", inla.nameunfix(labels[i]),"' is a " , x$model.random[i],"\n", sep=""))
            cat("\n")
        }
    } else {
        cat("The model has no random effects\n\n")
    }
}
