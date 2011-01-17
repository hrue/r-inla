
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
    
    for(ii in 1:length(x$family)) {
        cat(paste("Likelihood model",
                  inla.ifelse(length(x$family)>1, paste("[", ii, "]", sep=""), ""),
                  ": ", x$family[ii],"\n", sep=""))
        prop = inla.lmodel.properties(x$family[ii])
        if (prop$ntheta > 0)
            for (i in 1:prop$ntheta) {
                if(x$control.data[[ii]]$fixed[i])
                    cat(paste("\t", prop$theta[i], "in the", x$family[ii], "likelihood is fixed\n"))
                else 
                    cat(paste("\t", prop$theta[i], "in the", x$family[ii], "likelihood is random\n"))
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
