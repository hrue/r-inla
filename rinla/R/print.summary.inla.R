`print.summary.inla` = function(x, ...)
{
    digits = 4
    cat("\nCall:\n", inla.formula2character(x$call), "\n\n", sep = "")
    if (inla.is.element("cpu.used",  x)) {
        cat("Time used:\n")
        print(x$cpu.used)
        cat("\n")
    }
    
    if (inla.is.element("fixed", x)) {
        cat("Fixed effects:\n")
        print.default(x$fixed)
        cat("\n")
    } else {
        cat("The model has no fixed effects\n\n")
    }

    if (inla.is.element("lincomb.derived", x)) {
        cat("Linear combinations (derived):\n")
        print.default(x$lincomb.derived)
        cat("\n")
    }

    ## auto-tab error if not....
    if (inla.is.element("lincomb", x)) {
        cat("Linear combinations:\n")
        print.default(x$lincomb)
        cat("\n")
    }

    if (inla.is.element("random.names", x)) {
        cat("Random effects:\n")
        cat("Name\t ", "Model\n ")
        for(i in 1:length(x$random.names))
            cat(paste(inla.nameunfix(x$random.names[i])," ", x$random.model[i],"\n"))
        cat("\n")
    } else {
        cat("The model has no random effects\n\n")
    }

    if (inla.is.element("hyperpar", x)) {
        cat("Model hyperparameters:\n")
        print(format(x$hyperpar, digits=digits, nsmall=2), quote=FALSE)
        cat("\n")
    } else {
        cat("The model has no hyperparameters\n\n")
    }
    
    if (inla.is.element("neffp", x)) {
        cat("Expected number of effective parameters(std dev): ", format(x$neffp[1], digits=digits, nsmall=2),"(",
            format(x$neffp[2], digits=digits, nsmall=2),")\n", sep="")
        cat("Number of equivalent replicates :", format(x$neffp[3], digits=digits, nsmall=2),"\n")
        cat("\n")
    } else {
        cat("Expected number of effective parameters and Number of equivalent replicates not computed\n")
        cat("\n")
    }

    if (inla.is.element("dic", x)) {
        cat(paste("Deviance Information Criterion: ", format(x$dic$dic, digits=digits, nsmall=2), "\n", 
                  "Effective number of parameters: ", format(x$dic$p.eff, digits=digits, nsmall=2), "\n", sep=""))
        cat("\n")
    }
     
    if (inla.is.element("mlik", x)) {
        cat(paste("Marginal Likelihood: ", format(x$mlik[2], digits=digits, nsmall=2),"\n"))
    }

    if (inla.is.element("cpo", x)) {
        cat("CPO and PIT are computed\n")
        cat("\n")
    }

    if (inla.is.element("linear.predictor", x)) {
        cat("Posterior marginals for linear predictor and fitted values computed\n\n")
    }
}
