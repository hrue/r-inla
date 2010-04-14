### RCSId = "$Id: print.summary.inla.R,v 1.8 2010/03/16 12:15:48 hrue Exp $"

`print.summary.inla` <-
function(x,...)
{
    digits = 4
    cat("\nCall:\n", inla.formula2character(x$call), "\n\n", sep = "")
    if (!is.null(x$cpu.used)) {
        cat("Time used:\n")
        print(x$cpu.used)
        cat("\n")
    }
    
    if(!is.null(x$fixed)) {
        cat("Fixed effects:\n")
        print.default(x$fixed)
        cat("\n")
    } else {
        cat("The model has no fixed effects\n\n")
    }

    if(!is.null(x$lincomb)) {
        cat("Linear combinations:\n")
        print.default(x$lincomb)
        cat("\n")
    }

    if(!is.null(x$random.names)) {
        cat("Random effects:\n")
        cat("Name\t ","Model\t ", "\tMax KLD \n"  )
        for(i in 1:length(x$random.names))
            cat(paste(inla.nameunfix(x$random.names[i])," ",x$random.model[i]," ",
                      round(x$random.max.kld[i],5),"\n"))
        cat("\n")
    } else {
        cat("The model has no random effects\n\n")
    }

    if(!is.null(x$hyperpar)) {
        cat("Model hyperparameters:\n")
        print(format(x$hyperpar,digits=digits,nsmall=2),quote=FALSE)
        cat("\n")
    } else {
        cat("The model has no hyperparameters\n\n")
    }
    
    if(!is.null(x$neffp)) {
        cat("Expected number of effective parameters(std dev): ",format(x$neffp[1],digits=digits,nsmall=2),"(",
            format(x$neffp[2],digits=digits,nsmall=2),")\n",sep="")
        cat("Number of equivalent replicates :",format(x$neffp[3],digits=digits,nsmall=2),"\n")
        if(x$neffp[3]<2)
            cat("WARNING: The approximations could be not very accurate\n")
        cat("\n")
    }
    else
        cat("Expected number of effective parameters and Number of equivalent replicates not computed\n")
    cat("\n")

    if(!is.null(x$dic)) {
        cat(paste("Deviance Information Criterion: ",format(x$dic[4],digits=digits,nsmall=2),
                  " Effective number of parameters:",
                  format(x$dic[3],digits=digits,nsmall=2),"\n"))
        cat("\n")
    }
     
    if(!is.null(x$mlik)) {
        cat(paste("Marginal Likelihood: ",format(x$mlik[2],digits=digits,nsmall=2),"\n"))
        cat("Warning: Interpret the marginal likelihood with care if the prior model is improper.\n")
    }

    if(!is.null(x$pit)) {
        cat("CPO and PIT are computed\n")
        cat("\n")
    }
    if(!is.null(x$linear.predictor))
        cat("Posterior marginals for linear predictor and fitted values computed\n\n")
}

