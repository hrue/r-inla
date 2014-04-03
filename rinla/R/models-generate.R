## Nothing to export

## These functions are all just internal functions for building the
## documentation,  and generates Rd and tex for the models

`inla.models.generate.Rd` = function(file = NULL)
{
    ## this function is used to generate the man-page for the models

    if (!is.null(file))
        sink(file)
    
    tab = "\t"
    tab2 = paste(tab, "\t",  sep="")
    tab3 = paste(tab2, "\t",  sep="")
    tab4 = paste(tab3, "\t",  sep="")
    tab5 = paste(tab4, "\t",  sep="")
    tab6 = paste(tab5, "\t",  sep="")
    tab7 = paste(tab6, "\t",  sep="")
    tab8 = paste(tab7, "\t",  sep="")
    

    cat("
%% DO NOT EDIT!
%% This file is generated automatically from models.R
\\name{inla.models}
\\alias{inla.models}
\\title{Valid models in INLA}
\\description{
This page describe the models implemented in \\code{inla}, divided into sections: ")
    cat(inla.paste(names(inla.models()), sep=", "), ".\n}\n")

    cat("
\\usage{
inla.models()
}\n")
    
    cat("\\value{
", tab, "Valid sections are: ")
    cat(inla.paste(names(inla.models()),  sep=", "),  "\n")
    cat("
", tab, "\\describe{
")

    for(section in names(inla.models())) {
        cat(tab2, "\\item{Section `", section, "'.}{\n",sep="")
        cat(tab3, "Valid models in this section are:\n")
        cat(tab3, "\\describe{\n")
        
        for(model in names(inla.models()[[section]])) {
            cat(tab4, "\\item{Model `", model, "'.}{",  sep="")

            if (section != "prior") {
                nhyper = length(inla.models()[[section]][[model]]$hyper)
                cat("Number of hyperparmeters are ", nhyper, ".\n", sep="")
                if (nhyper > 0) {
                    cat(tab5, "\\describe{\n")
                    
                    h = inla.models()[[section]][[model]]$hyper
                    nhyper = length(h)
                    if (nhyper > 0) {
                        for(nh in 1:nhyper) {
                            nm = names(h)[nh]
                            cat(tab6, "\\item{Hyperparameter `", nm, "'}{\n", sep="")
                            cat(tab7, "\\describe{\n")

                            for(m in names(h[[nm]])) {
                                mval = h[[nm]][[m]]
                                if (is.null(mval))
                                    mval = "NULL"
                                if (is.function(mval)) {
                                    mval.src = attr(mval, "srcref")
                                    if (is.null(mval.src)) {
                                        mval.src = inla.paste(deparse(mval, control="keepInteger"))
                                    }
                                    mval = paste("\\code{", mval.src, "}",  sep="")
                                }
                                cat(tab7, "\\item{", m, " = }{`", inla.paste(mval), "'}\n", sep="")
                            }

                            cat(tab7, "}\n", sep="")
                            cat(tab6, "}\n", sep="")
                        }
                        cat(tab5, "}\n")
                    }

                    h = inla.models()[[section]][[model]]
                    ms = names(h)
                    ms = ms[ ms != "hyper" ]
                    
                    cat(tab5, "\\describe{\n")
                    cat(tab6, "\\item{Properties:}{\n", sep="")
                    cat(tab7, "\\describe{\n", sep="")
                    
                    for(m in ms) {
                        mval = h[[m]]
                        if (is.null(mval))
                            mval = "NULL"
                        cat(tab8, "\\item{", m, " = }{`", inla.paste(mval), "'}\n", sep="")
                    }
                    
                    cat(tab7, "}\n")
                    cat(tab6, "}\n")
                    cat(tab5, "}\n")
                }
            } else if (section == "prior") {
                np = inla.models()[[section]][[model]]$nparameters
                cat(tab5, "Number of parameters in the prior =", np, "\n")
            } else {
                stop("Should not happen.")
            }
                
            cat(tab4, "}\n")
        }
        cat(tab3, "}\n")
        cat(tab2, "}\n")
    }

    cat(tab, "} \n")
    cat("} \n")

    cat("\\examples{
## How to set hyperparameters to pass as the argument 'hyper'. This
## format is compatible with the old style (using 'initial', 'fixed',
## 'prior', 'param'), but the new style using 'hyper' take preceedence
## over the old style. The two styles can also be mixed. The old style
## might be removed from the code in the future...

## Only a subset need to be given
   hyper = list(theta = list(initial = 2))
## The `name' can be used instead of 'theta', or 'theta1', 'theta2',...
   hyper = list(precision = list(initial = 2))
   hyper = list(precision = list(prior = \"flat\", param = numeric(0)))
   hyper = list(theta2 = list(initial=3), theta1 = list(prior = \"gaussian\"))
## The 'short.name' can be used instead of 'name'
   hyper = list(rho = list(param = c(0,1)))
}\n")    
    if (!is.null(file))
        sink()

}

`inla.models.generate.tex` = function(a.list = inla.models(), file = NULL)
{
    my.doit.recursively = function(a.list, level = 0L)
    {
        if (level == 0L) {
            cat("%% DO NOT EDIT!\n")
            cat("%% This file is generated automatically from models.R\n")
            cat("\\begin{description}\n")
        }

        tab = inla.paste(rep("\t", level + 1L))

        for(nm in names(a.list)) {
            if (is.list(a.list[[nm]])) {
                if (length(a.list[[nm]]) > 0) {
                    cat(tab, "\\item[", nm, "]\\ ", "\n", sep="")
                    cat(tab, "\\begin{description}\n")
                    my.doit.recursively(a.list[[nm]], level + 1L)
                    cat(tab, "\\end{description}\n")
                } else {
                    cat(tab, "\\item[", nm, "]\\ ", "\n", sep="")
                }
            } else {
                val = a.list[[nm]]
                if (is.function(val)) {
                    val.src = attr(val, "srcref")
                    if (is.null(val.src)) {
                        val.src = inla.paste(deparse(val, control="keepInteger"))
                    }
                    val = paste("\\verb!", val.src, "!", sep="")
                }               
                cat(tab, "\\item[", nm, "]", " ", inla.paste(val), "\n", sep="")
            }
        }

        if (level == 0L)
            cat("\\end{description}\n")
    }

    ##
    ##
    
    if (!is.null(file))
        sink(file)
    
    my.doit.recursively(a.list)
    
    if (!is.null(file))
        sink()
}
