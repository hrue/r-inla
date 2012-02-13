
inla.cpo = function(...,  force.cpo.manual = FALSE, be.verbose = TRUE, recompute.mode = TRUE)
{
    ## evaluate inla() using the same arguments, but compute the
    ## cpo-values manually for those with a 'failure'

    if (TRUE) {
        lib = "INLA::"
    } else {
        print("***FIXME*** revert back to INLA::")
        lib = ""
    }

    arg.char = as.character(as.expression(match.call()))
    just.args = gsub(paste("^(", lib, ")?inla.cpo[(]", sep=""), "", arg.char)
    just.args = gsub("[)]$", "", just.args)

    ## remove the local arguments
    just.args = gsub(",[ ]*force\\.cpo\\.manual[ ]*=[^,]*", "", just.args)
    just.args = gsub(",[ ]*be\\.verbose[ ]*=[^,]*", "", just.args)
    just.args = gsub(",[ ]*recompute\\.mode[ ]*=[^,]*", "", just.args)

    ## call the standard inla()
    r = inla.eval(paste(lib, "inla(", just.args, ")", sep=""))

    ## if there is no cpo, then done
    if (is.null(r$cpo))
        return(r)

    ## loop over those with failure > 0
    if (!force.cpo.manual)
        idx.fail = which(r$cpo$failure > 0)
    else
        idx.fail = 1:length(r$cpo$failure)
    if (length(idx.fail) > 0) {
        cpo.old = r$cpo$cpo[idx.fail]
        pit.old = r$cpo$pit[idx.fail]
        k=1
        for(idx in idx.fail) {
            if (be.verbose) {
                cat("\r                                       ", inla.tictac(k))
                cat("\rCompute manually CPO/PIT for ", length(idx.fail), " nodes:\r", sep="")
                k = k+1
            }
            
            argument = paste(sep="", lib, "inla(", just.args, ",",
                    "control.expert = list(cpo.manual = TRUE, cpo.idx =", idx, "),",
                    "control.mode = list(result = r, restart=", inla.ifelse(recompute.mode, "TRUE", "FALSE"), "))")
            rr = inla.eval(argument)
            r$cpo$cpo[idx] = rr$cpo$cpo[idx]
            r$cpo$pit[idx] = rr$cpo$pit[idx]
            r$cpo$failure[idx] = 0
        }

        if (be.verbose)
            cat("\n")

        cpo.new = r$cpo$cpo[idx.fail]
        pit.new = r$cpo$pit[idx.fail]

        if (be.verbose) {
            print(cbind(index = idx.fail,
                        cpo.old = cpo.old, cpo.new = cpo.new,
                        pit.old = pit.old, pit.new = pit.new))
            
            cat("\nThe retured object contain the new values\n")
        }
    }

    return(r)
}
