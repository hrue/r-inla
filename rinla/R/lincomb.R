## utility-functions for creating linear combinations.

inla.make.lincomb = function(...)
{

    ## makes a lincomb-entry using the variables in the call, like
    ## > inla.make.lincomb(x=runif(2), b=1, c=3, Predictor=c(NA,2,3))
    ## [1] " x 1 0.736 2 0.231 b 1 1 c 1 3 Predictor 2 2 3 3"

    str = ""
    arg = match.call(expand.dots=TRUE)

    ## we might need to expand arguments?
    if (all(is.null(names(arg))))
        arg = eval.parent(parse(text=arg[2]))

    for(var in names(arg)) {
        if (var != "") {
            value = eval(inla.eval(paste("arg$", var, sep="")))
            idx = which(!is.na(value))
            if (length(idx) > 0) {
                str = paste(str, " ", var)
                for(i in idx)
                    str = paste(str, i, value[i])
            }
        }
    }
    return (str)
}

inla.make.lincombs = function(...)
{
    ## this is the more general version, which constructs one lincomb
    ## for each 'row', like
    ##
    ##> A
    ## [,1] [,2] [,3]
    ##[1,]    1    2    3
    ##[2,]    2    1    3
    ##> inla.make.lincombs(b=1:2,a=A,c=1:2)
    ##$lincomb.1
    ##[1] "   b 1 1   a 1 1 2 2 3 3   c 1 1"
    ##
    ##$lincomb.2
    ##[1] "   b 2 2   a 1 2 2 1 3 3   c 2 2"
    
    ## unfortunately, this gets quite complicated...

    arg = match.call(expand.dots=TRUE)

    ## we might need to expand arguments?
    if (all(is.null(names(arg))))
        arg = eval(parse(text=arg[2], 2))

    ## check that all 'rows' has the same lengths
    n = -1
    for(var in names(arg)) {
        if (var != "") {
            value = eval(inla.eval(paste("arg$", var, sep="")))
            if (is.matrix(value)) {
                stopifnot(nrow(value) == n || n < 0)
                n = nrow(value) 
            } else {
                stopifnot(length(value) == n || n < 0)
                n = length(value)
            }
        }
    }

    ## make the nameing so it is the correct order, like lc001, lc002,
    ## ... etc.
    name = paste("lc", inla.num(1:n, width = floor(log10(n))+1), sep="")

    ## create one lincomb for each row.
    lc = list()
    for(idx in 1:n) {
        f.arg = list()
        for(var in names(arg)) {
            if (var != "") {
                value = eval(inla.eval(paste("arg$", var, sep="")))
                if (is.matrix(value)) {
                    val = rep(NA, ncol(value))
                    ii = which( !is.na(value[idx,]) )
                    val[ii] = value[idx, ii]
                } else {
                    val = value[idx]
                }
                val = inla.eval(inla.paste(c("list(", var, "=", inla.2list(val), ")")))
                f.arg = c(f.arg, val)
            }
        }
        lincomb = inla.make.lincomb(f.arg)
        names(lincomb) = name[idx]
        lc = c(lc, lincomb)
    }

    return (lc)
}
