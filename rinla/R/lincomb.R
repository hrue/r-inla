## utility-functions for creating linear combinations.

##!\name{make.lincomb}
##!\alias{make.lincomb}
##!\alias{make.lincombs}
##!\alias{inla.make.lincomb}
##!\alias{inla.make.lincombs}
##!\description{Create a linear combination or several linear combinations, as input to \code{inla}}
##!\usage{
##!lincomb = inla.make.lincomb(...)
##!lincombs = inla.make.lincombs(...)
##!}
##!\arguments{
##!TODO
##!}
##!\value{
##TODO
##!}
##!\author{Havard Rue \email{hrue@math.ntnu.no}}
##!\seealso{TODO}
##!\examples{
##!See the worked out examples in the FAQ section on \code{r-inla.org}
##!}

`inla.make.lincomb` = function(...)
{
    arg = match.call(expand.dots=TRUE)

    ## we might need to expand arguments?
    if (all(is.null(names(arg))))
        arg = eval.parent(parse(text=arg[2]))

    f.arg = list()
    f.arg[[ length(arg) ]] = list
    for(k in 1:length(arg)) {
        var = names(arg)[k]
        if (var != "") {
            value = eval.parent(arg[[k]])
        
            if (length(value) == 1) {
                ff.arg = list(list(weight = value))
            } else {
                ii = which( !is.na(value) )
                ff.arg = list(list(idx = ii, weight = value[ii]))
            }
        
            names(ff.arg) = var
            f.arg[[k]] = ff.arg
        }
    }
    if (is.null(f.arg[[1]][[1]]))
        f.arg[[1]] = NULL

    lc = list(f.arg)
    names(lc) = "lc"
    
    return (lc)
}

`inla.make.lincombs` = function(...)
{
    ## this is the more general version, which constructs one lincomb
    ## for each row
    
    arg = match.call(expand.dots=TRUE)

    ## we might need to expand arguments?
    if (all(is.null(names(arg))))
        arg = eval(parse(text=arg[2], 2))

    ## check that all 'rows' has the same lengths. store the evalued
    ## args in values[[.]]
    values = list()
    is.m = numeric(length(arg))
    n = -1
    for(k in 1:length(arg)) {
        if (names(arg)[k] != "") {
            values[[k]] = eval.parent(arg[[k]])
            if (is.matrix(values[[k]]) || is(values[[k]], "dgTMatrix")) {
                stopifnot(nrow(values[[k]]) == n || n < 0)
                n = nrow(values[[k]]) 
                is.m[k] = TRUE
            } else {
                stopifnot(length(values[[k]]) == n || n < 0)
                n = length(values[[k]])
                is.m[k] = FALSE
            }
        }
    }
    ## make the nameing so it is the correct order, like lc001, lc002,
    ## ... etc.
    name = paste("lc", inla.num(1:n), sep="")

    ## create one lincomb for each row. faster if dynamic arrays are
    ## allocated upfront.
    lc = list()
    lc[[n]] = list

    if (TRUE) {

        ## this is the fast version, which is kind of ugly, therefore
        ## the ``slow'' code is below. we simply pass all additional
        ## arguments in the reminder of lapply. There is obviously
        ## something in the memory management of R that I don't
        ## understand....

        lc = lapply(1:n,
                function(idx, ...) {
                    ## preallocate
                    f.arg = list()
                    f.arg[[ length(arg) ]] = list
                    for(k in 1:length(arg)) {
                        if (names(arg)[k] != "") {
                            var = names(arg)[k]
                            if (is.m[k]) {
                                if (is(values[[k]], "dgTMatrix")) {
                                    row = inla.sparse.get(values[[k]], row=idx)
                                    ff.arg = list(list(idx = row$j, weight = row$values))
                                } else {
                                    value = values[[k]][idx,]
                                    ii = which( !is.na(value) )
                                    ff.arg = list(list(idx = ii, weight = value[ii]))
                                }
                            } else {
                                ff.arg = list(list(weight = values[[k]][idx]))
                            }
                            names(ff.arg) = var
                            f.arg[[k]] = ff.arg
                        }
                    }
                    ## the first might or might not be relevant
                    if (is.null(f.arg[[1]][[1]]))
                        f.arg[[1]] = NULL

                    return (f.arg)
                },
                arg=arg, is.m=is.m, values=values)
        names(lc) = name
    } else {

        ## this is the slow version for which the lapply-version is
        ## buildt upon.
        
        for(idx in 1:n) {

            f.arg = list()
            for(k in 1:length(arg)) {
                if (names(arg)[k] != "") {
                    var = names(arg)[k]
                    if (is.m[k]) {
                        value = values[[k]][idx,]
                        ii = which( !is.na(value) )
                        ff.arg = list(list(idx = ii, weight = value[ii]))
                    } else {
                        ff.arg = list(list(weight = values[[k]][idx]))
                    }
                    names(ff.arg) = var
                    f.arg = c(f.arg, ff.arg)
                }
                if (is.null(f.arg[[1]][[1]]))
                    f.arg[[1]] = NULL
            }
            lc[[idx]] = f.arg
            names(lc)[idx] = name[idx]
        }
    }

    return (lc)
}

`inla.uncbind` = function(A, name.prefix="col")
{
    ## given a matrix with colnames, expand this into
    ## list(a=c1, b=c2, ....)

    if (!is.matrix(A))
        return (NULL)
    
    result = list()

    ## it's easier to check names of the matrix upfront...
    rownames(A) = NULL
    if (is.null(colnames(A))) {
        colnames(A) = 
            paste(name.prefix,
                  inla.num(1:ncol(A), width = inla.numlen(ncol(A))),
                  sep="")
    } else {
        ii = which(is.na(colnames(A)))
        if (length(ii) > 0) {
            colnames(A)[ii] = paste(name.prefix,
                               inla.num(ii, width = inla.numlen(ncol(A))),
                               sep="")
        }
    }

    if (TRUE) {
        ## this is the fast version
        tmp.result = apply(A, 2, function(x) list(x))
        result = sapply(tmp.result, function(x) c(x))
    } else {
        ## and this is the slow one
        for(icol in 1:ncol(A)) {
            ## extract each column and give it its name or 'column001' if
            ## its NULL.
            x = list(A[, icol])
            name = colnames(A)[icol]
            names(x) = name
            result = c(result, x)
        }
    }

    return (result)
}

`inla.unrbind` = function(A, name.prefix = "row")
{
    return (inla.uncbind(t(A), name.prefix = name.prefix))
}

##########################################################################
 ### OLD VERSIONS GOES BELOW HERE! I THINK THE FORMAT HAS CHANGED... ###
##########################################################################

`inla.make.lincombs.OLD` = function(...)
{
    stop("might not work... please check format!")

    ## this is the more general version, which constructs one lincomb
    ## for each 'row', like
    ##
    ##> A
    ## [, 1] [, 2] [, 3]
    ##[1,]    1    2    3
    ##[2,]    2    1    3
    ##> inla.make.lincombs(b=1:2, a=A, c=1:2)
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
    name = paste("lc", inla.num(1:n), sep="")
            
    ## create one lincomb for each row.
    lc = list()
    for(idx in 1:n) {
        f.arg = list()
        for(var in names(arg)) {
            if (var != "") {
                value = eval(inla.eval(paste("arg$", var, sep="")))
                if (is.matrix(value)) {
                    ii = which( !is.na(value[idx,]) )
                    val = rep(NA, max(ii))
                    val[ii] = value[idx, ii]
                } else {
                    val = value[idx]
                }
                if (TRUE) {
                    val = list(c(val))
                    names(val) = var
                } else {
                    val = inla.eval(inla.paste(c("list(", var, "=", inla.2list(val), ")")))
                }
                f.arg = c(f.arg, val)
            }
        }
        lincomb = inla.make.lincomb.OLD(f.arg)
        names(lincomb) = name[idx]
        lc = c(lc, lincomb)
    }

    return (lc)
}

`inla.make.lincomb.OLD` = function(...)
{
    stop("might not work... please check format")

    ## makes a lincomb-entry using the variables in the call, like
    ## > inla.make.lincomb(x=runif(2), b=1, c=3, Predictor=c(NA, 2, 3))
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
                if (TRUE) {
                    str = paste(str, " ", var)
                    str = paste(str, paste(rbind(idx, value[idx]), collapse=" "))
                } else {
                    str = paste(str, " ", var)
                    for(i in idx)
                        str = paste(str, i, value[i])
                }
            }
        }
    }
    return (str)
}
