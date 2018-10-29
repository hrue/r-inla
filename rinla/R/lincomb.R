## Export: inla.make.lincomb inla.make.lincombs

##!\name{make.lincomb}
##!\alias{make.lincomb}
##!\alias{make.lincombs}
##!\alias{inla.make.lincomb}
##!\alias{inla.make.lincombs}
##!\title{Create linear combinations}
##!\description{Create a linear combination or several
##!linear combinations, as input to \code{inla(...,  lincomb = <lincomb>)}}
##!\usage{
##!inla.make.lincomb(...)
##!inla.make.lincombs(...)
##!}
##!\arguments{
##!   \item{...}{Arguments; see examples}
##!}
##!\value{
##! A structure to be passed on to \code{\link{inla}} argument \code{lincomb} 
##!}
##!\author{Havard Rue \email{hrue@r-inla.org}}
##!\seealso{TODO}
##!\examples{
##!##See the worked out examples and description in the FAQ
##!##section on {www.r-inla.org}
##!}

`inla.make.lincomb` = function(...)
{
    arg = match.call(expand.dots=TRUE)

    ## we might need to expand arguments?
    if (all(is.null(names(arg)))) {
        arg = eval(parse(text=arg[2]), envir = parent.frame(), enclos = environment(list(...)))
    }

    f.arg = list()
    f.arg[[ length(arg) ]] = list
    for(k in 1:length(arg)) {
        var = names(arg)[k]
        if (var != "") {
            value = eval(arg[[k]], envir = parent.frame(), enclos = environment(list(...)))

            if (length(value) == 1) {
                ff.arg = list(list(weight = value))
            } else {
                ii = which(!is.na(value) & (value != 0))
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
    if (all(is.null(names(arg)))) {
        arg = eval(parse(text=arg[2]), envir = parent.frame(),
                enclos = environment(list(...)))
    }

    ## check that all 'rows' has the same lengths. store the evalued
    ## args in values[[.]]
    values = list()
    is.m = numeric(length(arg))
    n = -1
    for(k in 1:length(arg)) {
        if (names(arg)[k] != "") {
            values[[k]] = eval(arg[[k]], envir = parent.frame(),  enclos = environment(list(...)))
            if (is.matrix(values[[k]]) || is(values[[k]], "Matrix")) {
                if (is(values[[k]], "Matrix")) {
                    ## Make sure we have a dgTMatrix with unique representation:
                    values[[k]] = inla.as.dgTMatrix(values[[k]])
                }
                stopifnot(nrow(values[[k]]) == n || n < 0)
                n = nrow(values[[k]])
                is.m[k] = TRUE
            } else {
                stopifnot(length(values[[k]]) == n || n < 0)
                values[[k]][ is.na(values[[k]]) ] = 0
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
                                    ii = which(!is.na(value) & (value != 0))
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
                        ii = which(!is.na(value) & (value != 0))
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
