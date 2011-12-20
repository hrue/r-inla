## a proof-of-concent for Finns idea of removing NA's in the data when
## optimising.

`inla2` = function(formula, family, data, control.inla, ...)
{
    stopifnot(!missing(formula))
    stopifnot(!missing(family))
    stopifnot(!missing(data))

    ## simplifies stuff later
    if (missing(control.inla)) {
        control.inla = list()
    }

    y.name = as.character(formula[2])
    data.new = NULL
    if (is.data.frame(data)) {
        idx = which(names(data) == y.name)
        idx.na = which(is.na(data[, idx]) == TRUE)
        if (length(idx.na) > 0L) {
            data.new = data[-idx.na, ]
        }
    } else {
        stop("data must be a data.frame")
    }

    if (is.null(data.new)) {
        return (inla(formula, family, data, control.inla, ...))
    } else {
        print("pass 1...")
        result = inla(formula, family, data.new,
                ## override the users choice
                control.inla = list(int.strategy = "eb",
                        strategy = "gaussian"),
                ...)
        print(paste("found mode", paste(result$mode$theta, collapse=", ")))
        print("pass 2")
        return (inla(formula, family, data,
                     control.mode = list(theta = result$mode$theta, 
                             restart = FALSE),
                     control.inla, 
                     ...))
    }
}

##

n = 100
N = n
Ntot = n + N
x = runif(Ntot)
y = 1 + x + rnorm(Ntot, sd=0.1)
y[(n+1):Ntot] = NA

i = 1:Ntot
formula = y ~ 1 + x + f(i)
result = inla2(formula, "gaussian", data.frame(y=y, x=x, i=i), verbose=F)
