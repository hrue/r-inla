##!\name{inla.group}
##!\alias{inla.group}
##!
##!\title{Group or cluster covariates}
##!
##!\description{\code{inla.group} group or cluster covariates so to reduce
##!  the number of unique values}
##!
##!\usage{
##!inla.group(x, n = 25, method = c("cut", "quantile"), idx.only = FALSE)
##!}
##!
##!\arguments{
##!
##!  \item{x}{The vector of covariates to group.}
##!
##!  \item{n}{Number of classes or bins to group into.}
##!  
##!  \item{method}{Group either using bins with equal length intervals
##!    (\code{method = "cut"}), or equal distance in the `probability'
##!    scale using the quantiles (\code{method = "quantile"}).}
##!
##!  \item{idx.only}{Option to return the index only and not the
##!    \code{method}.}
##!
##!}
##!
##!
##!\value{
##!  \code{inla.group} return the new grouped covariates where the classes
##!  are set to the median of all the covariates belonging to that group.
##!}
##!
##!\author{Havard Rue \email{hrue@math.ntnu.no}}
##!
##!\seealso{\code{\link{f}}}
##!
##!\examples{
##!## this gives groups 3 and 8
##!x = 1:10
##!x.group = inla.group(x, n = 2)
##!
##!## this is the intended use, to reduce the number of unique values in
##!## the of first argument of f()
##!n = 100
##!x = rnorm(n)
##!y = x + rnorm(n)
##!result = inla(y ~ f(inla.group(x, n = 20), model = "iid"), data=data.frame(y=y,x=x))
##!}

`inla.group` = function(x, n = 25, method = c("cut", "quantile"), idx.only = FALSE)
{
    `inla.group.core` = function(x, n = 25, method = c("cut", "quantile"), idx.only)
    {
        ## group covariates into N groups using method "quantile" or
        ## "cut", i.e., the functions quantile() or cut().  the cut
        ## use 'even length' wheras the 'quantile' use even length
        ## quantiles.

        ## I make the "cut" default, as then we have control over the
        ## minimum distance between each cell, whereas the "quantile"
        ## approach, we do not have a such control.

        if (n < 1)
            stop("Number of groups must be > 0")

        if (n == 1)
            return (rep(median(x), length(x)))

        method = match.arg(method)
        if (method == "cut") {
            ## use equal length
            a = cut(x, n)
        } else {
            ## use break-points corresponding to the quantiles
            aq = unique(quantile(x, probs = c(0, ppoints(n-1), 1)))
            a = cut(x, breaks = as.numeric(aq), include.lowest=TRUE)
        }
        ## the rest is then the same
        nlev = nlevels(a)
        xx = list()
        for(i in 1:nlev)
            xx[[i]] = list()

        for(i in 1:length(x))
            xx[[as.numeric(a[i])]] = c(unlist(xx[[as.numeric(a[i])]]), x[i])
        values = numeric(nlev)

        ff.local = function(xx) {
            if (length(xx) > 0)
                return (median(xx))
            else
                return (NA)
        }
        if (!idx.only) {
            values = unlist(sapply(xx, ff.local))
            return (as.numeric(values[as.numeric(a)]))
        } else {
            return (as.numeric(a))
        }
    }

    if (missing(x))
        return (NULL)

    if (any(is.na(x))) {
        idx.ok = !is.na(x)
        x[idx.ok] = inla.group.core(x[idx.ok], n, method, idx.only)

        return (x)
    } else {
        return (inla.group.core(x, n, method, idx.only))
    }
}

`inla.group.old` = function(x, n)
{
    ### old version
    cutpoints = seq(range(x)[1], range(x)[2], length.out=(n+1))
    lev = (cutpoints[1:(n)]+cutpoints[2:(n+1)])/2
    int = cut(x, breaks=cutpoints, include.lowest=TRUE, labels=FALSE)
    return(lev[int])
}
