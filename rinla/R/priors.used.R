## Export: inla.priors.used

##!\name{inla.priors.used}
##!\alias{inla.priors.used}
##!\alias{priors.used}
##!\title{Print priors used}
##!\description{Print the priors used for the hyperparameters}
##!\usage{
##!    inla.priors.used(result, digits=6L)
##!}
##!\arguments{
##!\item{result}{An \code{inla}-object, typically the output from an \code{inla()}-call}
##!\item{digits}{The \code{digits} argument to the function \code{format()}}
##!}
##!\details{
##!This function provides a more human-friendly output of \code{result$all.hyper}
##!of all the priors used for the hyperparameters.
##!Since not all information about the model is encoded in this object, more hyperparameters
##!than actually used, may be printed. In particular, \code{group.theta1} is printed
##!even though the argument \code{group} in \code{f()} is not used. Similarly for spde-models, 
##!but the user should know that, for example, only the two first ones
##!are actually used.
##!Hopefully, this issue will be fixed in the future.
##!}
##!\author{Havard Rue \email{hrue@r-inla.org}}
##!\examples{
##!r = inla(y ~ 1 + x, data = data.frame(y = 1:10, x = rep(1:5, 2)))
##!inla.priors.used(r)
##!}

`inla.priors.used` = function(result, digits = 6L) 
{
    p2char = function(p) {
        return (paste(as.character(format(p, digits = digits)), sep="", collapse=", "))
    }

    mvnorm.n = function(len) {
        for(n in 1:300) {
            if (len == n + n^2)
                return (n)
        }
        stop(paste0("This should not happen: Fail to find 'n' from len=", len))
    }

    mvnorm.prior.print = function(param){
        n = mvnorm.n(length(param))
        m = matrix(param[1:n], nrow = n, ncol = 1)
        P = matrix(param[-(1:n)], nrow = n, ncol = n)
        for(i in 1:n) {
            cat("\t\t\t[",  format(m[i], digits = digits), "]  [", sep="")
            for(j in 1:n)
                cat(format(P[i, j], digits = digits), " ", sep="")
            cat("]\n", sep="")
        }
    }

    h4cat = function(h4, pre = NULL) {
        if (h4$fixed == FALSE) {
            cat("\t\t",
                if (!is.null(pre)) paste0(pre, ".") else "",
                "theta", "", j, ":", "\n",
                "\t\t\t", "parameter=[", h4$name, "]", "\n",
                "\t\t\t", "prior=[", h4$prior, "]", "\n",
                "\t\t\t", "param=[", p2char(h4$param), "]", "\n", sep="")
            if (FALSE) {
                if (h4$prior %in% "mvnorm")
                    mvnorm.prior.print(h4$param)
            }
        }
    }

    stopifnot(class(result) == "inla")
    h = result$all.hyper
    ntheta = length(result$mode$theta)
    
    for(nm in names(h)) {
        h2 = h[[nm]]
        if (!(nm %in% c("predictor"))) {
            cat("section=[", nm, "]", "\n", sep="")
            for(i in seq_along(h2)) {
                h3 = h2[[i]]
                if (is.null(h3$label)) h3$label = h3$hyperid
                if (is.null(h3$hyperid)) h3$hyperid = unlist(h3$label)
                cat("\t", "tag=[", h3$hyperid, "] component=[", unlist(h3$label),
                    "]", "\n", sep="")
                if (nm %in% "fixed") {
                    ## this is a special case
                    h3$hyper = list(list(name = h3$label,
                                         prior = "normal",
                                         param = c(h3$prior.mean, h3$prior.prec),
                                         fixed = FALSE))
                }
                for(j in seq_along(h3$hyper)) {
                    if (j <= ntheta) h4cat(h3$hyper[[j]])
                }
                for(j in seq_along(h3$group.hyper)) {
                    if (j <= ntheta) h4cat(h3$group.hyper[[j]], "group")
                }
                for(j in seq_along(h3$link$hyper)) {
                    if (j <= ntheta) h4cat(h3$link$hyper[[j]], "link")
                }
                for(j in seq_along(h3$mix$hyper)) {
                    if (j <= ntheta) h4cat(h3$mix$hyper[[j]], "mix")
                }
            }
        }
    }
}

