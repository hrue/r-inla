## Export: inla.posterior.sample inla.posterior.sample.eval

## ! \name{inla.sample}
## ! \alias{inla.posterior.sample}
## ! \alias{posterior.sample}
## ! \alias{inla.posterior.sample.eval}
## ! \alias{posterior.sample.eval}
## !
## ! \title{Generate samples, and functions thereof, from an approximated posterior of a fitted model}
## !
## ! \description{This function generate samples, and functions of those,
## !              from an approximated posterior of a fitted model (an inla-object)}
## ! \usage{
## !     inla.posterior.sample(n = 1L, result, selection = list(),
## !                           intern = FALSE,
## !                           use.improved.mean = TRUE, skew.corr = TRUE,
## !                           add.names = TRUE, seed = 0L, num.threads = NULL,
## !                           parallel.configs = TRUE,  verbose=FALSE)
## !     inla.posterior.sample.eval(fun, samples, return.matrix = TRUE, ...)
## ! }
## !
## ! \arguments{
## !   \item{n}{Number of samples.}
## !   \item{result}{The inla-object, ie the output from an \code{inla}-call.
## !       The \code{inla}-object must be created with
## !       \code{control.compute=list(config=TRUE)}.}
## !   \item{selection}{Select what part of the sample to return. By default, the whole sample
## !       is returned. \code{selection} is a named list with the name of the components of
## !       the sample, and what indices of them to return. Names include \code{APredictor},
## !       \code{Predictor}, \code{(Intercept)},  and otherwise names in the formula.
## !       The values of the list, is interpreted as indices. If they
## !       are negative, they are interpreted as 'not', a zero is interpreted as 'all',  and
## !       positive indices are interpreted as 'only'. The names of elements of each samples
## !       refer to the indices in the full sample.
## !       DO NOT USE this feature together with \code{inla.posterior.sample.eval}.}
## !   \item{intern}{Logical. If \code{TRUE} then produce samples in the
## !        internal scale for the hyperparmater, if \code{FALSE} then produce
## !        samples in the user-scale. (For example log-precision (intern)
## !        and precision (user-scale))}
## !   \item{use.improved.mean}{Logical. If \code{TRUE} then use the
## !        marginal mean values when constructing samples. If \code{FALSE}
## !        then use the mean in the Gaussian approximations.}
## !   \item{skew.corr}{Logical. If \code{TRUE} then correct samples for skewness,
## !       if \code{FALSE},  do not correct samples for skewness (ie use the
## !       Gaussian).}
## !   \item{add.names}{Logical. If \code{TRUE} then add name for each elements of each
## !       sample. If \code{FALSE}, only add name for the first sample.
## !       (This save space.)}
## !   \item{seed}{See the same argument in \code{?inla.qsample} for further
## !               information. In order to produce reproducible results,  you
## !               ALSO need to make sure the RNG in R is in the same state,
## !               see example below.  When \code{seed} is non-zero,
## !               \code{num.threads} is forced to "1:1" and parallel.configs is
## !               set to \code{FALSE}, since parallel sampling would not produce
## !               a reproducible sequence of pseudo-random numbers.}
## !   \item{num.threads}{The number of threads to use in the format 'A:B' defining the number threads in the
## !                      outer (A) and inner (B) layer for nested parallelism. A '0' will be replaced
## !                      intelligently.
## !                      \code{seed!=0} requires serial comptuations.}
## !   \item{parallel.configs}{Logical. If \code{TRUE} and not on Windows,
## !                           then try to run each configuration in
## !                           parallel (not Windows) using \code{A} threads (see \code{num.threads}),
## !                           where each of them is using \code{B:0} threads.}
## !   \item{verbose}{Logical. Run in verbose mode or not.}
## !   \item{fun}{The function to evaluate for each sample. Upon entry, the variable names
## !              defined in the model are defined as the value of the sample.
## !              The list of names are defined in \code{result$misc$configs$contents} where
## !              \code{result} is an \code{inla}-object. This includes predefined names for
## !              for the linear predictor (\code{Predictor} and \code{APredictor}),  and the
## !              intercept (\code{(Intercept)} or \code{Intercept}).
## !              The hyperparameters are defined as \code{theta},  no matter if they are in the
## !              internal scale or not. The function \code{fun} can also return a vector.
## !              To simplify usage, \code{fun} can also be a vector character's. In this case
## !              \code{fun} it is interpreted as (strict) variable
## !              names, and a function is created that return these variables:
## !              if argument \code{fun} equals \code{c("Intercept", "a[1:2]")},  then this is equivalent to
## !              pass \code{function() return(c(get('Intercept'), get('a[1:2]')))}.}
## !   \item{samples}{\code{samples} is the output from \code{inla.posterior.sample()}}
## !   \item{return.matrix}{Logical. If \code{TRUE},  then return the samples of \code{fun}
## !                         as matrix,  otherwise,  as a list.}
## !   \item{...}{Additional arguments to \code{fun}}
## !}
## !\details{The hyperparameters are sampled from the configurations used to do the
## !       numerical integration, hence if you want a higher resolution, you need to
## !       to change the \code{int.stratey} variable and friends. The latent field is
## !       sampled from the Gaussian approximation conditioned on the hyperparameters,
## !       but with a correction for the mean (default),
## !        and optional (and by default) corrected for the estimated skewness.
## !
## !       The log.density report is only correct when there is no constraints.
## !       With constraints, it correct the Gaussian part of the sample for the constraints.
## !
## !       After the sample is (optional) skewness corrected, the log.density
## !       is is not exact for correcting for constraints, but the error is very
## !       small in most cases.
## !}
## !\value{\code{inla.posterior.sample} returns a list of the samples,
## !       where each sample is a list with
## !     names \code{hyperpar} and \code{latent}, and with their marginal
## !     densities in \code{logdens$hyperpar} and \code{logdens$latent}
## !     and the joint density is in \code{logdens$joint}.
## !     \code{inla.posterior.sample.eval} return a list or a matrix of
## !     \code{fun} applied to each sample.
## !}
## !\author{Havard Rue \email{hrue@r-inla.org} and Cristian Chiuchiolo \email{cristian.chiuchiolo@kaust.edu.sa}}
## !
## !\examples{
## !  r = inla(y ~ 1 ,data = data.frame(y=rnorm(1)), control.compute = list(config=TRUE))
## !  samples = inla.posterior.sample(2,r)
## !
## !  ## reproducible results:
## !  inla.seed = as.integer(runif(1)*.Machine$integer.max)
## !  set.seed(12345)
## !  x = inla.posterior.sample(10, r, seed = inla.seed, num.threads="1:1")
## !  set.seed(12345)
## !  xx = inla.posterior.sample(10, r, seed = inla.seed, num.threads="1.1")
## !  all.equal(x, xx)
## !
## ! set.seed(1234)
## ! n = 25
## ! xx = rnorm(n)
## ! yy = rev(xx)
## ! z = runif(n)
## ! y = rnorm(n)
## ! r = inla(y ~ 1 + z + f(xx) + f(yy, copy="xx"),
## !         data = data.frame(y, z, xx, yy),
## !         control.compute = list(config=TRUE),
## !         family = "gaussian")
## ! r.samples = inla.posterior.sample(10, r)
## !
## ! fun = function(...) {
## !     mean(xx) - mean(yy)
## ! }
## ! f1 = inla.posterior.sample.eval(fun, r.samples)
## !
## ! fun = function(...) {
## !     c(exp(Intercept), exp(Intercept + z))
## ! }
## ! f2 = inla.posterior.sample.eval(fun, r.samples)
## !
## ! fun = function(...) {
## !     return (theta[1]/(theta[1] + theta[2]))
## ! }
## ! f3 = inla.posterior.sample.eval(fun, r.samples)
## !
## ! ## Predicting nz new observations, and
## ! ## comparing the estimated one with the true one
## ! set.seed(1234)
## ! n = 100
## ! alpha = beta = s = 1
## ! z = rnorm(n)
## ! y = alpha + beta * z + rnorm(n, sd = s)
## ! r = inla(y ~ 1 + z,
## !         data = data.frame(y, z),
## !         control.compute = list(config=TRUE),
## !         family = "gaussian")
## ! r.samples = inla.posterior.sample(10^3, r)
## !
## ! ## just return samples of the intercept
## ! intercepts = inla.posterior.sample.eval("Intercept", r.samples)
## !
## ! nz = 3
## ! znew = rnorm(nz)
## ! fun = function(zz = NA) {
## !     ## theta[1] is the precision
## !     return (Intercept + z * zz +
## !             rnorm(length(zz), sd = sqrt(1/theta[1])))
## ! }
## ! par(mfrow=c(1, nz))
## ! f1 = inla.posterior.sample.eval(fun, r.samples, zz = znew)
## ! for(i in 1:nz) {
## !     hist(f1[i, ], n = 100, prob = TRUE)
## !     m = alpha + beta * znew[i]
## !     xx = seq(m-4*s, m+4*s, by = s/100)
## !     lines(xx, dnorm(xx, mean=m, sd = s), lwd=2)
## ! }
## !
## ! ## 
## ! ## Be aware that using non-clean variable names might be a little tricky
## ! ## 
## ! n <- 100
## ! X <- matrix(rnorm(n^2), n, 2)
## ! x <- X[, 1]
## ! xx <- X[, 2]
## ! xxx <- x*xx
## ! 
## ! y <- 1 + 2*x + 3*xx + 4*xxx + rnorm(n, sd = 0.01)
## ! 
## ! r <- inla(y ~ X[, 1]*X[, 2],
## !           data = list(y = y, X = X),
## !           control.compute = list(config = TRUE))
## ! print(round(dig = 4, r$summary.fixed[,"mean"]))
## ! 
## ! sam <- inla.posterior.sample(100, r)
## ! sam.extract <- inla.posterior.sample.eval(
## !     (function(...) {
## !         beta.1 <- get("X[, 1]")
## !         beta.2 <- get("X[, 2]")
## !         beta.12 <- get("X[, 1]:X[, 2]")
## !         return(c(Intercept, beta.1, beta.2, beta.12))
## !     }), sam)
## ! print(round(dig = 4, rowMeans(sam.extract)))
## ! 
## ! ## a simpler form can also be used here, and in the examples below
## ! sam.extract <- inla.posterior.sample.eval(
## !                c("Intercept", "X[, 1]", "X[, 2]", "X[, 1]:X[, 2]"), sam)
## ! print(round(dig = 4, rowMeans(sam.extract)))
## !
## ! r <- inla(y ~ x + xx + xxx,
## !           data = list(y = y, x = x, xx = xx, xxx = xxx), 
## !           control.compute = list(config = TRUE))
## ! 
## ! sam <- inla.posterior.sample(100, r)
## ! sam.extract <- inla.posterior.sample.eval(
## !     (function(...) {
## !         return(c(Intercept, x, xx, xxx))
## !     }), sam)
## ! print(round(dig = 4, rowMeans(sam.extract)))
## !
## ! sam.extract <- inla.posterior.sample.eval(c("Intercept", "x", "xx", "xxx"), sam)
## ! print(round(dig = 4, rowMeans(sam.extract)))
## !
## ! r <- inla(y ~ x*xx,
## !           data = list(y = y, x = x, xx = xx), 
## !           control.compute = list(config = TRUE))
## ! 
## ! sam <- inla.posterior.sample(100, r)
## ! sam.extract <- inla.posterior.sample.eval(
## !     (function(...) {
## !         return(c(Intercept, x, xx, get("x:xx")))
## !     }), sam)
## ! print(round(dig = 4, rowMeans(sam.extract)))
## !
## ! sam.extract <- inla.posterior.sample.eval(c("Intercept", "x", "xx", "x:xx"), sam)
## ! print(round(dig = 4, rowMeans(sam.extract)))
## ! }
 
 
## Comments to the code below (contributed by CC):
##
## In order to call the right interpolating functions related to the different skewness mapping
## contained in the object "fun.splines" below we use a small trick for defining the right index order
## considering ties as well ( which is common in our skewness context). The following example should
## be clear enough:
##
## > x = sample(round(runif(3), digits = 2), size = 10, replace = T)
##
## > x
## [1] 0.88 0.88 0.16 0.16 0.16 0.53 0.53 0.53 0.16 0.53
##
## > as.numeric(factor(rank(a)))
## [1] 3 2 3 3 2 3 2 1 1 2
 
inla.create.sn.cache <- function() {
    ## Faster function for qsn,psn and dsn functions (using interpolations from the true ones)
    speed.fsn <- function(s, x, skew, fun.splines, deriv = 0) {
        ind <- order(skew)
        sss <- skew[ind]
        r <- as.numeric(factor(rank(sss)))
        skew.unique <- unique(sss)
        ## index for stored interpolated skewness by matching
        skew.ind <- findInterval(skew.unique, s)
        if (is.matrix(x)) {
            x.sample <- x[ind, , drop = FALSE]
            res.fsn <- unlist(lapply(
                seq_along(skew.unique),
                function(i) {
                lapply(t(x.sample[which(r == i), ]), fun.splines[[skew.ind[i]]])
            }
            ))
            res.fsn <- matrix(res.fsn, nrow = length(skew), ncol = ncol(x.sample), byrow = T)
            res.fsn <- res.fsn[order(ind), , drop = FALSE]
        } else {
            x.sample <- x[ind]
            res.fsn <- unlist(lapply(seq_along(skew.unique), function(i) {
                lapply(x.sample[which(r == i)], fun.splines[[skew.ind[i]]], deriv = deriv)
            }))
            res.fsn <- res.fsn[order(ind)]
        }
        return(res.fsn)
    }

    ## This code section will run only once in order to create the local object 'sn.cache'
    envir <- inla.get.inlaEnv()
    inla.require("sn", stop.on.error = TRUE)

    ## The following quantities can be changed and affect the whole code
    dig <- 2 # skewness precision (better keep this one for fast computations)
    sub <- 0.99 # highest value for the mapping function
    slb <- -sub # lowest value for the mapping function
    step <- 0.01 # step for the skewness sequence to be computed
    s <- seq(slb, sub, by = step) # overall skewness to compute
    s.pos <- seq(0, sub, by = step) # partial skewness to compute
    points <- seq(-4, 4, len = 50) # 50 points is a good deal for speed/accuracy

    ## Needed functions to store the tabulated values for 'qsn','psn' and 'dsn' functions
    ## Skewness Mapping (function)
    sn.map <- function(s, mu, sigma) {
        ## delta <- alpha/sqrt(1+alpha^2)
        s[which(is.na(s))] <- 0
        delta <- sign(s) * sqrt((pi / 2) * (abs(s)^(2 / 3) / (((4 - pi) / 2)^(2 / 3) + abs(s)^(2 / 3))))
        alpha <- delta / sqrt(1 - delta^2)
        xi <- mu - delta * sqrt((2 * sigma^2) / (pi - 2 * delta^2))
        omega <- sqrt((pi * sigma^2) / (pi - 2 * delta^2))
        return(list(alpha = alpha, xi = xi, omega = omega))
    }

    ## Needed function to construct the interpolations
    skew.spline.al <- function(s, skew, skew.store, type = "qsn") {
        res <- 0
        f <- function(x) {
            xx <- splinefun(s, x)
            return(xx)
        }
        s.fun <- apply(skew.store, 2, f)
        ## values table from "skew.points.table"
        if (type == "qsn") {
            res <- sapply(s.fun, function(f) f(abs(skew)))
            if (is.vector(res)) res <- as.matrix(t(res))
            ## Negative skewness positions computed by anti-symmetry relation
            if (any(skew < 0)) {
                neg.pos <- which(skew < 0)
                ord <- order(-res[1, ])
                res[neg.pos, ] <- -1 * res[neg.pos, ord]
            }
        } else {
            res <- sapply(s.fun, function(f) f(skew))
            if (is.vector(res)) res <- as.matrix(t(res))
        }
        return(res)
    }

    skew.points.table <- function(s.pos, s, points) {
        skew.store.qsn <- matrix(NA, nrow = length(s.pos), ncol = length(points))
        skew.store.jac <- matrix(NA, nrow = length(s), ncol = length(points))
        s.ind <- sn.map(s = s, mu = rep(0, length(s)), sigma = rep(1, length(s)))
        counter <- 1
        for (i in seq_along(s)) {
            alpha.x <- s.ind$alpha[i]
            xi.x <- s.ind$xi[i]
            omega.x <- s.ind$omega[i]
            skew.store.jac[i, ] <- sn::psn(points,
                                           xi = xi.x,
                                           omega = omega.x, alpha = alpha.x
                                           )
            if (s[i] >= 0) {
                skew.store.qsn[counter, ] <- sn::qsn(pnorm(points),
                                                     xi = xi.x,
                                                     omega = omega.x, alpha = alpha.x
                                                     )
                counter <- counter + 1
            }
        }
        return(list(skew.store.qsn = skew.store.qsn, skew.store.jac = skew.store.jac))
    }

    ## overall tabulated values for qsn,psn,dsn functions
    store.list <- skew.points.table(s.pos = s.pos, s = s, points = points)

    ## tabulated values for 'qsn' function
    store.qsn <- store.list$skew.store.qsn

    ## tabulated values for 'psn','dsn' functions
    store.jac <- store.list$skew.store.jac
    f.spline <- function(ss, type) {
        if (type != "jac") {
            val <- splinefun(points, skew.spline.al(s = s.pos, skew = ss, skew.store = store.qsn))
        } else {
            val <- splinefun(points, skew.spline.al(s = s, skew = ss, skew.store = store.jac, type = "jac"))
        }
        return(val)
    }
    ## matrix version of the skewness vector
    skew.ind <- as.matrix(s)
    ## splines for 'qsn'
    skew.splines <- apply(skew.ind, 1, f.spline, type = "fun")
    ## splines for 'psn','dsn'
    skew.splines.jac <- apply(skew.ind, 1, f.spline, type = "jac")
    ## Storing permanently the needed objects in a local INLA environment
    envir$sn.cache <- list(
        skew.splines = skew.splines, skew.splines.jac = skew.splines.jac,
        s = s, dig = dig, speed.fsn = speed.fsn
    )
    return(invisible())
}


inla.posterior.sample <- function(n = 1L, result, selection = list(),
                                  intern = FALSE,
                                  use.improved.mean = TRUE, skew.corr = TRUE,
                                  add.names = TRUE, seed = 0L, num.threads = NULL,
                                  parallel.configs = TRUE, verbose = FALSE) {
    ## New inla.posterior.sample with skewness correction. contributed by CC.

    stopifnot(!missing(result) && inherits(result, "inla"))
    if (is.null(result$misc$configs)) {
        stop("You need an inla-object computed with option 
         'control.compute=list(config = TRUE)'.")
    }

    if (!is.null(result$misc$configs$.preopt) && result$misc$configs$.preopt) {
        rfake <- list(mlik = result$mlik, 
                      summary.hyperpar = result$summary.hyperpar,
                      internal.summary.hyperpar = result$internal.summary.hyperpar, 
                      misc = list(from.theta = result$misc$from.theta,
                                  to.theta = result$misc$to.theta,
                                  configs = result$misc$configs))
        rfake$misc$configs$.preopt <- FALSE
        class(rfake) <- "inla"

        ct <- result$misc$configs$contents
        for(nm in c("APredictor", "Predictor")) {
            if (ct$tag[1] == nm) {
                ct$tag <- ct$tag[-1]
                ct$start <- ct$start[-1] - ct$start[2] + 1
                ct$length <- ct$length[-1]
            }
        }
        rfake$misc$configs$contents <- ct

        ## yes, selection is not an argument in this call, but we account for it afterwards
        xx <- inla.posterior.sample(n, rfake, intern = intern,
                                    use.improved.mean = use.improved.mean,
                                    skew.corr = skew.corr,
                                    add.names = add.names, seed = seed,
                                    num.threads = num.threads,
                                    parallel.configs = parallel.configs, verbose = verbose)

        A <- rbind(result$misc$configs$pA %*% result$misc$configs$A,
                   result$misc$configs$A)
        sel <- inla.posterior.sample.interpret.selection(selection, result)

        if (nrow(result$misc$configs$pA) > 0) {
            pnam <- c(paste0("APredictor:", seq_len(nrow(result$misc$configs$pA))),
                      paste0("Predictor:", seq_len(nrow(result$misc$configs$A))))
        } else {
            pnam <- paste0("Predictor:", seq_len(nrow(result$misc$configs$A)))
        }

        cont.A <- NULL
        cont.pA <- NULL
        if (!is.null(result$misc$configs$A)) {
            cont.A <- result$misc$configs$A
        }
        if (!is.null(result$misc$configs$pA)) {
            cont.pA <- result$misc$configs$pA
        }
        attr(xx, ".contents") <- c(result$misc$configs$contents, A = cont.A, pA = cont.pA)

        off <- result$misc$configs$offsets
        for(i in seq_along(xx)) {
            nam <- rownames(xx[[i]]$latent)
            xx[[i]]$latent <- c(as.numeric(A %*% xx[[i]]$latent) + off, as.numeric(xx[[i]]$latent))[sel]
            if (!is.null(nam)) {
                names(xx[[i]]$latent) <- c(pnam, nam)[sel]
            }
            xx[[i]]$latent <- as.matrix(xx[[i]]$latent, ncol = 1)
        }
        return(xx)
    }

    if (is.null(num.threads)) {
        num.threads <- inla.getOption("num.threads")
    }
    if (seed != 0L) {
        num.threads.user <- inla.parse.num.threads(num.threads)
        num.threads <- inla.parse.num.threads("1:1")
        if (num.threads != num.threads.user) {
            warning("Since 'seed!=0', parallel model is disabled and serial model is selected, num.threads='1:1'")
        }
        ## Since parallel.configs aren't used on winows anyway, warning about
        ## overriding it on other systems is overkill and doesn't really help.
        # if (!missing(parallel.configs) && !isFALSE(parallel.configs)) {
        # warning("Since 'seed!=0', parallel model is disabled and serial model is selected, parallel.configs=FALSE")
        # }
        parallel.configs <- FALSE
    } else {
        num.threads <- inla.parse.num.threads(num.threads)
    }

    if (use.improved.mean == FALSE) {
        skew.corr <- FALSE
    }
    if (skew.corr == TRUE && !exists("sn.cache", envir = inla.get.inlaEnv())) {
        ## function call for creating a local object in INLA environment
        inla.create.sn.cache()
    }

    sel <- inla.posterior.sample.interpret.selection(selection, result)
    if (sum(sel) == 0) {
        return(matrix(NA, 0, 0))
    }

    sel.n <- sum(sel)
    sel.map <- which(sel)
    stopifnot(length(sel.map) == sel.n)
    names(sel.map) <- names(sel[sel.map])

    n <- as.integer(n)
    stopifnot(is.integer(n) && n > 0L)

    ## first sample the theta(-index)
    cs <- result$misc$configs
    ld <- numeric(cs$nconfig)
    for (i in 1:cs$nconfig) {
        ld[i] <- cs$config[[i]]$log.posterior
    }
    p <- exp(ld - max(ld))
    idx <- sort(sample(1:cs$nconfig, n, prob = p, replace = TRUE))
    n.idx <- numeric(cs$nconfig)
    n.idx[] <- 0
    for (i in 1:cs$nconfig) {
        n.idx[i] <- sum(idx == i)
    }
    con <- cs$contents
    all.samples <- rep(list(c()), n)
    i.sample <- 1L

    ## it is easier to parallise here,  just doing the samples in parallel. then we process them
    ## in serial afterwards, since that is much less of an job

    if (!inla.os("windows") && parallel.configs) {
        nt <- as.numeric(strsplit(num.threads, ":")[[1]])
        ncores <- parallel::detectCores(all.tests = TRUE, logical = FALSE)
        if (nt[1] == 0) {
            nt[1] <- max(1, min(cs$nconfig, ncores))
            nt[2] <- 1
        }
        if (nt[2] == 0) {
            nt[2] <- max(1, ncores %/% nt[1])
        }
        # In parallel mode, each thread needs it's own random sequence,
        # so seed must be 0; seed != 0 is handled above to prevent parallel
        # runs with nonzero seed.
        xx.list <- parallel::mclapply(1:cs$nconfig,
            (function(k) {
                if (n.idx[k] > 0) {
                    xx <- inla.qsample(
                        n = n.idx[k],
                        Q = cs$config[[k]]$Q,
                        mu = inla.ifelse(
                            use.improved.mean,
                            cs$config[[k]]$improved.mean,
                            cs$config[[k]]$mean
                        ),
                        constr = cs$constr,
                        logdens = TRUE,
                        seed = 0L,
                        num.threads = paste0(nt[2], ":1"),
                        selection = sel.map,
                        verbose = verbose
                    )
                    return(xx)
                } else {
                    return(NA)
                }
            }),
            mc.cores = nt[1],
            mc.preschedule = TRUE
        )
    } else {
        xx.list <- NULL
    }

    for (k in 1:cs$nconfig) {
        if (n.idx[k] > 0) {
            ## then the latent field
            if (is.null(xx.list)) {
                xx <- inla.qsample(
                    n = n.idx[k], Q = cs$config[[k]]$Q,
                    mu = inla.ifelse(use.improved.mean, cs$config[[k]]$improved.mean, cs$config[[k]]$mean),
                    constr = cs$constr, logdens = TRUE, seed = seed, num.threads = num.threads,
                    selection = sel.map, verbose = verbose
                )
            }
            else {
                xx <- xx.list[[k]]
            }
            ## if user set seed,  then just continue this rng-stream
            if (seed > 0L) seed <- -1L

            ld.theta <- cs$max.log.posterior + cs$config[[k]]$log.posterior
            nm <- names(sel.map)

            theta <- cs$config[[k]]$theta
            log.J <- 0.0
            if (!is.null(theta) && !intern) {
                for (j in 1:length(theta)) {
                    theta[j] <- do.call(result$misc$from.theta[[j]],
                        args = list(theta[j])
                    )
                }
                names(theta) <- inla.transform.names(result, names(theta))

                if (TRUE) {
                    ## new fancy code using the automatic differentiation feature in R
                    for (i in 1:length(theta)) {
                        arg.val <- formals(result$misc$from.theta[[i]])
                        arg <- names(arg.val)
                        if (length(arg) == 1L) {
                            deriv.func <- inla.eval(paste("function(", arg, ") {}"))
                        } else {
                            if (length(arg) == 2L) {
                                deriv.func <- inla.eval(paste(
                                    "function(", arg[1L], ",",
                                    arg[2L], "=", arg.val[2L], ") {}"
                                ))
                            }
                            else {
                                stopifnot(length(arg) == 3L)
                                deriv.func <- inla.eval(paste(
                                    "function(", arg[1L], ",",
                                    arg[2L], "=", arg.val[2L], ",",
                                    arg[3L], "=", arg.val[3L], ") {}"
                                ))
                            }
                        }
                        temptest <- try(D(body(result$misc$from.theta[[i]]), arg[1L]),
                            silent = TRUE
                        )
                        if (!inherits(temptest, "try-error")) {
                            body(deriv.func) <- temptest
                            log.J <- log.J - log(abs(deriv.func(cs$config[[k]]$theta[i])))
                            ## Yes, it's a minus...
                        }
                        else {
                            h <- .Machine$double.eps^0.25
                            theta.1 <- do.call(result$misc$from.theta[[i]],
                                args = list(cs$config[[k]]$theta[i] - h)
                            )
                            theta.2 <- do.call(result$misc$from.theta[[i]],
                                args = list(cs$config[[k]]$theta[i] + h)
                            )
                            log.J <- log.J - log(abs((theta.2 - theta.1) / (2.0 * h)))
                            ## Yes, it's a minus...
                        }
                    }
                    ## print(paste("logJ", log.J))
                } else {
                    ## old code using numerical differentiation
                    h <- .Machine$double.eps^0.25
                    for (i in 1:length(theta)) {
                        theta.1 <- do.call(result$misc$from.theta[[i]],
                            args = list(cs$config[[k]]$theta[i] - h)
                        )
                        theta.2 <- do.call(result$misc$from.theta[[i]],
                            args = list(cs$config[[k]]$theta[i] + h)
                        )
                        log.J <- log.J - log(abs((theta.2 - theta.1) / (2.0 * h)))
                        ## Yes, it's a minus...
                    }
                    ## print(paste("logJ", log.J))
                }
            }
            if (!skew.corr) {
                ## Default for mean correction only
                sample.corr <- xx$sample
                C.th <- xx$logdens
            } else {
                ## Skewness Correction applied on the samples
                sh.cache <- get("sn.cache", envir = inla.get.inlaEnv())
                dig <- sh.cache$dig
                s <- sh.cache$s
                skew.splines <- sh.cache$skew.splines
                skew.splines.jac <- sh.cache$skew.splines.jac
                speed.fsn <- sh.cache$speed.fsn
                mean.GA <- cs$config[[k]]$mean
                mean.SN <- cs$config[[k]]$improved.mean
                sigma.th <- sqrt(diag(cs$config[[k]]$Qinv))
                val.max <- (mean.GA - mean.SN) / sigma.th
                skew.val <- cs$config[[k]]$skewness
                ## if the skewness is NA then we simply put 0
                skew.val[which(is.na(skew.val))] <- 0
                ## we round skewness values to the second digit for matching the correct interpolations
                skew.val <- round(skew.val, digits = dig)
                s.new <- round(s, digits = dig)
                skew.table.max <- max(s.new)
                if (any(abs(skew.val) > skew.table.max)) {
                    skew.val <- pmax(-skew.table.max, pmin(skew.table.max, skew.val))
                    warning(paste0("One or more marginal skewness are too high. Coerced to be ", skew.table.max))
                }
                ## Positions for skewness different from 0
                zero.not <- which(skew.val != 0)
                if (length(zero.not) == 0 || !any(sel.map %in% zero.not)) {
                    sample.corr <- xx$sample
                    C.th <- xx$logdens
                } else {
                    skew.val.not <- skew.val[zero.not]
                    val.max.not <- val.max[zero.not]
                    ## 'dsn' evaluations by interpolation
                    dsn.not <- speed.fsn(s = s.new, x = val.max.not, skew = skew.val.not, fun.splines = skew.splines.jac, deriv = 1)
                    ## 'psn' evaluations by interpolation
                    psn.not <- speed.fsn(s = s.new, x = val.max.not, skew = skew.val.not, fun.splines = skew.splines.jac)
                    dsn.new <- dnorm(val.max)
                    dsn.new[zero.not] <- dsn.not
                    psn.new <- pnorm(val.max)
                    psn.new[zero.not] <- psn.not
                    if (length(sel.map) < length(skew.val)) {
                        mean.SN <- mean.SN[sel.map]
                        sigma.th <- sigma.th[sel.map]
                        skew.val <- skew.val[sel.map]
                        zero.not <- which(skew.val != 0)
                        if (length(zero.not) != 0) {
                            skew.val.not <- skew.val[zero.not]
                        }
                    }
                    x.sample.not <- xx$sample[zero.not, , drop = FALSE]
                    mean.SN.not <- mean.SN[zero.not]
                    sigma.th.not <- sigma.th[zero.not]
                    x.val.not <- (x.sample.not - mean.SN.not) / sigma.th.not
                    ## faster vectorized 'qsn' function call
                    fast.qsn <- speed.fsn(s = s.new, x = x.val.not, skew = skew.val.not, fun.splines = skew.splines)
                    sample.corr <- xx$sample
                    sample.store <- sigma.th.not * fast.qsn + mean.SN.not
                    sample.corr[zero.not, ] <- sample.store[1:nrow(sample.store), ]
                    C.th <- xx$logdens + sum(log(dsn.new)) - sum(log(dnorm(qnorm(psn.new))))
                }
            }
            for (i in 1:n.idx[k]) {
                if (is.null(theta)) {
                    a.sample <- list(
                        hyperpar = NULL,
                        latent = sample.corr[, i, drop = FALSE],
                        logdens = list(
                            hyperpar = NULL,
                            latent = as.numeric(C.th[i]),
                            joint = as.numeric(C.th[i])
                        )
                    )
                } else {
                    ld.h <- as.numeric(ld.theta - result$mlik[1, 1] + log.J)
                    a.sample <- list(
                        hyperpar = theta,
                        latent = sample.corr[, i, drop = FALSE],
                        logdens = list(
                            hyperpar = ld.h,
                            latent = as.numeric(C.th[i]),
                            joint = as.numeric(ld.h + C.th[i])
                        )
                    )
                }
                if (add.names || i.sample == 1L) {
                    n1 <- length(nm)
                    n2 <- length(a.sample$latent)
                    stopifnot(n2 >= n1) ## this must be true. just a check
                    if (n2 > n1) {
                        ## This is the case where lincomb.derived.only = FALSE, so these are
                        ## then added to the end. Should transfer the names of them all the way
                        ## to here, but...
                        xnm <- paste0("Lincomb:", inla.num(1:(n2 - n1)))
                        nm <- c(nm, xnm)

                        if (i.sample == 1L) {
                            ## add to the contents if needed
                            con$tag <- c(con$tag, "Lincomb")
                            con$start <- c(con$start, n1 + 1)
                            con$length <- c(con$length, n2 - n1)
                        }
                    }
                    rownames(a.sample$latent) <- nm
                } else {
                    rownames(a.sample$latent) <- NULL
                }
                all.samples[[i.sample]] <- a.sample
                i.sample <- i.sample + 1L
            }
        }
    }
    if (length(selection) == 0L) {
        attr(all.samples, ".contents") <- con
    } else {
        ## we use a selection, need to build a new 'contents' list
        con <- list(tag = names(selection), start = c(), length = c())
        for (nm in names(selection)) {
            ## from Hmisc::escapeRegex. Need it for the '(Intercept)'
            re <- paste0("^", gsub("([.|()\\^{}+$*?]|\\[|\\])", "\\\\\\1", nm), ":")
            m <- grep(re, names(sel.map))
            if (length(m) > 0) {
                con$start <- c(con$start, min(m))
                con$length <- c(con$length, diff(range(m)) + 1)
            }
        }
        attr(all.samples, ".contents") <- con
    }
    return(all.samples)
}


`inla.posterior.sample.orig` <- function(n = 1, result, selection = list(), intern = FALSE,
                                         use.improved.mean = TRUE, add.names = TRUE, seed = 0L,
                                         num.threads = NULL, verbose = FALSE) {
    ## this is the original version, before the skewness correction. keep it here for completeness.

    stopifnot(!missing(result) && inherits(result, "inla"))
    if (is.null(result$misc$configs)) {
        stop("You need an inla-object computed with option 'control.compute=list(config = TRUE)'.")
    }

    if (seed != 0L && is.null(num.threads)) {
        num.threads <- 1L
    }
    if (is.null(num.threads)) {
        num.threads <- inla.getOption("num.threads")
    }
    num.threads <- max(num.threads, 1L)
    if (num.threads > 1L && seed != 0L) {
        stop("num.threads > 1L require seed = 0L")
    }

    sel <- inla.posterior.sample.interpret.selection(selection, result)
    if (sum(sel) == 0) {
        return(matrix(NA, 0, 0))
    }

    sel.n <- sum(sel)
    sel.map <- which(sel)
    stopifnot(length(sel.map) == sel.n)
    names(sel.map) <- names(sel[sel.map])

    n <- as.integer(n)
    stopifnot(is.integer(n) && n > 0L)

    ## first sample the theta(-index)
    cs <- result$misc$configs
    ld <- numeric(cs$nconfig)
    for (i in 1:cs$nconfig) {
        ld[i] <- cs$config[[i]]$log.posterior
    }
    p <- exp(ld - max(ld))
    idx <- sort(sample(1:cs$nconfig, n, prob = p, replace = TRUE))
    n.idx <- numeric(cs$nconfig)
    n.idx[] <- 0
    for (i in 1:cs$nconfig) {
        n.idx[i] <- sum(idx == i)
    }

    con <- cs$contents
    all.samples <- rep(list(c()), n)
    i.sample <- 1L
    for (k in 1:cs$nconfig) {
        if (n.idx[k] > 0) {
            ## then the latent field
            xx <- inla.qsample(
                n = n.idx[k], Q = cs$config[[k]]$Q,
                mu = inla.ifelse(use.improved.mean, cs$config[[k]]$improved.mean, cs$config[[k]]$mean),
                constr = cs$constr, logdens = TRUE, seed = seed, num.threads = num.threads,
                selection = sel.map, verbose = verbose
            )
            ## if user set seed,  then just continue this rng-stream
            if (seed > 0L) seed <- -1L

            ld.theta <- cs$max.log.posterior + cs$config[[k]]$log.posterior
            nm <- names(sel.map)

            theta <- cs$config[[k]]$theta
            log.J <- 0.0
            if (!is.null(theta) && !intern) {
                for (j in 1:length(theta)) {
                    theta[j] <- do.call(result$misc$from.theta[[j]], args = list(theta[j]))
                }
                names(theta) <- inla.transform.names(result, names(theta))

                if (TRUE) {
                    ## new fancy code using the automatic differentiation feature in R
                    for (i in 1:length(theta)) {
                        arg.val <- formals(result$misc$from.theta[[i]])
                        arg <- names(arg.val)
                        if (length(arg) == 1L) {
                            deriv.func <- inla.eval(paste("function(", arg, ") {}"))
                        } else {
                            if (length(arg) == 2L) {
                                deriv.func <- inla.eval(paste("function(", arg[1L], ",", arg[2L], "=", arg.val[2L], ") {}"))
                            }
                            else {
                                stopifnot(length(arg) == 3L)
                                deriv.func <- inla.eval(paste("function(", arg[1L], ",", arg[2L], "=", arg.val[2L], ",", arg[3L], "=", arg.val[3L], ") {}"))
                            }
                        }
                        temptest <- try(D(body(result$misc$from.theta[[i]]), arg[1L]), silent = TRUE)
                        if (!inherits(temptest, "try-error")) {
                            body(deriv.func) <- temptest
                            log.J <- log.J - log(abs(deriv.func(cs$config[[k]]$theta[i]))) ## Yes, it's a minus...
                        }
                        else {
                            h <- .Machine$double.eps^0.25
                            theta.1 <- do.call(result$misc$from.theta[[i]], args = list(cs$config[[k]]$theta[i] - h))
                            theta.2 <- do.call(result$misc$from.theta[[i]], args = list(cs$config[[k]]$theta[i] + h))
                            log.J <- log.J - log(abs((theta.2 - theta.1) / (2.0 * h))) ## Yes, it's a minus...
                        }
                    }
                    ## print(paste("logJ", log.J))
                } else {
                    ## old code using numerical differentiation
                    h <- .Machine$double.eps^0.25
                    for (i in 1:length(theta)) {
                        theta.1 <- do.call(result$misc$from.theta[[i]], args = list(cs$config[[k]]$theta[i] - h))
                        theta.2 <- do.call(result$misc$from.theta[[i]], args = list(cs$config[[k]]$theta[i] + h))
                        log.J <- log.J - log(abs((theta.2 - theta.1) / (2.0 * h))) ## Yes, it's a minus...
                    }
                    ## print(paste("logJ", log.J))
                }
            }

            for (i in 1:n.idx[k]) {
                if (is.null(theta)) {
                    a.sample <- list(
                        hyperpar = NULL,
                        latent = xx$sample[, i, drop = FALSE],
                        logdens = list(
                            hyperpar = NULL,
                            latent = as.numeric(xx$logdens[i]),
                            joint = as.numeric(xx$logdens[i])
                        )
                    )
                } else {
                    ld.h <- as.numeric(ld.theta - result$mlik[1, 1] + log.J)
                    a.sample <- list(
                        hyperpar = theta,
                        latent = xx$sample[, i, drop = FALSE],
                        logdens = list(
                            hyperpar = ld.h,
                            latent = as.numeric(xx$logdens[i]),
                            joint = as.numeric(ld.h + xx$logdens[i])
                        )
                    )
                }
                if (add.names || i.sample == 1L) {
                    n1 <- length(nm)
                    n2 <- length(a.sample$latent)
                    stopifnot(n2 >= n1) ## this must be true. just a check
                    if (n2 > n1) {
                        ## This is the case where lincomb.derived.only = FALSE, so these are
                        ## then added to the end. Should transfer the names of them all the way
                        ## to here, but...
                        xnm <- paste0("Lincomb:", inla.num(1:(n2 - n1)))
                        nm <- c(nm, xnm)

                        if (i.sample == 1L) {
                            ## add to the contents if needed
                            con$tag <- c(con$tag, "Lincomb")
                            con$start <- c(con$start, n1 + 1)
                            con$length <- c(con$length, n2 - n1)
                        }
                    }
                    rownames(a.sample$latent) <- nm
                } else {
                    rownames(a.sample$latent) <- NULL
                }
                all.samples[[i.sample]] <- a.sample
                i.sample <- i.sample + 1L
            }
        }
    }
    if (length(selection) == 0L) {
        attr(all.samples, ".contents") <- con
    } else {
        ## we use a selection, need to build a new 'contents' list
        con <- list(tag = names(selection), start = c(), length = c())
        for (nm in names(selection)) {
            ## from Hmisc::escapeRegex. Need it for the '(Intercept)'
            re <- paste0("^", gsub("([.|()\\^{}+$*?]|\\[|\\])", "\\\\\\1", nm), ":")
            m <- grep(re, names(sel.map))
            if (length(m) > 0) {
                con$start <- c(con$start, min(m))
                con$length <- c(con$length, diff(range(m)) + 1)
            }
        }
        attr(all.samples, ".contents") <- con
    }
    return(all.samples)
}

`inla.posterior.sample.eval` <- function(fun, samples, return.matrix = TRUE, ...) {
    ## evaluate fun(...) over each sample

    contents <- attr(samples, which = ".contents", exact = TRUE)
    if (is.null(contents)) {
        stop("Argument 'samples' must be the output from 'inla.posterior.sample()'.")
    }

    ## special shorthand feature that is very useful. if this is vector of character's, then its
    ## interpreted as a function returning these names, like ....eval(c("a", "b[1:2]"), ...)
    ## will be converted into the function: function() return(c(get('a'), get('b[1:2]')))

    if (is.character(fun)) {
        arg <- NULL
        for(nm in fun) {
            arg <- paste0(if (is.null(arg)) "" else paste0(arg, ","), "get('", nm, "')")
        }
        fun <- inla.eval(paste0("function(...) return(c(", arg, "))"))
    } else {
        fun <- match.fun(fun)
    }

    my.fun <- function(a.sample, .contents, .fun, ...) {
        env <- new.env()
        theta <- as.vector(a.sample$hyperpar)
        assign("theta", theta, envir = env)
        if (!is.null(.contents$A)) {
            assign("A", .contents$A, envir = env)
        }
        if (!is.null(.contents$pA)) {
            assign("pA", .contents$pA, envir = env)
        }
        for (i in seq_along(.contents$tag)) {
            assign(.contents$tag[i],
                as.vector(a.sample$latent[.contents$start[i]:(.contents$start[i] +
                    .contents$length[i] - 1), 1]),
                envir = env
            )
        }
        ## this is special
        if (exists("(Intercept)", envir = env)) {
            assign("Intercept", get("(Intercept)", envir = env), envir = env)
        }

        parent.env(env) <- .GlobalEnv
        environment(.fun) <- env
        return(.fun(...))
    }

    ret <- inla.mclapply(samples, my.fun, .fun = fun, .contents = contents, ...)
    if (return.matrix) {
        ns <- length(ret)
        nm <- names(ret[[1]])
        ret <- matrix(unlist(ret), ncol = ns)
        colnames(ret) <- paste0("sample:", 1:ns)
        rownames(ret) <- if (!is.null(nm)) nm else paste0("fun[", 1:nrow(ret), "]")
    }

    return(ret)
}

`inla.posterior.sample.interpret.selection` <- function(selection = list(), result) {
    ## this function interpret a selection, of the form of a named list,
    ## list(NAME = idx's, ...),
    ## with the standard names 'APredictor', 'Predictor', '(Intercept)' as well. the idx can
    ## contains negative numbers for which will be interpreted as 'not'. if idx=0, then this is
    ## interpreted as the whole vector. the result is a named list of vector of logicals, that
    ## described which part of the sample to select

    cs <- result$misc$configs$contents
    nc <- length(cs$tag)
    n <- sum(cs$length)
    select <- rep(FALSE, n)
    nam <- rep(NA, n)

    ## is selection is NULL or an empty list, this means select all. just make a selection that
    ## do that
    if (is.null(selection) || length(selection) == 0) {
        selection <- as.list(rep(0, nc))
        names(selection) <- cs$tag
    }

    for (k in seq_along(cs$tag)) {
        tag <- cs$tag[k]
        start <- cs$start[k]
        end <- start + cs$length[k] - 1L
        len <- cs$length[k]

        if (inla.is.element(tag, selection)) {
            idx <- which(names(selection) == tag)
            sel <- selection[[idx]]
            selection[[idx]] <- NULL

            if (all(sel == 0)) {
                sel <- 1:len
            } else if (all(sel > 0)) {
                stopifnot(all(sel <= len))
            } else if (all(sel < 0)) {
                stopifnot(all(-sel <= len))
                sel <- (1:len)[sel]
            } else {
                stop(paste("This should not happen. Something wrong with the selection for tag=", tag))
            }

            select[start + sel - 1] <- TRUE
            nam[start + sel - 1] <- paste0(tag, ":", sel)
        }
    }

    if (length(selection) > 0) {
        stop(paste0(
            "Some selections are not used: ",
            paste(names(selection), collapse = ", ", sep = "")
        ))
    }
    names(select) <- nam

    return(select)
}

`inla.selection2lincombs` <- function(selection = list(), tag.base = ".selection.") {
    ## this function convert from selection format to lincombs format.
    lc <- list()
    n.lc <- 0
    nms <- names(selection)
    nnms <- c()
    for (i in seq_along(selection)) {
        a <- selection[[i]]
        for (j in seq_along(a)) {
            n.lc <- n.lc + 1
            llc <- list()
            llc[[1]] <- list(list(idx = a[j], weight = 1))
            names(llc[[1]]) <- nms[i]
            lc[[n.lc]] <- llc
            nnms <- c(nnms, paste0(tag.base, nms[i], ":", a[j]))
        }
    }
    names(lc) <- nnms
    return(lc)
}
