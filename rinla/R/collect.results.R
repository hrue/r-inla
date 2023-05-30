## Export: inla.collect.results

## ! \name{inla.collect.results}
## ! \alias{inla.collect.results}
## ! \alias{collect.results}
## ! \title{Collect results from a inla-call}
## ! \description{\code{inla.collect.results} collect results  from a inla-call}
## ! \usage{
## ! inla.collect.results(
## !              results.dir,
## !              debug=FALSE,
## !              only.hyperparam=FALSE,
## !              file.log = NULL,
## !              file.log2 = NULL)
## !}
## ! \arguments{

`inla.collect.results` <-
    function(
             ## ! \item{results.dir}{The directory where the results of the inla run are stored}
             results.dir,

             ## ! \item{debug}{Logical. If \code{TRUE} some debugging information are printed}
             debug = FALSE,

             ## ! \item{only.hyperparam}{Binary variable indicating wheather only the
             ## ! results for the hyperparameters should be collected}
             only.hyperparam = FALSE,

             ## ! \item{file.log}{Character. The filename, if any, of the logfile for
             ## ! the internal calculations}
             file.log = NULL,

             ## ! \item{file.log2}{Character. The filename, if any, of the logfile2 for
             ## ! the internal calculations}
             file.log2 = NULL,

             ## ! \item{silent}{Internal use only}
             silent = inla.getOption("silent")
             )
{
    ## ! }
    ## ! \value{ The function returns an object of class \code{"inla"}, see the
    ## ! help file for \code{inla} for details.}
    ## !
    ## ! \details{This function is mainly used inside \code{inla}
    ## !  to collect results after running the inla
    ## !  function. It can also be used to collect results into R after having
    ## !  runned a inla section outside R.  }

    if (is.na(file.info(results.dir)$isdir) ||
        !file.info(results.dir)$isdir) {
        stop(paste("This is not a directory: ", results.dir, "\n"))
    }

    filename <- paste0(results.dir, "/.ok")
    res.ok <- file.exists(filename)
    if (!res.ok) {
        ## try this one instead
        results.dir.new <- paste(results.dir, "/results.files", sep = "")
        filename <- paste(results.dir.new, "/.ok", sep = "")
        res.ok <- file.exists(filename)
        if (res.ok) {
            if (debug) {
                cat(paste("inla.collect.results: retry with directory", results.dir.new, "\n"))
            }
            return(inla.collect.results(results.dir.new,
                                        debug = debug,
                                        only.hyperparam = only.hyperparam,
                                        file.log = file.log,
                                        file.log2 = file.log2
                                        ))
        } else {
            ## neither directories contain the file /.ok, then we
            ## assume the inla-program has crashed
            return (if (silent != 2L) inla.inlaprogram.has.crashed() else NULL)
        }
    }

    filename <- paste0(results.dir, "/dryrun")
    if (file.exists(filename)) {
        res <- list(dryrun = readLines(filename))
        class(res) <- "inla.dryrun"
        return(res)
    }

    is.null.list <- function(alist) {
        return (all(sapply(alist, is.null)))
    }

    if (!only.hyperparam) {
        res.fixed <- inla.collect.fixed(results.dir, debug)
        res.lincomb <- inla.collect.lincomb(results.dir, debug, derived = FALSE)
        res.lincomb.derived <- inla.collect.lincomb(results.dir, debug, derived = TRUE)
        res.dic <- inla.collect.dic(results.dir, debug)
        res.cpo.pit <- inla.collect.cpo(results.dir, debug)
        res.gcpo <- inla.collect.gcpo(results.dir, debug)
        res.po <- inla.collect.po(results.dir, debug)
        if (!is.null.list(res.dic) && !is.null.list(res.po)) {
            res.deviance.residuals <- list(deviance.residuals = res.dic$residuals.deviance)
            res.dic$residuals.deviance <- NULL
        } else {
            res.deviance.residuals <- list(deviance.residuals = NULL)
        }
        res.waic <- inla.collect.waic(results.dir, debug)
        res.random <- inla.collect.random(results.dir, debug)
        res.predictor <- inla.collect.predictor(results.dir, debug)
        res.spde2.blc <- inla.collect.spde2.blc(results.dir, debug)
        res.spde3.blc <- inla.collect.spde3.blc(results.dir, debug)
    } else {
        res.fixed <- NULL
        res.lincomb <- NULL
        res.lincomb.derived <- NULL
        res.dic <- NULL
        res.cpo.pit <- NULL
        res.gcpo <- NULL
        res.po <- NULL
        res.deviance.residuals <- NULL
        res.waic <- NULL
        res.random <- NULL
        res.predictor <- NULL
        res.spde2.blc <- NULL
        res.spde3.blc <- NULL
    }
    res.mlik <- inla.collect.mlik(results.dir, debug)
    res.q <- inla.collect.q(results.dir, debug)
    res.graph <- inla.collect.graph(results.dir, debug)
    res.offset <- inla.collect.offset.linear.predictor(results.dir, debug)

    ## get the hyperparameters
    theta.mode <- inla.read.binary.file(paste(results.dir, .Platform$file.sep, ".theta_mode", sep = ""))[-1]
    x.mode <- inla.read.binary.file(paste(results.dir, .Platform$file.sep, ".x_mode", sep = ""))[-1]
    gitid <- readLines(paste(results.dir, .Platform$file.sep, ".gitid", sep = ""))

    lfn.fnm <- paste(results.dir, .Platform$file.sep, "linkfunctions.names", sep = "")
    if (file.exists(lfn.fnm)) {
        linkfunctions.names <- readLines(lfn.fnm)
        fp <- file(paste(results.dir, .Platform$file.sep, "linkfunctions.link", sep = ""), "rb")
        n <- readBin(fp, integer(), 1)
        idx <- readBin(fp, double(), n)
        ok <- which(!is.nan(idx))
        idx[ok] <- idx[ok] + 1
        close(fp)
        linkfunctions <- list(names = linkfunctions.names, link = as.integer(idx))
    } else {
        linkfunctions <- NULL
    }

    if (length(theta.mode) > 0) {
        res.hyper <- inla.collect.hyperpar(results.dir, debug)

        ## get the joint (if printed)
        alldir <- dir(results.dir)
        if (length(grep("joint.dat", alldir)) == 1) {
            if (debug) {
                print("inla.collect.joint hyperpar")
            }
            fnm <- paste(results.dir, "/joint.dat", sep = "")
            if (file.info(fnm)$size > 0) {
                joint.hyper <- read.table(fnm)
            } else {
                joint.hyper <- NULL
            }
        } else {
            joint.hyper <- NULL
        }
    } else {
        res.hyper <- NULL
        joint.hyper <- NULL
    }

    logfile <- list(logfile = c(
                        inla.collect.logfile(file.log, debug)$logfile,
                        "", paste(rep("*", 72), sep = "", collapse = ""), "",
                        inla.collect.logfile(file.log2, debug)$logfile
                    ))
    misc <- inla.collect.misc(results.dir, debug)
    theta.tags <- NULL
    mode.status <- NA
    if (!is.null(misc)) {
        ## put also theta.mode in here
        misc$theta.mode <- theta.mode
        ## we need theta.tags for later usage
        if (!is.null(misc$theta.tags)) {
            theta.tags <- misc$theta.tags
        }
        mode.status <- misc$mode.status

        if (!is.null(misc$lincomb.derived.correlation.matrix)) {
            if (!is.null(res.lincomb.derived)) {
                id <- res.lincomb.derived$summary.lincomb.derived$ID
                tag <- rownames(res.lincomb.derived$summary.lincomb.derived)
                R <- misc$lincomb.derived.correlation.matrix
                rownames(R) <- colnames(R) <- tag[id]
                misc$lincomb.derived.correlation.matrix <- R
            } else {
                misc$lincomb.derived.correlation.matrix <- NULL
            }
        }
        if (!is.null(misc$lincomb.derived.covariance.matrix)) {
            if (!is.null(res.lincomb.derived)) {
                id <- res.lincomb.derived$summary.lincomb.derived$ID
                tag <- rownames(res.lincomb.derived$summary.lincomb.derived)
                R <- misc$lincomb.derived.covariance.matrix
                rownames(R) <- colnames(R) <- tag[id]
                misc$lincomb.derived.covariance.matrix <- R
            } else {
                misc$lincomb.derived.covariance.matrix <- NULL
            }
        }

        ## also put the linkfunctions here
        misc$linkfunctions <- linkfunctions
        if (!is.null(linkfunctions)) {
            ## a better name
            misc$family <- linkfunctions$link
        }
    }

    ## add the names of the theta's here, as they are available.
    if (!is.null(misc) && !is.null(joint.hyper)) {
        colnames(joint.hyper) <- c(misc$theta.tags,
                                   "Log posterior density",
                                   "Total integration weight (log.dens included)")
    }

    names(theta.mode) <- theta.tags
    res <- c(res.fixed, res.lincomb, res.lincomb.derived, res.mlik,
             list(cpo = res.cpo.pit), list(gcpo = res.gcpo), list(po = res.po), list(waic = res.waic),
             list(residuals = res.deviance.residuals), 
             res.random, res.predictor, res.hyper,
             res.offset, res.spde2.blc, res.spde3.blc, logfile,
             list(
                 misc = misc,
                 dic = res.dic, mode = list(
                                    theta = theta.mode, x = x.mode,
                                    theta.tags = theta.tags, mode.status = mode.status,
                                    log.posterior.mode = misc$log.posterior.mode
                                ),
                 joint.hyper = joint.hyper, nhyper = length(theta.mode),
                 version = list(inla.call = gitid, R.INLA = inla.version("version"))
             ),
             list(Q = res.q),
             res.graph,
             ok = res.ok
             )
    class(res) <- "inla"
    return(res)
}

`inla.collect.misc` <- function(dir, debug = FALSE) 
{
    d <- paste(dir, "/misc", sep = "")
    d.info <- file.info(d)$isdir

    if (debug) {
        print(paste("collect misc from", d))
    }

    if (is.na(d.info) || (d.info == FALSE)) {
        return(NULL)
    }

    fnm <- paste(d, "/theta-tags", sep = "")
    if (file.exists(fnm)) {
        tags <- readLines(fnm)
    } else {
        tags <- NULL
    }

    fnm <- paste(d, "/theta-from", sep = "")
    if (file.exists(fnm)) {
        theta.from <- readLines(fnm)
        ## evaluate these as functions
        theta.from <- lapply(theta.from, inla.source2function)
        if (!is.null(tags)) {
            names(theta.from) <- tags
        }
    } else {
        theta.from <- NULL
    }

    fnm <- paste(d, "/theta-to", sep = "")
    if (file.exists(fnm)) {
        theta.to <- readLines(fnm)
        ## evaluate these as functions
        theta.to <- lapply(theta.to, inla.source2function)
        if (!is.null(tags)) {
            names(theta.to) <- tags
        }
    } else {
        theta.to <- NULL
    }

    fnm <- paste(d, "/covmat-hyper-internal.dat", sep = "")
    if (file.exists(fnm)) {
        siz <- inla.read.binary.file(fnm)
        n <- siz[1L]
        stopifnot(length(siz) == n^2L + 1L)

        cov.intern <- matrix(siz[-1L], n, n)
        dd <- diag(cov.intern)
        s <- matrix(0.0, n, n)
        diag(s) <- 1.0 / sqrt(dd)
        cor.intern <- s %*% cov.intern %*% s
        diag(cor.intern) <- 1.0
    } else {
        cov.intern <- NULL
        cor.intern <- NULL
    }

    fnm <- paste(d, "/covmat-eigenvectors.dat", sep = "")
    if (file.exists(fnm)) {
        siz <- inla.read.binary.file(fnm)
        n <- siz[1L]
        stopifnot(length(siz) == n^2L + 1L)
        cov.intern.eigenvectors <- matrix(siz[-1L], n, n)
    } else {
        cov.intern.eigenvectors <- NULL
    }

    fnm <- paste(d, "/covmat-eigenvalues.dat", sep = "")
    if (file.exists(fnm)) {
        siz <- inla.read.binary.file(fnm)
        n <- siz[1L]
        stopifnot(length(siz) == n + 1L)
        cov.intern.eigenvalues <- siz[-1L]
    } else {
        cov.intern.eigenvalues <- NULL
    }

    fnm <- paste(d, "/reordering.dat", sep = "")
    if (file.exists(fnm)) {
        r <- as.integer(inla.read.binary.file(fnm))
    } else {
        r <- NULL
    }

    fnm <- paste(d, "/stdev_corr_pos.dat", sep = "")
    if (file.exists(fnm)) {
        stdev.corr.positive <- as.numeric(inla.read.fmesher.file(fnm))
    } else {
        stdev.corr.positive <- NULL
    }

    fnm <- paste(d, "/stdev_corr_neg.dat", sep = "")
    if (file.exists(fnm)) {
        stdev.corr.negative <- as.numeric(inla.read.fmesher.file(fnm))
    } else {
        stdev.corr.negative <- NULL
    }

    fnm <- paste(d, "/lincomb_derived_correlation_matrix.dat", sep = "")
    if (file.exists(fnm)) {
        lincomb.derived.correlation.matrix <- inla.read.fmesher.file(fnm)
    } else {
        lincomb.derived.correlation.matrix <- NULL
    }

    fnm <- paste(d, "/lincomb_derived_covariance_matrix.dat", sep = "")
    if (file.exists(fnm)) {
        lincomb.derived.covariance.matrix <- inla.read.fmesher.file(fnm)
    } else {
        lincomb.derived.covariance.matrix <- NULL
    }

    fnm <- paste(d, "/opt_directions.dat", sep = "")
    if (file.exists(fnm)) {
        opt.directions <- inla.read.fmesher.file(fnm)
        n <- dim(opt.directions)[1]
        colnames(opt.directions) <- paste0("dir:", 1:n)
        rownames(opt.directions) <- paste0("theta:", 1:n)
    } else {
        opt.directions <- NULL
    }

    fnm <- paste(d, "/mode-status.dat", sep = "")
    if (file.exists(fnm)) {
        mode.status <- scan(fnm, quiet = TRUE)
    } else {
        mode.status <- NA
    }

    fnm <- paste(d, "/nfunc.dat", sep = "")
    if (file.exists(fnm)) {
        nfunc <- as.numeric(scan(fnm, quiet = TRUE))
    } else {
        nfunc <- NA
    }

    fnm <- paste(d, "/log-posterior-mode.dat", sep = "")
    if (file.exists(fnm)) {
        lpm <- scan(fnm, quiet = TRUE)
    } else {
        lpm <- NA
    }

    fnm <- paste(d, "/config/configs.dat", sep = "")
    if (file.exists(fnm)) {
        fp <- file(fnm, "rb")
        iarr <- readBin(fp, integer(), 3)
        configs <- list(
            .preopt = FALSE, 
            n = iarr[1],
            nz = iarr[2],
            ntheta = iarr[3]
        )
        configs.i <- readBin(fp, integer(), configs$nz) ## 0-based
        configs.j <- readBin(fp, integer(), configs$nz) ## 0-based
        configs$nconfig <- readBin(fp, integer(), 1)

        nc <- readBin(fp, integer(), 1)
        if (nc > 0) {
            A <- readBin(fp, numeric(), configs$n * nc)
            e <- readBin(fp, numeric(), nc)
            configs$constr <- list(
                nc = nc,
                A = matrix(A, nc, configs$n),
                e = e
            )
        } else {
            configs$constr <- NULL
        }

        theta.tag <- readLines(paste(d, "/config/theta-tag.dat", sep = ""))
        configs$contents <- list(
            tag = readLines(paste(d, "/config/tag.dat", sep = "")),
            start = as.integer(readLines(paste(d, "/config/start.dat", sep = ""))) + 1L,
            length = as.integer(readLines(paste(d, "/config/n.dat", sep = "")))
        )

        if (configs$nconfig > 0L) {
            configs$config[[configs$nconfig]] <- list()
            for (k in 1L:configs$nconfig) {
                log.post <- readBin(fp, numeric(), 1)
                log.post.orig <- readBin(fp, numeric(), 1)
                if (configs$ntheta > 0L) {
                    theta <- readBin(fp, numeric(), configs$ntheta)
                    names(theta) <- theta.tag
                } else {
                    theta <- NULL
                }
                mean <- readBin(fp, numeric(), configs$n)
                improved.mean <- readBin(fp, numeric(), configs$n)
                skewness <- readBin(fp, numeric(), configs$n)
                ## add the offsets to the mean here
                offsets <- readBin(fp, numeric(), configs$n)
                mean <- mean + offsets
                improved.mean <- improved.mean + offsets

                Q <- readBin(fp, numeric(), configs$nz)
                Qinv <- readBin(fp, numeric(), configs$nz)
                Qprior <- readBin(fp, numeric(), configs$n)
                dif <- which(configs$i != configs$j)
                if (length(dif) > 0L) {
                    iadd <- configs.j[dif] ## yes, its the transpose part
                    jadd <- configs.i[dif] ## yes, its the transpose part
                    Qadd <- Q[dif]
                    Qinvadd <- Qinv[dif]
                } else {
                    iadd <- c()
                    jadd <- c()
                    Qadd <- c()
                    Qinvadd <- c()
                }
                configs$config[[k]] <- list(
                    .preopt = FALSE, 
                    theta = theta,
                    log.posterior = log.post,
                    log.posterior.orig = log.post.orig,
                    mean = mean,
                    improved.mean = improved.mean,
                    skewness = skewness,
                    Q = sparseMatrix(
                        i = c(configs.i, iadd),
                        j = c(configs.j, jadd),
                        x = c(Q, Qadd),
                        dims = c(configs$n, configs$n),
                        index1 = FALSE,
                        repr = "C"
                    ),
                    Qinv = sparseMatrix(
                        i = c(configs.i, iadd),
                        j = c(configs.j, jadd),
                        x = c(Qinv, Qinvadd),
                        dims = c(configs$n, configs$n),
                        index1 = FALSE,
                        repr = "C"
                    ),
                    Qprior.diag = Qprior
                )
            }

            ## rescale the log.posteriors
            configs$max.log.posterior <- max(sapply(configs$config, function(x) x$log.posterior.orig))
            max.log.posterior <- max(sapply(configs$config, function(x) x$log.posterior))
            for (k in 1L:configs$nconfig) {
                ## with integration weight
                configs$config[[k]]$log.posterior <- configs$config[[k]]$log.posterior - max.log.posterior
                ## without integration weights
                configs$config[[k]]$log.posterior.orig <- configs$config[[k]]$log.posterior.orig - configs$max.log.posterior
            }
        } else {
            configs$config <- NULL
        }
        close(fp)
    } else {
        configs <- NULL
    }


    fnm <- paste(d, "/config_preopt/configs.dat", sep = "")
    if (file.exists(fnm)) {
        fp <- file(fnm, "rb")
        iarr <- readBin(fp, integer(), 8)
        configs <- list(
            .preopt = TRUE, 
            mpred = iarr[1], 
            npred = iarr[2], 
            mnpred = iarr[3], 
            Npred = iarr[4], 
            n = iarr[5], 
            nz = iarr[6], 
            prior_nz = iarr[7],
            ntheta = iarr[8])
        configs.i <- readBin(fp, integer(), configs$nz) ## 0-based
        configs.j <- readBin(fp, integer(), configs$nz) ## 0-based
        configs.iprior <- readBin(fp, integer(), configs$prior_nz) ## 0-based
        configs.jprior <- readBin(fp, integer(), configs$prior_nz) ## 0-based
        configs$nconfig <- readBin(fp, integer(), 1)

        nc <- readBin(fp, integer(), 1)
        if (nc > 0) {
            A <- readBin(fp, numeric(), configs$n * nc)
            e <- readBin(fp, numeric(), nc)
            configs$constr <- list(
                nc = nc,
                A = matrix(A, nc, configs$n),
                e = e
            )
        } else {
            configs$constr <- NULL
        }
        configs$offsets <- readBin(fp, numeric(), configs$mnpred)

        theta.tag <- readLines(paste(d, "/config_preopt/theta-tag.dat", sep = ""))
        configs$contents <- list(
            tag = readLines(paste(d, "/config_preopt/tag.dat", sep = "")),
            start = as.integer(readLines(paste(d, "/config_preopt/start.dat", sep = ""))) + 1L,
            length = as.integer(readLines(paste(d, "/config_preopt/n.dat", sep = "")))
        )

        configs$A <- inla.read.fmesher.file(paste0(d, "/config_preopt/A.dat"))
        fnm <- paste0(d, "/config_preopt/pA.dat")
        if (file.exists(fnm)) {
            configs$pA <- inla.read.fmesher.file(fnm)
        } else {
            configs$pA <- matrix(nrow = 0, ncol = nrow(configs$A))
        }

        if (configs$nconfig > 0L) {
            configs$config[[configs$nconfig]] <- list()
            for (k in 1L:configs$nconfig) {
                log.post <- readBin(fp, numeric(), 1)
                log.post.orig <- readBin(fp, numeric(), 1)
                if (configs$ntheta > 0L) {
                    theta <- readBin(fp, numeric(), configs$ntheta)
                    names(theta) <- theta.tag
                } else {
                    theta <- NULL
                }
                mean <- readBin(fp, numeric(), configs$n)
                improved.mean <- readBin(fp, numeric(), configs$n)

                Q <- readBin(fp, numeric(), configs$nz)
                Qinv <- readBin(fp, numeric(), configs$nz)
                Qprior <- readBin(fp, numeric(), configs$prior_nz)

                output <- readBin(fp, numeric(), 2)
                if (output[1]) {
                    cpodens.moments <- matrix(readBin(fp, numeric(), configs$Npred * 3), ncol = 3, byrow = TRUE)
                } else {
                    cpodens.moments <- matrix(1, nrow = 0, ncol = 3)
                }
                colnames(cpodens.moments) <- c("mean", "variance", "skewness")

                if (output[2]) {
                    gcpodens.moments <- matrix(readBin(fp, numeric(), configs$Npred * 3), ncol = 3, byrow = TRUE)
                } else {
                    gcpodens.moments <- matrix(1, nrow = 0, ncol = 3)
                }
                colnames(gcpodens.moments) <- c("mean", "variance", "log.theta.correction")
                
                arg.str <- NULL
                have.arg.str <- readBin(fp, numeric(), 1)
                if (have.arg.str > 0) {
                    arg.str <- vector('character', configs$Npred)
                    for (i in 1:configs$Npred) {
                        arg.str[i] <- readBin(fp, character(), 1)
                    }
                }

                ll.info <- NULL
                have.ll.info <- readBin(fp, numeric(), 1)
                if (have.ll.info > 0) {
                    ll.info <- readBin(fp, double(), 3 * configs$Npred)
                    ll.info[is.nan(ll.info)] <- NA
                    ll.info <- matrix(ll.info, configs$Npred, 3, byrow = TRUE)
                    colnames(ll.info) <- c("gradient", "hessian", "deriv3")
                }

                A.lpred.mean.variance <- matrix(0.0, nrow = 0, ncol = 2, dimnames = list(NULL, c("mean", "variance")))
                lpred.mean.variance <- matrix(0.0, nrow = 0, ncol = 2, dimnames = list(NULL, c("mean", "variance")))
                have.lpred <- readBin(fp, numeric(), 1)
                if (have.lpred > 0) {
                    ## need to add offsets here as its not added when its stored
                    lpred.mean <- readBin(fp, double(), configs$mnpred) + configs$offsets[1:configs$mnpred]
                    lpred.variance <- readBin(fp, double(), configs$mnpred)
                    lpred.mean[is.nan(lpred.mean)] <- NA
                    lpred.variance[is.nan(lpred.variance)] <- NA
                    offset <- 0
                    if (configs$mpred >0) {
                        offset <- configs$mpred
                        idx <- 1:configs$mpred
                        A.lpred.mean.variance <- cbind(mean = lpred.mean[idx], variance = lpred.variance[idx])
                    }
                    idx <- offset + 1:configs$npred
                    lpred.mean.variance <- cbind(mean = lpred.mean[idx], variance = lpred.variance[idx])
                }

                dif <- which(configs$i != configs$j)
                if (length(dif) > 0L) {
                    iadd <- configs.j[dif] ## yes, its the transpose part
                    jadd <- configs.i[dif] ## yes, its the transpose part
                    Qadd <- Q[dif]
                    Qinvadd <- Qinv[dif]
                } else {
                    iadd <- c()
                    jadd <- c()
                    Qadd <- c()
                    Qinvadd <- c()
                }
                dif <- which(configs$iprior != configs$jprior)
                if (length(dif) > 0L) {
                    iprioradd <- configs.jprior[dif] ## yes, its the transpose part
                    jprioradd <- configs.iprior[dif] ## yes, its the transpose part
                    Qprioradd <- Qprior[dif]
                } else {
                    iprioradd <- c()
                    jprioradd <- c()
                    Qprioradd <- c()
                }

                configs$config[[k]] <- list(
                    theta = theta,
                    log.posterior = log.post,
                    log.posterior.orig = log.post.orig,
                    mean = mean,
                    improved.mean = improved.mean,
                    skewness = rep(NA, configs$n), 
                    Q = sparseMatrix(
                        i = c(configs.i, iadd),
                        j = c(configs.j, jadd),
                        x = c(Q, Qadd),
                        dims = c(configs$n, configs$n),
                        index1 = FALSE,
                        repr = "C"
                    ),
                    Qinv = sparseMatrix(
                        i = c(configs.i, iadd),
                        j = c(configs.j, jadd),
                        x = c(Qinv, Qinvadd),
                        dims = c(configs$n, configs$n),
                        index1 = FALSE,
                        repr = "C"
                    ),
                    Qprior = sparseMatrix(
                        i = c(configs.iprior, iprioradd),
                        j = c(configs.jprior, jprioradd),
                        x = c(Qprior, Qprioradd),
                        dims = c(configs$n, configs$n),
                        index1 = FALSE,
                        repr = "C"
                    ),
                    cpodens.moments = cpodens.moments,
                    gcpodens.moments = gcpodens.moments,
                    arg.str = arg.str,
                    ll.info = ll.info,
                    APredictor = A.lpred.mean.variance, 
                    Predictor = lpred.mean.variance
                )
            }

            ## rescale the log.posteriors
            configs$max.log.posterior <- max(sapply(configs$config, function(x) x$log.posterior.orig))
            for (k in 1L:configs$nconfig) {
                configs$config[[k]]$log.posterior <- configs$config[[k]]$log.posterior - configs$max.log.posterior
                configs$config[[k]]$log.posterior.orig <- configs$config[[k]]$log.posterior.orig - configs$max.log.posterior
            }
        } else {
            configs$config <- NULL
        }
        close(fp)
    }

    fnm <- paste(d, "/warnings.txt", sep = "")
    if (file.exists(fnm)) {
        warn <- readLines(fnm)
    } else {
        warn <- NULL
    }

    if (debug) {
        print(paste("collect misc from", d, "...done"))
    }

    return(list(
        cov.intern = cov.intern, cor.intern = cor.intern,
        cov.intern.eigenvalues = cov.intern.eigenvalues, cov.intern.eigenvectors = cov.intern.eigenvectors,
        reordering = r, theta.tags = tags, log.posterior.mode = lpm,
        stdev.corr.negative = stdev.corr.negative, stdev.corr.positive = stdev.corr.positive,
        to.theta = theta.to, from.theta = theta.from, mode.status = mode.status,
        lincomb.derived.correlation.matrix = lincomb.derived.correlation.matrix,
        lincomb.derived.covariance.matrix = lincomb.derived.covariance.matrix,
        opt.directions = opt.directions,
        configs = configs, nfunc = nfunc, warnings = warn
    ))
}

`inla.collect.logfile` <- function(file.log = NULL, debug = FALSE)
{
    if (is.null(file.log)) {
        return(list(logfile = NULL))
    }
    if (debug) {
        print(paste("Read logfile", file.log))
    }
    if (file.exists(file.log)) {
        ## replace tab with spaces
        logfile = gsub("\t", "        ", readLines(file.log))
        logfile <- logfile[-which(nchar(logfile) == 0)]
        return(list(logfile = logfile))
    } else {
        return(list(logfile = NULL))
    }
}

`inla.collect.size` <- function(dir, debug = FALSE)
{
    fnm <- paste(dir, "/size.dat", sep = "")
    siz <- inla.read.binary.file(fnm)
    if (length(siz) != 5L) {
        return(rep(0L, 5))
        ## stop(paste("length of siz is not 5L: fnm=", fnm))
    }
    if (is.na(siz[1L]) || siz[1L] < 0L) stop("siz[1L] = NA")
    if (is.na(siz[2L]) || siz[2L] <= 0L) siz[2L] <- siz[1L]
    if (is.na(siz[3L]) || siz[3L] <= 0L) siz[3L] <- siz[2L]
    if (is.na(siz[4L]) || siz[4L] <= 0L) siz[4L] <- 1L
    if (is.na(siz[5L]) || siz[5L] <= 0L) siz[5L] <- 1L

    return(list(n = siz[1L], N = siz[2L], Ntotal = siz[3L], ngroup = siz[4L], nrep = siz[5L]))
}

`inla.collect.hyperid` <- function(dir, debug = FALSE)
{
    fnm <- paste(dir, "/hyperid.dat", sep = "")
    id <- readLines(fnm)
    return(id)
}

`inla.collect.fixed` <- function(results.dir, debug = FALSE)
{
    alldir <- dir(results.dir)
    if (debug) {
        print("collect fixed effects")
    }

    ## read FIXED EFFECTS
    fix <- alldir[grep("^fixed.effect", alldir)]
    fix <- c(fix, alldir[grep("^intercept$", alldir)])
    n.fix <- length(fix)

    ## read the names of the fixed effects
    if (n.fix > 0L) {
        names.fixed <- inla.trim(character(n.fix))
        for (i in 1L:n.fix) {
            tag <- paste(results.dir, .Platform$file.sep, fix[i], .Platform$file.sep, "TAG", sep = "")
            if (!file.exists(tag)) {
                names.fixed[i] <- "NameMissing"
            } else {
                names.fixed[i] <- readLines(tag, n = 1L)
            }
        }
        ## read summary the fixed effects
        if (debug) {
            print(names.fixed)
        }

        summary.fixed <- numeric()
        marginals.fixed <- list()
        marginals.fixed[[n.fix]] <- NA

        for (i in 1L:n.fix) {
            first.time <- (i == 1L)
            file <- paste(results.dir, .Platform$file.sep, fix[i], sep = "")
            dir.fix <- dir(file)
            if (length(dir.fix) > 3L) {
                summ <- inla.read.binary.file(paste(file, .Platform$file.sep, "summary.dat", sep = ""))[-1L]
                if (first.time) {
                    col.nam <- c("mean", "sd")
                }

                ## read quantiles if existing
                if (length(grep("^quantiles.dat$", dir.fix)) > 0L) {
                    qq <- inla.interpret.vector(inla.read.binary.file(paste(file, .Platform$file.sep, "quantiles.dat", sep = "")),
                                                debug = debug
                                                )
                    summ <- c(summ, qq[, 2L])
                    if (first.time) {
                        col.nam <- c(col.nam, paste(as.character(qq[, 1L]), "quant", sep = ""))
                    }
                }

                ## read mode if existing
                if (length(grep("^mode.dat$", dir.fix)) > 0L) {
                    mm <- inla.interpret.vector(inla.read.binary.file(paste(file, .Platform$file.sep, "mode.dat", sep = "")),
                                                debug = debug
                                                )
                    summ <- c(summ, mm[, 2L])
                    if (first.time) {
                        col.nam <- c(col.nam, "mode")
                    }
                }

                if (length(grep("^cdf.dat$", dir.fix)) > 0L) {
                    qq <- inla.interpret.vector(inla.read.binary.file(paste(file, .Platform$file.sep, "cdf.dat", sep = "")),
                                                debug = debug
                                                )
                    summ <- c(summ, qq[, 2L])
                    if (first.time) {
                        col.nam <- c(col.nam, paste(as.character(qq[, 1L]), "cdf", sep = ""))
                    }
                }

                ## read also kld distance
                kld.fixed <- inla.read.binary.file(paste(file, .Platform$file.sep, "symmetric-kld.dat", sep = ""))[-1L]
                summ <- c(summ, kld.fixed)
                if (first.time) {
                    col.nam <- c(col.nam, "kld")
                }
                summary.fixed <- rbind(summary.fixed, summ)

                ## read the marginals
                xx <- inla.interpret.vector(inla.read.binary.file(paste(file, .Platform$file.sep, "marginal-densities.dat", sep = "")),
                                            debug = debug
                                            )
                if (is.null(xx)) {
                    xx <- cbind(c(NA, NA, NA), c(NA, NA, NA))
                }
                colnames(xx) <- c("x", "y")
                marginals.fixed[[i]] <- xx
            } else {
                if (first.time) {
                    col.nam <- c("mean", "sd", "kld")
                }
                summary.fixed <- rbind(summary.fixed, c(NA, NA, NA))
                xx <- cbind(c(NA, NA, NA), c(NA, NA, NA))
                colnames(xx) <- c("x", "y")
                marginals.fixed[[i]] <- xx
            }
        }
        rownames(summary.fixed) <- names.fixed
        colnames(summary.fixed) <- col.nam
        if (length(marginals.fixed) > 0L) {
            names(marginals.fixed) <- names.fixed
        }
    } else {
        if (debug) {
            print("No fixed effects")
        }
        names.fixed <- NULL
        summary.fixed <- NULL
        marginals.fixed <- NULL
    }

    ret <- list(
        names.fixed = names.fixed,
        summary.fixed = as.data.frame(summary.fixed),
        marginals.fixed = marginals.fixed
    )
    return(ret)
}

`inla.collect.lincomb` <- function(results.dir, debug = FALSE, derived = TRUE)
{
    ## rewrite from collect.random
    alldir <- dir(results.dir)
    if (derived) {
        lincomb <- alldir[grep("^lincomb.*derived[.]all", alldir)]
    } else {
        lincomb1 <- alldir[grep("^lincomb.*derived[.]all", alldir)]
        lincomb2 <- alldir[grep("^lincomb", alldir)]
        lincomb <- setdiff(lincomb2, lincomb1)
        if (debug) {
            print(paste("lincomb", lincomb))
        }
    }
    n.lincomb <- length(lincomb)
    if (debug) {
        print("collect lincombs")
    }

    ## read the names and model of the lincomb effects
    if (n.lincomb > 0L) {
        names.lincomb <- character(n.lincomb)
        model.lincomb <- inla.trim(character(n.lincomb))

        summary.lincomb <- list()
        summary.lincomb[[n.lincomb]] <- NA
        marginals.lincomb <- list()
        marginals.lincomb[[n.lincomb]] <- NA
        size.lincomb <- list()
        size.lincomb[[n.lincomb]] <- NA

        for (i in 1L:n.lincomb) {
            if (debug) {
                print(paste("read lincomb ", i, " of ", n.lincomb))
            }

            ## read the summary
            file <- paste(results.dir, .Platform$file.sep, lincomb[i], sep = "")
            dir.lincomb <- dir(file)

            if (debug) {
                print(paste("read from dir ", file))
            }

            if (length(dir.lincomb) > 4L) {
                dd <- matrix(inla.read.binary.file(file = paste(file, .Platform$file.sep, "summary.dat", sep = "")),
                             ncol = 3L, byrow = TRUE
                             )
                col.nam <- c("ID", "mean", "sd")

                ## read quantiles if existing
                if (debug) {
                    cat("...quantiles.dat if any\n")
                }
                if (length(grep("^quantiles.dat$", dir.lincomb)) == 1L) {
                    xx <- inla.interpret.vector(inla.read.binary.file(paste(file, .Platform$file.sep,
                                                                            "quantiles.dat",
                                                                            sep = ""
                                                                            )), debug = debug)
                    len <- dim(xx)[2L]
                    qq <- xx[, seq(2L, len, by = 2L), drop = FALSE]
                    col.nam <- c(col.nam, paste(as.character(xx[, 1L]), "quant", sep = ""))
                    dd <- cbind(dd, t(qq))
                }

                ## read mode if existing
                if (length(grep("^mode.dat$", dir.lincomb)) > 0L) {
                    mm <- inla.interpret.vector(inla.read.binary.file(paste(file, .Platform$file.sep, "mode.dat", sep = "")),
                                                debug = debug
                                                )
                    len <- dim(mm)[2L]
                    qq <- mm[, seq(2L, len, by = 2L), drop = FALSE]
                    dd <- cbind(dd, t(qq))
                    col.nam <- c(col.nam, "mode")
                }

                ## read cdf if existing
                if (debug) {
                    cat("...cdf.dat if any\n")
                }
                if (length(grep("^cdf.dat$", dir.lincomb)) == 1L) {
                    xx <- inla.interpret.vector(inla.read.binary.file(paste(file, .Platform$file.sep, "cdf.dat", sep = "")),
                                                debug = debug
                                                )
                    len <- dim(xx)[2L]
                    qq <- xx[, seq(2L, len, by = 2L), drop = FALSE]
                    col.nam <- c(col.nam, paste(as.character(xx[, 1L]), " cdf", sep = ""))
                    dd <- cbind(dd, t(qq))
                }

                if (debug) {
                    cat("...NAMES if any\n")
                }
                if (length(grep("^NAMES$", dir.lincomb)) == 1L) {
                    row.names <- readLines(paste(file, .Platform$file.sep, "NAMES", sep = ""))
                    ## remove the prefix 'lincomb.' as we do not need it in the names.
                    row.names <- sapply(row.names, function(x) gsub("^lincomb[.]", "", x))
                    names(row.names) <- NULL
                } else {
                    row.names <- NULL
                }

                ## read kld
                if (debug) {
                    cat("...kld\n")
                }
                kld1 <- matrix(inla.read.binary.file(file = paste(file, .Platform$file.sep, "symmetric-kld.dat", sep = "")),
                               ncol = 2L, byrow = TRUE
                               )
                qq <- kld1[, 2L, drop = FALSE]
                dd <- cbind(dd, qq)
                if (debug) {
                    cat("...kld done\n")
                }

                col.nam <- c(col.nam, "kld")
                colnames(dd) <- col.nam
                summary.lincomb[[i]] <- as.data.frame(dd)
                if (!is.null(row.names)) {
                    rownames(summary.lincomb[[i]]) <- row.names
                }

                xx <- inla.read.binary.file(paste(file, .Platform$file.sep, "marginal-densities.dat", sep = ""))
                rr <- inla.interpret.vector.list(xx, debug = debug)
                rm(xx)
                if (!is.null(rr)) {
                    nd <- length(rr)
                    names(rr) <- paste("index.", as.character(1L:nd), sep = "")
                    for (j in 1L:nd) {
                        colnames(rr[[j]]) <- c("x", "y")
                    }
                }
                marginals.lincomb[[i]] <- rr

                if (!is.null(row.names) && (length(marginals.lincomb) > 0L)) {
                    names(marginals.lincomb[[i]]) <- row.names
                }
            } else {
                N.file <- paste(file, .Platform$file.sep, "N", sep = "")
                if (!file.exists(N.file)) {
                    N <- 0L
                } else {
                    N <- scan(file = N.file, what = numeric(0L), quiet = TRUE)
                }
                summary.lincomb[[i]] <- data.frame("mean" = rep(NA, N), "sd" = rep(NA, N), "kld" = rep(NA, N))
                marginals.lincomb <- NULL
            }
            size.lincomb[[i]] <- inla.collect.size(file)
        }
        names(summary.lincomb) <- names.lincomb

        ## could be that marginals.lincomb is a list of lists of NULL
        if (!is.null(marginals.lincomb)) {
            if (all(sapply(marginals.lincomb, is.null))) {
                marginals.lincomb <- NULL
            }
        }

        if (!is.null(marginals.lincomb) && (length(marginals.lincomb) > 0L)) {
            names(marginals.lincomb) <- names.lincomb
        }
    } else {
        if (debug) {
            cat("No lincomb effets\n")
        }
        summary.lincomb <- NULL
        marginals.lincomb <- NULL
        size.lincomb <- NULL
    }

    if (derived) {
        res <- list(
            summary.lincomb.derived = as.data.frame(summary.lincomb[[1L]]),
            marginals.lincomb.derived = inla.ifelse(length(marginals.lincomb) > 0L, marginals.lincomb[[1L]], NULL),
            size.lincomb.derived = size.lincomb[[1L]]
        )
    } else {
        res <- list(
            summary.lincomb = as.data.frame(summary.lincomb[[1L]]),
            marginals.lincomb = inla.ifelse(length(marginals.lincomb) > 0L, marginals.lincomb[[1L]], NULL),
            size.lincomb = size.lincomb[[1L]]
        )
    }
    return(res)
}

`inla.collect.gcpo` <- function(results.dir, debug = FALSE)
{
    alldir <- dir(results.dir)
    if (length(grep("^gcpo$", alldir)) == 1L) {
        if (debug) {
            cat(paste("collect gcpo\n", sep = ""))
        }

        xx <- inla.read.binary.file(file = paste(results.dir, .Platform$file.sep, "gcpo", .Platform$file.sep, "gcpo.dat", sep = ""))
        n <- xx[1L]
        xx <- xx[-1L]

        values <- xx[1:n]
        values[is.nan(values)] <- NA
        kld <- xx[n + 1:n]
        mean <- xx[2 * n + 1:n]
        sd <- xx[3 * n + 1:n]
        offset <- 4 * n
        groups <- rep(list(list()), n)
        for(i in 1:n) {
            nn <- xx[offset + 1]
            if (nn > 0) {
                gg <- xx[offset + 1 + seq_len(nn)]
                cor <- xx[offset + 1 + nn + seq_len(nn)]
            } else {
                gg <- NULL
                cor <- NULL
            }
            groups[[i]] <- list(idx = gg, corr = cor)
            offset <- offset + 1 + 2 * nn
        }
    } else {
        values <- NULL
        kld <- NULL
        mean <- NULL
        sd <- NULL
        groups <- NULL
    }

    return(list(gcpo = values, kld = kld, mean = mean, sd = sd, groups = groups))
}

`inla.collect.cpo` <- function(results.dir, debug = FALSE)
{
    alldir <- dir(results.dir)
    if (length(grep("^cpo$", alldir)) == 1L) {
        if (debug) {
            cat(paste("collect cpo\n", sep = ""))
        }

        xx <- inla.read.binary.file(file = paste(results.dir, .Platform$file.sep, "cpo", .Platform$file.sep, "cpo.dat", sep = ""))
        n <- xx[1L]
        xx <- xx[-1L]
        len <- length(xx)
        cpo.res <- numeric(n)
        cpo.res[1L:n] <- NA
        cpo.res[xx[seq(1L, len, by = 2L)] + 1L] <- xx[seq(2L, len, by = 2L)]

        xx <- inla.read.binary.file(file = paste(results.dir, .Platform$file.sep, "cpo", .Platform$file.sep, "pit.dat", sep = ""))
        n <- xx[1L]
        xx <- xx[-1L]
        len <- length(xx)
        pit.res <- numeric(n)
        pit.res[1L:n] <- NA
        pit.res[xx[seq(1L, len, by = 2L)] + 1L] <- xx[seq(2L, len, by = 2L)]

        fnm <- paste(results.dir, .Platform$file.sep, "cpo", .Platform$file.sep, "failure.dat", sep = "")
        if (file.exists(fnm)) {
            xx <- inla.read.binary.file(fnm)
            n <- xx[1L]
            xx <- xx[-1L]
            len <- length(xx)
            failure.res <- numeric(n)
            failure.res[1L:n] <- NA
            failure.res[xx[seq(1L, len, by = 2L)] + 1L] <- xx[seq(2L, len, by = 2L)]
        }
        else {
            failure.res <- NULL
        }
        rm(xx)
    } else {
        cpo.res <- NULL
        pit.res <- NULL
        failure.res <- NULL
    }

    ## want NA not NaN
    cpo.res[is.nan(cpo.res)] <- NA
    pit.res[is.nan(pit.res)] <- NA
    failure.res[is.nan(failure.res)] <- NA

    return(list(cpo = cpo.res, pit = pit.res, failure = failure.res))
}

`inla.collect.po` <- function(results.dir, debug = FALSE)
{
    alldir <- dir(results.dir)
    if (length(grep("^po$", alldir)) == 1L) {
        if (debug) {
            cat(paste("collect po\n", sep = ""))
        }
        xx <- inla.read.binary.file(file = paste(results.dir, .Platform$file.sep, "po", .Platform$file.sep, "po.dat", sep = ""))
        n <- xx[1L]
        xx <- xx[-1L]
        xx <- xx[-seq(3, length(xx), by = 3L)] 
        len <- length(xx)
        po.res <- numeric(n)
        po.sign <- numeric(n)
        po.res[1L:n] <- NA
        po.res[xx[seq(1L, len, by = 2L)] + 1L] <- xx[seq(2L, len, by = 2L)]
    } else {
        po.res <- NULL
    }

    ## want NA not NaN
    po.res[is.nan(po.res)] <- NA

    return(list(po = po.res))
}

`inla.collect.waic` <- function(results.dir, debug = FALSE)
{
    ## yes, here we use the po-results!!!!
    alldir <- dir(results.dir)
    if (length(grep("^po$", alldir)) == 1L) {
        if (debug) {
            cat(paste("collect waic from po-results\n", sep = ""))
        }

        xx <- inla.read.binary.file(file = paste(results.dir, .Platform$file.sep, "po", .Platform$file.sep, "po.dat", sep = ""))
        n <- xx[1L]
        xx <- xx[-1L]
        len <- length(xx)
        po.res <- numeric(n)
        po2.res <- numeric(n)
        po.res[1L:n] <- NA
        po.res[xx[seq(1L, len, by = 3L)] + 1L] <- xx[seq(2L, len, by = 3L)]
        po2.res[1L:n] <- NA
        po2.res[xx[seq(1L, len, by = 3L)] + 1L] <- xx[seq(3L, len, by = 3L)]

        ## want NA not NaN
        po.res[is.nan(po.res)] <- NA
        po2.res[is.nan(po2.res)] <- NA

        ## compute waic
        return(list(
            waic = -2 * (sum(log(po.res), na.rm = TRUE) - sum(po2.res, na.rm = TRUE)),
            p.eff = sum(po2.res, na.rm = TRUE),
            local.waic = -2 * (log(po.res) - po2.res),
            local.p.eff = po2.res
        ))
    } else {
        return(NULL)
    }
}

`inla.collect.dic` <- function(results.dir, debug = FALSE)
{
    alldir <- dir(results.dir)
    ## get dic (if exists)
    if (length(grep("^dic$", alldir)) == 1L) {
        if (debug) {
            cat(paste("collect dic\n", sep = ""))
        }
        file <- paste(results.dir, .Platform$file.sep, "dic", .Platform$file.sep, "dic.dat", sep = "")
        dic.values <- inla.read.binary.file(file)

        file <- paste(results.dir, .Platform$file.sep, "dic", .Platform$file.sep, "deviance_e.dat", sep = "")
        if (inla.is.fmesher.file(file)) {
            dev.e <- c(inla.read.fmesher.file(file))
            dev.e[is.nan(dev.e)] <- NA
        } else {
            dev.e <- NULL
        }

        file <- paste(results.dir, .Platform$file.sep, "dic", .Platform$file.sep, "deviance_e_sat.dat", sep = "")
        if (inla.is.fmesher.file(file)) {
            dev.e.sat <- c(inla.read.fmesher.file(file))
            dev.e.sat[is.nan(dev.e.sat)] <- NA
        } else {
            dev.e.sat <- NULL
        }

        file <- paste(results.dir, .Platform$file.sep, "dic", .Platform$file.sep, "e_deviance.dat", sep = "")
        if (inla.is.fmesher.file(file)) {
            e.dev <- c(inla.read.fmesher.file(file))
            e.dev[is.nan(e.dev)] <- NA
        } else {
            e.dev <- NULL
        }

        file <- paste(results.dir, .Platform$file.sep, "dic", .Platform$file.sep, "sign.dat", sep = "")
        if (inla.is.fmesher.file(file)) {
            sig <- c(inla.read.fmesher.file(file))
            sig[is.nan(sig)] <- NA
        } else {
            sig <- NULL
        }

        file <- paste(results.dir, .Platform$file.sep, "dic", .Platform$file.sep, "e_deviance_sat.dat", sep = "")
        if (inla.is.fmesher.file(file)) {
            e.dev.sat <- c(inla.read.fmesher.file(file))
            e.dev.sat[is.nan(e.dev.sat)] <- NA
        } else {
            e.dev.sat <- NULL
        }

        f.idx <- NULL
        file <- paste(results.dir, .Platform$file.sep, "dic", .Platform$file.sep, "family_idx.dat", sep = "")
        if (inla.is.fmesher.file(file)) {
            f.idx <- c(inla.read.fmesher.file(file)) + 1L ## convert to R-indexing
            f.idx[is.nan(f.idx)] <- NA
        }

        ## if there there is no data at all, then all dic'values are
        ## NA. the returned values are 0, so we override them here.
        if (!is.null(f.idx) && all(is.na(f.idx))) {
            dic.values[] <- NA
        }

        if (!is.null(sig)) {
            ## avoid warnings
            ee <- e.dev.sat
            idx <- which(is.na(ee))
            ee[idx] <- 0
            sig[idx] <- 0
            deviance.residuals <- sqrt(ee) * sig
            deviance.residuals[idx] <- NA
        } else {
            deviance.residuals <- NULL
        }
        
        local.dic <- 2.0 * e.dev - dev.e
        local.dic.sat <- 2.0 * e.dev.sat - dev.e.sat
        local.p.eff <- e.dev - dev.e
        fam.dic <- dic.values[4L]
        fam.p.eff <- dic.values[3L]

        if (!is.null(f.idx) && !all(is.na(f.idx))) {
            n.fam <- max(f.idx, na.rm = TRUE)
            fam.dic <- numeric(n.fam)
            fam.dic.sat <- numeric(n.fam)
            fam.p.eff <- numeric(n.fam)
            for (i in 1:n.fam) {
                idx <- which(f.idx == i)
                fam.dic[i] <- sum(local.dic[idx])
                fam.dic.sat[i] <- sum(local.dic.sat[idx])
                fam.p.eff[i] <- sum(local.p.eff[idx])
            }
        } else {
            fam.dic.sat <- rep(NA, length(fam.dic))
        }

        dic <- list(
            "dic" = dic.values[4L],
            "p.eff" = dic.values[3L],
            "mean.deviance" = dic.values[1L],
            "deviance.mean" = dic.values[2L],
            "dic.sat" = dic.values[4L + 4L],
            "mean.deviance.sat" = dic.values[4L + 1L],
            "deviance.mean.sat" = dic.values[4L + 2L],
            "family.dic" = fam.dic,
            "family.dic.sat" = fam.dic.sat,
            "family.p.eff" = fam.p.eff,
            "family" = f.idx,
            "local.dic" = local.dic,
            "local.dic.sat" = local.dic.sat,
            "local.p.eff" = local.p.eff,
            "residuals.deviance" = deviance.residuals
        )
    } else {
        dic <- NULL
    }
    return(dic)
}

`inla.collect.q` <- function(results.dir, debug = FALSE)
{
    my.read.pnm <- function(...) {
        args <- list(...)
        filename <- args[[1]]
        inla.require("pixmap", stop.on.error = TRUE)
        if (file.exists(filename)) {
            ## disable warnings
            warn <- getOption("warn")
            options(warn = -1L) ## disable...
            ret <- pixmap::read.pnm(...)
            do.call("options", args = list(warn = warn))
        } else {
            ret <- NULL
        }
        return(ret)
    }

    alldir <- dir(results.dir)
    if (length(grep("^Q$", alldir)) == 1L) {
        if (debug) {
            cat(paste("collect q\n", sep = ""))
        }

        file <- paste(results.dir, .Platform$file.sep, "Q/precision-matrix.pbm", sep = "")
        Q.matrix <- my.read.pnm(file)

        file <- paste(results.dir, .Platform$file.sep, "Q/precision-matrix-reordered.pbm", sep = "")
        Q.matrix.reorder <- my.read.pnm(file)

        file <- paste(results.dir, .Platform$file.sep, "Q/precision-matrix_L.pbm", sep = "")
        L <- my.read.pnm(file)

        if (is.null(Q.matrix) && is.null(Q.matrix.reorder) && is.null(L)) {
            q <- NULL
        } else {
            q <- list(Q = Q.matrix, Q.reorder = Q.matrix.reorder, L = L)
        }
    } else {
        q <- NULL
    }

    return(q)
}

`inla.collect.graph` <- function(results.dir, debug = FALSE)
{
    alldir <- dir(results.dir)
    if (length(grep("^graph.dat$", alldir)) == 1L) {
        if (debug) {
            cat(paste("collect graph\n", sep = ""))
        }
        file <- paste(results.dir, .Platform$file.sep, "graph.dat", sep = "")
        g <- inla.read.graph(file)
    } else {
        g <- NULL
    }

    return(list(graph = g))
}

`inla.collect.hyperpar` <- function(results.dir, debug = FALSE)
{
    alldir <- dir(results.dir)
    all.hyper <- alldir[grep("^hyperparameter", alldir)]
    hyper <- all.hyper[grep("user-scale$", all.hyper)]
    n.hyper <- length(hyper)
    if (n.hyper > 0L) {
        ## get names for hyperpar
        names.hyper <- character(n.hyper)
        for (i in 1L:n.hyper) {
            tag <- paste(results.dir, .Platform$file.sep, hyper[i], .Platform$file.sep, "TAG", sep = "")
            if (!file.exists(tag)) {
                names.hyper[i] <- "missing NAME"
            } else {
                names.hyper[i] <- readLines(tag, n = 1L)
            }
        }

        ## get summary and marginals
        summary.hyper <- numeric()
        marginal.hyper <- list()
        marginal.hyper[[n.hyper]] <- NA

        for (i in 1L:n.hyper) {
            first.time <- (i == 1L)
            dir.hyper <- paste(results.dir, .Platform$file.sep, hyper[i], sep = "")
            file <- paste(dir.hyper, .Platform$file.sep, "summary.dat", sep = "")
            hyperid <- inla.collect.hyperid(dir.hyper)
            dd <- inla.read.binary.file(file)[-1L]
            summ <- dd
            if (first.time) {
                col.nam <- c("mean", "sd")
            }
            if (length(grep("^quantiles.dat$", dir(dir.hyper))) > 0L) {
                qq <- inla.interpret.vector(inla.read.binary.file(paste(dir.hyper, .Platform$file.sep, "quantiles.dat", sep = "")),
                                            debug = debug
                                            )
                summ <- c(summ, qq[, 2L])
                if (first.time) {
                    col.nam <- c(col.nam, paste(as.character(qq[, 1L]), "quant", sep = ""))
                }
            }
            if (length(grep("^mode.dat$", dir(dir.hyper))) > 0L) {
                qq <- inla.interpret.vector(inla.read.binary.file(paste(dir.hyper, .Platform$file.sep, "mode.dat", sep = "")),
                                            debug = debug
                                            )
                summ <- c(summ, qq[, 2L])
                if (first.time) {
                    col.nam <- c(col.nam, "mode")
                }
            }

            if (length(grep("^cdf.dat$", dir(dir.hyper))) > 0L) {
                qq <- inla.interpret.vector(inla.read.binary.file(paste(dir.hyper, .Platform$file.sep, "cdf.dat", sep = "")),
                                            debug = debug
                                            )
                summ <- c(summ, qq[, 2L])
                if (first.time) {
                    col.nam <- c(col.nam, paste(as.character(qq[, 1L]), "cdf", sep = ""))
                }
            }
            summary.hyper <- rbind(summary.hyper, summ)
            file <- paste(results.dir, .Platform$file.sep, hyper[i], .Platform$file.sep, "marginal-densities.dat", sep = "")
            xx <- inla.read.binary.file(file)
            marg1 <- inla.interpret.vector(xx, debug = debug)
            attr(marg1, "hyperid") <- hyperid
            rm(xx)
            if (!is.null(marg1)) {
                colnames(marg1) <- c("x", "y")
            }
            marginal.hyper[[i]] <- marg1
        }
        names(marginal.hyper) <- names.hyper
        rownames(summary.hyper) <- names.hyper
        colnames(summary.hyper) <- col.nam
    } else {
        marginal.hyper <- NULL
        summary.hyper <- NULL
    }

    ## collect also the hyperparameters in the internal scale
    all.hyper <- alldir[grep("^hyperparameter", alldir)]
    hyper <- all.hyper[-grep("user-scale$", all.hyper)]
    n.hyper <- length(hyper)
    if (n.hyper > 0L) {
        ## get names for hyperpar
        names.hyper <- character(n.hyper)
        for (i in 1L:n.hyper) {
            tag <- paste(results.dir, .Platform$file.sep, hyper[i], .Platform$file.sep, "TAG", sep = "")
            if (!file.exists(tag)) {
                names.hyper[i] <- "missing NAME"
            } else {
                names.hyper[i] <- readLines(tag, n = 1L)
            }
        }

        ## get summary and marginals
        internal.summary.hyper <- numeric()
        internal.marginal.hyper <- list()
        internal.marginal.hyper[[n.hyper]] <- NA
        for (i in 1L:n.hyper) {
            first.time <- (i == 1L)
            dir.hyper <- paste(results.dir, .Platform$file.sep, hyper[i], sep = "")
            file <- paste(dir.hyper, .Platform$file.sep, "summary.dat", sep = "")
            hyperid <- inla.collect.hyperid(dir.hyper)
            dd <- inla.read.binary.file(file)[-1L]
            summ <- dd
            if (first.time) {
                col.nam <- c("mean", "sd")
            }
            if (length(grep("^quantiles.dat$", dir(dir.hyper))) > 0L) {
                qq <- inla.interpret.vector(inla.read.binary.file(paste(dir.hyper, .Platform$file.sep, "quantiles.dat", sep = "")),
                                            debug = debug
                                            )
                summ <- c(summ, qq[, 2L])
                if (first.time) {
                    col.nam <- c(col.nam, paste(as.character(qq[, 1L]), "quant", sep = ""))
                }
            }
            if (length(grep("^mode.dat$", dir(dir.hyper))) > 0L) {
                qq <- inla.interpret.vector(inla.read.binary.file(paste(dir.hyper, .Platform$file.sep, "mode.dat", sep = "")),
                                            debug = debug
                                            )
                summ <- c(summ, qq[, 2L])
                if (first.time) {
                    col.nam <- c(col.nam, "mode")
                }
            }
            if (length(grep("^cdf.dat$", dir(dir.hyper))) > 0L) {
                qq <- inla.interpret.vector(inla.read.binary.file(paste(dir.hyper, .Platform$file.sep, "cdf.dat", sep = "")),
                                            debug = debug
                                            )
                summ <- c(summ, qq[, 2L])
                if (first.time) {
                    col.nam <- c(col.nam, paste(as.character(qq[, 1L]), "cdf", sep = ""))
                }
            }
            if (first.time) {
                internal.summary.hyper <- matrix(NA, n.hyper, length(summ))
            }
            internal.summary.hyper[i, ] <- summ
            file <- paste(results.dir, .Platform$file.sep, hyper[i], .Platform$file.sep, "marginal-densities.dat", sep = "")
            xx <- inla.read.binary.file(file)
            marg1 <- inla.interpret.vector(xx, debug = debug)
            attr(marg1, "hyperid") <- hyperid
            rm(xx)
            if (!is.null(marg1)) {
                colnames(marg1) <- c("x", "y")
            }

            internal.marginal.hyper[[i]] <- marg1
        }
        names(internal.marginal.hyper) <- names.hyper
        rownames(internal.summary.hyper) <- names.hyper
        colnames(internal.summary.hyper) <- col.nam
    } else {
        internal.summary.hyper <- NULL
        internal.marginal.hyper <- NULL
    }

    ret <- list(
        summary.hyperpar = as.data.frame(summary.hyper),
        marginals.hyperpar = marginal.hyper,
        internal.summary.hyperpar = as.data.frame(internal.summary.hyper),
        internal.marginals.hyperpar = internal.marginal.hyper
    )
    return(ret)
}

`inla.collect.mlik` <- function(results.dir, debug = FALSE)
{
    alldir <- dir(results.dir)
    if (length(grep("^marginal-likelihood$", alldir)) == 1L) {
        if (debug) {
            cat(paste("collect mlik\n", sep = ""))
        }
        file <- paste(results.dir, .Platform$file.sep, "marginal-likelihood",
                      .Platform$file.sep, "marginal-likelihood.dat",
                      sep = ""
                      )
        mlik.res <- matrix(inla.read.binary.file(file), 2L, 1L)
        rownames(mlik.res) <- c(
            "log marginal-likelihood (integration)",
            "log marginal-likelihood (Gaussian)"
        )
    }
    else {
        mlik.res <- NULL
    }

    return(list(mlik = mlik.res))
}

`inla.collect.predictor` <- function(results.dir, debug = FALSE)
{
    alldir <- dir(results.dir)

    ## FIRST: get the linear predictor
    subdir <- paste(results.dir, .Platform$file.sep, "predictor", sep = "")

    if (length(dir(subdir)) > 3L) {
        if (debug) {
            cat(paste("collect linear predictor\n", sep = ""))
        }

        if (debug) {
            cat("...read summary.dat\n")
        }
        file <- paste(subdir, .Platform$file.sep, "summary.dat", sep = "")
        dd <- matrix(inla.read.binary.file(file = file), ncol = 3L, byrow = TRUE)[, -1L, drop = FALSE]
        col.nam <- c("mean", "sd")

        ## info about size
        size.info <- inla.collect.size(subdir)
        if (!is.null(size.info)) {
            A <- (size.info$nrep == 2)
            n <- size.info$n
            nA <- size.info$Ntotal - size.info$n
        } else {
            ## should not happen
            stop("This should not happen")
        }

        ## get quantiles if computed
        if (length(grep("^quantiles.dat$", dir(subdir))) == 1L) {
            if (debug) {
                cat("...read quantiles.dat\n")
            }
            file <- paste(subdir, .Platform$file.sep, "quantiles.dat", sep = "")
            xx <- inla.interpret.vector(inla.read.binary.file(file), debug = debug)
            len <- dim(xx)[2L]
            qq <- xx[, seq(2L, len, by = 2L), drop = FALSE]
            col.nam <- c(col.nam, paste(as.character(xx[, 1L]), "quant", sep = ""))
            dd <- cbind(dd, t(qq))
            rm(xx)
        }

        if (length(grep("^mode.dat$", dir(subdir))) == 1L) {
            if (debug) {
                cat("...read mode.dat\n")
            }
            file <- paste(subdir, .Platform$file.sep, "mode.dat", sep = "")
            xx <- inla.interpret.vector(inla.read.binary.file(file), debug = debug)
            len <- dim(xx)[2L]
            qq <- xx[, seq(2L, len, by = 2L)]
            col.nam <- c(col.nam, "mode")
            dd <- cbind(dd, qq)
            rm(xx)
        }

        ## get cdf if computed
        if (length(grep("^cdf.dat$", dir(subdir))) == 1L) {
            if (debug) {
                cat("...read cdf.dat\n")
            }
            file <- paste(subdir, .Platform$file.sep, "cdf.dat", sep = "")
            xx <- inla.interpret.vector(inla.read.binary.file(file), debug = debug)
            len <- dim(xx)[2L]
            qq <- xx[, seq(2L, len, by = 2L), drop = FALSE]
            col.nam <- c(col.nam, paste(as.character(xx[, 1L]), " cdf", sep = ""))
            dd <- cbind(dd, t(qq))
            rm(xx)
        } else {
            if (debug) {
                cat("... no cdf.dat\n")
            }
        }

        ## get kld
        if (debug) {
            cat("...read kld\n")
        }
        kld <- matrix(inla.read.binary.file(file = paste(subdir, .Platform$file.sep, "symmetric-kld.dat", sep = "")),
                      ncol = 2L, byrow = TRUE
                      )
        dd <- cbind(dd, kld[, 2L, drop = FALSE])
        col.nam <- c(col.nam, "kld")
        colnames(dd) <- col.nam
        summary.linear.predictor <- as.data.frame(dd)

        if (A) {
            rownames(summary.linear.predictor) <- c(
                paste("APredictor.", inla.num(1L:nA), sep = ""),
                paste("Predictor.", inla.num(1:n), sep = "")
            )
        } else {
            rownames(summary.linear.predictor) <- paste("Predictor.", inla.num(1L:size.info$Ntotal), sep = "")
        }

        if (debug) {
            cat("...read marginal-densities.dat\n")
        }
        file <- paste(subdir, .Platform$file.sep, "marginal-densities.dat", sep = "")
        xx <- inla.read.binary.file(file)
        rr <- inla.interpret.vector.list(xx, debug = debug)
        rm(xx)
        if (!is.null(rr)) {
            if (A) {
                names(rr) <- c(
                    paste("APredictor.", inla.num(1L:nA), sep = ""),
                    paste("Predictor.", inla.num(1L:n), sep = "")
                )
            } else {
                names(rr) <- paste("Predictor.", as.character(1L:length(rr)), sep = "")
            }
            names.rr <- names(rr)
            for (i in 1L:length(rr)) {
                colnames(rr[[i]]) <- c("x", "y")
            }
        }
        marginals.linear.predictor <- rr
    } else {
        summary.linear.predictor <- NULL
        marginals.linear.predictor <- NULL
        size.info <- NULL
    }

    summary.fitted.values <- NULL
    marginals.fitted.values <- NULL
    ## SECOND: get the inverse linear predictor(if computed)
    if (length(grep("^predictor-user-scale$", alldir)) == 1L) {
        subdir <- paste(results.dir, .Platform$file.sep, "predictor-user-scale", sep = "")
        if (length(dir(subdir)) > 3L) {
            if (debug) {
                cat(paste("collect fitted values\n", sep = ""))
            }

            file <- paste(subdir, .Platform$file.sep, "summary.dat", sep = "")
            dd <- matrix(inla.read.binary.file(file = file), ncol = 3L, byrow = TRUE)[, -1L, drop = FALSE]
            col.nam <- c("mean", "sd")

            ## get quantiles if computed
            if (length(grep("^quantiles.dat$", dir(subdir))) == 1L) {
                file <- paste(subdir, .Platform$file.sep, "quantiles.dat", sep = "")
                xx <- inla.interpret.vector(inla.read.binary.file(file), debug = debug)
                len <- dim(xx)[2L]
                qq <- xx[, seq(2L, len, by = 2L), drop = FALSE]
                col.nam <- c(col.nam, paste(as.character(xx[, 1L]), "quant", sep = ""))
                dd <- cbind(dd, t(qq))
                rm(xx)
            }
            if (length(grep("^mode.dat$", dir(subdir))) == 1L) {
                file <- paste(subdir, .Platform$file.sep, "mode.dat", sep = "")
                xx <- inla.interpret.vector(inla.read.binary.file(file), debug = debug)
                len <- dim(xx)[2L]
                qq <- xx[, seq(2L, len, by = 2L)]
                col.nam <- c(col.nam, "mode")
                dd <- cbind(dd, qq)
                rm(xx)
            }

            ## get cdf if computed
            if (length(grep("^cdf.dat$", dir(subdir))) == 1L) {
                file <- paste(subdir, .Platform$file.sep, "cdf.dat", sep = "")
                xx <- inla.interpret.vector(inla.read.binary.file(file), debug = debug)
                len <- dim(xx)[2L]
                qq <- xx[, seq(2L, len, by = 2L), drop = FALSE]
                col.nam <- c(col.nam, paste(as.character(xx[, 1L]), " cdf", sep = ""))
                dd <- cbind(dd, t(qq))
                rm(xx)
            }

            colnames(dd) <- col.nam
            if (A) {
                rownames(dd) <- c(
                    paste("fitted.APredictor.", inla.num(1L:nA), sep = ""),
                    paste("fitted.Predictor.", inla.num(1L:n), sep = "")
                )
            } else {
                rownames(dd) <- paste("fitted.Predictor.", inla.num(1L:n), sep = "")
            }
            summary.fitted.values <- as.data.frame(dd)

            file <- paste(subdir, .Platform$file.sep, "marginal-densities.dat", sep = "")
            xx <- inla.read.binary.file(file)
            rr <- inla.interpret.vector.list(xx, debug = debug)
            rm(xx)
            if (!is.null(rr)) {
                if (A) {
                    names(rr) <- c(
                        paste("fitted.APredictor.", inla.num(1L:nA), sep = ""),
                        paste("fitted.Predictor.", inla.num(1:n), sep = "")
                        )
                } else {
                    names(rr) <- paste("fitted.Predictor.", inla.num(1L:length(rr)), sep = "")
                }
                names.rr <- names(rr)
                for (i in 1L:length(rr)) {
                    colnames(rr[[i]]) <- c("x", "y")
                    }
            }
            marginals.fitted.values <- rr
        }
    } 

    res <- list(
        summary.linear.predictor = as.data.frame(summary.linear.predictor),
        marginals.linear.predictor = marginals.linear.predictor,
        summary.fitted.values = as.data.frame(summary.fitted.values),
        marginals.fitted.values = marginals.fitted.values,
        size.linear.predictor = size.info
    )

    return(res)
}

`inla.collect.random` <- function(results.dir, debug = FALSE) 
{
    alldir <- dir(results.dir)
    random <- alldir[grep("^random.effect", alldir)]
    n.random <- length(random)
    if (debug) {
        print("collect random effects")
    }

    ## read the names and model of the random effects
    if (n.random > 0L) {
        names.random <- character(n.random)
        model.random <- inla.trim(character(n.random))
        for (i in 1L:n.random) {
            tag <- paste(results.dir, .Platform$file.sep, random[i], .Platform$file.sep, "TAG", sep = "")
            if (!file.exists(tag)) {
                names.random[i] <- "missing NAME"
            } else {
                names.random[i] <- readLines(tag, n = 1L)
            }
            modelname <- inla.trim(paste(results.dir, .Platform$file.sep, random[i], .Platform$file.sep, "MODEL", sep = ""))
            if (!file.exists(modelname)) {
                model.random[i] <- "NoModelName"
            } else {
                model.random[i] <- inla.trim(readLines(modelname, n = 1L))
            }
        }

        summary.random <- list()
        summary.random[[n.random]] <- NA
        size.random <- list()
        size.random[[n.random]] <- NA

        marginals.random <- list()
        marginals.random[[n.random]] <- NA

        for (i in 1L:n.random) {
            if (debug) {
                print(paste("read random ", i, " of ", n.random))
            }
            ## read the summary
            file <- paste(results.dir, .Platform$file.sep, random[i], sep = "")
            dir.random <- dir(file)

            if (length(dir.random) > 5L) {
                dd <- matrix(inla.read.binary.file(file = paste(file, .Platform$file.sep, "summary.dat", sep = "")), ncol = 3L, byrow = TRUE)
                col.nam <- c("ID", "mean", "sd")
                ## read quantiles if existing
                if (debug) {
                    cat("...quantiles.dat if any\n")
                }
                if (length(grep("^quantiles.dat$", dir.random)) == 1L) {
                    xx <- inla.interpret.vector(inla.read.binary.file(paste(file, .Platform$file.sep, "quantiles.dat", sep = "")),
                                                debug = debug
                                                )
                    len <- dim(xx)[2L]
                    qq <- xx[, seq(2L, len, by = 2L), drop = FALSE]
                    col.nam <- c(col.nam, paste(as.character(xx[, 1L]), "quant", sep = ""))
                    dd <- cbind(dd, t(qq))
                }
                if (length(grep("^mode.dat$", dir.random)) == 1L) {
                    xx <- inla.interpret.vector(inla.read.binary.file(paste(file, .Platform$file.sep, "mode.dat", sep = "")),
                                                debug = debug
                                                )
                    len <- dim(xx)[2L]
                    qq <- xx[, seq(2L, len, by = 2L)]
                    col.nam <- c(col.nam, "mode")
                    dd <- cbind(dd, qq)
                }

                ## read cdf if existing
                if (debug) {
                    cat("...cdf.dat if any\n")
                }
                if (length(grep("^cdf.dat$", dir.random)) == 1L) {
                    xx <- inla.interpret.vector(inla.read.binary.file(paste(file, .Platform$file.sep, "cdf.dat", sep = "")),
                                                debug = debug
                                                )
                    len <- dim(xx)[2L]
                    qq <- xx[, seq(2L, len, by = 2L), drop = FALSE]
                    col.nam <- c(col.nam, paste(as.character(xx[, 1L]), " cdf", sep = ""))
                    dd <- cbind(dd, t(qq))
                }

                ## read kld
                if (debug) {
                    cat("...kld\n")
                }
                kld1 <- matrix(inla.read.binary.file(file = paste(file, .Platform$file.sep, "symmetric-kld.dat", sep = "")),
                               ncol = 2L, byrow = TRUE
                               )
                qq <- kld1[, 2L, drop = FALSE]
                dd <- cbind(dd, qq)
                if (debug) {
                    cat("...kld done\n")
                }


                col.nam <- c(col.nam, "kld")
                colnames(dd) <- col.nam
                summary.random[[i]] <- as.data.frame(dd)

                xx <- inla.read.binary.file(paste(file, .Platform$file.sep, "marginal-densities.dat", sep = ""))
                rr <- inla.interpret.vector.list(xx, debug = debug)
                rm(xx)
                if (!is.null(rr)) {
                    nd <- length(rr)
                    names(rr) <- paste("index.", as.character(1L:nd), sep = "")
                    names.rr <- names(rr)
                    for (j in 1L:nd) {
                        colnames(rr[[j]]) <- c("x", "y")
                        }
                }
                marginals.random[[i]] <- if (is.null(rr)) NA_real_ else rr

                ## if id.names are present,  override the default names
                id.names <- inla.readLines(paste(file, .Platform$file.sep, "id-names.dat", sep = ""))
                if (!is.null(id.names)) {
                    len.id.names <- length(id.names)
                    summary.random[[i]]$ID[1L:len.id.names] <- id.names
                    # Check for existing marginals.random[[i]] data:
                    # Cannot use is.na() directly, since the result must
                    # always be a single logical value.
                    if (length(marginals.random) >= i && !is.null(rr)) {
                        names(marginals.random[[i]][1L:len.id.names]) <- id.names
                    }
                }
            } else {
                N.file <- paste(file, .Platform$file.sep, "N", sep = "")
                if (!file.exists(N.file)) {
                    N <- 0L
                } else {
                    N <- scan(file = N.file, what = numeric(0L), quiet = TRUE)
                }
                summary.random[[i]] <- data.frame("mean" = rep(NA, N), "sd" = rep(NA, N), "kld" = rep(NA, N))
                marginals.random <- NULL
            }

            size.random[[i]] <- inla.collect.size(file)
        }
        names(summary.random) <- names.random

        ## could be that marginals.random is a list of lists of NULL or NA
        if (!is.null(marginals.random)) {
            if (all(sapply(marginals.random,
                           function(x) (
                               is.null(x) ||
                               (is.numeric(x) &&
                                (length(x) == 1) &&
                                is.na(x))
                               )))) {
                marginals.random <- NULL
            }
        }

        if (!is.null(marginals.random) && (length(marginals.random) > 0L)) {
            names(marginals.random) <- names.random
        }
    } else {
        if (debug) {
            cat("No random effets\n")
        }
        model.random <- NULL
        summary.random <- NULL
        marginals.random <- NULL
        size.random <- NULL
    }

    res <- list(
        model.random = model.random,
        summary.random = lapply(summary.random, as.data.frame),
        marginals.random = marginals.random,
        size.random = size.random
    )
    return(res)
}

`inla.collect.spde2.blc` <- function(results.dir, debug = FALSE) 
{
    ## a copy from collect.random
    alldir <- dir(results.dir)
    random <- alldir[grep("^spde2.blc", alldir)]
    n.random <- length(random)
    if (debug) {
        print("collect random effects")
    }

    ## read the names and model of the random effects
    if (n.random > 0L) {
        names.random <- character(n.random)
        model.random <- inla.trim(character(n.random))
        for (i in 1L:n.random) {
            tag <- paste(results.dir, .Platform$file.sep, random[i], .Platform$file.sep, "TAG", sep = "")
            if (!file.exists(tag)) {
                names.random[i] <- "missing NAME"
            } else {
                names.random[i] <- readLines(tag, n = 1L)
            }
            modelname <- inla.trim(paste(results.dir, .Platform$file.sep, random[i], .Platform$file.sep, "MODEL", sep = ""))
            if (!file.exists(modelname)) {
                model.random[i] <- "NoModelName"
            } else {
                model.random[i] <- inla.trim(readLines(modelname, n = 1L))
            }
        }


        summary.random <- list()
        summary.random[[n.random]] <- NA
        size.random <- list()
        size.random[[n.random]] <- NA

        marginals.random <- list()
        marginals.random[[n.random]] <- NA

        for (i in 1L:n.random) {
            if (debug) {
                print(paste("read random ", i, " of ", n.random))
            }
            ## read the summary
            file <- paste(results.dir, .Platform$file.sep, random[i], sep = "")
            dir.random <- dir(file)

            if (length(dir.random) > 4L) {
                dd <- matrix(inla.read.binary.file(file = paste(file, .Platform$file.sep, "summary.dat", sep = "")), ncol = 3L, byrow = TRUE)
                col.nam <- c("ID", "mean", "sd")
                ## read quantiles if existing
                if (debug) {
                    cat("...quantiles.dat if any\n")
                }
                if (length(grep("^quantiles.dat$", dir.random)) == 1L) {
                    xx <- inla.interpret.vector(inla.read.binary.file(paste(file, .Platform$file.sep, "quantiles.dat", sep = "")),
                                                debug = debug
                                                )
                    len <- dim(xx)[2L]
                    qq <- xx[, seq(2L, len, by = 2L), drop = FALSE]
                    col.nam <- c(col.nam, paste(as.character(xx[, 1L]), "quant", sep = ""))
                    dd <- cbind(dd, t(qq))
                }
                if (length(grep("^mode.dat$", dir.random)) == 1L) {
                    xx <- inla.interpret.vector(inla.read.binary.file(paste(file, .Platform$file.sep, "mode.dat", sep = "")),
                                                debug = debug
                                                )
                    len <- dim(xx)[2L]
                    qq <- xx[, seq(2L, len, by = 2L), drop = FALSE]
                    col.nam <- c(col.nam, "mode")
                    dd <- cbind(dd, t(qq))
                }

                ## read cdf if existing
                if (debug) {
                    cat("...cdf.dat if any\n")
                }
                if (length(grep("^cdf.dat$", dir.random)) == 1L) {
                    xx <- inla.interpret.vector(inla.read.binary.file(paste(file, .Platform$file.sep, "cdf.dat", sep = "")),
                                                debug = debug
                                                )
                    len <- dim(xx)[2L]
                    qq <- xx[, seq(2L, len, by = 2L), drop = FALSE]
                    col.nam <- c(col.nam, paste(as.character(xx[, 1L]), " cdf", sep = ""))
                    dd <- cbind(dd, t(qq))
                }

                ## read kld
                if (debug) {
                    cat("...kld\n")
                }
                kld1 <- matrix(inla.read.binary.file(file = paste(file, .Platform$file.sep, "symmetric-kld.dat", sep = "")),
                               ncol = 2L, byrow = TRUE
                               )
                qq <- kld1[, 2L, drop = FALSE]
                dd <- cbind(dd, qq)
                if (debug) {
                    cat("...kld done\n")
                }


                col.nam <- c(col.nam, "kld")
                colnames(dd) <- col.nam
                summary.random[[i]] <- as.data.frame(dd)

                xx <- inla.read.binary.file(paste(file, .Platform$file.sep, "marginal-densities.dat", sep = ""))
                rr <- inla.interpret.vector.list(xx, debug = debug)
                rm(xx)
                if (!is.null(rr)) {
                    nd <- length(rr)
                    names(rr) <- paste("index.", as.character(1L:nd), sep = "")
                    names.rr <- names(rr)
                    for (j in 1L:nd) {
                        colnames(rr[[j]]) <- c("x", "y")
                        }
                }
                marginals.random[[i]] <- rr
            } else {
                N.file <- paste(file, .Platform$file.sep, "N", sep = "")
                if (!file.exists(N.file)) {
                    N <- 0L
                } else {
                    N <- scan(file = N.file, what = numeric(0L), quiet = TRUE)
                }
                summary.random[[i]] <- data.frame("mean" = rep(NA, N), "sd" = rep(NA, N), "kld" = rep(NA, N))
                marginals.random <- NULL
            }

            size.random[[i]] <- inla.collect.size(file)
        }
        names(summary.random) <- names.random

        ## could be that marginals.random is a list of lists of NULL or NA
        if (!is.null(marginals.random)) {
            if (all(sapply(marginals.random,
                           function(x) (is.null(x) ||
                                        (is.numeric(x) &&
                                         (length(x) == 1) &&
                                         is.na(x)))))) {
                marginals.random <- NULL
            }
        }

        if (!is.null(marginals.random) && (length(marginals.random) > 0L)) {
            names(marginals.random) <- names.random
        }
    } else {
        if (debug) {
            cat("No random effets\n")
        }
        model.random <- NULL
        summary.random <- NULL
        marginals.random <- NULL
        size.random <- NULL
    }

    res <- list(
        model.spde2.blc = model.random,
        summary.spde2.blc = lapply(summary.random, as.data.frame),
        marginals.spde2.blc = marginals.random,
        size.spde2.blc = size.random
    )
    return(res)
}

`inla.collect.spde3.blc` <- function(results.dir, debug = FALSE)
{
    ## a copy from collect.random
    alldir <- dir(results.dir)
    random <- alldir[grep("^spde3.blc", alldir)]
    n.random <- length(random)
    if (debug) {
        print("collect random effects")
    }

    ## read the names and model of the random effects
    if (n.random > 0L) {
        names.random <- character(n.random)
        model.random <- inla.trim(character(n.random))
        for (i in 1L:n.random) {
            tag <- paste(results.dir, .Platform$file.sep, random[i], .Platform$file.sep, "TAG", sep = "")
            if (!file.exists(tag)) {
                names.random[i] <- "missing NAME"
            } else {
                names.random[i] <- readLines(tag, n = 1L)
            }
            modelname <- inla.trim(paste(results.dir, .Platform$file.sep, random[i], .Platform$file.sep, "MODEL", sep = ""))
            if (!file.exists(modelname)) {
                model.random[i] <- "NoModelName"
            } else {
                model.random[i] <- inla.trim(readLines(modelname, n = 1L))
            }
        }


        summary.random <- list()
        summary.random[[n.random]] <- NA
        size.random <- list()
        size.random[[n.random]] <- NA

        marginals.random <- list()
        marginals.random[[n.random]] <- NA

        for (i in 1L:n.random) {
            if (debug) {
                print(paste("read random ", i, " of ", n.random))
            }
            ## read the summary
            file <- paste(results.dir, .Platform$file.sep, random[i], sep = "")
            dir.random <- dir(file)

            if (length(dir.random) > 4L) {
                dd <- matrix(inla.read.binary.file(file = paste(file, .Platform$file.sep, "summary.dat", sep = "")), ncol = 3L, byrow = TRUE)
                col.nam <- c("ID", "mean", "sd")
                ## read quantiles if existing
                if (debug) {
                    cat("...quantiles.dat if any\n")
                }
                if (length(grep("^quantiles.dat$", dir.random)) == 1L) {
                    xx <- inla.interpret.vector(inla.read.binary.file(paste(file, .Platform$file.sep, "quantiles.dat", sep = "")),
                                                debug = debug
                                                )
                    len <- dim(xx)[2L]
                    qq <- xx[, seq(2L, len, by = 2L), drop = FALSE]
                    col.nam <- c(col.nam, paste(as.character(xx[, 1L]), "quant", sep = ""))
                    dd <- cbind(dd, t(qq))
                }
                if (length(grep("^mode.dat$", dir.random)) == 1L) {
                    xx <- inla.interpret.vector(inla.read.binary.file(paste(file, .Platform$file.sep, "mode.dat", sep = "")),
                                                debug = debug
                                                )
                    len <- dim(xx)[2L]
                    qq <- xx[, seq(2L, len, by = 2L), drop = FALSE]
                    col.nam <- c(col.nam, "mode")
                    dd <- cbind(dd, t(qq))
                }

                ## read cdf if existing
                if (debug) {
                    cat("...cdf.dat if any\n")
                }
                if (length(grep("^cdf.dat$", dir.random)) == 1L) {
                    xx <- inla.interpret.vector(inla.read.binary.file(paste(file, .Platform$file.sep, "cdf.dat", sep = "")),
                                                debug = debug
                                                )
                    len <- dim(xx)[2L]
                    qq <- xx[, seq(2L, len, by = 2L), drop = FALSE]
                    col.nam <- c(col.nam, paste(as.character(xx[, 1L]), " cdf", sep = ""))
                    dd <- cbind(dd, t(qq))
                }

                ## read kld
                if (debug) {
                    cat("...kld\n")
                }
                kld1 <- matrix(inla.read.binary.file(file = paste(file, .Platform$file.sep, "symmetric-kld.dat", sep = "")),
                               ncol = 2L, byrow = TRUE
                               )
                qq <- kld1[, 2L, drop = FALSE]
                dd <- cbind(dd, qq)
                if (debug) {
                    cat("...kld done\n")
                }

                col.nam <- c(col.nam, "kld")
                colnames(dd) <- col.nam
                summary.random[[i]] <- as.data.frame(dd)

                xx <- inla.read.binary.file(paste(file, .Platform$file.sep, "marginal-densities.dat", sep = ""))
                rr <- inla.interpret.vector.list(xx, debug = debug)
                rm(xx)
                if (!is.null(rr)) {
                    nd <- length(rr)
                    names(rr) <- paste("index.", as.character(1L:nd), sep = "")
                    names.rr <- names(rr)
                    for (j in 1L:nd) {
                        colnames(rr[[j]]) <- c("x", "y")
                        }
                }
                marginals.random[[i]] <- rr
            } else {
                N.file <- paste(file, .Platform$file.sep, "N", sep = "")
                if (!file.exists(N.file)) {
                    N <- 0L
                } else {
                    N <- scan(file = N.file, what = numeric(0L), quiet = TRUE)
                }
                summary.random[[i]] <- data.frame("mean" = rep(NA, N), "sd" = rep(NA, N), "kld" = rep(NA, N))
                marginals.random <- NULL
            }

            size.random[[i]] <- inla.collect.size(file)
        }
        names(summary.random) <- names.random

        ## could be that marginals.random is a list of lists of NULL or NA
        if (!is.null(marginals.random)) {
            if (all(sapply(marginals.random,
                           function(x) (is.null(x) ||
                                        (is.numeric(x) &&
                                         (length(x) == 1) &&
                                         is.na(x)))))) {
                marginals.random <- NULL
            }
        }

        if (!is.null(marginals.random) && (length(marginals.random) > 0L)) {
            names(marginals.random) <- names.random
        }
    } else {
        if (debug) {
            cat("No random effets\n")
        }
        model.random <- NULL
        summary.random <- NULL
        marginals.random <- NULL
        size.random <- NULL
    }

    res <- list(
        model.spde3.blc = model.random,
        summary.spde3.blc = lapply(summary.random, as.data.frame),
        marginals.spde3.blc = marginals.random,
        size.spde3.blc = size.random
    )
    return(res)
}

`inla.image.reduce` <- function(im, image.dim = 512)
{
    ## reduce image IM to image.dim IMAGE.DIM and return the image as a matrix.
    ## order the indices so the output can be plotted by image()

    if (!inherits(im, "pixmapGrey") || (im@size[1L] != im@size[2L])) {
        return(im)
    } else {
        return(im@grey)
    }

    ## do not need this anymore as we do this in GMRFLib.
    if (FALSE) {
        if (image.dim >= im@size[1L]) {
            n <- as.integer(im@size[1L])
            x <- matrix(NA, n, n)
            for (j in 1L:n) {
                x[j, n - (1L:n) + 1L] <- im@grey[1L:n, j]
            }
            return(x)
        }
        block <- ceiling(im@size[1L] / image.dim)
        n <- floor(im@size[1L] / block)
        ii <- jj <- 0L
        x <- matrix(NA, n, n)
        for (i in seq(1L, im@size[1L] - block + 1L, by = block)) {
            ii <- ii + 1L
            jj <- 0L
            for (j in seq(1L, im@size[1L] - block + 1L, by = block)) {
                jj <- jj + 1L
                x[jj, n - ii + 1L] <- min(im@grey[i:(i + block - 1L), j:(j + block - 1L)])
            }
        }
        return(x)
    }
}

`inla.collect.offset.linear.predictor` <- function(results.dir, debug = FALSE)
{
    filename <- paste(results.dir, "/totaloffset/totaloffset.dat", sep = "")
    stopifnot(file.exists(filename))

    xx <- inla.read.binary.file(filename)
    return(list(offset.linear.predictor = xx))
}
