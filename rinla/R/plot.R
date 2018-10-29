## Export: plot!inla

##! \name{plot.inla}
##! \alias{plot.inla}
##! \alias{inla.plot}
##! \title{Default INLA plotting}
##! \description{
##!   Takes am \code{inla} object produced by \code{inla} and plot the results
##! }
##! \usage{
##! \method{plot}{inla}(x,
##!              plot.fixed.effects = TRUE,
##!              plot.lincomb = TRUE,
##!              plot.random.effects = TRUE,
##!              plot.hyperparameters = TRUE,
##!              plot.predictor = TRUE,
##!              plot.q = TRUE,
##!              plot.cpo = TRUE,
##!              plot.prior = FALSE, 
##!              single = FALSE,
##!              postscript = FALSE,
##!              pdf = FALSE,
##!              prefix = "inla.plots/figure-",
##!              intern = FALSE, 
##!              debug = FALSE, 
##!              ...)
##! }
##! \arguments{
##!   \item{x}{A fitted  \code{inla} object produced by \code{inla} }
##!   \item{plot.fixed.effects}{Boolean indicating if posterior marginals
##!     for the fixed effects in the model should be plotted }
##!   \item{plot.lincomb}{Boolean indicating if posterior marginals
##!     for the linear combinations should be plotted }
##!   \item{plot.random.effects}{Boolean indicating if posterior mean and quantiles
##!     for the random effects in the model should be plotted  }
##!   \item{plot.hyperparameters}{Boolean indicating if posterior marginals
##!     for the hyperparameters in the model should be plotted }
##!   \item{plot.predictor}{Boolean indicating if posterior mean and quantiles
##!     for the linear predictor in the model should be plotted }
##!   \item{plot.q}{Boolean indicating if precision matrix should be displayed}
##!   \item{plot.cpo}{Boolean indicating if CPO/PIT valuesshould be plotted}
##!   \item{plot.prior}{Plot also the prior density for the hyperparameters}
##!   \item{single}{Boolean indicating if there should be more than one plot per page
##!                 (FALSE) or just one (TRUE)}
##!   \item{postscript}{Boolean indicating if postscript files should be produced instead}
##!   \item{pdf}{Boolean indicating if PDF files should be produced instead}
##!   \item{prefix}{The prefix for the created files. Additional numbering and suffix is added.}
##!   \item{intern}{Plot also the hyperparameters in its internal scale.}
##!   \item{debug}{Write some debug information}
##!   \item{...}{Additional arguments to \code{postscript()}, \code{pdf()} or \code{dev.new()}.}
##! }
##! \value{The return value is a list of the files created (if any).}
##! \author{Havard Rue \email{hrue@r-inla.org} }
##! \seealso{\code{\link{inla}}}
##! \examples{
##!\dontrun{
##!result = inla(...)
##!plot(result)
##!plot(result, single=TRUE)
##!plot(result, single=TRUE, pdf=TRUE, paper = "a4")
##!   }
##! }
##! \keyword{plot}

`plot.inla` =
    function(x,
             plot.fixed.effects = TRUE,
             plot.lincomb = TRUE,
             plot.random.effects = TRUE,
             plot.hyperparameters = TRUE,
             plot.predictor = TRUE,
             plot.q = TRUE,
             plot.cpo = TRUE,
             plot.prior = FALSE, 
             single = FALSE,
             postscript = FALSE,
             pdf = FALSE,
             prefix = "inla.plots/figure-",
             intern = FALSE, 
             debug = FALSE, 
             ...)
{
    figure.count = 1L
    figures = c()

    if (plot.prior) {
        all.hyper = inla.all.hyper.postprocess(x$all.hyper)
    }

    if (postscript && pdf) {
        stop("Only one of 'postscript' and 'pdf' can be generated at the time.")
    }

    initiate.plot = function(...)
    {
        if (!postscript && !pdf) {
            inla.dev.new(...)
        } else {
            dir = dirname(prefix)
            if (!file.exists(dir) && nchar(dir) > 0L) {
                dir.create(dir, recursive=TRUE)
            } else {
                stopifnot(file.info(dir)$isdir)
            }
        }
        return (invisible())
    }

    new.plot = function(...)
    {
        if (!postscript && !pdf) {
            inla.dev.new(...)
        } else {
            found = FALSE
            while(!found) {
                filename = paste(prefix,
                        inla.ifelse(regexpr("/$", prefix), "", "/"),
                        figure.count,
                        inla.ifelse(postscript, ".eps", ".pdf"), sep="")
                if (file.exists(filename)) {
                    figure.count <<- figure.count + 1L ## YES
                } else {
                    found = TRUE
                }
            }
            if (postscript) {
                postscript(file = filename, ...)
            } else if (pdf) {
                pdf(file = filename, ...)
            } else {
                stop("This should not happen")
            }
            figures <<- c(figures, filename)
        }
        return (invisible())
    }

    close.plot = function(...)
    {
        if (postscript || pdf) {
            if (names(dev.cur()) != "null device") {
                dev.off()
            }
        }
        return (invisible())
    }

    close.and.new.plot = function(...)
    {
        close.plot(...)
        new.plot(...)
        return (invisible())
    }

    ##
    initiate.plot(...)

    if (plot.fixed.effects) {
        ## plot marginals for the fixed effects
        fix = x$marginals.fixed
        labels.fix = names(x$marginals.fixed)
        nf = length(labels.fix)
        if (nf>0) {
            if (nf == 1 || single) {
                plot.layout = c(1, 1)
            } else if (nf == 2) {
                plot.layout = c(2, 1)
            } else {
                plot.layout = c(3, 3)
            }
            np = prod(plot.layout)
            ip = 0
            for(i in 1:nf) {
                if (ip%%np == 0) {
                    close.and.new.plot(...)
                    par(mfrow=c(plot.layout[1], plot.layout[2]))
                }
                ip = ip + 1

                if (!all(is.na(fix[[i]]))) {
                    ss = x$summary.fixed[i,]
                    sub=paste("Mean = ", round(ss[names(ss)=="mean"], 3)," SD = ", round(ss[names(ss)=="sd"], 3), sep="")
                    m = inla.smarginal(fix[[i]])
                    plot(m, type="l", main=paste("PostDens [", inla.nameunfix(labels.fix[i]),"]", sep=""),
                         sub=sub, xlab="", ylab="")

                    if (plot.prior) {
                        xy = (inla.get.prior.xy(section = "fixed", hyperid = labels.fix[i],
                                                all.hyper = all.hyper, range = range(m$x), debug = debug))
                        lines(xy, lwd = 1, col = "blue")
                    }

                }
            }
        }
    }

    if (plot.lincomb) {
        ##plot marginals for the derived lincombs
        fix = x$marginals.lincomb.derived
        labels.fix = names(x$marginals.lincomb.derived)
        nf = length(labels.fix)
        if (nf>0) {
            if (nf == 1 || single) {
                plot.layout = c(1, 1)
            } else if (nf == 2) {
                plot.layout = c(2, 1)
            } else {
                plot.layout = c(3, 3)
            }
            np = prod(plot.layout)
            ip = 0
            for(i in 1:nf) {
                if (ip%%np == 0) {
                    close.and.new.plot(...)
                    par(mfrow=c(plot.layout[1], plot.layout[2]))
                }
                ip = ip + 1

                if (!all(is.na(fix[[i]]))) {
                    ss = x$summary.lincomb.derived[i,]
                    sub=paste("Mean = ", round(ss[names(ss)=="mean"], 3)," SD = ", round(ss[names(ss)=="sd"], 3), sep="")
                    plot(inla.smarginal(fix[[i]]), type="l", main=paste("PostDens [", inla.nameunfix(labels.fix[i]),"] (derived)", sep=""),
                         sub=sub, xlab="", ylab="")
                }
            }
        }
    }

    if (plot.lincomb) {
        ## plot marginals for the lincombs
        if (inla.is.element("marginals.lincomb", x)) {
            fix = x$marginals.lincomb
            labels.fix = names(x$marginals.lincomb)
        } else {
            fix = NULL
            labels.fix = NULL
        }
        nf = length(labels.fix)
        if (nf>0) {
            if (nf == 1 || single) {
                plot.layout = c(1, 1)
            } else if (nf == 2) {
                plot.layout = c(2, 1)
            } else {
                plot.layout = c(3, 3)
            }
            np = prod(plot.layout)
            ip = 0
            for(i in 1:nf) {
                if (ip%%np == 0) {
                    close.and.new.plot(...)
                    par(mfrow=c(plot.layout[1], plot.layout[2]))
                }
                ip = ip + 1

                if (!all(is.na(fix[[i]]))) {
                    ss = x$summary.lincomb[i,]
                    sub=paste("Mean = ", round(ss[names(ss)=="mean"], 3)," SD = ", round(ss[names(ss)=="sd"], 3), sep="")
                    plot(inla.smarginal(fix[[i]]), type="l", main=paste("PostDens [", inla.nameunfix(labels.fix[i]),"]", sep=""),
                         sub=sub, xlab="", ylab="")
                }
            }
        }
    }

    if (plot.random.effects) {
        rand = x$summary.random
        labels.random = names(x$summary.random)
        nr = length(labels.random)
        if (nr>0) {
            for(i in 1:nr) {
                if (!all(is.na(rand[[i]]))) {
                    rr = rand[[i]]
                    dim1 = dim(rr)[1]
                    dim2 = dim(rr)[2]
                    tp = ifelse(labels.random[i] == "baseline.hazard", "s", "l")

                    r.n = x$size.random[[i]]$n
                    r.N = x$size.random[[i]]$N
                    r.Ntotal = x$size.random[[i]]$Ntotal
                    nrep = x$size.random[[i]]$nrep
                    ngroup = x$size.random[[i]]$ngroup
                    r.n.orig = r.n %/% ngroup
                    r.N.orig = r.N %/% ngroup

                    ## determine the plot-layout here
                    nrr = ngroup*nrep
                    if (nrr == 1 || single) {
                        plot.layout = c(1, 1)
                    } else if (nrr == 2) {
                        plot.layout = c(2, 1)
                    } else {
                        plot.layout = c(3, 3)
                    }
                    np = prod(plot.layout)
                    ip = 0

                    if (r.n > 1) {
                        for (r.rep in 1:nrep) {
                            for (r.group in 1:ngroup) {
                                for(ii in 1:inla.ifelse(r.N > r.n, r.N %/% r.n, 1)) {

                                    if (ip%%np == 0) {
                                        close.and.new.plot(...)
                                        par(mfrow=c(plot.layout[1], plot.layout[2]))
                                    }
                                    ip = ip + 1

                                    rep.txt = ""
                                    if (nrep > 1) {
                                        rep.txt = paste(rep.txt, " rep:", r.rep, sep="")
                                    }
                                    if (ngroup > 1) {
                                        rep.txt = paste(rep.txt, " group:", r.group, sep="")
                                    }
                                    if (r.N > r.n) {
                                        rep.txt = paste(rep.txt, " part:", ii, sep="")
                                    }

                                    idx = (r.rep-1)*r.N +  (r.group-1)*r.N.orig + (ii-1)*r.n.orig + (1:r.n.orig)

                                    ## if the dimension is > 1, then plot the means++
                                    ##
                                    xval = NULL
                                    if (tp == "s") ## baseline.hazard
                                    {
                                        xval = rr[, colnames(rr)=="ID"][idx]
                                        yval = rr[, colnames(rr)=="mean"][idx]
                                        plot(xval, yval,
                                             ylim=range(rr[, setdiff(colnames(rr), c("ID", "sd", "kld"))]),
                                             xlim=range(xval),
                                             axes=TRUE, ylab="", xlab="", type=tp, lwd=2)
                                        if (!is.null(x$.args$data$baseline.hazard.strata.coding)) {
                                            rep.txt = inla.paste(c(rep.txt, "[",
                                                    x$.args$data$baseline.hazard.strata.coding[r.rep], "]"), sep="")
                                        }
                                    } else {
                                        xval = rr[, colnames(rr)=="ID"][idx]
                                        yval = rr[, colnames(rr)=="mean"][idx]
                                        if (is.numeric(xval)) {
                                            plot(xval, yval,
                                                 ylim=range(rr[, setdiff(colnames(rr), c("ID", "sd", "kld"))]),
                                                 xlim=range(xval),
                                                 axes=FALSE, ylab="", xlab="", type=tp, lwd=2)
                                            axis(1)
                                            axis(2)
                                            box()
                                        } else {
                                            plot(as.factor(xval), yval,
                                                 ylim=range(rr[, setdiff(colnames(rr), c("ID", "sd", "kld"))]),
                                                 axes=TRUE, ylab="", xlab="", type=tp, lwd=2)
                                        }
                                    }

                                    lq = grep("quan", colnames(rr))
                                    main=inla.nameunfix(labels.random[i])

                                    if (length(lq)>0) {
                                        qq = rr[, lq]
                                        dq = dim(qq)[2]
                                        sub = paste("PostMean ")
                                        for(j in 1:dq) {
                                            yval = qq[, j][idx]
                                            ## this is the baseline.hazard case
                                            if (length(yval)+1 == length(xval)) {
                                                yval = c(yval, yval[ length(yval) ])
                                            }
                                            if (is.numeric(xval)) {
                                                points(xval, yval, type=tp, lty=2)
                                            } else {
                                                points(as.factor(xval), yval, pch=19)
                                            }
                                            sub = gsub("quant", "%", paste(sub, colnames(qq)[j]))
                                        }
                                        title(main=inla.nameunfix(main), sub=paste(inla.nameunfix(sub), rep.txt))
                                    }
                                    else
                                        title(main=inla.nameunfix(main), sub=paste("Posterior mean", rep.txt))
                                }
                            }
                        }
                    } else {
                        for (r.rep in 1:nrep) {
                            for(r.group in 1:ngroup) {
                                if (ip%%np == 0) {
                                    close.and.new.plot(...)
                                    par(mfrow=c(plot.layout[1], plot.layout[2]))
                                }
                                ip = ip + 1

                                rep.txt = ""
                                if (nrep > 1) {
                                    rep.txt = paste(rep.txt, ", replicate", r.rep)
                                }
                                if (ngroup > 1) {
                                    rep.txt = paste(rep.txt, ", group", r.group)
                                }

                                idx = (r.rep-1)*ngroup + r.group

                                ## if the dimension is 1, the plot the marginals
                                ##
                                if (!is.null(x$marginals.random[[i]])) {
                                    zz = x$marginals.random[[i]][[r.rep]]
                                    plot(inla.smarginal(zz), type="l",
                                         main=paste("PostDens [", inla.nameunfix(labels.random[i]),"]", " ", rep.txt, sep=""),
                                         xlab=inla.nameunfix(labels.random[i]), ylab="")
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    if (plot.hyperparameters) {
        hyper = x$marginals.hyperpar
        if (!is.null(hyper)) {
            nhyper = length(hyper)

            if (nhyper == 1 || single) {
                plot.layout = c(1, 1)
            } else if (nhyper == 2) {
                plot.layout = c(2, 1)
            } else {
                plot.layout = c(3, 3)
            }
            np = prod(plot.layout)

            ip = 0
            for(i in 1:nhyper) {
                if (ip%%np == 0) {
                    close.and.new.plot(...)
                    par(mfrow=c(plot.layout[1], plot.layout[2]))
                }
                ip = ip + 1

                hh = hyper[[i]]
                if (!is.null(hh)) {
                    label = inla.nameunfix(names(hyper)[i])
                    m = inla.smarginal(hh)
                    plot(m, type="l", ylab="", xlab="")
                    title(main=paste("PostDens [", label, "]", sep=""))

                    if (plot.prior) {
                        id = unlist(strsplit(attr(hyper[[i]], "hyperid"), "\\|"))
                        if (length(id) > 0) {
                            xy = (inla.get.prior.xy(section = tolower(id[2]), hyperid = id[1],
                                                    all.hyper = all.hyper, range = range(m$x), intern = FALSE,
                                                    debug = debug))
                            lines(xy, lwd = 1, col = "blue")
                        }
                    }
                }
            }
        }
    }

    if (plot.hyperparameters && intern) {
        hyper = x$internal.marginals.hyperpar
        if (!is.null(hyper)) {
            nhyper = length(hyper)

            if (nhyper == 1 || single) {
                plot.layout = c(1, 1)
            } else if (nhyper == 2) {
                plot.layout = c(2, 1)
            } else {
                plot.layout = c(3, 3)
            }
            np = prod(plot.layout)

            ip = 0
            for(i in 1:nhyper) {
                if (ip%%np == 0) {
                    close.and.new.plot(...)
                    par(mfrow=c(plot.layout[1], plot.layout[2]))
                }
                ip = ip + 1

                hh = hyper[[i]]
                if (!is.null(hh)) {
                    label = inla.nameunfix(names(hyper)[i])
                    m = inla.smarginal(hh)
                    plot(m, type="l", ylab="", xlab="")
                    title(main=paste("PostDens [", label, "]", sep=""))

                    if (plot.prior) {
                        id = unlist(strsplit(attr(hyper[[i]], "hyperid"), "\\|"))
                        if (length(id) > 0) {
                            xy = (inla.get.prior.xy(section = tolower(id[2]), hyperid = id[1], 
                                                    all.hyper = all.hyper, range = range(m$x), intern = TRUE,
                                                    debug = debug))
                            lines(xy, lwd = 1, col = "blue")
                        }
                    }
                }
            }
        }
    }

    if (plot.predictor) {

        ## linear predictor and fitted values
        n = x$size.linear.predictor$n
        if (!is.null(n)) {
            lp = x$summary.linear.predictor
            fv = x$summary.fitted.values
            ## remove the 'kld' column, we do not need it
            lp = lp[, names(lp) != "kld"]
            fv = fv[, names(fv) != "kld"]
            A = (x$size.linear.predictor$nrep == 2)

            if (A) {
                nm = dim(lp)[1] - n
            } else {
                nm = n
            }

            for(m in inla.ifelse(A, 1:2, 1)) {
                if (m == 1) {
                    idx = 1:nm
                    plot.idx = idx
                } else if (m == 2) {
                    idx = 1:n + nm
                    plot.idx = 1:n
                } else {
                    stop("This should not happen")
                }
                if (A) {
                    if (m == 1) {
                        msg = inla.ifelse(A, "(A*lin.pred)", "")
                    } else {
                        msg = "(orig.)"
                    }
                } else {
                    msg = ""
                }

                if (!is.null(lp)) {
                    close.and.new.plot(...)
                    if (!is.null(fv)) {
                        if (single) {
                            par(mfrow=c(1, 1))
                        } else {
                            par(mfrow=c(2, 1))
                        }
                    }
                    plot(plot.idx, lp[idx, colnames(lp)=="mean"], ylim=range(lp[idx, names(lp) != "sd"]), ylab="", xlab="Index", type="l", lwd=2, ...)
                    lq = grep("quan", colnames(lp))
                    if (length(lq)>0) {
                        qq = lp[, lq]
                        dq = dim(qq)[2]
                        sub = paste("Posterior mean together with ")
                        for(j in 1:dq) {
                            points(plot.idx, qq[idx, j], type="l", lty=2)
                            sub = paste(sub, colnames(qq)[j])
                        }
                        title(main=paste("Linear Predictor", msg), sub= inla.nameunfix(sub))
                    }
                    else
                        title(main=paste("Linear Predictor ", msg, inla.nameunfix(labels.random[i])), sub="Posterior mean")

                    if (single) {
                        close.and.new.plot()
                    }

                    if (!is.null(fv)) {
                        plot(plot.idx, fv[idx , colnames(fv)=="mean"], ylim=range(fv[idx, names(fv) != "sd"]), ylab="", xlab="Index", type="l", lwd=2, ...)
                        lq = grep("quan", colnames(fv))
                        if (length(lq)>0) {
                            qq = fv[, lq]
                            dq = dim(qq)[2]
                            sub = paste("Posterior mean together with ")
                            for(j in 1:dq) {
                                points(plot.idx, qq[idx, j], type="l", lty=2)
                                sub = paste(sub, colnames(qq)[j])
                            }
                            title(main=paste("Fitted values (inv.link(lin.pred))", msg), sub = inla.nameunfix(sub))
                        }
                        else
                            title(main=paste("Fitted values (inv.link(lin.pred))", msg, inla.nameunfix(labels.random[i])))
                    }
                }
            }
        }
    }

    if (plot.q && !is.null(x$Q)) {
        close.and.new.plot(...)
        if (single) {
            par(mfrow = c(1, 1))
        } else {
            par(mfrow = c(2, 2))
        }
        if (!is.null(x$Q$Q)) {
            plot(x$Q$Q, main = "The precision matrix", ...)
            nx = x$Q$Q@size[1L]
            lines(c(0, nx, nx, 0, 0), c(0, 0, nx, nx, 0), lwd=2)
        }
        if (!is.null(x$Q$Q.reorder)) {
            plot(x$Q$Q.reorder, main = "The reordered precision matrix", ...)
            nx = x$Q$Q.reorder@size[1L]
            lines(c(0, nx, nx, 0, 0), c(0, 0, nx, nx, 0), lwd=2)
        }
        if (!is.null(x$Q$L)) {
            plot(x$Q$L, main = "The Cholesky triangle", ...)
            nx = x$Q$L@size[1L]
            lines(c(0, nx, nx, 0, 0), c(0, 0, nx, nx, 0), lwd=2)
        }
    }

    if (plot.cpo && !is.null(x$cpo)) {
        if (!(is.null(x$cpo$pit) || length(x$cpo$pit) == 0L) || !(is.null(x$cpo$cpo) || length(x$cpo$cpo) == 0L)) {
            close.and.new.plot(...)
            if (single) {
                par(mfrow=c(1, 1))
            } else {
                par(mfrow=c(2, 2))
            }
        }

        ok = (!is.na(x$cpo$failure) & !is.na(x$cpo$pit) & !is.na(x$cpo$cpo))
        failure = x$cpo$failure[ok]
        pit = x$cpo$pit[ok]
        cpo = x$cpo$cpo[ok]

        if (!(is.null(pit) || length(pit) == 0L)) {
            ## if the observational model is discrete then do some
            ## adjustments: define the modified pit as Prob(y < y_obs)
            ## + 0.5*Prob(y = y_obs).
            ##
            n.fail = sum(failure > 0.0)
            if (!is.null(x$family) && inla.model.properties(x$family, "likelihood")$discrete) {
                m = "The (modified) PIT-values"
                p = pit - 0.5*cpo
            } else {
                m = "The PIT-values"
                p = pit
            }
            plot(p, main = paste(m, ", n.fail", n.fail, sep=""), ylab = "Probability", xlab = "index", ...)
            if (n.fail > 0) {
                points(p[ failure > 0 ], pch=20)
            }

            hist(p, main = paste(m, ", n.fail", n.fail, sep=""), xlab = "Probability",
                 n = max(20, min(round(length(pit)/10), 100)))
        }

        if (!(is.null(cpo) || length(cpo) == 0L)) {
            n.fail = sum(failure != 0.0)
            plot(cpo, main = paste("The CPO-values", ", nfail", n.fail, sep=""), ylab = "Probability", xlab = "index", ...)
            if (n.fail > 0) {
                points(cpo[ failure > 0 ], pch=20L)
            }
            hist(cpo, main = paste("Histogram of the CPO-values", ", nfail", n.fail, sep=""), xlab = "Probability",
                 n = max(20L, min(round(length(pit)/10L), 100L)))
        }
    }

    close.plot(...)
    return (invisible(figures))
}

inla.extract.prior = function(section = NULL, hyperid = NULL, all.hyper, debug=FALSE)
{
    str.trunc = function(..., max.len = 32)
    {
        str = paste(..., sep="", collapse = " ")
        str = substr(str, 1, min(max.len, nchar(str)))
        return(str)
    }

    output = function(..., warning=FALSE, force = FALSE)
    {
        msg = paste("*** inla.extract.prior: ", ..., sep="", collapse=" ")
        if (warning) warning(msg)
        if (debug || force) print(msg)
        return (invisible())
    }
    
    if (section == "fixed") {
        ## hyperid = name of fixed effect
        output("enter section [fixed]")
        output("searching for hyperid = ", hyperid)
        h = all.hyper$fixed
        found = FALSE
        for(i in seq_along(h)) {
            ## output("search for label = ", h[[i]]$label)
            if (inla.strcasecmp(h[[i]]$label, hyperid)) {
                found = TRUE
                break
            }
        }
        if (found) {
            prior = "gaussian"
            param = c(h[[i]]$prior.mean, h[[i]]$prior.prec)
            from.theta = function(x) x
            to.theta = function(x) x
        } else {
            ## try section 'linear' instead
            output("enter section [linear]")
            output("searching for hyperid = ", hyperid)
            h = all.hyper$linear
            found = FALSE
            for(i in seq_along(h)) {
                ## output("search for label = ", h[[i]]$label)
                if (inla.strcasecmp(h[[i]]$label, hyperid)) {
                    found = TRUE
                    break
                }
            }
            if (!found) {
                output(paste("cannot find hyperid '", hyperid,
                             "' in sections 'fixed' and 'linear'.",  sep=""), warning = TRUE)
                return (NA)
            }
            prior = "gaussian"
            param = c(h[[i]]$prior.mean, h[[i]]$prior.prec)
            from.theta = function(x) x
            to.theta = function(x) x
        }
    } else if (section == "predictor") {
        output("section predictor")
        h = all.hyper$predictor
        stopifnot(length(h) == 1)
        stopifnot(inla.strcasecmp(as.character(h[[1]]$theta$hyperid), hyperid))
        prior = h[[1]]$theta$prior
        param = h[[1]]$theta$param
        from.theta = h[[1]]$theta$from.theta
        to.theta = h[[1]]$theta$to.theta
    } else if (length(grep("^inla.data[0-9]+$", section)) > 0) {
        ## likelihood
        output("request for likelihood ", section, " with hyperid ",  hyperid)
        h = all.hyper$family
        found = FALSE
        for (idx.family in seq_along(h)) {
            if (inla.strcasecmp(section, h[[idx.family]]$hyperid)) {
                found = TRUE
                break
            }
        }
        if (!found) {
            output("likelihood ",  section, " with hyperid ", hyperid,  " is not found", warning=TRUE)
            return (NA)
        } else {
            output("likelihood ",  section, " with hyperid ", hyperid,  " is found with idx.family=", idx.family)
        }

        ## now we need to find the theta
        found = FALSE
        for (idx.theta in seq_along(h[[idx.family]]$hyper)) {
            if (inla.strcasecmp(hyperid, h[[idx.family]]$hyper[[idx.theta]]$hyperid)) {
                found = TRUE
                break
            }
        }
        if (found) {
            output("likelihood ",  section, " with hyperid ", hyperid,  ". theta is found with idx.theta=", idx.theta)
            prior = h[[idx.family]]$hyper[[idx.theta]]$prior
            param = h[[idx.family]]$hyper[[idx.theta]]$param
            from.theta = h[[idx.family]]$hyper[[idx.theta]]$from.theta
            to.theta = h[[idx.family]]$hyper[[idx.theta]]$to.theta
        } else {
            ## look in the link-section for theta
            found = FALSE
            for (idx.theta in seq_along(h[[idx.family]]$link$hyper)) {
                if (inla.strcasecmp(hyperid, h[[idx.family]]$link$hyper[[idx.theta]]$hyperid)) {
                    found = TRUE
                    break
                }
            }
            if (found) {
                output("likelihood ",  section, " with hyperid ", hyperid,  ". theta is found with idx.theta=", idx.theta)
                prior = h[[idx.family]]$link$hyper[[idx.theta]]$prior
                param = h[[idx.family]]$link$hyper[[idx.theta]]$param
                from.theta = h[[idx.family]]$link$hyper[[idx.theta]]$from.theta
                to.theta = h[[idx.family]]$link$hyper[[idx.theta]]$to.theta
            } else {
                output("likelihood ",  section, " with hyperid ", hyperid,  ". theta is not found",
                       warning = TRUE)
                return (NA)
            }
        }
    } else {
        ## random
        output("request for random ", section, " with hyperid ",  hyperid)
        h = all.hyper$random
        found = FALSE
        for (idx.random in seq_along(h)) {
            if (inla.strcasecmp(section, h[[idx.random]]$hyperid)) {
                found = TRUE
                break
            }
        }
        if (!found) {
            output("random ",  section, " with hyperid ", hyperid,  " is not found", warning = TRUE)
            return (NA)
        } else {
            output("random ",  section, " with hyperid ", hyperid,  " is found with idx.random=", idx.random)
        }

        ## now we need to find the theta
        found = FALSE
        for (idx.theta in seq_along(h[[idx.random]]$hyper)) {
            if (inla.strcasecmp(hyperid, h[[idx.random]]$hyper[[idx.theta]]$hyperid)) {
                found = TRUE
                break
            }
        }
        hyper = h[[idx.random]]$hyper
        if (!found) {
            for (idx.theta in seq_along(h[[idx.random]]$group.hyper)) {
                if (inla.strcasecmp(hyperid, h[[idx.random]]$group.hyper[[idx.theta]]$hyperid)) {
                    found = TRUE
                    break
                }
            }
            hyper = h[[idx.random]]$group.hyper
            if (!found) {
                output("random ",  section, " with hyperid ", hyperid,  ". theta is not found",
                       warning = TRUE)
                return (NA)
            } else {
                output("random ",  section, " with hyperid ", hyperid,  ". theta is found with (group) idx.theta=", idx.theta)
            }
        } else {
            output("random ",  section, " with hyperid ", hyperid,  ". theta is found with idx.theta=", idx.theta)
        }
        prior = hyper[[idx.theta]]$prior
        param = hyper[[idx.theta]]$param
        from.theta = hyper[[idx.theta]]$from.theta
        to.theta = hyper[[idx.theta]]$to.theta
    }

    output("prior ",  str.trunc(prior), " param ",  str.trunc(param))
    return (list(prior = prior, param = param, from.theta = from.theta, to.theta = to.theta))

}

inla.get.prior.xy = function(section = NULL, hyperid = NULL, all.hyper, debug=FALSE,
    len = 1000, range, intern = FALSE)
{
    str.trunc = function(..., max.len = 32)
    {
        str = paste(..., sep="", collapse = " ")
        str = substr(str, 1, min(max.len, nchar(str)))
        return(str)
    }
        
    output = function(..., stop=FALSE, force = FALSE)
    {
        msg = paste("*** inla.get.prior.xy: ", ..., sep="", collapse=" ")
        if (stop) stop(msg)
        if (debug || force) print(msg)
        return (invisible())
    }

    ## add priors here. the format is
    ##
    ## my.<NameOfPrior> = function(theta, param, log=FALSE)
    ##
    ## and should return the density (log=F) or log-density (log=T) for the internal parameter
    ## THETA (which is a vector), and the prior has parameters PARAM (which is the same ones as
    ## specified from within the R-interface and argument 'hyper'. the conversion to the prior
    ## density for the user-scale is done automatically.

    my.pc.gamma = function(theta, param, log=FALSE) 
    {
        ## see ?inla.pc.dgamma. this is the same prior, but for theta where x=exp(theta)
        x = exp(theta)
        ld = inla.pc.dgamma(x, lambda = param[1], log=TRUE) + theta
        return (if (log) ld else exp(ld))
    }

    my.pc.mgamma = function(theta, param, log=FALSE) 
    {
        ## see ?inla.pc.dgamma. this is the same prior, but for theta where x=exp(-theta)
        return (my.pc.gamma(-theta, param, log=log))
    }

    my.pc.gammacount = function(theta, param, log=FALSE) 
    {
        ## see ?inla.pc.dgammacount. this is the same prior, but for theta where x=exp(theta)
        x = exp(theta)
        ld = inla.pc.dgammacount(x, lambda = param[1], log=TRUE) + theta
        return (if (log) ld else exp(ld))
    }

    my.pc.cor0 = function(theta, param, log = FALSE)
    {
        e.theta = exp(theta)
        rho = 2.0 * e.theta/(1 + e.theta) - 1.0
        ld = (inla.pc.dcor0(rho, u = param[1], alpha = param[2], log=TRUE) +
              log(2) + theta - 2.0*log(1+e.theta))
        return (if (log) ld else exp(ld))
    }

    my.pc.rho0 = function(theta, param, log = FALSE)
    {
        return (my.pc.cor0(theta, param, log))
    }

    my.pc.cor1 = function(theta, param, log = FALSE)
    {
        e.theta = exp(theta)
        rho = 2.0 * e.theta/(1 + e.theta) - 1.0
        ld = (inla.pc.dcor1(rho, u = param[1], alpha = param[2], log=TRUE) +
              log(2) + theta - 2.0*log(1+e.theta))
        return (if (log) ld else exp(ld))
    }

    my.pc.rho1 = function(theta, param, log = FALSE)
    {
        return (my.pc.cor1(theta, param, log))
    }

    my.logitbeta = function(theta, param, log=FALSE)
    {
        e.theta = exp(theta)
        p = e.theta/(1.0 + e.theta)
        ld = (dbeta(p, shape1 = param[1], shape2 = param[2], log=TRUE) +
              theta - 2.0 * log(1 + e.theta))
        return (if (log) ld else exp(ld))
    }

    my.betacorrelation = function(theta, param, log=FALSE)
    {
        print(param)
        e.theta = exp(theta)
        p = e.theta/(1 + e.theta)
        ld = dbeta(p, shape1 = param[1], shape2 = param[2], log=TRUE) +
            theta - 2.0 * log(1.0 + e.theta)
        return (if (log) ld else exp(ld))
    }

    my.pcfgnh = function(theta, param, log=FALSE) 
    {
        ## we compute the PC-prior on the fly using these two packages. Its somewhat quick.
        inla.require("HKprocess")
        
        to.theta = inla.models()$latent$fgn$hyper$theta2$to.theta
        from.theta = inla.models()$latent$fgn$hyper$theta2$from.theta

        logdet.FGN = function(H, n) {
            ans = c()
            Hseq = H
            for(H in Hseq) {
                r = inla.acvfFGN(H, n-1)
                res = HKprocess::ltzc(r, rep(0, n))
                ans = c(ans,  as.numeric(res[2]))
            }
            return (ans)
        }
        d = function(H, n) return (sqrt(-logdet.FGN(H, n)))
        H.intern = seq(-10, 19, by = 0.1)
        dist = d(from.theta(H.intern), 100)
        log.d = log(dist/max(dist))
        res = cbind(H.intern = H.intern, log.d = log.d)
        log.d.spline = splinefun(res[, 1], res[, 2])
        d.spline = function(H.intern, ...) return (exp(log.d.spline(H.intern, ...)))
        d.spline.deriv = function(H.intern, ...) return (d.spline(H.intern, ...) * log.d.spline(H.intern, deriv=1L))

        U = param[1]
        alpha = param[2]
        lambda = -log(alpha)/d.spline(to.theta(U))
        ld = dexp(d.spline(theta), rate = lambda, log=TRUE) + log(abs(d.spline.deriv(theta)))
        ##print(cbind(theta = theta,  ld = ld))
        
        return (if (log) ld else exp(ld))
    }

    my.pc.prec = function(theta, param, log=FALSE)
    {
        ## sigma = exp(-theta/2)
        ## pi(theta) = lambda * exp(-lambda*sigma) * sigma/2
        ## param = c(U, alpha)
        lambda = -log(param[2])/param[1]
        sigma = exp(-theta/2.0)
        ld = log(lambda) + (-theta/2.0) - log(2.0) - lambda * sigma
        return (if (log) ld else exp(ld))
    }

    my.loggamma = function(theta, param, log=FALSE)
    {
        ## the log-of-a-gamma(a, b) density
        ld = dgamma(exp(theta), shape = param[1], rate = param[2], log=TRUE) + theta
        return (if (log) ld else exp(ld))
    }

    my.gaussian = function(theta, param, log=FALSE)
    {
        if (param[2] == 0) ## improper prior
            param[2] = 1.0/.Machine$double.eps
        ld = dnorm(theta, mean = param[1], sd = sqrt(1/param[2]), log=TRUE)
        return (if (log) ld else exp(ld))
    }
            
    my.normal = function(theta, param, log=FALSE)
    {
        return (my.gaussian(theta, param, log))
    }

    my.table = function(theta, param, log=FALSE)
    {
        fun = splinefun(param[, 1], param[, 2])
        ld = fun(theta)
        return (if (log) ld else exp(ld))
    }
        
    ## end of prior functions

    stopifnot(!missing(range))
    prior = inla.extract.prior(section, hyperid, all.hyper, debug)

    if (length(prior) == 1 && is.na(prior)) {
        return (list(x = NA, y = NA))
    } 

    if (length(grep("table:", prior$prior)) > 0) {
        tab = substr(prior$prior, nchar("table:")+1, nchar(prior$prior))
        xy = as.numeric(unlist(strsplit(tab, "[ \t\n\r]+")))
        xy = xy[!is.na(xy)]
        nxy = length(xy) %/% 2L
        xx = xy[1:nxy]
        yy = xy[1:nxy + nxy]
        xy = cbind(xx, yy)
        prior$param = xy
        prior$prior = "table"
    }
    if (length(grep("expression:", prior$prior)) > 0) {
        ## not available yet.
        return (list(x = NA, y = NA))
    }

    myp = paste("my.", prior$prior, sep="")
    if (!exists(myp) || !is.function(eval(parse(text = myp)))) {
        output("internal prior-function not found: ", myp, ", NEEDS TO BE IMPLEMENTED", force=TRUE)
        return (list(x = NA, y=NA))
    }
    output("prior: ", str.trunc(prior$prior), " param: ", str.trunc(prior$param))

    if (intern) {
        ## use a linear scale. 'x' is in the linear scale
        x = seq(range[1], range[2], len = len)
        y = do.call(myp, list(theta=x, param = prior$param))
    } else {
        ## 'x' is in the user-scale.
        range.theta = prior$to.theta(range)
        theta = seq(range.theta[1], range.theta[2], len = len)
        x = prior$from.theta(theta)
        ld = do.call(myp, args = list(theta = theta, param = prior$param, log=TRUE))
        fun = splinefun(x, theta)
        y = exp(ld + log(abs(fun(x, deriv=1))))
    }
    return (list(x = x, y = y))
}

`inla.all.hyper.postprocess` = function(all.hyper)
{
    ## postprocess all.hyper, by converting and replacing prior = 'mvnorm' into its p
    ## marginals. this is for the spde-models

    len.n = function(param, max.dim = 10000)
    {
        len = function(n) {
            return (n + n^2)
        }

        len.target = length(param)
        for(n in 1:max.dim) {
            if (len(n) == len.target) {
                return (n)
            }
        }
        stop(paste("length(param) is wrong:", len.target))
    }

    get.mvnorm.marginals = function(param)
    {
        n = len.n(param)
        mu = param[1:n]
        Q = matrix(param[-(1:n)], n, n)
        Sigma = solve((Q + t(Q))/2.)
        return (list(mean = mu, prec = 1/diag(Sigma)))
    }
    
    for (i in seq_along(all.hyper$random)) {
        for(j in seq_along(all.hyper$random[[i]]$hyper)) {

            if (all.hyper$random[[i]]$hyper[[j]]$prior == "mvnorm") {
                ## replace this one, and the p-following ones, with its marginals
                m = get.mvnorm.marginals(all.hyper$random[[i]]$hyper[[j]]$param)
                for(k in 1:length(m$mean)) {
                    kk = j + k - 1
                    all.hyper$random[[i]]$hyper[[kk]]$prior = "normal"
                    all.hyper$random[[i]]$hyper[[kk]]$param = c(m$mean[k], m$prec[k])
                }
            }
        }
    }

    return (all.hyper)
}
