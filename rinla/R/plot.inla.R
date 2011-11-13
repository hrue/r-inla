
`plot.inla` =
    function(x,
             plot.fixed.effects = TRUE,
             plot.lincomb = TRUE,
             plot.random.effects = TRUE,
             plot.hyperparameters = TRUE,
             plot.predictor = TRUE,
             plot.q = TRUE,
             plot.cpo = TRUE,
             single = FALSE, 
             ...)
{
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
            par(mfrow=c(plot.layout[1], plot.layout[2]))

            ip = 0
            for(i in 1:nf) {
                if (ip%%np == 0) {
                    dev.new()
                    par(mfrow=c(plot.layout[1], plot.layout[2]))
                }
                ip = ip + 1
                
                if (!all(is.na(fix[[i]]))) {
                    ss = x$summary.fixed[i,]
                    sub=paste("Mean = ", round(ss[names(ss)=="mean"], 3)," SD = ", round(ss[names(ss)=="sd"], 3), sep="")
                    plot(inla.smarginal(fix[[i]]), type="l", main=paste("PostDens [", inla.nameunfix(labels.fix[i]),"]", sep=""),
                         sub=sub, xlab="", ylab="", ...)
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
            par(mfrow=c(plot.layout[1], plot.layout[2]))

            ip = 0
            for(i in 1:nf) {
                if (ip%%np == 0) {
                    dev.new()
                    par(mfrow=c(plot.layout[1], plot.layout[2]))
                }
                ip = ip + 1
                
                if (!all(is.na(fix[[i]]))) {
                    ss = x$summary.lincomb.derived[i,]
                    sub=paste("Mean = ", round(ss[names(ss)=="mean"], 3)," SD = ", round(ss[names(ss)=="sd"], 3), sep="")
                    plot(inla.smarginal(fix[[i]]), type="l", main=paste("PostDens [", inla.nameunfix(labels.fix[i]),"] (derived)", sep=""),
                         sub=sub, xlab="", ylab="", ...)
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
            par(mfrow=c(plot.layout[1], plot.layout[2]))

            ip = 0
            for(i in 1:nf) {
                if (ip%%np == 0) {
                    dev.new()
                    par(mfrow=c(plot.layout[1], plot.layout[2]))
                }
                ip = ip + 1
                
                if (!all(is.na(fix[[i]]))) {
                    ss = x$summary.lincomb[i,]
                    sub=paste("Mean = ", round(ss[names(ss)=="mean"], 3)," SD = ", round(ss[names(ss)=="sd"], 3), sep="")
                    plot(inla.smarginal(fix[[i]]), type="l", main=paste("PostDens [", inla.nameunfix(labels.fix[i]),"]", sep=""),
                         sub=sub, xlab="", ylab="", ...)
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
                                        dev.new()
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
                                        xval= x$.internal$baseline.hazard.cutpoints
                                        yval = rr[, colnames(rr)=="mean"][idx]
                                        yval = c(yval, yval[length(yval)])
                                        plot(xval, yval,
                                             ylim=range(rr[, setdiff(colnames(rr), c("ID", "sd", "kld"))]),
                                             xlim=range(xval),
                                             axes=TRUE, ylab="", xlab="", type=tp, lwd=2, ...)
                                        if (!is.null(x$.internal$baseline.hazard.strata.coding)) {
                                            rep.txt = inla.paste(c(rep.txt, "[",
                                                    x$.internal$baseline.hazard.strata.coding[r.rep], "]"), sep="")
                                        }
                                    } else {
                                        xval = rr[, colnames(rr)=="ID"][idx]
                                        yval = rr[, colnames(rr)=="mean"][idx]
                                        if (is.numeric(xval)) {
                                            plot(xval, yval,
                                                 ylim=range(rr[, setdiff(colnames(rr), c("ID", "sd", "kld"))]),
                                                 xlim=range(xval),
                                                 axes=FALSE, ylab="", xlab="", type=tp, lwd=2, ...)
                                            axis(1)
                                            axis(2)
                                            box()
                                        } else {
                                            plot(as.factor(xval), yval,
                                                 ylim=range(rr[, setdiff(colnames(rr), c("ID", "sd", "kld"))]),
                                                 axes=TRUE, ylab="", xlab="", type=tp, lwd=2, ...)
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
                                    dev.new()
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
                                         xlab=inla.nameunfix(labels.random[i]), ylab="", ...)
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
                    dev.new()
                    par(mfrow=c(plot.layout[1], plot.layout[2]))
                }
                ip = ip + 1

                hh = hyper[[i]]
                if (!is.null(hh)) {
                    label = inla.nameunfix(names(hyper)[i])
                    plot(inla.smarginal(hh), type="l", ylab="", xlab="", ...)
                    title(main=paste("PostDens [", label, "]", sep=""))
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
                    dev.new()
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
    if (plot.q && !is.null(x$Q.matrix)) {
        dev.new()
        if (single) {
            par(mfrow = c(1, 1))
        } else {
            par(mfrow = c(2, 2))
        }
        if (!is.null(x$Q.matrix)) {
            if (class(x$Q.matrix) == "matrix") {
                image(x$Q.matrix, main = "The precision matrix")
            } else {
                ## old style format
                xx = x$Q.matrix
                xx@grey = 1-xx@grey
                plot(xx, main = "The precision matrix",
                     axes = TRUE, asp = 1, xlim = c(1, xx@size[1]), ylim = c(1, xx@size[2]), ...)
                box()
                rm(xx)
            }
        }
        if (!is.null(x$Q.matrix.reorder)) {
            if (class(x$Q.matrix.reorder) == "matrix") {
                image(x$Q.matrix.reorder, main = "The reordered precision matrix")
            } else {
                ## old style format
                xx = x$Q.matrix.reorder
                xx@grey = 1-xx@grey
                plot(xx, main = "The reordered precision matrix",
                     axes = TRUE, asp = 1, xlim = c(1, xx@size[1]), ylim = c(1, xx@size[2]), ...)
                box()
                rm(xx)
            }
        }
        if (!is.null(x$L)) {
            if (class(x$L) == "matrix") {
                image(x$L, main = "The Cholesky triangle")
            } else {
                ## old style format
                xx = x$L
                xx@grey = 1-xx@grey
                plot(xx, main = "The Cholesky triangle",
                     axes = TRUE, asp = 1, xlim = c(1, xx@size[1]), ylim = c(1, xx@size[2]), ...)
                box()
                rm(xx)
            }
        }
    }
    if (plot.cpo) {
        if (!is.null(x$pit) || !is.null(x$cpo)) {
            dev.new()
            if (single) {
                par(mfrow=c(1, 1))
            } else {
                par(mfrow=c(2, 2))
            }
        }
        if (!is.null(x$pit)) {
            ## if the observational model is discrete then do some
            ## adjustments: define the modified pit as Prob(y < y_obs)
            ## + 0.5*Prob(y = y_obs).
            ##

            oks = !is.nan(x$failure)
            x$failure = x$failure[oks]
            x$cpo = x$cpo[oks]
            x$pit = x$pit[oks]

            n.fail = sum(x$failure > 0.0)
            if (!is.null(x$family) && inla.model.properties(x$family, "likelihood")$discrete) {
                m = "The (modified) PIT-values"
                pit = x$pit - 0.5*x$cpo
            } else {
                m = "The PIT-values"
                pit = x$pit
            }
            plot(pit, main = paste(m, ", n.fail", n.fail, sep=""), ylab = "Probability", xlab = "index", ...)
            if (n.fail > 0) {
                points(pit[ x$failure > 0 ], pch=20)
            }

            hist(pit, main = paste(m, ", n.fail", n.fail, sep=""), xlab = "Probability",
                 n = max(20, min(round(length(x$pit)/10), 100)))
        }
        if (!is.null(x$cpo)) {
            n.fail = sum(x$failure != 0.0)
            plot(x$cpo, main = paste("The CPO-values", ", n.fail", n.fail, sep=""), ylab = "Probability", xlab = "index", ...)
            if (n.fail > 0) {
                points(x$cpo[ x$failure > 0 ], pch=20)
            }
            hist(x$cpo, main = paste("Histogram of the CPO-values", ", n.fail", n.fail, sep=""), xlab = "Probability",
                 n = max(20, min(round(length(x$pit)/10), 100)))
        }
    }
}

