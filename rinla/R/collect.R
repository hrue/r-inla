
### AUXILIARY FUNCTIONS TO COLLECT RESULTS FROM THE DIRECTORY
### "inla.model/results.files" or any directory where the results from a
### inla run are stored

## temporary only...
inla.internal.experimental.mode = FALSE


`inla.collect.misc` = function(dir, debug = FALSE)
{
    d = paste(dir,"/misc", sep="")
    d.info = file.info(d)$isdir

    if (debug)
        print(paste("collect misc from", d))

    if (is.na(d.info) || (d.info == FALSE))
        return (NULL)
    
    fnm = paste(d, "/theta-tags", sep="")
    if (file.exists(fnm)) {
        tags = readLines(fnm)
    } else {
        tags = NULL
    }
        
    fnm = paste(d, "/covmat-hyper-internal.dat", sep="")
    if (file.exists(fnm)) {
        siz = inla.read.binary.file(fnm)
        n = siz[1]
        stopifnot(length(siz) == n^2+1)

        cov.intern = matrix(siz[-1], n, n)
        dd = diag(cov.intern)
        s = matrix(0, n, n)
        diag(s) = 1/sqrt(dd)
        cor.intern = s %*% cov.intern %*% s
        diag(cor.intern) = 1
    } else {
        cov.intern = NULL
        cor.intern = NULL
    }

    fnm = paste(d, "/reordering.dat", sep="")
    if (file.exists(fnm)) {
        r = as.integer(inla.read.binary.file(fnm))
    } else {
        r = NULL
    }

    fnm = paste(d, "/stdev_corr_pos.dat", sep="")
    if (file.exists(fnm)) {
        stdev.corr.positive = as.numeric(inla.read.fmesher.file(fnm))
    } else {
        stdev.corr.positive = NULL
    }

    fnm = paste(d, "/stdev_corr_neg.dat", sep="")
    if (file.exists(fnm)) {
        stdev.corr.negative = as.numeric(inla.read.fmesher.file(fnm))
    } else {
        stdev.corr.negative = NULL
    }
    
    if (debug)
        print(paste("collect misc from", d, "...done"))

    return (list(cov.intern = cov.intern, cor.intern = cor.intern, reordering = r, theta.tags = tags,
                 stdev.corr.negative = stdev.corr.negative, stdev.corr.positive = stdev.corr.positive))
}

`inla.collect.size` = function(dir, debug = FALSE)
{
    fnm = paste(dir, "/size.dat", sep="")
    siz = inla.read.binary.file(fnm)
    if (length(siz) != 5)
        stop(paste("length of siz is not 5: fnm=", fnm))

    if (is.na(siz[1]) || siz[1] < 0) stop("siz[1] = NA")
    if (is.na(siz[2]) || siz[2] <= 0) siz[2] = siz[1]
    if (is.na(siz[3]) || siz[3] <= 0) siz[3] = siz[2]
    if (is.na(siz[4]) || siz[4] <= 0) siz[4] = 1
    if (is.na(siz[5]) || siz[5] <= 0) siz[5] = 1

    return (list(n=siz[1], N = siz[2], Ntotal = siz[3], ngroup = siz[4], nrep=siz[5]))
}

`inla.collect.fixed` = function(results.dir, debug = FALSE)
{
    alldir=dir(results.dir)
    if (debug)
        print("collect fixed effects")
    
    ## read FIXED EFFECTS
    fix = alldir[grep("^fixed.effect", alldir)]
    fix = c(fix, alldir[grep("^intercept$", alldir)])
    n.fix = length(fix)

    ##read the names of the fixed effects
    if (n.fix>0) {
        names.fixed = inla.trim(character(n.fix))
        for(i in 1:n.fix) {
            tag = paste(results.dir, .Platform$file.sep, fix[i], .Platform$file.sep,"TAG", sep="")
            if (!file.exists(tag))
                names.fixed[i] = "missing NAME"
            else
                names.fixed[i] = inla.namefix(readLines(tag, n=1))
        }
        ##read summary the fixed effects
        if (debug)
            print(names.fixed)
        
        summary.fixed = numeric()
        marginals.fixed = list()
        for(i in 1:n.fix) {
            file =  paste(results.dir, .Platform$file.sep, fix[i], sep="")
            dir.fix = dir(file)
            if (length(dir.fix) > 3) {
                sum = inla.read.binary.file(paste(file, .Platform$file.sep,"summary.dat", sep=""))[-1]
                col.nam = c("mean","sd")
            
                ##read quantiles if existing
                if (length(grep("^quantiles.dat$", dir.fix))>0) {
                    qq = inla.interpret.vector(inla.read.binary.file(paste(file, .Platform$file.sep, "quantiles.dat", sep="")),
                            debug=debug)
                    sum = c(sum, qq[, 2])
                    col.nam = c(col.nam, paste(as.character(qq[, 1]),"quant", sep=""))
                }

                ##read quantiles if existing
                if (length(grep("^cdf.dat$", dir.fix))>0) {
                    qq = inla.interpret.vector(inla.read.binary.file(paste(file, .Platform$file.sep, "cdf.dat", sep="")),
                            debug=debug)
                    sum = c(sum, qq[, 2])
                    col.nam = c(col.nam, paste(as.character(qq[, 1]),"cdf", sep=""))
                }
            
                ##read also kld distance
                kld.fixed = inla.read.binary.file(paste(file, .Platform$file.sep,"symmetric-kld.dat", sep=""))[-1]
                sum = c(sum, kld.fixed)
                col.nam = c(col.nam, "kld")
                summary.fixed = rbind(summary.fixed, sum)

                ##read the marginals
                xx = inla.interpret.vector(inla.read.binary.file(paste(file, .Platform$file.sep,"marginal-densities.dat", sep="")),
                        debug=debug)
                if (is.null(xx))
                    xx = cbind(c(NA, NA, NA), c(NA, NA, NA))
                colnames(xx) = c("x", "y")
                marginals.fixed[[i]] = xx

                if (inla.internal.experimental.mode) {
                    class(marginals.fixed[[i]]) = "inla.marginal"
                    attr(marginals.fixed[[i]], "inla.tag") = paste("marginal fixed", names.fixed[i])
                }
            } else {
                col.nam = c("mean", "sd", "kld")
                summary.fixed = rbind(summary.fixed, c(NA, NA, NA))
                xx = cbind(c(NA, NA, NA), c(NA, NA, NA))
                colnames(xx) = c("x", "y")
                marginals.fixed[[i]] = xx

                if (inla.internal.experimental.mode) {
                    class(marginals.fixed[[i]]) = "inla.marginal"
                    attr(marginals.fixed[[i]], "inla.tag") = paste("marginal fixed", names.fixed[i])
                }
            }
        }    
        rownames(summary.fixed) = inla.namefix(names.fixed)
        colnames(summary.fixed) = inla.namefix(col.nam)
        if (length(marginals.fixed) > 0) {
            names(marginals.fixed) = inla.namefix(names.fixed)
        }
    }
    else {
        if (debug)
            print("No fixed effects")
        names.fixed=NULL
        summary.fixed=NULL
        marginals.fixed=NULL
    }
    
    if (inla.internal.experimental.mode) {
        class(marginals.fixed) = "inla.marginals"
        attr(marginals.fixed,  "inla.tag", "marginals fixed")
    }

    ret = list(names.fixed=names.fixed, summary.fixed=summary.fixed, marginals.fixed=marginals.fixed)
    return(ret)
}

`inla.collect.lincomb` =
    function(results.dir,
             debug = FALSE,
             derived = TRUE)
{
    ## rewrite from collect.random
    alldir = dir(results.dir)
    if (derived) {
        lincomb = alldir[grep("^lincomb.*derived[.]all", alldir)]
    } else {
        lincomb1 = alldir[grep("^lincomb.*derived[.]all", alldir)]
        lincomb2 = alldir[grep("^lincomb", alldir)]
        lincomb = setdiff(lincomb2, lincomb1)
        if (debug)
            print(paste("lincomb",  lincomb))
    }        
    n.lincomb = length(lincomb)
    if (debug)
        print("collect lincombs")

    ##read the names and model of the lincomb effects
    if (n.lincomb > 0) {
        names.lincomb = inla.namefix(character(n.lincomb))
        model.lincomb = inla.trim(character(n.lincomb))

        summary.lincomb = list()
        marginals.lincomb = list()
        size.lincomb = list()
        
        for(i in 1:n.lincomb) {
            if (debug)
                print(paste("read lincomb ", i , " of ", n.lincomb))

            ##read the summary
            file= paste(results.dir, .Platform$file.sep, lincomb[i], sep="")
            dir.lincomb = dir(file)

            if (debug)
                print(paste("read from dir ",  file))

            if (length(dir.lincomb) > 4) {
                dd = matrix(inla.read.binary.file(file=paste(file, .Platform$file.sep,"summary.dat", sep="")),
                        ncol=3, byrow=TRUE)
                col.nam = c("ID","mean","sd")
            
                ##read quantiles if existing
                if (debug)
                    cat("...quantiles.dat if any\n")
                if (length(grep("^quantiles.dat$", dir.lincomb))==1) {
                    xx = inla.interpret.vector(inla.read.binary.file(paste(file, .Platform$file.sep,
                            "quantiles.dat", sep="")), debug=debug)
                    len = dim(xx)[2]
                    qq = xx[, seq(2, len, by=2), drop=FALSE]
                    col.nam = c(col.nam, paste(as.character(xx[, 1]),"quant", sep=""))
                    dd = cbind(dd, t(qq))
                }

                ##read cdf if existing
                if (debug)
                    cat("...cdf.dat if any\n")
                if (length(grep("^cdf.dat$", dir.lincomb))==1) {
                    xx = inla.interpret.vector(inla.read.binary.file(paste(file, .Platform$file.sep,"cdf.dat", sep="")),
                            debug=debug)
                    len = dim(xx)[2]
                    qq = xx[, seq(2, len, by=2), drop=FALSE]
                    col.nam = c(col.nam, paste(as.character(xx[, 1])," cdf", sep=""))
                    dd = cbind(dd, t(qq))
                }

                if (debug)
                    cat("...NAMES if any\n")
                if (length(grep("^NAMES$", dir.lincomb))==1) {
                    row.names = readLines(paste(file, .Platform$file.sep,"NAMES", sep=""))
                    ## remove the prefix 'lincomb.' as we do not need it in the names.
                    row.names = sapply(row.names, function(x) gsub("^lincomb[.]", "", x))
                    names(row.names) = NULL
                } else {
                    row.names = NULL
                }
                
                ##read kld
                if (debug)
                    cat("...kld\n")
                kld1 = matrix(inla.read.binary.file(file=paste(file, .Platform$file.sep,"symmetric-kld.dat", sep="")),
                        ncol=2, byrow=TRUE)
                qq = kld1[, 2, drop=FALSE]
                dd = cbind(dd, qq)
                if (debug)
                    cat("...kld done\n")
            
                col.nam = c(col.nam, "kld")
                colnames(dd) = inla.namefix(col.nam)
                summary.lincomb[[i]] = as.data.frame(dd)
                if (!is.null(row.names)) {
                    rownames(summary.lincomb[[i]]) = row.names
                }

                xx = inla.read.binary.file(paste(file, .Platform$file.sep,"marginal-densities.dat", sep=""))
                rr = inla.interpret.vector.list(xx, debug=debug)
                rm(xx)
                if (!is.null(rr)) {
                    nd = length(rr)
                    names(rr) = inla.namefix(paste("index.", as.character(1:nd), sep=""))
                    for(j in 1:nd) {
                        colnames(rr[[j]]) = inla.namefix(c("x", "y"))
                        if (inla.internal.experimental.mode) {
                            class(rr[[j]]) = "inla.marginal"
                            if (derived) {
                                attr(rr[[j]], "inla.tag") = paste("marginal lincomb derived", names(rr)[j])
                            } else {
                                attr(rr[[j]], "inla.tag") = paste("marginal lincomb", names(rr)[j])
                            }
                        }
                    }
                }
                marginals.lincomb[[i]] = rr
                
                if (!is.null(row.names) && (length(marginals.lincomb)>0)) {
                    names(marginals.lincomb[[i]]) = row.names
                }
            } else {
                N.file = paste(file, .Platform$file.sep,"N", sep="")
                if (!file.exists(N.file)) {
                    N = 0
                } else {
                    N = scan(file=N.file, what = numeric(0), quiet=TRUE)
                }
                summary.lincomb[[i]] = data.frame("mean" = rep(NA, N), "sd" = rep(NA, N), "kld" = rep(NA, N))
                marginals.lincomb = NULL
            }
            size.lincomb[[i]] = inla.collect.size(file)

            if (inla.internal.experimental.mode) {
                if (!is.null(marginals.lincomb)) {
                    class(marginals.lincomb[[i]]) = "inla.marginals"
                    if (derived) {
                        attr(marginals.lincomb[[i]], "inla.tag") = "marginal lincomb derived"
                    } else {
                        attr(marginals.lincomb[[i]], "inla.tag") = "marginal lincomb"
                    }                    
                }
            }
        }
        names(summary.lincomb) = inla.namefix(names.lincomb)
        if (!is.null(marginals.lincomb) && (length(marginals.lincomb) > 0))
            names(marginals.lincomb) = inla.namefix(names.lincomb)
    } else {
        if (debug)
            cat("No lincomb effets\n")
        summary.lincomb=NULL
        marginals.lincomb=NULL
        size.lincomb = NULL
    }

    if (derived) {
        res = list(summary.lincomb.derived = summary.lincomb[[1]],
                marginals.lincomb.derived = inla.ifelse(length(marginals.lincomb) > 0, marginals.lincomb[[1]], NULL), 
                size.lincomb.derived = size.lincomb[[1]])
    } else {
        res = list(summary.lincomb = summary.lincomb[[1]],
                marginals.lincomb = inla.ifelse(length(marginals.lincomb)>0, marginals.lincomb[[1]], NULL), 
                size.lincomb = size.lincomb[[1]])
    }
    return(res)
}

`inla.collect.cpo` =
    function(results.dir,
             debug = FALSE)
{
    alldir = dir(results.dir)
    if (length(grep("^cpo$", alldir))==1) {
        if (debug)
            cat(paste("collect cpo\n", sep=""))
      
        xx = inla.read.binary.file(file=paste(results.dir, .Platform$file.sep,"cpo", .Platform$file.sep,"cpo.dat", sep=""))
        n = xx[1]
        xx = xx[-1]
        len = length(xx)
        cpo.res=numeric(n)
        cpo.res[1:n] = NA
        cpo.res[xx[seq(1, len, by=2)] +1] = xx[seq(2, len, by=2)]
        
        xx = inla.read.binary.file(file=paste(results.dir, .Platform$file.sep,"cpo", .Platform$file.sep,"pit.dat", sep=""))
        n = xx[1]
        xx = xx[-1]
        len = length(xx)
        pit.res = numeric(n)
        pit.res[1:n] = NA
        pit.res[xx[seq(1, len, by=2)] +1] = xx[seq(2, len, by=2)]

        fnm=paste(results.dir, .Platform$file.sep,"cpo", .Platform$file.sep,"failure.dat", sep="")
        if (file.exists(fnm)) {
            xx = inla.read.binary.file(fnm)
            n = xx[1]
            xx = xx[-1]
            len = length(xx)
            failure.res = numeric(n)
            failure.res[1:n] = NA
            failure.res[xx[seq(1, len, by=2)] +1] = xx[seq(2, len, by=2)]
        }
        else
            failure.res = NULL
        rm(xx)
    }
    else {
        cpo.res = NULL
        pit.res = NULL
        failure.res = NULL
    }
    return(list(cpo=cpo.res, pit=pit.res, failure=failure.res))
}

`inla.collect.dic` =
    function(results.dir,
             debug = FALSE)
{
    alldir = dir(results.dir)
    ## get dic (if exists)
    if (length(grep("^dic$", alldir))==1) {
        if (debug)
            cat(paste("collect dic\n", sep=""))
        file=paste(results.dir, .Platform$file.sep,"dic", .Platform$file.sep,"dic.dat", sep="")
        dic.values = inla.read.binary.file(file)

        file=paste(results.dir, .Platform$file.sep,"dic", .Platform$file.sep,"deviance_e.dat", sep="")
        if (inla.is.fmesher.file(file)) {
            dev.e = c(inla.read.fmesher.file(file))
        } else {
            dev.e = NULL
        }

        file=paste(results.dir, .Platform$file.sep,"dic", .Platform$file.sep,"e_deviance.dat", sep="")
        if (inla.is.fmesher.file(file)) {
            e.dev = c(inla.read.fmesher.file(file))
        } else {
            e.dev = NULL
        }
        dic = list(
                "dic" = dic.values[4],
                "p.eff"= dic.values[3],
                "local.dic" = 2.0*e.dev - dev.e, 
                "local.p.eff" = e.dev - dev.e, 
                "mean.deviance" = dic.values[1],
                "deviance.mean" = dic.values[2])
    } else {
        dic = NULL
    }
    return(dic)     
}

`inla.collect.q` =
    function(results.dir,
             image.dim = 256,
             debug = FALSE)
{
    alldir = dir(results.dir)
    if (length(grep("^Q$", alldir))==1) {
        w = getOption("warn")
        options(warn = -1)
        pixm = require("pixmap", quietly = TRUE)
        options(warn = w)
        
        if (debug)
            cat(paste("collect q\n", sep=""))
        
        file=paste(results.dir, .Platform$file.sep,"Q/precision-matrix.pbm", sep="")
        if (file.exists(file) && pixm)
            Q.matrix = inla.image.reduce(read.pnm(file), image.dim = image.dim)
        else
            Q.matrix = NULL
        
        file=paste(results.dir, .Platform$file.sep,"Q/precision-matrix-reordered.pbm", sep="")
        if (file.exists(file) && pixm)
            Q.matrix.reorder = inla.image.reduce(read.pnm(file), image.dim = image.dim)
        else
            Q.matrix.reorder = NULL
        
        file=paste(results.dir, .Platform$file.sep,"Q/precision-matrix_L.pbm", sep="")
        if (file.exists(file) && pixm)
            L = inla.image.reduce(read.pnm(file), image.dim = image.dim)
        else
            L = NULL

        if (is.null(Q.matrix) && is.null(Q.matrix.reorder) && is.null(L))
            q = NULL
        else
            q = list(Q.matrix = Q.matrix, Q.matrix.reorder = Q.matrix.reorder, L = L)
    }
    return(q)     
}

`inla.collect.hyperpar` =
    function(results.dir,
             debug=FALSE)
{
    alldir = dir(results.dir)
    all.hyper = alldir[grep("^hyperparameter", alldir)]
    hyper = all.hyper[grep("user-scale$", all.hyper)]
    n.hyper = length(hyper)
    if (n.hyper>0) {
        ## get names for hyperpar
        names.hyper = character(n.hyper)
        for(i in 1:n.hyper) {
            tag = paste(results.dir, .Platform$file.sep, hyper[i], .Platform$file.sep,"TAG", sep="")
            if (!file.exists(tag))
                names.hyper[i] = "missing NAME"
            else
                names.hyper[i] = inla.namefix(readLines(tag, n=1))
        }

        ## get summary and marginals
        summary.hyper = numeric()
        marginal.hyper = list()
        
        for(i in 1:n.hyper) {
            dir.hyper =  paste(results.dir, .Platform$file.sep, hyper[i], sep="")
            file = paste(dir.hyper, .Platform$file.sep,"summary.dat", sep="")
            dd = inla.read.binary.file(file)[-1]
            sum = dd
            col.nam = c("mean","sd")
            if (length(grep("^quantiles.dat$", dir(dir.hyper)))>0) {
                qq = inla.interpret.vector(inla.read.binary.file(paste(dir.hyper, .Platform$file.sep, "quantiles.dat", sep="")),
                        debug=debug)
                sum = c(sum, qq[, 2])
                col.nam = c(col.nam, paste(as.character(qq[, 1]),"quant", sep=""))
            }
            if (length(grep("^cdf.dat$", dir(dir.hyper)))>0) {
                qq = inla.interpret.vector(inla.read.binary.file(paste(dir.hyper, .Platform$file.sep, "cdf.dat", sep="")),
                        debug=debug)
                sum = c(sum, qq[, 2])
                col.nam = c(col.nam, paste(as.character(qq[, 1]),"cdf", sep=""))
            }
            summary.hyper = rbind(summary.hyper, sum)
            file =paste(results.dir, .Platform$file.sep, hyper[i], .Platform$file.sep,"marginal-densities.dat", sep="")
            xx = inla.read.binary.file(file)
            marg1 = inla.interpret.vector(xx, debug=debug)
            rm(xx)
            if (!is.null(marg1)) {
                colnames(marg1) = c("x","y")
            }

            if (inla.internal.experimental.mode) {
                class(marg1) = "inla.marginal"
                attr(marg1, "inla.tag") = paste("marginal hyper", names.hyper[i])
            }
            
            marginal.hyper = c(marginal.hyper, list(marg1))
        }
        names(marginal.hyper) = inla.namefix(names.hyper)
        rownames(summary.hyper) = inla.namefix(names.hyper)
        colnames(summary.hyper) = inla.namefix(col.nam)
    } else {
        marginal.hyper=NULL
        summary.hyper=NULL
    }

    if (inla.internal.experimental.mode) {
        if (!is.null(marginal.hyper)) {
            class(marginal.hyper) = "inla.marginals"
            attr(marginal.hyper, "inla.tag") = "marginal hyper"
        }
    }
    
    ## collect also the hyperparameters in the internal scale
    all.hyper = alldir[grep("^hyperparameter", alldir)]
    hyper = all.hyper[-grep("user-scale$", all.hyper)]
    n.hyper = length(hyper)
    if (n.hyper>0) {
        ## get names for hyperpar
        names.hyper = character(n.hyper)
        for(i in 1:n.hyper) {
            tag = paste(results.dir, .Platform$file.sep, hyper[i], .Platform$file.sep,"TAG", sep="")
            if (!file.exists(tag))
                names.hyper[i] = "missing NAME"
            else
                names.hyper[i] = inla.namefix(readLines(tag, n=1))
        }

        ## get summary and marginals
        internal.summary.hyper = numeric()
        internal.marginal.hyper = list()
        
        for(i in 1:n.hyper) {
            dir.hyper =  paste(results.dir, .Platform$file.sep, hyper[i], sep="")
            file = paste(dir.hyper, .Platform$file.sep,"summary.dat", sep="")
            dd = inla.read.binary.file(file)[-1]
            sum = dd
            col.nam = c("mean","sd")
            if (length(grep("^quantiles.dat$", dir(dir.hyper)))>0) {
                qq = inla.interpret.vector(inla.read.binary.file(paste(dir.hyper, .Platform$file.sep, "quantiles.dat", sep="")),
                        debug=debug)
                sum = c(sum, qq[, 2])
                col.nam = c(col.nam, paste(as.character(qq[, 1]),"quant", sep=""))
            }
            if (length(grep("^cdf.dat$", dir(dir.hyper)))>0) {
                qq = inla.interpret.vector(inla.read.binary.file(paste(dir.hyper, .Platform$file.sep, "cdf.dat", sep="")),
                        debug=debug)
                sum = c(sum, qq[, 2])
                col.nam = c(col.nam, paste(as.character(qq[, 1]),"cdf", sep=""))
            }
            internal.summary.hyper = rbind(internal.summary.hyper, sum)
            file =paste(results.dir, .Platform$file.sep, hyper[i], .Platform$file.sep,"marginal-densities.dat", sep="")
            xx = inla.read.binary.file(file)
            marg1 = inla.interpret.vector(xx, debug=debug)
            rm(xx)
            if (!is.null(marg1))
                colnames(marg1) = c("x","y")

            if (inla.internal.experimental.mode) {
                class(marg1) = "inla.marginal"
                attr(marg1, "inla.tag") = paste("marginal hyper internal", names.hyper[i])
            }
            
            internal.marginal.hyper = c(internal.marginal.hyper, list(marg1))
        }
        names(internal.marginal.hyper) = inla.namefix(names.hyper)
        rownames(internal.summary.hyper) = inla.namefix(names.hyper)
        colnames(internal.summary.hyper) = inla.namefix(col.nam)
    }
    else {
        internal.summary.hyper=NULL
        internal.marginal.hyper=NULL
    }
    
    if (inla.internal.experimental.mode) {
        if (!is.null(internal.marginal.hyper)) {
            class(internal.marginal.hyper) = "inla.marginals"
            attr(internal.marginal.hyper, "inla.tag") = "marginal hyper internal"
        }
    }
    
    ret=list(summary.hyperpar=summary.hyper,
        marginals.hyperpar=marginal.hyper,
        internal.summary.hyperpar = internal.summary.hyper,
        internal.marginals.hyperpar = internal.marginal.hyper)
    return(ret)
}

`inla.collect.mlik` =
    function(results.dir,
             debug = FALSE)
{
    alldir = dir(results.dir)
    if (length(grep("^marginal-likelihood$", alldir))==1) {
        if (debug)
            cat(paste("collect mlik\n", sep=""))
        file=paste(results.dir, .Platform$file.sep,"marginal-likelihood",
            .Platform$file.sep,"marginal-likelihood.dat", sep="")
        mlik.res = matrix(inla.read.binary.file(file), 2, 1)
        rownames(mlik.res) = inla.namefix(c("log marginal-likelihood (integration)",
                    "log marginal-likelihood (Gaussian)"))
    }
    else
        mlik.res = NULL

    return(list(mlik=mlik.res))
}

`inla.collect.predictor` =
    function(results.dir,
             return.marginals.predictor = TRUE,
             debug = FALSE)
{
    alldir = dir(results.dir)

    ##FIRST: get the linear predictor
    subdir=paste(results.dir, .Platform$file.sep,"predictor", sep="")

    if (length(dir(subdir))>3) {
        if (debug)
            cat(paste("collect linear predictor\n", sep=""))
        
        if (debug)
            cat("...read summary.dat\n")
        file=paste(subdir, .Platform$file.sep,"summary.dat", sep="")
        dd = matrix(inla.read.binary.file(file=file), ncol=3, byrow=TRUE)[,-1, drop=FALSE]
        col.nam = c("mean","sd")

        ##get quantiles if computed
        if (length(grep("^quantiles.dat$", dir(subdir)))==1) {
            if (debug)
                cat("...read quantiles.dat\n")
            file=paste(subdir, .Platform$file.sep,"quantiles.dat", sep="")
            xx = inla.interpret.vector(inla.read.binary.file(file), debug=debug)
            len = dim(xx)[2]
            qq = xx[, seq(2, len, by=2), drop=FALSE]
            col.nam = c(col.nam, paste(as.character(xx[, 1]),"quant", sep=""))
            dd = cbind(dd, t(qq))
            rm(xx)
        }
        else {
            if (debug)
                cat("predictor: no quantiles.dat\n")
        }

        ##get cdf if computed
        if (length(grep("^cdf.dat$", dir(subdir)))==1) {
            if (debug)
                cat("...read cdf.dat\n")
            file=paste(subdir, .Platform$file.sep,"cdf.dat", sep="")
            xx = inla.interpret.vector(inla.read.binary.file(file), debug=debug)
            len = dim(xx)[2]
            qq = xx[, seq(2, len, by=2), drop=FALSE]
            col.nam = c(col.nam, paste(as.character(xx[, 1])," cdf", sep=""))
            dd = cbind(dd, t(qq))
            rm(xx)
        }
        else {
            if (debug)
                cat("... no cdf.dat\n")
        }
        
        ##get kld
        if (debug)
            cat("...read kld\n")
        kld =  matrix(inla.read.binary.file(file=paste(subdir, .Platform$file.sep,"symmetric-kld.dat", sep="")),
            ncol=2, byrow=TRUE)
        dd = cbind(dd, kld[, 2, drop=FALSE])
        col.nam = c(col.nam, "kld")
        colnames(dd) = inla.namefix(col.nam)
        summary.linear.predictor = as.data.frame(dd)
        
        if (return.marginals.predictor) {
            if (debug)
                cat("...read marginal-densities.dat\n")
            file=paste(subdir, .Platform$file.sep,"marginal-densities.dat", sep="")
            xx = inla.read.binary.file(file)
            rr = inla.interpret.vector.list(xx, debug=debug)
            rm(xx)
            if (!is.null(rr)) {
                names(rr) = inla.namefix(paste("index.", as.character(1:length(rr)), sep=""))
                names.rr = names(rr)
                for(i in 1:length(rr)) {
                    colnames(rr[[i]]) = inla.namefix(c("x", "y"))

                    if (inla.internal.experimental.mode) {
                        class(rr[[i]]) = "inla.marginal"
                        attr(rr[[i]], "inla.tag") = paste("marginal linear predictor", names.rr[i])
                    }
                }
            }

            if (inla.internal.experimental.mode) {
                class(rr) = "inla.marginals"
                attr(rr, "inla.tag") = "marginals linear predictor"
            }
            marginals.linear.predictor = rr
        } else {
            marginals.linear.predictor = NULL
        }

        ## info about size
        size.info = inla.collect.size(subdir)
    }
    else {
        summary.linear.predictor = NULL
        marginals.linear.predictor = NULL
        size.info = NULL
    }

    ##SECOND: get the inverse linear predictor(if computed)
    if (length(grep("^predictor-user-scale$", alldir))==1) {
        subdir=paste(results.dir, .Platform$file.sep,"predictor-user-scale", sep="")
        if (length(dir(subdir))>3) {
            if (debug)
                cat(paste("collect fitted values\n", sep=""))
            
            file=paste(subdir, .Platform$file.sep,"summary.dat", sep="")
            dd = matrix(inla.read.binary.file(file=file), ncol=3, byrow=TRUE)[,-1, drop=FALSE]
            col.nam = c("mean","sd")

            ##get quantiles if computed
            if (length(grep("^quantiles.dat$", dir(subdir)))==1) {
                file=paste(subdir, .Platform$file.sep,"quantiles.dat", sep="")
                xx = inla.interpret.vector(inla.read.binary.file(file), debug=debug)
                len = dim(xx)[2]
                qq = xx[, seq(2, len, by=2), drop=FALSE]
                col.nam = c(col.nam, paste(as.character(xx[, 1]),"quant", sep=""))
                dd = cbind(dd, t(qq))
                rm(xx)
            }

            ##get cdf if computed
            if (length(grep("^cdf.dat$", dir(subdir)))==1) {
                file=paste(subdir, .Platform$file.sep,"cdf.dat", sep="")
                xx = inla.interpret.vector(inla.read.binary.file(file), debug=debug)
                len = dim(xx)[2]
                qq = xx[, seq(2, len, by=2), drop=FALSE]
                col.nam = c(col.nam, paste(as.character(xx[, 1])," cdf", sep=""))
                dd = cbind(dd, t(qq))
                rm(xx)
            }
       
            colnames(dd) = inla.namefix(col.nam)
            summary.fitted.values = as.data.frame(dd)
            if (return.marginals.predictor) {
                file=paste(subdir, .Platform$file.sep,"marginal-densities.dat", sep="")
                xx = inla.read.binary.file(file)
                rr = inla.interpret.vector.list(xx, debug=debug)
                rm(xx)
                if (!is.null(rr)) {
                    names(rr) = inla.namefix(paste("index.", as.character(1:length(rr)), sep=""))
                    names.rr = names(rr)
                    for(i in 1:length(rr)) {
                        colnames(rr[[i]]) = inla.namefix(c("x", "y"))
                        if (inla.internal.experimental.mode) {
                            class(rr[[i]]) = "inla.marginal"
                            attr(rr[[i]], "inla.tag") = paste("marginal fitted values", names.rr[i])
                        }
                    }
                }

                if (inla.internal.experimental.mode) {
                    class(rr) = "inla.marginals"
                    attr(rr, "inla.tag") = "marginals fitted values"
                }
                marginals.fitted.values = rr
            } else {
                marginals.fitted.values = NULL
            }
        } else {
            summary.fitted.values = NULL
            marginals.fitted.values = NULL
        }
    } else {
        summary.fitted.values = NULL
        marginals.fitted.values = NULL
    }

    res = list(summary.linear.predictor=summary.linear.predictor,
            marginals.linear.predictor=marginals.linear.predictor,
            summary.fitted.values=summary.fitted.values,
            marginals.fitted.values=marginals.fitted.values,
            size.linear.predictor = size.info)

    return(res)
}

`inla.collect.random` =
    function(results.dir,
             return.marginals.random,
             debug = FALSE)
{
    alldir = dir(results.dir)
    random = alldir[grep("^random.effect", alldir)]
    n.random = length(random)
    if (debug)
        print("collect random effects")

    ##read the names and model of the random effects
    if (n.random>0) {
        names.random = inla.namefix(character(n.random))
        model.random = inla.trim(character(n.random))
        for(i in 1:n.random) {
            tag = paste(results.dir, .Platform$file.sep, random[i], .Platform$file.sep,"TAG", sep="")
            if (!file.exists(tag))
                names.random[i] = "missing NAME"
            else
                names.random[i] = inla.namefix(readLines(tag, n=1))
            modelname = inla.trim(paste(results.dir, .Platform$file.sep, random[i], .Platform$file.sep,"MODEL", sep=""))
            if (!file.exists(modelname))
                model.random[i] = "NoModelName"
            else
                model.random[i] = inla.trim(readLines(modelname, n=1))
        }
        
        summary.random = list()
        marginals.random = list()
        size.random = list()
        
        for(i in 1:n.random) {
            if (debug)
                print(paste("read random ", i , " of ", n.random))
            ##read the summary
            file= paste(results.dir, .Platform$file.sep, random[i], sep="")
            dir.random = dir(file)

            if (length(dir.random) > 4) {
                dd = matrix(inla.read.binary.file(file=paste(file, .Platform$file.sep,"summary.dat", sep="")), ncol=3, byrow=TRUE)
                col.nam = c("ID","mean","sd")
            
                ##read quantiles if existing
                if (debug)
                    cat("...quantiles.dat if any\n")
                if (length(grep("^quantiles.dat$", dir.random))==1) {
                    xx = inla.interpret.vector(inla.read.binary.file(paste(file, .Platform$file.sep,"quantiles.dat", sep="")),
                            debug=debug)
                    len = dim(xx)[2]
                    qq = xx[, seq(2, len, by=2), drop=FALSE]
                    col.nam = c(col.nam, paste(as.character(xx[, 1]),"quant", sep=""))
                    dd = cbind(dd, t(qq))
                }

                ##read cdf if existing
                if (debug)
                    cat("...cdf.dat if any\n")
                if (length(grep("^cdf.dat$", dir.random))==1) {
                    xx = inla.interpret.vector(inla.read.binary.file(paste(file, .Platform$file.sep,"cdf.dat", sep="")),
                            debug=debug)
                    len = dim(xx)[2]
                    qq = xx[, seq(2, len, by=2), drop=FALSE]
                    col.nam = c(col.nam, paste(as.character(xx[, 1])," cdf", sep=""))
                    dd = cbind(dd, t(qq))
                }

                ##read kld
                if (debug)
                    cat("...kld\n")
                kld1 = matrix(inla.read.binary.file(file=paste(file, .Platform$file.sep,"symmetric-kld.dat", sep="")),
                    ncol=2, byrow=TRUE)
                qq = kld1[, 2, drop=FALSE]
                dd = cbind(dd, qq)
                if (debug)
                    cat("...kld done\n")

            
                col.nam = c(col.nam, "kld")
                colnames(dd) = inla.namefix(col.nam)
                summary.random[[i]] = as.data.frame(dd)

                if (return.marginals.random) {
                    xx = inla.read.binary.file(paste(file, .Platform$file.sep,"marginal-densities.dat", sep=""))
                    rr = inla.interpret.vector.list(xx, debug=debug)
                    rm(xx)
                    if (!is.null(rr)) {
                        nd = length(rr)
                        names(rr) = inla.namefix(paste("index.", as.character(1:nd), sep=""))
                        names.rr = names(rr)
                        for(j in 1:nd) {
                            colnames(rr[[j]]) = inla.namefix(c("x", "y"))
                            if (inla.internal.experimental.mode) {
                                class(rr[[j]]) = "inla.marginal"
                                attr(rr[[j]], "inla.tag") = paste("marginal random", names.random[i], names.rr[j])
                            }
                        }
                    }

                    if (inla.internal.experimental.mode) {
                        class(rr) = "inla.marginals"
                        attr(rr, "inla.tag") = paste("marginals random",  names.random[i])
                    }
                    marginals.random[[i]] = rr
                } else {
                    marginals.random=NULL
                }
            } else {
                N.file = paste(file, .Platform$file.sep,"N", sep="")
                if (!file.exists(N.file))
                    N = 0
                else
                    N = scan(file=N.file, what = numeric(0), quiet=TRUE)
                summary.random[[i]] = data.frame("mean" = rep(NA, N), "sd" = rep(NA, N), "kld" = rep(NA, N))
                marginals.random = NULL
            }

            size.random[[i]] = inla.collect.size(file)
        }
        names(summary.random) = inla.namefix(names.random)
        if (!is.null(marginals.random) && (length(marginals.random) > 0))
            names(marginals.random) = inla.namefix(names.random)
    } else {
        if (debug)
            cat("No random effets\n")
        model.random=NULL
        summary.random=NULL
        marginals.random=NULL
        size.random = NULL
    }
    
    res = list(model.random=model.random, summary.random=summary.random, marginals.random=marginals.random, size.random = size.random)
    return(res)
}

`inla.image.reduce` = function(im, image.dim=256)
{
    ## reduce image IM to image.dim IMAGE.DIM and return the image as a matrix.
    ## order the indices so the output can be plotted by image()
    
    if ((class(im) != "pixmapGrey") || (im@size[1] != im@size[2])) {
        return (im)
    } else {
        return (im@grey)
    }

    ## do not need this anymore as we do this in GMRFLib.
    if (FALSE) {
        if (image.dim >= im@size[1]) {
            n = as.integer(im@size[1])
            x  = matrix(NA, n, n)
            for(j in 1L:n)
                x[j, n-(1L:n)+1L] = im@grey[1L:n, j]
            return (x)
        }
        block = ceiling(im@size[1]/image.dim)
        n = floor(im@size[1]/block)
        ii = jj = 0
        x = matrix(NA, n, n)
        for(i in seq(1, im@size[1]-block+1, by=block)) {
            ii = ii + 1
            jj = 0
            for(j in seq(1, im@size[1]-block+1, by=block)) {
                jj = jj + 1
                x[jj, n-ii+1] = min(im@grey[i:(i+block-1), j:(j+block-1)])
            }
        }
        return (x)
    }
}

`inla.collect.configurations` = function (results.dir, debug=FALSE) 
{
    d = paste(results.dir, "/si", sep="")
    if (file.exists(d)) {
        if (debug) {
            print(paste("collect configurations for si from dir=", d))
        }
        for (f in list.files(d, pattern="^configuration.*[.]R$", full.names=TRUE)) {
            if (debug) {
                print(paste("...read", f))
            }
            source(f)
        }
        return (list(si = inla.si.configuration))
    } else {
        if (debug) {
            print(paste("no si configurations to collect from dir=", d))
        }
        return (list(si = NULL))
    }
}

`inla.collect.offset.linear.predictor` = function(results.dir, debug = FALSE)
{
    filename = paste(results.dir, "/totaloffset/totaloffset.dat", sep="")
    stopifnot(file.exists(filename))

    xx = inla.read.binary.file(filename)
    return (list(total.offset = xx))
}
