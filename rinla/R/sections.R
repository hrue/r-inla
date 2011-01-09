### Functions to write the different sections in the .ini-file

`inla.write.boolean.field` =
    function(tag, val, file)
{
    ## write tag = 1 or tag = 0 depending on val. if val is NULL do not write
    if (!is.null(val)) {
        if (val) {
            cat(tag, " = 1\n", sep = " ", file = file, append = TRUE)
        } else {
            cat(tag, " = 0\n", sep = " ", file = file, append = TRUE)
        }
    }
}

`inla.data.section` =
    function(file, family, file.data, control, i.family="")
{
    cat("[INLA.Data", i.family, "]\n", sep = "", file = file,  append = TRUE)
    cat("type = data\n", sep = " ", file = file,  append = TRUE)
    cat("likelihood = ",family,"\n", sep = " ", file = file,  append = TRUE)
    cat("filename = ", file.data,"\n", sep = " ", file = file,  append = TRUE)

    if (inla.one.of(family, "laplace")) {
        ## two parameters, alpha and epsilon is require for LAPLACE
        cat("alpha = ", inla.ifelse(is.null(control$alpha), 0.5, control$alpha), "\n",
            sep="", file=file, append=TRUE)
        cat("epsilon = ", inla.ifelse(is.null(control$epsilon), 0.01, control$epsilon), "\n",
            sep="", file=file, append=TRUE)
        cat("gamma = ", inla.ifelse(is.null(control$gamma), 1.0, control$gamma), "\n",
            sep="", file=file, append=TRUE)
    }

    if (inla.one.of(family,c("sn", "skewnormal"))) {
        cat("sn.shape.max = ", inla.ifelse(is.null(control$sn.shape.max), 5.0, control$sn.shape.max), "\n",
            sep="", file=file, append=TRUE)
    }

    if (inla.one.of(family, "gev")) {
        cat("gev.scale.xi = ", inla.ifelse(is.null(control$gev.scale.xi), 0.01, control$gev.scale.xi), "\n",
            sep="", file=file, append=TRUE)
    }
    
    prop = inla.lmodel.properties(family)
    if (prop$ntheta > 0) {
        k = 1
        for(j in 1:prop$ntheta) {
            jj = j-1
            if (!is.null(control$fixed[j])) {
                cat("fixed", inla.ifelse(prop$ntheta==1, "", as.character(jj)), " = ",
                    inla.ifelse(control$fixed[j], 1, 0),"\n", sep = "", file = file,  append = TRUE)
            }
            if (!is.null(control$initial[j])) {
                cat("initial", inla.ifelse(prop$ntheta==1, "", as.character(jj)), " = ",
                    control$initial[j],"\n", sep = "", file = file,  append = TRUE)
            }
        }
    }
    if (prop$npriors > 0) {
        k = 1
        for(j in 1:prop$npriors) {
            jj = j-1
            if (!is.null(control$prior) && !is.null(control$prior[j])) {
                cat("prior", inla.ifelse(prop$npriors==1, "", as.character(jj)), " = ",
                    control$prior[j],"\n", sep = "", file = file,  append = TRUE)
            }
            if (!is.null(control$param)) {
                cat("parameters", inla.ifelse(prop$npriors == 1, "", as.character(jj)), " = ", sep="", file = file, append=TRUE)
                ## divide equally
                off = prop$nparameters %/% prop$npriors
                for(kk in 1:off) {
                    cat(" ", control$param[k], sep = " ", file = file,  append = TRUE)
                    k = k + 1
                }
                cat("\n", file = file,  append = TRUE)
            }
        }
    }
    if (!is.null(control$dof.max))
        cat("dof.max = ",control$dof.max,"\n", sep = " ", file = file, append = TRUE)
    cat("\n", sep = " ", file = file,  append = TRUE)
}

`inla.ffield.section` =
    function(file, file.loc, n, nrep, ngroup, 
             file.cov, file.extraconstr, file.weights, random.spec, results.dir, only.hyperparam, data.dir)
{
    ## yes, I need to do this....
    if (inla.one.of(random.spec$model, "positive")) {
        inla.offset.section(file, file.loc, n, nrep, file.cov, file.extraconstr, file.weights,
                            random.spec, results.dir, only.hyperparam, data.dir)
    } else {
        ## as it should be
        
        label= inla.namefix(random.spec$term)
        prop = inla.model.properties(random.spec$model, stop.on.error=TRUE)
    
        cat("[",label,"]\n", sep = "", file = file,  append = TRUE)
        cat("dir = ",results.dir,"\n", sep = " ", file = file,  append = TRUE)
        cat("type = ffield\n", sep = " ", file = file,  append = TRUE)
        cat("model = ",random.spec$model,"\n", sep = " ", file = file,  append = TRUE)
        if (!is.null(random.spec$same.as)) {
            cat("same.as = ",random.spec$same.as,"\n", sep = " ", file = file,  append = TRUE)
        }
        cat("covariates = ", file.cov,"\n", sep = " ", file = file,  append = TRUE)
        if (!is.null(random.spec$diagonal)) {
            cat("diagonal =", random.spec$diagonal,"\n", sep = " ", file = file,  append = TRUE)
        }
        inla.write.boolean.field("constraint", random.spec$constr, file)
        inla.write.boolean.field("si", random.spec$si, file)
        if (!is.null(file.extraconstr)) {
            cat("extraconstraint =", file.extraconstr,"\n", sep = " ", file = file,  append = TRUE)
        }
        if (!is.null(random.spec$weights)) {
            cat("weights =", file.weights,"\n", sep = " ", file = file,  append = TRUE)
        }
        if (!is.null(random.spec$spde.prefix)) {
            ## need a special one, as spde.prefix is not a file or a directory...
            fnm = inla.copy.dir.for.section.spde(random.spec$spde.prefix, data.dir)
            cat("spde.prefix =", fnm, "\n", sep = " ", file = file,  append = TRUE)
        }
        if (!is.null(random.spec$of)) {
            cat("of =", random.spec$of, "\n", sep = " ", file = file,  append = TRUE)
        }
        if (!is.null(random.spec$precision)) {
            cat("precision =", random.spec$precision, "\n", sep = " ", file = file,  append = TRUE)
        }

        if (!is.null(random.spec$range)) {
            cat("range.low  =", random.spec$range[1], "\n", sep = " ", file = file, append = TRUE)
            cat("range.high =", random.spec$range[2], "\n", sep = " ", file = file, append = TRUE)
        }

        if (prop$ntheta > 0) {
            k = 1
            for(j in 1:prop$ntheta) {
                jj = j-1
                if (!is.null(random.spec$fixed[j])) {
                    cat("fixed", inla.ifelse(prop$ntheta==1, "", as.character(jj)), " = ",
                        inla.ifelse(random.spec$fixed[j], 1, 0),"\n", sep = "", file = file,  append = TRUE)
                }
                if (!is.null(random.spec$initial[j])) {
                    cat("initial", inla.ifelse(prop$ntheta==1, "", as.character(jj)), " = ",
                        random.spec$initial[j],"\n", sep = "", file = file,  append = TRUE)
                }
            }
        }
        if (prop$npriors > 0) {
            k = 1
            for(j in 1:prop$npriors) {
                jj = j-1
                if (!is.null(random.spec$prior) && !is.null(random.spec$prior[j])) {
                    cat("prior", inla.ifelse(prop$npriors==1, "", as.character(jj)), " = ",
                        random.spec$prior[j], "\n", sep = "", file = file,  append = TRUE)
                }
                if (!is.null(random.spec$param)) {
                    cat("parameters", inla.ifelse(prop$npriors == 1, "", as.character(jj)), " = ",
                        sep="", file = file, append=TRUE)
                    ## divide equally
                    off = prop$nparameters %/% prop$npriors
                    for(kk in 1:off) {
                        cat(" ", control$param[k], sep = " ", file = file,  append = TRUE)
                        k = k + 1
                    }
                    cat("\n", file = file,  append = TRUE)
                }
            }
        }

        if (inla.model.properties(random.spec$model)$nrow.ncol) {
            cat("nrow = ", random.spec$nrow, "\n", sep = " ", file = file,  append = TRUE)
            cat("ncol = ", random.spec$ncol, "\n", sep = " ", file = file,  append = TRUE)

            if (!is.null(random.spec$bvalue)) {
                cat("bvalue = ", random.spec$bvalue, "\n", sep = " ", file = file,  append = TRUE)
            }
            if (inla.one.of(random.spec$model, c("matern2d", "matern2dx2part0", "matern2dx2p1"))) {
                if (!is.null(random.spec$nu)) {
                    cat("nu = ", random.spec$nu, "\n", sep = " ", file = file,  append = TRUE)
                }
            }
        } else {
            cat("n = ",n,"\n", sep = " ", file = file,  append = TRUE)
        }
        cat("nrep = ", inla.ifelse(is.null(nrep), 1, nrep), "\n", sep = " ", file = file,  append = TRUE)

        if (!is.null(ngroup) && ngroup > 1) {
            cat("ngroup = ", ngroup, "\n", sep = " ", file = file,  append = TRUE)
            if (!is.null(random.spec$control.group$prior)) {
                cat("group.prior = ", random.spec$control.group$prior, "\n", sep = " ", file = file,  append = TRUE)
            }
            if (!is.null(random.spec$control.group$param)) {
                cat("group.parameters = ", random.spec$control.group$param, "\n", sep = " ", file = file,  append = TRUE)
            }
            if (!is.null(random.spec$control.group$initial)) {
                cat("group.initial = ", random.spec$control.group$initial, "\n", sep = " ", file = file,  append = TRUE)
            }
            if (!is.null(random.spec$control.group$fixed)) {
                if (random.spec$control.group$fixed) {
                    cat("group.fixed = ", 1, "\n", sep = " ", file = file,  append = TRUE)
                } else {
                    cat("group.fixed = ", 0, "\n", sep = " ", file = file,  append = TRUE)
                }
            }
            if (!is.null(random.spec$control.group$model)) {
                cat("group.model = ", random.spec$control.group$model, "\n", sep = " ", file = file,  append = TRUE)
            }
        }
        
        if (!is.null(random.spec$cyclic)) {
            cat("cyclic = ",as.numeric(random.spec$cyclic),"\n", sep = " ", file = file,  append = TRUE)
        }
        if (!is.null(random.spec$season.length)) {
            cat("season = ",random.spec$season.length,"\n", sep = " ", file = file,  append = TRUE)
        }
        if (!is.null(random.spec$graph.file)) {
            fnm = inla.copy.file.for.section(random.spec$graph.file, data.dir)
            cat("graph = ",fnm, "\n", sep = " ", file = file,  append = TRUE)
        }
        if (!is.null(file.loc)) {
            cat("locations = ", file.loc,"\n", sep = " ", file = file,  append = TRUE)
        }

        if (!is.null(random.spec$Cmatrix)) {
            if (is.character(random.spec$Cmatrix)) {
                fnm = inla.copy.file.for.section(random.spec$Cmatrix, data.dir)
                cat("Cmatrix = ", fnm, "\n", append=TRUE, sep = " ", file = file)
            } else {
                file.C = inla.tempfile(tmpdir=data.dir)
                inla.sparse2file(random.spec$Cmatrix, file.C, c.indexing = TRUE, symmetric = TRUE)
                file.C = gsub(data.dir, "$inladatadir", file.C, fixed=TRUE)
                cat("Cmatrix = ",file.C, "\n",append=TRUE, sep = " ", file = file)
            }
        }
        if (!is.null(random.spec$rankdef)) {
            cat("rankdef = ",random.spec$rankdef,"\n",append=TRUE, sep = " ", file = file)
        }
        if (!is.null(random.spec$cdf)) {
            cat("cdf = ",random.spec$cdf, "\n", sep = " ", file = file,  append = TRUE)
        }
        if (!is.null(random.spec$quantiles)) {
            cat("quantiles = ",random.spec$quantiles, "\n", sep = " ", file = file,  append = TRUE)
        }
        if (only.hyperparam || !random.spec$compute) {
            cat("compute = 0\n", sep = " ", file = file,  append = TRUE)
        } else {
            cat("compute = 1\n", sep = " ", file = file,  append = TRUE)
        }
        cat("\n", sep = " ", file = file,  append = TRUE)
    }
}

`inla.offset.section` =
    function(file, file.loc, n, nrep, file.cov, file.extraconstr, file.weights,
             random.spec, results.dir, only.hyperparam, data.dir)
{
    ## special version of inla.ffield.section() but for type = offset
    ##
    label= inla.namefix(random.spec$term)
    prop = inla.model.properties(random.spec$model, stop.on.error=TRUE)
    
    cat("[",label,"]\n", sep = "", file = file,  append = TRUE)
    cat("dir = ",results.dir,"\n", sep = " ", file = file,  append = TRUE)
    cat("type = offset\n", sep = " ", file = file,  append = TRUE)
    cat("model = ",random.spec$model,"\n", sep = " ", file = file,  append = TRUE)
    if (!is.null(random.spec$weights))
        cat("weights =", file.weights,"\n", sep = " ", file = file,  append = TRUE)

    if (prop$ntheta > 0) {
        k = 1
        for(j in 1:prop$ntheta) {
            jj = j-1
            if (!is.null(random.spec$fixed[j])) {
                cat("fixed", inla.ifelse(prop$ntheta==1, "", as.character(jj)), " = ",
                    inla.ifelse(random.spec$fixed[j], 1, 0),"\n", sep = "", file = file,  append = TRUE)
            }
            if (!is.null(random.spec$initial[j])) {
                cat("initial", inla.ifelse(prop$ntheta==1, "", as.character(jj)), " = ",
                    random.spec$initial[j],"\n", sep = "", file = file,  append = TRUE)
            }
        }
    }
    if (prop$npriors > 0) {
        k = 1
        for(j in 1:prop$npriors) {
            jj = j-1
            if (!is.null(random.spec$prior) && !is.null(random.spec$prior[j])) {
                cat("prior", inla.ifelse(prop$npriors==1, "", as.character(jj)), " = ",
                    random.spec$prior[j], "\n", sep = "", file = file,  append = TRUE)
            }
            if (!is.null(random.spec$param)) {
                cat("parameters", inla.ifelse(prop$npriors == 1, "", as.character(jj)), " = ",
                    sep="", file = file, append=TRUE)
                ## divide equally
                off = prop$nparameters %/% prop$npriors
                for(kk in 1:off) {
                    cat(" ", random.spec$param[k], sep = " ", file = file,  append = TRUE)
                    k = k + 1
                }
                cat("\n", file = file,  append = TRUE)
            }
        }
    }

    ##cat("n = ",n,"\n", sep = " ", file = file,  append = TRUE)

    if (!is.null(random.spec$cdf)) {
        cat("cdf = ",random.spec$cdf, "\n", sep = " ", file = file,  append = TRUE)
    }
    if (!is.null(random.spec$quantiles)) {
        cat("quantiles = ",random.spec$quantiles, "\n", sep = " ", file = file,  append = TRUE)
    }
    
    cat("\n", sep = " ", file = file,  append = TRUE)
}

`inla.inla.section` =
    function(file,inla.spec)
{
    cat("[INLA.Parameters]\n", sep = " ", file = file,  append = TRUE)
    cat("type = inla\n", sep = " ", file = file,  append = TRUE)

    if (!is.null(inla.spec$int.strategy)) {
        cat("int.strategy = ",inla.spec$int.strategy,"\n", sep = " ", file = file,  append = TRUE)
    }
    if (!is.null(inla.spec$strategy)) {
        cat("strategy = ",inla.spec$strategy,"\n", sep = " ", file = file,  append = TRUE)
    }
    if (!is.null(inla.spec$linear.correction)) {
        cat("linear.correction = ",inla.spec$linear.correction,"\n", sep = " ", file = file,  append = TRUE)
    }
    if (!is.null(inla.spec$h)) {
        cat("h = ",inla.spec$h,"\n", sep = " ", file = file,  append = TRUE)
    }
    if (!is.null(inla.spec$dz)) {
        cat("dz = ",inla.spec$dz,"\n", sep = " ", file = file,  append = TRUE)
    }
    if (!is.null(inla.spec$interpolator)) {
        cat("interpolator = ", inla.spec$interpolator,"\n", sep = " ", file = file,  append = TRUE)
    }
    if (!is.null(inla.spec$diff.logdens)) {
        cat("diff.log.dens = ",inla.spec$diff.logdens,"\n", sep = " ", file = file,  append = TRUE)
    }
    if (!is.null(inla.spec$print.joint.hyper)) {
        cat("fp.hyperparam = $inlaresdir/joint.dat\n", sep = "", file = file,  append = TRUE)
    }
    if (!is.null(inla.spec$tolerance)) {
        ## tolerance is a generic option, which sets 'all' tolerance to the given value
        cat("epsg = ", inla.spec$tol,"\n", sep = " ", file = file,  append = TRUE)
        cat("epsf = ", inla.spec$tol,"\n", sep = " ", file = file,  append = TRUE)
        cat("epsx = ", inla.spec$tol,"\n", sep = " ", file = file,  append = TRUE)
        cat("abserr.func = ", inla.spec$tol,"\n", sep = " ", file = file,  append = TRUE)
        cat("abserr.step = ", inla.spec$tol,"\n", sep = " ", file = file,  append = TRUE)
    }
    
    inla.write.boolean.field("hessian.force.diagonal", inla.spec$force.diagonal, file)
    inla.write.boolean.field("skip.configurations", inla.spec$skip.configurations, file)
    inla.write.boolean.field("mode.known", inla.spec$mode.known.conf, file)
    inla.write.boolean.field("adjust.weights", inla.spec$adjust.weights, file)
    inla.write.boolean.field("derived.only", inla.spec$derived.only, file)

    if (!is.null(inla.spec$restart) && inla.spec$restart >= 0) {
        cat("restart = ", as.integer(inla.spec$restart), "\n", file = file, sep = " ", append = TRUE)
    }

    if (!is.null(inla.spec$optimiser)) {
        cat("optimiser = ", inla.spec$optimiser,"\n", sep = " ", file = file,  append = TRUE)
    }
    if (!is.null(inla.spec$verbose) && inla.spec$verbose) {
        cat("optpar.fp = stdout\n", sep = " ", file = file,  append = TRUE)
    }

    if (!is.null(inla.spec$reordering)) {
        ## reordering could be a number -1, 0, ....  or a string.
        r = inla.spec$reordering
        if (is.character(r)) {
            r.idx = which(tolower(r) == tolower(c("AUTO", "DEFAULT", "IDENTITY", "BAND", "METIS", "GENMMD", "AMD", "MD")))
            if (length(r.idx) > 0) {
                r.idx = r.idx - 2
            } else {
                r.idx = -1
            }
        } else {
            r.idx = r
        }
        cat("reordering = ", r.idx, "\n", sep = " ", file = file,  append = TRUE)
    }

    if (!is.null(inla.spec$cpo.diff)) {
        cat("cpo.diff = ", inla.spec$cpo.diff, "\n", sep = " ", file = file,  append = TRUE)
    }
    if (!is.null(inla.spec$npoints)) {
        cat("npoints = ", inla.spec$npoints, "\n", sep = " ", file = file,  append = TRUE)
    }
    inla.write.boolean.field("adapt.hessian.mode", inla.spec$adapt.hessian.mode, file)

    if (!is.null(inla.spec$adapt.hessian.max.trials) && inla.spec$adapt.hessian.max.trials >= 0) {
        cat("adapt.hessian.max.trials = ", inla.spec$adapt.hessian.max.trials, "\n", sep = " ", file = file,  append = TRUE)
    }
    if (!is.null(inla.spec$adapt.hessian.scale) && inla.spec$adapt.hessian.scale >= 1) {
        cat("adapt.hessian.scale = ", inla.spec$adapt.hessian.scale, "\n", sep = " ", file = file,  append = TRUE)
    }
    if (!is.null(inla.spec$step.len)) {
        cat("step.len = ", inla.spec$step.len, "\n", sep = " ", file = file,  append = TRUE)
    }
    if (!is.null(inla.spec$diagonal) && inla.spec$diagonal >= 0.0) {
        cat("diagonal = ", inla.spec$diagonal, "\n", sep = " ", file = file,  append = TRUE)
    }
    if (!is.null(inla.spec$numint.maxfeval)) {
        cat("numint.maxfeval = ", as.integer(inla.spec$numint.maxfeval), "\n", file = file, append = TRUE)
    }
    if (!is.null(inla.spec$numint.relerr)) {
        cat("numint.relerr = ", inla.spec$numint.relerr, "\n", file = file, append = TRUE)
    }
    if (!is.null(inla.spec$numint.abserr)) {
        cat("numint.abserr = ", inla.spec$numint.abserr, "\n", file = file, append = TRUE)
    }

    cat("\n", sep = " ", file = file,  append = TRUE)
}

`inla.predictor.section` =
    function(file, n, predictor.spec, file.offset, data.dir)
{
    cat("[Predictor]\n", sep = " ", file = file,  append = TRUE)
    cat("type = predictor\n", sep = " ", file = file,  append = TRUE)
    cat("dir = predictor\n", sep = " ", file = file, append = TRUE)
    cat("n = ",n, "\n", sep = " ", file = file,  append = TRUE)

    inla.write.boolean.field("fixed", predictor.spec$fixed, file)
    inla.write.boolean.field("compute", predictor.spec$compute, file)
    
    if (!is.null(predictor.spec$cdf)) {
        cat("cdf = ",predictor.spec$cdf, "\n", sep = " ", file = file,  append = TRUE)
    }
    if (!is.null(predictor.spec$quantiles)) {
        cat("cdf = ",predictor.spec$quantiles, "\n", sep = " ", file = file,  append = TRUE)
    }
    if (!is.null(file.offset)) {
        cat("offset = ", file.offset,"\n", sep = " ", file = file, append=TRUE)
    }
    if (!is.null(predictor.spec$prior)) {
        cat("prior = ", predictor.spec$prior, sep = " ", file = file, append=TRUE)
    }
    ## hard-coded to have two parameters... FIX LATER
    if (!is.null(predictor.spec$param)) {
        cat("parameters = ",predictor.spec$param[1]," ",predictor.spec$param[2],"\n", sep = " ", file = file, append=TRUE)
    }
    if (!is.null(predictor.spec$initial)) {
        cat("initial = ",predictor.spec$initial,"\n", sep = " ", file = file, append=TRUE)
    } else {
        if (predictor.spec$fixed) {
            cat("initial = 10\n", sep = " ", file = file,  append = TRUE)
        }
    }
    if (!is.null(predictor.spec$predictor.usermap)){
        cat("predictor.usermap=", predictor.spec$predictor.usermap, "\n", sep=" ", file = file, append = TRUE)
    }
    if (!is.null(predictor.spec$cross) && length(predictor.spec$cross) > 0) {
        if (length(predictor.spec$cross) != n) {
            stop(paste("Length of cross does not match length of predictor", length(predictor.spec$cross), "!=", n))
        }
        file.cross = inla.tempfile(tmpdir=data.dir)
        predictor.spec$cross[is.na(predictor.spec$cross)] = 0
        write(predictor.spec$cross, ncol=1, file=file.cross)
        fnm = gsub(data.dir, "$inladatadir", file.cross, fixed=TRUE)
        cat("cross.constraint =", fnm, "\n", file=file, append = TRUE)
    }

    if (!is.null(predictor.spec$A)) {
        ## Since we will expand the A into [ I, -A; -A^T, A^T A
        ## ], we have to read the matrix it if its given in a
        ## file. The fileformat is given as a three colums file, with
        ## columns i, j and values.
        if (is.character(predictor.spec$A)) {
            A = read.table(predictor.spec$A, col.names = c("i", "j", "x"))
            A = sparseMatrix(i = A$i, j = A$j, x = A$x, index1=TRUE)
        } else {
            A = predictor.spec$A
        }
        if (!is(A, "dgTMatrix")) {
            A = inla.sparse.check(A)
        }

        ## Now we will build the extended Matrix, which is
        ##
        ## Aextended = [ I, -A; -A^T, A^T A ] (2n x 2x)
        ##
        ## This matrix is the one that is needed for input to inla. We
        ## will only build the upper triangular part!
        
        ## Note that only a subset of the matrix has to be
        ## given, and the rest is filled with zeros:
        ## > inla.sparse2dgTMatrix(list(i=2, j=2, values=2), dims=c(3,3))
        ##   3 x 3 sparse Matrix of class "dgTMatrix"
        ##   [1,] . . .
        ##   [2,] . 2 .
        ##   [3,] . . .
        ##

        ## if A is incomplete, just fill it with unit diagonals. NOTE:
        ## If this takes time, then we may have to require that A is
        ## completely specified!
        nA = dim(A)[1]
        if (nA < n) {
            A = inla.as.dgTMatrix(
                    sparseMatrix(i = c(A@i+1L, (nA+1):n), j = c(A@j+1L, (nA+1):n),
                    x = c(A@x, rep(1.0, n-nA)), index1 = TRUE, dims = c(n,n))
                    )
            stopifnot(dim(A)[1] == n)
        }

        ## The `I'
        Aext = list(i = 1L:n, j = 1L:n, x = rep(1.0,n))

        ## add -A. Ooops; the internal storage @i etc, are zero-based indexing.
        Aext$i = c(Aext$i, (A@i+1L))
        Aext$j = c(Aext$j, (A@j+1L) + n)
        Aext$x = c(Aext$x, -A@x)

        ## add -A^T. Ooops; the internal storage @i etc, are zero-based indexing.
        Aext$i = c(Aext$i, (A@j+1L) + n)
        Aext$j = c(Aext$j, (A@i+1L))
        Aext$x = c(Aext$x, -A@x)

        ## add A^T A. Ooops; the internal storage @i etc, are zero-based indexing.
        ATA = inla.as.dgTMatrix( t(A) %*% A )
        Aext$i = c(Aext$i, (ATA@i+1L) + n)
        Aext$j = c(Aext$j, (ATA@j+1L) + n)
        Aext$x = c(Aext$x, ATA@x)

        stopifnot(length(Aext$i) == length(Aext$j))
        stopifnot(length(Aext$i) == length(Aext$x))
        
        file.A=inla.tempfile(tmpdir=data.dir)
        Aext = sparseMatrix(i = Aext$i, j = Aext$j, x = Aext$x)
        inla.sparse2file(Aext, filename = file.A, c.indexing = TRUE, symmetric = TRUE)
        
        file.A = gsub(data.dir, "$inladatadir", file.A, fixed=TRUE)
        cat("Aext = ", file.A, "\n", append=TRUE, sep = " ", file = file)
        cat("precision = ", predictor.spec$precision, "\n", append=TRUE, sep = " ", file = file)
    }

    cat("\n", sep = " ", file = file,  append = TRUE)
}

`inla.problem.section` =
    function(file , data.dir, result.dir, hyperpar, return.marginals, dic, cpo, mlik, quantiles, smtp, q, strategy)
{
    cat("", sep = "", file = file, append=FALSE)
    cat("inladatadir = ",data.dir, "\n", sep = "", file = file,  append = TRUE)
    cat("inlaresdir = ",result.dir, "\n", sep = "", file = file,  append = TRUE)
    cat("#inladatadir = ",gsub("^.*/","",data.dir), "\n", sep = "", file = file,  append = TRUE) #
    cat("#inlaresdir = ", gsub("^.*/","",result.dir), "-%d\n", sep = "", file = file,  append = TRUE) #

    cat("\n", sep = " ", file = file,  append = TRUE)
    cat("[INLA.Model]\n", sep = " ", file = file,  append = TRUE)
    cat("type = problem\n", sep = " ", file = file,  append = TRUE)
    cat("dir = $inlaresdir\n", sep = " ", file = file,  append = TRUE)
    cat("strategy = ", strategy, "\n", sep = " ", file = file,  append = TRUE)
    inla.write.boolean.field("return.marginals", return.marginals, file)
    inla.write.boolean.field("hyperparameters", hyperpar, file)
    inla.write.boolean.field("cpo", cpo, file)
    inla.write.boolean.field("dic", dic, file)
    inla.write.boolean.field("mlik", mlik, file)
    inla.write.boolean.field("q", q, file)

    if (!is.null(smtp)) {
        cat("smtp = ",smtp, "\n", sep = " ", file = file,  append = TRUE)
    }
    if (!is.null(quantiles)) {
        cat("quantiles = ",quantiles, "\n", sep = " ", file = file,  append = TRUE)
    }
    cat("\n", sep = " ", file = file,  append = TRUE)
}

`inla.parse.fixed.prior` = function(name, prior)
{
    if (is.null(prior)) {
        return (NULL)
    } else if (is.numeric(prior)) {
        return (prior)
    } else {
        if (is.null(inla.eval(paste("prior$", name, sep="")))) {
            if (!is.null(prior$default)) {
                return (prior$default)
            } else {
                return (NULL)
            }
        } else {
            return (inla.eval(paste("prior$", name, sep="")))
        }
    }
}

`inla.linear.section` =
    function(file, file.fixed, label, results.dir, control.fixed, only.hyperparam)
{
    cat("[", inla.namefix(label), "]\n", sep = "", file = file,  append = TRUE)
    cat("dir = ",results.dir,"\n", sep = " ", file = file,  append = TRUE)
    cat("type = linear\n", sep = " ", file = file,  append = TRUE)
    cat("covariates = ", file.fixed,"\n", sep = " ", file = file,  append = TRUE)
    if (only.hyperparam || !control.fixed$compute) {
        cat("compute = 0\n", sep = " ", file = file,  append = TRUE)
    } else {
        cat("compute = 1\n", sep = " ", file = file,  append = TRUE)
    }
    if (!is.null(control.fixed$cdf)) {
        cat("cdf = ",control.fixed$cdf, "\n", sep = " ", file = file,  append = TRUE)
    }
    if (!is.null(control.fixed$quantiles)) {
        cat("cdf = ",control.fixed$quantiles, "\n", sep = " ", file = file,  append = TRUE)
    }
    if (length(grep("^[(]Intercept[)]$", inla.trim(label))) == 1) {
        prec = control.fixed$prec.intercept
        mean = control.fixed$mean.intercept
        if (!is.null(mean)) {
            cat("mean = ", mean, "\n", sep = " ", file = file,  append = TRUE)
        }
        if (!is.null(prec)) {
            cat("precision = ", prec, "\n", sep = " ", file = file,  append = TRUE)
        }
    } else {
        prec = inla.parse.fixed.prior(label, control.fixed$prec)
        mean = inla.parse.fixed.prior(label, control.fixed$mean)
        if (!is.null(mean)) {
            cat("mean = ", mean, "\n", sep = " ", file = file,  append = TRUE)
        }
        if (!is.null(prec)) {
            cat("precision = ", prec, "\n", sep = " ", file = file,  append = TRUE)
        }
    }
    cat("\n", sep = " ", file = file,  append = TRUE)
}

`inla.mode.section` =
    function(file, args, data.dir)
{
    if (!is.null(args$result) || !is.null(args$theta) || !is.null(args$x)) {

        cat("[INLA.Control.Mode]\n", sep = " ", file = file,  append = TRUE)
        cat("type = mode\n", sep = " ", file = file,  append = TRUE)

        ## use default the mode in result if given
        if (is.null(args$theta) && !is.null(args$result)) {
            args$theta = args$result$mode$theta
        }
        
        if (!is.null(args$theta)) {
            file.theta = inla.tempfile(tmpdir=data.dir)
            ##cat("theta = ", inla.paste(as.character(args$theta)), "\n", sep = " ", file = file,  append = TRUE)
            fp.binary = file(file.theta, "wb")
            ## convert from character to a vector of doubles
            if (is.character(args$theta)) {
                args$theta = as.numeric(strsplit(args$theta, "[ \t]+")[[1]])
            }
            args$theta = args$theta[!is.na(args$theta)]
            writeBin(as.integer(length(args$theta)), fp.binary)
            writeBin(args$theta, fp.binary)
            close(fp.binary)
            fnm = gsub(data.dir, "$inladatadir", file.theta, fixed=TRUE)
            cat("theta =", fnm, "\n", file=file, append = TRUE)
        }
        ## use default the mode in result if given
        if (is.null(args$x) && !is.null(args$result)) {
            args$x = args$result$mode$x
        }
        if (!is.null(args$x)) {
            file.x = inla.tempfile(tmpdir=data.dir)
            ##write(args$x, ncol=1, file=file.x)
            fp.binary = file(file.x, "wb")
            writeBin(as.integer(length(args$x)), fp.binary)
            writeBin(args$x, fp.binary)
            close(fp.binary)
            fnm = gsub(data.dir, "$inladatadir", file.x, fixed=TRUE)
            cat("x =", fnm, "\n", file=file, append = TRUE)
        }
        inla.write.boolean.field("restart", args$restart, file)
    }
}

`inla.expert.section` =
    function(file, args)
{
    if (!is.null(args$cpo.manual) && args$cpo.manual) {
        cat("\n## If you edit this section it is assumed you know what you're doing ;-)\n", file=file, append=TRUE) #
        cat("[INLA.Expert]\n", sep = " ", file = file,  append = TRUE)
        cat("type = expert\n", sep = " ", file = file,  append = TRUE)
        inla.write.boolean.field("cpo.manual", args$cpo.manual, file)
        ## recall to convert to 0-based index'ing
        cat("cpo.idx = ", args$cpo.idx -1,"\n", sep = " ", file = file,  append = TRUE)
    }
}

`inla.lincomb.section` =
    function(file, data.dir, contr, lincomb)
{
    ## this one write binary format files...

    ## format is either
    ##
    ##     list("lc1" = "a 1 1 b 2 1 3 2...", ...)
    ##
    ## or
    ##
    ##     list("lc1" = list( "a" = list(idx=1, weight=1), "b" = list(idx=c(2,3), weight = c(1,2)), ...), ...)
    ##       
    ## use the functions 'inla.make.lincomb()' and 'inla.make.lincombs()'
    
    ## if use.one.file = TRUE, then use one file for all lincombs and
    ## the 'ENTRY keyword', otherwise, use one file for each lincomb.

    if (!is.null(lincomb)) {

        fnm = inla.tempfile(tmpdir=data.dir)
        file.create(fnm)
        fp.binary = file(fnm, "wb")
        stopifnot(!is.null(fnm))
        stopifnot(!is.null(fp.binary))

        numlen = inla.numlen(length(lincomb))

        for(i in 1:length(lincomb)) {
            
            if (is.null(names(lincomb[i])) || is.na(names(lincomb[i]))) {
                secname = paste("lincomb.", inla.num(i, width=numlen), sep="")
                lc = lincomb[[i]]
            } else if (names(lincomb[i]) == "") {
                secname = paste("lincomb.", inla.num(i, width=numlen), sep="")
                lc = lincomb[[i]]
            } else {
                secname = paste("lincomb.", names(lincomb[i])[1], sep="")
                lc = lincomb[[i]]
            }
                
            cat("\n[", secname, "]\n", sep = "", file = file,  append = TRUE)
            cat("type = lincomb\n", sep = " ", file = file,  append = TRUE)
            if (!is.null(contr$precision)) {
                cat("precision = ", contr$precision,"\n", sep = " ", file = file,  append = TRUE)
            }
            if (!is.null(contr$usermap)) {
                cat("usermap = ", contr$usermap,"\n", sep = " ", file = file,  append = TRUE)
            }
            inla.write.boolean.field("verbose", contr$verbose, file)

            cat("file.offset = ", as.integer(seek(fp.binary, where=NA)), "\n", sep="", file = file, append = TRUE)
            
            ## number of entries
            writeBin(as.integer(length(lc)), fp.binary)

            for(j in 1:length(lc)) {

                lc.j.name = as.character( names(lc[[j]]) )
                lc.j = lc[[j]][[1]]

                ## if $idx is not there, its default 1.
                if (is.null(lc.j$idx) && length(lc.j$weight) == 1) {
                    lc.j$idx = 1
                }

                ## NA's are allowed; but we just remove them.
                idx = lc.j$idx[ !is.na(lc.j$idx) ]
                weight = lc.j$weight[ !is.na(lc.j$weight) ]
                if (length(idx) == 0 || length(weight) == 0)
                    stop(paste("lincomb", secname, "has only zero entries. This is not allowed"))
                stopifnot(length(idx) == length(weight))

                ## this the old code:
                ## cat(c( lc.j.name, c(rbind( idx, weight))), "\n", file=fnm, append=TRUE)
                writeBin(as.integer(nchar(lc.j.name)), fp.binary)
                writeBin(as.character(lc.j.name), fp.binary)
                writeBin(as.integer(length(idx)), fp.binary) ## number of pairs
                writeBin(as.integer(idx), fp.binary)
                writeBin(as.double(weight), fp.binary)

                ##print(paste(" stop writing at position ", as.integer(seek(fp.binary, where=NA))))
            }

            fnm.new = gsub(data.dir, "$inladatadir", fnm, fixed=TRUE)
            cat("filename = ", fnm.new, "\n", sep = " ", file = file, append = TRUE)
        }

        close(fp.binary)
    }
}

`inla.lincomb.section.OLDFORMAT` =
    function(file, data.dir, contr, lincomb, use.one.file = TRUE)
{
    ## format is either
    ##
    ##     list("lc1" = "a 1 1 b 2 1 3 2...", ...)
    ##
    ## or
    ##
    ##     list("lc1" = list( "a" = list(idx=1, weight=1), "b" = list(idx=c(2,3), weight = c(1,2)), ...), ...)
    ##       

    ## if use.one.file = TRUE, then use one file for all lincombs and
    ## the 'ENTRY keyword', otherwise, use one file for each lincomb.

    if (!is.null(lincomb)) {

        ## I need this to be set upfront
        fnm = NULL

        numlen = inla.numlen(length(lincomb))

        for(i in 1:length(lincomb)) {
            
            if (is.null(names(lincomb[i])) || is.na(names(lincomb[i]))) {
                secname = paste("lincomb.", inla.num(i, width=numlen), sep="")
                lc = lincomb[[i]]
            } else if (names(lincomb[i]) == "") {
                secname = paste("lincomb.", inla.num(i, width=numlen), sep="")
                lc = lincomb[[i]]
            } else {
                secname = paste("lincomb.", names(lincomb[i])[1], sep="")
                lc = lincomb[[i]]
            }
                
            cat("\n[", secname, "]\n", sep = "", file = file,  append = TRUE)
            cat("type = lincomb\n", sep = " ", file = file,  append = TRUE)
            if (!is.null(contr$precision)) {
                cat("precision = ", contr$precision,"\n", sep = " ", file = file,  append = TRUE)
            }
            if (!is.null(contr$usermap)) {
                cat("usermap = ", contr$usermap,"\n", sep = " ", file = file,  append = TRUE)
            }

            if (use.one.file) {
                ## create file if we havn't done already
                if (is.null(fnm)) {
                    fnm = inla.tempfile(tmpdir=data.dir)
                    file.create(fnm)
                }
                ## use the ENTRY option, set option entry and add ENTRY to the file
                cat("entry = ", secname, "\n", sep = "", file = file, append = TRUE)
                cat("entryoffset = ", file.info(fnm)$size, "\n", sep="", file = file, append = TRUE)
                cat("\nENTRY ", secname, "\n", sep= "", file = fnm, append = TRUE)
            } else {
                ## need new file for each lincomb
                fnm = inla.tempfile(tmpdir=data.dir)
                file.create(fnm)
            }
            
            if (is.character(lc)) {
                write(lc, file=fnm)
            } else {
                for(j in 1:length(lc)) {
                    lc.j.name = names(lc[[j]])
                    lc.j = lc[[j]][[1]]

                    ## if $idx is not there, its default 1.
                    if (is.null(lc.j$idx) && length(lc.j$weight) == 1) {
                        cat(c( lc.j.name, 1, lc.j$weight), "\n", file=fnm, append=TRUE)
                    } else {
                        idx = lc.j$idx[ !is.na(lc.j$idx) ]
                        weight = lc.j$weight[ !is.na(lc.j$weight) ]
                        stopifnot(length(idx) == length(weight))
                        cat(c( lc.j.name, c(rbind( idx, weight))), "\n", file=fnm, append=TRUE)
                    }
                }
            }

            fnm.new = gsub(data.dir, "$inladatadir", fnm, fixed=TRUE)
            cat("filename = ", fnm.new, "\n", sep = " ", file = file, append = TRUE)
        }
    }
}

`inla.copy.file.for.section` =
    function(filename, data.dir)
{
    if (missing(filename)) {
        return (NULL)
    }
    if (missing(data.dir)) {
        stop("data.dir required")
    }

    fnm = inla.tempfile(tmpdir=data.dir)
    file.copy(filename, fnm, overwrite=TRUE)
    return (gsub(data.dir, "$inladatadir", fnm, fixed=TRUE))
}

`inla.copy.dir.for.section` =
    function(dir.name, data.dir)
{
    d.fnm = inla.tempfile(tmpdir=data.dir)
    inla.dir.create(d.fnm)
    files.to.copy = paste(dir.name, "/", dir(dir.name, recursive=TRUE), sep="")
    file.copy(files.to.copy, d.fnm, recursive=TRUE)
    return (gsub(data.dir, "$inladatadir", d.fnm, fixed=TRUE))
}

`inla.copy.dir.for.section.spde` =
    function(prefix, data.dir)
{
    dir.name = dirname(prefix)
    file.prefix = basename(prefix)

    d.fnm = inla.tempfile(tmpdir=data.dir)
    inla.dir.create(d.fnm)
    files.to.copy = paste(dir.name, "/",
            dir(dir.name, pattern = paste("^", file.prefix, sep=""), recursive=TRUE), sep="")
    file.copy(files.to.copy, d.fnm, recursive=TRUE)
    rdir = gsub(data.dir, "$inladatadir", d.fnm, fixed=TRUE)
    rprefix = paste(rdir, "/", file.prefix, sep="")
    return (rprefix)
}

`inla.z.section` =
    function(file, random.spec, data.dir, results.dir, only.hyperparam, k.off)
{
    label= inla.namefix(random.spec$term)
    if (!is.matrix(random.spec$Z)) {
        stop("Argument Z in model=[z]has to be a matrix.")
    }
    n = dim(random.spec$Z)[1]
    m = dim(random.spec$Z)[2]
    ind = 0:(n-1)
    
    for(k in 1:m) {
        kk = k + k.off
        cat("[",label,".", k, "]\n", sep = "", file = file,  append = TRUE)
        cat("type = ffield\n", sep = " ", file = file,  append = TRUE)
        cat("dir = ", "random.effect", inla.num(kk), "\n", sep = "", file = file, append = TRUE)
        cat("model = ", inla.ifelse(k == 1, "z", "zadd"),"\n", sep = " ", file = file,  append = TRUE)
        cat("n = 1\n", file=file, append = TRUE)
        if (k == 1) {
            if (!is.null(random.spec$param)) {
                cat("parameters = ",random.spec$param[1]," ",random.spec$param[2],"\n", sep = " ", file = file,
                    append = TRUE)
            }
            cat("prior = ",random.spec$prior.logprec,"\n", sep = " ", file = file,  append = TRUE)
        }

        file.cov=inla.tempfile(tmpdir=data.dir)
        file.create(file.cov)
        write(t(cbind(ind,rep(0,n))),ncolumns=2,file=file.cov,append=FALSE)
        file.cov = gsub(data.dir, "$inladatadir", file.cov, fixed=TRUE)
        cat("covariates = ", file.cov,"\n", sep = " ", file = file,  append = TRUE)

        file.w=inla.tempfile(tmpdir=data.dir)
        file.create(file.w)
        write(t(cbind(ind,random.spec$Z[,k])),ncolumns=2,file=file.w,append=FALSE)
        file.w = gsub(data.dir, "$inladatadir", file.w, fixed=TRUE)
        cat("weights = ", file.w,"\n", sep = " ", file = file,  append = TRUE)

        if (only.hyperparam) {
            cat("compute = 0\n", sep = " ", file = file,  append = TRUE)
        }
        cat("\n", sep = " ", file = file,  append = TRUE)
    }
}

        

    
