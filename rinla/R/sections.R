### Functions to write the different sections in the .ini-file

`inla.write.hyper` = function(hyper, file, prefix="")
{
    stopifnot(!missing(hyper))
    stopifnot(!missing(file))
    
    len = length(hyper)
    if (len == 0L) {
        return ()
    }

    for(k in 1L:len) {
        if (len == 1L) {
            suff = ""
        } else {
            suff = as.character(k-1L)
        }
        cat(prefix, "initial",    suff, " = ", hyper[[k]]$initial, "\n", file = file, append = TRUE, sep="")
        cat(prefix, "fixed",      suff, " = ", as.numeric(hyper[[k]]$fixed), "\n", file = file, append = TRUE, sep="")

        ## these are for "expression:"...
        ## if there are newlines,  remove them
        tmp.prior = gsub("\n", "", hyper[[k]]$prior)
        ## remove preceding spaces
        tmp.prior = gsub("^[ \t]+", "", tmp.prior)
        ## if the expression ends with a ";" with or without spaces, remove it
        tmp.prior = gsub(";*[ \t]*$", "", tmp.prior)
        ## for all priors except the "expression:" one,  then trim the name
        if (length(grep("^expression[ \t]*:", tolower(tmp.prior))) == 0L) {
            tmp.prior = inla.trim.family(tmp.prior)
        }
        cat(prefix, "prior",      suff, " = ", tmp.prior, "\n", file = file, append = TRUE, sep="")

        cat(prefix, "parameters", suff, " = ", inla.paste(hyper[[k]]$param), "\n", file = file, append = TRUE, sep="")
        cat(prefix, "to.theta",   suff, " = ", inla.function2source(hyper[[k]]$to.theta), "\n", file = file, append = TRUE, sep="")
        cat(prefix, "from.theta",   suff, " = ", inla.function2source(hyper[[k]]$from.theta), "\n", file = file, append = TRUE, sep="")
    }

    return ()
}

`inla.write.boolean.field` = function(tag, val, file)
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

`inla.family.section` = function(...) {
    ## this is just a wrapper to make the naming better
    return (inla.data.section(...))
}
    
`inla.data.section` = function(file, family, file.data, file.weights, control, i.family="")
{
    ## this function is called from 'inla.family.section' only.
    cat("[INLA.Data", i.family, "]\n", sep = "", file = file,  append = TRUE)
    cat("type = data\n", sep = " ", file = file,  append = TRUE)
    cat("likelihood = ", family,"\n", sep = " ", file = file,  append = TRUE)
    cat("filename = ", file.data,"\n", sep = " ", file = file,  append = TRUE)
    cat("weights = ", file.weights,"\n", sep = " ", file = file,  append = TRUE)
    cat("link = ", inla.model.validate.link.function(family, control$link), "\n", 
        file = file,  append = TRUE)

    cat("variant = ",
        inla.ifelse(is.null(control$variant), 0L, as.integer(control$variant)),
        "\n", file = file,  append = TRUE)

    if (inla.one.of(family, "laplace")) {
        ## two parameters, alpha and epsilon is require for LAPLACE
        cat("alpha = ", inla.ifelse(is.null(control$alpha), 0.5, control$alpha), "\n",
            sep="", file=file, append=TRUE)
        cat("epsilon = ", inla.ifelse(is.null(control$epsilon), 0.01, control$epsilon), "\n",
            sep="", file=file, append=TRUE)
        cat("gamma = ", inla.ifelse(is.null(control$gamma), 1.0, control$gamma), "\n",
            sep="", file=file, append=TRUE)
    }

    if (inla.one.of(family, c("sn", "skewnormal"))) {
        cat("sn.shape.max = ", inla.ifelse(is.null(control$sn.shape.max), 5.0, control$sn.shape.max), "\n",
            sep="", file=file, append=TRUE)
    }

    if (inla.one.of(family, "gev")) {
        cat("gev.scale.xi = ", inla.ifelse(is.null(control$gev.scale.xi), 0.01, control$gev.scale.xi), "\n",
            sep="", file=file, append=TRUE)
    }
    
    inla.write.hyper(control$hyper, file)
    
    cat("\n", sep = " ", file = file,  append = TRUE)
}

`inla.ffield.section` = function(file, file.loc, file.cov, file.id.names = NULL,  n, nrep, ngroup, 
        file.extraconstr, file.weights, random.spec, results.dir, only.hyperparam, data.dir)
{
    label= inla.namefix(random.spec$term)
    prop = inla.model.properties(random.spec$model, "latent", stop.on.error=TRUE)
    
    cat("[", label,"]\n", sep = "", file = file,  append = TRUE)
    cat("dir = ", results.dir,"\n", sep = " ", file = file,  append = TRUE)
    cat("type = ffield\n", sep = " ", file = file,  append = TRUE)
    cat("model = ", random.spec$model,"\n", sep = " ", file = file,  append = TRUE)
    if (!is.null(random.spec$same.as)) {
        cat("same.as = ", random.spec$same.as,"\n", sep = " ", file = file,  append = TRUE)
    }
    cat("covariates = ", file.cov,"\n", sep = " ", file = file,  append = TRUE)
    if (!is.null(file.id.names)) {
        cat("id.names =", file.id.names,"\n", sep = " ", file = file,  append = TRUE)
    }
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
    if (!is.null(random.spec$spde2.prefix)) {
        ## need a special one, as spde2.prefix is not a file or a directory...
        fnm = inla.copy.dir.for.section.spde(random.spec$spde2.prefix, data.dir)
        cat("spde2.prefix =", fnm, "\n", sep = " ", file = file,  append = TRUE)
        cat("spde2.transform =", random.spec$spde2.transform, "\n", sep = " ", file = file,  append = TRUE)
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

    inla.write.hyper(random.spec$hyper, file)

    if (inla.model.properties(random.spec$model, "latent")$nrow.ncol) {
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
        cat("n = ", n,"\n", sep = " ", file = file,  append = TRUE)
    }
    cat("nrep = ", inla.ifelse(is.null(nrep), 1, nrep), "\n", sep = " ", file = file,  append = TRUE)

    if (!is.null(ngroup) && ngroup > 1) {
        cat("ngroup = ", ngroup, "\n", sep = " ", file = file,  append = TRUE)
        if (!is.null(random.spec$control.group$model))
            cat("group.model = ", random.spec$control.group$model, "\n", sep = " ", file = file,  append = TRUE)
        inla.write.hyper(random.spec$control.group$hyper, file = file,  prefix = "group.")
    }
        
    if (!is.null(random.spec$cyclic)) {
        cat("cyclic = ", as.numeric(random.spec$cyclic),"\n", sep = " ", file = file,  append = TRUE)
    }
    if (!is.null(random.spec$season.length)) {
        cat("season = ", random.spec$season.length,"\n", sep = " ", file = file,  append = TRUE)
    }
    if (!is.null(random.spec$graph)) {
        gfile = inla.write.graph(random.spec$graph, filename = inla.tempfile())
        fnm = inla.copy.file.for.section(gfile, data.dir)
        unlink(gfile)
        cat("graph = ", fnm, "\n", sep = " ", file = file,  append = TRUE)
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
            inla.write.fmesher.file(random.spec$Cmatrix, filename = file.C)
            file.C = gsub(data.dir, "$inladatadir", file.C, fixed=TRUE)
            cat("Cmatrix = ", file.C, "\n", append=TRUE, sep = " ", file = file)
        }
    }
    if (!is.null(random.spec$rankdef)) {
        cat("rankdef = ", random.spec$rankdef,"\n", append=TRUE, sep = " ", file = file)
    }
    if (!is.null(random.spec$cdf)) {
        cat("cdf = ", random.spec$cdf, "\n", sep = " ", file = file,  append = TRUE)
    }
    if (!is.null(random.spec$quantiles)) {
        cat("quantiles = ", random.spec$quantiles, "\n", sep = " ", file = file,  append = TRUE)
    }
    if (only.hyperparam || !random.spec$compute) {
        cat("compute = 0\n", sep = " ", file = file,  append = TRUE)
    } else {
        cat("compute = 1\n", sep = " ", file = file,  append = TRUE)
    }

    if (random.spec$model == "me") {
        ## possible scale-variable
        if (!is.null(random.spec$scale)) {
            file.scale=inla.tempfile(tmpdir=data.dir)
            ns = length(random.spec$scale)
            if (is.null(random.spec$values.order)) {
                idxs = 1:ns
            } else {
                idxs = random.spec$values.order
                stopifnot(length(random.spec$values.order) == ns)
            }
            inla.write.fmesher.file(as.matrix(cbind(idxs -1L, random.spec$scale)), filename=file.scale, debug = FALSE)
            ##print(cbind(idxs, random.spec$scale))
            file.scale = gsub(data.dir, "$inladatadir", file.scale, fixed=TRUE)
            cat("scale =", file.scale,"\n", sep = " ", file = file,  append = TRUE)
        }
    }

    if (random.spec$model == "rgeneric") {
        cat("rgeneric.Id = ", random.spec$rgeneric$Id, "\n", append=TRUE, sep = " ", file = file)
        stopifnot(!file.exists(random.spec$rgeneric$fifo$R2c))
        stopifnot(!file.exists(random.spec$rgeneric$fifo$c2R))
        cat("rgeneric.fifo.R2c = ", random.spec$rgeneric$fifo$R2c, "\n", append=TRUE, sep = " ", file = file)
        cat("rgeneric.fifo.c2R = ", random.spec$rgeneric$fifo$c2R, "\n", append=TRUE, sep = " ", file = file)
    }
            
    cat("\n", sep = " ", file = file,  append = TRUE)
}

`inla.inla.section` = function(file, inla.spec)
{
    cat("[INLA.Parameters]\n", sep = " ", file = file,  append = TRUE)
    cat("type = inla\n", sep = " ", file = file,  append = TRUE)

    if (!is.null(inla.spec$int.strategy)) {
        cat("int.strategy = ", inla.spec$int.strategy,"\n", sep = " ", file = file,  append = TRUE)
    }
    if (!is.null(inla.spec$strategy)) {
        cat("strategy = ", inla.spec$strategy,"\n", sep = " ", file = file,  append = TRUE)
    }
    if (!is.null(inla.spec$linear.correction)) {
        cat("linear.correction = ", inla.spec$linear.correction,"\n", sep = " ", file = file,  append = TRUE)
    }
    if (!is.null(inla.spec$h)) {
        cat("h = ", inla.spec$h,"\n", sep = " ", file = file,  append = TRUE)
    }
    if (!is.null(inla.spec$dz)) {
        cat("dz = ", inla.spec$dz,"\n", sep = " ", file = file,  append = TRUE)
    }
    if (!is.null(inla.spec$interpolator)) {
        cat("interpolator = ", inla.spec$interpolator,"\n", sep = " ", file = file,  append = TRUE)
    }
    if (!is.null(inla.spec$diff.logdens)) {
        cat("diff.log.dens = ", inla.spec$diff.logdens,"\n", sep = " ", file = file,  append = TRUE)
    }
    if (!is.null(inla.spec$print.joint.hyper)) {
        cat("fp.hyperparam = $inlaresdir/joint.dat\n", sep = "", file = file,  append = TRUE)
    }

    if (is.null(inla.spec$tolerance) || is.na(inla.spec$tolerance)) {
        ## we need a default value
        inla.spec$tolerance = 0.005
    }

    if (is.null(inla.spec$tolerance.f) || is.na(inla.spec$tolerance.f)) {
        inla.spec$tolerance.f = inla.spec$tolerance^(3/2) ## yes. 
    }
    cat("tolerance.f = ", inla.spec$tolerance.f,"\n", sep = " ", file = file,  append = TRUE)

    if (is.null(inla.spec$tolerance.g) || is.na(inla.spec$tolerance.g)) {
        inla.spec$tolerance.g = inla.spec$tolerance
    }
    cat("tolerance.g = ", inla.spec$tolerance.g,"\n", sep = " ", file = file,  append = TRUE)

    if (is.null(inla.spec$tolerance.x) || is.na(inla.spec$tolerance.x)) {
        inla.spec$tolerance.x = inla.spec$tolerance
    }
    cat("tolerance.x = ", inla.spec$tolerance.x,"\n", sep = " ", file = file,  append = TRUE)

    inla.write.boolean.field("hessian.force.diagonal", inla.spec$force.diagonal, file)
    inla.write.boolean.field("skip.configurations", inla.spec$skip.configurations, file)
    inla.write.boolean.field("mode.known", inla.spec$mode.known.conf, file)
    inla.write.boolean.field("adjust.weights", inla.spec$adjust.weights, file)
    inla.write.boolean.field("lincomb.derived.only", inla.spec$lincomb.derived.only, file)
    inla.write.boolean.field("lincomb.derived.correlation.matrix", inla.spec$lincomb.derived.correlation.matrix, file)

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
        ## reordering could be a number -1, 0, ....  or a string,  or the output from inla.qreordering()
        r = inla.spec$reordering
        if (is.list(r)) {
            ## output from inla.qreordering
            r = r$name
        }
        if (is.character(r)) {
            r.code = inla.reorderings.name2code(r)
        } else if (is.integer(r)) {
            ## this will fail is code is wrong
            dummy = inla.reorderings.code2name(r) 
            r.code = r
        } else {
            stop("This should not happen.")
        }
        cat("reordering = ", r.code, "\n", sep = " ", file = file,  append = TRUE)
    }

    if (!is.null(inla.spec$cpo.diff)) {
        cat("cpo.diff = ", inla.spec$cpo.diff, "\n", sep = " ", file = file,  append = TRUE)
    }
    if (!is.null(inla.spec$npoints)) {
        cat("n.points = ", inla.spec$npoints, "\n", sep = " ", file = file,  append = TRUE)
    }
    if (!is.null(inla.spec$cutoff)) {
        cat("cutoff = ", inla.spec$cutoff, "\n", sep = " ", file = file,  append = TRUE)
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
    if (!is.null(inla.spec$cmin)) {
        cat("cmin = ", inla.spec$cmin, "\n", file = file, append = TRUE)
    }
    if (!is.null(inla.spec$step.factor)) {
        cat("nr.step.factor = ", inla.spec$step.factor, "\n", file = file, append = TRUE)
    }
    if (!is.null(inla.spec$global.node.factor)) {
        cat("global.node.factor = ", inla.spec$global.node.factor, "\n", file = file, append = TRUE)
    }
    if (!is.null(inla.spec$global.node.degree)) {
        cat("global.node.degree = ", inla.spec$global.node.degree, "\n", file = file, append = TRUE)
    }

    ## options related to 'stupid search'.
    inla.write.boolean.field("stupid.search", inla.spec$stupid.search, file)
    if (!is.null(inla.spec$stupid.search.max.iter)) {
        cat("stupid.search.max.iter = ", as.integer(inla.spec$stupid.search.max.iter), "\n", file = file,  append = TRUE)
    }
    if (!is.null(inla.spec$stupid.search.factor)) {
        fac = as.numeric(inla.spec$stupid.search.factor)
        stopifnot(fac >= 1.0)
        cat("stupid.search.factor = ", fac, "\n", file = file,  append = TRUE)
    }

    cat("\n", sep = " ", file = file,  append = TRUE)
}

`inla.predictor.section` = function(file, n, m, predictor.spec, file.offset, data.dir, file.link.fitted.values)
{
    ## n = NPredictor
    ## m = MPredictor

    cat("[Predictor]\n", sep = " ", file = file,  append = TRUE)
    cat("type = predictor\n", sep = " ", file = file,  append = TRUE)
    cat("dir = predictor\n", sep = " ", file = file, append = TRUE)
    cat("n = ", n, "\n", sep = " ", file = file,  append = TRUE)
    cat("m = ", m, "\n", sep = " ", file = file,  append = TRUE)

    inla.write.boolean.field("fixed", predictor.spec$fixed, file)
    inla.write.boolean.field("compute", predictor.spec$compute, file)
    
    if (!is.null(predictor.spec$cdf)) {
        cat("cdf = ", predictor.spec$cdf, "\n", sep = " ", file = file,  append = TRUE)
    }
    if (!is.null(predictor.spec$quantiles)) {
        cat("cdf = ", predictor.spec$quantiles, "\n", sep = " ", file = file,  append = TRUE)
    }
    if (!is.null(file.offset)) {
        cat("offset = ", file.offset,"\n", sep = " ", file = file, append=TRUE)
    }
    if (!is.null(file.link.fitted.values)) {
        cat("link.fitted.values = ", file.link.fitted.values,"\n", sep = " ", file = file, append=TRUE)
    }

    inla.write.hyper(predictor.spec$hyper, file)

    if (!is.null(predictor.spec$cross) && length(predictor.spec$cross) > 0) {
        if (length(predictor.spec$cross) != n + m) {
            stop(paste("Length of cross does not match the total length of predictor", length(predictor.spec$cross), "!=", n+m))
        }
        file.cross = inla.tempfile(tmpdir=data.dir)
        ## better to go through factor to get levels 1...ncross. 
        cross = as.factor(predictor.spec$cross)
        cross = as.integer(cross)
        cross[is.na(cross)] = 0L ## means not in use
        if (inla.getOption("internal.binary.mode")) {
            inla.write.fmesher.file(as.matrix(cross, ncol=1), filename=file.cross)
        } else {
            write(cross, ncol=1, file=file.cross)
        }

        fnm = gsub(data.dir, "$inladatadir", file.cross, fixed=TRUE)
        cat("cross.constraint =", fnm, "\n", file=file, append = TRUE)
    }

    if (!is.null(predictor.spec$A)) {
        ## Now we will build the extended Matrix, which is
        ##
        ## Aextended = [ I, -A; -A^T, A^T A ] ((n+m) x (n+m))
        ##
        ## This matrix is the one that is needed for input to inla. 
    
        if (is.character(predictor.spec$A)) {
            A = read.table(predictor.spec$A, col.names = c("i", "j", "x"))
            A = sparseMatrix(i = A$i, j = A$j, x = A$x, index1=TRUE)
        } else {
            A = predictor.spec$A
        }
        A = inla.sparse.check(A, must.be.squared=FALSE)

        ## check dimensions
        stopifnot(dim(A)[1] == m)
        stopifnot(dim(A)[2] == n)

        ## replace NA's with zeros.
        A[ is.na(A) ] = 0.0

        ## Aext = [ I, -A; -A^T, A^T A ] ((n+m) x (n+m))
        Aext = rBind(cBind(Diagonal(m), -A), cBind(-t(A), t(A) %*% A))
        stopifnot(dim(Aext)[1] == m+n)
        stopifnot(dim(Aext)[2] == m+n)

        file.A=inla.tempfile(tmpdir=data.dir)
        inla.write.fmesher.file(Aext, filename = file.A)

        file.A = gsub(data.dir, "$inladatadir", file.A, fixed=TRUE)
        cat("Aext = ", file.A, "\n", append=TRUE, sep = " ", file = file)
        cat("AextPrecision = ", predictor.spec$precision, "\n", append=TRUE, sep = " ", file = file)
    }

    cat("\n", sep = " ", file = file,  append = TRUE)
}

`inla.problem.section` = function(file , data.dir, result.dir, hyperpar, return.marginals, dic,
        cpo, mlik, quantiles, smtp, q, strategy, graph)
{
    cat("", sep = "", file = file, append=FALSE)
    cat("###  ", inla.version("hgid"), "\n", sep = "", file = file,  append = TRUE) 
    cat("###  ", inla.paste(Sys.info()), "\n", sep = "", file = file,  append = TRUE) 
    cat("###  ", inla.os.type(), "-", inla.os.32or64bit(), "bit", " ", date(), "\n", sep = "", file = file,  append = TRUE) 
    cat("inladatadir = ", data.dir, "\n", sep = "", file = file,  append = TRUE)
    cat("inlaresdir = ", result.dir, "\n", sep = "", file = file,  append = TRUE)
    cat("#inladatadir = ", gsub("^.*/","", data.dir), "\n", sep = "", file = file,  append = TRUE) #
    cat("#inlaresdir = ", gsub("^.*/","", result.dir), "-%d\n", sep = "", file = file,  append = TRUE) #
 

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
    inla.write.boolean.field("graph", graph, file)

    if (!is.null(smtp)) {
        cat("smtp = ", smtp, "\n", sep = " ", file = file,  append = TRUE)
    }
    if (!is.null(quantiles)) {
        cat("quantiles = ", quantiles, "\n", sep = " ", file = file,  append = TRUE)
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

`inla.linear.section` = function(file, file.fixed, label, results.dir, control.fixed, only.hyperparam)
{
    cat("[", inla.namefix(label), "]\n", sep = "", file = file,  append = TRUE)
    cat("dir = ", results.dir,"\n", sep = " ", file = file,  append = TRUE)
    cat("type = linear\n", sep = " ", file = file,  append = TRUE)
    cat("covariates = ", file.fixed,"\n", sep = " ", file = file,  append = TRUE)
    if (only.hyperparam || !control.fixed$compute) {
        cat("compute = 0\n", sep = " ", file = file,  append = TRUE)
    } else {
        cat("compute = 1\n", sep = " ", file = file,  append = TRUE)
    }
    if (!is.null(control.fixed$cdf)) {
        cat("cdf = ", control.fixed$cdf, "\n", sep = " ", file = file,  append = TRUE)
    }
    if (!is.null(control.fixed$quantiles)) {
        cat("cdf = ", control.fixed$quantiles, "\n", sep = " ", file = file,  append = TRUE)
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

`inla.mode.section` = function(file, args, data.dir)
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
        inla.write.boolean.field("fixed", args$fixed, file)
    }
}

`inla.expert.section` = function(file, args)
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

`inla.lincomb.section` = function(file, data.dir, contr, lincomb)
{
    ## this one write binary format files...

    ## format is either
    ##
    ##     list("lc1" = "a 1 1 b 2 1 3 2...", ...)
    ##
    ## or
    ##
    ##     list("lc1" = list( "a" = list(idx=1, weight=1), "b" = list(idx=c(2, 3), weight = c(1, 2)), ...), ...)
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

`inla.copy.file.for.section` = function(filename, data.dir)
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

`inla.copy.dir.for.section` = function(dir.name, data.dir)
{
    d.fnm = inla.tempfile(tmpdir=data.dir)
    inla.dir.create(d.fnm)
    files.to.copy = paste(dir.name, "/", dir(dir.name, recursive=TRUE), sep="")
    file.copy(files.to.copy, d.fnm, recursive=TRUE)
    return (gsub(data.dir, "$inladatadir", d.fnm, fixed=TRUE))
}

`inla.copy.dir.for.section.spde` = function(prefix, data.dir)
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

`inla.z.section` = function(file, random.spec, data.dir, results.dir, only.hyperparam, k.off)
{
    ## binary io is not yet implemented. I think this model is on its way out...

    label= inla.namefix(random.spec$term)
    if (!is.matrix(random.spec$Z)) {
        stop("Argument Z in model=[z]has to be a matrix.")
    }
    n = dim(random.spec$Z)[1L]
    m = dim(random.spec$Z)[2L]
    ind = 0L:(n-1L)
    
    for(k in 1L:m) {
        kk = k + k.off
        cat("[", label,".", k, "]\n", sep = "", file = file,  append = TRUE)
        cat("type = ffield\n", sep = " ", file = file,  append = TRUE)
        cat("dir = ", "random.effect", inla.num(kk), "\n", sep = "", file = file, append = TRUE)
        cat("model = ", inla.ifelse(k == 1, "z", "zadd"),"\n", sep = " ", file = file,  append = TRUE)
        cat("n = 1\n", file=file, append = TRUE)
        if (k == 1L) {
            inla.write.hyper(random.spec$hyper, file)
        }

        file.cov=inla.tempfile(tmpdir=data.dir)
        file.create(file.cov)
        write(t(cbind(ind, rep(0, n))), ncolumns=2, file=file.cov, append=FALSE)
        file.cov = gsub(data.dir, "$inladatadir", file.cov, fixed=TRUE)
        cat("covariates = ", file.cov,"\n", sep = " ", file = file,  append = TRUE)

        file.w=inla.tempfile(tmpdir=data.dir)
        file.create(file.w)
        write(t(cbind(ind, random.spec$Z[, k])), ncolumns=2, file=file.w, append=FALSE)
        file.w = gsub(data.dir, "$inladatadir", file.w, fixed=TRUE)
        cat("weights = ", file.w,"\n", sep = " ", file = file,  append = TRUE)

        if (only.hyperparam) {
            cat("compute = 0\n", sep = " ", file = file,  append = TRUE)
        }
        cat("\n", sep = " ", file = file,  append = TRUE)
    }
}
