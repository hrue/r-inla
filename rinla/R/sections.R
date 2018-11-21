## Nothing to export

### Functions to write the different sections in the .ini-file

`inla.secsep` = function(secname) 
{
    sep = "!"
    if (missing(secname)) {
        return (sep)
    } else {
        return (paste0(sep, inla.namefix(secname), sep))
    }
}
    
`inla.write.hyper` = function(hyper, file, prefix="", data.dir, ngroup = -1L, low = -Inf, high = Inf)
{
    stopifnot(!missing(hyper))
    stopifnot(!missing(file))

    if (is.null(hyper) || length(hyper) == 0L) {
        return ()
    }
    len = length(hyper)
    for(k in 1L:len) {
        if (len == 1L) {
            suff = ""
        } else {
            suff = as.character(k-1L)
        }
        cat(prefix, "initial", suff, " = ", hyper[[k]]$initial, "\n", file = file, append = TRUE, sep="")
        cat(prefix, "fixed", suff, " = ", as.numeric(hyper[[k]]$fixed), "\n", file = file, append = TRUE, sep="")
        cat(prefix, "hyperid", suff, " = ", as.numeric(hyper[[k]]$hyperid), "\n", file = file, append = TRUE, sep="")

        ## these are for "expression:"...
        ## if there are newlines,  remove them
        tmp.prior = gsub("\n", "", hyper[[k]]$prior)
        ## remove preceding spaces
        tmp.prior = gsub("^[ \t]+", "", tmp.prior)
        ## if the expression ends with a ";" with or without spaces, remove it
        tmp.prior = gsub(";*[ \t]*$", "", tmp.prior)
        ## if there is a 'return (' replace it with 'return('
        tmp.prior = gsub("return[ \t]+\\(", "return(", tmp.prior)
        ## for all priors except the "expression:" one,  then trim the name
        if (length(grep("^(expression|table)[ \t]*:", tolower(tmp.prior))) == 0L) {
            tmp.prior = inla.trim.family(tmp.prior)
        }

        ## table: is now stored in a file
        if (length(grep("^table:",  tmp.prior)) > 0) {
            tab = substr(tmp.prior, nchar("table:")+1, nchar(tmp.prior))
            xy = as.numeric(unlist(strsplit(tab, "[ \t\n\r]+")))
            xy = xy[!is.na(xy)]
            nxy = length(xy) %/% 2L
            xx = xy[1:nxy]
            yy = xy[1:nxy + nxy]
            xy = cbind(xx, yy)
            file.xy = inla.tempfile(tmpdir=data.dir)
            inla.write.fmesher.file(xy, filename = file.xy)
            file.xy = gsub(data.dir, "$inladatadir", file.xy, fixed=TRUE)
            cat(prefix, "prior", suff, " = table: ", file.xy, "\n", append=TRUE, sep = "", file = file)
        } else {
            cat(prefix, "prior", suff, " = ", tmp.prior, "\n", file = file, append = TRUE, sep="")
        }

        cat(prefix, "parameters", suff, " = ", inla.paste(hyper[[k]]$param), "\n", file = file, append = TRUE, sep="")
        to.t = gsub("REPLACE.ME.ngroup", paste("ngroup=", as.integer(ngroup), sep=""), inla.function2source(hyper[[k]]$to.theta))
        from.t = gsub("REPLACE.ME.ngroup", paste("ngroup=", as.integer(ngroup), sep=""), inla.function2source(hyper[[k]]$from.theta))
        to.t = gsub("REPLACE.ME.low", paste("low=", as.numeric(low), sep=""), to.t)
        from.t = gsub("REPLACE.ME.low", paste("low=", as.numeric(low), sep=""), from.t)
        to.t = gsub("REPLACE.ME.high", paste("high=", as.numeric(high), sep=""), to.t)
        from.t = gsub("REPLACE.ME.high", paste("high=", as.numeric(high), sep=""), from.t)
        cat(prefix, "to.theta", suff, " = ", to.t, "\n", file = file, append = TRUE, sep="")
        cat(prefix, "from.theta", suff, " = ", from.t, "\n", file = file, append = TRUE, sep="")

        ## do a second replacement so that we replace the functions with actual functions after
        ## the replacement of REPLACE.ME.....
        hyper[[k]]$from.theta = eval(parse(text = gsub("REPLACE.ME.ngroup", paste("ngroup=", as.integer(ngroup), sep=""),
                                               inla.function2source(hyper[[k]]$from.theta, newline = "
"))))
        hyper[[k]]$from.theta = eval(parse(text = gsub("REPLACE.ME.low", paste("low=", as.numeric(low), sep=""),
                                               inla.function2source(hyper[[k]]$from.theta, newline = "
"))))
        hyper[[k]]$from.theta = eval(parse(text = gsub("REPLACE.ME.high", paste("high=", as.numeric(high), sep=""),
                                               inla.function2source(hyper[[k]]$from.theta, newline = "
"))))
        hyper[[k]]$to.theta = eval(parse(text= gsub("REPLACE.ME.low", paste("low=", as.numeric(low), sep=""),
                                             inla.function2source(hyper[[k]]$to.theta, newline = "
"))))
        hyper[[k]]$to.theta = eval(parse(text= gsub("REPLACE.ME.high", paste("high=", as.numeric(high), sep=""),
                                             inla.function2source(hyper[[k]]$to.theta, newline = "
"))))
        hyper[[k]]$to.theta = eval(parse(text= gsub("REPLACE.ME.ngroup", paste("ngroup=", as.integer(ngroup), sep=""),
                                             inla.function2source(hyper[[k]]$to.theta, newline = "
"))))
    }

    return (hyper)
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

`inla.family.section` = function(...)
{
    ## this is just a wrapper to make the naming better
    return (inla.data.section(...))
}

`inla.data.section` = function(
        file, family, file.data, file.weights, control, i.family="",
        link.covariates = link.covariates, data.dir)
{
    ## this function is called from 'inla.family.section' only.
    cat(inla.secsep(), "INLA.Data", i.family, inla.secsep(), "\n", sep = "", file = file,  append = TRUE)
    cat("type = data\n", sep = " ", file = file,  append = TRUE)
    cat("likelihood = ", family,"\n", sep = " ", file = file,  append = TRUE)
    cat("filename = ", file.data,"\n", sep = " ", file = file,  append = TRUE)
    cat("weights = ", file.weights,"\n", sep = " ", file = file,  append = TRUE)

    cat("variant = ",
        inla.ifelse(is.null(control$variant), 0L, as.integer(control$variant)),
        "\n", file = file,  append = TRUE)

    if (inla.one.of(family, "cenpoisson")) {
        if (is.null(control$cenpoisson.I)) {
            interval = inla.set.control.family.default()$cenpoisson.I
        } else {
            ## must be integers
            interval = round(control$cenpoisson.I)
        }
        if (length(interval) != 2L) {
            stop(paste("cenpoisson.I: Must be a vector of length 2.", length(interval)))
        }
        if (interval[1] > interval[2]) {
            stop(paste("cenpoisson.I = c(Low, High): Low > High!", interval[1],  interval[2]))
        }
        cat("cenpoisson.I = ", interval[1], " ",  interval[2], "\n", sep="", file=file, append=TRUE)
    }
    
    if (TRUE) {
        if (!is.null(control$quantile))
            stop("control.family=list(quantile=...) is disabled. Use control.family=list(control.link=list(quantile=...)) instead")
        quantile = control$control.link$quantile
        if (is.numeric(quantile)) {
            if ((quantile <= 0.0) || (quantile >= 1.0)) {
                stop(paste("quantile: Must be a numeric in the interval (0, 1)"))
            }
        } else {
            quantile = -1  ## so we get an error if used.
        }
        cat("quantile = ", quantile, "\n", sep="", file=file, append=TRUE)
    }

    if (inla.one.of(family, c("sn", "skewnormal"))) {
        cat("sn.shape.max = ", inla.ifelse(is.null(control$sn.shape.max), 5.0, control$sn.shape.max), "\n",
            sep="", file=file, append=TRUE)
    }

    if (inla.one.of(family, "gev")) {
        cat("gev.scale.xi = ", inla.ifelse(is.null(control$gev.scale.xi), 0.01, control$gev.scale.xi), "\n",
            sep="", file=file, append=TRUE)
    }

    inla.write.hyper(control$hyper, file, data.dir = data.dir)

    ## the link-part. first make it backward-compatible...
    if (!(is.null(control$link) || inla.strcasecmp(control$link, "default"))) {
        ## control$link is set,  use that if not control.link$model is set
        if (!(is.null(control$control.link$model) || inla.strcasecmp(control$control.link$model, "default"))) {
            stop("Both control.family$link (OBSOLETE) and control.family$control.link$model is set. Please use the latter only.")
        }
        control$control.link$model = control$link
    }
    lmod = control$control.link$model = inla.model.validate.link.function(family, control$control.link$model)
    ord = control$control.link$order
    cat("link.model = ", lmod, "\n", file = file,  append = TRUE)
    if (inla.one.of(lmod, c("special1"))) {
        ## for these models, the argument order is required
        if (is.null(ord)) {
            stop(paste("For link-model:", lmod, ", 'order' must be spesified."))
        }
        cat("link.order = ", ord, "\n", file = file,  append = TRUE)
        ## special
        if (ord > 1L && length(control$control.link$hyper$theta2$param) == 2L) {
            par = control$control.link$hyper$theta2$param
            control$control.link$hyper$theta2$param = c(rep(par[1], ord), c(diag(rep(par[2], ord))))
        }
    } else {
        if (!is.null(ord)) {
            stop(paste("For link-model:", lmod, ", the argument 'order' is not used and must be NULL."))
        }
    }

    variant = control$control.link$variant
    if (inla.one.of(lmod, c("logoffset"))) {
        ## for these models, the argument 'variant' is required
        if (## general
            is.null(variant) ||
            ## model-spesific
            (inla.one.of(lmod, c("logoffset")) && all(variant != c(0, 1)))) {
            stop(paste("For link-model:", lmod, ", the argument variant must be 0 or 1,  not", variant))
        }
        cat("link.variant = ", variant, "\n", file = file,  append = TRUE)
    } else {
        if (!is.null(variant)) {
            stop(paste("For link-model:", lmod, ", the argument 'variant' is not used and must be NULL."))
        }
    }

    inla.write.hyper(control$control.link$hyper, file, prefix = "link.", data.dir = dirname(file))
    if (!is.null(link.covariates)) {
        if (!is.matrix(link.covariates)) {
            link.covariates = matrix(c(link.covariates), ncol = 1)
        }
        file.link.cov = inla.tempfile(tmpdir=data.dir)
        inla.write.fmesher.file(link.covariates, filename = file.link.cov)
        file.link.cov = gsub(data.dir, "$inladatadir", file.link.cov, fixed=TRUE)
        cat("link.covariates = ", file.link.cov, "\n", append=TRUE, sep = " ", file = file)
    }

    ## the mix-part
    inla.write.boolean.field("mix.use", !is.null(control$control.mix$model), file)
    if (!is.null(control$control.mix$model)) {
        cat("mix.model = ", control$control.mix$model, "\n", sep="", file=file, append=TRUE)
        nq = as.integer(control$control.mix$nq)
        stopifnot(nq >= 5L)
        cat("mix.nq = ", nq, "\n", sep="", file=file, append=TRUE)
        integrator = match.arg(control$control.mix$integrator, c("default", "gq", "simpson"))
        cat("mix.integrator = ", integrator, "\n", sep="", file=file, append=TRUE)
        inla.write.hyper(control$control.mix$hyper, file, prefix = "mix.", data.dir = dirname(file))
    }

    cat("\n", sep = " ", file = file,  append = TRUE)
}

`inla.ffield.section` = function(file, file.loc, file.cov, file.id.names = NULL,  n, nrep, ngroup,
        file.extraconstr, file.weights, random.spec, results.dir, only.hyperparam, data.dir)
{
    label= random.spec$term
    prop = inla.model.properties(random.spec$model, "latent", stop.on.error=TRUE)

    cat(inla.secsep(label), "\n", sep = "", file = file,  append = TRUE)
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
    if (!is.null(random.spec$spde3.prefix)) {
        ## need a special one, as spde3.prefix is not a file or a directory...
        fnm = inla.copy.dir.for.section.spde(random.spec$spde3.prefix, data.dir)
        cat("spde3.prefix =", fnm, "\n", sep = " ", file = file,  append = TRUE)
        cat("spde3.transform =", random.spec$spde3.transform, "\n", sep = " ", file = file,  append = TRUE)
    }
    if (inla.one.of(random.spec$model, "copy")) {
        if (!is.null(random.spec$of)) {
            cat("of =", random.spec$of, "\n", sep = " ", file = file,  append = TRUE)
        }
    }
    if (inla.one.of(random.spec$model, c("copy", "sigm", "revsigm", "log1exp", "fgn"))) {
        if (!is.null(random.spec$precision)) {
            cat("precision =", random.spec$precision, "\n", sep = " ", file = file,  append = TRUE)
        }
    }

    if (inla.one.of(random.spec$model, c("clinear", "copy", "mec", "meb"))) {
        if (is.null(random.spec$range)) {
            ## default is the identity mapping
            random.spec$range = c(0, 0)
        }
        cat("range.low  =", random.spec$range[1], "\n", sep = " ", file = file, append = TRUE)
        cat("range.high =", random.spec$range[2], "\n", sep = " ", file = file, append = TRUE)
    }
    if (inla.one.of(random.spec$model, c("rw1", "rw2", "besag", "bym", "bym2", "besag2", "rw2d", "rw2diid", "seasonal"))) {
        if (is.null(random.spec$scale.model)) {
            random.spec$scale.model = inla.getOption("scale.model.default")
        }
        cat("scale.model = ", as.numeric(random.spec$scale.model), "\n", sep = " ", file = file,  append = TRUE)
    }
    if (inla.one.of(random.spec$model, c("besag", "bym", "bym2", "besag2"))) {
        cat("adjust.for.con.comp = ", as.numeric(random.spec$adjust.for.con.comp), "\n", sep = " ", file = file,  append = TRUE)
    }

    if (FALSE) {
        ## this is only for the mvnorm prior,  which we do not use anymore (hrue/16/02/2015)
        if (inla.one.of(random.spec$model, "ar")) {
            ## set a default prior for order > 1 if the param is given only for p=1
            par = random.spec$hyper$theta2$param
            if (length(par) == 2L && random.spec$order > 1L) {
                random.spec$hyper$theta2$param = c(rep(par[1], random.spec$order), par[2]*diag(random.spec$order))
            }
        }
    }

    ## possible adaptive priors
    if (inla.one.of(random.spec$model, c("bym2", "rw2diid"))) {
        stopifnot(random.spec$hyper$theta2$short.name == "phi") ## just a check...
        if (random.spec$hyper$theta2$prior == "pc") {
            ## this is the pc-prior which is adaptive for this
            ## model. compute the prior here.
            if (inla.one.of(random.spec$model, "bym2")) {
                random.spec$hyper$theta2$prior = inla.pc.bym.phi(
                    graph = random.spec$graph,
                    ## have to do this automatic
                    ## rankdef = random.spec$rankdef,
                    u = random.spec$hyper$theta2$param[1L],
                    alpha = random.spec$hyper$theta2$param[2L],
                    scale.model = TRUE,
                    return.as.table = TRUE,
                    adjust.for.con.comp = as.numeric(random.spec$adjust.for.con.comp))
                random.spec$hyper$theta2$param = numeric(0)
            } else if (inla.one.of(random.spec$model, "rw2diid")) {
                random.spec$hyper$theta2$prior = inla.pc.rw2diid.phi(
                    size = c(random.spec$nrow,  random.spec$ncol),
                    u = random.spec$hyper$theta2$param[1L],
                    alpha = random.spec$hyper$theta2$param[2L],
                    return.as.table = TRUE)
                random.spec$hyper$theta2$param = numeric(0)
            } else {
                stop("This should not happpen.")
            }
        }
    }

    if (is.null(random.spec$range)) {
        low = -Inf
        high = Inf
    } else {
        low = random.spec$range[1]
        high = random.spec$range[2]
    }
    random.spec$hyper = inla.write.hyper(random.spec$hyper, file, data.dir = data.dir, ngroup = ngroup, low = low, high = high)

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

        if (!inla.one.of(random.spec$model, "copy")) {
            cat("group.model = ", random.spec$control.group$model, "\n", sep = " ", file = file,  append = TRUE)
            inla.write.boolean.field("group.cyclic", random.spec$control.group$cyclic, file)
            if (is.null(random.spec$control.group$scale.model)) {
                random.spec$control.group$scale.model = inla.getOption("scale.model.default")
            }
            inla.write.boolean.field("group.scale.model", random.spec$control.group$scale.model, file)
            inla.write.boolean.field("group.adjust.for.con.comp", random.spec$control.group$adjust.for.con.comp, file)
            if (inla.one.of(random.spec$control.group$model, "ar")) {
                ## 'order' is only used for model=ar
                p = inla.ifelse(is.null(random.spec$control.group$order), 0, as.integer(random.spec$control.group$order))
                cat("group.order = ", p, "\n", sep = " ", file = file,  append = TRUE)
            }
            if (inla.one.of(random.spec$control.group$model, "besag")) {
                stopifnot(!is.null(random.spec$control.group$graph))
                gfile = inla.write.graph(random.spec$control.group$graph, filename = inla.tempfile())
                fnm = inla.copy.file.for.section(gfile, data.dir)
                unlink(gfile)
                cat("group.graph = ", fnm, "\n", sep = " ", file = file,  append = TRUE)
            } else {
                stopifnot(is.null(random.spec$control.group$graph))
            }
            random.spec$control.group$hyper = (inla.write.hyper(random.spec$control.group$hyper,
                                                                file = file,  prefix = "group.",
                                                                data.dir = data.dir, ngroup = ngroup))
        }
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

    if (inla.one.of(random.spec$model, "z")) {
        ## This is for the random-effect Z*z, where Z is a n x m
        ## matrix and z are N_m(0, prec*C). We rewrite this as a model
        ## for zz = (v, z), where v is N(Z*z, high.precision). This
        ## gives the precision matrix for zz as A+prec*B, where A and
        ## B are defined below.
        Z = inla.as.sparse(random.spec$Z)
        tZ = t(Z)
        Z.n = dim(Z)[1]
        Z.m = dim(Z)[2]
        A = inla.as.sparse(random.spec$precision * cbind(rbind(Diagonal(Z.n), -tZ), rbind(-Z, tZ %*% Z)))
        if (is.null(random.spec$Cmatrix)) {
            Cm = inla.as.sparse(Diagonal(Z.m))
        } else {
            Cm = inla.as.sparse(random.spec$Cmatrix)
        }
        stopifnot(all(Z.m == dim(Cm)))
        B = inla.as.sparse(cbind(
            rbind(Diagonal(Z.n, 0.0),                                     # n x n zero-matrix
                  sparseMatrix(dims = c(Z.m, Z.n), i = 1, j = 1, x = 0)), # m x n zero-matrix
            rbind(sparseMatrix(dims = c(Z.n, Z.m), i = 1, j = 1, x = 0),  # n x m zero-matrix
                  Cm)))

        ## dimensions
        cat("z.n = ", Z.n,"\n", append=TRUE, sep = " ", file = file)
        cat("z.m = ", Z.m,"\n", append=TRUE, sep = " ", file = file)
        ## matrix A
        file.A = inla.tempfile(tmpdir=data.dir)
        inla.write.fmesher.file(A, filename = file.A)
        file.A = gsub(data.dir, "$inladatadir", file.A, fixed=TRUE)
        cat("z.Amatrix = ", file.A, "\n", append=TRUE, sep = " ", file = file)
        ## matrix B
        file.B = inla.tempfile(tmpdir=data.dir)
        inla.write.fmesher.file(B, filename = file.B)
        file.B = gsub(data.dir, "$inladatadir", file.B, fixed=TRUE)
        cat("z.Bmatrix = ", file.B, "\n", append=TRUE, sep = " ", file = file)
    }

    if (inla.one.of(random.spec$model, "dmatern")) {
        ## need the matrix of locations
        file.loc = inla.tempfile(tmpdir=data.dir)
        inla.write.fmesher.file(random.spec$locations, filename = file.loc)
        file.loc = gsub(data.dir, "$inladatadir", file.loc, fixed=TRUE)
        cat("dmatern.locations = ", file.loc, "\n", append=TRUE, sep = " ", file = file)
    }

    if (inla.one.of(random.spec$model, "generic3")) {
        ## For this model, Cmatrix is a list of matrices. Error checking have be done already in
        ## f()
        nC = length(random.spec$Cmatrix)
        stopifnot(nC > 0L)
        cat("generic3.n = ", random.spec$n, "\n", append=TRUE, sep = "", file = file)
        cat("generic3.m = ", nC, "\n", append=TRUE, sep = "", file = file)
        for(k in 1L:nC) {
            file.A = inla.tempfile(tmpdir=data.dir)
            inla.write.fmesher.file(inla.as.sparse(random.spec$Cmatrix[[k]]), filename = file.A)
            file.A = gsub(data.dir, "$inladatadir", file.A, fixed=TRUE)
            cat("generic3.Cmatrix.", as.integer(k-1L), " = ", file.A, "\n", append=TRUE, sep = "", file = file)
        }
    }

    ## if the Cmatrix is defined we need to process it except if its
    ## the z/generic3-model for which this has already been done.
    if (!inla.one.of(random.spec$model, c("z", "generic3")) && !is.null(random.spec$Cmatrix)) {
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

    if (inla.one.of(random.spec$model, "slm")) {
        ## This is the spatial-lag-model (SLM),  with
        ## N_C(b,  Q = kappa*A1 + A2 + rho*kappa*B + rho^2*kappa*C)
        X = random.spec$args.slm$X
        W = random.spec$args.slm$W
        Q = random.spec$args.slm$Q.beta
        slm.n = dim(X)[1L]
        slm.m = dim(X)[2L]

        cat("slm.n = ", slm.n,"\n", append=TRUE, sep = " ", file = file)
        cat("slm.m = ", slm.m,"\n", append=TRUE, sep = " ", file = file)
        cat("slm.rho.min = ", random.spec$args.slm$rho.min,"\n", append=TRUE, sep = " ", file = file)
        cat("slm.rho.max = ", random.spec$args.slm$rho.max,"\n", append=TRUE, sep = " ", file = file)

        ## matrix A1
        A1 = cbind(
            rbind(Diagonal(slm.n),
                  -t(X)),
            rbind(-X,
                  t(X) %*% X))
        file.A1 = inla.tempfile(tmpdir=data.dir)
        inla.write.fmesher.file(A1, filename = file.A1)
        file.A1 = gsub(data.dir, "$inladatadir", file.A1, fixed=TRUE)
        cat("slm.A1matrix = ", file.A1, "\n", append=TRUE, sep = " ", file = file)

        ## matrix A2
        A2 = cbind(
            rbind(Matrix(0, slm.n, slm.n),
                  Matrix(0, slm.m, slm.n)),
            rbind(Matrix(0, slm.n, slm.m),
                  Q))
        file.A2 = inla.tempfile(tmpdir=data.dir)
        inla.write.fmesher.file(A2, filename = file.A2)
        file.A2 = gsub(data.dir, "$inladatadir", file.A2, fixed=TRUE)
        cat("slm.A2matrix = ", file.A2, "\n", append=TRUE, sep = " ", file = file)

        ## matrix B
        B = cbind(
            rbind(-(t(W) + W),
                  t(X) %*% W),
            rbind(t(W) %*% X,
                  Matrix(0, slm.m, slm.m)))
        file.B = inla.tempfile(tmpdir=data.dir)
        inla.write.fmesher.file(B, filename = file.B)
        file.B = gsub(data.dir, "$inladatadir", file.B, fixed=TRUE)
        cat("slm.Bmatrix = ", file.B, "\n", append=TRUE, sep = " ", file = file)

        ## matrix C
        C = cbind(
            rbind(t(W) %*% W,
                  Matrix(0, slm.m, slm.n)),
            rbind(Matrix(0, slm.n, slm.m),
                  Matrix(0, slm.m, slm.m)))
        file.C = inla.tempfile(tmpdir=data.dir)
        inla.write.fmesher.file(C, filename = file.C)
        file.C = gsub(data.dir, "$inladatadir", file.C, fixed=TRUE)
        cat("slm.Cmatrix = ", file.C, "\n", append=TRUE, sep = " ", file = file)
    }

    if (inla.one.of(random.spec$model, "ar1c")) {
        Z = random.spec$args.ar1c$Z
        Q.beta = random.spec$args.ar1c$Q.beta
        ar1c.n = dim(Z)[1L]
        ar1c.m = dim(Z)[2L]

        cat("ar1c.n = ", ar1c.n,"\n", append=TRUE, sep = " ", file = file)
        cat("ar1c.m = ", ar1c.m,"\n", append=TRUE, sep = " ", file = file)

        ## matrix Z
        Z[ar1c.n, ] = 0  ## not used
        file.Z = inla.tempfile(tmpdir=data.dir)
        inla.write.fmesher.file(Z, filename = file.Z)
        file.Z = gsub(data.dir, "$inladatadir", file.Z, fixed=TRUE)
        cat("ar1c.Z = ", file.Z, "\n", append=TRUE, sep = " ", file = file)

        ## matrix t(Z) %*% Z (makes my life simpler)
        ZZ = t(Z) %*% Z
        file.ZZ = inla.tempfile(tmpdir=data.dir)
        inla.write.fmesher.file(ZZ, filename = file.ZZ)
        file.ZZ = gsub(data.dir, "$inladatadir", file.ZZ, fixed=TRUE)
        cat("ar1c.ZZ = ", file.ZZ, "\n", append=TRUE, sep = " ", file = file)

        file.Qbeta = inla.tempfile(tmpdir=data.dir)
        inla.write.fmesher.file(Q.beta, filename = file.Qbeta)
        file.Qbeta = gsub(data.dir, "$inladatadir", file.Qbeta, fixed=TRUE)
        cat("ar1c.Qbeta = ", file.Qbeta, "\n", append=TRUE, sep = " ", file = file)
    }

    ## the z-model for which this has already been done.
    if (!inla.one.of(random.spec$model, c("z", "generic3")) && !is.null(random.spec$Cmatrix)) {
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

    if (inla.one.of(random.spec$model, c("mec", "meb", "iid"))) {
        ## possible scale-variable
        if (!is.null(random.spec$scale)) {
            file.scale=inla.tempfile(tmpdir=data.dir)
            ns = length(random.spec$scale)
            if (is.null(random.spec$values)) {
                idxs = 1:ns
            } else {
                idxs = 1:ns
                stopifnot(length(random.spec$values) == ns)
                stopifnot(!is.null(idxs))
            }
            inla.write.fmesher.file(as.matrix(cbind(idxs -1L, random.spec$scale)), filename=file.scale, debug = FALSE)
            ##print(cbind(idxs, random.spec$scale))
            file.scale = gsub(data.dir, "$inladatadir", file.scale, fixed=TRUE)
            cat("scale =", file.scale,"\n", sep = " ", file = file,  append = TRUE)
        }
    } else {
        if (!is.null(random.spec$scale)) {
            stop(paste("Section [",  label, "]: option 'scale' is not NULL but not used.",  sep=""))
        }
    }

    if (random.spec$model == "rgeneric") {
        file.rgeneric = inla.tempfile(tmpdir=data.dir)
        ## this must be the same name as R_GENERIC_MODEL in inla.h
        model = paste(".inla.rgeneric.model", ".", random.spec$rgeneric$Id, sep="")
        assign(model, random.spec$rgeneric$model)
        ## save model, or the object that 'model' expands to
        inla.eval(paste("save(", model,
                        ", file = ", "\"", file.rgeneric, "\"",
                        ", ascii = FALSE, compress = TRUE)",  sep=""))
        fnm = gsub(data.dir, "$inladatadir", file.rgeneric, fixed=TRUE)
        cat("rgeneric.file =", fnm, "\n", file=file, append = TRUE)
        cat("rgeneric.model =", model, "\n", file=file, append = TRUE)
        rm(model) ## do not need it anymore
    }

    if (inla.one.of(random.spec$model, c("ar", "fgn", "fgn2"))) {
        cat("order = ", random.spec$order, "\n", append=TRUE, sep = " ", file = file)
    }

    if (is.null(random.spec$correct)) {
        random.spec$correct = -1L  ## code for ``make the default choice''
    }
    cat("correct = ", as.numeric(random.spec$correct), "\n", append=TRUE, sep = "", file = file)
    cat("\n", sep = " ", file = file,  append = TRUE)

    ## need to store the updated one
    return (random.spec)
}

`inla.inla.section` = function(file, inla.spec, data.dir)
{
    cat(inla.secsep("INLA.Parameters"), "\n", sep = " ", file = file,  append = TRUE)
    cat("type = inla\n", sep = " ", file = file,  append = TRUE)

    if (!is.null(inla.spec$int.strategy)) {
        cat("int.strategy = ", inla.spec$int.strategy,"\n", sep = " ", file = file,  append = TRUE)
    }

    if (inla.one.of(inla.spec$int.strategy, c("user", "user.std"))) {
        if (is.null(inla.spec$int.design)) {
            stop(paste0("int.strategy = 'user' or 'user.std' require the integration design in 'int.design'"))
        }
        file.A = inla.tempfile(tmpdir=data.dir)
        inla.write.fmesher.file(as.matrix(inla.spec$int.design), filename = file.A)
        file.A = gsub(data.dir, "$inladatadir", file.A, fixed=TRUE)
        cat("int.design = ", file.A, "\n", sep = " ", file = file,  append = TRUE)
    }

    if (!is.null(inla.spec$strategy)) {
        cat("strategy = ", inla.spec$strategy,"\n", sep = " ", file = file,  append = TRUE)
    }
    if (!is.null(inla.spec$adaptive.max)) {
        cat("adaptive.max = ", as.integer(inla.spec$adaptive.max),"\n", sep = " ", file = file,  append = TRUE)
    }
    inla.write.boolean.field("fast", inla.spec$fast, file)
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

    if (!(is.null(inla.spec$tolerance.step) || is.na(inla.spec$tolerance.step))) {
        cat("tolerance.step = ", inla.spec$tolerance.step,"\n", sep = " ", file = file, append = TRUE)
    }
        
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
    if (!is.null(inla.spec$stencil)) {
        stopifnot(inla.spec$stencil %in% c(3, 5, 7, 9))
        cat("stencil = ", inla.spec$stencil, "\n", sep = " ", file = file,  append = TRUE)
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

    inla.write.boolean.field("correct", inla.spec$correct, file)
    inla.write.boolean.field("correct.verbose", inla.spec$correct.verbose, file)
    if (!is.null(inla.spec$correct.factor)) {
        stopifnot(inla.spec$correct.factor > 0);
        cat("correct.factor = ", inla.spec$correct.factor, "\n", file = file, append = TRUE)
    }
    if (!is.null(inla.spec$correct.strategy)) {
        cat("correct.strategy = ", inla.spec$correct.strategy,"\n", sep = " ", file = file,  append = TRUE)
    }

    cat("\n", sep = " ", file = file,  append = TRUE)
}

`inla.predictor.section` = function(file, n, m, predictor.spec, file.offset, data.dir, file.link.fitted.values)
{
    ## n = NPredictor
    ## m = MPredictor

    cat(inla.secsep("Predictor"), "\n", sep = " ", file = file,  append = TRUE)
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
        cat("quantiles = ", predictor.spec$quantiles, "\n", sep = " ", file = file,  append = TRUE)
    }
    if (!is.null(file.offset)) {
        cat("offset = ", file.offset,"\n", sep = " ", file = file, append=TRUE)
    }
    if (!is.null(file.link.fitted.values)) {
        cat("link.fitted.values = ", file.link.fitted.values,"\n", sep = " ", file = file, append=TRUE)
    }

    inla.write.hyper(predictor.spec$hyper, file, data.dir = data.dir)

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
            inla.write.fmesher.file(as.matrix(cross, ncol=1L), filename=file.cross)
        } else {
            write(cross, ncolumns=1L, file=file.cross)
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

        ## replace NA's with zeros. (This is now done already in inla.R)
        ## A[ is.na(A) ] = 0.0

        ## Aext = [ I, -A; -A^T, A^T A ] ((n+m) x (n+m))
        Aext = rbind(cbind(Diagonal(m), -A), cbind(-t(A), t(A) %*% A))
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
        cpo, po, mlik, quantiles, smtp, q, openmp.strategy, graph, config, gdensity)
{
    cat("", sep = "", file = file, append=FALSE)
    cat("###  ", inla.version("hgid"), "\n", sep = "", file = file,  append = TRUE)
    cat("###  ", inla.paste(Sys.info()), "\n", sep = "", file = file,  append = TRUE)
    cat("###  ", inla.os.type(), "-", inla.os.32or64bit(), "bit", " ", date(), "\n", sep = "", file = file,  append = TRUE)

    cat("\n### [[[start of output from sessionInfo()]]]", "\n",  file = file, append = TRUE)
    s = paste("###   ", capture.output(sessionInfo()))
    for(i in seq_along(s)) {
        cat(s[i], "\n",  file = file, append = TRUE)
    }
    cat("### [[[end of output from sessionInfo()]]]\n", "\n",  file = file, append = TRUE)

    cat("inladatadir = ", data.dir, "\n", sep = "", file = file,  append = TRUE)
    cat("inlaresdir = ", result.dir, "\n", sep = "", file = file,  append = TRUE)
    cat("##inladatadir = ", gsub("^.*/","", data.dir), "\n", sep = "", file = file,  append = TRUE) #
    cat("##inlaresdir = ", gsub("^.*/","", result.dir), "-%d\n", sep = "", file = file,  append = TRUE) #

    ## libR-stuff
    cat("\n", sep = " ", file = file,  append = TRUE)
    cat(inla.secsep("INLA.libR"), "\n", sep = " ", file = file,  append = TRUE)
    cat("type = libR\n", sep = " ", file = file,  append = TRUE)
    cat("R_HOME = ", Sys.getenv("R_HOME"), "\n", sep = "", file = file,  append = TRUE)

    cat("\n", sep = " ", file = file,  append = TRUE)
    cat(inla.secsep("INLA.Model"), "\n", sep = " ", file = file,  append = TRUE)
    cat("type = problem\n", sep = " ", file = file,  append = TRUE)
    cat("dir = $inlaresdir\n", sep = " ", file = file,  append = TRUE)
    inla.write.boolean.field("return.marginals", return.marginals, file)
    inla.write.boolean.field("hyperparameters", hyperpar, file)
    inla.write.boolean.field("cpo", cpo, file)
    inla.write.boolean.field("po", po, file)
    inla.write.boolean.field("dic", dic, file)
    inla.write.boolean.field("mlik", mlik, file)
    inla.write.boolean.field("q", q, file)
    inla.write.boolean.field("graph", graph, file)
    inla.write.boolean.field("config", config, file)
    inla.write.boolean.field("gdensity", gdensity, file)

    if (is.null(smtp) || !(is.character(smtp) && (nchar(smtp) > 0))) {
        smtp = inla.getOption("smtp")
    }
    smtp = match.arg(tolower(smtp), c("band", "taucs", "pardiso", "default"))
    cat("smtp = ", smtp, "\n", sep = " ", file = file,  append = TRUE)

    if (is.null(openmp.strategy) || !(is.character(openmp.strategy) && (nchar(openmp.strategy) > 0))) {
        openmp.strategy = "default"
    }
    openmp.strategy = match.arg(tolower(openmp.strategy),
                                c("default", "small", "medium", "large", "huge", "pardiso.serial", "pardiso.parallel"))
    cat("openmp.strategy = ", openmp.strategy, "\n", sep = " ", file = file,  append = TRUE)

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
    } else if (any(name == names(prior))) {
        return (prior[[which(name == names(prior))]])
    } else if (any("default" == names(prior))) {
        return (prior[[which("default" == names(prior))]])
    } 
    return (NULL)
}

`inla.linear.section` = function(file, file.fixed, label, results.dir, control.fixed, only.hyperparam)
{
    cat(inla.secsep(label), "\n", sep = "", file = file,  append = TRUE)
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
        cat("quantiles = ", control.fixed$quantiles, "\n", sep = " ", file = file,  append = TRUE)
    }
    if (length(grep("^[(]Intercept[)]$", inla.trim(label))) == 1) {
        prec = control.fixed$prec.intercept
        mean = control.fixed$mean.intercept
        if (is.null(mean)) {
            mean =inla.set.control.fixed.default()$mean.intercept
        }
        if (is.null(prec)) {
            prec =inla.set.control.fixed.default()$prec.intercept
        }
        stopifnot(!is.null(mean) || !is.null(prec))
    } else {
        mean = inla.parse.fixed.prior(label, control.fixed$mean)
        prec = inla.parse.fixed.prior(label, control.fixed$prec)
        if (is.null(mean)) {
            mean = inla.set.control.fixed.default()$mean
        }
        if (is.null(prec)) {
            prec = inla.set.control.fixed.default()$prec
        }
        stopifnot(!is.null(mean) || !is.null(prec))
    }
    cat("mean = ", mean, "\n", sep = " ", file = file,  append = TRUE)
    cat("precision = ", prec, "\n", sep = " ", file = file,  append = TRUE)
    cat("\n", sep = " ", file = file,  append = TRUE)

    return (list(label = label, prior.mean = mean, prior.prec = prec))
}

`inla.mode.section` = function(file, args, data.dir)
{
    if (!is.null(args$result) || !is.null(args$theta) || !is.null(args$x)) {

        cat(inla.secsep("INLA.Control.Mode"), "\n", sep = " ", file = file,  append = TRUE)
        cat("type = mode\n", sep = " ", file = file,  append = TRUE)

        ## use default the mode in result if given
        if (is.null(args$theta) && !is.null(args$result)) {
            args$theta = args$result$mode$theta
        }

        if (!is.null(args$theta)) {
            ## set non-finite's to 0
            args$theta[!is.finite(args$theta)] = 0
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
            ## set non-finite's to 0
            args$x[!is.finite(args$x)] = 0
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

`inla.expert.section` = function(file, args, data.dir)
{
    cat(inla.secsep("INLA.Expert"), "\n", sep = " ", file = file,  append = TRUE)
    cat("type = expert\n", sep = " ", file = file,  append = TRUE)

    if (!is.null(args$cpo.manual) && args$cpo.manual) {
        inla.write.boolean.field("cpo.manual", args$cpo.manual, file)
        ## recall to convert to 0-based index'ing
        cat("cpo.idx = ", args$cpo.idx -1,"\n", sep = " ", file = file,  append = TRUE)
    }

    if (!is.null(args$jp)) {
        stopifnot(inherits(args$jp, "inla.jp"))
        file.jp = inla.tempfile(tmpdir=data.dir)
        ## this must be the same name as R_JP_MODEL in inla.h
        model = ".inla.jp.model"
        assign(model, args$jp$model)
        ## save model, or the object that 'model' expands to
        inla.eval(paste("save(", model,
                        ", file = ", "\"", file.jp, "\"",
                        ", ascii = FALSE, compress = TRUE)",  sep=""))
        fnm = gsub(data.dir, "$inladatadir", file.jp, fixed=TRUE)
        cat("jp.file =", fnm, "\n", file=file, append = TRUE)
        cat("jp.model =", model, "\n", file=file, append = TRUE)
        rm(model) ## do not need it anymore
    }

    if (is.null(args$disable.gaussian.check)) {
        args$disable.gaussian.check = FALSE
    }
    inla.write.boolean.field("DISABLE.GAUSSIAN.CHECK", args$disable.gaussian.check, file)
}

`inla.update.section` = function(file, data.dir, contr)
{
    if (!is.null(contr$result) && length(contr$result$mode$theta) > 0L)
    {
        file.update = inla.tempfile(tmpdir=data.dir)
        cat(inla.secsep("INLA.update"), "\n", sep = " ", file = file,  append = TRUE)
        cat("type = update\n", sep = " ", file = file,  append = TRUE)
        x = c(
            length(contr$result$mode$theta),
            contr$result$mode$theta,
            contr$result$misc$stdev.corr.positive,
            contr$result$misc$stdev.corr.negative,
            1/sqrt(contr$result$misc$cov.intern.eigenvalues), ## this is what is required
            c(contr$result$misc$cov.intern.eigenvectors)      ## column-wise storage
            )
        x = as.matrix(x, ncol = 1L)
        inla.write.fmesher.file(x, filename = file.update)
        file.update = gsub(data.dir, "$inladatadir", file.update, fixed=TRUE)
        cat("filename = ", file.update, "\n", sep = " ", file = file,  append = TRUE)
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
        prev.secnames = c()

        for(i in 1:length(lincomb)) {

            nam = names(lincomb[i])
            if (is.null(nam) || is.na(nam)) {
                secname = paste("lincomb.", inla.num(i, width=numlen), sep="")
                lc = lincomb[[i]]
            } else if (nam == "") {
                secname = paste("lincomb.", inla.num(i, width=numlen), sep="")
                lc = lincomb[[i]]
            } else {
                secname = paste("lincomb.", nam[1], sep="")
                lc = lincomb[[i]]
            }

            ## check if the name is used previously, if so, stop.
            if (secname %in% prev.secnames) {
                stop(paste("Duplicated name [", secname, "] in 'lincomb'; need unique names or NA or ''.",
                           sep=""))
            } else {
                prev.secnames = c(secname, prev.secnames)
            }

            cat("\n", inla.secsep(secname), "\n", sep = "", file = file,  append = TRUE)
            cat("type = lincomb\n", sep = " ", file = file,  append = TRUE)
            cat("lincomb.order = ", i, "\n", sep = " ", file = file,  append = TRUE)
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
    ## They original files were created by f(), so remove them:
    unlink(files.to.copy)
    return (rprefix)
}
