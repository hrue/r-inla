### Functions to write the different sections in the .ini-file

`inla.write.boolean.field` =
    function(tag, val, file)
{
    ## write tag = 1 or tag = 0 depending on val. if val is NULL do not write
    if (!is.null(val)) {
        if (val)
            cat(tag, " = 1\n", sep = " ", file = file, append = TRUE)
        else
            cat(tag, " = 0\n", sep = " ", file = file, append = TRUE)
    }
}

`inla.data.section` =
    function(file, family, file.data, control, i.family="")
{
    cat("[Data", i.family, "]\n", sep = "", file = file,  append = TRUE)
    cat("type = data\n", sep = " ", file = file,  append = TRUE)
    cat("likelihood = ",family,"\n", sep = " ", file = file,  append = TRUE)
    cat("filename = ", file.data,"\n", sep = " ", file = file,  append = TRUE)

    if (inla.one.of(family,"laplace")) {
        ## two parameters, alpha and epsilon is require for LAPLACE
        cat("alpha = ", ifelse(is.null(control$alpha), 0.5, control$alpha), "\n",
            sep="", file=file, append=TRUE)
        cat("epsilon = ", ifelse(is.null(control$epsilon), 0.01, control$epsilon), "\n",
            sep="", file=file, append=TRUE)
        cat("gamma = ", ifelse(is.null(control$gamma), 1.0, control$gamma), "\n",
            sep="", file=file, append=TRUE)
    }

    prop = inla.lmodel.properties(family)
    if (prop$ntheta > 0) {
        k = 1
        for(j in 1:prop$ntheta) {
            jj = j-1
            if (!is.null(control$fixed[j]))
                cat("fixed", ifelse(prop$ntheta==1, "", as.character(jj)), " = ",
                    ifelse(control$fixed[j], 1, 0),"\n", sep = "", file = file,  append = TRUE)
            if (!is.null(control$initial[j]))
                cat("initial", ifelse(prop$ntheta==1, "", as.character(jj)), " = ",
                    control$initial[j],"\n", sep = "", file = file,  append = TRUE)
        }
    }
    if (prop$npriors > 0) {
        k = 1
        for(j in 1:prop$npriors) {
            jj = j-1
            if (!is.null(control$prior) && !is.null(control$prior[j]))
                cat("prior", ifelse(prop$npriors==1, "", as.character(jj)), " = ",
                    control$prior[j],"\n", sep = "", file = file,  append = TRUE)
            if(!is.null(control$param)) {
                cat("parameters", ifelse(prop$npriors == 1, "", as.character(jj)), " = ", sep="", file = file, append=TRUE)
                ## divide equally
                for(k in 1:as.integer(round(prop$nparameters/prop$npriors))) {
                    cat(" ", control$param[k], sep = " ", file = file,  append = TRUE)
                    k = k + 1
                }
                cat("\n", file = file,  append = TRUE)
            }
        }
    }
    if(!is.null(control$dof.max))
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
        cat("covariates = ", file.cov,"\n", sep = " ", file = file,  append = TRUE)
        if (!is.null(random.spec$diagonal))
            cat("diagonal =", random.spec$diagonal,"\n", sep = " ", file = file,  append = TRUE)
        inla.write.boolean.field("constraint", random.spec$constr, file)
        inla.write.boolean.field("si", random.spec$si, file)
        if(!is.null(file.extraconstr))
            cat("extraconstraint =", file.extraconstr,"\n", sep = " ", file = file,  append = TRUE)
        if(!is.null(random.spec$weights)) 
            cat("weights =", file.weights,"\n", sep = " ", file = file,  append = TRUE)
        if (!is.null(random.spec$sphere.dir)) {
            fnm = inla.copy.dir.for.section(random.spec$sphere.dir, data.dir)
            cat("sphere.dir =", fnm, "\n", sep = " ", file = file,  append = TRUE)
        }
        if (!is.null(random.spec$T.order))
            cat("T.order =", random.spec$T.order, "\n", sep = " ", file = file,  append = TRUE)
        if (!is.null(random.spec$T.model))
            cat("T.model =", random.spec$T.model, "\n", sep = " ", file = file,  append = TRUE)
        if (!is.null(random.spec$K.order))
            cat("K.order =", random.spec$K.order, "\n", sep = " ", file = file,  append = TRUE)
        if (!is.null(random.spec$K.model))
            cat("K.model =", random.spec$K.model, "\n", sep = " ", file = file,  append = TRUE)

        if (!is.null(random.spec$of))
            cat("of =", random.spec$of, "\n", sep = " ", file = file,  append = TRUE)
        if (!is.null(random.spec$precision))
            cat("precision =", random.spec$precision, "\n", sep = " ", file = file,  append = TRUE)


        if (prop$ntheta > 0) {
            k = 1
            for(j in 1:prop$ntheta) {
                jj = j-1
                if (!is.null(random.spec$fixed[j]))
                    cat("fixed", ifelse(prop$ntheta==1, "", as.character(jj)), " = ",
                        ifelse(random.spec$fixed[j], 1, 0),"\n", sep = "", file = file,  append = TRUE)
                if (!is.null(random.spec$initial[j]))
                    cat("initial", ifelse(prop$ntheta==1, "", as.character(jj)), " = ",
                        random.spec$initial[j],"\n", sep = "", file = file,  append = TRUE)
            }
        }
        if (prop$npriors > 0) {
            k = 1
            for(j in 1:prop$npriors) {
                jj = j-1
                if (!is.null(random.spec$prior) && !is.null(random.spec$prior[j]))
                    cat("prior", ifelse(prop$npriors==1, "", as.character(jj)), " = ",
                        random.spec$prior[j], "\n", sep = "", file = file,  append = TRUE)
                if(!is.null(random.spec$param)) {
                    cat("parameters", ifelse(prop$npriors == 1, "", as.character(jj)), " = ", sep="", file = file, append=TRUE)
                    ## divide equally
                    for(kk in 1:as.integer(round(prop$nparameters/prop$npriors))) {
                        cat(" ", random.spec$param[k], sep = " ", file = file,  append = TRUE)
                        k = k + 1
                    }
                    cat("\n", file = file,  append = TRUE)
                }
            }
        }

        if (inla.model.properties(random.spec$model)$nrow.ncol) {
            cat("nrow = ", random.spec$nrow, "\n", sep = " ", file = file,  append = TRUE)
            cat("ncol = ", random.spec$ncol, "\n", sep = " ", file = file,  append = TRUE)

            if (!is.null(random.spec$bvalue))
                cat("bvalue = ", random.spec$bvalue, "\n", sep = " ", file = file,  append = TRUE)
            
            if (inla.one.of(random.spec$model, c("matern2d", "matern2dx2part0", "matern2dx2p1"))) {
                if (!is.null(random.spec$nu))
                    cat("nu = ", random.spec$nu, "\n", sep = " ", file = file,  append = TRUE)
            }
        }
        else
            cat("n = ",n,"\n", sep = " ", file = file,  append = TRUE)

        cat("nrep = ", inla.ifelse(is.null(nrep), 1, nrep), "\n", sep = " ", file = file,  append = TRUE)

        if (!is.null(ngroup) && ngroup > 1) {
            cat("ngroup = ", ngroup, "\n", sep = " ", file = file,  append = TRUE)
            if (!is.null(random.spec$control.group$prior))
                cat("group.prior = ", random.spec$control.group$prior, "\n", sep = " ", file = file,  append = TRUE)
            if (!is.null(random.spec$control.group$param))
                cat("group.parameters = ", random.spec$control.group$param, "\n", sep = " ", file = file,  append = TRUE)
            if (!is.null(random.spec$control.group$initial))
                cat("group.initial = ", random.spec$control.group$initial, "\n", sep = " ", file = file,  append = TRUE)
            if (!is.null(random.spec$control.group$fixed)) {
                if (random.spec$control.group$fixed)
                    cat("group.fixed = ", 1, "\n", sep = " ", file = file,  append = TRUE)
                else 
                    cat("group.fixed = ", 0, "\n", sep = " ", file = file,  append = TRUE)
            }
            if (!is.null(random.spec$control.group$model)) {
                cat("group.model = ", random.spec$control.group$model, "\n", sep = " ", file = file,  append = TRUE)
            }
        }
        
        if(!is.null(random.spec$cyclic))
            cat("cyclic = ",as.numeric(random.spec$cyclic),"\n", sep = " ", file = file,  append = TRUE)
        if(!is.null(random.spec$season.length))
            cat("season = ",random.spec$season.length,"\n", sep = " ", file = file,  append = TRUE)
        if(!is.null(random.spec$graph.file)) {
            fnm = inla.copy.file.for.section(random.spec$graph.file, data.dir)
            cat("graph = ",fnm, "\n", sep = " ", file = file,  append = TRUE)
        }
        if (!is.null(file.loc))
            cat("locations = ", file.loc,"\n", sep = " ", file = file,  append = TRUE)

        if (!is.null(random.spec$Cmatrix)) {
            if (is.character(random.spec$Cmatrix)) {
                cat("Cmatrix = ",random.spec$Cmatrix,"\n",append=TRUE, sep = " ", file = file)
            } else {
                file.C=inla.tempfile(tmpdir=data.dir)
                file.create(file.C)

                ii = sort(unique(random.spec$Cmatrix$i))
                ni = length(ii)
                if (!all(ii == 1:ni))
                    stop(paste("Cmatrix$i is not the numbers 1, 2, ...,", ni))

                jj = sort(unique(random.spec$Cmatrix$j))
                nj = length(jj)
                if (!all(jj == 1:nj))
                    stop(paste("Cmatrix$j is not the numbers 1, 2, ...,", nj))

                if (ni != nj || length(random.spec$Cmatrix$i) != length(random.spec$Cmatrix$Cij) ||
                    length(random.spec$Cmatrix$j) != length(random.spec$Cmatrix$i)) {
                    stop(paste("Wrong dimensions:",
                               "length(sort(unique(random.spec$Cmatrix$i)))", ni,
                               ", length(sort(unique(random.spec$Cmatrix$j)))", nj,
                               ", length(random.spec$Cmatrix$i)", length(random.spec$Cmatrix$i),
                               ", length(random.spec$Cmatrix$j)", length(random.spec$Cmatrix$j),
                               ", length(random.spec$Cmatrix$Cij)", length(random.spec$Cmatrix$Cij)))
                }
                write(t(cbind(random.spec$Cmatrix$i -1 , random.spec$Cmatrix$j - 1, random.spec$Cmatrix$Cij)),
                      ncolumns=3,file=file.C,append=FALSE)
                file.C = gsub(data.dir, "$DATADIR", file.C, fixed=TRUE)
                cat("Cmatrix = ",file.C, "\n",append=TRUE, sep = " ", file = file)
            }
        }
        if(!is.null(random.spec$rankdef))
            cat("rankdef = ",random.spec$rankdef,"\n",append=TRUE, sep = " ", file = file)
        if(!is.null(random.spec$cdf))
            cat("cdf = ",random.spec$cdf, "\n", sep = " ", file = file,  append = TRUE)
        if(!is.null(random.spec$quantiles))
            cat("quantiles = ",random.spec$quantiles, "\n", sep = " ", file = file,  append = TRUE)
        if(only.hyperparam)
            cat("compute=0\n", sep = " ", file = file,  append = TRUE)
        cat("\n", sep = " ", file = file,  append = TRUE)
    }
}

`inla.offset.section` =
    function(file, file.loc, n, nrep, file.cov, file.extraconstr, file.weights, random.spec, results.dir, only.hyperparam, data.dir)
{
    ## special version of inla.ffield.section() but for type = offset
    ##
    label= inla.namefix(random.spec$term)
    prop = inla.model.properties(random.spec$model, stop.on.error=TRUE)
    
    cat("[",label,"]\n", sep = "", file = file,  append = TRUE)
    cat("dir = ",results.dir,"\n", sep = " ", file = file,  append = TRUE)
    cat("type = offset\n", sep = " ", file = file,  append = TRUE)
    cat("model = ",random.spec$model,"\n", sep = " ", file = file,  append = TRUE)
    if(!is.null(random.spec$weights))
        cat("weights =", file.weights,"\n", sep = " ", file = file,  append = TRUE)

    if (prop$ntheta > 0) {
        k = 1
        for(j in 1:prop$ntheta) {
            jj = j-1
            if (!is.null(random.spec$fixed[j]))
                cat("fixed", ifelse(prop$ntheta==1, "", as.character(jj)), " = ",
                    ifelse(random.spec$fixed[j], 1, 0),"\n", sep = "", file = file,  append = TRUE)
            if (!is.null(random.spec$initial[j]))
                cat("initial", ifelse(prop$ntheta==1, "", as.character(jj)), " = ",
                    random.spec$initial[j],"\n", sep = "", file = file,  append = TRUE)
        }
    }
    if (prop$npriors > 0) {
        k = 1
        for(j in 1:prop$npriors)
         
            jj = j-1
        if (!is.null(random.spec$prior) && !is.null(random.spec$prior[j]))
            cat("prior", ifelse(prop$npriors==1, "", as.character(jj)), " = ",
                random.spec$prior[j], "\n", sep = "", file = file,  append = TRUE)
        if(!is.null(random.spec$param)) {
            cat("parameters", ifelse(prop$npriors == 1, "", as.character(jj)), " = ", sep="", file = file, append=TRUE)
            ## divide equally
            for(kk in 1:as.integer(round(prop$nparameters/prop$npriors))) {
                cat(" ", random.spec$param[k], sep = " ", file = file,  append = TRUE)
                k = k + 1
            }
            cat("\n", file = file,  append = TRUE)
        }
    }

    ##cat("n = ",n,"\n", sep = " ", file = file,  append = TRUE)

    if(!is.null(random.spec$cdf))
        cat("cdf = ",random.spec$cdf, "\n", sep = " ", file = file,  append = TRUE)
    if(!is.null(random.spec$quantiles))
        cat("quantiles = ",random.spec$quantiles, "\n", sep = " ", file = file,  append = TRUE)
    
    cat("\n", sep = " ", file = file,  append = TRUE)
}

`inla.inla.section` =
    function(file,inla.spec)
{
    cat("[INLA param]\n", sep = " ", file = file,  append = TRUE)
    cat("type = inla\n", sep = " ", file = file,  append = TRUE)

    if(!is.null(inla.spec$int.strategy)) 
        cat("int.strategy = ",inla.spec$int.strategy,"\n", sep = " ", file = file,  append = TRUE)
    if(!is.null(inla.spec$strategy)) 
        cat("strategy = ",inla.spec$strategy,"\n", sep = " ", file = file,  append = TRUE)
    if(!is.null(inla.spec$linear.correction)) 
        cat("linear.correction = ",inla.spec$linear.correction,"\n", sep = " ", file = file,  append = TRUE)
    if(!is.null(inla.spec$h)) 
        cat("h = ",inla.spec$h,"\n", sep = " ", file = file,  append = TRUE)
    if(!is.null(inla.spec$dz)) 
        cat("dz = ",inla.spec$dz,"\n", sep = " ", file = file,  append = TRUE)
    if(!is.null(inla.spec$interpolator)) 
        cat("interpolator = ", inla.spec$interpolator,"\n", sep = " ", file = file,  append = TRUE)
    if(!is.null(inla.spec$diff.logdens)) 
        cat("diff.log.dens = ",inla.spec$diff.logdens,"\n", sep = " ", file = file,  append = TRUE)
    if(inla.spec$print.joint.hyper)
        cat("fp.hyperparam = $RESDIR/joint.dat\n", sep = " ", file = file,  append = TRUE)
    if(!is.null(inla.spec$tolerance)) {
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
    inla.write.boolean.field("huge", inla.spec$huge, file)
    inla.write.boolean.field("derived.only", inla.spec$derived.only, file)

    if (!is.null(inla.spec$restart) && inla.spec$restart >= 0)
        cat("restart = ", as.integer(inla.spec$restart), "\n", file = file, sep = " ", append = TRUE)

    if (!is.null(inla.spec$optimiser))
        cat("optimiser = ", inla.spec$optimiser,"\n", sep = " ", file = file,  append = TRUE)
    if (!is.null(inla.spec$verbose) && inla.spec$verbose)
        cat("optpar.fp = stdout\n", sep = " ", file = file,  append = TRUE)

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

    if (!is.null(inla.spec$cpo.diff)) 
        cat("cpo.diff = ", inla.spec$cpo.diff, "\n", sep = " ", file = file,  append = TRUE)

    if (!is.null(inla.spec$npoints))
        cat("npoints = ", inla.spec$npoints, "\n", sep = " ", file = file,  append = TRUE)

    inla.write.boolean.field("adapt.hessian.mode", inla.spec$adapt.hessian.mode, file)

    if (!is.null(inla.spec$adapt.hessian.max.trials) && inla.spec$adapt.hessian.max.trials >= 0)
        cat("adapt.hessian.max.trials = ", inla.spec$adapt.hessian.max.trials, "\n", sep = " ", file = file,  append = TRUE)

    if (!is.null(inla.spec$adapt.hessian.scale) && inla.spec$adapt.hessian.scale >= 1)
        cat("adapt.hessian.scale = ", inla.spec$adapt.hessian.scale, "\n", sep = " ", file = file,  append = TRUE)

    if (!is.null(inla.spec$step.len) && inla.spec$step.len > 0)
        cat("step.len = ", inla.spec$step.len, "\n", sep = " ", file = file,  append = TRUE)

    if (!is.null(inla.spec$diagonal) && inla.spec$diagonal >= 0.0)
        cat("diagonal = ", inla.spec$diagonal, "\n", sep = " ", file = file,  append = TRUE)
        
    cat("\n", sep = " ", file = file,  append = TRUE)
}

`inla.predictor.section` =
    function(file, n, predictor.spec, file.offset, data.dir)
{
    cat("[Predictor]\n", sep = " ", file = file,  append = TRUE)
    cat("type = predictor\n", sep = " ", file = file,  append = TRUE)
    cat("dir = predictor\n", sep = " ", file = file, append = TRUE)
    cat("n = ",n,"\n", sep = " ", file = file,  append = TRUE)

    inla.write.boolean.field("fixed", predictor.spec$fixed, file)
    inla.write.boolean.field("compute", predictor.spec$compute, file)
    
    if(!is.null(predictor.spec$cdf))
        cat("cdf = ",predictor.spec$cdf, "\n", sep = " ", file = file,  append = TRUE)
    if(!is.null(predictor.spec$quantiles))
        cat("cdf = ",predictor.spec$quantiles, "\n", sep = " ", file = file,  append = TRUE)
    if(!is.null(file.offset))
        cat("offset = ", file.offset,"\n", sep = " ", file = file, append=TRUE)
    if(!is.null(predictor.spec$param))
        cat("parameters = ",predictor.spec$param[1]," ",predictor.spec$param[2],"\n", sep = " ", file = file, append=TRUE)
    if(!is.null(predictor.spec$initial)) {
        cat("initial = ",predictor.spec$initial,"\n", sep = " ", file = file, append=TRUE)
    } else {
        if(predictor.spec$fixed) 
            cat("initial = 12\n", sep = " ", file = file,  append = TRUE)
    }
    if(!is.null(predictor.spec$predictor.usermap)){
        cat("predictor.usermap=", predictor.spec$predictor.usermap, "\n", sep=" ", file = file, append = TRUE)
    }
    if (!is.null(predictor.spec$cross) && length(predictor.spec$cross) > 0) {
        if (length(predictor.spec$cross) != n)
            stop(paste("Length of cross does not match length of predictor", length(predictor.spec$cross), "!=", n))
        file.cross = inla.tempfile(tmpdir=data.dir)
        write(predictor.spec$cross, ncol=1, file=file.cross)
        fnm = gsub(data.dir, "$DATADIR", file.cross, fixed=TRUE)
        cat("cross.constraint =", fnm, "\n", file=file, append = TRUE)
    }

    cat("\n", sep = " ", file = file,  append = TRUE)
}

`inla.problem.section` =
    function(file , data.dir, result.dir, hyperpar, dic, cpo, mlik, quantiles, smtp, q)
{
    cat("", sep = "", file = file, append=FALSE)
    cat("DATADIR = ",data.dir, "\n", sep = "", file = file,  append = TRUE)
    cat("RESDIR = ",result.dir, "\n", sep = "", file = file,  append = TRUE)
    cat("#DATADIR = ",gsub("^.*/","",data.dir), "\n", sep = "", file = file,  append = TRUE) #
    cat("#RESDIR = ", gsub("^.*/","",result.dir), "-%d\n", sep = "", file = file,  append = TRUE) #

    cat("\n", sep = " ", file = file,  append = TRUE)
    cat("[Model]\n", sep = " ", file = file,  append = TRUE)
    cat("type = problem\n", sep = " ", file = file,  append = TRUE)
    cat("dir = $RESDIR\n", sep = " ", file = file,  append = TRUE)
    inla.write.boolean.field("hyperparameters",hyperpar, file)
    inla.write.boolean.field("cpo", cpo, file)
    inla.write.boolean.field("dic", dic, file)
    inla.write.boolean.field("mlik", mlik, file)
    inla.write.boolean.field("q", q, file)

    if(!is.null(smtp))
        cat("smtp = ",smtp, "\n", sep = " ", file = file,  append = TRUE)
    if(!is.null(quantiles))
        cat("quantiles = ",quantiles, "\n", sep = " ", file = file,  append = TRUE)
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
    function(file, file.fixed,label,results.dir,control.fixed,only.hyperparam)
{
    cat("[",inla.namefix(label),"]\n", sep = "", file = file,  append = TRUE)
    cat("dir = ",results.dir,"\n", sep = " ", file = file,  append = TRUE)
    cat("type = linear\n", sep = " ", file = file,  append = TRUE)
    cat("covariates = ", file.fixed,"\n", sep = " ", file = file,  append = TRUE)
    if(only.hyperparam) 
        cat("compute = 0\n", sep = " ", file = file,  append = TRUE)
    if(!is.null(control.fixed$cdf))
        cat("cdf = ",control.fixed$cdf, "\n", sep = " ", file = file,  append = TRUE)
    if(!is.null(control.fixed$quantiles))
        cat("cdf = ",control.fixed$quantiles, "\n", sep = " ", file = file,  append = TRUE)
    if (length(grep("^[(]Intercept[)]$", inla.trim(label))) == 1) {
        prec = control.fixed$prec.intercept
        mean = control.fixed$mean.intercept
        if(!is.null(mean))
            cat("mean = ", mean, "\n", sep = " ", file = file,  append = TRUE)
        if(!is.null(prec))
            cat("precision = ", prec, "\n", sep = " ", file = file,  append = TRUE)
    } else {
        prec = inla.parse.fixed.prior(label, control.fixed$prec)
        mean = inla.parse.fixed.prior(label, control.fixed$mean)
        if(!is.null(mean))
            cat("mean = ", mean, "\n", sep = " ", file = file,  append = TRUE)
        if(!is.null(prec))
            cat("precision = ", prec, "\n", sep = " ", file = file,  append = TRUE)
    }
    cat("\n", sep = " ", file = file,  append = TRUE)
}

`inla.mode.section` =
    function(file, args, data.dir)
{
    if (!is.null(args$mode) && length(args$mode) > 0) {
        cat("[Mode]\n", sep = " ", file = file,  append = TRUE)
        cat("type = mode\n", sep = " ", file = file,  append = TRUE)

        cat("mode = ", inla.paste(as.character(args$mode)), "\n", sep = " ", file = file,  append = TRUE)
        if (!is.null(args$x.mode) && length(args$x.mode) > 0) {
            file.x.mode = inla.tempfile(tmpdir=data.dir)
            write(args$x.mode, ncol=1, file=file.x.mode)
            fnm = gsub(data.dir, "$DATADIR", file.x.mode, fixed=TRUE)
            cat("xmode =", fnm, "\n", file=file, append = TRUE)
        }
        inla.write.boolean.field("restart", args$restart, file)
    }
}

`inla.expert.section` =
    function(file, args)
{
    if (!is.null(args$cpo.manual) && args$cpo.manual) {
        cat("\n## If you edit this section it is assumed you know what you're doing ;-)\n", file=file, append=TRUE) #
        cat("[Expert]\n", sep = " ", file = file,  append = TRUE)
        cat("type = expert\n", sep = " ", file = file,  append = TRUE)
        inla.write.boolean.field("cpo.manual", args$cpo.manual, file)
        ## recall to convert to 0-based index'ing
        cat("cpo.idx = ", args$cpo.idx -1,"\n", sep = " ", file = file,  append = TRUE)
    }
}

`inla.lincomb.section` =
    function(file, data.dir, contr, lincomb)
{
    if (!is.null(lincomb)) {
        for(i in 1:length(lincomb)) {
            if (is.null(names(lincomb[i])) || is.na(names(lincomb[i]))) {
                secname = paste("lincomb.", i, sep="")
                lc = lincomb[i]
            } else if (names(lincomb[i]) == "") {
                secname = paste("lincomb.", i, sep="")
                lc = lincomb[[i]]
            } else {
                secname = paste("lincomb.", names(lincomb[i])[1], sep="")
                lc = lincomb[[i]]
            }
                
            cat("\n[", secname, "]\n", sep = "", file = file,  append = TRUE)
            cat("type = lincomb\n", sep = " ", file = file,  append = TRUE)
            inla.write.boolean.field("c.indexing", contr$c.indexing, file)
            if (!is.null(contr$precision))
                cat("precision = ", contr$precision,"\n", sep = " ", file = file,  append = TRUE)
            if (!is.null(contr$usermap))
                cat("usermap = ", contr$usermap,"\n", sep = " ", file = file,  append = TRUE)
            fnm = inla.tempfile(tmpdir=data.dir)
            file.create(fnm)
            write(lc, file=fnm)
            fnm = gsub(data.dir, "$DATADIR", fnm, fixed=TRUE)
            cat("filename = ", fnm, "\n", sep = " ", file = file, append = TRUE)
        }
    }
}

`inla.copy.file.for.section` =
    function(filename, data.dir)
{
    if (missing(filename))
        return (NULL)
    if (missing(data.dir))
        stop("data.dir required")

    fnm = inla.tempfile(tmpdir=data.dir)
    file.copy(filename, fnm, overwrite=TRUE)
    return (gsub(data.dir, "$DATADIR", fnm, fixed=TRUE))
}

`inla.copy.dir.for.section` =
    function(dirname, data.dir)
{
    if (missing(dirname))
        return (NULL)
    if (missing(data.dir))
        stop("data.dir required")

    d.fnm = inla.tempfile(tmpdir=data.dir)
    dir.create(d.fnm, recursive=TRUE)
    files.to.copy = paste(dirname, "/", dir(dirname, recursive=TRUE), sep="")
    file.copy(files.to.copy, d.fnm, recursive=TRUE)
    return (gsub(data.dir, "$DATADIR", d.fnm, fixed=TRUE))
}

`inla.z.section` =
    function(file, random.spec, data.dir, results.dir, only.hyperparam, k.off)
{
    label= inla.namefix(random.spec$term)
    if (!is.matrix(random.spec$Z))
        stop("Argument Z in model=[z]has to be a matrix.")
    n = dim(random.spec$Z)[1]
    m = dim(random.spec$Z)[2]
    ind = 0:(n-1)
    
    for(k in 1:m) {
        kk = k + k.off
        cat("[",label,".", k, "]\n", sep = "", file = file,  append = TRUE)
        cat("type = ffield\n", sep = " ", file = file,  append = TRUE)
        cat("dir = ", "random.effect", inla.num(kk), "\n", sep = "", file = file, append = TRUE)
        cat("model = ", ifelse(k == 1, "z", "zadd"),"\n", sep = " ", file = file,  append = TRUE)
        cat("n = 1\n", file=file, append = TRUE)
        if (k == 1) {
            if(!is.null(random.spec$param))
                cat("parameters = ",random.spec$param[1]," ",random.spec$param[2],"\n", sep = " ", file = file,
                    append = TRUE)
            cat("prior = ",random.spec$prior.logprec,"\n", sep = " ", file = file,  append = TRUE)
        }

        file.cov=inla.tempfile(tmpdir=data.dir)
        file.create(file.cov)
        write(t(cbind(ind,rep(0,n))),ncolumns=2,file=file.cov,append=FALSE)
        file.cov = gsub(data.dir, "$DATADIR", file.cov, fixed=TRUE)
        cat("covariates = ", file.cov,"\n", sep = " ", file = file,  append = TRUE)

        file.w=inla.tempfile(tmpdir=data.dir)
        file.create(file.w)
        write(t(cbind(ind,random.spec$Z[,k])),ncolumns=2,file=file.w,append=FALSE)
        file.w = gsub(data.dir, "$DATADIR", file.w, fixed=TRUE)
        cat("weights = ", file.w,"\n", sep = " ", file = file,  append = TRUE)

        if(only.hyperparam) 
            cat("compute = 0\n", sep = " ", file = file,  append = TRUE)
        cat("\n", sep = " ", file = file,  append = TRUE)
    }
}

        

    
