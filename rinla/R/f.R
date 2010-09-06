
`f` =
function(...,
         model = "iid",
         copy=NULL,
         n=NULL,
         nrep = NULL,
         replicate = NULL,
         ngroup = NULL,
         group = NULL,
         control.group = inla.set.control.group.default(),
         initial=NULL,
         season.length=NULL,
         prior=NULL,
         param = NULL,
         fixed = NULL,
         constr = NULL,
         extraconstr=NULL,
         values=NULL,
         cyclic = NULL,
         diagonal = NULL, 
         graph.file=NULL,
         cdf=NULL,
         quantiles=NULL,
         Qmatrix=NULL,  ### not used anymore
         Cmatrix=NULL,
         rankdef=NULL,
         Z = NULL,
         nrow = NULL,
         ncol = NULL,
         nu = NULL,
         bvalue = NULL,
         sphere.dir = NULL,
         spde.prefix = NULL,
         T.order=NULL,
         T.model=NULL,
         K.order=NULL,
         K.model=NULL,
         mean.linear=NULL,
         prec.linear=NULL,
         of=NULL,
         precision=NULL,
         si=NULL)
{
    debug = FALSE

    ## this is a nice trick
    if (!is.null(copy)) {
        if (!missing(model))
            warning(paste("Ignored argument model=`", model,
                          "' in f() due to copy=`", copy, "'", sep=""))
        if (!is.null(of))
            stop("Argument `of=NULL' is required when `copy=...' is used.")
        model = "copy"
        of = copy
        copy = NULL
    } 

    if (is.null(model))
        stop("No model is specified.")
    inla.is.model(model,stop.on.error=TRUE)
    
    inla.check.control(control.group)
    cont.group = inla.set.control.group.default()
    cont.group[(namc = names(control.group))] = control.group
    
    ## CHECK ARGUMENTS.
    ## This is a bit tricky. We want to check if there are arguments
    ## in `...' which is of type `name = value' which are invalid. If
    ## so, this lead to an obscoure error later on...
    ##
    ## we first collect all arguments of type `name = value'
    args.eq = c() 
    for (arg in unlist(strsplit(as.character(as.expression(match.call(expand.dots=TRUE))), ",")))
        if (length(grep("=", arg)) > 0)
            args.eq = c(args.eq, gsub(" ", "", unlist(strsplit(arg, "="))[1]))
    ##
    ## then we compare these with the legal ones in INLA::f(), and
    ## flag an error its not among the legal ones.  OOPS: Need to add
    ## some dummy arguments which are those inside the extraconstr and
    ## Cmatrix argument, and inla.group() as well.
    arguments = c(names(formals(INLA::f)), "A", "e", "i", "j", "values", "method", "Cij")
    arguments = arguments[-grep("^[.][.][.]$", arguments)]
    for(elm in args.eq) {
        if (!is.element(elm, arguments)) {
            f.call = as.character(as.expression(match.call(expand.dots=TRUE)))
            valid.args = inla.paste(sort(arguments), sep="\n\t")
            stop(paste("Argument `", elm, "' in formula specification\n\n\t\t",
                       f.call, "\n\n  is invalid. Valid arguments are:\n\n\t", valid.args, sep=""))
        }
    }

    if (!is.null(Qmatrix))
        stop("Argument Qmatrix in f(), has changed name to Cmatrix; please fix...")
    
    ## check that the Q matrix is defined if and only if hte model
    ## is generic. same with the Cmatrix
    if (inla.one.of(model, c("generic", "generic0","generic1", "generic2"))) {
        if (is.null(Cmatrix))
            stop("For generic models the Cmatrix has to be provided")

        if (!is.character(Cmatrix) && !is.list(Cmatrix))
            stop("Argument `Cmatrix' is not of type `character' or `list'")
            
        if (is.character(Cmatrix)) {
            if (!file.exists(Cmatrix))
                stop("Filename defined in argument `Cmatrix' does not exists.")
        }
        else {
            inla.sparse.check(Cmatrix)
            if (is.null(Cmatrix$i) || is.null(Cmatrix$j) || is.null(Cmatrix$values))
                stop("List defined in argument `Cmatrix' is not of type `Cmatrix = list(i=c(), j=c(), values=c())'")
            if (length(Cmatrix$i) != length(Cmatrix$j) || length(Cmatrix$i) != length(Cmatrix$values))
                stop("Entries in the list `Cmatrix' has not equal length")
        }

        if (is.null(n)) {
            ## find `n' from the Cmatrix

            if (is.character(Cmatrix)) {
                cm = read.table(Cmatrix, col.names=c("i", "j", "values"))
                n = max(cm$i, cm$j)
            } else {
                n = max(Cmatrix$i, Cmatrix$j)
            }
        }
    } else {
        if (!is.null(Cmatrix))
            stop("Cmatrix is only used for generic models")
    }

    ## chech that the graph.file is provided, if required. Set 'n' from the graph.file.
    if(inla.one.of(model, c("besag", "bym"))) {
        if (is.null(graph.file))
            stop(paste("The graph.file has to be provided for model", model))
        if (!file.exists(graph.file))
            stop(paste("Cannot find graph.file", graph.file))
        ## read n from the graph
        n.from.graph = scan(graph.file, n=1, what = integer(0), comment.char = "#") #
        if (n.from.graph <= 0)
            stop(paste("Argument 'n from graph.file' is void:", n.from.graph))
        if (!is.null(n) && n != n.from.graph)
            stop(paste("Argument 'n' and 'n from graph.file' does not match", n, n.from.graph))
        n = n.from.graph
    }
    if(inla.one.of(model, c("besag2"))) {
        if (is.null(graph.file))
            stop(paste("The graph.file has to be provided for model", model))
        if (!file.exists(graph.file))
            stop(paste("Cannot find graph.file", graph.file))
        ## read n from the graph
        n.from.graph = 2*scan(graph.file, n=1, what = integer(0), comment.char = "#") #
        if (n.from.graph <= 0)
            stop(paste("Argument 'n from graph.file' is void:", n.from.graph))
        if (!is.null(n) && n != n.from.graph)
            stop(paste("Argument 'n' and 2*'n from graph.file' does not match", n, n.from.graph))
        n = n.from.graph
    }

    ## is N required?
    if (is.null(n) && (!is.null(inla.model.properties(model)$n.required) && inla.model.properties(model)$n.required)) {
        stop(paste("Argument `n' in f() is required for model:", model))
    }
    
    ## special N required?
    if ((!is.null(inla.model.properties(model)$n.div.by) && inla.model.properties(model)$n.div.by) && !is.null(n)) {
        if (!inla.divisible(n, inla.model.properties(model)$n.div.by))
            stop(paste("Argument `n'", n, "is not divisible by", inla.model.properties(model)$n.div.by))
    }

    ## set default 'values'?
    if (!is.null(n) && is.null(values) &&
        (!is.null(inla.model.properties(model)$set.default.values) && inla.model.properties(model)$set.default.values)) {
        values = 1:n
    }

    if(inla.one.of(model, "linear")) {
        vars = as.list(substitute(list(...)))[-1]
        d = length(vars)
        term = deparse(vars[[1]], backtick = TRUE, width.cutoff = 500) 

        if(d!=1)
            stop("Too many or no term specified ")
        term = attr(terms(reformulate(term)),"term.labels")
        ret = list(d=d, term=term, model=model, mean.linear=mean.linear, prec.linear=prec.linear, label=term,
                cdf=cdf, quantiles = quantiles)
    }
    
    if (inla.one.of(model, "positive")) {
        if (!is.null(n))
            stop(paste("model = positive require n = 1, not n =",n))
        n=1
    }
            
    if(!is.null(prec.linear) | !is.null(mean.linear))
        stop("'mean.linear' and 'prec.linear' defined only for model='linear'")

    if (inla.one.of(model, "sphere")) {
        if (is.null(sphere.dir))
            stop("Argument sphere.dir=NULL is required for model = sphere")
        if (!(file.exists(sphere.dir) && file.info(sphere.dir)$isdir))
            stop(paste("Argument sphere.dir=", sphere.dir, "must be an existing directory."))

        ## set the Occilation parameter default to fixed, unless its set already
        if (missing(fixed)) {
            fixed = c(FALSE, FALSE, FALSE, TRUE)
        }
        if (missing(initial)) {
            initial = c(NA, NA, NA, -20)
        }
    } 

    if (inla.one.of(model, "spde")) {
        if (is.null(spde.prefix))
            stop("Argument spde.prefix=NULL is required for model = spde")
        ## file ``PREFIXs'' must exists... test this one
        if (!(file.exists(paste(spde.prefix, "s", sep=""))))
            stop(paste("Argument spde.prefix=", spde.prefix, "does not seems to be valid (no file `PREFIXs')"))

        ## set the Occilation parameter default to fixed, unless its set already
        if (missing(fixed)) {
            fixed = c(FALSE, FALSE, FALSE, TRUE)
        }
        if (missing(initial)) {
            initial = c(NA,NA,NA,-20)
        }
    } 
        
    ## in ... is the name of the covariate  and possibly the location of the weights
    ## like f(covariate,weights)
        
    vars = as.list(substitute(list(...)))[-1]
    d = length(vars)
    term = deparse(vars[[1]], backtick = TRUE, width.cutoff = 500)
        
    if (debug)
        print(vars)
        
    ##the second term in ... is the (possible) weights for the selected covariate!
    if(d==1) 
        weights = NULL
    else if(d==2)
        weights = deparse(vars[[2]], backtick = TRUE, width.cutoff = 500)         
    else if(d>2) 
        stop(paste("To many variables included in f():", inla.paste(vars)))
    else if (d==0) 
        stop("At least one variable in f() needs to be defined")
        
    if(inla.one.of(model,"seasonal") &&
       is.null(season.length)) stop("The length of the season has to
        be provided in season.length")
        
    ## cyclic is only valid for rw1, rw2 and rw2d-models
    if(!is.null(cyclic) && cyclic &&
       !inla.one.of(model, c("rw1", "rw2", "rw2d", "rw1c2", "rw2c2")))
        stop("Cyclic defined only for rw1, rw1c2, rw2, rw2c2 and rw2d models")
        

    need.nrow.ncol = inla.model.properties(model)$nrow.ncol
    ## nrow/ncol
    if ((!is.null(nrow) || !is.null(ncol)) && !need.nrow.ncol)
        stop(paste("nrow and ncol are not needed for model = ", model))
    if (need.nrow.ncol) {
        if (is.null(nrow) || is.null(ncol))
            stop(paste("nrow and ncol must be specified for model", model))
        if (nrow <= 0 || ncol <= 0 || trunc(nrow) != nrow || trunc(ncol) != ncol)
            stop("nrow and ncol must be positive intergers.")
        ## then it's ok...
        if (!is.null(values))
            stop(paste("values are not used for model = ", model))
    }

    ## check the NU parameter
    if (model != "matern2d" && model != "matern2dx2part0" && model != "matern2dx2p1") {
        if (!is.null(nu))
            stop("Argument NU is only used for matern2d/matern2dx2(part0/p1)-model.")
    }
    else {
        if (!is.null(nu) && !is.element(nu, c(0,1,2,3)))
            stop("For matern2d/matern2dx2(part0/p1)-model, the NU-parameter must be 0,1,2 or 3.")
    }
        
    ## get the weights
    term = attr(terms(reformulate(term)),"term.labels")
    if(d > 1)
        weigths = attr(terms(reformulate(weights)),"weights.labels")

    ##for all instrinsic model the constraint has to be ON...
    ##...except if the rw is cyclic!!!!!
    if(is.null(constr)) {
        constr = inla.model.properties(model)$constr
        if(!is.null(cyclic) && cyclic)
            constr=FALSE
    }   
        
    ## if diagonal is not set, the set this depending on the constr
    if (is.null(diagonal)) {
        if (constr)
            diagonal = inla.set.f.default()$diagonal
        else
            diagonal = 0
    }
        
    if (inla.one.of(model, "z") && is.null(Z))
        stop("With model [z] then covariate-matrix Z is required. Example: f(ind, Z=Z, model=\"z\")")

    if(!is.null(extraconstr)) {
        ##check
        A=extraconstr$A
        e=extraconstr$e
        if(!is.matrix(A))
            stop("A(extraconstraint) has to be a matrix")
        else
            if(nrow(A)!=length(e))
                stop("Dimension of A and e do not correspond")
        ##print.extraconstr = paste("list(A=", deparse(extraconstr$A, backtick = TRUE, width.cutoff = 500),",",
        ##    "e=",deparse(extraconstr$e, backtick = TRUE, width.cutoff = 500),")")
    }

    prop = inla.model.properties(model, stop.on.error=TRUE)
    if (!is.null(param))
        if (length(param) != prop$nparameters)
            stop(paste("The length of `param' in ", model, " has to be ", prop$nparameters))
        
    ## FIXED and no INITIAL
    if (is.null(initial) && !is.null(fixed) && !all(fixed == FALSE))
        stop("The fixed hyperparameters in model", model,"needs initialisation; use initial=c(....)")
        
    ## exceptions to default fixed = 0...
    if (model == "copy" && is.null(fixed)) {
        fixed = 1
        if (is.null(initial))
            initial = 1
    }

    if (!is.null(initial))
        if (length(initial) != prop$ntheta && !is.na(prop$ntheta))
            stop(paste("The length of `initial' on model", model,"has to be ", prop$ntheta))
        
    if (!is.null(fixed)) {
        if (length(fixed) != prop$ntheta && !is.na(prop$ntheta))
            stop(paste("The length of `fixed' on model", model,"has to be ", prop$ntheta))
    }
    else {
        if (prop$ntheta && !is.na(prop$ntheta))
            fixed = rep(0, prop$ntheta)
        else
            fixed = NULL
    }
    
    if (!is.null(fixed))
        for(j in 1:prop$ntheta)
            if (fixed[j] && !is.numeric(initial[j]))
                stop(paste("Model ", model, ", parameter no ", j, ", is fixed but has no intial value", sep=""))

    if (!is.null(prior)) {
        if (length(prior) != prop$npriors)
            stop(paste("Prior specification for model", model,"must contains", prop$npriors, "terms"))
        for(i in 1:length(prior))
            if (!is.null(prior[i]))
                inla.is.prior(prior[i], stop.on.error=TRUE)
    }

    ret=list(d=d, term=term, weights=weights, n=n, nrep = nrep, replicate = replicate,
            ngroup = ngroup, group = group, control.group = cont.group,
            Z=Z, model=model, prior=prior, 
            initial=initial, diagonal = diagonal, param = param,  fixed = fixed,
            cyclic=cyclic, season.length=season.length, 
            constr = constr, label=term, graph.file=graph.file, cdf=cdf, quantiles = quantiles,
            Cmatrix = Cmatrix, rankdef=rankdef, extraconstr=extraconstr, values=values,
            nrow = nrow, ncol = ncol, nu = nu, bvalue = bvalue,
            sphere.dir = sphere.dir, T.order = T.order, T.model = T.model, K.order = K.order, K.model = K.model,
            of = of, precision = precision, si = si,
            spde.prefix = spde.prefix )

    return (ret)
}
