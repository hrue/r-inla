##!\name{f}
##!\alias{f}
##!\title{ Defining smooth and spatial effects in INLA formulae }
##!\description{
##!
##!  Function used in definition of smooth and spatial terms within \code{inla} model
##!  formulae. The function does not evaluate a anything - it
##!  exists purely to help set up a model. The function specifies one
##!  smooth function in the linear predictor (see \code{\link{inla.models}}) as
##!  \deqn{w\ f(x)}{weight*f(var)}
##!
##!}
##!\usage{
##!`f` =
##!function(..., model=NULL,...)
##!}
##!\arguments{
`f` = function(
        ##!\item{...}{ Name of the covariate and, possibly of the
        ##!weights vector. NB: order counts!!!! The first specified
        ##!term is the covariate and the second one is the vector of
        ##!weights.}
        ...,

        ##!\item{model}{ A string indicating the choosen model. The
        ##! default is \code{iid}. See
        ##! \code{names(inla.models()$latent)} for a list of possible
        ##! alternatives.}
        model = "iid",

        ##!\item{copy}{}
        copy=NULL,

        ##!\item{same.as}{}
        same.as = NULL,

        ##!\item{n}{An optional argument which defines the dimension
        ##!of the model if this is different from
        ##!\code{length(sort(unique(covariate)))}}
        n=NULL,

        ##!\item{nrep}{}
        nrep = NULL,

        ##!\item{replicate}{}
        replicate = NULL,

        ##!\item{group}{}
        ngroup = NULL,

        ##!\item{group}{}
        group = NULL,

        ##!\item{control.group}{}
        control.group = inla.set.control.group.default(),

        ##!\item{hyper}{Spesification of the hyperparameter, fixed or
        ##!random, initial values, priors and its parameters. See
        ##!\code{?inla.models} for the list of hyparameters for each
        ##!model and its default options.}
        hyper = NULL,

        ##!\item{initial}{THIS OPTION IS OBSOLETE; use
        ##!\code{hyper}!!! Vector indicating the starting values for
        ##!the optimization algorithm. The length of the vector
        ##!depends on the number of hyperparamters in the choosen
        ##!\code{model}. If \code{fixed=T} the value at which the
        ##!parameters are fixed is determines through \code{initial}.
        ##!See \code{inla.models()$latent$'model name'} to have info
        ##!about the choosen model.}
        initial=NULL,

        ##!\item{prior}{THIS OPTION IS OBSOLETE; use \code{hyper}!!!
        ##!Prior distribution(s) for the hyperparameters of the
        ##!!random model. The default value depends on the type of
        ##!model, see !\url{www.r-inla.org} for a detailed
        ##!description of the models. See
        ##!\code{names(inla.models()$priors)} for possible prior
        ##!choices}
        prior=NULL,

        ##!\item{param}{THIS OPTION IS OBSOLETE; use \code{hyper}!!!
        ##!Vector indicating the parameters \eqn{a}{a} and \eqn{b}{b}
        ##!of the prior distribution for the hyperparameters. The
        ##!length of the vector depends on the choosen \code{model}.
        ##!See \code{inla.models()$latent$'model name'} to have info
        ##!about the choosen model.}
        param = NULL,

        ##!\item{fixed}{THIS OPTION IS OBSOLETE; use \code{hyper}!!!
        ##!Vector of boolean variables indicating wheater the
        ##!hyperparameters of the model are fixed or random. The
        ##!length of the vector depends on the choosen \code{model}
        ##!See \code{inla.models()$latent$'model name'} to have info
        ##!about the choosen model.}
        fixed = NULL,

        ##!\item{season.length}{Lenght of the seasonal compoment
        ##!(ONLY if \code{model="seasonal"}) }
        season.length=NULL,

        ##!\item{constr}{A boolean variable indicating whater to set
        ##!a sum to 0 constraint on the term. By default the sum to 0
        ##!constraint is imposed on all intrinsic models
        ##!("iid","rw1","rw1","besag", etc..).}
        constr = NULL,

        ##!\item{extraconstr}{This argument defines extra linear
        ##!constraints. The argument is a list with two elements, a
        ##!matrix \code{A} and a vector \code{e}, which defines the
        ##!extra constraint \code{Ax = e}; for example
        ##!\code{extraconstr = list(A = A, e=e)}. The number of
        ##!columns of \code{A} must correspond to the length of this
        ##!\code{f}-model.  Note that this constraint comes
        ##!additional to the sum-to-zero constraint defined if
        ##!\code{constr = TRUE}.}
        extraconstr=NULL,

        ##!\item{values}{An optional vector giving all values
        ##!assumed by the covariate for which we want estimated the
        ##!effect. It must be a numeric vector, a vector of factors
        ##!or \code{NULL}.}
        values=NULL,

        ##!\item{cyclic}{A boolean specifying wheather the model is
        ##!cyclical. Only valid for "rw1" and "rw2" models, is
        ##!cyclic=T then the sum to 0 constraint is removed. For the
        ##!correct form of the grah file see \cite{Martinoand Rue
        ##!(2008)}.}
        cyclic = NULL,

        ##!\item{diagonal}{An extra constant added to the diagonal of
        ##!the precision matrix.}
        diagonal = NULL,

        ##!\item{graph}{Defines the graph-object either as a file with
        ##!a graph-description, an \code{inla.graph}-object, or as a
        ##!(sparse) symmetric matrix.}
        graph=NULL,

        ##!\item{graph.file}{THIS OPTION IS OBSOLETE AND REPLACED BY
        ##!THE MORE GENERAL ARGUMENT \code{graph}. PLEASE CHANGE YOUR
        ##!CODE.
        ##!Name of the file containing the graph
        ##!of the model; see
        ##!\url{http://www.r-inla.org/help/faq}.}
        graph.file=NULL,

        ##!\item{cdf}{A vector of maximum 10 values between 0 and 1
        ##!\eqn{x(0), x(1),\ldots}{x(0), x(1),\ldots}. The function
        ##!returns, for each posterior marginal the probabilities
        ##!\deqn{\mbox{Prob}(X<x(p))}{Prob(X<x(p))} }
        cdf=NULL,

        ##!\item{quantiles}{A vector of maximum 10 quantiles,
        ##!\eqn{p(0), p(1),\dots}{p(0), p(1),\ldots} to compute for
        ##!each posterior marginal. The function returns, for each
        ##!posterior marginal, the values
        ##!\eqn{x(0), x(1),\dots}{x(0), x(1),\ldots} such that
        ##!\deqn{\mbox{Prob}(X<x(p))=p}{Prob(X<x)=p}}
        quantiles=NULL,

        ##!\item{Cmatrix}{The specification of the precision matrix
        ##!for the generic models (up to a scaling constant).
        ##!\code{Cmatrix} is either a
        ##!(dense) matrix, a matrix create using
        ##!\code{Matrix::sparseMatrix()}, or a filename which stores the
        ##!non-zero elements of \code{Cmatrix}, in three columns:
        ##!\code{i}, \code{j} and \code{Qij}.}
        Cmatrix=NULL,

        ##!\item{rankdef}{A number \bold{defining} the rank
        ##!deficiency of the model, with sum-to-zero constraint and
        ##!possible extra-constraints taken into account. See
        ##!details.}
        rankdef=NULL,

        ##!\item{Z}{}
        Z = NULL,

        ##!\item{nrow}{Number of rows for 2d-models}
        nrow = NULL,

        ##!\item{ncol}{Number of columns for 2d-models}
        ncol = NULL,

        ##!\item{nu}{Smoothing parameter for the Matern2d-model,
        ##!possible values are \code{c(0, 1, 2, 3)}}
        nu = NULL,

        ##!\item{bvalue}{}
        bvalue = NULL,

        ##!\item{spde.prefix}{}
        spde.prefix = NULL,

        ##!\item{spde2.prefix}{}
        spde2.prefix = NULL,

        ##!\item{spde2.transform}{}
        spde2.transform = c("logit", "log", "identity"),

        ##!\item{mean.linear}{Prior mean for the linear component,
        ##!only used if \code{model="linear"}}
        mean.linear = inla.set.control.fixed.default()$mean, 

        ##!\item{prec.linear}{Prior precision for the linear
        ##!component, only used if \code{model="linear"}}
        prec.linear = inla.set.control.fixed.default()$prec, 

        ##!\item{si}{}
        si=FALSE,

        ##!\item{compute}{ A boolean variable indicating wheather the
        ##! marginal posterior distribution for the nodes in the
        ##! \code{f()} model should be computed or not. This is
        ##! usefull for large models where we are only interested in
        ##! some posterior marginals.}
        compute = TRUE,

        ##!\item{of}{}
        of=NULL,

        ##!\item{precision}{The precision for the artifical noise added when creating a copy of a model.}
        precision=NULL,

        ##!\item{range}{A vector of size two giving the lower and
        ##!upper range for the scaling parameter \code{beta} in the
        ##!model \code{COPY}. If low = high then the identity mapping
        ##!is used.}
        range = NULL,

        ##!\item{adjust.for.con.comp}{If TRUE (default), adjust some
        ##!of the models (currently: besag, bym and besag2) if the
        ##!number of connected components in graph is larger than
        ##!1. If FALSE, do nothing.}
        adjust.for.con.comp = TRUE,

        ##!\item{scale}{A scaling vector. Its meaning depends on the model.}
        scale = NULL, 

        ##!\item{strata}{A stratum vector. It meaning depends on the model.}
        strata = NULL, 

        ##!\item{rgeneric}{A object of class \code{inla-rgeneric} which defines the model. (EXPERIMENTAL!)}
        rgeneric = NULL, 

        ## local debug-flag
        debug = FALSE)
{
    ##!}
    ##!\value{}

    ##!\details{There is no default value for \code{rankdef}, if it
    ##!is not defined by the user then it is computed by the rank
    ##!deficiency of the prior model (for the generic model, the
    ##!default is zero), plus 1 for the sum-to-zero constraint if the
    ##!prior model is proper, plus the number of extra
    ##!constraints. \bold{Oops:} This can be wrong, and then the user
    ##!must define the \code{rankdef} explicitely.}
    ##!\author{Havard Rue \email{hrue@math.ntnu.no}}
    ##!\seealso{\code{\link{inla}}, \code{\link{hyperpar.inla}}}

    ## this is required. the hyper.defaults can only be changed in the
    ## model=model.object
    hyper.default = NULL

    ## if model is a particular class, then use this to set default
    ## arguments to all names(formals(f)) (except the "...").
    ## THIS FEATURE IS EXPERIMENTAL FOR THE MOMENT!  OOPS: the classes
    ## are defined in the function inla.model.object.classes() as this
    ## is also used in the inla() function itself.
    if (any(inherits(model, inla.model.object.classes()))) {
        arguments = names(formals(INLA::f))
        arguments = arguments[ -grep("^[.][.][.]$",arguments) ]
        ## add this one manually
        arguments = c(arguments, "hyper.default")

        ## evaluate arguments in 'model' and set those that are in
        ## 'arguments'. However, due to the 'missing || is.null', the
        ## argument 'model' does not pass the test and must be set
        ## manually afterwards.
        for(nm in names(model$f)) {
            if (inla.one.of(nm, arguments)) {
                inla.eval(paste("if (missing(",nm,") || is.null(", nm,")) ",
                                nm, "=", "model$f$", nm,  sep=""))
            } else {
                ## this should not happen
                stop(paste("Argument `", nm, "' is not an argument in f().",
                           " This is likely not what you want!", sep=""))
            }
        }

        ## the 'model' argument is required!
        stopifnot(is.character(model$f$model))
        model = model$f$model
    }

    ## this is a nice trick
    if (!is.null(copy)) {
        if (!missing(model)) {
            warning(paste("Ignored argument model=`", model,
                          "' in f() due to copy=`", copy, "'", sep=""))
        }
        if (!is.null(of)) {
            stop("Argument `of=NULL' is required when `copy=...' is used.")
        }
        model = "copy"
        of = copy
        copy = NULL
    }

    if (is.null(model)) {
        stop("No model is specified.")
    }
    inla.is.model(model, "latent", stop.on.error=TRUE)

    ## in ... is the name of the covariate  and possibly the location of the weights
    ## like f(covariate, weights)
    vars = as.list(substitute(list(...)))[-1]
    d = length(vars)
    if (d == 0L)
        stop(paste("Missing covariate in f(...) for model=", model))
    term = deparse(vars[[1]], backtick = TRUE, width.cutoff = 500)
    if (debug) {
        print(vars)
    }
    ## the second term in ... is the (possible) weights for the selected covariate!
    if (d==1) {
        weights = NULL
    } else if (d==2) {
        weights = deparse(vars[[2]], backtick = TRUE, width.cutoff = 500)
    } else if (d>2) {
        stop(paste("To many variables included in f():", inla.paste(vars)))
    } else if (d==0) {
        stop("At least one variable in f() needs to be defined")
    }
    ## get the weights
    term = attr(terms(reformulate(term)),"term.labels")
    if (d > 1) {
        weigths = attr(terms(reformulate(weights)),"weights.labels")
    }

    ## set the hyperparameters
    hyper = inla.set.hyper(model = model,  section = "latent",
            hyper = hyper, hyper.default = hyper.default,
            initial = initial, fixed = fixed,  prior = prior,  param = param)

    ## for model = copy, its is not allowed to define constr or extraconstr
    if (inla.one.of(model, "copy")) {
        stopifnot(missing(constr))
        stopifnot(missing(extraconstr))
    } else {
        ## this is only used for model = copy.
        stopifnot(missing(range))
    }

    if (!missing(same.as) && !is.null(same.as) && model != "copy") {
        stop(paste("Argument 'same.as' must be NULL is model != 'copy' :", model))
    }

    inla.check.control(control.group)
    cont.group = inla.set.control.group.default()
    cont.group[(namc = names(control.group))] = control.group
    cont.group$hyper = inla.set.hyper(cont.group$model, "group", cont.group$hyper,
            cont.group$initial, cont.group$fixed, cont.group$prior, cont.group$param)

    ## CHECK ARGUMENTS.
    ## This is a bit tricky. We want to check if there are arguments
    ## in `...' which is of type `name = value' which are invalid. If
    ## so, this lead to an obscoure error later on...
    ##
    ## we first collect all arguments of type `name = value'
    if (TRUE) {
        ## New code
        args.eq = names(match.call(expand.dots = TRUE))
        args.eq = args.eq[ args.eq != "" ]
    } else {
        args.eq = c()
        for (arg in unlist(strsplit(as.character(as.expression(match.call(expand.dots=TRUE))), ","))) {
            if (length(grep("=", arg)) > 0) {
                args.eq = c(args.eq, gsub(" ", "", unlist(strsplit(arg, "="))[1]))
            }
        }
    }

    ## then we compare these with the legal ones in f(), and
    ## flag an error its not among the legal ones.  OOPS: Need to add
    ## some dummy arguments which are those inside the extraconstr and
    ## Cmatrix argument, and inla.group() as well.
    arguments = c(names(formals(INLA::f)), "A", "e")
    arguments = arguments[-grep("^[.][.][.]$", arguments)]
    for(elm in args.eq) {
        if (!is.element(elm, arguments)) {
            f.call = as.character(as.expression(match.call(expand.dots=TRUE)))
            valid.args = inla.paste(sort(arguments), sep="\n\t")
            stop(paste("Argument `", elm, "' in formula specification\n\n\t\t",
                       f.call, "\n\n  is invalid. Valid arguments are:\n\n\t", valid.args, sep=""))
        }
    }

    ## check that the Q matrix is defined if and only if the model is
    ## generic. same with the Cmatrix
    if (inla.one.of(model, c("generic", "generic0","generic1", "generic2"))) {
        if (is.null(Cmatrix)) {
            stop("For generic models the Cmatrix has to be provided")
        }
        Cmatrix = inla.sparse.check(Cmatrix)
        if (is.null(n)) {
            n = inla.sparse.dim(Cmatrix)[1]
        }
    } else {
        if (!is.null(Cmatrix)) {
            stop("Cmatrix is only used for generic models")
        }
    }

    if (!is.null(graph) && !is.null(graph.file)) {
        stop("Don't know what to do, as 'graph' and 'graph.file' is set; please use only argument 'graph'.")
    }
    if (is.null(graph) && !is.null(graph.file)) {
        ## then it is easy and we can just go ahead with 'graph'
        graph = graph.file

        ## show a warning,  if this hasn't been shown before.
        if (inla.getOption("show.warning.graph.file")) {
            warning("Argument 'graph.file' in 'f()' is obsolete; please use the more general argument 'graph' instead.")
            ## disable this warning from now on.
            inla.setOption(show.warning.graph.file=FALSE)
        }
    }

    ## chech that the graph is provided, if required. Set 'n' from the graph.
    if (inla.one.of(model, c("besag", "bym", "besagproper"))) {
        if (is.null(graph)) {
            stop(paste("The 'graph' has to be provided for model", model))
        }
        n.from.graph = inla.read.graph(graph)$n
        if (n.from.graph <= 0) {
            stop(paste("Argument 'n from graph' is void:", n.from.graph))
        }
        if (!is.null(n) && n != n.from.graph) {
            stop(paste("Argument 'n' and 'n from graph' does not match", n, n.from.graph))
        }
        n = n.from.graph
    }
    if (inla.one.of(model, c("besag2"))) {
        if (is.null(graph)) {
            stop(paste("The graph has to be provided for model", model))
        }
        ## read n from the graph
        n.from.graph = 2L*inla.read.graph(graph)$n
        if (n.from.graph <= 0) {
            stop(paste("Argument 'n from graph' is void:", n.from.graph))
        }
        if (!is.null(n) && n != n.from.graph) {
            stop(paste("Argument 'n' and 2*'n from graph' does not match", n, n.from.graph))
        }
        n = n.from.graph

        if (!is.null(constr) && constr) {
            stop(paste("'constr=TRUE' does not make sense for model 'besag2'"))
        }
    }

    ## is N required?
    if (is.null(n) && (!is.null(inla.model.properties(model, "latent")$n.required)
                       && inla.model.properties(model, "latent")$n.required)) {
        stop(paste("Argument `n' in f() is required for model:", model))
    }

    ## special N required?
    if ((!is.null(inla.model.properties(model, "latent")$n.div.by)
         && inla.model.properties(model, "latent")$n.div.by) && !is.null(n)) {
        if (!inla.divisible(n, inla.model.properties(model, "latent")$n.div.by)) {
            stop(paste("Argument `n'", n, "is not divisible by", inla.model.properties(model, "latent")$n.div.by))
        }
    }

    ## set default 'values'?  Do the check for values for nrow.ncol
    ## models further below.
    if (!is.null(n) && is.null(values) &&
        (!is.null(inla.model.properties(model, "latent")$set.default.values)
         && inla.model.properties(model, "latent")$set.default.values)) {
        values = 1:n
    }

    if (inla.one.of(model, "linear")) {
        if (d != 1L) {
            stop("Model = 'linear' do not accept weights. Just set 'z.new = z * weights' as the covariates.")
        }
        ret = list(d=d, term=term, model=model, mean.linear=mean.linear, prec.linear=prec.linear, label=term,
                cdf=cdf, quantiles = quantiles, compute = compute)
        ## return here!
        return (ret)
    }
    if (!missing(prec.linear) || !missing(mean.linear)) {
        stop("Arguments 'mean.linear' and 'prec.linear' are defined only for model='linear'.")
    }

    if (inla.one.of(model, "positive")) {
        if (!is.null(n)) {
            stop(paste("model = positive require n = 1, not n =", n))
        }
        n=1
    }

    if (inla.one.of(model, c("spde"))) {
        if (is.null(spde.prefix)) {
            stop("Argument spde.prefix=NULL is required for model = spde")
        }
        ## file ``PREFIXs'' must exists... test this one
        if (!(file.exists(paste(spde.prefix, "s", sep="")))) {
            stop(paste("Argument spde.prefix=", spde.prefix, "does not seems to be valid (no file `PREFIXs')"))
        }
    }

    if (inla.one.of(model, c("spde2"))) {
        if (is.null(spde2.prefix)) {
            stop("Argument spde2.prefix=NULL is required for model = spde2")
        }
    }
    spde2.transform = match.arg(spde2.transform, several.ok = FALSE)

    if (inla.one.of(model,"seasonal") &&
        is.null(season.length)) {
        stop("The length of the season has to be provided in season.length")
    }

    ## cyclic is only valid for rw1, rw2 and rw2d-models
    if (!is.null(cyclic) && cyclic &&
        !inla.one.of(model, c("rw1", "rw2", "rw2d", "rw1c2", "rw2c2"))) {
        stop("Cyclic defined only for rw1, rw1c2, rw2, rw2c2 and rw2d models")
    }

    need.nrow.ncol = inla.model.properties(model, "latent")$nrow.ncol
    ## nrow/ncol
    if ((!is.null(nrow) || !is.null(ncol)) && !need.nrow.ncol) {
        stop(paste("nrow and ncol are not needed for model = ", model))
    }
    if (need.nrow.ncol) {
        if (is.null(nrow) || is.null(ncol)) {
            stop(paste("nrow and ncol must be specified for model", model))
        }
        if (nrow <= 0 || ncol <= 0 || trunc(nrow) != nrow || trunc(ncol) != ncol) {
            stop("nrow and ncol must be positive intergers.")
        }

        ## set n as well, makes it easier.
        if (!missing(n)) {
            stopifnot(n == nrow*ncol)
        } else {
            n = nrow * ncol
        }

        ## and set default values, if required
        if (inla.model.properties(model, "latent")$set.default.values) {
            if (missing(values)) {
                values = 1:n
            } else {
                stopifnot(length(values) == n)
            }
        }
    }

    ## check the NU parameter
    if (model != "matern2d" && model != "matern2dx2part0" && model != "matern2dx2p1") {
        if (!is.null(nu)) {
            stop("Argument NU is only used for matern2d/matern2dx2(part0/p1)-model.")
        }
    } else {
        if (!is.null(nu) && !is.element(nu, c(0, 1, 2, 3))) {
            stop("For matern2d/matern2dx2(part0/p1)-model, the NU-parameter must be 0, 1, 2 or 3.")
        }
    }

    ##for all instrinsic model the constraint has to be ON...
    ##...except if the rw is cyclic!!!!!
    if (is.null(constr)) {
        constr = inla.model.properties(model, "latent")$constr
        if (!is.null(cyclic) && cyclic) {
            constr=FALSE
        }
    }

    ## if diagonal is not set, the set this depending on the constr
    if (is.null(diagonal)) {
        if (constr) {
            diagonal = inla.set.f.default()$diagonal
        } else {
            diagonal = 0
        }
    }

    if (inla.one.of(model, "z") && is.null(Z)) {
        stop("With model [z] then covariate-matrix Z is required. Example: f(ind, Z=Z, model=\"z\")")
    }

    if (!is.null(extraconstr)) {
        ##check
        A=extraconstr$A
        e=extraconstr$e
        if (!is.matrix(A)) {
            stop("A(extraconstraint) has to be a matrix")
        } else {
            if (nrow(A)!=length(e)) {
                stop("Dimension of A and e do not correspond")
            }
        }

        ##print.extraconstr = paste("list(A=", deparse(extraconstr$A, backtick = TRUE, width.cutoff = 500),",",
        ##    "e=", deparse(extraconstr$e, backtick = TRUE, width.cutoff = 500),")")
    }

    if (adjust.for.con.comp && inla.one.of(model, c("besag", "besag2", "bym"))) {

        ## I am not sure if this is the best place to do this, but for
        ## the time being I'll do this here.

        ## For the model = "besag", "bym" and "besag2", we need to
        ## check if the graph contains more than 1 connected
        ## components. If so, we need to modify the meaning of
        ## constr=TRUE, and set the correct value of rankdef.

        ## However, although this is the correct way to do this, it is
        ## not common to do it like this, so therefore, I issue a
        ## warning.

        g = inla.read.graph(graph)
        if (g$cc$n == 1) {
            ## hole graph is just one connected component. all is
            ## fine, no need to do anything
        } else {

            ## issue a warning, as the model has 'changed' compared to earlier versions.
            warning(paste("The graph for the model", model, "has", g$cc$n,
                          "connected components!!! Model is revised accordingly.", sep= " "))

            cc.n = sapply(g$cc$nodes, length)
            cc.n2 = sum(cc.n >= 2L)
            dimA = 0

            if (debug) {
                print(paste("modify model", model))
                print(paste("number of connected components", inla.paste(cc.n)))
                print(paste("number of connected components of size >= 2L", cc.n2))
            }

            if (constr) {
                ## need to redefine the meaning of constr = TRUE to mean
                ## constr=TRUE for all connected components with size > 1

                ## like bym place the constr on the second half
                m = inla.model.properties(model, "latent")
                if (m$augmented) {
                    N = m$aug.factor * n
                    offset = (m$aug.constr -1L) * n
                    stopifnot(length(offset) == 1)
                } else {
                    N = n
                    offset = 0
                }
                constr = FALSE
                AA = matrix(0, cc.n2, N)
                ee = rep(0, cc.n2)

                if (debug) {
                    print(paste("add new extraconstr, dim = ", cc.n2, "x", n))
                }

                k = 1L
                for(i in 1L:length(cc.n)) {
                    if (cc.n[i] >= 2L) {
                        AA[k, offset + g$cc$nodes[[i]]] = 1
                        k = k + 1L
                    }
                }
                stopifnot(k-1L == cc.n2)

                if (!is.null(extraconstr)) {
                    dimA = dim(extraconstr$A)[1] ## need to remember this one for the 'rankdef'
                    extraconstr$A = rbind(AA, extraconstr$A)
                    extraconstr$e = c(ee, extraconstr$e)
                } else {
                    extraconstr = list(A = AA, e = ee)
                }
            }

            ## set correct rankdef if not set manually. well, this is
            ## not failsafe without detailed analysis, but then we
            ## asssume that the user knows the correct rankdef (s)he
            ## will set it correctly. Note that regions with no
            ## neigbours only reduce 'n' as they do not contribute to the rankdef.
            if (is.null(rankdef)) {
                cc.n1 = sum(cc.n == 1)
                rankdef = cc.n1 + dimA + (g$n - sum(cc.n[ cc.n >= 2L ] -1L))
                if (debug) {
                    print(paste("redefine rankdef = ", rankdef))
                }
            }

            ## need this, either if there is extraconstr or if there
            ## are regions without neighbours.
            if (is.null(diagonal) && (!is.null(extraconstr) || any(cc.n == 1))) {
                diagonal = inla.set.f.default()$diagonal
                if (debug) {
                    print(paste("set diagonal = ", diagonal))
                }
            }
        }
    }

    if (model %in% "rgeneric") {
        stopifnot(inherits(rgeneric, "inla-rgeneric"))
        if (inla.os("windows")) {
            stop("Model 'rgeneric' is not available for Windows; please use Linux or MacOSX. (No, there is no quick fix.)")
        }
    }

    ret=list(
            Cmatrix = Cmatrix,
            Z=Z,
            bvalue = bvalue,
            cdf=cdf,
            compute = compute,
            constr = constr,
            control.group = cont.group,
            cyclic=cyclic,
            d=d,
            diagonal = diagonal,
            extraconstr=extraconstr,
            graph=graph,
            group = group,
            hyper = hyper,
            label = term,
            model=model,
            n=n,
            ncol = ncol,
            ngroup = ngroup,
            nrep = nrep,
            nrow = nrow,
            nu = nu,
            of = of,
            precision = precision,
            quantiles = quantiles,
            range = range,
            rankdef=rankdef,
            replicate = replicate,
            same.as = same.as,
            season.length = season.length,
            si = si,
            spde.prefix = spde.prefix,
            spde2.prefix = spde2.prefix,
            spde2.transform = spde2.transform,
            term=term,
            values=values,
            weights=weights,
            scale = scale,
            strata = strata,
            rgeneric = rgeneric
            )

    return (ret)
}

## "inla.model.class" is a generic model class that can be inherited
## from to mark a class as an inla model object class.
`inla.model.object.classes` = function()
{
    return (c("inla.model.class", "inla.wrapper.model",
              "inla.spde", "inla.spde1", "inla.spde2"))
}
