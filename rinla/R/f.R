## Export: f

##!\name{f}
##!\alias{f}
##!\title{Define general Gaussian models in the INLA formula }
##!\description{
##!
##!  Function used for defining of smooth and spatial terms within \code{inla} model
##!  formulae. The function does not evaluate anything - it
##!  exists purely to help set up a model. The function specifies one
##!  smooth function in the linear predictor (see \code{\link{inla.list.models}}) as
##!  \deqn{w\ f(x)}{weight*f(var)}
##!
##!}
##!\usage{
##!    f(...,
##!         model = "iid",
##!         copy=NULL,
##!         same.as = NULL,
##!         n=NULL,
##!         nrep = NULL,
##!         replicate = NULL,
##!         ngroup = NULL,
##!         group = NULL,
##!         control.group = inla.set.control.group.default(),
##!         hyper = NULL,
##!         initial=NULL,
##!         prior=NULL,
##!         param = NULL,
##!         fixed = NULL,
##!         season.length=NULL,
##!         constr = NULL,
##!         extraconstr=list(A=NULL, e=NULL),
##!         values=NULL,
##!         cyclic = NULL,
##!         diagonal = NULL,
##!         graph=NULL,
##!         graph.file=NULL,
##!         cdf=NULL,
##!         quantiles=NULL,
##!         Cmatrix=NULL,
##!         rankdef=NULL,
##!         Z = NULL,
##!         nrow = NULL,
##!         ncol = NULL,
##!         nu = NULL,
##!         bvalue = NULL,
##!         spde.prefix = NULL,
##!         spde2.prefix = NULL,
##!         spde2.transform = c("logit", "log", "identity"),
##!         spde3.prefix = NULL,
##!         spde3.transform = c("logit", "log", "identity"),
##!         mean.linear = inla.set.control.fixed.default()$mean,
##!         prec.linear = inla.set.control.fixed.default()$prec,
##!         compute = TRUE,
##!         of=NULL,
##!         precision = exp(14),
##!         range = NULL,
##!         adjust.for.con.comp = TRUE,
##!         order = NULL,
##!         scale = NULL,
##!         strata = NULL,
##!         rgeneric = NULL,
##!         scale.model = NULL,
##!         args.slm = list(rho.min = NULL, rho.max = NULL, 
##!                         X = NULL, W = NULL, Q.beta = NULL),
##!         args.ar1c = list(Z = NULL, Q.beta = NULL),
##!         correct = NULL,
##!         locations = NULL, 
##!         debug = FALSE)
##!}
##!\arguments{
`f` = function(
    ##!\item{...}{ Name of the covariate and, possibly of the
    ##!weights vector. NB: order counts!!!! The first specified
    ##!term is the covariate and the second one is the vector of
    ##!weights (which can be negative).}
    ...,

    ##!\item{model}{ A string indicating the choosen model. The
    ##! default is \code{iid}. See
    ##! \code{names(inla.models()$latent)} for a list of possible
    ##! alternatives and \code{\link{inla.doc}} for detailed docs.}
    model = "iid",

    ##!\item{copy}{TODO}
    copy=NULL,

    ##!\item{same.as}{TODO}
    same.as = NULL,

    ##!\item{n}{An optional argument which defines the dimension
    ##!of the model if this is different from
    ##!\code{length(sort(unique(covariate)))}}
    n=NULL,

    ##!\item{nrep}{TODO}

    nrep = NULL,

    ##!\item{replicate}{We need to write documentation here}
    replicate = NULL,

    ##!\item{ngroup}{TODO}
    ngroup = NULL,

    ##!\item{group}{TODO}
    group = NULL,

    ##!\item{control.group}{TODO}
    control.group = inla.set.control.group.default(),

    ##!\item{hyper}{Specification of the hyperparameter, fixed or
    ##!random, initial values, priors and its parameters. See
    ##!\code{?inla.models} for the list of hyparameters for each
    ##!model and its default options or
	##!use \code{inla.doc()} for
	##!detailed info on the family and
	##!supported prior distributions.}
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
    extraconstr=list(A=NULL, e=NULL),

    ##!\item{values}{An optional vector giving all values
    ##!assumed by the covariate for which we want estimated the
    ##!effect. It must be a numeric vector, a vector of factors
    ##!or \code{NULL}.}
    values=NULL,

    ##!\item{cyclic}{A boolean specifying wheather the model is
    ##!cyclical. Only valid for "rw1" and "rw2" models, is
    ##!cyclic=T then the sum to 0 constraint is removed. For the
    ##!correct form of the grah file see \cite{Martino and Rue
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
    ##!\url{www.r-inla.org/faq}.}
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
    ##!for the generic,  generic3 or z models (up to a scaling constant).
    ##!\code{Cmatrix} is either a
    ##!(dense) matrix, a matrix created using
    ##!\code{Matrix::sparseMatrix()}, or a filename which stores the
    ##!non-zero elements of \code{Cmatrix}, in three columns:
    ##!\code{i}, \code{j} and \code{Qij}. In case of the generic3 model,
    ##!it is a list of such specifications.}
    Cmatrix=NULL,

    ##!\item{rankdef}{A number \bold{defining} the rank
    ##!deficiency of the model, with sum-to-zero constraint and
    ##!possible extra-constraints taken into account. See
    ##!details.}
    rankdef=NULL,

    ##!\item{Z}{The matrix for the z-model}
    Z = NULL,

    ##!\item{nrow}{Number of rows for 2d-models}
    nrow = NULL,

    ##!\item{ncol}{Number of columns for 2d-models}
    ncol = NULL,

    ##!\item{nu}{Smoothing parameter for the Matern2d-model,
    ##!possible values are \code{c(0, 1, 2, 3)}}
    nu = NULL,

    ##!\item{bvalue}{TODO}
    bvalue = NULL,

    ##!\item{spde.prefix}{TODO}
    spde.prefix = NULL,

    ##!\item{spde2.prefix}{TODO}
    spde2.prefix = NULL,

    ##!\item{spde2.transform}{TODO}
    spde2.transform = c("logit", "log", "identity"),

    ##!\item{spde3.prefix}{TODO}
    spde3.prefix = NULL,

    ##!\item{spde3.transform}{TODO}
    spde3.transform = c("logit", "log", "identity"),

    ##!\item{mean.linear}{Prior mean for the linear component,
    ##!only used if \code{model="linear"}}
    mean.linear = inla.set.control.fixed.default()$mean,

    ##!\item{prec.linear}{Prior precision for the linear
    ##!component, only used if \code{model="linear"}}
    prec.linear = inla.set.control.fixed.default()$prec,

    ##!\item{compute}{ A boolean variable indicating wheather the
    ##! marginal posterior distribution for the nodes in the
    ##! \code{f()} model should be computed or not. This is
    ##! usefull for large models where we are only interested in
    ##! some posterior marginals.}
    compute = TRUE,

    ##!\item{of}{TODO}
    of=NULL,

    ##!\item{precision}{The precision for the artifical noise added when creating a copy of a model and others.}
    precision = exp(14),

    ##!\item{range}{A vector of size two giving the lower and
    ##!upper range for the scaling parameter \code{beta} in the
    ##!model \code{COPY}, \code{CLINEAR}, \code{MEC} and \code{MEB}.
    ##!If \code{low = high} then the identity mapping
    ##!is used.}
    range = NULL,

    ##!\item{adjust.for.con.comp}{If TRUE (default), adjust some
    ##!of the models (currently: besag, bym, bym2 and besag2) if the
    ##!number of connected components in graph is larger than
    ##!1. If FALSE, do nothing.}
    adjust.for.con.comp = TRUE,

    ##!\item{order}{Defines the \code{order} of the model: for
    ##!model \code{ar} this defines the order p, in AR(p). Not
    ##!used for other models at the time being.}
    order = NULL,

    ##!\item{scale}{A scaling vector. Its meaning depends on the model.}
    scale = NULL,

    ##!\item{strata}{A stratum vector. It meaning depends on the model.}
    strata = NULL,

    ##!\item{rgeneric}{A object of class \code{inla.rgeneric} which defines the model. (EXPERIMENTAL!)}
    rgeneric = NULL,

    ##!\item{scale.model}{Logical. If \code{TRUE} then scale the RW1 and RW2 and BESAG and BYM and BESAG2 and RW2D models so the their (generlized) variance is 1. Default value is \code{inla.getOption("scale.model.default")}}
    scale.model = NULL,

    ##!\item{args.slm}{Required arguments to the model="slm"; see the documentation for further details.},
    args.slm = list(rho.min = NULL, rho.max = NULL, X = NULL, W = NULL, Q.beta = NULL),

    ##!\item{args.ar1c}{Required arguments to the model="ar1c"; see the documentation for further details.},
    args.ar1c = list(Z = NULL, Q.beta = NULL),

    ##!\item{correct}{Add this model component to the list of variables to be used in the corrected Laplace approximation? If \code{NULL} use default choice,  otherwise correct if \code{TRUE} and do not if \code{FALSE}. (This option is currently experimental.)},
    correct = NULL,

    ##!\item{locations}{A matrix with locations for the model \code{dmatern}. This also defines \code{n}.}
    locations = NULL,
    
    ##!\item{debug}{Enable local debug output}
    debug = FALSE)
{
    ##!}
    ##!\value{TODO}

    ##!\details{There is no default value for \code{rankdef}, if it
    ##!is not defined by the user then it is computed by the rank
    ##!deficiency of the prior model (for the generic model, the
    ##!default is zero), plus 1 for the sum-to-zero constraint if the
    ##!prior model is proper, plus the number of extra
    ##!constraints. \bold{Oops:} This can be wrong, and then the user
    ##!must define the \code{rankdef} explicitely.}
    ##!\author{Havard Rue \email{hrue@r-inla.org}}
    ##!\seealso{\code{\link{inla}}, \code{\link{hyperpar.inla}}}

    ## this is required. the hyper.defaults can only be changed in the
    ## model=model.object
    hyper.default = NULL

    empty.extraconstr = function(ec)
    {
        return (is.null(ec) || all(sapply(ec, is.null)))
    }

    ## if model is a particular class, then use this to set default
    ## arguments to all names(formals(f)) (except the "...").
    ## THIS FEATURE IS EXPERIMENTAL FOR THE MOMENT!  OOPS: the classes
    ## are defined in the function inla.model.object.classes() as this
    ## is also used in the inla() function itself.
    if (any(inherits(model, inla.model.object.classes()))) {
        if (any(inherits(model, inla.spde.object.classes()))) {
            ## Write specific spde models to disk. These files are
            ## removed later in inla.copy.dir.for.section.spde, in section.R
            my.spde.prefix <- inla.fmesher.make.prefix(NULL, NULL)
            if (any(inherits(model, "inla.spde1"))) {
                spde.prefix <- my.spde.prefix
                ## inla checks PREFIX valididy by looking for "s"
                model$param.inla <- list(
                    s=model$mesh$loc,
                    c0=model$internal$c0,
                    g1=model$internal$g1,
                    g2=model$internal$g2,
                    basis.T=model$internal$basis.T,
                    basis.K=model$internal$basis.K)
                spde.matrices <- names(model$param.inla)
            } else if (any(inherits(model, "inla.spde2"))) {
                spde2.prefix <- my.spde.prefix
                spde.matrices <- c("M0", "M1", "M2",
                                   "B0", "B1", "B2", "BLC")
            } else if (any(inherits(model, "inla.spde3"))) {
                spde3.prefix <- my.spde.prefix
                spde.matrices <- c("M0", "M1", "M2", "M3",
                                   "B0", "B1", "B2", "B3", "BLC")
            }
            for (spde.matrix.name in spde.matrices) {
                ## Only write matrix if it is non-empty (only happens for BLC)
                if (nrow(model$param.inla[[spde.matrix.name]]) > 0) {
                    fmesher.write(inla.affirm.double(
                        model$param.inla[[spde.matrix.name]]),
                        my.spde.prefix, spde.matrix.name)
                }
            }
        }

        atmp = paste(unlist(lapply(as.list(formals(INLA::f)), function(x) names(x))))
        arguments = unique(sort(c(names(formals(INLA::f)), atmp[-which(nchar(atmp) == 0L)])))
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

    if (!is.null(constr)) {
        if (inla.one.of(model, c("spde3", "spde2", "spde")) && constr) {
            stop("Option 'constr=TRUE' is disabled for model='spde2' and 'spde' and 'spde3'; please refer to the spde-tutorial.")
        }
    }

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
    }

    if (!inla.one.of(model, c("copy", "clinear", "mec", "meb"))) {
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
    ## flag an error its not among the legal ones.
    atmp = paste(unlist(lapply(as.list(formals(INLA::f)), function(x) names(x))))
    arguments = unique(sort(c(names(formals(INLA::f)), atmp[-which(nchar(atmp) == 0L)])))
    arguments = arguments[-grep("^[.][.][.]$", arguments)]
    for(elm in args.eq) {
        if (!is.element(elm, arguments)) {
            f.call = as.character(as.expression(match.call(expand.dots=TRUE)))
            valid.args = inla.paste(sort(arguments), sep="\n\t")
            stop(paste("Argument `", elm, "' in formula specification\n\n\t\t",
                       f.call, "\n\n  is invalid. Valid arguments are:\n\n\t", valid.args, sep=""))
        }
    }

    ## check that 'order' is defined
    if (inla.one.of(model, c("ar"))) {
        if (is.null(order) || missing(order) || order < 1L) {
            stop("Model 'ar' needs 'order' to be defined as an integer > 0.")
        }
        order = as.integer(order)
        max.order = length(inla.models()$latent$ar$hyper) -1L
        if (order > max.order) {
            stop(paste("Model 'ar': order=", order, ", is to large. max.order =", max.order, sep=""))
        }
    }
    if (inla.one.of(model, c("fgn", "fgn2"))) {
        if (is.null(order) || missing(order)) {
            order = inla.models()$latent$fgn$order.default
        } else {
            order = as.integer(order)
        }
        stopifnot(any(order == inla.models()$latent$fgn$order.defined))
    }

    ## Check that the Cmatrix is defined for those models needing it, and oposite.
    if (!inla.one.of(model, "z")) {
        if (inla.one.of(model, c("generic", "generic0","generic1", "generic2"))) {
            if (is.null(Cmatrix)) {
                stop("For generic models the Cmatrix has to be provided")
            }
            Cmatrix = inla.sparse.check(Cmatrix)
            if (is.null(n)) {
                n = inla.sparse.dim(Cmatrix)[1]
            }
        } else if (inla.one.of(model, "generic3")) {
            if (!is.list(Cmatrix) || (is.list(Cmatrix) && length(Cmatrix) == 0L)) {
                stop("For model='generic3' then argument 'Cmatrix' must be a non-empty list of matrices.")
            }
            stopifnot(length(Cmatrix) > 0L)
            ## check each matrix
            Cdim = numeric(length(Cmatrix))
            for(i in 1:length(Cmatrix)) {
                Cmatrix[[i]] = inla.sparse.check(Cmatrix[[i]], must.be.squared = TRUE)
                Cdim[i] = inla.sparse.dim(Cmatrix[[i]])[1L]
            }
            stopifnot(all(Cdim == Cdim[1L]))
            if (is.null(n)) {
                n = Cdim[1L]
            }
        } else {
            if (!is.null(Cmatrix)) {
                stop(paste("Cmatrix is not used for this model:",  model))
            }
        }
    } else {
        ## for the z-model, then Cmatrix is required but Cmatrix=NULL is allowed; so we do not
        ## need to do anything at this stage.
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
    if (inla.one.of(model, c("besag", "bym", "bym2", "besagproper", "besagproper2"))) {
        if (is.null(graph)) {
            stop(paste("The 'graph' has to be provided for model", model))
        }
        n.from.graph = inla.graph.size(graph)
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
        n.from.graph = 2L*inla.graph.size(graph)
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

    if (inla.one.of(model, c("z"))) {
        if (is.null(Z)) {
            stop("With model [z] then covariate-matrix Z is required. Example: f(ind, Z=Z, model=\"z\")")
        }
        if (!is.null(rankdef) && rankdef != 0) {
            stop("Option rankdef!=0, is not allowed for model='z'; use 'constr' or 'extraconstr' to define intrinsic models from proper ones.")
        }
        if (is.null(n)) {
            n = sum(dim(Z))
        }
        if (!is.null(constr) && constr == TRUE) {
            ## let constr=TRUE be defined as sum(z)=0 only.
            constr=FALSE
            zn = dim(Z)[1L]
            zm = dim(Z)[2L]
            z.row = c(rep(0, zn), rep(1, zm))
            if (empty.extraconstr(extraconstr)) {
                extraconstr = list(A = matrix(z.row, 1, zn+zm), e=0)
            } else {
                extraconstr$A = rbind(extraconstr$A, z.row)
                extraconstr$e = c(extraconstr$e, 0)
            }
        }
    }

    if (inla.one.of(model, c("slm"))) {
        stopifnot(!is.null(args.slm))
        stopifnot(!is.null(args.slm$rho.min) && length(args.slm$rho.min) == 1L)
        stopifnot(!is.null(args.slm$rho.max) && length(args.slm$rho.max) == 1L)
        stopifnot(args.slm$rho.max > args.slm$rho.min)
        stopifnot(!is.null(args.slm$X) && inla.is.matrix(args.slm$X))
        stopifnot(!is.null(args.slm$W) && inla.is.matrix(args.slm$W))
        stopifnot(!is.null(args.slm$Q.beta) && inla.is.matrix(args.slm$Q.beta))

        ## make sure they are sparse
        args.slm$X = inla.as.sparse(args.slm$X)
        args.slm$W = inla.as.sparse(args.slm$W)
        args.slm$Q.beta = inla.as.sparse(args.slm$Q.beta)

        slm.n = dim(args.slm$X)[1L]
        slm.m = dim(args.slm$X)[2L]
        stopifnot(all(dim(args.slm$X) == c(slm.n, slm.m)))
        stopifnot(all(dim(args.slm$W) == c(slm.n, slm.n)))
        stopifnot(all(dim(args.slm$Q.beta) == c(slm.m, slm.m)))

        if (missing(n) || is.null(n)) {
            n = slm.n + slm.m
        } else {
            stopifnot(n == slm.n + slm.m)
        }
    }

    if (inla.one.of(model, c("ar1c"))) {
        stopifnot(!is.null(args.ar1c))
        stopifnot(!is.null(args.ar1c$Z) && inla.is.matrix(args.ar1c$Z))
        stopifnot(!is.null(args.ar1c$Q.beta) && inla.is.matrix(args.ar1c$Q.beta))

        args.ar1c$Z = as.matrix(args.ar1c$Z)           ## is dense
        args.ar1c$Q.beta = as.matrix(args.ar1c$Q.beta) ## is dense

        ar1c.n = dim(args.ar1c$Z)[1L]
        ar1c.m = dim(args.ar1c$Z)[2L]
        stopifnot(all(dim(args.ar1c$Z) == c(ar1c.n, ar1c.m)))
        stopifnot(all(dim(args.ar1c$Q.beta) == c(ar1c.m, ar1c.m)))

        if (missing(n) || is.null(n)) {
            n = ar1c.n + ar1c.m
        } else {
            stopifnot(n == ar1c.n + ar1c.m)
        }
    }

    if (inla.one.of(model, c("dmatern"))) {
        stopifnot(!missing(locations) && !is.null(locations))
        if (is.vector(locations)) {
            locations = matrix(locations, ncol = 1)
        }
        stopifnot(is.matrix(locations))
        stopifnot(nrow(locations) > 1)
        stopifnot(!any(is.na(locations)))
        n = nrow(locations)
    } else {
        stopifnot(missing(locations) || is.null(locations))
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
        ## If spde2.prefix is given, use it.
        if (is.null(spde2.prefix)) {
            stop("Argument spde2.prefix=NULL is required for model = spde2")
        }
    }
    spde2.transform = match.arg(spde2.transform, several.ok = FALSE)

    if (inla.one.of(model, c("spde3"))) {
        if (is.null(spde3.prefix)) {
            stop("Argument spde3.prefix=NULL is required for model = spde3")
        }
    }
    spde3.transform = match.arg(spde3.transform, several.ok = FALSE)

    if (inla.one.of(model,"seasonal") &&
        is.null(season.length)) {
        stop("The length of the season has to be provided in season.length")
    }

    ## cyclic is only valid for rw1, rw2 and rw2d-models
    if (!is.null(cyclic) && cyclic &&
        !inla.one.of(model, c("rw1", "rw2", "rw2d", "rw2diid", "ar1"))) {
        stop("Cyclic defined only for rw1, rw1c2, rw2, rw2d, rw2diid and ar1 models")
    }

    need.nrow.ncol = inla.model.properties(model, "latent")$nrow.ncol
    ## nrow/ncol
    if ((!is.null(nrow) || !is.null(ncol)) && !need.nrow.ncol) {
        stop(paste("'nrow' and 'ncol' are not needed for model = ", model))
    }
    if (need.nrow.ncol) {
        if (is.null(nrow) || is.null(ncol)) {
            stop(paste("'nrow' and 'ncol' must be specified for model", model))
        }
        if (nrow <= 0 || ncol <= 0 || trunc(nrow) != nrow || trunc(ncol) != ncol) {
            stop("'nrow' and 'ncol' must be positive intergers.")
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
        if (constr || !empty.extraconstr(extraconstr)) {
            diagonal = inla.set.f.default()$diagonal
        } else {
            diagonal = 0.0
        }
    }

    if (!empty.extraconstr(extraconstr)) {
        ##check
        A=extraconstr$A
        e=extraconstr$e
        if (!is.matrix(A)) {
            stop("A(extraconstr) has to be a matrix")
        } else {
            if (nrow(A)!=length(e)) {
                stop("Dimension of A and e do not correspond")
            }
        }
    }

    if (!missing(scale.model) && !inla.one.of(model, c("rw1", "rw2", "besag", "bym", "bym2", "besag2", "rw2d", "rw2diid", "seasonal"))) {
        stop(paste("Option 'scale.model' is not used for model:", model))
    }
    if (missing(scale.model) || is.null(scale.model)) {
        ## must doit like this otherwise we run into problems when
        ## compiling the package
        scale.model = inla.getOption("scale.model.default")
        if (inla.one.of(model, c("bym2", "rw2diid"))) {
            scale.model = TRUE
        }
    }
    if (inla.one.of(model, c("bym2", "rw2diid")) && !scale.model) {
        stop("Model 'bym2' and 'rw2diid' require scale.model=TRUE or 'missing(scale.model)'")
    }

    if (inla.one.of(model, c("besag", "besag2", "bym", "bym2"))) {
        ## this is a somewhat complicated case
        g = inla.read.graph(graph)
        cc.n = sapply(g$cc$nodes, length)
        cc.n1 = sum(cc.n == 1L)
        cc.n2 = sum(cc.n >= 2L)

        if (is.null(rankdef)) {
            if (scale.model) {
                rankdef = cc.n2
            } else {
                rankdef = cc.n2 + cc.n1
            }
            if (!(empty.extraconstr(extraconstr))) {
                rankdef = rankdef + dim(extraconstr$A)[1]
            }
        }
            
        if (adjust.for.con.comp) {
            ## we need to check if the graph contains more than 1 connected components. 
            if (g$cc$n == 1) {
                ## nothing to do 
            } else {
                if (debug) {
                    print(paste("modify model", model))
                    print(paste("number of connected components", inla.paste(cc.n)))
                    print(paste("number of connected components of size >= 2L", cc.n2))
                }
                if (constr) {
                    ## need to redefine the meaning of constr = TRUE to mean
                    ## constr=TRUE for all connected components with size >= 2L

                    ## like bym/bym2 place the constr on the second half
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

                    if (!empty.extraconstr(extraconstr)) {
                        extraconstr$A = rbind(AA, extraconstr$A)
                        extraconstr$e = c(ee, extraconstr$e)
                    } else {
                        extraconstr = list(A = AA, e = ee)
                    }

                    if (debug) {
                        print(AA)
                        print(ee)
                    }
                } else {
                    if (!missing(adjust.for.con.comp)) {
                        stop("The option 'adjust.for.con.comp' is only used for models 'besag', 'besag2', 'bym' and 'bym2'.")
                    }
                }
                if (is.null(diagonal)) {
                    diagonal = inla.set.f.default()$diagonal
                    if (debug) {
                        print(paste("set diagonal = ", diagonal))
                    }
                }
            }
        } else {
            ## nothing to do
        }
    }

    if (model %in% "rgeneric") {
        if (inla.is.element("f", rgeneric) && inla.is.element("rgeneric", rgeneric$f)) {
            rgeneric = rgeneric$f$rgeneric
        }
        R.init = rgeneric$R.init
        stopifnot(inherits(rgeneric, "inla.rgeneric"))
        ## add an 'Id' so we know who we are
        rgeneric = list(model = rgeneric, Id = vars[[1]], R.init = R.init)
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
        extraconstr= inla.ifelse(empty.extraconstr(extraconstr), NULL,  extraconstr),
        graph=graph,
        group = group,
        hyper = hyper,
        label = term,
        model=model,
        n = n,
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
        spde.prefix = spde.prefix,
        spde2.prefix = spde2.prefix,
        spde2.transform = spde2.transform,
        spde3.prefix = spde3.prefix,
        spde3.transform = spde3.transform,
        term=term,
        values=values,
        order = order,
        weights=weights,
        scale = scale,
        strata = strata,
        rgeneric = rgeneric,
        scale.model = as.logical(scale.model),
        adjust.for.con.comp = as.logical(adjust.for.con.comp),
        args.slm = args.slm,
        args.ar1c = args.ar1c,
        correct = correct,
        locations = locations
        )

    if (debug) print(ret)
    return (ret)
}

## "inla.model.class" is a generic model class that can be inherited
## from to mark a class as an inla model object class.
`inla.model.object.classes` = function()
{
    return (c("inla.model.class", "inla.wrapper.model",
              "inla.spde", "inla.spde1", "inla.spde2", "inla.spde3",
              "inla.rgeneric"))
}

## List spde classes that f() can write to disk when needed.
`inla.spde.object.classes` = function()
{
    return (c("inla.spde1", "inla.spde2", "inla.spde3"))
}
