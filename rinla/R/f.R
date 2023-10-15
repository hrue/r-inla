#' @export
#' @title Define general Gaussian models in the INLA formula
#' @description
#'
#'   Function used for defining of smooth and spatial terms within `inla` model
#'   formulae. The function does not evaluate anything - it
#'   exists purely to help set up a model. The function specifies one
#'   smooth function in the linear predictor (see [inla.list.models()]) as
#'   \deqn{w\ f(x)}{weight*f(var)}

`f` <- function(
                #' @param ... Name of the covariate and, possibly of the
                #' weights vector. NB: order counts!!!! The first specified
                #' term is the covariate and the second one is the vector of
                #' weights (which can be negative).
                ...,

                #' @param model A string indicating the chosen model. The
                #'  default is `iid`. See
                #'  `names(inla.models()$latent)` for a list of possible
                #'  alternatives and [inla.doc()] for detailed docs.
                model = "iid",

                #' @param copy The name of the model-component to copy
                copy = NULL,

                #' @param scopy The name of the model-component to smooth-copy (where the
                #' copy-function is a spline)
                scopy = NULL,

                #' @param same.as Can be used with `copy=".."`. `same.as="A"` says
                #' that this `copy` should use the same scaling parameter as another
                #' `copy` "A"
                same.as = NULL,

                #' @param n An optional argument which defines the dimension
                #' of the model if this is different from
                #' `length(sort(unique(covariate)))`
                n = NULL,

                #' @param nrep Number of replications, if not given,  then `nrep=max(replications)`
                nrep = NULL,

                #' @param replicate A vector of which replications to use.
                replicate = NULL,

                #' @param ngroup Number of groups,  if not given,  then `ngroup=max(group)`
                ngroup = NULL,

                #' @param group A vector of which groups to use.
                group = NULL,

                #' @param control.group Controls the use of `group`
                control.group = inla.set.control.group.default(),

                #' @param control.scopy Controls the use of `scopy`
                control.scopy = inla.set.control.scopy.default(),

                #' @param hyper Specification of the hyperparameter, fixed or
                #' random, initial values, priors and its parameters. See
                #' `?inla.models` for the list of hyparameters for each
                #' model and its default options or
                #' use `inla.doc()` for
                #' detailed info on the family and
                #' supported prior distributions.
                hyper = NULL,

                #' @param initial THIS OPTION IS OBSOLETE,  DO NOT USE
                initial = NULL,

                #' @param prior THIS OPTION IS OBSOLETE,  DO NOT USE
                prior = NULL,

                #' @param param THIS OPTION IS OBSOLETE,  DO NOT USE
                param = NULL,

                #' @param fixed THIS OPTION IS OBSOLETE; DO NOT USE
                fixed = NULL,

                #' @param season.length Length of the seasonal component for `model="seasonal"`
                season.length = NULL,

                #' @param constr  A boolean variable indicating whater to set
                #' a sum to 0 constraint on the term. By default the sum to 0
                #' constraint is imposed on all intrinsic models
                #' ("iid","rw1","rw1","besag", etc..).
                constr = NULL,

                #' @param extraconstr This argument defines extra linear
                #' constraints. The argument is a list with two elements, a
                #' matrix `A` and a vector `e`, which defines the
                #' extra constraint `Ax = e`; for example
                #' `extraconstr = list(A = A, e=e)`. The number of
                #' columns of `A` must correspond to the length of this
                #' `f`-model.  Note that this constraint comes
                #' additional to the sum-to-zero constraint defined if
                #' `constr = TRUE`.
                extraconstr = list(A = NULL, e = NULL),

                #' @param values An optional vector giving all values
                #' assumed by the covariate for which we want estimated the
                #' effect. It must be a numeric vector, a vector of factors
                #' or `NULL`.
                values = NULL,

                #' @param cyclic A boolean specifying wheather the model is
                #' cyclical. Only valid for "rw1" and "rw2" models, is
                #' cyclic=T then the sum to 0 constraint is removed. For the
                #' correct form of the grah file see \cite{Martino and Rue
                #' (2008)}.
                cyclic = NULL,

                #' @param diagonal An extra constant added to the diagonal of
                #' the precision matrix to prevent numerical issues.
                diagonal = NULL,

                #' @param graph Defines the graph-object either as a file with
                #' a graph-description, an `inla.graph`-object, or as a
                #' (sparse) symmetric matrix .
                graph = NULL,

                #' @param graph.file THIS OPTION IS OBSOLETE,  DO NOT USE
                graph.file = NULL,

                #' @param cdf THIS OPTION IS OBSOLETE,  DO NOT USE
                cdf = NULL,

                #' @param quantiles A vector of maximum 10 quantiles,
                #' \eqn{p(0), p(1),\dots}{p(0), p(1),\ldots} to compute for
                #' each posterior marginal. The function returns, for each
                #' posterior marginal, the values
                #' \eqn{x(0), x(1),\dots}{x(0), x(1),\ldots} such that
                #' \deqn{\mbox{Prob}(X<x(p))=p}{Prob(X<x)=p}
                quantiles = NULL,

                #' @param Cmatrix The specification of the precision matrix
                #' for the generic,  generic3 or z models (up to a scaling constant).
                #' `Cmatrix` is either a
                #' (dense) matrix, a matrix created using
                #' `Matrix::sparseMatrix()`, or a filename which stores the
                #' non-zero elements of `Cmatrix`, in three columns:
                #' `i`, `j` and `Qij`. In case of the generic3 model,
                #' it is a list of such specifications.
                Cmatrix = NULL,

                #' @param rankdef A number **defining** the rank
                #' deficiency of the model, with sum-to-zero constraint and
                #' possible extra-constraints taken into account. See
                #' details.
                rankdef = NULL,

                #' @param Z The matrix for the z-model
                Z = NULL,

                #' @param nrow Number of rows for 2d-models
                nrow = NULL,

                #' @param ncol Number of columns for 2d-models
                ncol = NULL,

                #' @param nu Smoothing parameter for the Matern2d-model,
                #' possible values are `c(0, 1, 2, 3)`
                nu = NULL,

                #' @param bvalue The boundary conditions for model `rw2d`,  0 means use
                #' the correct subspace (default),  while 1 means condition on 0's outside
                bvalue = NULL,

                #' @param spde.prefix Internal use only
                spde.prefix = NULL,

                #' @param spde2.prefix Internal use only
                spde2.prefix = NULL,

                #' @param spde2.transform Internal use only
                spde2.transform = c("logit", "log", "identity"),

                #' @param spde3.prefix Internal use only
                spde3.prefix = NULL,

                #' @param spde3.transform Internal use only
                spde3.transform = c("logit", "log", "identity"),

                #' @param mean.linear Prior mean for `model="linear"`
                mean.linear = inla.set.control.fixed.default()$mean,

                #' @param prec.linear Prior precision for `model="linear"`
                prec.linear = inla.set.control.fixed.default()$prec,

                #' @param compute  A boolean variable indicating whether the
                #'  marginal posterior distribution for the nodes in the
                #'  `f()` model should be computed or not. This is
                #'  usefull for large models where we are only interested in
                #'  some posterior marginals.
                compute = TRUE,

                #' @param of Internal use only
                of = NULL,

                #' @param precision The precision for the artificial noise added when creating a copy of a model and others.
                precision = 10^8,

                #' @param range A vector of size two giving the lower and
                #' upper range for the scaling parameter `beta` in the
                #' model `COPY`, `CLINEAR`, `MEC` and `MEB`.
                #' If `low = high` then the identity mapping
                #' is used.
                range = NULL,

                #' @param adjust.for.con.comp If TRUE (default), adjust some
                #' of the models (currently: besag, bym, bym2 and besag2) if the
                #' number of connected components in graph is larger than 1.
                #' If FALSE, do nothing.
                adjust.for.con.comp = TRUE,

                #' @param order Defines the `order` of the model: for
                #' model `ar` this defines the order p, in AR(p). Not
                #' used for other models at the time being.
                order = NULL,

                #' @param scale A scaling vector. Its meaning depends on the model.
                scale = NULL,

                #' @param rgeneric A object of class `inla.rgeneric` which defines the model. (EXPERIMENTAL!)
                rgeneric = NULL,

                #' @param cgeneric A object of class `inla.cgeneric` which defines the model. (EXPERIMENTAL!)
                cgeneric = NULL,

                #' @param scale.model Logical. If `TRUE` then scale the RW1 and RW2 and BESAG and BYM and BESAG2 and RW2D models so the their (generlized) variance is 1. Default value is `inla.getOption("scale.model.default")`
                scale.model = NULL,

                #' @param args.slm Required arguments to the model="slm"; see the documentation for further details.
                args.slm = list(rho.min = NULL, rho.max = NULL, X = NULL, W = NULL, Q.beta = NULL),

                #' @param args.ar1c Required arguments to the model="ar1c"; see the documentation for further details.
                args.ar1c = list(Z = NULL, Q.beta = NULL),

                #' @param args.intslope A list with the `subject` (factor),  `strata` (factor) and `covariates` (numeric) for the `intslope` model; see the documentation for further details,
                args.intslope = list(subject = NULL, strata = NULL, covariates = NULL),

                #' @param vb.correct Add this model component to the list of nodes to be used for the (potential) vb correction?  If `TRUE` do,  and do not if `FALSE`. Can also be a vector of nodes to add in the correction-set.
                vb.correct = TRUE,

                #' @param locations A matrix with locations for the model `dmatern`. This also defines `n`.
                locations = NULL,

                #' @param debug Enable local debug output
                debug = FALSE) {
    #'
    #' @returns TODO

    #' @details There is no default value for `rankdef`, if it
    #' is not defined by the user then it is computed by the rank
    #' deficiency of the prior model (for the generic model, the
    #' default is zero), plus 1 for the sum-to-zero constraint if the
    #' prior model is proper, plus the number of extra
    #' constraints. **Oops:** This can be wrong, and then the user
    #' must define the `rankdef` explicitly.
    #' @author Havard Rue \email{hrue@@r-inla.org}
    #' @seealso [inla()], [hyperpar.inla()]

    ## this is required. the hyper.defaults can only be changed in the
    ## model=model.object
    hyper.default <- NULL

    empty.extraconstr <- function(ec) {
        return(is.null(ec) || all(sapply(ec, is.null)))
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
                    s = model$mesh$loc,
                    c0 = model$internal$c0,
                    g1 = model$internal$g1,
                    g2 = model$internal$g2,
                    basis.T = model$internal$basis.T,
                    basis.K = model$internal$basis.K
                )
                spde.matrices <- names(model$param.inla)
            } else if (any(inherits(model, "inla.spde2"))) {
                spde2.prefix <- my.spde.prefix
                spde.matrices <- c(
                    "M0", "M1", "M2",
                    "B0", "B1", "B2", "BLC"
                )
            } else if (any(inherits(model, "inla.spde3"))) {
                spde3.prefix <- my.spde.prefix
                spde.matrices <- c(
                    "M0", "M1", "M2", "M3",
                    "B0", "B1", "B2", "B3", "BLC"
                )
            }
            for (spde.matrix.name in spde.matrices) {
                ## Only write matrix if it is non-empty (only happens for BLC)
                if (nrow(model$param.inla[[spde.matrix.name]]) > 0) {
                    fmesher.write(
                        inla.affirm.double(
                            model$param.inla[[spde.matrix.name]]
                        ),
                        my.spde.prefix, spde.matrix.name
                    )
                }
            }
        }

        atmp <- paste(unlist(lapply(as.list(formals(INLA::f)), function(x) names(x))))
        arguments <- unique(sort(c(names(formals(INLA::f)), atmp[-which(nchar(atmp) == 0L)])))
        arguments <- arguments[-grep("^[.][.][.]$", arguments)]
        ## add this one manually
        arguments <- c(arguments, "hyper.default")

        ## evaluate arguments in 'model' and set those that are in
        ## 'arguments'. However, due to the 'missing || is.null', the
        ## argument 'model' does not pass the test and must be set
        ## manually afterwards.
        for (nm in names(model$f)) {
            if (inla.one.of(nm, arguments)) {
                inla.eval(paste("if (missing(", nm, ") || is.null(", nm, ")) ",
                                nm, "=", "model$f$", nm,
                                sep = ""
                                ))
            } else {
                ## this should not happen
                stop(paste("Argument `", nm, "' is not an argument in f().",
                           " This is likely not what you want!",
                           sep = ""
                           ))
            }
        }

        ## the 'model' argument is required!
        stopifnot(is.character(model$f$model))
        model <- model$f$model
    }

    stopifnot(!(!is.null(copy) && !is.null(scopy)))
    if (!is.null(cdf)) {
        stop("The argument 'cdf' is no longer in use. Please fix.")
    }

    ## this is a nice trick
    if (!is.null(copy)) {
        if (!missing(model)) {
            warning(paste("Ignored argument model=`", model,
                          "' in f() due to copy=`", copy, "'",
                          sep = ""
                          ))
        }
        if (!is.null(of)) {
            stop("Argument `of=NULL' is required when `copy=...' is used.")
        }
        model <- "copy"
        of <- copy
        copy <- NULL
    }
    if (!is.null(scopy)) {
        if (!missing(model)) {
            warning(paste("Ignored argument model=`", model,
                          "' in f() due to scopy=`", scopy, "'",
                          sep = ""
                          ))
        }
        if (!is.null(of)) {
            stop("Argument `of=NULL' is required when `scopy=...' is used.")
        }
        model <- "scopy"
        of <- scopy
        scopy <- NULL
    }

    if (is.null(model)) {
        stop("No model is specified.")
    }
    inla.is.model(model, "latent", stop.on.error = TRUE)

    if (!is.null(constr)) {
        if (inla.one.of(model, c("spde3", "spde2", "spde")) && constr) {
            stop("Option 'constr=TRUE' is disabled for model='spde2' and 'spde' and 'spde3'; please refer to the spde-tutorial.")
        }
    }

    ## in ... is the name of the covariate  and possibly the location of the weights
    ## like f(covariate, weights)
    vars <- as.list(substitute(list(...)))[-1]
    d <- length(vars)
    if (d == 0L) {
        stop(paste("Missing covariate in f(...) for model=", model))
    }
    term <- deparse(vars[[1]], backtick = TRUE, width.cutoff = 500)
    if (debug) {
        print(vars)
    }
    ## the second term in ... is the (possible) weights for the selected covariate!
    if (d == 1) {
        weights <- NULL
    } else if (d == 2) {
        weights <- deparse(vars[[2]], backtick = TRUE, width.cutoff = 500)
    } else if (d > 2) {
        stop(paste("To many variables included in f():", inla.paste(vars)))
    } else if (d == 0) {
        stop("At least one variable in f() needs to be defined")
    }
    ## get the weights
    term <- attr(terms(reformulate(term)), "term.labels")
    if (d > 1) {
        weigths <- attr(terms(reformulate(weights)), "weights.labels")
    }

    ## set the hyperparameters
    hyper <- inla.set.hyper(
        model = model, section = "latent",
        hyper = hyper, hyper.default = hyper.default,
        initial = initial, fixed = fixed, prior = prior, param = param
    )

    ## for model = copy, its is not allowed to define constr or extraconstr
    if (inla.one.of(model, c("copy", "scopy"))) {
        stopifnot(missing(constr))
        stopifnot(missing(extraconstr))
        ## unless stated otherwise, do not vb.correct copies
        if (is.null(vb.correct)) {
            vb.correct <- FALSE
        }
    }

    if (!inla.one.of(model, c("copy", "clinear", "mec", "meb"))) {
        stopifnot(missing(range))
    }

    if (!missing(same.as) && !is.null(same.as) && model != "copy") {
        stop(paste("Argument 'same.as' must be NULL is model != 'copy' :", model))
    }

    inla.check.control(control.group)
    cont.group <- inla.set.control.group.default()
    cont.group[(namc <- names(control.group))] <- control.group
    cont.group$hyper <- inla.set.hyper(
        cont.group$model, "group", cont.group$hyper,
        cont.group$initial, cont.group$fixed, cont.group$prior, cont.group$param
    )

    inla.check.control(control.scopy)
    cont.scopy <- inla.set.control.scopy.default()
    cont.scopy[(namc <- names(control.scopy))] <- control.scopy
    cont.scopy$hyper <- inla.set.hyper(
        cont.scopy$model, "scopy", cont.scopy$hyper,
        cont.scopy$initial, cont.scopy$fixed, cont.scopy$prior, cont.scopy$param
    )

    ## CHECK ARGUMENTS.
    ## This is a bit tricky. We want to check if there are arguments
    ## in `...' which is of type `name = value' which are invalid. If
    ## so, this lead to an obscoure error later on...
    ##
    ## we first collect all arguments of type `name = value'
    if (TRUE) {
        ## New code
        args.eq <- names(match.call(expand.dots = TRUE))
        args.eq <- args.eq[args.eq != ""]
    } else {
        args.eq <- c()
        for (arg in unlist(strsplit(as.character(as.expression(match.call(expand.dots = TRUE))), ","))) {
            if (length(grep("=", arg)) > 0) {
                args.eq <- c(args.eq, gsub(" ", "", unlist(strsplit(arg, "="))[1]))
            }
        }
    }

    ## then we compare these with the legal ones in f(), and
    ## flag an error its not among the legal ones.
    atmp <- paste(unlist(lapply(as.list(formals(INLA::f)), function(x) names(x))))
    arguments <- unique(sort(c(names(formals(INLA::f)), atmp[-which(nchar(atmp) == 0L)])))
    arguments <- arguments[-grep("^[.][.][.]$", arguments)]
    for (elm in args.eq) {
        if (!is.element(elm, arguments)) {
            f.call <- as.character(as.expression(match.call(expand.dots = TRUE)))
            valid.args <- inla.paste(sort(arguments), sep = "\n\t")
            stop(paste("Argument `", elm, "' in formula specification\n\n\t\t",
                       f.call, "\n\n  is invalid. Valid arguments are:\n\n\t", valid.args,
                       sep = ""
                       ))
        }
    }

    ## check that 'order' is defined
    if (inla.one.of(model, c("ar"))) {
        if (is.null(order) || missing(order) || order < 1L) {
            stop("Model 'ar' needs 'order' to be defined as an integer > 0.")
        }
        order <- as.integer(order)
        max.order <- length(inla.models()$latent$ar$hyper) - 1L
        if (order > max.order) {
            stop(paste("Model 'ar': order=", order, ", is to large. max.order =", max.order, sep = ""))
        }
    }
    if (inla.one.of(model, c("fgn", "fgn2"))) {
        if (is.null(order) || missing(order)) {
            order <- inla.models()$latent$fgn$order.default
        } else {
            order <- as.integer(order)
        }
        stopifnot(any(order == inla.models()$latent$fgn$order.defined))
    }
    if (inla.one.of(model, c("iidkd"))) {
        if (is.null(order) || missing(order)) {
            stop("Model 'iidkd' needs 'order' to be defined")
        } else {
            order <- as.integer(order)
            if (!(order >= 2 && order <= 20)) {
                stop("Model 'iidkd' needs 'order' to be in the interval 2 to 20")
            }
        }
    }

    ## Check that the Cmatrix is defined for those models needing it, and oposite.
    if (!inla.one.of(model, "z")) {
        if (inla.one.of(model, c("generic", "generic0", "generic1", "generic2"))) {
            if (is.null(Cmatrix)) {
                stop("For generic models the Cmatrix has to be provided")
            }
            Cmatrix <- inla.sparse.check(Cmatrix)
            if (is.null(n)) {
                n <- inla.sparse.dim(Cmatrix)[1]
            }
        } else if (inla.one.of(model, "generic3")) {
            if (!is.list(Cmatrix) || (is.list(Cmatrix) && length(Cmatrix) == 0L)) {
                stop("For model='generic3' then argument 'Cmatrix' must be a non-empty list of matrices.")
            }
            stopifnot(length(Cmatrix) > 0L)
            ## check each matrix
            Cdim <- numeric(length(Cmatrix))
            for (i in 1:length(Cmatrix)) {
                Cmatrix[[i]] <- inla.sparse.check(Cmatrix[[i]], must.be.squared = TRUE)
                Cdim[i] <- inla.sparse.dim(Cmatrix[[i]])[1L]
            }
            stopifnot(all(Cdim == Cdim[1L]))
            if (is.null(n)) {
                n <- Cdim[1L]
            }
        } else {
            if (!is.null(Cmatrix)) {
                stop(paste("Cmatrix is not used for this model:", model))
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
        graph <- graph.file

        ## show a warning,  if this hasn't been shown before.
        if (inla.getOption("show.warning.graph.file")) {
            warning("Argument 'graph.file' in 'f()' is obsolete; please use the more general argument 'graph' instead.")
            ## disable this warning from now on.
            inla.setOption(show.warning.graph.file = FALSE)
        }
    }

    ## chech that the graph is provided, if required. Set 'n' from the graph.
    if (inla.one.of(model, c("besag", "bym", "bym2", "besagproper", "besagproper2"))) {
        if (is.null(graph)) {
            stop(paste("The 'graph' has to be provided for model", model))
        }
        n.from.graph <- inla.graph.size(graph)
        if (n.from.graph <= 0) {
            stop(paste("Argument 'n from graph' is void:", n.from.graph))
        }
        if (!is.null(n) && n != n.from.graph) {
            stop(paste("Argument 'n' and 'n from graph' does not match", n, n.from.graph))
        }
        n <- n.from.graph
    }
    if (inla.one.of(model, c("besag2"))) {
        if (is.null(graph)) {
            stop(paste("The graph has to be provided for model", model))
        }
        ## read n from the graph
        nn <- inla.graph.size(graph)
        n.from.graph <- 2L * nn
        if (n.from.graph <= 0) {
            stop(paste("Argument 'n from graph' is void:", n.from.graph))
        }
        if (!is.null(n) && n != n.from.graph) {
            stop(paste("Argument 'n' and 2*'n from graph' does not match", n, n.from.graph))
        }
        n <- n.from.graph

        ## this is a special case
        if (is.null(constr) || constr) {
            constr <- FALSE
            con <- list(
                A = rbind(
                    rep(1:0, each = nn),
                    rep(0:1, each = nn)
                ),
                e = c(0, 0)
            )
            extraconstr$A <- rbind(extraconstr$A, con$A)
            extraconstr$e <- c(extraconstr$e, con$e)
            if (missing(rankdef) || is.null(rankdef)) {
                rankdef <- nrow(extraconstr$A) %/% 2L
            }
        }
        if (is.null(diagonal)) {
            ## need a larger value...
            diagonal <- inla.set.f.default()$diagonal * 10
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
            n <- sum(dim(Z))
        }
        if (!is.null(constr) && constr == TRUE) {
            ## let constr=TRUE be defined as sum(z)=0 only.
            constr <- FALSE
            zn <- dim(Z)[1L]
            zm <- dim(Z)[2L]
            z.row <- c(rep(0, zn), rep(1, zm))
            if (empty.extraconstr(extraconstr)) {
                extraconstr <- list(A = matrix(z.row, 1, zn + zm), e = 0)
            } else {
                extraconstr$A <- rbind(extraconstr$A, z.row)
                extraconstr$e <- c(extraconstr$e, 0)
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
        args.slm$X <- inla.as.sparse(args.slm$X)
        args.slm$W <- inla.as.sparse(args.slm$W)
        args.slm$Q.beta <- inla.as.sparse(args.slm$Q.beta)

        slm.n <- dim(args.slm$X)[1L]
        slm.m <- dim(args.slm$X)[2L]
        stopifnot(all(dim(args.slm$X) == c(slm.n, slm.m)))
        stopifnot(all(dim(args.slm$W) == c(slm.n, slm.n)))
        stopifnot(all(dim(args.slm$Q.beta) == c(slm.m, slm.m)))

        if (missing(n) || is.null(n)) {
            n <- slm.n + slm.m
        } else {
            stopifnot(n == slm.n + slm.m)
        }
    }

    if (inla.one.of(model, c("ar1c"))) {
        stopifnot(!is.null(args.ar1c))
        stopifnot(!is.null(args.ar1c$Z) && inla.is.matrix(args.ar1c$Z))
        stopifnot(!is.null(args.ar1c$Q.beta) && inla.is.matrix(args.ar1c$Q.beta))

        args.ar1c$Z <- as.matrix(args.ar1c$Z) ## is dense
        args.ar1c$Q.beta <- as.matrix(args.ar1c$Q.beta) ## is dense

        ar1c.n <- dim(args.ar1c$Z)[1L]
        ar1c.m <- dim(args.ar1c$Z)[2L]
        stopifnot(all(dim(args.ar1c$Z) == c(ar1c.n, ar1c.m)))
        stopifnot(all(dim(args.ar1c$Q.beta) == c(ar1c.m, ar1c.m)))

        if (missing(n) || is.null(n)) {
            n <- ar1c.n + ar1c.m
        } else {
            stopifnot(n == ar1c.n + ar1c.m)
        }
    }

    if (inla.one.of(model, c("dmatern"))) {
        stopifnot(!missing(locations) && !is.null(locations))
        if (is.vector(locations)) {
            locations <- matrix(locations, ncol = 1)
        }
        stopifnot(is.matrix(locations))
        stopifnot(nrow(locations) > 1)
        stopifnot(!any(is.na(locations)))
        n <- nrow(locations)
    } else {
        stopifnot(missing(locations) || is.null(locations))
    }

    if (inla.one.of(model, "intslope")) {
        stopifnot(!is.null(args.intslope$subject))
        stopifnot(!is.null(args.intslope$strata))
        stopifnot(!is.null(args.intslope$covariates))
        zz <- args.intslope$covariates
        zz[is.na(zz)] <- 0
        args.intslope <- list(
            subject = as.numeric(as.factor(args.intslope$subject)),
            strata = as.numeric(as.factor(args.intslope$strata)),
            covariates = zz
        )
        if (is.null(n)) {
            n <- length(args.intslope$subject)
        } else {
            stopifnot(n == length(args.intslope$subject))
        }
        stopifnot(length(args.intslope$subject) == length(args.intslope$strata))
        stopifnot(length(args.intslope$subject) == length(args.intslope$covariates))
    }

    ## is N required?
    if (is.null(n) && (!is.null(inla.model.properties(model, "latent")$n.required)
        && inla.model.properties(model, "latent")$n.required)) {
        stop(paste("Argument `n' in f() is required for model:", model))
    }

    ## special N required?
    if (inla.one.of(model, "iidkd")) {
        n.div.by <- order
        if (!inla.divisible(n, n.div.by)) {
            stop(paste("Argument `n'", n, "is not divisible by", n.div.by))
        }
    } else if ((!is.null(inla.model.properties(model, "latent")$n.div.by)
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
        values <- 1:n
    }

    if (inla.one.of(model, "linear")) {
        if (d != 1L) {
            stop("Model = 'linear' do not accept weights. Just set 'z.new = z * weights' as the covariates.")
        }
        if (!is.null(group)) {
            stop("Model = 'linear' do not accept argument 'group'. Please use another model component.")
        }
        if (!is.null(replicate)) {
            stop("Model = 'linear' do not accept argument 'replicate'. Please use another model component.")
        }
        ret <- list(
            d = d, term = term, model = model, mean.linear = mean.linear, prec.linear = prec.linear, label = term,
            cdf = cdf, quantiles = quantiles, compute = compute
        )
        ## return here!
        return(ret)
    }
    if (!missing(prec.linear) || !missing(mean.linear)) {
        stop("Arguments 'mean.linear' and 'prec.linear' are defined only for model='linear'.")
    }

    if (inla.one.of(model, c("spde"))) {
        if (is.null(spde.prefix)) {
            stop("Argument spde.prefix=NULL is required for model = spde")
        }
        ## file ``PREFIXs'' must exists... test this one
        if (!(file.exists(paste(spde.prefix, "s", sep = "")))) {
            stop(paste("Argument spde.prefix=", spde.prefix, "does not seems to be valid (no file `PREFIXs')"))
        }
    }

    if (inla.one.of(model, c("spde2"))) {
        ## If spde2.prefix is given, use it.
        if (is.null(spde2.prefix)) {
            stop("Argument spde2.prefix=NULL is required for model = spde2")
        }
    }
    spde2.transform <- match.arg(spde2.transform, several.ok = FALSE)

    if (inla.one.of(model, c("spde3"))) {
        if (is.null(spde3.prefix)) {
            stop("Argument spde3.prefix=NULL is required for model = spde3")
        }
    }
    spde3.transform <- match.arg(spde3.transform, several.ok = FALSE)

    if (inla.one.of(model, "seasonal") &&
        is.null(season.length)) {
        stop("The length of the season has to be provided in season.length")
    }

    ## cyclic is only valid for rw1, rw2 and rw2d-models
    if (!is.null(cyclic) && cyclic &&
        !inla.one.of(model, c("rw1", "rw2", "rw2d", "rw2diid", "ar1"))) {
        stop("Cyclic defined only for rw1, rw1c2, rw2, rw2d, rw2diid and ar1 models")
    }

    need.nrow.ncol <- inla.model.properties(model, "latent")$nrow.ncol
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
            stopifnot(n == nrow * ncol)
        } else {
            n <- nrow * ncol
        }

        ## and set default values, if required
        if (inla.model.properties(model, "latent")$set.default.values) {
            if (missing(values)) {
                values <- 1:n
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

    ## for all instrinsic model the constraint has to be default ON
    if (is.null(constr)) {
        constr <- inla.model.properties(model, "latent")$constr
    }

    ## if diagonal is not set, the set this depending on the constr
    if (is.null(diagonal)) {
        if (constr || !empty.extraconstr(extraconstr)) {
            diagonal <- inla.set.f.default()$diagonal
        } else {
            diagonal <- 0.0
        }
    }

    if (!empty.extraconstr(extraconstr)) {
        ## check
        A <- extraconstr$A
        e <- extraconstr$e
        if ((!is.matrix(A) && !is(A, "Matrix"))) {
            stop("extraconstr$A has to be a dense or sparse matrix")
        } else {
            if (nrow(A) != length(e)) {
                stop("Dimension of A and e do not correspond")
            }
        }
    }

    if (!missing(scale.model) && !inla.one.of(model, c(
                                                         "rw1", "rw2", "besag", "bym",
                                                         "bym2", "besag2", "rw2d",
                                                         "rw2diid", "seasonal"
                                                     ))) {
        stop(paste("Option 'scale.model' is not used for model:", model))
    }
    if (missing(scale.model) || is.null(scale.model)) {
        ## must doit like this otherwise we run into problems when
        ## compiling the package
        scale.model <- inla.getOption("scale.model.default")
        if (inla.one.of(model, c("bym2", "rw2diid"))) {
            scale.model <- TRUE
        }
    }
    if (inla.one.of(model, c("bym2", "rw2diid")) && !scale.model) {
        stop("Model 'bym2' and 'rw2diid' require scale.model=TRUE or 'missing(scale.model)'")
    }

    if (inla.one.of(model, c("besag", "besag2", "bym", "bym2"))) {
        ## this is a somewhat complicated case
        g <- inla.read.graph(graph)
        cc.n <- sapply(g$cc$nodes, length)
        cc.n1 <- sum(cc.n == 1L)
        cc.n2 <- sum(cc.n >= 2L)

        if (is.null(rankdef)) {
            if (scale.model) {
                rankdef <- cc.n2
            } else {
                rankdef <- cc.n2 + cc.n1
            }
            if (!(empty.extraconstr(extraconstr))) {
                rankdef <- rankdef + dim(extraconstr$A)[1]
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
                    m <- inla.model.properties(model, "latent")
                    if (m$augmented) {
                        N <- m$aug.factor * n
                        offset <- (m$aug.constr - 1L) * n
                        stopifnot(length(offset) == 1)
                    } else {
                        N <- n
                        offset <- 0
                    }
                    constr <- FALSE
                    AA <- matrix(0, cc.n2, N)
                    ee <- rep(0, cc.n2)

                    if (debug) {
                        print(paste("add new extraconstr, dim = ", cc.n2, "x", n))
                    }

                    k <- 1L
                    for (i in 1L:length(cc.n)) {
                        if (cc.n[i] >= 2L) {
                            AA[k, offset + g$cc$nodes[[i]]] <- 1
                            k <- k + 1L
                        }
                    }
                    stopifnot(k - 1L == cc.n2)

                    if (!empty.extraconstr(extraconstr)) {
                        extraconstr$A <- rbind(AA, extraconstr$A)
                        extraconstr$e <- c(ee, extraconstr$e)
                    } else {
                        extraconstr <- list(A = AA, e = ee)
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
                    diagonal <- inla.set.f.default()$diagonal
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
            rgeneric <- rgeneric$f$rgeneric
        }
        R.init <- rgeneric$R.init
        stopifnot(inherits(rgeneric, "inla.rgeneric"))
        ## add an 'Id' so we know who we are
        rgeneric <- list(model = rgeneric, Id = vars[[1]], R.init = R.init)
    }

    if (model %in% "cgeneric") {
        if (inla.is.element("f", cgeneric) && inla.is.element("cgeneric", cgeneric$f)) {
            cgeneric <- cgeneric$f$cgeneric
        }
        stopifnot(inherits(cgeneric, "inla.cgeneric"))
        ## add an 'Id' so we know who we are
        cgeneric <- list(model = cgeneric, Id = vars[[1]])
    }

    ret <- list(
        Cmatrix = Cmatrix,
        Z = Z,
        bvalue = bvalue,
        cdf = cdf,
        compute = compute,
        constr = constr,
        control.group = cont.group,
        control.scopy = cont.scopy,
        cyclic = cyclic,
        d = d,
        diagonal = diagonal,
        extraconstr = inla.ifelse(empty.extraconstr(extraconstr), NULL, extraconstr),
        graph = graph,
        group = group,
        hyper = hyper,
        label = term,
        model = model,
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
        rankdef = rankdef,
        replicate = replicate,
        same.as = same.as,
        season.length = season.length,
        spde.prefix = spde.prefix,
        spde2.prefix = spde2.prefix,
        spde2.transform = spde2.transform,
        spde3.prefix = spde3.prefix,
        spde3.transform = spde3.transform,
        term = term,
        values = values,
        order = order,
        weights = weights,
        scale = scale,
        rgeneric = rgeneric,
        cgeneric = cgeneric,
        scale.model = as.logical(scale.model),
        adjust.for.con.comp = as.logical(adjust.for.con.comp),
        args.intslope = args.intslope,
        args.slm = args.slm,
        args.ar1c = args.ar1c,
        vb.correct = vb.correct,
        locations = locations
    )

    if (debug) print(ret)
    return(ret)
}

## "inla.model.class" is a generic model class that can be inherited
## from to mark a class as an inla model object class.
`inla.model.object.classes` <- function() {
    return(c(
        "inla.model.class", "inla.wrapper.model",
        "inla.spde", "inla.spde1", "inla.spde2", "inla.spde3",
        "inla.rgeneric", "inla.cgeneric"
    ))
}

## List spde classes that f() can write to disk when needed.
`inla.spde.object.classes` <- function() {
    return(c("inla.spde1", "inla.spde2", "inla.spde3"))
}
