## Export: inla

##! \name{inla}
##! \alias{inla}
##!
##! \title{Bayesian analysis of
##! structured additive models} \description{\code{inla} performs a
##! full Bayesian analysis of additive models using Integrated Nested
##! Laplace approximation
##! }
##!
##! \usage{
##!inla(
##!    formula,
##!    family = "gaussian", 
##!    contrasts = NULL,
##!    data,
##!    quantiles=c(0.025, 0.5, 0.975),
##!    E = NULL,
##!    offset=NULL,
##!    scale = NULL,
##!    weights = NULL,
##!    Ntrials = NULL,
##!    strata = NULL,
##!    link.covariates = NULL,
##!    verbose = FALSE,
##!    lincomb = NULL,
##!    control.compute = list(),
##!    control.predictor = list(),
##!    control.family = list(),
##!    control.inla = list(),
##!    control.results = list(),
##!    control.fixed = list(),
##!    control.mode = list(),
##!    control.expert = list(),
##!    control.hazard = list(),
##!    control.lincomb = list(),
##!    control.update = list(), 
##!    only.hyperparam = FALSE,
##!    inla.call = inla.getOption("inla.call"),
##!    inla.arg = inla.getOption("inla.arg"),
##!    num.threads = inla.getOption("num.threads"),
##!    blas.num.threads = inla.getOption("blas.num.threads"),
##!    keep = inla.getOption("keep"),
##!    working.directory = inla.getOption("working.directory"),
##!    silent = inla.getOption("silent"),
##!    debug = inla.getOption("debug"), 
##!    .parent.frame = parent.frame()
##!    )
##!
##! }
##! \arguments{
    
`inla` =
    function (##! \item{formula}{ A \code{inla} formula like \code{y
              ##!~1 + z + f(ind, model="iid")} + f(ind2,
              ##!weights, model="ar1") This is much like the formula
              ##!for a \code{glm} except that smooth or spatial terms
              ##!can be added to the right hand side of the formula.
              ##!See \code{\link{f}} for full details and the web site
              ##!\url{www.r-inla.org} for several worked out
              ##!examples. Each smooth or spatial term specified
              ##!through \code{f} should correspond to separate column
              ##!of the data frame \code{data}.
              ##!The response
              ##!variable, \code{y} can be a univariate response
              ##!variable, a list or the output of the function
              ##!\code{inla.surf} for survival analysis models.}
              formula,
              
              ##!\item{family}{ A string indicating the likelihood
              ##! family. The default is \code{gaussian} with identity
              ##! link. See \code{names(inla.models()$likelihood)} for a
              ##! list of possible alternatives and use \code{\link{inla.doc}}
              ##! for detailed docs for individual families.}
              family = "gaussian", 
              
              ##!\item{contrasts}{Optional contrasts for the fixed
              ##!effects; see \code{?lm} or \code{?glm} for details.}
              contrasts = NULL,
              
              ##!\item{data}{ A data frame or list containing the
              ##!variables in the model.  The data frame MUST be
              ##!provided}
              data,
              
              ##!\item{quantiles}{ A vector of quantiles,
              ##!\eqn{p(0), p(1),\dots}{p(0), p(1),\ldots} to compute
              ##!for each posterior marginal. The function returns,
              ##!for each posterior marginal, the values
              ##!\eqn{x(0), x(1),\dots}{x(0), x(1),\ldots} such that
              ##!\deqn{\mbox{Prob}(X<x(p))=p}{Prob(X<x)=p} }
              quantiles=c(0.025, 0.5, 0.975),
              
              ##!\item{E}{ Known component in the mean for the Poisson
              ##!likelihoods defined as \deqn{E_i\exp(\eta_i)}{E
              ##!exp(eta)} where \deqn{\eta_i}{eta} is the linear
              ##!predictor. If not provided it is set to \code{rep(1, n.data)}.}
              E = NULL,
              
              ##!\item{offset}{This argument is used to specify an
              ##!a-priori known and fixed component to be included in
              ##!the linear predictor during fitting.  This should be
              ##!\code{NULL} or a numeric vector of length either one
              ##!or equal to the number of cases. One or more
              ##!\code{offset()} terms can be included in the formula
              ##!instead or as well, and if both are used, they are
              ##!combined into a common offset.  If the
              ##!\code{A}-matrix is used in the linear predictor
              ##!statement \code{control.predictor}, then the
              ##!\code{offset} given in this argument is added to
              ##!\code{eta*}, the linear predictor related to the
              ##!observations, as \code{eta* = A eta + offset},
              ##!whereas an offset in the formula is added to
              ##!\code{eta}, the linear predictor related to the
              ##!formula, as \code{eta = ... + offset.formula}. So in
              ##!this case, the offset defined here and in the formula
              ##!has a different meaning and usage.}
              offset=NULL,
              
              ##!\item{scale}{ Fixed (optional) scale parameters of
              ##!the precision for Gaussian and Student-T response
              ##!models. Default value is rep(1, n.data).}
              scale = NULL,
              
              ##!\item{weights}{ Fixed (optional) weights parameters of
              ##!the likelihood, so the log-likelihood[i] is changed into
              ##!weights[i]*log-likelihood[i]. Default value is rep(1,
              ##!n.data). Due to the danger of mis-interpreting the results (see below), this option is DISABLED
              ##!by default. You can enable this option for the rest of your \code{R} session,
              ##!doing \code{inla.setOption(enable.inla.argument.weights=TRUE)}.
              ##!WARNING: The normalizing constant for the likelihood is NOT recomputed, so
              ##!ALL marginals (and the marginal likelihood) must be interpreted with great care.
              ##!Possibly, you may want to set the prior for the hyperparameters to \code{"uniform"}
              ##!and the integration strategy to \code{"eb"} to mimic a maximum-likelihood approach.}
              weights = NULL,
              
              ##!\item{Ntrials}{A vector containing the number of trials for the \code{binomial} 
              ##!likelihood and variantes, or the number of required successes for the
              ##!\code{nbinomial2} likelihood. Default value is \code{rep(1, n.data)}.}
              Ntrials = NULL,
              
              ##!\item{strata}{Fixed (optional) strata indicators 
              ##!for tstrata likelihood model.}
              strata = NULL,
              
              ##!\item{link.covariates}{A vector or matrix with covariates for link functions}
              link.covariates = NULL,

              ##!\item{verbose}{
              ##!Boolean indicating if the \code{inla}-program should
              ##!run in a verbose mode (default \code{FALSE}).}
              verbose = FALSE,
              
              ##!\item{lincomb}{ Used to define linear combination of
              ##!nodes in the latent field. The posterior distribution
              ##!of such linear combination is computed by the
              ##!\code{inla} function. See
              ##!\url{www.r-inla.org/faq} for examples of
              ##!how to define such linear combinations.}
              lincomb = NULL,
              
              ##!\item{control.compute}{ See \code{?control.compute}}
              control.compute = list(),
              
              ##!\item{control.predictor}{ See
              ##!\code{?control.predictor}}
              control.predictor = list(),
              
              ##!\item{control.family}{ See \code{?control.family}}
              control.family = list(),
              
              ##!\item{control.inla}{ See \code{?control.inla}}
              control.inla = list(),
              
              ##!\item{control.results}{ See \code{?control.result}}
              control.results = list(),
              
              ##!\item{control.fixed}{ See \code{?control.fixed}}
              control.fixed = list(),
              
              ##!\item{control.mode}{ See \code{?control.mode}}
              control.mode = list(),
              
              ##!\item{control.expert}{ See \code{?control.expert}}
              control.expert = list(),
              
              ##!\item{control.hazard}{ See \code{?control.hazard}}
              control.hazard = list(),
              
              ##!\item{control.lincomb}{ See \code{?control.lincomb}}
              control.lincomb = list(),
              
              ##!\item{control.update}{ See \code{?control.update}}
              control.update = list(),
              
              ##!\item{only.hyperparam}{ A boolean variable saying if
              ##!only the hyperparameters should be computed. This option is mainly used
              ##!internally. (TODO: This option should not be located here,  change it!)}
              only.hyperparam = FALSE,
              
              ##!\item{inla.call}{ The path to, or the name of, the
              ##!\code{inla}-program. This is program is installed
              ##!together with the \code{R}-package, but, for example,
              ##!a native compiled version can be used instead to
              ##!improve the performance.}
              inla.call = inla.getOption("inla.call"),
              
              ##!\item{inla.arg}{ A string indicating ALL arguments to
              ##!the 'inla' program and do not include default
              ##!arguments. (OOPS: This is an expert option!)}
              inla.arg = inla.getOption("inla.arg"),
              
              ##!\item{num.threads}{ Maximum number of threads the
              ##!\code{inla}-program will use}
              num.threads = inla.getOption("num.threads"),
              
              ##!\item{blas.num.threads}{The absolute value of \code{blas.num.threads} is the maximum
              ##!number of threads the the \code{openblas}/\code{mklblas} will use (if available). If
              ##!\code{blas.num.threads} > 0, then the environment variables
              ##!\code{OPENBLAS_NUM_THREADS} and \code{MKL_NUM_THREADS} will be assigned. If
              ##!\code{blas.num.threads} < 0, then the environment variables
              ##!\code{OPENBLAS_NUM_THREADS} and \code{MKL_NUM_THREADS} will be assigned unless they
              ##!are already defined. If \code{blas.num.threads} = 0, then variables
              ##!\code{OPENBLAS_NUM_THREADS} and \code{MKL_NUM_THREADS} will be removed.}
              blas.num.threads = inla.getOption("blas.num.threads"),
              
              ##!\item{keep}{ A boolean variable indicating that the
              ##!working files (ini file, data files and results
              ##!files) should be kept. If TRUE and no
              ##!\code{working.directory} is specified the working
              ##!files are stored in a directory called "inla".  }
              keep = inla.getOption("keep"),
              
              ##!\item{working.directory}{ A string giving the name
              ##!of an non-existing directory where to store the
              ##!working files.}
              working.directory = inla.getOption("working.directory"),
              
              ##!\item{silent}{If equal to 1L or TRUE, then the
              ##!\code{inla}-program would be ``silent''. If equal to
              ##!2L, then supress also error messages from the
              ##!\code{inla}-program.}
              silent = inla.getOption("silent"),
              
              ##!\item{debug}{ If \code{TRUE}, then enable some debug
              ##!output.  }
              debug = inla.getOption("debug"),

              ##!\item{.parent.frame}{Internal use only}
              .parent.frame = parent.frame()
              )
    ##!}

    ##!\value{%%
    
    ##!\code{inla} returns an object of class \code{"inla"}. This is a list
    ##!containing at least the following arguments:
    
    ##!\item{summary.fixed}{Matrix containing the mean and standard
    ##!deviation (plus, possibly quantiles and cdf) of the the fixed
    ##!effects of the model.}
    
    ##!\item{marginals.fixed}{
    ##!A list containing the posterior marginal
    ##!densities of the fixed effects of the model.}

    ##!\item{summary.random}{List of matrices containing the mean and
    ##!standard deviation (plus, possibly quantiles and cdf) of the
    ##!the smooth or spatial effects defined through \code{f()}.}
  
    ##!\item{marginals.random}{If
    ##!\code{return.marginals.random}=\code{TRUE} in
    ##!\code{control.results} (default), a list containing the
    ##!posterior marginal densities of the random effects defined
    ##!through \code{f}.}

    ##!\item{summary.hyperpar}{
    ##!A matrix containing the mean and sd
    ##!(plus, possibly quantiles and cdf) of the hyperparameters of
    ##!the model }

    ##!\item{marginals.hyperpar}{
    ##!A list containing the posterior marginal
    ##!densities of the hyperparameters of the model.} 
    

    ##!\item{summary.linear.predictor}{
    ##! A matrix containing the mean and sd
    ##! (plus, possibly quantiles and cdf) of the linear predictors
    ##! \eqn{\eta} in
    ##! the model
    ##! }

    ##!\item{marginals.linear.predictor}{
    ##! If \code{compute=TRUE} in
    ##! \code{control.predictor}, a list containing the posterior
    ##! marginals of the linear predictors \eqn{\eta} in the model.
    ##! }

    ##!\item{summary.fitted.values}{ A matrix containing the mean and
    ##! sd (plus, possibly quantiles and cdf) of the fitted values
    ##! \eqn{g^{-1}(\eta)} obtained by transforming the linear
    ##! predictors by the inverse of the link function. This quantity
    ##! is only computed if \code{marginals.fitted.values} is
    ##! computed. Note that if an observation is \code{NA} then the
    ##! identity link is used. You can manually transform a marginal
    ##! using \code{inla.marginal.transform()} or set the argument
    ##! \code{link} in the \code{control.predictor}-list;
    ##! see \code{?control.predictor}
    ##! }

    ##!\item{marginals.fitted.values}{
    ##! If \code{compute=TRUE} in
    ##! \code{control.predictor}, a list containing the posterior
    ##! marginals  of the fitted values
    ##! \eqn{g^{-1}(\eta)} obtained by
    ##! transforming the linear predictors by the inverse of the link
    ##! function.
    ##! Note that if an observation is \code{NA} then the
    ##! identity link is used. You can manually transform a marginal
    ##! using \code{inla.marginal.transform()} or set the argument
    ##! \code{link} in the \code{control.predictor}-list;
    ##! see \code{?control.predictor}
    ##! }
    
    ##!\item{summary.lincomb}{
    ##!If \code{lincomb != NULL} a list of
    ##!matrices containing the mean and sd (plus, possibly quantiles
    ##!and cdf) of all linear combinations defined.  }
    
    ##!\item{marginals.lincomb}{
    ##!If \code{lincomb != NULL} a list of
    ##! posterior marginals of all linear combinations defined.  } 
    

    ##!\item{joint.hyper}{
    ##!A matrix containing the joint density of
    ##!the hyperparameters (in the internal scale) }

    ##!\item{dic}{
    ##!If \code{dic}=\code{TRUE} in \code{control.compute}, the
    ##!deviance information criteria and effective number of parameters,
    ##!otherwise \code{NULL}
    ##!}

    ##!\item{cpo}{
    ##!If \code{cpo}=\code{TRUE} in \code{control.compute}, a list
    ##!of three elements: \code{cpo$cpo} are the values of the conditional 
    ##!predictive ordinate (CPO), \code{cpo$pit} are the values of the 
    ##!probability integral transform (PIT) and \code{cpo$failure} 
    ##!indicates whether some assumptions are violated. In short, if 
    ##!cpo$failure[i] > 0 then some assumption is violated, the higher the 
    ##!value (maximum 1) the more seriously.
    ##!}

    ##!\item{po}{
    ##!If \code{po}=\code{TRUE} in \code{control.compute}, a list
    ##!of one elements: \code{po$po} are the values of the 
    ##!predictive ordinate (CPO) (\code{pi(yi|y)})
    ##!}

    ##!\item{waic}{
    ##!If \code{waic}=\code{TRUE} in \code{control.compute}, a list
    ##!of two elements: \code{waic$waic} is the Watanabe-Akaike information criteria,  and 
    ##!\code{waic$p.eff} is the estimated effective number of parameters
    ##!}

    ##!\item{mlik}{
    ##!If \code{mlik}=\code{TRUE} in \code{control.compute}, the
    ##! log marginal likelihood of the model (using two different estimates), otherwise \code{NULL}
    ##!}
    
    ##! \item{neffp}{
    ##!Expected effective number of parameters in the model. The
    ##!standard deviation of the expected number of parameters and the
    ##!number of replicas for parameter are also returned}
    
    ##!\item{mode}{
    ##!A list of two elements: \code{mode$theta} is the
    ##!computed mode of the hyperparameters and \code{mode$x} is the
    ##!mode of the latent field given the modal value of the
    ##!hyperparamters.
    ##!}

    ##!\item{call}{
    ##!The matched call.}                   

    ##!\item{formula}{
    ##!The formula supplied}

    ##!\item{nhyper}{
    ##!The number of hyperparameters in the model}

    ##!\item{cpu.used}{
    ##!The cpu time used by the \code{inla} function}
    ##!}

    ##!\references{
    ##!Rue, H. and Martino, S. and Chopin, N. (2009)
    ##!\emph{Approximate Bayesian Inference for latent Gaussian models
    ##!using Integrated Nested Laplace Approximations, JRSS-series B
    ##!(with discussion)}, vol 71, no 2, pp 319-392.
    ##!Rue, H and Held, L. (2005) \emph{Gaussian Markov Random Fields
    ##!- Theory and Applications} Chapman and Hall}
    ##!\author{Havard Rue \email{hrue@r-inla.org} and Sara Martino}
    ##!\seealso{\code{\link{f}}, 
    ##!\code{\link{inla.hyperpar}} }
    ##!\examples{
    ##!\dontrun{
    ##!##See the web page \url{www.r-inla.org} for a series of worked out examples
    ##!}
    ##!}
{
    ## This will prevent values of 'OutDec' not '.' to cause error, as we create the Model.ini
    ## file with cat().
    old.options = options()
    on.exit(options(old.options))
    options(OutDec=".")
    
    my.time.used = numeric(4)
    my.time.used[1] = Sys.time()

    ## keep track of all priors and keep the original data
    all.hyper = list()
    data.orig = data
    
    if (missing(formula)) {
        stop("Usage: inla(formula, family, data, ...); see ?inla\n")
    }

    if (missing(data)) {
        stop("Missing data.frame/list `data'. Leaving `data' empty might lead to\n\t\tuncontrolled behaviour, therefore is it required.")
    }
    if (!is.data.frame(data) && !is.list(data)) {
        stop("\n\tArgument `data' must be a data.frame or a list.")
    }

    if (!missing(weights) && !is.null(weights)) {
        if (!inla.getOption("enable.inla.argument.weights")) {
            stop(paste("Argument 'weights' must be enabled before use due to the risk of mis-interpreting the results.\n",
                       "\tUse 'inla.setOption(\"enable.inla.argument.weights\", TRUE)' to enable it; see ?inla"))
        }
    }


    ## if data is a list, then it can contain elements that defines a
    ## model, like f(idx, model = model.objects). These objects crash
    ## the formula routines in R since then the data cannot be cast
    ## into a data.frame. to solve this, we remove such objects from
    ## data, and create a second data-object, data.model, which hold
    ## these.
    data.model = NULL
    if (is.list(data) && length(data) > 0L) {

        ## first check if all elements have names
        if (any(nchar(names(data)) == 0L)) {
            stop(paste("Elements in the 'data'-list has no name: data[[k]], for k=c(", inla.paste(which(nchar(names(data)) == 0L), sep=","), ")."))
        }
            
        i.remove = c()
        for(i in 1:length(data)) {
            ## these are the objects which we want to remove:
            ## inla.model.object.classes()
            if (any(inherits(data[[i]], inla.model.object.classes()))) {
                add.to = list(data[[i]])
                names(add.to) = names(data)[i]
                data.model = c(data.model,  add.to)
                i.remove = c(i.remove, i)
            }
        }
        names.remove = names(data)[i.remove]
        for(nm in names.remove) {
            idx = which(names(data) == nm)
            data[[idx]] = NULL
        }
    }

    ## check all control.xx arguments here. do the assign as variable
    ## expansion might occur.
    control.compute = inla.check.control(control.compute, data)
    control.predictor = inla.check.control(control.predictor, data)
    ## I need to check for NA's already here.
    if (!is.null(control.predictor$A)) {
        control.predictor$A = inla.as.sparse(control.predictor$A, na.rm=TRUE, zeros.rm=TRUE)
    }
    ## do not check control.family here, as we need to know n.family
    control.inla = inla.check.control(control.inla, data)
    control.results = inla.check.control(control.results, data)
    control.fixed = inla.check.control(control.fixed, data)
    control.mode = inla.check.control(control.mode, data)
    control.expert = inla.check.control(control.expert, data)
    control.hazard = inla.check.control(control.hazard, data)
    control.lincomb = inla.check.control(control.lincomb, data)
    control.update = inla.check.control(control.update, data)

    n.family = length(family)
    for(i in 1:n.family) {
        family[i] = inla.trim.family(family[i])
        ## check if this family is legal.
        if (is.null(inla.model.properties(family[i], "likelihood"))) {
            stop(paste("Unknown family: ", family[i], ". Do `names(inla.models()$likelihood)' to list available families.", sep=""))
        }
    }

    ## if the user specify inla.call="remote" or "inla.remote" then
    ## use the internal inlaprogram
    remote = FALSE
    submit = FALSE
    ownfun = FALSE
    submit.id = ""

    if (is.function(inla.call)) {
        ## in this case, the responsibility is with the user
        ownfun = TRUE
    } else if (inla.strcasecmp(inla.call, "remote") ||
        inla.strcasecmp(inla.call, "inla.remote") ||
        length(grep("/inla.remote$", inla.call)) > 0 ||
        length(grep("/inla.remote.cygwin$", inla.call)) > 0) {
        remote = TRUE
        inla.call = system.file("bin/remote/inla.remote", package="INLA")
        if (inla.os("windows")) {
            inla.call = paste(inla.call, ".cygwin", sep="")
        }
    } else if (inla.strcasecmp(inla.call, "submit") ||
               inla.strcasecmp(inla.call, "inla.submit") ||
               length(grep("/inla.submit$", inla.call)) > 0) {
        remote = TRUE
        submit = TRUE
        submit.id = paste(gsub("[ :]", "-", date()), "---", as.integer(runif(1,min=1E8,max=1E9-1)), sep="")
        inla.call = system.file("bin/remote/inla.submit", package="INLA")
        if (inla.os("windows")) {
            inla.call = paste(inla.call, ".cygwin", sep="")
        } 
    }
    
    ## Need to do this here.
    cont.fixed = inla.set.control.fixed.default()
    cont.fixed[names(control.fixed)] = control.fixed

    ##
    ## check for survival model with a baseline-hazard. if so, then
    ## expand the data-frame and call inla() again.
    ##
    have.surv = FALSE
    cont.hazard = NULL
    for(i in 1:n.family) {
        have.surv = have.surv || inla.model.properties(family[i], "likelihood")$survival
    }

    if (have.surv && (inla.one.of(family, c("coxph")))) {
        ## This is not supported yet. 
        stopifnot(is.null(control.predictor$A))

        cph = inla.coxph(formula, data, control.hazard, debug = debug)
        result = inla(
            cph$formula,
            family = cph$family,
            data = c(as.list(cph$data), cph$data.list), 
            contrasts = contrasts, 
            quantiles=quantiles,
            E = cph$E,
            offset= offset, 
            scale= scale, 
            weights= weights, 
            Ntrials = NULL,             # Not used for the poisson
            strata = NULL,              # Not used for the poisson
            lincomb = lincomb,
            verbose = verbose,
            control.compute = control.compute,
            control.predictor = control.predictor,
            control.family = control.family,
            control.inla = control.inla,
            control.results = control.results,
            control.fixed = control.fixed,
            control.mode = control.mode,
            control.expert = control.expert,
            control.hazard = control.hazard,
            control.lincomb = control.lincomb,
            control.update = control.update,
            only.hyperparam = only.hyperparam,
            inla.call = inla.call,
            inla.arg = inla.arg,
            num.threads = num.threads,
            blas.num.threads = blas.num.threads,
            keep = keep,
            working.directory = working.directory,
            silent = silent,
            debug = debug)

        ## replace the argument so it can be reused, if...
        result$call.orig = deparse(match.call())
        return (result)
    }
    
    ## this is nice hack ;-) we keep the original response. then we
    ## delete it from 'data' keeping a copy of the original one
    y...orig = eval(parse(text=formula[2]), data)

    if (n.family > 1) {
        y...orig = inla.as.list.of.lists(y...orig)
        ny = max(sapply(y...orig, function(xx) if (is.list(xx)) max(sapply(xx, length)) else length(xx)))
        nc = length(y...orig)
        if (n.family != nc) {
            stop(paste("Number of families", n.family,
                       "does not match number of response variables", nc))
        }
    } else {
        nc = NULL ## not in use
        if (inherits(y...orig, "inla.surv")) {
            class(y...orig) = NULL
            ny = max(sapply(y...orig, length))
        } else if (inherits(y...orig, "inla.mdata")) {
            class(y...orig) = NULL
            ny = max(sapply(y...orig, length))
        } else {
            if (length(dim(y...orig)) == 2) {
                ## some matrix type, could be a sparse matrix
                stopifnot(dim(y...orig)[2] == 1)
                y...orig = c(drop(as.matrix(y...orig)))
            }
            ny = length(y...orig)
        }
    }

    data = inla.remove(as.character(formula[2]), data) 
    if (!is.null(control.predictor$A)) {
        MN = inla.sparse.dim(control.predictor$A)
        MPredictor = MN[1]
        NPredictor = MN[2]
        stopifnot(MPredictor == ny)
    } else {
        MPredictor = 0L
        NPredictor = ny
    }
    NData = inla.ifelse(MPredictor > 0,  MPredictor, NPredictor)
    stopifnot(NData > 0)
    stopifnot(NPredictor > 0)
    if (debug) {
        print(paste("MPredictor", MPredictor))
        print(paste("NPredictor", NPredictor))
        print(paste("NData", NPredictor))
    }

    if (n.family == 1) {
        y...fake = c(rep(Inf, NData))
        if (debug)
            print(paste("y...fake has length", length(y...fake)))
        
    } else {
        y...fake = c(rep(Inf, NData))
        if (debug)
            print(paste("y...fake has length", length(y...fake)))
    }

    if (FALSE) {
        if (n.family != 1)
            stop(paste("length(family) are", n.family, "but the response has only one column!"))
        
        y...fake = c(y...orig)
        y...fake[is.na(y...fake)] = Inf ## otherwise model.matrix() fails below
    }
    if (debug) {
        cat("n.family", n.family, "\n")
    }

    ## ...and then we add the fake to the data-frame after removing the original response
    if (is.data.frame(data)) {
        if (MPredictor > 0) {
            if (MPredictor != NPredictor) {
                stop(paste("It can be ``dangerous'' to use a 'data.frame' as data, when the\n", 
                           "\tA-matrix is not a square matrix; please pass the data\n", 
                           "\tas a list: data = list(...)."))
            }
        }
        data = as.data.frame(c(as.data.frame(data), list(y...fake=y...fake)))
    } else {
        data = c(data, list(y...fake=y...fake))
        if (MPredictor > 0 && (MPredictor > NPredictor)) {
            data$y...fake = data$y...fake[1:NPredictor]
        }
        if (MPredictor > 0 && (MPredictor < NPredictor)) {
            data$y...fake = c(data$y...fake, rep(Inf, NPredictor - MPredictor))
        }
    }

    data.same.len = inla.fix.data(data, NPredictor)
    if (debug) {
        print(c("Entries with same length:", names(data.same.len)))
    }

    ## creates a new version of ``formula'' and keep the old one for reference
    formula.orig = formula
    ##inla.eval(paste("formula = y...fake ~ ", inla.formula2character(formula.orig[3])))
    formula = update.formula(formula, y...fake ~ .)
    ## parse the formula
    gp = inla.interpret.formula(formula, data.same.len=data.same.len, data=data,
        data.model = data.model, parent.frame = .parent.frame)
    call = deparse(match.call())

    ## issue a warning if the intercept is spesified while the
    ## control.predictor A matrix is used.
    if (length(grep("\\+ ?1($| )", gp$fixf)) && !is.null(control.predictor$A)) {
        warning("The A-matrix in the predictor (see ?control.predictor) is specified
  but an intercept is in the formula. This will likely result
  in the intercept being applied multiple times in the model, and is likely
  not what you want. See ?inla.stack for more information.
  You can remove the intercept adding ``-1'' to the formula.")
    }

    if (gp$n.fix > 0) {
        ## use y...fake also here (which is the same as in gp$fixf[2])
        ## and construct the model.matrix().
        new.fix.formula = gp$fixf
        ##inla.eval(paste("new.fix.formula = y...fake ~ ", inla.formula2character(gp$fixf[3])))
        new.fix.formula = update.formula(new.fix.formula, y...fake ~ .) 

        ## replace NA's in covariates with 0. Ignore factors. NA's in
        ## factors is done afterwards.
        inla.na.action = function(x, ...) {
            for(k in seq_along(x)) {
                if ((is.numeric(x[[k]]) || inla.is.matrix(x[[k]])) && !is.factor(x[[k]])) {
                    x[[k]][is.na(x[[k]])] = 0
                }
            }
            return (na.pass(x))
        }
        if (inla.one.of(cont.fixed$expand.factor.strategy, "inla")) {
            ## expand all factors as matrices, exclude possible
            ## factors in the f()'s. it is ok to include all terms in
            ## the f()'s as duplicated names are not allowed.
            exclude.names = c()
            for (k in seq_along(gp$random.spec)) {
                exclude.names = c(exclude.names, gp$random.spec[[k]]$label)
            }
            data.same.len = inla.expand.factors(data.same.len, exclude.names)
        } else if (inla.one.of(cont.fixed$expand.factor.strategy, "model.matrix")) {
            ## do nothing
        } else {
            stop(paste("Unknown value for flag 'expand.factor.strategy' in 'control.fixed':",
                       cont.fixed$expand.factor.strategy))
        }
        if (inla.require("MatrixModels")) {
            gp$model.matrix = MatrixModels::model.Matrix(
                new.fix.formula,
                data = model.frame(new.fix.formula, data.same.len, na.action=inla.na.action),
                contrasts.arg = contrasts, sparse=TRUE)
        } else {
            gp$model.matrix = model.matrix(
                new.fix.formula,
                data = model.frame(new.fix.formula, data.same.len, na.action=inla.na.action),
                contrasts.arg = contrasts)
        }
        ## as NA's in factors are not set to zero in
        ## 'inla.na.action'. Do that here if the strategy is 'inla',
        ## otherwise signal an error.
        if (any(is.na(gp$model.matrix))) {
            if (inla.one.of(cont.fixed$expand.factor.strategy, "inla")) {
                gp$model.matrix[is.na(gp$model.matrix)] = 0
            } else {
                stop(paste("With control.fixed = list(expand.factor.strategy='model.matrix'),", 
                           "then NA's in factor are not allowd. Please use strategy 'inla' instead."))
            }
        }
        ## this have to match
        stopifnot(dim(gp$model.matrix)[1L] == NPredictor)

        if (any(is.infinite(gp$model.matrix))) {
            n.infs = sum(is.infinite(gp$model.matrix))
            stop(paste("There are", n.infs, "Inf's in the model.matrix. This is not allowed."))
        }

        ## n.fix can have been changed here due to a `-1'
        gp$n.fix = dim(gp$model.matrix)[2L]

        ## if reqested, then setup the fixed effects as linear
        ## combinations, so we can compute their posterior correlation
        ## matrix. so the lincombs are 'inla.make.lincomb(z=1)'
        ## etc. we need the quotes \"..\" for the "(Intercept)".
        if (!is.null(control.fixed$correlation.matrix) &&
            control.fixed$correlation.matrix && !is.null(gp$model.matrix)) {
            fix.names = colnames(gp$model.matrix)
            lc.all.fix = c()
            for(fix.name in fix.names) {
                lc.fix = inla.eval(paste("inla.make.lincomb(\"", fix.name,"\"=1)", sep=""))
                names(lc.fix) = fix.name
                lc.all.fix = c(lc.all.fix, lc.fix)
            }
            ## if we have defined lincomb's also in the inla call,
            ## then simply append the fixed-effects at the end,
            ## otherwise set the 'lincomb' argument. however. for this
            ## to work with repeated calls (like inla.hyperpar()) we
            ## have to turn off control.fixed$correlation.matrix as we
            ## have moved that part into 'lincomb'. also, we need to
            ## make sure that we actually compute the correlation
            ## matrix for the linear combinations in control.inla$...
            lincomb = c(lincomb, lc.all.fix)
            control.fixed$correlation.matrix = FALSE
            control.inla$lincomb.derived.correlation.matrix = TRUE
        }
        
    } else {
        gp$model.matrix = NULL
    }
    
    ## control what should be computed
    cont.compute = inla.set.control.compute.default()
    cont.compute[names(control.compute)] = control.compute
    if (only.hyperparam) {
        cont.compute$hyperpar = TRUE
        cont.compute$dic = cont.compute$cpo = cont.compute$po = cont.compute$waic = FALSE 
    } 
    
    ## control predictor section
    cont.predictor = inla.set.control.predictor.default()
    cont.predictor[names(control.predictor)] = control.predictor
    cont.predictor$hyper = inla.set.hyper("predictor", "predictor",
        cont.predictor$hyper, cont.predictor$initial,
        cont.predictor$fixed, cont.predictor$prior, cont.predictor$param)
    all.hyper$predictor$hyper = cont.predictor$hyper
    if (cont.compute$cpo || cont.compute$dic || cont.compute$po || cont.compute$waic || !is.null(cont.predictor$link))
        cont.predictor$compute=TRUE
    if (only.hyperparam) {
        cont.predictor$compute = cont.predictor$return.marginals = FALSE
        cont.predictor$cdf = cont.predictor$quantiles = NULL
    }
    
    ## control inla
    cont.inla =inla.set.control.inla.default(family=family)
    cont.inla[names(control.inla)] = control.inla

    ## control.family
    control.family.orig = control.family
    if (n.family == 1) {
        if (!missing(control.family) && (inla.is.list.of.lists(control.family) && length(control.family) > 1L))
            stop(paste("Argument 'control.family' does not match length(family)=", n.family))
    } else {
        if (!missing(control.family) && !(inla.is.list.of.lists(control.family) && length(control.family) == n.family))
            stop(paste("Argument 'control.family' does not match length(family)=", n.family))
    }
    
    if (missing(control.family)) {
        tt = list(list())
        for(i in 1:n.family)
            tt[[i]] = control.family
        control.family = tt
    } else if (n.family == 1) {
        control.family = list(control.family)
    }
    
    if (length(control.family) == 1 && n.family > 1)
        stop(paste("length(control.family) = 1 while length(family) > 1."))

    ## finally, do the check. Note that the name of the argument is
    ## also used in this function, so we need to borrow the name.
    control.family.save = control.family
    for(ii in 1:n.family) {
        control.family = control.family.save[[ii]]
        ## need to be able to say when n.family =  1
        ## control.family = list(
        ##    list(hyper = list(), control.link = list())))
        while (is.list(control.family) &&
               length(control.family) > 0 &&
               is.null(names(control.family))) {
            control.family = control.family[[1]]
        }
        control.family.save[[ii]] = inla.check.control(control.family, data)
    }
    control.family = control.family.save
    cont.family = list(list())
    for(i.family in 1:n.family) {
        cont.family[[i.family]] = inla.set.control.family.default()
        cont.family[[i.family]]$control.mix = inla.set.control.mix.default()
        cont.family[[i.family]]$control.link = inla.set.control.link.default()

        ## need to take option 'control.mix' and 'control.link' out and process it seperately
        c.mix = control.family[[i.family]]$control.mix
        c.link = control.family[[i.family]]$control.link
        control.family[[i.family]]$control.mix = NULL
        control.family[[i.family]]$control.link = NULL

        cont.family[[i.family]][names(control.family[[i.family]])] = control.family[[i.family]]
        cont.family[[i.family]]$hyper = inla.set.hyper(
                                   family[i.family],
                                   "likelihood",
                                   cont.family[[i.family]]$hyper, 
                                   cont.family[[i.family]]$initial, 
                                   cont.family[[i.family]]$fixed,
                                   cont.family[[i.family]]$prior,
                                   cont.family[[i.family]]$param)
        all.hyper$family[[i.family]] = list(
                            hyperid = paste("INLA.Data", i.family, sep=""), 
                            label = family[i.family],
                            hyper = cont.family[[i.family]]$hyper)
        
        cont.family[[i.family]]$control.mix[names(c.mix)] = c.mix
        cont.family[[i.family]]$control.link[names(c.link)] = c.link
        if (!is.null(cont.family[[i.family]]$control.mix$model)) {
            cont.family[[i.family]]$control.mix$hyper = inla.set.hyper(
                                       cont.family[[i.family]]$control.mix$model,
                                       "mix",
                                       cont.family[[i.family]]$control.mix$hyper, 
                                       cont.family[[i.family]]$control.mix$initial, 
                                       cont.family[[i.family]]$control.mix$fixed,
                                       cont.family[[i.family]]$control.mix$prior,
                                       cont.family[[i.family]]$control.mix$param)
            all.hyper$family[[i.family]]$mix$hyper= cont.family[[i.family]]$control.mix$hyper
        }
        cont.family[[i.family]]$control.link$hyper = inla.set.hyper(
                                   cont.family[[i.family]]$control.link$model,
                                   "link",
                                   cont.family[[i.family]]$control.link$hyper, 
                                   cont.family[[i.family]]$control.link$initial, 
                                   cont.family[[i.family]]$control.link$fixed,
                                   cont.family[[i.family]]$control.link$prior,
                                   cont.family[[i.family]]$control.link$param)
        all.hyper$family[[i.family]]$link$hyper = cont.family[[i.family]]$control.link$hyper
    }
    
    ## control results
    cont.results = inla.set.control.results.default()
    cont.results[names(control.results)] = control.results
    
    ## Create the directory where to store Model.ini and data.files
    ## and results.file
    if (!is.null(working.directory))
        keep=TRUE
    if (keep) {
        ##create the directory locally or whereever specified by the user
        if (is.null(working.directory)) {
            working.directory="inla.model"
            working.directory.start="inla.model"
        } else {
            working.directory.start= working.directory
        }

        ## if this directory already exists, try the numbered versions
        kk=1
        ans=file.exists(working.directory)
        while(ans) {
            working.directory=paste(working.directory.start,"-", kk, sep="")
            kk=kk+1
            ans=file.exists(working.directory)
        }
        inla.dir=working.directory
        tmp = inla.dir.create(inla.dir, StopOnError = FALSE) ## return NULL if fail
        if (is.null(tmp)) {
            ## this is the failsafe-option
            inla.dir = paste(inla.dir, "-", substring(as.character(runif(1)), 3), sep="")
            tmp = inla.dir.create(inla.dir, StopOnError = FALSE)
            if (is.null(tmp)) {
                stop(paste("Cannot create directory [", inla.dir,
                           "] even after trying a random dirname. I give up.", sep=""))
            }
        }
        if (verbose) {
            cat("Model and results are stored in working directory [", inla.dir,"]\n", sep="")
        }
    } else {
        ##create a temporary directory
        inla.dir=inla.tempfile()
        inla.dir=gsub("\\\\", "/", inla.dir)
        inla.dir.create(inla.dir)
    }
    ## Create a directory where to store data and results
    inla.dir = normalizePath(inla.dir)
    data.dir=paste(inla.dir, "/data.files", sep="")
    results.dir = paste(inla.dir, "/results.files", sep="")
    inla.dir.create(data.dir)

    ## create the .file.ini and make the problem.section
    file.ini = paste(inla.dir, "/Model.ini", sep="")
    file.log = paste(inla.dir, "/Logfile.txt", sep="")
    file.log2 = paste(inla.dir, "/Logfile2.txt", sep="")

    ## problem section
    if (debug) 
        print("prepare problem section")

    inla.problem.section(file = file.ini, data.dir = data.dir, result.dir = results.dir,
                         hyperpar = cont.compute$hyperpar, return.marginals = cont.compute$return.marginals,
                         dic = cont.compute$dic, mlik = cont.compute$mlik,
                         cpo = cont.compute$cpo,
                         ## these two are merged together as they are compute together
                         po = (cont.compute$po || cont.compute$waic), 
                         quantiles = quantiles, smtp = cont.compute$smtp, q = cont.compute$q,
                         openmp.strategy = cont.compute$openmp.strategy, graph = cont.compute$graph,
                         config = cont.compute$config, gdensity = cont.compute$gdensity)

    ## PREPARE RESPONSE AND FIXED EFFECTS
    if (debug)
        cat("Prepare inla file.....")

    ## copy the argument-lists
    mf = match.call(expand.dots = FALSE)
    mf$family = NULL; mf$quantiles=NULL; 
    mf$verbose = NULL; mf$control.compute = NULL; mf$control.predictor = NULL;
    mf$silent = NULL; mf$control.hazard=NULL;
    mf$control.family = NULL;  mf$control.update = NULL;
    mf$control.inla = NULL; mf$control.results = NULL; mf$control.fixed = NULL; mf$control.lincomb=NULL;
    mf$control.mode = NULL; mf$control.expert = NULL; mf$inla.call = NULL;
    mf$num.threads = NULL; mf$blas.num.threads = NULL; mf$keep = NULL;
    mf$working.directory = NULL; mf$only.hyperparam = NULL; mf$debug = NULL; mf$contrasts = NULL; 
    mf$inla.arg = NULL; mf$lincomb=NULL; mf$.parent.frame = NULL;
    mf$data = data.same.len

    if (gp$n.fix > 0)
        mf$formula = gp$fixf
    else if (gp$n.random > 0)
        mf$formula = gp$randf
    else
        mf$formula = y...fake ~ 1
    
    ## these we need
    mf$na.action = na.pass
    mf[[1]] = as.name("model.frame")

    if (gp$n.random > 0) {
        rf = mf ## for later use
        rf$weights = rf$scale = rf$Ntrials = rf$offset = rf$E =  rf$strata = rf$link.covariates = NULL ## these we do not need
        rf$formula = gp$randf
        rf$data = data.same.len
        rf = eval.parent(rf)
    } else {
        rf = NULL
    }
    
    if (gp$n.weights > 0) {
        wf = mf
        wf$weights = wf$scale = wf$Ntrials = wf$offset = wf$E =  wf$strata = wf$link.covariates = NULL ## these we do not need
        wf$formula = gp$weightf
        wf$data = data.same.len
        wf = eval.parent(wf)
    } else {
        wf = NULL
    }

    ## We set these here, instead of::
    ## scale = model.extract(mf, "scale")
    ## Ntrials = model.extract(mf, "Ntrials")
    ## E = model.extract(mf, "E")
    ## offset = as.vector(model.extract(mf, "offset"))

    for (nm in c("scale", "weights", "Ntrials", "offset", "E", "strata", "link.covariates")) {
        inla.eval(paste("tmp = try(eval(mf$", nm, ", data, enclos = parent.frame()), silent=TRUE)", sep=""))
        if (!is.null(tmp) && !inherits(tmp, "try-error")) {
            inla.eval(paste("mf$", nm, " = NULL", sep=""))
            inla.eval(paste(nm, " = tmp"))
        } else {
            inla.eval(paste("tmp = try(eval.parent(mf$", nm, "), silent=TRUE)", sep=""))
            if (!is.null(tmp) && !inherits(tmp, "try-error")) {
                inla.eval(paste("mf$", nm, " = NULL", sep=""))
                inla.eval(paste(nm, " = tmp"))
            } else {
                ## this *is* defined by default,  as all variables are default NULL
                inla.eval(paste("mf$", nm, " = NULL", sep=""))
                inla.eval(paste(nm, "= inla.eval(nm)", sep=""))
            }
        }
    }
    
    ## as there are functions as well with this name....
    if (is.function(offset))
        offset = NULL
    if (is.function(scale))
        scale = NULL

    ## ## ## ## ## ## ## ## ## ## ## ##
    ## ## ## ## ## ## ## ## ## ## ## ##

    mf = eval.parent(mf)
    indN = seq(0L, NPredictor-1L)
    indM = seq(0L, MPredictor-1L)
    indD = seq(0L, NData-1)
    
    ## this takes care of the offset: `offset' is the argument,
    ## `offset.formula' is in the formula. Note that 'offset' goes
    ## into eta* = A %*% eta + offset, whereas, eta = .... +
    ## offset.formula
    if (!is.null(gp$offset)) {
        ## there can be more offsets
        offset.formula = 0
        for(i in 1:length(gp$offset)) {
            offset.formula = offset.formula + as.vector(eval(parse(text=gp$offset[i]), data))
        }
    } else {
        offset.formula = rep(0, NPredictor)
    }

    offset.len = inla.ifelse(MPredictor > 0, MPredictor, NPredictor)
    offset.formula.len = NPredictor
    if (is.null(offset)) {
        offset = 0
    }
    if (length(offset) == 1L) {
        offset = rep(offset,  offset.len)
    }
    if (length(offset) != offset.len) {
        stop(paste("Length of argument 'offset' is wrong:", length(offset), "!=", offset.len))
    } 
    if (length(offset.formula) != offset.formula.len) {
        stop(paste("Length of 'offset(...)' in the formula is wrong:", length(offset.formula), "!=", offset.formula.len))
    } 

    if (length(family) == 1)
        family = rep(family, n.family)
    
    ## Create a file with the response, each for each 'family'
    for(i.family in 1:n.family) {
        if (n.family == 1)
            yy = y...orig
        else
            yy = y...orig[[i.family]]
        
        if (MPredictor > 0) {
            if (is.list(yy)) {
                stopifnot(max(sapply(yy, length)) == MPredictor)
            } else if (inla.is.matrix(yy)) {
                stopifnot(dim(yy)[1L] == MPredictor)
            } else {
                stopifnot(length(yy) == MPredictor)
            }
        }

        if (!is.null(yy) && !(is.numeric(yy) || is.list(yy) || is.matrix(yy) || inla.is.matrix(yy) || all(is.na(yy)))) {
            stop(paste("The response for family[", i.family, "] is not of type 'numeric|list|matrix'; don't know what to do.", sep=""))
        }

        ## we can check for 'survial' here, either if data is and family isn't,  or the oposite
        if (is.list(yy) &&
            all(is.element(names(yy), names(inla.surv(1, 1))))) {
            ## check if family is of survival type. if not, check if appending 'surv' is.
            if (!inla.model.properties(family[i.family], "likelihood")$survival) {
                new.fam = paste0(family[i.family], "surv")
                ok = inla.model.properties(new.fam, "likelihood", stop.on.error=FALSE)$survival
                if (!is.null(ok) && ok) {
                    warning(paste0("*** WARNING *** Input family is '", family[i.family],
                                   "' but input data is of 'inla.surv(...)' type\n",
                                   "  *** WARNING *** Do you ment to use family '", new.fam,
                                   "' ?"))
                }
            }
        }
        if (inla.model.properties(family[i.family], "likelihood")$survival && !is.list(yy)) {
            ## check if removing 'surv' makes a valid family
            new.fam = gsub("surv$", "", family[i.family])
            ok = inla.model.properties(new.fam, "likelihood",
                                       stop.on.error=FALSE)$survival
            if (!is.null(ok) && !ok) {
                warning(paste0("*** WARNING *** Input family is '", family[i.family],
                               "' and require input of type 'inla.surv(...)'\n",
                               "  *** WARNING *** Do you ment to use family '", new.fam,
                               "' ?"))
            }
        }

        files = inla.create.data.file(y.orig= yy, mf=mf, E=E, scale=scale, 
            weights=weights, Ntrials=Ntrials, strata=strata, 
            family=family[i.family], data.dir=data.dir, file=file.ini, debug=debug)
        
        ## add a section to the file.ini
        prop = inla.model.properties(family[i.family], "likelihood", stop.on.error=TRUE)

        if (debug) 
            print("prepare data section")

        ##....then create the new section 
        inla.family.section(file=file.ini, family=family[[i.family]], file.data=files$file.data, file.weights=files$file.weights,
                            control=cont.family[[i.family]], i.family=i.family, link.covariates = link.covariates, data.dir = data.dir)
    }

    ##create the PREDICTOR section. 
    if (debug) 
        print("prepare predictor section")

    stopifnot(!is.null(offset.formula) && !is.null(offset)) ## must be zeros if not used. this makes it easier
    offset.formula[is.na(offset.formula)] = 0
    offset[is.na(offset)] = 0
    if (!is.null(control.predictor$A)) {
        off = cbind(c(indM, MPredictor + indN), c(as.vector(control.predictor$A %*% offset.formula + offset), offset.formula))
    } else {
        off = cbind(indN, offset + offset.formula)
    }
    file.offset = inla.tempfile(tmpdir=data.dir)
    if (inla.getOption("internal.binary.mode")) {
        inla.write.fmesher.file(as.matrix(off), filename=file.offset, debug = debug)
    } else {
        file.create(file.offset)
        write(t(off), ncolumns=2L, file=file.offset, append=FALSE)
    }
    file.offset = gsub(data.dir, "$inladatadir", file.offset, fixed=TRUE)

    if (!is.null(cont.predictor$link)) {
        not.na = which(!is.na(cont.predictor$link))
        if (any(cont.predictor$link[not.na] != as.integer(cont.predictor$link[not.na])) || !all(cont.predictor$link[not.na] %in% 1L:n.family))
            stop(paste("Argument 'control.predictor$link' must consist of integers from the set 1:", n.family, " or NA's.", sep=""))
        if (!all(is.na(cont.predictor$link))) {
            if (max(cont.predictor$link[not.na]) > n.family) {
                stop(paste("n.family =", n.family,  "while argument 'control.predictor$link' has max=", max(cont.predictor$link[not.na])))
            }
        }
        if (!is.null(control.predictor$A)) {
            if (length(cont.predictor$link) == 1L) {
                ## shortcut
                cont.predictor$link = rep(cont.predictor$link, MPredictor)
            }
            if (!(length(cont.predictor$link) == MPredictor || length(cont.predictor$link) == MPredictor + NPredictor)) {
                stop(paste("Length of argument 'control.predictor$link' is wrong: length(link) = ",
                           length(cont.predictor$link), ". Length must be equal to ",
                           "length(A%*%eta)=", MPredictor, " or length(c(A%*%eta,eta))=", MPredictor + NPredictor,
                           ".", sep=""))
            }
            if (length(cont.predictor$link) == MPredictor) {
                tlink = cbind(0L:(MPredictor + NPredictor -1L), c(cont.predictor$link -1L, rep(NA, NPredictor)))
            } else {
                tlink = cbind(0L:(MPredictor + NPredictor -1L), cont.predictor$link - 1L)
            }
        } else {
            if (length(cont.predictor$link) == 1L) {
                ## shortcut
                cont.predictor$link = rep(cont.predictor$link, NPredictor)
            }
            if (!(length(cont.predictor$link) == NPredictor)) {
                stop(paste("Length of argument 'control.predictor$link' is wrong: length(link)=", length(control.predictor$link),
                           ". Length must be equal to ",
                           "length(eta)=", NPredictor, ".", sep=""))
            }
            stopifnot(length(cont.predictor$link) == NPredictor)
            tlink = cbind(indN, as.numeric(cont.predictor$link) -1L)
        }
        file.link.fitted.values = inla.tempfile(tmpdir=data.dir)
        if (inla.getOption("internal.binary.mode")) {
            inla.write.fmesher.file(as.matrix(tlink), filename=file.link.fitted.values, debug = debug)
        } else {
            file.create(file.link.fitted.values)
            write(t(tlink), ncolumns=2L, file=file.offset, append=FALSE)
        }
        file.link.fitted.values = gsub(data.dir, "$inladatadir", file.link.fitted.values, fixed=TRUE)
    } else {
        file.link.fitted.values = NULL
    }

    inla.predictor.section(file=file.ini, n=NPredictor, m=MPredictor,
                           predictor.spec=cont.predictor, file.offset=file.offset, data.dir = data.dir,
                           file.link.fitted.values = file.link.fitted.values)

    all.labels = character(0)
    if (!is.null(cont.predictor$compute) && cont.predictor$compute)
        all.labels = c(all.labels,"predictor")

    ## FIXED EFFECTS
    if (gp$n.fix >0) {
        nc = ncol(gp$model.matrix)
        labels = colnames(gp$model.matrix)
        for(i in 1:nc) {
            if (debug)
                cat("write label[", labels[i],"]\n")

            fixed.eff=cbind(indN, as.numeric(gp$model.matrix[, i]))
            ##remove lines with NA
            fixed.eff = fixed.eff[!is.na(fixed.eff[, 2]),, drop=FALSE]

            file.fixed=inla.tempfile(tmpdir=data.dir)
            if (inla.getOption("internal.binary.mode")) {
                inla.write.fmesher.file(as.matrix(fixed.eff), filename = file.fixed, debug=debug)
            } else {
                file.create(file.fixed)
                write(t(fixed.eff), ncolumns=ncol(fixed.eff), file=file.fixed, append=FALSE)
            }
            file.fixed = gsub(data.dir, "$inladatadir", file.fixed, fixed=TRUE)

            if (debug)
                cat("file fixed", file.fixed,"\n")

            all.hyper$fixed[[i]] =
                inla.linear.section(file=file.ini, file.fixed=file.fixed, label=labels[i],
                                    results.dir=paste("fixed.effect", inla.num(i), sep=""),
                                    control.fixed=cont.fixed, only.hyperparam=only.hyperparam)
        }
    }

    ##RANDOM EFFECT OR LINEAR EFFCT DEFINED VIA f()
    nr=gp$n.random
    n.weights=0
    j=0
    extra.fixed=0

    if (nr>0) {
        name.random.dir=c()
        if (nr!=(ncol(rf)-1))
            stop("\n\tSOMETHING STRANGE: You probably used one variable as a covariate more than once.")
        
        ## make the covariate numeric if its all NA (for which it is of type 'logical')
        stopifnot(ncol(rf) >= 2L)
        for(r in 1:(ncol(rf)-1)) { ## YES, '-1' is correct
            idx = rf[, r+1]
            if (is.logical(idx) && all(is.na(idx))) {
                rf[, r+1] = idx = as.numeric(idx)
            }
            if (all(is.na(idx))) {
                ## check that either 'n' or 'values' are given, to prevent errors from happening
                ## below.
                if (is.null(gp$random.spec[[r]]$values) && is.null(gp$random.spec[[r]]$n)) {
                    stop(paste("Model component 'f(", 
                               gp$random.spec[[r]]$term,
                               ", ...)' have only NA values in '",
                               gp$random.spec[[r]]$term, "'. ", 
                               "In this case, ",
                               "either argument 'n' or 'values' must be spesified.", sep=""))
                }
            }
        }

        location = list()
        covariate = list()

        count.linear = 0
        count.random = 0
        
        ## sort the models so copies are at the end
        if (nr > 1) {
            swapped = TRUE
            while (swapped) {
                swapped = FALSE
                for(r in 1:(nr-1)) {
                    if (
                        ##
                        ## not both are copy, then copy is at the end
                        ##
                        (gp$random.spec[[r]]$model == "copy" && gp$random.spec[[r+1]]$model != "copy") ||
                        ##
                        ## both are copy, then 'same.as' is at the end
                        ##
                        (gp$random.spec[[r]]$model == "copy" && gp$random.spec[[r+1]]$model == "copy" &&
                         !is.null(gp$random.spec[[r]]$same.as) && is.null(gp$random.spec[[r+1]]$same.as))) {

                        tmp = gp$random.spec[[r]]
                        gp$random.spec[[r]] = gp$random.spec[[r+1]]
                        gp$random.spec[[r+1]] = tmp

                        tmp = rf[, r+1]
                        rf[, r+1] = rf[, r+2]
                        rf[, r+2] = tmp
                        
                        swapped = TRUE
                    }
                }
            }
        }            
        if (debug) {
            for(r in 1:nr)
                cat(paste(r, gp$random.spec[[r]]$model, gp$random.spec[[r]]$of, gp$random.spec[[r]]$same.as), "\n")
        }
        
        all.terms = c()
        for(r in 1:nr)
            all.terms = c(all.terms, gp$random.spec[[r]]$term)
        
        for (r in 1:nr) {
            n = nrep = ngroup = N = NULL
            
            if (gp$random.spec[[r]]$model != "linear") {
                ##in this case we have to add a FFIELD section.........
                count.random = count.random+1
                xx=rf[, r +1]

                ## Do some checking for 'copy', which might lead
                ## to some strange errors if not met.  Its OK to
                ## check these here as we already have sorted the
                ## fields.
                
                ## if copy, then make sure that some critical variables are the same
                if (gp$random.spec[[r]]$model == "copy") {
                    r.idx = which(all.terms == gp$random.spec[[r]]$of)
                    if (length(r.idx) == 0)
                        stop(paste("Something wrong: cannot find the match of copy='", gp$random.spec[[r]]$of, "'", sep=""))
                    
                    gp$random.spec[[r]]$ngroup = gp$random.spec[[r.idx]]$ngroup
                    gp$random.spec[[r]]$nrep = gp$random.spec[[r.idx]]$nrep
                    gp$random.spec[[r]]$n = gp$random.spec[[r.idx]]$n
                    gp$random.spec[[r]]$values = gp$random.spec[[r.idx]]$values
                }

                if (!is.null(gp$random.spec[[r]]$nrep)) {
                    nrep = as.integer(gp$random.spec[[r]]$nrep)
                } else {
                    if (!is.null(gp$random.spec[[r]]$replicate)) {
                        nrep = as.integer(max(gp$random.spec[[r]]$replicate, na.rm=TRUE))
                    } else {
                        nrep = 1
                    }
                }

                if (!is.null(gp$random.spec[[r]]$ngroup)) {
                    ngroup = as.integer(gp$random.spec[[r]]$ngroup)
                } else {
                    if (!is.null(gp$random.spec[[r]]$group)) {
                        ngroup = as.integer(max(gp$random.spec[[r]]$group, na.rm = TRUE))
                    } else {
                        ngroup = 1
                    }
                }
                
                if (nrep > 1 || ngroup > 1) {
                    if (is.null(gp$random.spec[[r]]$n)) {
                        if (is.factor(xx)) {
                            n = length(unique(xx[!is.na(xx)]))
                        } else {
                            if (!is.null(gp$random.spec[[r]]$values))
                                n = length(unique(gp$random.spec[[r]]$values[ !is.na(gp$random.spec[[r]]$values) ]))
                            else
                                n = length(unique(xx[!is.na(xx)]))
                        }
                        gp$random.spec[[r]]$n = n
                    } else {
                        n = gp$random.spec[[r]]$n
                    }
                    ##warning(paste("inla: I do not think `n' is required to give...; compute n=", n))

                    N = inla.model.properties(gp$random.spec[[r]]$model, "latent")$aug.factor * n
                    replicate = gp$random.spec[[r]]$replicate
                    if (is.null(replicate))
                        replicate = rep(1, length(xx))
                    
                    group = gp$random.spec[[r]]$group
                    if (is.null(group))
                        group = rep(1, length(xx))
                    
                    if (nrep > 1) {
                        ## if replicate is a single number, then reuse
                        if (length(replicate) == 1)
                            replicate = rep(replicate, length(xx))

                        if (length(replicate) != length(xx))
                            stop(paste("length(replicate) != length(xx)", length(replicate), "!=",  length(xx)))
                        if (!all(is.element(replicate, 1:nrep) | is.na(replicate)))
                            stop(paste("Error in the values of `replicate'; not in [1,...,", nrep,"]", sep=""))
                        if (any(is.na(replicate[ !is.na(xx) ])))
                            stop(paste("There are one or more NA's in 'replicate' where 'idx' in f(idx,...) is not NA: idx = \'",
                                       gp$random.spec[[r]]$term, "\'", sep=""))
                        replicate[ is.na(xx) ] = 1
                    }

                    if (ngroup > 1) {
                        ## if group is a single number, then reuse
                        if (length(group) == 1)
                            group = rep(group, length(xx))

                        if (length(group) != length(xx))
                            stop(paste("length(group) != length(xx)", length(group), "!=",  length(xx)))
                        if (!all(is.element(group, 1:ngroup) | is.na(group)))
                            stop(paste("Error in the values of `group'; not in [1,...,", ngroup,"]", sep=""))
                        if (any(is.na( group[ !is.na(xx) ])))
                            stop(paste("There are one or more NA's in 'group' where 'idx' in f(idx,...) is not NA: idx = \'",
                                       gp$random.spec[[r]]$term, "\'", sep=""))
                        group[ is.na(xx) ] = 1

                        ## issue a WARNING,  if there are to many unused groups
                        g.used = unique(sort(group))
                        g.unused = setdiff(1:ngroup, g.used)
                        ng.used = length(g.used)
                        ng.unused = length(g.unused)

                        txt = paste("f(", gp$random.spec[[r]]$term, ", ...)",  sep="")
                        if (!is.element(1, g.used)) {
                            warning(paste(txt, ": ", 
                                          "There is no indices where group[]=1, this is *usually* a misspesification"),
                                    immediate. = TRUE)
                        }
                        if (ng.unused >= ng.used) {
                            warning(paste(txt, ": ", 
                                          "Number of unused groups >= the number of groups used: ", ng.unused, " >= ", ng.used, 
                                        ", this is *usually* a misspesification", sep=""),
                                    immediate. = TRUE)
                        }
                    }
                } else {
                    N = NULL
                }
                NG = N*ngroup

                ## might be used for plotting etc...
                gp$random.spec[[r]]$N = N
                gp$random.spec[[r]]$NG = NG
                gp$random.spec[[r]]$NREP = NREP = nrep

                ## make sure these are set so they will be used by the `copy'-feature
                if (is.null(gp$random.spec[[r]]$n)) {
                    gp$random.spec[[r]]$n = n
                }
                gp$random.spec[[r]]$ngroup = ngroup
                gp$random.spec[[r]]$nrep = nrep

                if (nrep > 1 || ngroup > 1) {
                    if (is.null(replicate))
                        stop("replicate is NULL")
                    if (is.null(group))
                        stop("group is NULL")
                }

                if (!is.null(gp$random.spec[[r]]$values)) {

                    ## in any case they must be unique:
                    if (length(gp$random.spec[[r]]$values) != length(unique(gp$random.spec[[r]]$values))) {
                        stop(paste("The 'values' given in f(",
                                   gp$random.spec[[r]]$term, ",..., values=...) are not unique:\n",
                                   "    length(values)=", length(gp$random.spec[[r]]$values),
                                   " but length(unique(values))=", length(unique(gp$random.spec[[r]]$values)), 
                                   sep=""))
                    }

                    ## values are given. then either both must be numeric or a factor,  or factor x character

                    if (is.numeric(xx) && is.numeric(gp$random.spec[[r]]$values)) {
                        ## case 1: both numeric

                        ## no sort for values here, since they are given as they should be !!!!
                        location[[r]] = unique(gp$random.spec[[r]]$values)
                        cov = match(xx, location[[r]])-1L + inla.ifelse(nrep > 1L || ngroup > 1L,  (replicate-1L)*NG + (group-1)*N, 0L)
                        
                        ## check that not do anything we will regret
                        i.cov = is.na(match(xx, location[[r]]))
                        i.xx = is.na(xx)
                        if (!all(i.cov == i.xx)) {
                            stop(paste("f(", gp$random.spec[[r]]$term, "). Covariate does not match 'values' ",
                                       sum(i.cov != i.xx), " times. Indexes for mismatch:", inla.paste(which(i.cov != i.xx)),
                                       ". This is not what you want. Use NA values in the covariate!",
                                       sep=""))
                        }
                        cov[is.na(cov)] = -1L
                        covariate[[r]] = cov
                    } else if ((is.factor(xx) || is.character(xx))
                               &&
                               (is.factor(gp$random.spec[[r]]$values) || is.character(gp$random.spec[[r]]$values))) {
                        ## case 2: both factors or character's

                        if (is.character(xx)) {
                            ## convert to factor (easier)
                            xx = factor(xx, unique(xx))
                        }
                        if (is.character(gp$random.spec[[r]]$values)) {
                            ## convert it to a factor
                            gp$random.spec[[r]]$values = factor(gp$random.spec[[r]]$values, gp$random.spec[[r]]$values)
                        }

                        gp$random.spec[[r]]$id.names = levels(gp$random.spec[[r]]$values)
                        location[[r]] = 1L:length(levels(gp$random.spec[[r]]$values))

                        ## let us first check if there are levels in
                        ## xx that is not present in values. this
                        ## should not happen as NA's should be used
                        ## instead

                        if (length(setdiff(levels(xx), gp$random.spec[[r]]$values)) != 0L) {
                            dif = setdiff(levels(xx), levels(gp$random.spec[[r]]$values))
                            stop(paste("In f(", gp$random.spec[[r]]$term, "): Covariate does not match 'values' ",
                                       length(dif), " times. Levels that mismatch:", inla.paste(c(dif),  sep=" "), 
                                       ".\n  This is not what you want. Use NA values in the covariate!",
                                       sep=""))
                        }

                        cov = match(xx, levels(gp$random.spec[[r]]$values)) -1L
                        cov[is.na(cov)] = -1L
                        covariate[[r]] = cov
                    } else {
                        ## all combinations not covered above is not allowed. 
                        stop(paste("In f(", gp$random.spec[[r]]$term, "): 'covariate' must match 'values'",
                                   ",  and both must either be 'numeric', or 'factor'/'character'.",  sep=""))
                    }
                } else {
                    ## values are not given. then 'values' depends on the type of the covariate.
                    if (is.factor(xx)) {
                        gp$random.spec[[r]]$id.names = levels(xx)
                        location[[r]] = 1L:length(levels(xx))
                        xx = as.numeric(xx)
                    } else if (is.character(xx)) {
                        gp$random.spec[[r]]$id.names = levels(as.factor(xx))
                        location[[r]] = 1L:length(levels(as.factor(xx)))
                        xx = as.numeric(sapply(xx, function(xx.i, id.names) which(xx.i == id.names), id.names = gp$random.spec[[r]]$id.names))
                    } else if (is.numeric(xx)) {
                        gp$random.spec[[r]]$id.names = NULL
                        location[[r]] = sort(unique(xx))
                        ## need to store the mapping for later use
                    } else {
                        stop(paste("f(", gp$random.spec[[r]]$term, "). Covariate is not of type 'factor', 'character' or 'numeric'.",  sep=""))
                    }

                    cov = match(xx, location[[r]])-1L + inla.ifelse(nrep > 1L || ngroup > 1L,  (replicate-1L)*NG + (group-1L)*N, 0L)

                    ## check that we do not do anything we will regret
                    i.cov = is.na(match(xx, location[[r]]))
                    i.xx = is.na(xx)
                    if (!all(i.cov == i.xx)) {
                        stop(paste("f(", gp$random.spec[[r]]$term, "). Covariate does not match 'values' ",
                                   sum(i.cov != i.xx), " times. Indexes for mismatch:", inla.paste(which(i.cov != i.xx)),
                                   ". This is not what you want. Use NA values in the covariate!",
                                   sep=""))
                    }

                    cov[is.na(cov)] = -1L
                    covariate[[r]] = cov
                }

                ## just make sure this is all set
                if (is.null(gp$random.spec[[r]]$values)) {
                    gp$random.spec[[r]]$values = location[[r]]
                }

                ## create a location and covariate file
                file.loc=inla.tempfile(tmpdir=data.dir)
                if (inla.getOption("internal.binary.mode")) {
                    inla.write.fmesher.file(as.matrix(as.numeric(location[[r]]), ncol = 1),  filename = file.loc, debug = debug)
                    ## prevent some numerical instabilities for models rw1, rw2, crw2, etc...
                    inla.check.location(location[[r]], term = gp$random.spec[[r]]$term,
                                        model = gp$random.spec[[r]]$model, section = "latent")
                } else {
                    file.create(file.loc)
                    write(as.numeric(location[[r]]), ncolumns=1, file=file.loc, append=FALSE)
                }
                file.loc = gsub(data.dir, "$inladatadir", file.loc, fixed=TRUE)
                
                ## this have to match
                stopifnot(length(covariate[[r]]) == NPredictor)
                file.cov=inla.tempfile(tmpdir=data.dir)
                if (inla.getOption("internal.binary.mode")) {
                    inla.write.fmesher.file(as.matrix(cbind(indN, covariate[[r]])), filename=file.cov, debug = debug)
                } else {
                    file.create(file.cov)
                    write(t(cbind(indN, covariate[[r]])), ncolumns=2, file=file.cov, append=FALSE)
                }
                file.cov = gsub(data.dir, "$inladatadir", file.cov, fixed=TRUE)

                ## name of the 'names of the values'
                if (!is.null(gp$random.spec[[r]]$id.names)) {
                    file.id.names = inla.tempfile(tmpdir=data.dir)
                    inla.writeLines(file.id.names, gp$random.spec[[r]]$id.names)
                    file.id.names = gsub(data.dir, "$inladatadir", file.id.names, fixed=TRUE)
                } else {
                    file.id.names = NULL
                }

                file.cov=inla.tempfile(tmpdir=data.dir)
                if (inla.getOption("internal.binary.mode")) {
                    inla.write.fmesher.file(as.matrix(cbind(indN, covariate[[r]])), filename=file.cov, debug = debug)
                } else {
                    file.create(file.cov)
                    write(t(cbind(indN, covariate[[r]])), ncolumns=2, file=file.cov, append=FALSE)
                }
                file.cov = gsub(data.dir, "$inladatadir", file.cov, fixed=TRUE)


                if (nrep == 1 && ngroup == 1) {
                    n = inla.ifelse(is.null(gp$random.spec[[r]]$n), length(location[[r]][!is.na(location[[r]])]),
                        gp$random.spec[[r]]$n)
                } else {
                    if (is.null(n))
                        stop("n is NULL!!!!")
                }
                if (is.null(gp$random.spec[[r]]$n))
                    gp$random.spec[[r]]$n = n

                if (debug) {
                    print(paste("n", n))
                    print(paste("nrep", nrep))
                }

                n.div.by = inla.model.properties(gp$random.spec[[r]]$model, "latent")$n.div.by
                if (!is.null(n.div.by)) {
                    if (!inla.divisible(n, by=n.div.by))
                        stop(paste("Model [", gp$random.spec[[r]]$model,"] require `n'", n, "divisible by", n.div.by))
                }

                if (nrep < 1 || nrep != round(nrep))
                    stop(paste("\nArgument NREP is void:", nrep))
                if (ngroup < 1 || ngroup != round(ngroup))
                    stop(paste("\nArgument NGROUP is void:", ngroup))

                ## for the CRW2/BYM (or agumented) models, we have to
                ## implement the constr=TRUE using the argument
                ## EXTRACONSTR

                if (inla.model.properties(gp$random.spec[[r]]$model, "latent")$augmented) {
                    fac = inla.model.properties(gp$random.spec[[r]]$model, "latent")$aug.factor
                    con = inla.model.properties(gp$random.spec[[r]]$model, "latent")$aug.constr

                    if (fac > 1) {
                        if ((max(con) > fac || min(con) < 1) || is.null(con))
                            stop(paste("INTERNAL ERROR: aug.constr = ", con))
                        
                        ## this case is a very special case
                        if (gp$random.spec[[r]]$constr) {
                            if (is.null(gp$random.spec[[r]]$extraconstr)) {
                                A = matrix(0, 1, fac*n)
                                for(con.elm in con)
                                    A[1, (con.elm-1)*n + 1:n] = 1
                                e = 0
                            } else {
                                nrow = dim(gp$random.spec[[r]]$extraconstr$A)[1]
                                ncol = dim(gp$random.spec[[r]]$extraconstr$A)[2]

                                if (ncol != fac*n)
                                    stop(paste("Wrong dimension for the extraconstr: ncol", ncol, "n", n))
                                
                                A = matrix(0, nrow+1, ncol)
                                e = c(gp$random.spec[[r]]$extraconstr$e, 0)
                                A[1:nrow, 1:ncol] = gp$random.spec[[r]]$extraconstr$A

                                for(con.elm in con)
                                    A[nrow+1, (con.elm-1)*n + 1:n] = 1
                            }
                            gp$random.spec[[r]]$extraconstr = list(A=A, e=e)
                            gp$random.spec[[r]]$constr = FALSE
                        }
                    } else {
                        stopifnot(fac == 1)

                        ## this is the case fac = 1. this case is a
                        ## very special case, for iid2d, iid3d, etc

                        n.small = n %/% max(con)

                        if (gp$random.spec[[r]]$constr) {
                            if (is.null(gp$random.spec[[r]]$extraconstr)) {
                                A = matrix(0, length(con), n)
                                k = 1
                                for(con.elm in con) {
                                    A[k, (con.elm-1)*n.small + 1:n.small] = 1
                                    k = k + 1
                                }
                                e = rep(0, length(con))
                            } else {
                                nrow = dim(gp$random.spec[[r]]$extraconstr$A)[1]
                                ncol = dim(gp$random.spec[[r]]$extraconstr$A)[2]

                                if (ncol != n)
                                    stop(paste("Wrong dimension for the extraconstr: ncol", ncol, "n", n))
                                
                                A = matrix(0, nrow+length(con), ncol)
                                e = c(gp$random.spec[[r]]$extraconstr$e, rep(0, length(con)))
                                A[1:nrow, 1:ncol] = gp$random.spec[[r]]$extraconstr$A

                                k = 1
                                for(con.elm in con) {
                                    A[nrow+k, (con.elm-1)*n.small + 1:n.small] = 1
                                    k = k + 1
                                }
                            }
                            gp$random.spec[[r]]$extraconstr = list(A=A, e=e)
                            gp$random.spec[[r]]$constr = FALSE
                        }
                    }
                } 

                ##print(gp$random.spec[[r]]$extraconstr$A)
                ##print(gp$random.spec[[r]]$extraconstr$e)
                
                ##and in case a file for the extraconstr
                if (!is.null(gp$random.spec[[r]]$extraconstr)) {
                    A=gp$random.spec[[r]]$extraconstr$A
                    e=gp$random.spec[[r]]$extraconstr$e

                    if ((gp$random.spec[[r]]$model != "rgeneric") &&
                        (ncol(A) != inla.model.properties(gp$random.spec[[r]]$model, "latent")$aug.factor*n))
                        stop(paste("\n\tncol in matrix A(extraconstr) does not correspont to the length of f:",
                                   ncol(A),
                                   inla.model.properties(gp$random.spec[[r]]$model, "latent")$aug.factor*n))

                    file.extraconstr=inla.tempfile(tmpdir=data.dir)
                    if (inla.getOption("internal.binary.mode")) {
                        inla.write.fmesher.file(as.matrix(c(as.vector(t(A)), e), ncol=1), filename=file.extraconstr, debug=debug)
                    } else {
                        file.create(file.extraconstr)
                        write(c(as.vector(t(A)), e), ncolumns=1, file=file.extraconstr, append=FALSE)
                    }
                    file.extraconstr = gsub(data.dir, "$inladatadir", file.extraconstr, fixed=TRUE)
                } else {
                    file.extraconstr = NULL
                }
                
                ##....also if necessary a file for the weights (not to be confused with argument 'weights' in the inla() call...)
                www = NULL
                if (!is.null(gp$random.spec[[r]]$weights)) {
                    ## $weights is the name
                    www = wf[, gp$random.spec[[r]]$weights ]
                    www[is.na(www)] = 0

                    ##create a file for the weights
                    file.weights=inla.tempfile(tmpdir=data.dir)
                    if (inla.getOption("internal.binary.mode")) {
                        inla.write.fmesher.file(as.matrix(cbind(indN, www)), filename=file.weights, debug = debug)
                    } else {
                        file.create(file.weights)
                        write(t(cbind(indN, www)), ncolumns=2, file=file.weights, append=FALSE)
                    }
                    file.weights = gsub(data.dir, "$inladatadir", file.weights, fixed=TRUE)

                    n.weights = n.weights+1
                }

                ## for some models, the priors are computed in this function. in these cases, we
                ## need to updated all.hyper with these priors, not the ones that goes into this
                ## function...
                rs.updated = (inla.ffield.section(file=file.ini,
                                                  file.loc=file.loc, file.cov=file.cov,
                                                  file.id.names = file.id.names, 
                                                  file.extraconstr=file.extraconstr, 
                                                  file.weights=file.weights, n=n, nrep = nrep, ngroup = ngroup,
                                                  random.spec=gp$random.spec[[r]], 
                                                  results.dir=paste("random.effect", inla.num(count.random), sep=""), 
                                                  only.hyperparam= only.hyperparam,
                                                  data.dir=data.dir, www = www))
                all.hyper$random[[r]] = (list(hyperid = inla.namefix(gp$random.spec[[r]]$term),
                                              hyper = rs.updated$hyper,
                                              group.hyper = rs.updated$control.group$hyper))


            } else if (inla.one.of(gp$random.spec[[r]]$model, "linear")) {
                ##....while here we have to add a LINEAR section
                count.linear = count.linear+1
                xx=rf[, r +1]
                xx[is.na(xx)] = 0
                file.linear = inla.tempfile(tmpdir=data.dir)
                if (inla.getOption("internal.binary.mode")) {
                    inla.write.fmesher.file(as.matrix(cbind(indN, xx)), filename=file.linear, debug=debug)
                } else {
                    file.create(file.linear)
                    write(t(cbind(indN, xx)), ncolumns=2, file=file.linear, append=FALSE)
                }
                file.linear = gsub(data.dir, "$inladatadir", file.linear, fixed=TRUE)

                cont = list(cdf=gp$random.spec[[r]]$cdf,
                    quantiles=gp$random.spec[[r]]$quantiles,
                    prec=gp$random.spec[[r]]$prec.linear,
                    mean=gp$random.spec[[r]]$mean.linear,
                    compute=gp$random.spec[[r]]$compute)
                
                if (is.null(all.hyper$linear)) {
                    lin.count = 1L
                } else {
                    lin.count = length(all.hyper$linear) + 1L
                }
                all.hyper$linear[[lin.count]] = inla.linear.section(
                                    file=file.ini, file.fixed=file.linear, label=gp$random.spec[[r]]$term,
                                    results.dir=paste("fixed.effect", inla.num(gp$n.fix+count.linear), sep=""),
                                    control.fixed = cont, only.hyperparam=only.hyperparam)
            } else {
                stop("This should not happen.")
            }
        }
    }
    
    if (n.weights != gp$n.weights) 
        stop("\n\tSomething strange with weights in the covariate...")

    ## the inla section
    inla.inla.section(file=file.ini, inla.spec=cont.inla, data.dir)

    ## create mode section
    cont.mode = inla.set.control.mode.default()
    cont.mode[names(control.mode)] = control.mode
    inla.mode.section(file=file.ini, cont.mode, data.dir)

    ## create expert section
    cont.expert = inla.set.control.expert.default()
    cont.expert[names(control.expert)] = control.expert
    inla.expert.section(file=file.ini, cont.expert, data.dir = data.dir)
    
    ## create lincomb section
    cont.lincomb = inla.set.control.lincomb.default()
    cont.lincomb[names(control.lincomb)] = control.lincomb
    inla.lincomb.section(file=file.ini, data.dir=data.dir, contr=cont.lincomb, lincomb=lincomb)
    
    ## create update section
    cont.update = inla.set.control.update.default()
    cont.update[names(control.update)] = control.update
    inla.update.section(file=file.ini, data.dir=data.dir, contr=cont.update)
    
    ## now, do the job
    if (debug) {
        cat("...done\n")
        cat("Run inla...")
    }

    ## inla.arg override all default arguments including `-b' !!!
    if (is.null(inla.arg)) {
        arg.arg = ""
        arg.nt = inla.ifelse(is.numeric(num.threads), paste(" -t ", num.threads, " ", sep=""), "")

        ## due to the weird behaviour,  we will do the verbose-mode differently now
        if (inla.os("linux") || inla.os("mac")) {
            arg.v = inla.ifelse(verbose, "-v", "-v")
        } else {
            arg.v = inla.ifelse(verbose, "-v", "-v")
        }

        arg.s = inla.ifelse(silent, "-s", "")
        arg.b = "-b"
    } else {
        arg.arg = inla.arg
        arg.nt = ""
        arg.v = ""
        arg.s = ""
        arg.b = ""
    }
    
    ## collect all. we might add '-p' later if inla.call="submit"
    all.args = paste(arg.arg, arg.b, arg.s, arg.v, arg.nt, sep=" ")

    ## define some environment variables for remote computing
    vars = list(INLA_PATH = system.file("bin", package="INLA"),
                INLA_OS = inla.os.type(), 
                INLA_HGVERSION = inla.version("hgid"), 
                INLA_RVERSION = paste0(R.Version()$major, ".",
                                       strsplit(R.Version()$minor,"[.]")[[1]][1]),
                INLA_RHOME = Sys.getenv("R_HOME"))
    do.call("Sys.setenv", vars)
    inla.set.sparselib.env(inla.dir, blas.num.threads = blas.num.threads)

    vars = NULL
    if (debug) {
        vars = c(vars, INLA_DEBUG=1)
    }
    if (remote || submit) {
        if (submit) {
            all.args = paste(all.args,  "-p") ## need this option
            vars = c(vars,
                     INLA_SUBMIT_ID = submit.id)
        }
        if (inla.os("windows")) {
            vars = c(vars, 
                     INLA_SSH_AUTH_SOCK = inla.getOption("ssh.auth.sock"), 
                     INLA_CYGWIN_HOME = inla.getOption("cygwin.home"), 
                     INLA_HOME = inla.cygwin.map.filename(gsub("\\\\", "/", inla.get.HOME())))
        } else {
            vars = c(vars,
                     INLA_HOME = inla.get.HOME())
            if (Sys.getenv("SSH_AUTH_SOCK") == "") {
                vars = c(vars,
                         INLA_SSH_AUTH_SOCK = inla.getOption("ssh.auth.sock"))
            }
        }
    }
    if (!is.null(vars))
        do.call("Sys.setenv", as.list(vars))

    ## write the list of environment variables set, so they can be reset if needed
    env = Sys.getenv()
    env.n = names(env)
    idx = grep("^(INLA_|(OPENBLAS|MKL)_NUM_THREADS|PARDISO)", env.n)
    env.list = env[idx]
    file.env = paste0(inla.dir, "/environment")
    cat(file=file.env)
    for(i in seq_along(env.list)) {
        cat(names(env.list[i]), "=\"", env.list[i], "\"\n", sep="", file=file.env, append=TRUE)
    }

    my.time.used[2] = Sys.time()
    ## ...meaning that if inla.call = "" then just build the files (optionally...)
    if (ownfun || nchar(inla.call) > 0) {
        if (ownfun) {
            ## undocumented feature for PB
            echoc = inla.call(file.ini = file.ini,
                              file.log = if (verbose) NULL else file.log,
                              results.dir = results.dir,
                              inla.call.args = all.args)
        } else if (inla.os("linux") || inla.os("mac")) {
            if (verbose) {
                echoc = system(paste(shQuote(inla.call), all.args, shQuote(file.ini)))
            } else {
                echoc = system(paste(shQuote(inla.call), all.args, shQuote(file.ini), " > ", shQuote(file.log),
                    inla.ifelse(silent == 2L, " 2>/dev/null", "")))
            }
        } else if (inla.os("windows")) {
            if (!remote && !submit) {
                if (verbose) {
                    echoc = try(system2(inla.call, args=paste(all.args, shQuote(file.ini)), stdout="", stderr="", wait=TRUE))
                } else {
                    if (FALSE) {
                        ## old .bat-solution
                        bat.file = paste(tempfile(), ".BAT",  sep="")
                        cat("@echo off\n",  file=bat.file, append=FALSE)
                        cat(paste(shQuote(inla.call), all.args, "-v", shQuote(file.ini), ">", shQuote(file.log),
                                  inla.ifelse(silent == 2L, "2>NUL", "")), file=bat.file, append=TRUE)
                        echoc = try(system2(bat.file, wait=TRUE), silent=FALSE)
                        unlink(bat.file)
                    } else {
                        ## new try
                        echoc = try(system2(inla.call, args=paste(all.args, shQuote(file.ini)), stdout=file.log, stderr=file.log2, wait=TRUE))
                    }
                }
                if (echoc != 0L) {
                    if (!verbose && (silent != 2L)) {
                        inla.inlaprogram.has.crashed()
                    }
                }
            } else {
                ## remote || submit
                echoc = try(inla.cygwin.run.command(
                    paste(inla.cygwin.map.filename(inla.call),
                          all.args,
                          inla.cygwin.map.filename(file.ini)),
                    file.log = inla.ifelse(verbose, NULL, inla.cygwin.map.filename(file.log))), silent=TRUE)
                ## echoc = 0L
            }
        } else {
            stop("\n\tNot supported architecture.")
        }

        if (debug) {
            cat("..done\n")
        }

        my.time.used[3] = Sys.time()

        if (echoc == 0L) {
            if (!submit) {
                ret = try(inla.collect.results(results.dir, control.results=cont.results, debug=debug,
                    only.hyperparam=only.hyperparam, file.log = file.log, file.log2=file.log2), silent=FALSE)
                if (!is.list(ret)) {
                    ret = list()
                }
            } else {
                ret = list()
            }
            
            my.time.used[4] = Sys.time()
            cpu.used = c(
                "Pre" = diff(my.time.used)[1],
                "Running" = diff(my.time.used)[2],
                "Post" = diff(my.time.used)[3],
                "Total" = my.time.used[4] - my.time.used[1])

            ret$cpu.used = cpu.used
            ## store all arguments; replacing 'control.xxx' with 'cont.xxx'
            the.args = list()
            for (nm in names(formals(inla))) {
                nnm = nm
                nnm = gsub("^control\\.", "cont.", nnm) ## these are the processed ones
                nnm = gsub("^data$", "data.orig", nnm) 
                nnm = gsub("^formula$", "formula.orig", nnm) 
                nnm = gsub("^cont(rol)?\\.family$", "control.family.orig", nnm)
                inla.eval(paste("the.args$", nm, " = ", nnm, sep=""))
            }
            ## remove the .Evironment attribute, as it will fail if
            ## its rerun if .Environment is not there.
            attr(the.args$formula, ".Environment") = NULL
            
            ## further fix for $control.family until the c(param = numeric(0)) issue is solved. 
            if (n.family > 1L) {
                ## we need to fix the case where there this option is
                ## not set. otherwise, we will get list(), instead of
                ## list(list(), list()), say, for n.family=2
                if (is.list(the.args$control.family) && length(the.args$control.family) == 0L) {
                    the.args$control.family = lapply(1:n.family, function(x) list())
                }
            }
            ## OLD CODE: if (n.family == 1) the.args$control.family = the.args$control.family[[1L]]
            ret$all.hyper = all.hyper
            ret$.args = the.args
            ret$call = call
            ret$model.matrix = gp$model.matrix
            class(ret) = "inla"
        } else {
            ## do this instead
            if (silent != 2L) {
                inla.inlaprogram.has.crashed()
            } else {
                ## with a crash, try to collect the logfile only
                ret1 = try(inla.collect.logfile(file.log, debug), silent = TRUE)
                ret2 = try(inla.collect.logfile(file.log2, debug), silent = TRUE)
                if (inherits(ret1, "try-error")) { ## yes,  its 'ret1'
                    ret = NULL
                } else {
                    ret = list(logfile = c(ret1$logfile,
                                           "", paste(rep("*",72),sep="",collapse=""), "", 
                                           ret2$logfile))
                    class(ret) = "inla"
                }
            }
        }

        if (debug && !keep) {
            cat("clean up\n")
        }
        if (!keep && !submit) {
            try(unlink(inla.dir, recursive=TRUE), silent = TRUE)
        }
    } else {
        ret = NULL
    }
        
    ##
    if (submit) {
        ret$misc$inla.dir = inla.dir
        ret.sub = list(ret = ret, id = submit.id)
    } else {
        return (ret)
    }
}


`inla.fix.data` = function(data, n, revert = FALSE)
{
    ## extract all entries in 'data' with length='n'. if 'revert',  do the oposite
    if (is.data.frame(data)) {
        if (dim(data)[1L] == n) {
            return (inla.ifelse(revert, data.frame(), data))
        } else {
            return (inla.ifelse(!revert, data.frame(), data))
        }
    }
    if (is.list(data) && length(data) > 0L) {
        idx = which(sapply(
            data, 
            function(a, n) {
                if (inla.is.matrix(a) && (length(dim(a)) > 1L) && (dim(a)[1L] == n)) {
                    return (TRUE)
                } else if (length(a) != n || inla.is.matrix(a)) {
                    return (FALSE)
                } else if (length(a) == n) {
                    return (TRUE)
                } else {
                    stop("This should not happen.")
                }
            },
            n = n))
        if (revert) {
            idx = (1:length(data))[-idx]
        }
        if (length(idx) > 0L) {
            return (data[idx])
        } else {
            return (inla.ifelse(revert || is.list(data), list(), data.frame()))
        }
    } else {
        return (data)
    }
    stop("Should not happen.")
}

`inla.expand.factors` = function(data, exclude.names = c())
{
    ## replace factors in 'data' with their expanded version we get
    ## from model.matrix without any intercept.
    for(k in seq_along(data)) {
        if (is.factor(data[[k]])) {
            if (!(names(data)[k] %in% exclude.names)) {
                formula = as.formula(paste("~ -1 + ",  names(data)[k]))
                tmp = model.matrix(formula, model.frame(formula,  data, na.action = na.pass))
                colnames(tmp) = paste0(names(data)[k], levels(data[[k]]))
                data[[k]] = tmp
            }
        }
    }
    return (data)
}

`inla.set.sparselib.env` = function(inla.dir = NULL, blas.num.threads = 1L) 
{
    ## environment variables for sparse libraries
    if (is.null(inla.dir)) {
        inla.dir = inla.tempdir()
    }

    lic.filename = "pardiso.lic" ## do not change
    lic.filename.dir = paste0(inla.dir, "/", lic.filename)
    file.create(lic.filename.dir)
    
    if (!is.null(inla.getOption("pardiso.license"))) {
        lic.file = normalizePath(inla.getOption("pardiso.license"))
        lic.path = NA
        if (file.exists(lic.file)) {
            info = file.info(lic.file)
            if (!is.na(info$isdir)) {
                if (info$isdir) {
                    lic.path = lic.file
                } else if (!is.null(inla.dir)) {
                    file.copy(lic.file, lic.filename.dir, overwrite=TRUE)
                    lic.path = inla.dir
                } else {
                    stop("This should not happen")
                }
            } else {
                lic.path = lic.file
            }
        } else {
            lic.path = lic.file
        }
        do.call("Sys.setenv", list(PARDISO_LIC_PATH = normalizePath(lic.path), 
                                   INLA_LOAD_PARDISO = 1))
    } else {
        Sys.unsetenv("INLA_LOAD_PARDISO")
    }

    if (Sys.getenv("PARDISOLICMESSAGE") == "")
        Sys.setenv(PARDISOLICMESSAGE=1)
    
    blas.num.threads = as.integer(blas.num.threads)
    if (blas.num.threads == 0) {
        Sys.unsetenv("OPENBLAS_NUM_THREADS")
        Sys.unsetenv("MKL_NUM_THREADS")
    } else if (blas.num.threads > 0) {
        Sys.setenv(OPENBLAS_NUM_THREADS = blas.num.threads)
        Sys.setenv(MKL_NUM_THREADS = blas.num.threads)
    } else {
        if (Sys.getenv("OPENBLAS_NUM_THREADS") == "")
            Sys.setenv(OPENBLAS_NUM_THREADS = abs(blas.num.threads))
        if (Sys.getenv("MKL_NUM_THREADS") == "")
            Sys.setenv(MKL_NUM_THREADS = abs(blas.num.threads))
    }
        
    return (invisible())
}    
