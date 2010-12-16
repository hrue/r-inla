##! \name{inla}
##! \alias{inla}
##! \alias{xinla}

##! \title{Bayesian analysis of
##! structured additive models} \description{\code{inla} performs a
##! full Bayesian analysis of additive models using Integrated Nested
##! Laplace approximation
##! }

##! \usage{
##! inla = function (formula,
##!              family = "gaussian", 
##!              data = data.frame(),...)
##! }

##! \arguments{
    
`inla` =
    function (
              ##! \item{formula}{ A \code{inla} formula like \code{y
              ##!~1 + z + f(ind, model="iid")} + f(ind2,
              ##!weights,model="ar1") This is much like the formula
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
              ##! link. See \code{names(inla.models()$lmodels)} for a
              ##! list of possible alternatives.}
              family = "gaussian", 

              ##!\item{data}{ A data frame or list containing the
              ##!variables in the model.  The data frame MUST be
              ##!provided}
              data = NULL,

              ##!\item{quantiles}{ A vector of quantiles,
              ##!\eqn{p(0),p(1),\dots}{p(0),p(1),\ldots} to compute
              ##!for each posterior marginal. The function returns,
              ##!for each posterior marginal, the values
              ##!\eqn{x(0),x(1),\dots}{x(0),x(1),\ldots} such that
              ##!\deqn{\mbox{Prob}(X<x(p))=p}{Prob(X<x)=p} }
              quantiles=c(0.025, 0.5, 0.975),

              ##!\item{E}{ Known component in the mean for the Poisson
              ##!likelihoods defined as \deqn{E_i\exp(\eta_i)}{E
              ##!exp(eta)} where \deqn{\eta_i}{eta} is the linear
              ##!predictor. If not provided it is set to 1 }
              E = NULL,

              ##!\item{offset}{ This can be used to specify an
              ##!a-priori known component to be included in the linear
              ##!predictor during fitting.  This should be \code{NULL}
              ##!or a numeric vector of length either one or equal to
              ##!the number of cases. One or more \code{offset()}
              ##!terms can be included in the formula instead or as
              ##!well, and if both are specified their sum is used.}
              offset=NULL,

              ##!\item{scale}{ Fixed (optional) scale parameters of
              ##!the precision for Gaussian and Student-T response
              ##!models. Default value is 1.}
              scale = NULL,

              ##!\item{Ntrials}{ A vector containing the number of
              ##!trials for the \code{binomial} likelihood. Default
              ##!value is 1}
              Ntrials = NULL,

              ##!\item{verbose}{
              ##!Boolean indicating if the \code{inla}-program should
              ##!run in a verbose mode (default \code{FALSE}).}
              verbose = FALSE,

              ##!\item{lincomb}{ Used to define linear combination of
              ##!nodes in the latent field. The posterior distribution
              ##!of such linear combination is computed by the
              ##!\code{inla} function. See
              ##!\url{http://www.r-inla.org/help/faq} for examples of
              ##!how to define such linear combinations.}
              lincomb = NULL,

              ##!\item{control.compute}{ See \code{?control.compute}}
              control.compute = list(),

              ##!\item{control.predictor}{ See
              ##!\code{?control.predictor}}
              control.predictor = list(),

              ##!\item{control.data}{ See \code{?control.data}}
              control.data = list(),

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

              ##!\item{only.hyperparam}{ A boolean variable saying if
              ##!only the hyperparameters are to be computed. Mainly
              ##!for internal use.}
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
              ##!\code{inla}-program will use. xFor Windows this
              ##!defaults to 1, otherwise its defaults to \code{NULL}
              ##!(for which the system takes over control).}
              num.threads = inla.getOption("num.threads"),

              ##!\item{keep}{ A boolean variable indicating that the
              ##!working files (ini file, data files and results
              ##!files) should be kept. If TRUE and no
              ##!\code{working.directory} is specified the working
              ##!files are stored in a directory called "inla".  }
              keep = inla.getOption("keep"),

              ##!\item{working.directory}{ A string giving the name
              ##!of an alternative directory where to store the
              ##!working files }
              working.directory = inla.getOption("working.directory"),

              ##!\item{silent}{ A boolean variable defining how the
              ##!\code{inla}-program should be ``silent''.  }
              silent = inla.getOption("silent"),

              ##!\item{debug}{ If \code{TRUE}, then enable some debug
              ##!output.  }
              debug = inla.getOption("debug"),

              ##!\item{user.hook}{ This defines an optional
              ##!user-defined function, which can be called just after
              ##!the .ini-file is created, usually, to add extra
              ##!information to the .ini-file. See the function
              ##!\code{inla.user.hook} for an example of a such
              ##!function, and the arguments to it.  }
              user.hook = NULL,

              ##!\item{user.hook.arg}{ This defines an optional
              ##!argument to \code{user.hook} }
              user.hook.arg = NULL,
              ##
              ## these are ``internal'' options, used to transfer info
              ## from expansions
              ##
              .internal = list()
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
    ##!(plus, possibly quantiles and cdf) of the linear predictors
    ##!\eqn{\eta} in
    ##!the model }

    ##!\item{marginals.linear.predictor}{
    ##! If \code{compute=TRUE} in
    ##!\code{control.predictor}, a list containing the posterior
    ##!marginals of the linear predictors \eqn{\eta} in the model }

    ##!\item{summary.fitted.values}{
    ##! A matrix containing the mean and sd
    ##!(plus, possibly quantiles and cdf) of the fitted values
    ##!\eqn{g^{-1}(\eta)} obtained by
    ##!transforming the linear predictors by the inverse of the link
    ##!function.
    ##! }

    ##!\item{marginals.fitted.values}{
    ##! If \code{compute=TRUE} in
    ##!\code{control.predictor}, a list containing the posterior
    ##!marginals  of the fitted values
    ##!\eqn{g^{-1}(\eta)} obtained by
    ##!transforming the linear predictors by the inverse of the link
    ##!function.
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
    ##!If \code{cpo}=\code{TRUE} in \code{control.compute}, the
    ##!values of conditional predictive ordinate (CPO) otherwise \code{NULL}.
    ##!}

    ##!\item{pit}{
    ##!If \code{cpo}=\code{TRUE} in \code{control.compute}, the
    ##!values of the probability
    ##!integral transform (PIT) for the model, otherwise \code{NULL}.
    ##!}

    ##!\item{mlik}{
    ##!If \code{mlik}=\code{TRUE} in \code{control.compute}, the
    ##! marginal likelihood of the model, otherwise \code{NULL}
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
    ##!
    
    ##!Rue, H and Held, L. (2005) \emph{Gaussian Markov Random Fields
    ##!- Theory and Applications} Chapman and Hall
    ##!
    
    ##!Martino, S. and Rue, H. (2008) \emph{Implementing Approximate
    ##!Bayesian Inference using Integrated Nested Laplace
    ##!Approximation: a manual for the inla program} Preprint N.2 from
    ##!Dep. of Mathematical Sciences (NTNU Norway)
    ##!}

    ##!\author{Sara Martino and Havard Rue \email{hrue@math.ntnu.no} }
    
    ##!\seealso{\code{\link{f}}, 
    ##!\code{\link{inla.hyperpar}} }

    ##!\examples{
    ##! \dontrun{See the web page \url{www.r-inla.org} for a series of
    ##!worked out examples}}


{
    my.time.used = numeric(4)
    my.time.used[1] = Sys.time()
    
    if (nargs() == 0) {
        cat("\tUsage: inla(formula, family, data, other.arguments...); see ?inla\n")
        return (NULL)
    }

    if (is.null(data))
        stop("\t\tMissing data.frame argument `data'. Leaving `data' empty might lead to\n\t\tuncontrolled behaviour, therefore is it required.")

    ## check all control.xx arguments here...
    inla.check.control(control.compute)
    inla.check.control(control.predictor)
    ## do not check control.data here, as we need to know n.family
    ## first. do this later: inla.check.control(control.data)
    inla.check.control(control.inla)
    inla.check.control(control.results)
    inla.check.control(control.fixed)
    inla.check.control(control.mode)
    inla.check.control(control.expert)
    inla.check.control(control.hazard)
    inla.check.control(control.lincomb)

    n.family = length(family)
    for(i. in 1:n.family)
        family[i.] = inla.trim.family(family[i.])

    ## if the user specify inla.call="remote" or "inla.remote" then
    ## use the internal one
    remote = FALSE
    if (inla.strcasecmp(inla.call, "remote") ||
        inla.strcasecmp(inla.call, "inla.remote") ||
        length(grep("/inla.remote$", inla.call)) > 0 ||
        length(grep("/inla.remote.cygwin$", inla.call)) > 0) {
        remote = TRUE
        inla.call = system.file("bin/remote/inla.remote", package="INLA")
        if (inla.os("windows"))
            inla.call = paste(inla.call, ".cygwin", sep="")
    }
    
    if (!is.data.frame(data) && !is.list(data))
        stop("\n\tArgument `data' must be a data.frame or a list")
    
    ##
    ## check for survival model with a baseline-hazard. if so, then
    ## expand the data-frame and call inla() again.
    ##
    have.surv = FALSE
    for(i in 1:n.family)
        have.surv = have.surv || inla.lmodel.properties(family[i])$survival

    if (have.surv && (inla.one.of(family,c("coxph")))) {
        ## in this case, we expand the data-frame into a sequence of
        ## Poisson observations, and call inla() again.

        cont.hazard = inla.set.control.hazard.default()
        cont.hazard[names(control.hazard)] = control.hazard
      
        y.surv = eval(parse(text=formula[2]),data)
        if (class(y.surv) != "inla.surv")
            stop(paste("For survival models, then the reponse has to be of class `inla.surv'; you have `",
                       class(y.surv), "'", sep=""))
        data.orig = data
        data = inla.remove(as.character(formula[2]), data)
        if (is.null(y.surv$subject)) {
            res = inla.expand.dataframe.1(y.surv, as.data.frame(data), control.hazard = cont.hazard)
            new.data = res$data
            .internal$baseline.hazard.cutpoints = res$cutpoints
        } else {
            res = inla.expand.dataframe.2(y.surv, as.data.frame(data), control.hazard = cont.hazard)
            new.data = res$data
            .internal$baseline.hazard.cutpoints = res$cutpoints
        }
        data = data.orig
        
        if (cont.hazard$model != "rw1" && cont.hazard$model != "rw2")
            stop(paste("Valid models for the hazard is `rw1' and `rw2'; you have", cont.hazard$model))
        
        ## This is the strata-part, making the baseline hazard
        ## replicates according to the strata.
        strata.var = NULL
        if (!is.null(cont.hazard$strata.name)) {
            if (is.character(cont.hazard$strata.name) && length(cont.hazard$strata.name)==1) {
                ## strata = "x"
                strata.var = cont.hazard$strata
            } else {
                stop("Argument to `strata.name' must be the name of a variable in the data.frame.")
            }
        }
        if (debug)
            print(paste("strata.var", strata.var))
        if (!is.null(strata.var)) {
            if (!is.element(strata.var, names(new.data))) {
                stop(inla.paste(c("Variable `", strata.var,
                                  "' in control.hazard=list(strata=...) needs to be in the data.frame: names(data) = ",
                                  names(new.data))))
            }

            if (debug) print("apply inla.strata() on strata.var")
            inla.eval(paste("strata.tmp = inla.strata(new.data$", strata.var, ")", sep=""))
            inla.eval(paste("new.data$", strata.var, " = strata.tmp$strata"))
            inla.eval(paste(".internal$baseline.hazard.strata.coding = strata.tmp$coding"))
        }
        
        ## make sure that the dimension of the baseline-hazard is
        ## correct, but setting 'values' correct.
        if (!is.null(cont.hazard$cutpoints)) {
            baseline.hazard.values = seq(1, length(cont.hazard$cutpoints)-1)
        } else if (!is.null(cont.hazard$n.intervals)) {
            baseline.hazard.values = seq(1, cont.hazard$n.intervals)
        } else {
            baseline.hazard.values = NULL
        }
        
        f.hazard = paste(
                "+ f(baseline.hazard, model=\"", cont.hazard$model,"\"",
                inla.ifelse(!is.null(baseline.hazard.values),
                            inla.paste(c(", values = ", inla.2list(baseline.hazard.values))), ""),
                ", fixed = ", cont.hazard$fixed,
                inla.ifelse(is.null(cont.hazard$initial), ", initial = NULL",
                            paste(", initial = ", cont.hazard$initial, "", sep="")),
                ", constr = ", cont.hazard$constr,
                inla.ifelse(is.null(cont.hazard$prior), ", prior = NULL",
                            paste(", prior = \"", cont.hazard$prior, "\"", sep="")),
                inla.ifelse(is.null(cont.hazard$param), ", param = NULL",
                            paste(", param = ", inla.2list(cont.hazard$param), sep="")),
                ", si = ", inla.ifelse(cont.hazard$si, "TRUE", "FALSE"),
                inla.ifelse(is.null(strata.var), "", paste(", replicate=", strata.var)),
                ")", sep="")
        
        inla.eval(paste("surv.formula = .y.surv ~ ", inla.formula2character(formula[3]), f.hazard))

        if (debug) {
            print(f.hazard)
            print(inla.formula2character(surv.formula[3]))
            print(paste("surv.formula = y.surv ~ ", inla.formula2character(formula[3]), f.hazard))
        }

        result = inla(surv.formula,
                     family = "poisson",
                     data = new.data, 
                     quantiles=quantiles,
                     E = .E,
                     offset= offset,
                     scale = scale,
                     Ntrials = NULL,
                     lincomb = lincomb,
                     verbose = verbose,
                     control.compute = control.compute,
                     control.predictor = control.predictor,
                     control.data = control.data,
                     control.inla = control.inla,
                     control.results = control.results,
                     control.fixed = control.fixed,
                     control.mode = control.mode,
                     control.expert = control.expert,
                     control.lincomb = control.lincomb,
                     control.hazard = control.hazard,
                     only.hyperparam = only.hyperparam,
                     inla.call = inla.call,
                     inla.arg = inla.arg,
                     num.threads = num.threads,
                     keep = keep,
                     working.directory = working.directory,
                     silent = silent,
                     debug = debug,
                     user.hook = user.hook,
                     user.hook.arg = user.hook.arg,
                     ## internal options used to transfer data after expansions
                     .internal = .internal)

        ## replace the argument so it can be reused, if...
        result$call.orig = deparse(match.call())
        return (result)
    }

    ## this is nice hack ;-) we keep the original response. then we
    ## delete it from 'data' keeping a copy of the original one
    y...orig = eval(parse(text=formula[2]),data)
    if (!inherits(y...orig, "inla.surv"))
        class(y...orig) = NULL
    
    if (n.family > 1) {
        y...orig = inla.as.list.of.lists(y...orig)
        ny = max(sapply(y...orig, length))
        nc = length(y...orig)
        if (n.family != nc)
            stop(paste("Number of families", n.family, "does not match number of responce variables", nc))
    } else {
        if (inherits(y...orig, "inla.surv"))
            ny = max(sapply(y...orig, length))
        else
            ny = length(y...orig)
        nc = 1
    }

    data.orig = data
    data = inla.remove(as.character(formula[2]), data) 

    if (n.family == 1) {
        y...fake = c(rep(Inf, ny))
        if (debug)
            print(paste("y...fake has length", length(y...fake)))
        
    } else {
        if (FALSE) {
            ## number of NA's in one row must be >= nc-1
            n.na = apply(y...orig, 1, function(x) sum(is.na(x)))
            if (sum(!(n.na >= nc-1)) > 0)
                stop(paste("\t",sum(!(n.na >= nc-1)), " of the rows of the responce ``", formula[2],
                           "'', has more than one real responce.\n",
                           "\tOnly one (or zero) real responce is allowed for each row", sep=""))
        }
        y...fake = c(rep(Inf, ny))
        if (debug)
            print(paste("y...fake has length", length(y...fake)))
    }

    if (FALSE) {
        if (n.family != 1)
            stop(paste("length(family) are", n.family, "but the responce has only one column!"))
        
        y...fake = c(y...orig)
        y...fake[is.na(y...fake)] = Inf ## otherwise model.matrix() fails below
    }

    if (debug)
        cat("n.family", n.family, "\n")

    
    ## ...and then we add the fake to the data-frame after removing the original responce
    if (is.data.frame(data))
        data = as.data.frame(c(as.data.frame(data), list(y...fake=y...fake)))
    else
        data = c(data, list(y...fake=y...fake))

    ## creates a new version of ``formula'' and keep the old one for reference
    formula.orig = formula
    ##inla.eval(paste("formula = y...fake ~ ", inla.formula2character(formula.orig[3])))
    formula = update.formula(formula, y...fake ~ .)
    
    ## parse the formula
    gp = inla.interpret.formula(formula, data=data)
    call = deparse(match.call())
     
    if (gp$n.fix > 0) {
        inla.na.action = function(x, ...) {
            ## set fixed effects that are NA to 0
            x[is.na(x)] = 0
            return (x)
        }

        ## use y...fake also here (which is the same as in gp$fixf[2])
        ## and construct the model.matrix().
        new.fix.formula = gp$fixf
        ##inla.eval(paste("new.fix.formula = y...fake ~ ", inla.formula2character(gp$fixf[3])))
        new.fix.formula = update.formula(new.fix.formula, y...fake ~ .)
        gp$model.matrix = model.matrix(new.fix.formula, data=model.frame(new.fix.formula, data, na.action=inla.na.action))

        ## n.fix can have been changed here due to a `-1'
        gp$n.fix = dim(gp$model.matrix)[2]
    } else {
        gp$model.matrix = NULL
    }
    
    ## control what should be computed
    cont.compute = inla.set.control.compute.default()
    cont.compute[names(control.compute)] = control.compute
    if (only.hyperparam) {
        cont.compute$hyperpar = TRUE
        cont.compute$dic = cont.compute$cpo = FALSE 
    } 
    
    ## control predictor section
    cont.pred = inla.set.control.predictor.default()
    cont.pred[names(control.predictor)] = control.predictor
    if (cont.compute$cpo || cont.compute$dic) 
        cont.pred$compute=TRUE
    if (only.hyperparam) {
        cont.pred$compute = cont.pred$return.marginals = FALSE
        cont.pred$cdf = cont.pred$quantiles = NULL
    }
    
    cont.inla =inla.set.control.inla.default(family)
    cont.inla[names(control.inla)] = control.inla

    ## control fixed
    cont.fixed = inla.set.control.fixed.default()
    cont.fixed[names(control.fixed)] = control.fixed

    ## control.data
    if (n.family == 1) {
        if (!missing(control.data) && (inla.is.list.of.lists(control.data) && length(control.data) > 1))
            stop(paste("Argument control.data does not match length(family)=", n.family))
    } else {
        if (!missing(control.data) && !(inla.is.list.of.lists(control.data) && length(control.data) == n.family))
            stop(paste("Argument control.data does not match lengt)h(family)=", n.family))
    }
    
    if (missing(control.data)) {
        tt = list(list())
        for(i in 1:n.family)
            tt[[i]] = control.data
        control.data = tt
    } else if (n.family == 1) {
        if (!inla.is.list.of.lists(control.data))
            control.data = list(control.data)
    }
    
    if (length(control.data) == 1 && n.family > 1)
        stop(paste("length(control.data) = 1 while length(family) > 1."))

    ## finally, do the check. Note that the name of the argument is
    ## also used in this function, so we need to borrow the name.
    control.data.save = control.data
    for(ii in 1:n.family) {
        control.data = control.data.save[[ii]]
        inla.check.control(control.data)
    }
    control.data = control.data.save
    
    cont.data = list(list())
    for(i.family in 1:n.family) {
        cont.data[[i.family]] = inla.set.control.data.default()
        cont.data[[i.family]][names(control.data[[i.family]])] = control.data[[i.family]]
    }

    ## control results
    cont.result = inla.set.control.results.default()
    cont.result[names(control.results)] = control.results
    
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

        ## if this already exists then create one more 
        kk=1
        ans=file.exists(working.directory)
        while(ans) {
            working.directory=paste(working.directory.start,"-", kk, sep="")
            kk=kk+1
            ans=file.exists(working.directory)
        }
        inla.dir=working.directory
        inla.dir.create(inla.dir)
        cat("Model and results are stored in working directory [",inla.dir,"]\n", sep="")
    } else {
        ##create a temporary directory
        inla.dir=inla.tempfile()
        inla.dir=gsub("\\\\", "/", inla.dir)
        inla.dir.create(inla.dir)
    }
    ## Create a directory where to store data and results
    data.dir=paste(inla.dir, "/data.files",sep="")
    results.dir = paste(inla.dir, "/results.files",sep="")
    inla.dir.create(data.dir)

    ## create the .file.ini and make the problem.section
    file.ini = paste(inla.dir, "/Model.ini", sep="")

    ## problem section
    if (debug) 
        print("prepare problem section")

    inla.problem.section(file = file.ini, data.dir = data.dir, result.dir = results.dir,
                         hyperpar = cont.compute$hyperpar, return.marginals = cont.compute$return.marginals,
                         dic = cont.compute$dic, mlik = cont.compute$mlik, cpo = cont.compute$cpo,
                         quantiles = quantiles, smtp = cont.compute$smtp, q = cont.compute$q,
                         strategy = cont.compute$strategy)

    ## PREPARE RESPONSE AND FIXED EFFECTS
    if (debug)
        cat("Prepare inla file.....")

    ## copy the argument-lists
    mf = match.call(expand.dots = FALSE)
    ##mf = mf[names(mf) != "user.hook"]
    ##mf = mf[names(mf) != "user.hook.arg"]
    mf$family = NULL; mf$quantiles=NULL; 
    mf$verbose = NULL; mf$control.compute = NULL; mf$control.predictor = NULL; mf$silent = NULL; mf$control.hazard=NULL;
    mf$control.data = NULL; mf$control.inla = NULL; mf$control.results = NULL; mf$control.fixed = NULL; mf$control.lincomb=NULL;
    mf$control.mode = NULL; mf$control.expert = NULL; mf$inla.call = NULL; mf$num.threads = NULL; mf$keep = NULL;
    mf$working.directory = NULL; mf$only.hyperparam = NULL; mf$debug = NULL; 
    mf$user.hook = NULL; mf$user.hook.arg = NULL; mf$inla.arg = NULL; mf$lincomb=NULL;
    mf$.internal = NULL;
    mf$data = data

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
        rf$scale = rf$Ntrials = rf$offset = rf$E =  NULL ## these we do not need
        rf$formula = gp$randf
        rf = eval.parent(rf)
    } else {
        rf = NULL
    }
        
    if (gp$n.weights > 0) {
        wf = mf
        wf$scale = wf$Ntrials = wf$offset = wf$E =  NULL ## these we do not need
        wf$formula = gp$weightf
        wf = eval.parent(wf)
    } else {
        wf = NULL
    }

    ## I don't know why eval.parent(mf) fails here otherwise
    ## sometimes when called from inla.hyperpar()

    ## We set these here, instead of::
    ## scale = model.extract(mf, "scale")
    ## Ntrials = model.extract(mf, "Ntrials")
    ## E = model.extract(mf, "E")
    ## offset = as.vector(model.extract(mf, "offset"))

    for (nm in c("scale", "Ntrials", "offset", "E")) {
        inla.eval(paste("tmp = try(eval(mf$", nm, ", data), silent=TRUE)", sep=""))
        if (!is.null(tmp) && !inherits(tmp, "try-error")) {
            inla.eval(paste("mf$", nm, " = NULL", sep=""))
            inla.eval(paste(nm, " = tmp"))
        } else {
            inla.eval(paste("tmp = try(eval.parent(mf$", nm, "), silent=TRUE)", sep=""))
            if (!is.null(tmp) && !inherits(tmp, "try-error")) {
                inla.eval(paste("mf$", nm, " = NULL", sep=""))
                inla.eval(paste(nm, " = tmp"))
            } else {
                inla.eval(paste(nm, "= mf$", nm, " = NULL", sep=""))
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
    tot.data=length(mf[,1])
    ind=seq(0,tot.data-1)

    ## this takes care of the offset: `offset' is the argument,
    ## `offset.formula' is in the formula and `offset.sum' is their
    ## sum
    if (!is.null(gp$offset)) {
        ## there can be more offsets
        offset.formula = 0
        for(i in 1:length(gp$offset))
            offset.formula = offset.formula + as.vector(eval(parse(text=gp$offset[i]), data))
    } else {
        offset.formula = NULL
    }

    if (!is.null(offset.formula) && is.null(offset))
        offset.sum = offset.formula
    else if (is.null(offset.formula) && !is.null(offset))
        offset.sum = offset
    else if (!is.null(offset.formula) && !is.null(offset.formula)) {
        if (length(offset.formula) == length(offset)) {
            offset.sum = offset.formula + offset
        } else {
            stop("\n\tThe offset defined in formula and in argument has different length.")
        }
    } else {
        offset.sum = NULL
    }

    ## cat("offset.formula ", offset.formula, "\n")
    ## cat("offset         ", offset, "\n")
    ## cat("offset.sum     ", offset.sum, "\n")
    
    if (length(family) == 1)
        family = rep(family, n.family)
    
    ## Create a file with the response, each for each 'family'
    for(i.family in 1:n.family) {
        if (n.family == 1)
            yy = y...orig
        else
            yy = y...orig[[i.family]]
        
        file.data = inla.create.data.file(y.orig= yy, mf=mf, E=E, scale=scale, Ntrials=Ntrials, family=family[i.family],
                data.dir=data.dir, file=file.ini, debug=debug)
    
        ## add a section to the file.ini
        prop = inla.lmodel.properties(family[i.family], stop.on.error=TRUE)
        if (!is.null(cont.data[[i.family]]$param))
            if (length(cont.data[[i.family]]$param) != prop$nparameters)
                stop(paste("\n\tThe length of `param' on `control.data[[",i.family,"]]' has to be ", prop$nparameters))

        ## FIXED and no INITIAL
        if (is.null(cont.data[[i.family]]$initial) && !is.null(cont.data[[i.family]]$fixed) &&
            !all(cont.data[[i.family]]$fixed == FALSE))
            stop(paste("\n\tThe fixed hyperparameters in the likelihood needs initialisation;",
                       "use control.data[[",i.family,"]]=list(initial=c(...))"))

        if (!is.null(cont.data[[i.family]]$initial))
            if (length(cont.data[[i.family]]$initial) != prop$ntheta)
                stop(paste("\n\tThe length of `initial' in `control.data[[", i.family,"]]' has to be ", prop$ntheta))

        if (!is.null(cont.data[[i.family]]$fixed)) {
            if (length(cont.data[[i.family]]$fixed) != prop$ntheta)
                stop(paste("\n\tThe length of `fixed' in `control.data[[", i.family, "]]' has to be ", prop$ntheta))
        } else {
            if (prop$ntheta)
                cont.data[[i.family]]$fixed = rep(0, prop$ntheta)
            else
                cont.data[[i.family]]$fixed = NULL
        }

        if (!is.null(cont.data[[i.family]]$fixed))
            for(j. in 1:prop$ntheta)
                if (cont.data[[i.family]]$fixed[j.] && !is.numeric(cont.data[[i.family]]$initial[j.]))
                    stop(paste("\n\tLikelihood model ", i.family," parameter no ", j., " is fixed but has no intial value", sep=""))

        if (debug) 
            print("prepare data section")

        ##....then create the new section
        inla.data.section(file=file.ini,family=family[[i.family]],file.data=file.data,control=cont.data[[i.family]],
                          i.family=i.family)
    }

    ##create the PREDICTOR section. if necessary create a file with
    ##the offset for all likelihood
    if (debug) 
        print("prepare predictor section")

    if (!is.null(offset.sum)) {
        if (sum(is.na(offset.sum))>0) 
            stop("\n\tNo NA values allowed in the offset vector!")
        os = cbind(ind,offset.sum)
        file.offset = inla.tempfile(tmpdir=data.dir)
        file.create(file.offset)
        write(t(os),ncolumns=2,file=file.offset,append=FALSE)
        file.offset = gsub(data.dir, "$inladatadir", file.offset, fixed=TRUE)
    } else {
        file.offset = NULL
    }
    inla.predictor.section(file=file.ini, n=tot.data, predictor.spec=cont.pred, file.offset=file.offset, data.dir = data.dir)

    ##
    all.labels = character(0)
    if (!is.null(cont.pred$compute) && cont.pred$compute)
        all.labels = c(all.labels,"predictor")

    ##FIXED EFFECTS
    if (gp$n.fix >0) {
        nc = ncol(gp$model.matrix)
        labels = colnames(gp$model.matrix)
        for(i in 1:nc) {
            if (debug)
                cat("write label[", labels[i],"]\n")

            fixed.eff=cbind(ind, as.numeric(gp$model.matrix[,i]))
            ##remove lines with NA
            fixed.eff = fixed.eff[!is.na(fixed.eff[,2]),,drop=FALSE]

            file.fixed=inla.tempfile(tmpdir=data.dir)
            file.create(file.fixed)
            write(t(fixed.eff),ncolumns=ncol(fixed.eff),file=file.fixed,append=FALSE)
            file.fixed = gsub(data.dir, "$inladatadir", file.fixed, fixed=TRUE)

            if (debug)
                cat("file fixed", file.fixed,"\n")

            inla.linear.section(file=file.ini, file.fixed=file.fixed, label=labels[i],
                                results.dir=paste("fixed.effect",inla.num(i),sep=""),
                                control=cont.fixed, only.hyperparam=only.hyperparam)
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
                print(paste(r, gp$random.spec[[r]]$model, gp$random.spec[[r]]$of))
        }
            
        all.terms = c()
        for(r in 1:nr)
            all.terms = c(all.terms, gp$random.spec[[r]]$term)
        
        for(r in 1:nr) {
            n = nrep = ngroup = N = NULL
            
            if (gp$random.spec[[r]]$model != "linear" && gp$random.spec[[r]]$model != "z") {
                ##in this case we have to add a FFIELD section.........
                count.random = count.random+1
                xx=rf[,r +1]

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

                    N = inla.model.properties(gp$random.spec[[r]]$model)$aug.factor * n
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
                gp$random.spec[[r]]$n = n
                gp$random.spec[[r]]$ngroup = ngroup
                gp$random.spec[[r]]$nrep = nrep

                if (nrep > 1 || ngroup > 1) {
                    if (is.null(replicate))
                        stop("replicate is NULL")
                    if (is.null(group))
                        stop("group is NULL")
                }

                if (!is.null(gp$random.spec[[r]]$values)) {
                    ## no sort for values here, since they are given as they should be !!!!
                    location[[r]] = unique(gp$random.spec[[r]]$values)
                    cov = match(xx,location[[r]])-1 + inla.ifelse(nrep > 1 || ngroup > 1,  (replicate-1)*NG + (group-1)*N, 0)
                    cov[is.na(cov)] = -1
                    covariate[[r]] = cov
                } else {
                    location[[r]] = sort(unique(xx))
                    cov = match(xx,location[[r]])-1 + inla.ifelse(nrep > 1 || ngroup > 1,  (replicate-1)*NG + (group-1)*N, 0)
                    cov[is.na(cov)] = -1
                    covariate[[r]] = cov
                }

                ## this case is special.
                if (is.factor(xx) && is.factor(gp$random.spec[[r]]$values)) {
                    ## then we need to redefine the values to
                    ## something numerical. perhaps we need to fix
                    ## this later on, as then the plotting should take
                    ## into account the "labels" ???
                    if (debug)
                        print("redefine values and locations since both idx and values are factors")
                    gp$random.spec[[r]]$values = location[[r]]= 1:length(location[[r]])
                }
                
                if (is.null(gp$random.spec[[r]]$values))
                    gp$random.spec[[r]]$values = location[[r]]

                ##create a location and covariate file
                file.loc=inla.tempfile(tmpdir=data.dir)
                file.create(file.loc)
                write(location[[r]],ncolumns=1,file=file.loc,append=FALSE)
                file.loc = gsub(data.dir, "$inladatadir", file.loc, fixed=TRUE)
                
                file.cov=inla.tempfile(tmpdir=data.dir)
                file.create(file.cov)
                write(t(cbind(ind,covariate[[r]])),ncolumns=2,file=file.cov,append=FALSE)
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

                n.div.by = inla.model.properties(gp$random.spec[[r]]$model)$n.div.by
                if (!is.null(n.div.by)) {
                    if (!inla.divisible(n,by=n.div.by))
                        stop(paste("Model [", gp$random.spec[[r]]$model,"] require `n'", n, "divisible by", n.div.by))
                }

                if (nrep < 1 || nrep != round(nrep))
                    stop(paste("\nArgument NREP is void:", nrep))
                if (ngroup < 1 || ngroup != round(ngroup))
                    stop(paste("\nArgument NGROUP is void:", ngroup))

                ## for the CRW2/BYM (or agumented) models, we have to
                ## implement the constr=TRUE using the argument
                ## EXTRACONSTR

                if (inla.model.properties(gp$random.spec[[r]]$model)$augmented) {
                    fac = inla.model.properties(gp$random.spec[[r]]$model)$aug.factor
                    con = inla.model.properties(gp$random.spec[[r]]$model)$aug.constr

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
                                    stop(paste("Wrong dimension for the extraconstraint: ncol", ncol, "n", n))
                            
                                A = matrix(0,nrow+1,ncol)
                                e = c(gp$random.spec[[r]]$extraconstr$e, 0)
                                A[1:nrow,1:ncol] = gp$random.spec[[r]]$extraconstr$A

                                for(con.elm in con)
                                    A[nrow+1, (con.elm-1)*n + 1:n] = 1
                            }
                            gp$random.spec[[r]]$extraconstr = list(A=A,e=e)
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
                                    stop(paste("Wrong dimension for the extraconstraint: ncol", ncol, "n", n))
                            
                                A = matrix(0,nrow+length(con),ncol)
                                e = c(gp$random.spec[[r]]$extraconstr$e, rep(0,length(con)))
                                A[1:nrow,1:ncol] = gp$random.spec[[r]]$extraconstr$A

                                k = 1
                                for(con.elm in con) {
                                    A[nrow+k, (con.elm-1)*n.small + 1:n.small] = 1
                                    k = k + 1
                                }
                            }
                            gp$random.spec[[r]]$extraconstr = list(A=A,e=e)
                            gp$random.spec[[r]]$constr = FALSE
                        }
                    }
                } 

                ##print(gp$random.spec[[r]]$extraconstr$A)
                ##print(gp$random.spec[[r]]$extraconstr$e)
                ##browser()
                
                ##and in case a file for the extraconstraint
                if (!is.null(gp$random.spec[[r]]$extraconstr)) {
                    A=gp$random.spec[[r]]$extraconstr$A
                    e=gp$random.spec[[r]]$extraconstr$e
                    
                    if (ncol(A) != inla.model.properties(gp$random.spec[[r]]$model)$aug.factor*n)
                        stop(paste("\n\tncol in matrix A(extraconstraint) does not correspont to the length of f:",
                                   ncol(A),
                                   inla.model.properties(gp$random.spec[[r]]$model)$aug.factor*n))

                    file.extraconstr=inla.tempfile(tmpdir=data.dir)
                    file.create(file.extraconstr)
                    write(c(as.vector(t(A)),e),ncolumns=1,file=file.extraconstr,append=FALSE)
                    file.extraconstr = gsub(data.dir, "$inladatadir", file.extraconstr, fixed=TRUE)
                } else {
                    file.extraconstr = NULL
                }
                
                ##....also if necessary a file for the weights
                if (!is.null(gp$random.spec[[r]]$weights)) {
                    ## $weights is the name
                    www = wf[, gp$random.spec[[r]]$weights ]
                    www[is.na(www)] = 0

                    ##create a file for the weights
                    file.weights=inla.tempfile(tmpdir=data.dir)
                    file.create(file.weights)
                    write(t(cbind(ind,www)),ncolumns=2,file=file.weights,append=FALSE)
                    file.weights = gsub(data.dir, "$inladatadir", file.weights, fixed=TRUE)

                    n.weights = n.weights+1
                }
                ##create a FFIELD section
                inla.ffield.section(file=file.ini, file.loc=file.loc, file.cov=file.cov,
                                    file.extraconstr=file.extraconstr, 
                                    file.weights=file.weights, n=n, nrep = nrep, ngroup = ngroup,
                                    random.spec=gp$random.spec[[r]], 
                                    results.dir=paste("random.effect", inla.num(count.random), sep=""), 
                                    only.hyperparam= only.hyperparam,
                                    data.dir=data.dir)
            }
            else if (inla.one.of(gp$random.spec[[r]]$model, "linear")) {
                ##....while here we have to add a LINEAR section
                count.linear = count.linear+1
                xx=rf[,r +1]
                file.linear = inla.tempfile(tmpdir=data.dir)
                file.create(file.linear)
                write(t(cbind(ind,xx)),ncolumns=2,file=file.linear,append=FALSE)
                file.linear = gsub(data.dir, "$inladatadir", file.linear, fixed=TRUE)

                cont = list(cdf=gp$random.spec[[r]]$cdf,
                        quantiles=gp$random.spec[[r]]$quantiles,
                        prec=gp$random.spec[[r]]$prec.linear,
                        mean=gp$random.spec[[r]]$mean.linear)
                inla.linear.section(file=file.ini, file.fixed=file.linear, label=gp$random.spec[[r]]$term,
                                    results.dir=paste("fixed.effect",inla.num(gp$n.fix.revised+count.linear),sep=""),
                                    control=cont, only.hyperparam=only.hyperparam)
            }
            else if (inla.one.of(gp$random.spec[[r]]$model, "z")) {
                if (dim(gp$random.spec[[r]]$Z)[1] != tot.data)
                    stop(paste("\n\tNumber of data is", tot.data, "but dimension of Z is", dim(gp$random.spec[[r]]$Z)))

                inla.z.section(file=file.ini, random.spec = gp$random.spec[[r]], data.dir = data.dir,
                               results.dir = results.dir, only.hyperparam = only.hyperparam, k.off = count.random)
                count.random  = count.random + dim(gp$random.spec[[r]]$Z)[2]
            } else {
                stop("This should not happen.")
            }
        }
    }
    
    if (n.weights != gp$n.weights) 
        stop("\n\tSomething strange with weights in the covariate...")

    ## the inla section
    inla.inla.section(file=file.ini,inla.spec=cont.inla)

    ## create mode section
    cont.mode = inla.set.control.mode.default()
    cont.mode[names(control.mode)] = control.mode
    inla.mode.section(file=file.ini, cont.mode, data.dir)

    ## create expert section
    cont.expert = inla.set.control.expert.default()
    cont.expert[names(control.expert)] = control.expert
    inla.expert.section(file=file.ini, cont.expert)
    
    ## create lincomb section
    cont.lincomb = inla.set.control.lincomb.default()
    cont.lincomb[names(control.lincomb)] = control.lincomb
    inla.lincomb.section(file=file.ini, data.dir=data.dir, contr=cont.lincomb, lincomb=lincomb)
    
    if (!is.null(user.hook) && class(user.hook) == "function")
        user.hook(file.ini = file.ini, data.dir = data.dir, results.dir = results.dir,
                  formula = formula.orig, data = data.orig, arg = user.hook.arg)

    ## now, do the job
    if (debug) {
        cat("...done\n")
        cat("Run inla...")
    }

    ## inla.arg override all default arguments including `-b' !!!
    if (is.null(inla.arg)) {
        arg.arg = ""
        arg.nt = inla.ifelse(is.numeric(num.threads), paste(" -t ", num.threads, " ", sep=""), "")
        arg.v = inla.ifelse(verbose, "-v", "")
        arg.s = inla.ifelse(silent, "-s", "")
        arg.b = "-b"
    } else {
        arg.arg = inla.arg
        arg.nt = ""
        arg.v = ""
        arg.s = ""
        arg.b = ""
    }
    
    all.args = paste(arg.arg, arg.b, arg.s, arg.v, arg.nt, sep=" ")

    ## define some environment variables for remote computing
    inla.eval(paste("Sys.setenv(", "\"INLA_PATH\"", "=\"", system.file("bin", package="INLA"), "\"", ")", sep=""))
    inla.eval(paste("Sys.setenv(", "\"INLA_OS\"", "=\"", inla.os.type() , "\"", ")", sep=""))
    if (remote && inla.os("windows")) {
        inla.eval(paste("Sys.setenv(", "\"INLA_SSH_AUTH_SOCK\"", "=\"", inla.getOption("ssh.auth.sock"), "\"", ")", sep=""))
        inla.eval(paste("Sys.setenv(", "\"INLA_CYGWIN_HOME\"", "=\"", inla.getOption("cygwin.home"), "\"", ")", sep=""))
        inla.eval(paste("Sys.setenv(", "\"INLA_HOME\"", "=\"",
                        inla.cygwin.map.filename(gsub("\\\\", "/", inla.get.HOME())), "\"", ")", sep=""))
    } else {
        inla.eval(paste("Sys.setenv(", "\"INLA_HOME\"", "=\"", inla.get.HOME(), "\"", ")", sep=""))

        ## if SSH_AUTH_SOCK is not set, then we can pass it to the remote computing script
        if (Sys.getenv("SSH_AUTH_SOCK") == "") {
            inla.eval(paste("Sys.setenv(", "\"INLA_SSH_AUTH_SOCK\"", "=\"", inla.getOption("ssh.auth.sock"), "\"", ")", sep=""))
        }
    }

    my.time.used[2] = Sys.time()
    
    ## ...meaning that if inla.call = "" then just build the files (optionally...)
    if (nchar(inla.call) > 0) {
        if (inla.os("linux") || inla.os("mac")) {
            echoc = system(paste(shQuote(inla.call), all.args, shQuote(file.ini)))
        }
        else if (inla.os("windows")) {
            if (!remote) {
                echoc = try(system(paste(shQuote(inla.call), all.args, shQuote(file.ini))), silent=TRUE)
                echoc = 0
            } else {
                echoc = try(inla.cygwin.run.command(
                        paste(inla.cygwin.map.filename(inla.call),
                              all.args,
                              inla.cygwin.map.filename(file.ini))), silent=TRUE)
                echoc = 0
            }
        }
        else
            stop("\n\tNot supported architecture.")

        if (debug)
            cat("..done\n")

        my.time.used[3] = Sys.time()

        if (echoc==0) {
            ret = try(inla.collect.results(results.dir,control.results=cont.result, debug=debug,
                    only.hyperparam=only.hyperparam), silent=TRUE)
            if (!is.list(ret))
                ret = list()
            
            my.time.used[4] = Sys.time()
            cpu.used = c(
                    "Pre-processing"  = diff(my.time.used)[1],
                    "Running inla"    = diff(my.time.used)[2],
                    "Post-processing" = diff(my.time.used)[3],
                    "Total" = my.time.used[4] - my.time.used[1])

            ret$cpu.used = cpu.used
            ret$lincomb = lincomb

            ret$control.compute=cont.compute
            ret$control.predictor=cont.pred
            ret$control.lincomb = control.lincomb
            ret$control.data=cont.data
            ret$control.inla=cont.inla
            ret$control.results = cont.result
            ret$control.fixed = cont.fixed
            ret$control.mode = cont.mode
            ret$control.expert = cont.expert
            if (exists("cont.hazard"))
                ret$control.hazard = cont.hazard
            ret$control.lincomb = cont.lincomb

            ret$call=call
            ret$family=family
            ret$data=data.orig
            ret$offset=offset
            ret$Ntrials=Ntrials
            ret$E=E
            ret$scale=scale
            ret$formula=formula.orig
            ret$control.fixed=control.fixed
            ret$inla.call = inla.call
            ret$silent = silent
            ret$num.threads = num.threads
            ret$model.matrix = gp$model.matrix
            ret$user.hook = user.hook
            ret$user.hook.arg = user.hook.arg
            
            ## store also the default control-parameters; just in case
            ret$.control.defaults = list(
                    control.compute = inla.set.control.compute.default(),
                    control.predictor = inla.set.control.predictor.default(),
                    control.data = inla.set.control.data.default(),
                    control.inla = inla.set.control.inla.default(),
                    control.results = inla.set.control.results.default(),
                    control.fixed = inla.set.control.fixed.default(),
                    control.mode = inla.set.control.mode.default(),
                    control.expert = inla.set.control.expert.default(),
                    control.hazard = inla.set.control.hazard.default(),
                    control.lincomb = inla.set.control.lincomb.default())
            ##
            ret$.internal = .internal
            ##
            class(ret) = "inla"
        } else {
            ret = NULL
        }

        if (debug && !keep) cat("clean up\n")
        if (!keep) unlink(inla.dir, recursive=TRUE)
    }
    else
        ret = NULL

    ##
    return (ret)
}

`inla.check.inla.call` =
    function(inla.call = inla.getOption("inla.call"))
{
    ##
    ## Signal an error if the inla-program cannot be started
    ##

    ok = TRUE
    
    ## this is a special option: inla.call = ""
    if (is.character(inla.call) && nchar(inla.call) == 0)
        return (ok)

    if (inla.os("windows")) {
        ret = try(system(paste(shQuote(inla.call), " --ping"), intern = TRUE, ignore.stderr = TRUE,
                wait = TRUE, input = NULL, show.output.on.console = FALSE, minimized = TRUE,
                invisible = TRUE),silent = TRUE)
        if (class(ret) == "try-error")
            ok = FALSE
    }
    else if (inla.os("linux") || inla.os("mac")) {
        ret = try(system(paste(shQuote(inla.call), " --ping"), intern = TRUE, ignore.stderr = TRUE,
                wait = TRUE, input = NULL),silent = TRUE)
        ok =  (length(grep("ALIVE", ret)) > 0 || length(grep("unknown option", ret)) > 0)
    } else {
        stop("Not supported OS")
    }
    
    if (!ok)
        stop(paste("\n\n***ERROR*** Cannot find the inla-program[", inla.call, "].\n",
                   "            Please modify the [inla.call] argument\n\n", sep=""))

    return (ok)
}

`inla.user.hook` =
    function(file.ini = NULL, data.dir = NULL, results.dir = NULL, formula = NULL, data = NULL, arg = NULL)
{
    cat("\n\ninla.user.hook called with \n")

    ## This is the ini-file that is used
    cat("\tfile.ini = ", file.ini, "\n")

    ## data.dir also coded in the .ini-file as the variable inladatadir,
    ## so files in the data.dir should be specified in the .ini-file
    ## as $inladatadir/file
    cat("\tdata.dir = ", data.dir, "\n")

    ## results.dir also coded in the .ini-file as the variable inlaresdir,
    ## is where the result-files are stored. Note that the argument
    ## 'dir' in the ini-file, like
    ##
    ## dir = fixed.effect001
    ##
    ## spesify that the results are stored in inlaresdir/fixed.effect001
    ##
    cat("\tresults.dir = ", results.dir, "\n")

    ## the formula (of type 'formula') and the data.frame/list that
    ## was passed to inla().
    cat("\tformula = ", inla.formula2character(formula), "\n")
    cat("\tdata = ", names(data), "\n")

    ## arg
    cat("\targ = ")
    dput(arg)

    cat("\n")
}

