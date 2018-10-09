## Export: inla.set.control.lincomb.default
## Export: inla.set.control.update.default
## Export: inla.set.control.group.default
## Export: inla.set.control.mix.default
## Export: inla.set.control.link.default
## Export: inla.set.control.expert.default
## Export: inla.set.control.compute.default
## Export: inla.set.control.family.default
## Export: inla.set.control.fixed.default
## Export: inla.set.control.inla.default
## Export: inla.set.control.predictor.default
## Export: inla.set.control.results.default
## Export: inla.set.control.mode.default
## Export: inla.set.control.hazard.default

## Export: control.lincomb
## Export: control.update
## Export: control.group
## Export: control.mix
## Export: control.link
## Export: control.expert
## Export: control.compute
## Export: control.family
## Export: control.fixed
## Export: control.inla
## Export: control.predictor
## Export: control.results
## Export: control.mode
## Export: control.hazard



### Defines default arguments

`inla.set.control.update.default` =
    function(...)
{
    ##:EXTRA: 
    ##:NAME: control.update
    list(
        ##:ARGUMENT: result Update the joint posterior for the hyperparameters from result
        result = NULL
        )
    
    ##:SEEALSO: inla
}

`inla.set.control.lincomb.default` =
    function(...)
{
    ##:EXTRA: 
    ##:NAME: control.lincomb
    list(
        ##:ARGUMENT: precision The precision for the artificial tiny noise. Default 1e09.
        precision = 10^9,

        ##:ARGUMENT: verbose Use verbose mode for linear combinations if verbose model is set globally. (Default TRUE)
        verbose = TRUE)

    ##:SEEALSO: inla
}

`inla.set.control.group.default` =
    function(...)
{
    ##:EXTRA: 
    ##:NAME: control.group
    list(
        ##:ARGUMENT: model Group model (one of 'exchangable', 'exchangablepos', 'ar1',  'ar', 'rw1', 'rw2', 'besag', or 'iid')
        model = "exchangeable",
        
        ##:ARGUMENT: order Defines the \code{order} of the model: for model \code{ar} this defines the order p, in AR(p). Not used for other models at the time being.
        order = NULL, 

        ##:ARGUMENT: cyclic Make the group model cyclic? (Only applies to models 'ar1',  'rw1' and 'rw2')
        cyclic = FALSE,
         
        ##:ARGUMENT: graph The graph spesification (Only applies to model 'besag')
        graph = NULL, 

        ##:ARGUMENT: scale.model Scale the intrinsic model (RW1, RW2, BESAG) so the generalized variance is 1. (Default \code{TRUE})
        scale.model = TRUE, 

        ##:ARGUMENT: adjust.for.con.comp Adjust for connected components when \code{scale.model=TRUE}? (default \code{TRUE})
        adjust.for.con.comp = TRUE, 

        ##:ARGUMENT: hyper Definition of the hyperparameter(s)
        hyper = NULL,

        ##:ARGUMENT: initial (OBSOLETE!) The initial value for the group correlation or precision in the internal scale.
        initial = NULL,

        ##:ARGUMENT: fixed (OBSOLETE!) A boolean variable if the group correction or precision is assumed to be fixed or random.
        fixed = NULL,

        ##:ARGUMENT: prior (OBSOLETE!) The name of the prior distribution for the group correlation or precision in the internal scale
        prior = NULL,

        ##:ARGUMENT: param (OBSOLETE!) Prior parameters
        param = NULL)

    ##:SEEALSO: inla
}


`inla.set.control.mix.default` =
    function(...)
{
    ##:EXTRA: The \code{control.mix} -list is set within the corresponding \code{control.family}-list a the mixture of the likelihood is likelihood spesific. (This option is EXPERIMENTAL.)
    ##:NAME: control.mix
    list(
        ##:ARGUMENT: model The model for the random effect. Currently, only \code{model='gaussian'} is implemented
        model = NULL, 

        ##:ARGUMENT: hyper Definition of the hyperparameter(s) for the random effect model chosen
        hyper = NULL,

        ##:ARGUMENT: initial (OBSOLETE!) The initial value(s) for the hyperparameter(s)
        initial = NULL,

        ##:ARGUMENT: fixed (OBSOLETE!) A boolean variable if hyperparmater(s) is/are fixed or random
        fixed = NULL,

        ##:ARGUMENT: prior (OBSOLETE!) The name of the prior distribution(s) for the hyperparmater(s)
        prior = NULL,

        ##:ARGUMENT: param (OBSOLETE!) The parameters for the prior distribution(s) for the hyperparmater(s)
        param = NULL)

    ##:SEEALSO: inla
}

`inla.set.control.link.default` =
    function(...)
{
    ##:EXTRA: The \code{control.link}-list is set within the corresponding \code{control.family}-list as the link is likelihood-familiy spesific.
    ##:NAME: control.link
    list(
        ##:ARGUMENT: model The name of the link function/model
        model = "default",

        ##:ARGUMENT: order The \code{order} of the link function, where the interpretation of \code{order} is model-dependent.
        order = NULL, 

        ##:ARGUMENT: variant The \code{variant} of the link function, where the interpretation of \code{variant} is model-dependent.
        variant = NULL, 

        ##:ARGUMENT: nq Number of quadrature-points used to do the numerical integration (default 15)
        nq = 15, 

        ##:ARGUMENT: hyper Definition of the hyperparameter(s) for the link model chosen
        hyper = NULL,

        ##:ARGUMENT: quantile The quantile for quantile link function
        quantile = NULL,

        ##:ARGUMENT: initial (OBSOLETE!) The initial value(s) for the hyperparameter(s)
        initial = NULL,

        ##:ARGUMENT: fixed (OBSOLETE!) A boolean variable if hyperparmater(s) is/are fixed or random
        fixed = NULL,

        ##:ARGUMENT: prior (OBSOLETE!) The name of the prior distribution(s) for the hyperparmater(s)
        prior = NULL,

        ##:ARGUMENT: param (OBSOLETE!) The parameters for the prior distribution(s) for the hyperparmater(s)
        param = NULL)

    ##:SEEALSO: inla
}


`inla.set.f.default` =
    function(...)
{
    list(diagonal = .Machine$double.eps^0.319) ## almost 1e-5 on my computer
}


`inla.set.control.expert.default` =
    function(...)
{
    ##:EXTRA: 
    ##:NAME: control.expert
    list(
        ##:ARGUMENT: cpo.manual A boolean variable to decide if the inla-program is to be runned in a manual-cpo-mode. (EXPERT OPTION: DO NOT USE)
        cpo.manual = FALSE,

        ##:ARGUMENT: cpo.idx  The index/indices of the data point(s) to remove. (EXPERT OPTION: DO NOT USE)
        cpo.idx = -1,

        ##:ARGUMENT: disable.gaussian.check Disable the check for fast computations with a Gaussian likelihood and identity link (default \code{FALSE})
        disable.gaussian.check = FALSE, 

        ##:ARGUMENT: jp An object of class \code{inla.jp} defining a joint prior
        jp = NULL
    )

    ##:SEEALSO: inla
}


`inla.set.control.compute.default`=
    function(...)
{
    ##:EXTRA: 
    ##:NAME: control.compute
    list(
        ##:ARGUMENT: openmp.strategy The computational strategy to use: 'small', 'medium', 'large', 'huge' and 'default'. There are also two options for the pardiso solver: 'pardiso.serial' and 'pardiso.parallel'. The difference is how the parallelisation is done, and is tuned for 'small'-sized models, 'medium'-sized models, etc. The default option tries to make an educated guess, but this allows to overide this selection. Default is 'default'
        openmp.strategy = "default", 

        ##:ARGUMENT: hyperpar A boolean variable if the marginal for the hyperparameters should be computed. Default TRUE.
        hyperpar=TRUE,

        ##:ARGUMENT: return.marginals A boolean variable if the marginals for the latent field should be returned (although it is computed). Default TRUE
        return.marginals=TRUE,

        ##:ARGUMENT: dic A boolean variable if the DIC-value should be computed. Default FALSE.
        dic=FALSE,

        ##:ARGUMENT: mlik A boolean variable if the marginal likelihood should be computed. Default \code{TRUE}.
        mlik=TRUE,

        ##:ARGUMENT: cpo A boolean variable if the cross-validated predictive measures (cpo, pit) should be computed (default \code{FALSE})
        cpo=FALSE,

        ##:ARGUMENT: po A boolean variable if the predictive ordinate should be computed (default \code{FALSE})
        po=FALSE,
        
        ##:ARGUMENT: waic A boolean variable if the Watanabe-Akaike information criteria should be computed (default \code{FALSE})
        waic=FALSE,
        
        ##:ARGUMENT: q A boolean variable if binary images of the precision matrix, the reordered precision matrix and the Cholesky triangle should be generated. (Default FALSE.)
        q=FALSE,

        ##:ARGUMENT: config A boolean variable if the internal GMRF approximations be stored. (Default FALSE. EXPERIMENTAL)
        config=FALSE,

        ##:ARGUMENT: smtp The sparse-matrix solver, one of 'default', 'taucs', 'band' or 'pardiso' (default \code{inla.getoption("smtp")})
        smtp = NULL,

        ##:ARGUMENT: graph A boolean variable if the graph itself should be returned. (Default FALSE.)
        graph = FALSE, 
        
        ##:ARGUMENT: gdensity A boolean variable if the Gaussian-densities itself should be returned. (Default FALSE.)
        gdensity = FALSE)
        
    ##:SEEALSO: inla
}

`inla.set.control.family.default`=
    function(...)
{
    ##:EXTRA: 
    ##:NAME: control.family
    list(
        ##:ARGUMENT: hyper Definition of the hyperparameters
        hyper = NULL,

        ##:ARGUMENT: initial (OBSOLETE!) Initial value for the hyperparameter(s) of the likelihood in the internal scale.
        initial=NULL,

        ##:ARGUMENT: prior (OBSOLETE!) The name of the prior distribution(s) for othe hyperparameter(s).
        prior=NULL,

        ##:ARGUMENT: param (OBSOLETE!) The parameters for the prior distribution
        param=NULL,

        ##:ARGUMENT: fixed (OBSOLETE!) Boolean variable(s) to say if the hyperparameter(s) is fixed or random.
        fixed=NULL,

        ##:ARGUMENT: link (OBSOLETE! Use \code{control.link=list(model=)} instead.) The link function to use.
        link= "default",

        ##:ARGUMENT: sn.shape.max Maximum value for the shape-parameter for Skew Normal observations (default 5.0)
        sn.shape.max = 5.0,

        ##:ARGUMENT: gev.scale.xi The internal scaling of the shape-parameter for the GEV distribution. (default 0.01)
        gev.scale.xi = 0.01,

        ##:ARGUMENT: cenpoisson.I The censoring interval for the censored Poisson
        cenpoisson.I = c(-1L, -1L),

        ##:ARGUMENT: variant This variable is used to give options for various variants of the likelihood,  like chosing different parameterisations for example. See the relevant likelihood documentations for options (does only apply to some likelihoods).
        variant = 0L,

        ##:ARGUMENT: control.mix See \code{?control.mix}
        control.mix = NULL, 

        ##:ARGUMENT: control.link See \code{?control.link}
        control.link = NULL
        )

    ##:SEEALSO: inla
}

`inla.set.control.fixed.default`=
    function(...)
{
    ##:EXTRA: 
    ##:NAME: control.fixed
    list(
        ##:ARGUMENT: cdf  A list of values to compute the CDF for, for all fixed effects
        cdf=NULL,

        ##:ARGUMENT: quantiles  A list of quantiles to compute for all fixed effects
        quantiles = NULL,

        ##:ARGUMENT: expand.factor.strategy The strategy used to expand factors into fixed effects based on their levels. The default strategy is us use the \code{model.matrix}-function for which NA's are not allowed (\code{expand.factor.strategy="model.matrix"}) and levels are possible removed. The alternative option (\code{expand.factor.strategy="inla"}) use an \code{inla}-spesific expansion which expand a factor into one fixed effects for each level, do allow for NA's and all levels are present in the model. In this case, factors MUST BE factors in the data.frame/list and NOT added as \code{.+factor(x1)+.} in the formula only.
        expand.factor.strategy = "model.matrix", 

        ##:ARGUMENT: mean Prior mean for all fixed effects except the intercept. Alternatively, a named list with specific means where name=default applies to unmatched names. For example \code{control.fixed=list(mean=list(a=1, b=2, default=0))} assign 'mean=1' to fixed effect 'a' , 'mean=2' to effect 'b' and 'mean=0' to all others. (default 0.0)
        mean = 0.0,

        ##:ARGUMENT: mean.intercept Prior mean for the intercept (default 0.0)
        mean.intercept = 0.0,

        ##:ARGUMENT: prec  Default precision for all fixed effects except the intercept. Alternatively, a named list with specific means where name=default applies to unmatched names.  For example \code{control.fixed=list(prec=list(a=1, b=2, default=0.01))} assign 'prec=1' to fixed effect 'a' , 'prec=2' to effect 'b' and 'prec=0.01' to all others. (default 0.001)
        prec= 0.001,

        ##:ARGUMENT: prec.intercept  Default precision the intercept (default 0.0)
        prec.intercept = 0.0,

        ##:ARGUMENT: compute Compute marginals for the fixed effects ? (default TRUE)
        compute = TRUE,

        ##:ARGUMENT: correlation.matrix Compute the posterior correlation matrix for all fixed effects? (default FALSE) OOPS: This option will set up appropriate linear combinations and the results are shown as the posterior correlation matrix of the linear combinations. This option will imply \code{control.inla=list(lincomb.derived.correlation.matrix=TRUE)}.
         
        correlation.matrix = FALSE)

    ##:SEEALSO: inla
}

`inla.set.control.inla.default`=
    function(...)
{
    family = "gaussian"
    xx = list(...)[1]
    if (!is.null(xx$family)) {
        family = xx$family
    }

    ##:EXTRA: 
    ##:NAME: control.inla
    ans = list(
        ##:ARGUMENT: strategy  Character The strategy to use for the approximations; one of 'gaussian', 'simplified.laplace' (default), 'laplace' or 'adaptive'
        strategy="simplified.laplace",

        ##:ARGUMENT: int.strategy  Character The integration strategy to use; one of 'auto' (default),  'ccd', 'grid', 'eb' (empirical bayes),  'user' or 'user.std'
        int.strategy="auto",

        ##:ARGUMENT: int.design  Matrix Matrix of user-defined integration points and weights. Each row consists theta values and the integration weight. (EXPERIMENTAL!)
        int.design=NULL,

        ##:ARGUMENT: interpolator  Character The interpolator used to compute the marginals for the hyperparameters. One of 'auto', 'nearest', 'quadratic', 'weighted.distance', 'ccd', 'ccdintegrate', 'gridsum', 'gaussian'. Default is 'auto'.
        interpolator="auto",

        ##:ARGUMENT: fast Logical If TRUE, then replace conditional modes in the Laplace approximation with conditional expectation (default TRUE)
        fast = TRUE,

        ##:ARGUMENT: linear.correction Logical Default TRUE for the 'strategy = laplace' option.
        linear.correction=NULL,

        ##:ARGUMENT: h Numerical The step-length for the gradient calculations for the hyperparameters. Default 0.01.
        h=0.01,

        ##:ARGUMENT: dz Numerical The step-length in the standarised scale for the integration of the hyperparameters. Default 0.75.
        dz=0.75,

        ##:ARGUMENT: diff.logdens Numerical The difference of the log.density for the hyperpameters to stop numerical integration using int.strategy='grid'. Default 6.
        diff.logdens=6,

        ##:ARGUMENT: print.joint.hyper Logical If TRUE, the store also the joint distribution of the hyperparameters (without any costs). Default TRUE.
        print.joint.hyper=TRUE,

        ##:ARGUMENT: force.diagonal Logical If TRUE, then force the Hessian to be diagonal. (Default \code{FALSE})
        force.diagonal=FALSE,

        ##:ARGUMENT: skip.configurations Logical Skip configurations if the values at the main axis are to small. (Default \code{TRUE})
        skip.configurations=TRUE,

        ##:ARGUMENT: mode.known Logical If TRUE then no optimisation is done. (Default FALSE.)
        mode.known=FALSE,

        ##:ARGUMENT: adjust.weights Logical If TRUE then just more accurate integration weights. (Default TRUE.)
        adjust.weights=TRUE,

        ##:ARGUMENT: tolerance Numerical The tolerance for the optimisation of the hyperparameters. If set, this is the default value for for 'tolerance.f^(2/3)',  'tolerance.g' and  'tolerance.x'; see below.
        tolerance = 0.005,

        ##:ARGUMENT: tolerance.f Numerical The tolerance for the absolute change in the log posterior in the optimisation of the hyperparameters.
        tolerance.f = NULL,

        ##:ARGUMENT: tolerance.g Numerical The tolerance for the absolute change in the gradient of the log posterior in the optimisation of the hyperparameters.
        tolerance.g = NULL,

        ##:ARGUMENT: tolerance.x Numerical The tolerance for the change in the hyperparameters (root-mean-square) in the optimisation of the hyperparameters.
        tolerance.x = NULL, 

        ##:ARGUMENT: tolerance.step Numerical The tolerance for the change in root-mean_squre in the inner Newton-like optimisation of the latent field.
        tolerance.step = 0.0005, 

        ##:ARGUMENT: restart Numerical To improve the optimisation, the optimiser is restarted at the found optimum 'restart' number of times.
        restart = 0L,

        ##:ARGUMENT: optimiser Character The optimiser to use; one of 'gsl' or 'default'.
        optimiser = "default",

        ##:ARGUMENT: verbose Logical Run in verbose mode? (Default FALSE)
        verbose = NULL,

        ##:ARGUMENT: reordering Character Type of reordering to use. (EXPERT OPTION; one of "AUTO", "DEFAULT", "IDENTITY", "REVERSEIDENTITY",  "BAND", "METIS", "GENMMD", "AMD", "MD", "MMD", "AMDBAR", "AMDC", "AMDBARC",  or the output from \code{inla.qreordering}. Default is 'auto'.)
        reordering = "auto",

        ##:ARGUMENT: cpo.diff Numerical Threshold to define when the cpo-calculations are inaccurate. (EXPERT OPTION.)
        cpo.diff = NULL,

        ##:ARGUMENT: npoints Numerical Number of points to use in the 'stratey=laplace' approximation (default 9)
        npoints = 9,

        ##:ARGUMENT: cutoff Numerical The cutoff used in the 'stratey=laplace' approximation. (Smaller value is more accurate and more slow.) (default 1e-4)
        cutoff = 1e-4,

        ##:ARGUMENT: adapt.hessian.mode Logical Should optimisation be continued if the Hessian estimate is void? (Default TRUE)
        adapt.hessian.mode = NULL,

        ##:ARGUMENT: adapt.hessian.max.trials Numerical Number of steps in the adaptive Hessian optimisation
        adapt.hessian.max.trials = NULL,

        ##:ARGUMENT: adapt.hessian.scale Numerical The scaling of the 'h' after each trial.
        adapt.hessian.scale = NULL, 

        ##:ARGUMENT: adaptive.max Selecting \code{strategy="adaptive"} will chose the default strategy for all fixed effects and model components with length less or equal to \code{adaptive.max}, for others, the gaussian strategy will be applied.
        adaptive.max = 10L,
        
        ##:ARGUMENT: huge Logical If TRUE then try to do some of the internal parallisations differently. Hopefully this will be of benefite for 'HUGE' models. (Default FALSE.) [THIS OPTION IS OBSOLETE AND NOT USED!]
        huge = FALSE,

        ##:ARGUMENT: step.len Numerical The step-length used to compute numerical derivaties of the log-likelihood
        step.len = .Machine$double.eps^(1.0/3.9134),

        ##:ARGUMENT: stencil Numerical Number of points in the stencil used to compute the numerical derivaties of the log-likelihood (3, 5, 7 or 9). (default 5)
        stencil = 5L, 

        ##:ARGUMENT: lincomb.derived.only Logical If TRUE the only compute the marginals for the derived linear combinations and if FALSE, the and also the linear combinations to the graph (Default TRUE)
        lincomb.derived.only = TRUE,

        ##:ARGUMENT: lincomb.derived.correlation.matrix Logical If TRUE compute also the correlations for the derived linear combinations, if FALSE do not (Default FALSE)
        lincomb.derived.correlation.matrix = FALSE,

        ##:ARGUMENT: diagonal Numerical Expert use only! Add a this value on the diagonal of the joint precision matrix. (default 0.0)
        diagonal = 0.0,

        ##:ARGUMENT: numint.maxfeval Numerical Maximum number of function evaluations in the the numerical integration for the hyperparameters. (Default 100000.)
        numint.maxfeval = 100000,

        ##:ARGUMENT: numint.relerr Numerical Relative error requirement in the the numerical integration for the hyperparameters. (Default 1e-5)
        numint.relerr = 1e-5,

        ##:ARGUMENT: numint.abserr Numerical Absolute error requirement in the the numerical integration for the hyperparameters. (Default 1e-6)
        numint.abserr = 1e-6,

        ##:ARGUMENT: cmin Numerical The minimum value for the negative Hessian from the likelihood. Increasing this value will stabalise the optimisation but can introduce bias in some estimates unless \code{-Inf} is used.  (Default -Inf)
        cmin = -Inf,

        ##:ARGUMENT: step.factor Numerical The step factor in the Newton-Raphson algorithm saying how large step to take (Default 1.0)
        ## YES! setting this to a negative values means = 1,  EXCEPT the first time (for each thread) where |step.factor| is used.
        ## This is an hidden option.
        step.factor = -0.1,

        ##:ARGUMENT: global.node.factor Numerical The factor which defines the degree required (how many neighbors), as a fraction of n-1, that is required to be classified as a global node and numbered last (whatever the reordering routine says). Here,  n,  is the size of the graph. (Disabled if larger than 1.) (default 2.0)
        global.node.factor = 2.0, 

        ##:ARGUMENT: global.node.degree Numerical The degree required (number of neighbors) to be classified as a global node and numbered last (whatever the reordering routine says). (default \code{.Machine$integer.max})
        global.node.degree = .Machine$integer.max,

        ##:ARGUMENT: stupid.search Logical Enable or disable the stupid-search-algorithm, if the Hessian calculations reveals that the mode is not found. (Default \code{TRUE}.)
        stupid.search = TRUE,
            
        ##:ARGUMENT: stupid.search.max.iter Numerical Maximum number of iterations allowed for the stupid-search-algorithm. (default 1000)
        stupid.search.max.iter = 1000L, 

        ##:ARGUMENT: stupid.search.factor Numerical Factor (>=1) to increase the step-length with after each new interation. (default 1.05)
        stupid.search.factor = 1.05,
        
        ##:ARGUMENT: correct Logical Add correction for the Laplace approximation. (default \code{FALSE})
        correct = FALSE,

        ##:ARGUMENT: correct.factor Numerical Factor used in adjusting the correction factor (default=10) if correct=TRUE
        correct.factor = 10.0,

        ##:ARGUMENT: correct.strategy  Character The strategy used to compute the correction; one of 'simplified.laplace' (default) or 'laplace'
        correct.strategy = "simplified.laplace", 

        ##:ARGUMENT: correct.verbose  Logical Be verbose when computing the correction? (default \code{FALSE})
        correct.verbose = FALSE)
 
    ##:SEEALSO: inla

    return (ans)
}


`inla.set.control.predictor.default`=
    function(...)
{
    ##:EXTRA: 
    ##:NAME: control.predictor
    list(
        ##:ARGUMENT: hyper Definition of the hyperparameters.
        hyper = NULL,

        ##:ARGUMENT: fixed (OBSOLETE!) If the precision for the artificial noise is fixed or not (defualt TRUE)
        fixed=NULL,

        ##:ARGUMENT: prior (OBSOLETE!) The prior for the artificial noise
        prior=NULL,

        ##:ARGUMENT: param (OBSOLETE!) Prior parameters for the artificial noise
        param=NULL,

        ##:ARGUMENT: initial (OBSOLETE!) The value of the log precision of the artificial noise
        initial=NULL,

        ##:ARGUMENT: compute A boolean variable; should the marginals for the linear predictor be computed? (Default FALSE.)
        compute=FALSE,

        ##:ARGUMENT: cdf A list of values to compute the CDF for the linear predictor
        cdf=NULL,

        ##:ARGUMENT: quantiles A list of quantiles to compute for the linear predictor
        quantiles = NULL,

        ##:ARGUMENT: cross Cross-sum-to-zero constraints with the linear predictor. All linear predictors with the same level of 'cross' are constrained to have sum zero. Use 'NA' for no contribution. 'Cross' has the same length as the linear predictor (including the 'A' matrix extention). (THIS IS AN EXPERIMENTAL OPTION, CHANGES MAY APPEAR.)
        cross=NULL,

        ##:ARGUMENT: A The observation matrix (matrix or Matrix::sparseMatrix).
        A = NULL,

        ##:ARGUMENT: precision The precision for eta* - A*eta, (default \code{exp(15)})
        precision = exp(15),

        ##:ARGUMENT: link Define the family-connection for unobserved observations (\code{NA}). \code{link} is integer values which defines the family connection; \code{family[link[idx]]} unless \code{is.na(link[idx])} for which the identity-link is used. The \code{link}-argument only influence the \code{fitted.values} in the \code{result}-object. If \code{is.null(link)} (default) then the identity-link is used for all missing observations. If the length of \code{link} is 1, then this value is replicated with the length of the responce vector. If an element of the responce vector is \code{!NA} then the corresponding entry in \code{link} is not used (but must still be a legal value). Setting this variable implies \code{compute=TRUE}.
        link = NULL)

    ##:SEEALSO: inla
}

`inla.set.control.results.default`=
    function(...)
{
    ##:EXTRA: 
    ##:NAME: control.results
    list(
        ##:ARGUMENT: return.marginals.random A boolean variable; read the marginals for the fterms? (Default TRUE)
        return.marginals.random=TRUE,

        ##:ARGUMENT: return.marginals.predictor A boolean variable; read the marginals for the linear predictor? (Default TRUE)
        return.marginals.predictor=TRUE)
    ##:SEEALSO: inla
}

`inla.set.control.mode.default`=
    function(...)
{
    ## this is internal use only...
    ##:EXTRA: 
    ##:NAME: control.mode
    list(
        ##:ARGUMENT: result Prevous result from inla(). Use the theta- and x-mode from this run.
        result = NULL,

        ##:ARGUMENT: theta The theta-mode/initial values for theta. This option has preference over result$mode$theta.
        theta = NULL,

        ##:ARGUMENT: x The x-mode/intitial values for x. This option has preference over result$mode$x.
        x = NULL,

        ##:ARGUMENT: restart A boolean variable; should we restart the optimisation from this configuration or fix the mode at this configuration? (Default FALSE.)
        restart = FALSE,

        ##:ARGUMENT: fixed A boolean variable. If TRUE then treat all thetas as known and fixed, and if FALSE then treat all thetas as unknown and random (default).
        fixed = FALSE)
    ##:SEEALSO: inla
}

`inla.set.control.hazard.default` =
    function(...)
{
    ##:EXTRA:
    ##:NAME: control.hazard
    list(
        ##:ARGUMENT: model The model for the baseline hazard model. One of 'rw1' or 'rw2'. (Default 'rw1'.)
        model = "rw1",

        ##:ARGUMENT: hyper The definition of the hyperparameters.
        hyper = NULL,

        ##:ARGUMENT: fixed (OBSOLETE!) A boolean variable; is the precision for 'model' fixed? (Default FALSE.)
        fixed = FALSE,

        ##:ARGUMENT: initial (OBSOLETE!) The initial value for the precision.
        initial = NULL,

        ##:ARGUMENT: prior (OBSOLETE!) The prior distribution for the precision for 'model'
        prior = NULL,

        ##:ARGUMENT: param (OBSOLETE!) The parameters in the prior distribution
        param = NULL,

        ##:ARGUMENT: constr A boolean variable; shall the  'model' be constrained to sum to zero?
        constr = TRUE,

        ##:ARGUMENT: diagonal An extra constant added to the diagonal of the precision matrix
        diagonal = NULL,

        ##:ARGUMENT: n.intervals Number of intervals in the baseline hazard. (Default 15)
        n.intervals = 15,

        ##:ARGUMENT: cutpoints The cutpoints to use. If not specified the they are compute from 'n.intervals' and the maximum length of the interval. (Default NULL)
        cutpoints = NULL,

        ##:ARGUMENT: strata.name The name of the stratefication variable for the baseline hazard in the data.frame
        strata.name = NULL,

        ##:ARGUMENT: scale.model Scale the baseline hazard model (RW1, RW2) so the generalized variance is 1. (Default \code{inla.getOption("scale.model.default")}.)
        scale.model = NULL)
    ##:SEEALSO: inla
}

## check control-arguments

`inla.check.control` = function(contr, data = NULL)
{
    ## This function will signal an error if the arguments in CONTR
    ## does not match the ones in the corresponding
    ## `inla.set.XX.default()' routine.  EXAMPLE: contr is
    ## `control.inla' and default arguments is found in
    ## `inla.set.control.inla.default()'
    
    ## Will expand unexpanded names from the names in 'data' first
    contr = local({
        name = paste("inla.tmp.env", as.character(runif(1)), sep="")
        attach(data, name = name, warn.conflicts = FALSE)
        ccontr = contr
        detach(name, character.only = TRUE)
        ccontr
    })

    stopifnot(!missing(contr))
    stopifnot(is.list(contr))
    if (length(contr) == 0) {
        return(contr)
    }

    nm = paste(sys.call()[2])
    f = paste("inla.set.", nm, ".default()", sep="")
    elms = names(inla.eval(f))

    if (is.null(names(contr))) {
        stop(inla.paste(c("Named elements in in control-argument `", nm, "', is required: ", contr,
                          "\n\n  Valid ones are:\n\t",
                          inla.paste(sort(elms), sep="\n\t")), sep=""))
    }
        
    for(elm in names(contr)) {
        if (!is.element(elm, elms)) {
            stop(inla.paste(c("Name `", elm,"' in control-argument `", nm, "', is void.\n\n  Valid ones are:\n\t",
                              inla.paste(sort(elms), sep="\n\t")), sep=""))
        }
    }

    return(contr)
}



## test-implementation
##`control.lincomb` = function(precision, verbose)
##{
##    aa = match.call()[-1]
##    ret = list()
##    for(a in names(aa)) {
##        if (!missing(a)) {
##            xx = get(a)
##            names(xx) = a
##            ret = c(ret, xx)
##        }
##    }
##    return (ret)
##}

inla.make.completion.function = function(...)
{
    my.eval = function(command, envir = parent.frame(),
        enclos = if (is.list(envir) || is.pairlist(envir)) parent.frame() else baseenv()) 
    {
        return(eval(parse(text = command), envir, enclos))
    }

    xx = sort(list(...)[[1L]])
    my.eval(paste("function(", paste(xx, sep="", collapse=", "), ") {
    aa = match.call()[-1L]
    ret = list()
    for(a in names(aa)) {
        if (!missing(a)) {
            xx = get(a)
            names(xx) = a
            ret = c(ret, xx)
        }
    }
    return (ret)
}"))
}

control.update = inla.make.completion.function(names(inla.set.control.update.default()))
control.lincomb = inla.make.completion.function(names(inla.set.control.lincomb.default()))
control.group = inla.make.completion.function(names(inla.set.control.group.default()))
control.mix = inla.make.completion.function(names(inla.set.control.mix.default()))
control.link = inla.make.completion.function(names(inla.set.control.link.default()))
control.expert = inla.make.completion.function(names(inla.set.control.expert.default()))
control.compute = inla.make.completion.function(names(inla.set.control.compute.default()))
control.family = inla.make.completion.function(names(inla.set.control.family.default()))
control.fixed = inla.make.completion.function(names(inla.set.control.fixed.default()))
control.inla = inla.make.completion.function(names(inla.set.control.inla.default()))
control.predictor = inla.make.completion.function(names(inla.set.control.predictor.default()))
control.results = inla.make.completion.function(names(inla.set.control.results.default()))
control.mode = inla.make.completion.function(names(inla.set.control.mode.default()))
control.hazard = inla.make.completion.function(names(inla.set.control.hazard.default()))
