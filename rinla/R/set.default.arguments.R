### Defines default arguments

## The ... documentation, seealso, and description sections are inherited
## by all other control.* method documentation.

#' @title control.update
#' 
#' @description
#' Control variables in `control.*` for use with [inla()].
#' The functions can be used to TAB-complete arguments, and
#' returns a list of the default control arguments, unless overridden by
#' specific input arguments.
#' 
#' @param ... Named arguments passed on to the main function
#' @seealso [inla()]
#' @family control
#' @export
`control.update` <-
    function(
        #' @param result Update the joint posterior for the hyperparameters from result
        result = NULL
    ) {
        as.list(environment())
    }

#' @title control.lincomb
#' @inherit control.update params description seealso 
#' @family control
#' @export
`control.lincomb` <-
    function(
        #' @param verbose Use verbose mode for linear combinations if verbose model is
        #' set globally. (Default FALSE). This option is only available for
        #' the default `inla.mode` (`inla.mode="compact"`).
        verbose = FALSE
    ) {
        as.list(environment())
    }

    
#' @title control.group
#' @inherit control.update params description seealso 
#' @family control
#' @export
`control.group` <-
    function(
        #' @param model Group model (one of 'exchangable', 'exchangablepos', 'ar1',
        #' 'ar', 'rw1', 'rw2', 'besag', or 'iid')
        model = "exchangeable",
        
        #' @param order Defines the `order` of the model: for model `ar` this
        #' defines the order p, in AR(p). Not used for other models at the time
        #' being.
        order = NULL,
        
        #' @param cyclic Make the group model cyclic? (Only applies to models 'ar1',
        #' 'rw1' and 'rw2')
        cyclic = FALSE,
        
        #' @param graph The graph specification (Only applies to model 'besag')
        graph = NULL,
        
        #' @param scale.model Scale the intrinsic model (RW1, RW2, BESAG) so the
        #' generalized variance is 1. (Default `TRUE`)
        scale.model = TRUE,
        
        #' @param adjust.for.con.comp Adjust for connected components when
        #' `scale.model=TRUE`? (default `TRUE`)
        adjust.for.con.comp = TRUE,
        
        #' @param hyper Definition of the hyperparameter(s)
        hyper = NULL,
        
        #' @param initial (OBSOLETE!) The initial value for the group correlation or
        #' precision in the internal scale.
        initial = NULL,
        
        #' @param fixed (OBSOLETE!) A boolean variable if the group correction or
        #' precision is assumed to be fixed or random.
        fixed = NULL,
        
        #' @param prior (OBSOLETE!) The name of the prior distribution for the group
        #' correlation or precision in the internal scale
        prior = NULL,
        
        #' @param param (OBSOLETE!) Prior parameters
        param = NULL
    ) {
        as.list(environment())
    }

    
#' @title control.scopy
#' @inherit control.update params description seealso 
#' @family control
#' @export
`control.scopy` <-
    function(
            #' @param covariate The covariate for the scopy function
            covariate = NULL,

            #' @param n Number of betas
            n = 5,

            #' @param model scopy model (one of 'rw1' and 'rw2')
            model = "rw2",

            #' @param mean The prior mean for mean(betas)
            mean = 1.0,

            #' @param prec.mean The prior precision for mean(betas)
            prec.mean =  1.0, 

            #' @param prec.betas The prior precision prec(betas-mean(betas))
            prec.betas =  10.0
    ) {
        as.list(environment())
    }

#' @title control.mix
#' @inherit control.update params description seealso 
#' @details The `control.mix` list is set within the corresponding `control.family`-list a the mixture of the likelihood is likelihood specific. (This option is EXPERIMENTAL.)
#' @family control
#' @export
`control.mix` <-
    function(
            #' @param model The model for the random effect. Currently, only
            #' `model='gaussian'` is implemented
            model = NULL,

            #' @param hyper Definition of the hyperparameter(s) for the random effect model
            #' chosen
            hyper = NULL,

            #' @param initial (OBSOLETE!) The initial value(s) for the hyperparameter(s)
            initial = NULL,

            #' @param fixed (OBSOLETE!) A boolean variable if hyperparmater(s) is/are fixed
            #' or random
            fixed = NULL,

            #' @param prior (OBSOLETE!) The name of the prior distribution(s) for the
            #' hyperparmater(s)
            prior = NULL,

            #' @param param (OBSOLETE!) The parameters for the prior distribution(s) for the
            #' hyperparmater(s)
            param = NULL,

            #' @param npoints Number of points used to do the numerical integration (default
            #' 101)
            npoints = 101,

            #' @param integrator The integration scheme to use (`default`,
            #' `quadrature`, `simpson`)
            integrator = "default"
    ) {
        as.list(environment())
    }

#' @title control.pom
#' @inherit control.update params description seealso 
## @details The \code{control.pom} -list is set within the corresponding \code{control.family}
#' @family control
#' @export
`control.pom` <-
    function(
            #' @param cdf character The cdf to use, "logit" (default) or "probit"
            cdf  = "logit",

            #' @param fast Logical Use a faster but approximate form for the probit cdf
            #' (default `FALSE`)?
            fast = FALSE
        ) {
            as.list(environment())
        }

#' @title control.link
#' @inherit control.update params description seealso 
#' @family control
#' @export
`control.link` <-
    function(
        #' @details The `control.link`-list is set within the corresponding `control.family`-list as the link is likelihood-family specific.
            #' @param model The name of the link function/model
            model = "default",

            #' @param order The `order` of the link function, where the interpretation
            #' of `order` is model-dependent.
            order = NULL,

            #' @param variant The `variant` of the link function, where the
            #' interpretation of `variant` is model-dependent.
            variant = NULL,

            #' @param hyper Definition of the hyperparameter(s) for the link model chosen
            hyper = NULL,

            #' @param quantile The quantile for quantile link function
            quantile = NULL,

            #' @param a The parameter `a` in the LOGa link
            a = 1.0,

            #' @param initial (OBSOLETE!) The initial value(s) for the hyperparameter(s)
            initial = NULL,

            #' @param fixed (OBSOLETE!) A boolean variable if hyperparmater(s) is/are fixed
            #' or random
            fixed = NULL,

            #' @param prior (OBSOLETE!) The name of the prior distribution(s) for the
            #' hyperparmater(s)
            prior = NULL,

            #' @param param (OBSOLETE!) The parameters for the prior distribution(s) for
            #' the hyperparmater(s)
            param = NULL
    ) {
        as.list(environment())
    }

`inla.set.f.default` <-
    function(...) {
        list(diagonal = 1.0E-4)
    }

#' @title control.expert
#' @inherit control.update params description seealso 
#' @family control
#' @export
`control.expert` <-
    function(
            #' @param cpo.manual A boolean variable to decide if the inla-program is to be
            #' runned in a manual-cpo-mode. (EXPERT OPTION: DO NOT USE)
            cpo.manual = FALSE,

            #' @param cpo.idx  The index/indices of the data point(s) to remove. (EXPERT
            #' OPTION: DO NOT USE)
            cpo.idx = -1,

            #' @param disable.gaussian.check Disable the check for fast computations with a
            #' Gaussian likelihood and identity link (default `FALSE`)
            disable.gaussian.check = FALSE,

            #' @param jp An object of class `inla.jp` defining a joint prior
            jp = NULL, 

            #' @param dot.product.gain Output the gain in
            #' optimizing dot-products? (Default `FALSE`)
            dot.product.gain = FALSE,

            #' @param globalconstr Add a global constraint (see `?f` and argument
            #' `extraconstr`). Note that a global constraint does NOT
            #' correct the normalisation constant. 
            globalconstr = list(A = NULL, e = NULL)
        ) {
            as.list(environment())
        }


#' @title control.gcpo
#' @inherit control.update params description seealso 
#' @family control
#' @export
#' @details (For experts only!) Set control variables for the gcpo in [control.compute].
#' The intended use is to use `inla.group.cv`.
#' Refer to `?inla.group.cv` and the vignette for details.
control.gcpo <-
    function(
        #' @param enable TODO
        enable = FALSE,
        #' @param num.level.sets TODO
        num.level.sets = -1,
        #' @param size.max TODO
        size.max = 32,
        #' @param strategy TODO
        strategy = c("posterior", "prior"),
        #' @param groups TODO
        groups = NULL,
        #' @param selection TODO
        selection = NULL,
        #' @param friends TODO
        friends = NULL, 
        #' @param verbose TODO
        verbose = FALSE,
        #' @param epsilon TODO
        epsilon = 0.005,
        #' @param prior.diagonal TODO
        prior.diagonal = 1e-4, 
        #' @param correct.hyperpar TODO
        correct.hyperpar = TRUE,
        #' @param keep TODO
        keep = NULL,
        #' @param remove TODO
        remove = NULL,
        #' @param remove.fixed TODO
        remove.fixed = TRUE
    ) {
        as.list(environment())
    }

#' @title control.compute
#' @inherit control.update params description seealso 
#' @family control
#' @export
`control.compute` <-
    function(
            #' @param openmp.strategy The computational strategy to use: 'small', 'medium',
            #' 'large', 'huge', 'default' and 'pardiso'.
            openmp.strategy = "default",

            #' @param hyperpar A boolean variable if the marginal for the hyperparameters
            #' should be computed. Default TRUE.
            hyperpar = TRUE,

            #' @param return.marginals A boolean variable if the marginals for the latent
            #' field should be returned (although it is computed). Default TRUE
            return.marginals = TRUE,

            #' @param return.marginals.predictor A boolean variable if the marginals for
            #' the linear predictor should be returned (although it is computed).
            #' Default FALSE
            return.marginals.predictor = FALSE,

            #' @param dic A boolean variable if the DIC-value should be computed. Default
            #' FALSE.
            dic = FALSE,

            #' @param mlik A boolean variable if the marginal likelihood should be computed.
            #' Default `TRUE`.
            mlik = TRUE,

            #' @param cpo A boolean variable if the cross-validated predictive measures
            #' (cpo, pit) should be computed (default `FALSE`)
            cpo = FALSE,

            #' @param po A boolean variable if the predictive ordinate should be computed
            #' (default `FALSE`)
            po = FALSE,

            #' @param waic A boolean variable if the Watanabe-Akaike information criteria
            #' should be computed (default `FALSE`)
            waic = FALSE,

            #' @param residuals Provide estimates of `residuals` (whatever we mean by that).
            #' (default `FALSE`)
            #' Currently only residuals base on expected (saturated) deviance are available.
            #' The sign of the residuals are only `very likely` correct.
            #' These residuals are not properly justified from a Bayesian point of
            #' view, hence must be used with caution. It is provided in the hope
            #' they would be useful. This feature is EXPERIMENTAL for the moment,
            #' so changes can happen at any time.
            residuals = FALSE,

            #' @param q A boolean variable if binary images of the precision matrix, the
            #' reordered precision matrix and the Cholesky triangle should be generated.
            #' (Default FALSE.)
            q = FALSE,

            #' @param config A boolean variable if the internal GMRF approximations be
            #' stored. (Default FALSE. EXPERIMENTAL)
            config = FALSE,

            #' @param likelihood.info A boolean variable to store likelihood-information or not.
            #' This option requires `config=TRUE` (Default FALSE. EXPERIMENTAL)
            likelihood.info = FALSE,

            #' @param smtp The sparse-matrix solver, one of 'default', 'taucs', 'band' or
            #' 'pardiso' (default `inla.getOption("smtp")`). `smtp='pardiso'` implies
            #' `openmp.strategy='pardiso'`.
            smtp = NULL,

            #' @param graph A boolean variable if the graph itself should be returned.
            #' (Default FALSE.)
            graph = FALSE,

            #' @param internal.opt A boolean variable, if to do internal online
            #' optimisations or not. (Default TRUE.)
            internal.opt = TRUE,

            #' @param save.memory A boolean variable, make choices which
            #' saves memory over accuracy. (Default 'inla.getOption("save.memory")')
            save.memory = NULL, 

            #' @param control.gcpo (For experts only!) Set control variables for the gcpo.
            #' The intended use is to use `inla.group.cv`.
            #' Refer to [control.gcpo], `?inla.group.cv` and the vignette for details.
            control.gcpo = INLA::control.gcpo()
            ## INLA:: needed to avoid name ambiguity with the parameter itself
            ## and to allow calling the function without INLA in the namespace.
            ## During development, use
            ##   control.compute = list(control.gcpo = control.gcpo())
            ## to test updated settings.
        ) {
            as.list(environment())
        }

#' @title control.lp.scale
#' @inherit control.update params description seealso 
#' @family control
#' @export
`control.lp.scale` <-
    function(
            #' @param hyper Definition of the hyperparameter(s)
            hyper = NULL
    ) {
        as.list(environment())
    }

#' @title control.pardiso
#' @inherit control.update params description seealso 
#' @details Extra options controlling the PARDISO library
#' @family control
#' @export
`control.pardiso` <-
    function(
            #' @param verbose Show detailed output (default FALSE)
            verbose = FALSE,
        
            #' @param debug Show internal debug output (default FALSE)
            debug = FALSE,
            
            #' @param parallel.reordering Do reordering in parallel (default TRUE)
            parallel.reordering = TRUE,
            
            #' @param nrhs Number of right-hand sides to solve for in parallel (`-1`
            #' will determine this adapative)
            nrhs = -1
    ) {
        as.list(environment())
    }

#' @title control.bgev
#' @inherit control.update params description seealso 
#' @family control
#' @export
`control.bgev` <-
    function(
        #' @details The `control.bgev`-list is set within the corresponding `control.family`-list as control parameters to the `family="bgev"`
            #' @param q.location The quantile level for the location parameter
            q.location = 0.5,

            #' @param q.spread The quantile level for the spread parameter (must be < 0.5)
            q.spread = 0.25,

            #' @param q.mix The lower and upper quantile level for the mixing function
            q.mix = c(0.10, 0.20),

            #' @param beta.ab The parameters a and b in the Beta mixing function
            beta.ab = 5L
        ) {
            as.list(environment())
        }

#' @title control.family
#' @inherit control.update params description seealso 
#' @family control
#' @export
`control.family` <-
    function(
            #' @param dummy A dummy argument that can be used as a workaround
            dummy = 0,

            #' @param hyper Definition of the hyperparameters
            hyper = NULL,

            #' @param initial (OBSOLETE!) Initial value for the hyperparameter(s) of the
            #' likelihood in the internal scale.
            initial = NULL,

            #' @param prior (OBSOLETE!) The name of the prior distribution(s) for othe
            #' hyperparameter(s).
            prior = NULL,

            #' @param param (OBSOLETE!) The parameters for the prior distribution
            param = NULL,

            #' @param fixed (OBSOLETE!) Boolean variable(s) to say if the hyperparameter(s)
            #' is fixed or random.
            fixed = NULL,

            #' @param link (OBSOLETE! Use `control.link=list(model=)` instead.) The
            #' link function to use.
            link = "default",

            #' @param sn.shape.max Maximum value for the shape-parameter for Skew Normal
            #' observations (default 5.0)
            sn.shape.max = 5.0,

            #' @param gev.scale.xi (Expert option, do not use unless you know what you are
            #' doing.) The internal scaling of the shape-parameter for the GEV distribution.
            #' (default 0.1)
            gev.scale.xi = 0.1,

            #' @param control.bgev  See `?control.bgev`
            control.bgev = NULL,

            #' @param cenpoisson.I The censoring interval for the censored Poisson
            cenpoisson.I = c(-1L, -1L),

            #' @param beta.censor.value The censor value for the Beta-likelihood `(0
            #' <= beta.censor.value < 1/2)`
            beta.censor.value = 0.0,

            #' @param variant This variable is used to give options for various variants of
            #' the likelihood,  like chosing different parameterisations for example. See the
            #' relevant likelihood documentations for options (does only apply to some
            #' likelihoods).
            variant = 0L,

            #' @param control.mix See `?control.mix`
            control.mix = NULL,

            #' @param control.pom See `?control.pom`
            control.pom = NULL,

            #' @param control.link See `?control.link`
            control.link = NULL,

            #' @param link.simple See `inla.doc("0inflated")`
            link.simple = "default"
        ) {
            as.list(environment())
        }

#' @title control.fixed
#' @inherit control.update params description seealso 
#' @family control
#' @export
`control.fixed` <-
    function(
            #' @param cdf  A list of values to compute the CDF for, for all fixed effects
            cdf = NULL,

            #' @param quantiles  A list of quantiles to compute for all fixed effects
            quantiles = NULL,

            #' @param expand.factor.strategy The strategy used to expand factors into fixed
            #' effects based on their levels. The default strategy is us use the
            #' `model.matrix`-function for which NA's are not allowed
            #' (`expand.factor.strategy="model.matrix"`) and levels are possible removed.
            #' The alternative option (`expand.factor.strategy="inla"`) use an
            #' `inla`-specific expansion which expand a factor into one fixed effects for
            #' each level, do allow for NA's and all levels are present in the model. In this
            #' case, factors MUST BE factors in the data.frame/list and NOT added as
            #' `.+factor(x1)+.` in the formula only.
            expand.factor.strategy = "model.matrix",

            #' @param mean Prior mean for all fixed effects except the intercept.
            #'  Alternatively, a named list with specific means where name=default applies to
            #'  unmatched names. For example `control.fixed=list(mean=list(a=1, b=2,
            #'  default=0))` assign 'mean=1' to fixed effect 'a' , 'mean=2' to effect 'b' and
            #' 'mean=0' to all others. (default 0.0)
            mean = 0.0,

            #' @param mean.intercept Prior mean for the intercept (default 0.0)
            mean.intercept = 0.0,

            #' @param prec  Default precision for all fixed effects except the intercept.
            #' Alternatively, a named list with specific means where name=default applies to
            #' unmatched names.  For example `control.fixed=list(prec=list(a=1, b=2,
            #' default=0.01))` assign 'prec=1' to fixed effect 'a' , 'prec=2' to effect 'b' and
            #' 'prec=0.01' to all others. (default 0.001)
            prec = 0.001,

            #' @param prec.intercept  Default precision the intercept (default 0.0)
            prec.intercept = 0.0,

            #' @param compute Compute marginals for the fixed effects ? (default TRUE)
            compute = TRUE,

            #' @param correlation.matrix Compute the posterior correlation matrix for all
            #'  fixed effects? (default FALSE) OOPS: This option will set up appropriate linear
            #'  combinations and the results are shown as the posterior correlation matrix of the
            #'  linear combinations. This option will imply
            #' `control.inla=list(lincomb.derived.correlation.matrix=TRUE)`.
            correlation.matrix = FALSE,

            #' @param remove.names A vector of names of expanded fixed effects to remove
            #' from the model-matrix. This is an expert option, and should only be used if you
            #' know what you are doing.
            remove.names = NULL
        ) {
            as.list(environment())
        }

#' @title control.fixed
#' @inherit control.update params description seealso 
#' @family control
#' @export
#' @details control.vb List of arguments for various VB corrections.
#' Used for [control.inla] `control.vb` specifications.
control.vb <- function(
        #' @param enable Logical/Character Use this feature? If `"auto"` this will be selected automatically. 
        enable = "auto",
        #' @param strategy Character What to correct, either "mean" or "variance".
        strategy = c("mean", "variance"),
        #' @param verbose Logical Be verbose or not. 
        verbose = TRUE,
        #' @param iter.max Integer Maximum number of iterations. 
        iter.max = 25,
        #' @param emergency Numeric If the standardized correction for the mean is larger
        #'           than this value, then call the vb.correction off and issue a warning
        emergency = 25,
        #' @param f.enable.limit Vector of length 2. The size limit to correct for a `f()`.
        #'         First element is for `strategy="mean"`.
        #'         Second element is for `strategy="mean"`.
        f.enable.limit = c(30, 25),
        #' @param hessian.update How many times the Hessian is updated for each
        #'         correction (`strategy="variance"` only).
        hessian.update = 2,
        #' @param hessian.strategy Select strategy for computing the Hessian
        #'         matrix for `strategy="variance"`, one of `"full"`, 
        #'         `"diagonal"`, `"partial"` and `"default"`.
        hessian.strategy = c("default", "full", "partial", "diagonal")
) {
    as.list(environment())
}


#' @title control.inla
#' @inherit control.update params description seealso 
#' @family control
#' @export
`control.inla` <-
    function(
        #' @param strategy  Character The strategy to use for the approximations; one
        #' of 'auto' (default), 'gaussian', 'simplified.laplace', 'laplace' or 'adaptive'.
        strategy = "auto",
        
        #' @param int.strategy  Character The integration strategy to use; one of
        #' 'auto' (default),  'ccd', 'grid', 'eb' (empirical bayes),  'user' or 'user.std'.
        #' For the experimental mode,  then 'grid' equal 'ccd' for more than two
        #' hyperparameters.
        int.strategy = "auto",
        
        #' @param int.design  Matrix Matrix of user-defined integration points and
        #' weights. Each row consists theta values and the integration weight.
        #' (EXPERIMENTAL!).
        int.design = NULL,
        
        #' @param interpolator  Character The interpolator used to compute the
        #' marginals for the hyperparameters. One of 'auto', 'nearest', 'quadratic',
        #' 'weighted.distance', 'ccd', 'ccdintegrate', 'gridsum', 'gaussian'. Default is
        #' 'auto'.
        interpolator = "auto",
        
        #' @param fast Logical If TRUE, then replace conditional modes in the Laplace
        #' approximation with conditional expectation (default TRUE).
        fast = TRUE,
        
        #' @param linear.correction Logical Default TRUE for the 'strategy = laplace'
        #' option.
        linear.correction = NULL,
        
        #' @param h Numerical The step-length for the gradient calculations for the
        #' hyperparameters. Default 0.005.
        h = 0.005,
        
        #' @param dz Numerical The step-length in the standarised scale for the
        #' integration of the hyperparameters. Default 0.75.
        dz = 0.75,
        
        #' @param diff.logdens Numerical The difference of the log.density for the
        #' hyperpameters to stop numerical integration using int.strategy='grid'. Default 6.
        diff.logdens = 6,
        
        #' @param print.joint.hyper Logical If TRUE, the store also the joint
        #' distribution of the hyperparameters (without any costs). Default TRUE.
        print.joint.hyper = TRUE,
        
        #' @param force.diagonal Logical If TRUE, then force the Hessian to be
        #' diagonal. (Default `FALSE`)
        force.diagonal = FALSE,
        
        #' @param skip.configurations Logical Skip configurations if the values at the
        #' main axis are to small. (Default `TRUE`)
        skip.configurations = TRUE,
        
        #' @param mode.known Logical If TRUE then no optimisation is done. (Default
        #' FALSE.)
        mode.known = FALSE,
        
        #' @param adjust.weights Logical If TRUE then just more accurate integration
        #' weights. (Default TRUE.)
        adjust.weights = TRUE,
        
        #' @param tolerance Numerical The tolerance for the optimisation of the
        #' hyperparameters. If set, this is the default value for for '2.5*tolerance.f',
        #' 'tolerance.g' and  '5*tolerance.x'; see below.
        tolerance = 0.005,
        
        #' @param tolerance.f Numerical The tolerance for the absolute change in the
        #' log posterior in the optimisation of the hyperparameters.
        tolerance.f = NULL,
        
        #' @param tolerance.g Numerical The tolerance for the absolute change in the
        #' gradient of the log posterior in the optimisation of the hyperparameters.
        tolerance.g = NULL,
        
        #' @param tolerance.x Numerical The tolerance for the change in the
        #' hyperparameters (root-mean-square) in the optimisation of the hyperparameters.
        tolerance.x = NULL,
        
        #' @param tolerance.step Numerical The tolerance for the change in
        #' root-mean_squre in the inner Newton-like optimisation of the latent field.
        tolerance.step = 0.001,
        
        #' @param restart Numerical To improve the optimisation, the optimiser is
        #' restarted at the found optimum 'restart' number of times.
        restart = 0L,
        
        #' @param optimiser Character The optimiser to use; one of 'gsl' or 'default'.
        optimiser = "default",
        
        #' @param verbose Logical Run in verbose mode? (Default FALSE)
        verbose = NULL,
        
        #' @param reordering Character Type of reordering to use. (EXPERT OPTION; one
        #' of "AUTO", "DEFAULT", "IDENTITY", "REVERSEIDENTITY",  "BAND", "METIS", "GENMMD",
        #' "AMD", "MD", "MMD", "AMDBAR", "AMDC", "AMDBARC",  or the output from
        #' `inla.qreordering`. Default is 'auto'.)
        reordering = "auto",
        
        #' @param cpo.diff Numerical Threshold to define when the cpo-calculations are
        #' inaccurate. (EXPERT OPTION.)
        cpo.diff = NULL,
        
        #' @param npoints Numerical Number of points to use in the 'stratey=laplace'
        #' approximation (default 9)
        npoints = 9,
        
        #' @param cutoff Numerical The cutoff used in the 'stratey=laplace'
        #' approximation. (Smaller value is more accurate and more slow.) (default 1e-4)
        cutoff = 1e-4,
        
        #' @param adapt.hessian.mode Logical Should optimisation be continued if the
        #' Hessian estimate is void? (Default TRUE)
        adapt.hessian.mode = NULL,
        
        #' @param adapt.hessian.max.trials Numerical Number of steps in the adaptive
        #' Hessian optimisation
        adapt.hessian.max.trials = NULL,
        
        #' @param adapt.hessian.scale Numerical The scaling of the 'h' after each
        #' trial.
        adapt.hessian.scale = NULL,
        
        #' @param adaptive.max Selecting `strategy="adaptive"` will chose the
        #' default strategy for all fixed effects and model components with length less or
        #' equal to `adaptive.max`, for others, the gaussian strategy will be applied.
        adaptive.max = 25L,
        
        #' @param huge Logical If TRUE then try to do some of the internal
        #' parallelisations differently. Hopefully this will be of benefit for 'HUGE' models.
        #' (Default FALSE.) THIS OPTION IS OBSOLETE AND NOT USED!
        huge = FALSE,
        
        #' @param step.len Numerical The step-length used to compute numerical
        #' derivaties of the log-likelihood (0 means `default` which
        #' depends on `stencil`)
        step.len = 0.0, 
        
        #' @param stencil Numerical Number of points in the stencil used to compute the
        #' numerical derivaties of the log-likelihood (5, 7 or 9). (default 5)
        stencil = 5L,
        
        #' @param lincomb.derived.correlation.matrix Logical If TRUE compute also the
        #' correlations for the derived linear combinations, if FALSE do not (Default FALSE)
        lincomb.derived.correlation.matrix = FALSE,
        
        #' @param diagonal Numerical Expert use only! Add a this value on the diagonal
        #' of the joint precision matrix. (default 0.0)
        diagonal = 0.0,
        
        #' @param numint.maxfeval Numerical Maximum number of function evaluations in
        #' the the numerical integration for the hyperparameters. (Default 100000.)
        numint.maxfeval = 100000,
        
        #' @param numint.relerr Numerical Relative error requirement in the the
        #' numerical integration for the hyperparameters. (Default 1e-5)
        numint.relerr = 1e-5,
        
        #' @param numint.abserr Numerical Absolute error requirement in the the
        #' numerical integration for the hyperparameters. (Default 1e-6)
        numint.abserr = 1e-6,
        
        #' @param cmin Numerical The minimum value for the negative Hessian from the
        #' likelihood. Increasing this value will stabalise the optimisation but can
        #' introduce bias.  (Default -Inf)
        cmin = -Inf,
        
        #' @param b.strategy Character If `cmin` is used, either keep the linear
        #'     term (with `b.strategy="keep"`) or skip the contribution by
        #'     setting the linear term to zero (`b.strategy="skip"`). The
        #'     default value is `"skip"`
        b.strategy = "skip",
        
        #' @param step.factor Numerical The step factor in the Newton-Raphson algorithm
        #' saying how large step to take (Default 1.0) YES! setting this to a
        #' negative values means = 1,  EXCEPT the first time (for each thread)
        #' where |step.factor| is used. 
        step.factor = -0.1,
        
        #' @param global.node.factor Numerical The factor which defines the degree
        #'     required (how many neighbors), as a fraction of `n-1`,
        #'     that is required to be classified as a global node and numbered
        #'     last (whatever the reordering routine says). Here, `n`,
        #'     is the size of the graph. (Disabled if larger than 1,  default 2)
        global.node.factor = 2.0,
        
        #' @param global.node.degree Numerical The degree required (number of
        #'     neighbors) to be classified as a global node and numbered
        #'     last (whatever the reordering routine says).
        #'     (default `.Machine$integer.max`)
        global.node.degree = .Machine$integer.max,
        
        #' @param stupid.search Logical Enable or disable the stupid-search-algorithm,
        #' if the Hessian calculations reveals that the mode is not found.
        #' (Default `TRUE`.)
        stupid.search = TRUE,
        
        #' @param stupid.search.max.iter Numerical Maximum number of iterations allowed
        #' for the stupid-search-algorithm. (default 1000)
        stupid.search.max.iter = 1000L,
        
        #' @param stupid.search.factor Numerical Factor (>=1) to increase the
        #' step-length with after each new iteration. (default 1.05)
        stupid.search.factor = 1.05,
        
        #' @param control.vb list of arguments for various VB corrections.
        #' See [control.vb()] for details.
        control.vb = INLA::control.vb(),
        ## INLA:: needed to avoid name ambiguity with the parameter itself
        ## and to allow calling the function without INLA in the namespace.
        ## During development, use
        ##   control.inla = list(control.vb = control.vb())
        ## to test updated settings.
        
        #' @param num.gradient Character Set the numerical scheme to compute the
        #' gradient,  one of `"forward"` or `"central"` (default).
        num.gradient = "central",
        
        #' @param num.hessian Character Set the numerical scheme to compute the
        #' Hessian,  one of `"forward"` or `"central"` (default).
        num.hessian = "central",
        
        #' @param optimise.strategy Character THIS OPTION IS EXPERIMENTAL. Chose the
        #' optimiser strategy,  one of `"plain"` or `"smart"` (default)
        optimise.strategy = "smart",
        
        #' @param use.directions THIS OPTION IS EXPERIMENTAL. Unless `FALSE` or
        #' `NULL`,  use directions for computing gradient and Hessian, initialised with
        #' `use.directions` if a matrix.
        use.directions = TRUE,
        
        #' @param constr.marginal.diagonal Add stability to `AQ^-1A^T` by adding a
        #' small diagonal term. (default `epsilon^0.5`)
        constr.marginal.diagonal = sqrt(.Machine$double.eps),
        
        #' @param improved.simplified.laplace If `TRUE` use an experimental
        #' improved variant, otherwise, use the standard one.
        improved.simplified.laplace = FALSE,
        
        #' @param parallel.linesearch Use serial (default) or parallel line-search
        #' (highly experimental for the moment)
        parallel.linesearch = FALSE,
        
        #' @param compute.initial.values Compute initial values for the latent field or not.
        #' (experimental-mode only)
        compute.initial.values = TRUE,
        
        #' @param hessian.correct.skewness.only If TRUE then correct only
        #' skewness in the Hessian, for the hyperparameters. If FALSE (default),
        #' correct also variance (experimental-mode only)
        hessian.correct.skewness.only = FALSE
    ) {
        as.list(environment())
    }


#' @title control.predictor
#' @inherit control.update params description seealso 
#' @family control
#' @export
`control.predictor` <-
    function(
            #' @param hyper Definition of the hyperparameters.
            hyper = NULL,

            #' @param fixed (OBSOLETE!) If the precision for the artificial noise is fixed
            #' or not (default TRUE)
            fixed = NULL,

            #' @param prior (OBSOLETE!) The prior for the artificial noise
            prior = NULL,

            #' @param param (OBSOLETE!) Prior parameters for the artificial noise
            param = NULL,

            #' @param initial (OBSOLETE!) The value of the log precision of the artificial
            #' noise
            initial = NULL,

            #' @param compute A boolean variable; should the marginals for the linear
            #' predictor be computed? (Default FALSE.)
            compute = FALSE,

            #' @param cdf A list of values to compute the CDF for the linear predictor
            cdf = NULL,

            #' @param quantiles A list of quantiles to compute for the linear predictor
            quantiles = NULL,

            #' @param cross Cross-sum-to-zero constraints with the linear predictor. All
            #' linear predictors with the same level of 'cross' are constrained to have sum
            #' zero. Use 'NA' for no contribution. 'Cross' has the same length as the linear
            #' predictor (including the 'A' matrix extention). (THIS IS AN EXPERIMENTAL OPTION,
            #' CHANGES MAY APPEAR.)
            cross = NULL,

            #' @param A The observation matrix (matrix or Matrix::sparseMatrix).
            A = NULL,

            #' @param precision The precision for eta* - A*eta, (default `exp(15)`)
            precision = exp(15),

            #' @param link Define the family-connection for unobserved observations
            #' (`NA`). `link` is integer values which defines the family connection;
            #' `family[link[idx]]` unless `is.na(link[idx])` for which the
            #' identity-link is used. The `link`-argument only influence the
            #' `fitted.values` in the `result`-object. If `is.null(link)`
            #' (default) then the identity-link is used for all missing observations. If the
            #' length of `link` is 1, then this value is replicated with the length of the
            #' responce vector. If an element of the responce vector is `!NA` then the
            #' corresponding entry in `link` is not used (but must still be a legal value).
            #' Setting this variable implies `compute=TRUE`.
            link = NULL
        ) {
            as.list(environment())
        }

#' @title control.mode
#' @inherit control.update params description seealso 
#' @family control
#' @export
#' @details For internal use and for algorithms built on to of INLA.
`control.mode` <-
    function(
            #' @param result Previous result from inla(). Use the theta- and x-mode from
            #' this run.
            result = NULL,

            #' @param theta The theta-mode/initial values for theta. This option has
            #' preference over result$mode$theta.
            theta = NULL,

            #' @param x The x-mode/initial values for x. This option has preference over
            #' result$mode$x.
            x = NULL,

            #' @param restart A boolean variable; should we restart the optimisation from
            #' this configuration or fix the mode at this configuration? (Default FALSE.)
            restart = FALSE,

            #' @param fixed A boolean variable. If TRUE then treat all thetas as known and
            #' fixed, and if FALSE then treat all thetas as unknown and random (default).
            fixed = FALSE
        ) {
            as.list(environment())
        }

#' @title control.hazard
#' @inherit control.update params description seealso 
#' @family control
#' @export
`control.hazard` <-
    function(
            #' @param model The model for the baseline hazard model. One of 'rw1', 'rw2' or
            #' 'iid'. (Default 'rw1'.)
            model = "rw1",

            #' @param hyper The definition of the hyperparameters.
            hyper = NULL,

            #' @param fixed (OBSOLETE!) A boolean variable; is the precision for 'model'
            #' fixed? (Default FALSE.)
            fixed = FALSE,

            #' @param initial (OBSOLETE!) The initial value for the precision.
            initial = NULL,

            #' @param prior (OBSOLETE!) The prior distribution for the precision for
            #' 'model'
            prior = NULL,

            #' @param param (OBSOLETE!) The parameters in the prior distribution
            param = NULL,

            #' @param constr A boolean variable; shall the  'model' be constrained to sum
            #' to zero?
            constr = TRUE,

            #' @param diagonal An extra constant added to the diagonal of the precision
            #' matrix
            diagonal = NULL,

            #' @param n.intervals Number of intervals in the baseline hazard. (Default 15)
            n.intervals = 15,

            #' @param cutpoints The cutpoints to use. If not specified the they are compute
            #' from 'n.intervals' and the maximum length of the interval. (Default NULL)
            cutpoints = NULL,

            #' @param strata.name The name of the stratefication variable for the baseline
            #' hazard in the data.frame
            strata.name = NULL,

            #' @param scale.model Scale the baseline hazard model (RW1, RW2) so the
            #' generalized variance is 1. (Default
            #' `inla.getOption("scale.model.default")`.)
            scale.model = NULL
        ) {
            as.list(environment())
        }

## check control-arguments

`inla.check.control` <- function(contr, data = NULL) {
    ## This function will signal an error if the arguments in CONTR
    ## does not match the ones in the corresponding
    ## `inla.set.XX.default()' routine.  EXAMPLE: contr is
    ## `control.inla' and default arguments is found in
    ## `inla.set.control.inla.default()'

    ## Will expand unexpanded names from the names in 'data' first
    contr <- local({
        name <- paste("inla.tmp.env", as.character(runif(1)), sep = "")
        attach(data, name = name, warn.conflicts = FALSE)
        ccontr <- contr
        detach(name, character.only = TRUE)
        ccontr
    })

    stopifnot(!missing(contr))
    if (length(contr) == 0) {
        return(contr)
    }

    nm <- paste(sys.call()[2])
    f <- paste("inla.set.", nm, ".default()", sep = "")
    elms <- names(inla.eval(f))

    if (is.null(names(contr))) {
        stop(inla.paste(c(
            "Named elements in in control-argument `", nm, "', is required: ", contr,
            "\n\n  Valid ones are:\n\t",
            inla.paste(sort(elms), sep = "\n\t")
        ), sep = ""))
    }

    for (elm in names(contr)) {
        if (!is.element(elm, elms)) {
            stop(inla.paste(c(
                "Name `", elm, "' in control-argument `", nm, "', is void.\n\n  Valid ones are:\n\t",
                inla.paste(sort(elms), sep = "\n\t")
            ), sep = ""))
        }
    }

    return(contr)
}



## test-implementation
## `control.lincomb` = function(precision, verbose)
## {
## aa = match.call()[-1]
## ret = list()
## for(a in names(aa)) {
## if (!missing(a)) {
## xx = get(a)
## names(xx) = a
## ret = c(ret, xx)
## }
## }
## return (ret)
## }

#inla.make.completion.function <- function(...) {
#    my.eval <- function(command, envir = parent.frame(),
#                        enclos = if (is.list(envir) || is.pairlist(envir)) parent.frame() else baseenv()) {
#        return(eval(parse(text = command), envir, enclos))
#    }
#    
#    xx <- sort(list(...)[[1L]])
#    my.eval(paste("function(", paste(xx, sep = "", collapse = ", "), ") {
#    aa = match.call()[-1L]
#    ret = list()
#    for(a in names(aa)) {
#        if (!missing(a)) {
#            xx = get(a)
#            names(xx) = a
#            ret = c(ret, xx)
#        }
#    }
#    return (ret)
#}"))
#}


#' @rdname control.update
#' @export
`inla.set.control.update.default` <- function(...) { control.update(...) }

#' @rdname control.lincomb
#' @export
`inla.set.control.lincomb.default` <- function(...) { control.lincomb(...) }

#' @rdname control.group
#' @export
`inla.set.control.group.default` <- function(...) { control.group(...) }

#' @rdname control.scopy
#' @export
`inla.set.control.scopy.default` <- function(...) { control.scopy(...) }

#' @rdname control.mix
#' @export
`inla.set.control.mix.default` <- function(...) { control.mix(...) }

#' @rdname control.pom
#' @export
`inla.set.control.pom.default` <- function(...) { control.pom(...) }

#' @rdname control.link
#' @export
`inla.set.control.link.default` <- function(...) { control.link(...) }

#' @rdname control.expert
#' @export
`inla.set.control.expert.default` <- function(...) { control.expert(...) }

#' @rdname control.gcpo
#' @export
`inla.set.control.gcpo.default` <- function(...) { control.gcpo(...) }

#' @rdname control.compute
#' @export
`inla.set.control.compute.default` <- function(...) { control.compute(...) }

#' @rdname control.lp.scale
#' @export
`inla.set.control.lp.scale.default` <- function(...) { control.lp.scale(...) }

#' @rdname control.pardiso
#' @export
`inla.set.control.pardiso.default` <- function(...) { control.pardiso(...) }

#' @rdname control.bgev
#' @export
`inla.set.control.bgev.default` <- function(...) { control.bgev(...) }

#' @rdname control.family
#' @export
`inla.set.control.family.default` <- function(...) { control.family(...) }

#' @rdname control.fixed
#' @export
`inla.set.control.fixed.default` <- function(...) { control.fixed(...) }

#' @rdname control.vb
#' @export
`inla.set.control.vb.default` <- function(...) { control.vb(...) }

#' @rdname control.inla
#' @export
`inla.set.control.inla.default` <- function(...) { control.inla(...) }

#' @rdname control.predictor
#' @export
`inla.set.control.predictor.default` <- function(...) { control.predictor(...) }

#' @rdname control.mode
#' @export
`inla.set.control.mode.default` <- function(...) { control.mode(...) }

#' @rdname control.hazard
#' @export
`inla.set.control.hazard.default` <- function(...) { control.hazard(...) }
