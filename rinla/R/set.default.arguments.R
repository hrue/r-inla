## Export: inla.set.control.bgev.default
## Export: inla.set.control.compute.default
## Export: inla.set.control.expert.default
## Export: inla.set.control.family.default
## Export: inla.set.control.fixed.default
## Export: inla.set.control.group.default
## Export: inla.set.control.hazard.default
## Export: inla.set.control.inla.default
## Export: inla.set.control.lincomb.default
## Export: inla.set.control.link.default
## Export: inla.set.control.lp.scale.default
## Export: inla.set.control.mix.default
## Export: inla.set.control.mode.default
## Export: inla.set.control.pardiso.default
## Export: inla.set.control.pom.default
## Export: inla.set.control.predictor.default
## Export: inla.set.control.scopy.default
## Export: inla.set.control.update.default

## Export: control.bgev
## Export: control.compute
## Export: control.expert
## Export: control.family
## Export: control.fixed
## Export: control.gcpo
## Export: control.group
## Export: control.hazard
## Export: control.inla
## Export: control.lincomb
## Export: control.link
## Export: control.lp.scale
## Export: control.mix
## Export: control.mode
## Export: control.pardiso
## Export: control.pom
## Export: control.predictor
## Export: control.scopy
## Export: control.update
## Export: control.vb


### Defines default arguments

`inla.set.control.update.default` <-
    function(...) {
        ## :EXTRA:
        ## :NAME: control.update
        list(
            ## :ARGUMENT: result Update the joint posterior for the hyperparameters from result
            result = NULL
        )

        ## :SEEALSO: inla
    }

`inla.set.control.lincomb.default` <-
    function(...) {
        ## :EXTRA:
        ## :NAME: control.lincomb
        list(
            ## :ARGUMENT+: verbose Use verbose mode for linear combinations if verbose model is
            ## :ARGUMENT+: set globally. (Default TRUE). This option is only available for
            ## :ARGUMENT+: the default \code{inla.mode} (\code{inla.mode="compact"}).
            verbose = FALSE
        )

        ## :SEEALSO: inla
    }

`inla.set.control.group.default` <-
    function(...) {
        ## :EXTRA:
        ## :NAME: control.group
        list(
            ## :ARGUMENT+: model Group model (one of 'exchangable', 'exchangablepos', 'ar1',
            ## :ARGUMENT+: 'ar', 'rw1', 'rw2', 'besag', or 'iid')
            model = "exchangeable",

            ## :ARGUMENT+: order Defines the \code{order} of the model: for model \code{ar} this
            ## :ARGUMENT+: defines the order p, in AR(p). Not used for other models at the time
            ## :ARGUMENT+: being.
            order = NULL,

            ## :ARGUMENT+: cyclic Make the group model cyclic? (Only applies to models 'ar1',
            ## :ARGUMENT+: 'rw1' and 'rw2')
            cyclic = FALSE,

            ## :ARGUMENT: graph The graph spesification (Only applies to model 'besag')
            graph = NULL,

            ## :ARGUMENT+: scale.model Scale the intrinsic model (RW1, RW2, BESAG) so the
            ## :ARGUMENT+: generalized variance is 1. (Default \code{TRUE})
            scale.model = TRUE,

            ## :ARGUMENT+: adjust.for.con.comp Adjust for connected components when
            ## :ARGUMENT+: \code{scale.model=TRUE}? (default \code{TRUE})
            adjust.for.con.comp = TRUE,

            ## :ARGUMENT: hyper Definition of the hyperparameter(s)
            hyper = NULL,

            ## :ARGUMENT+: initial (OBSOLETE!) The initial value for the group correlation or
            ## :ARGUMENT+: precision in the internal scale.
            initial = NULL,

            ## :ARGUMENT+: fixed (OBSOLETE!) A boolean variable if the group correction or
            ## :ARGUMENT+: precision is assumed to be fixed or random.
            fixed = NULL,

            ## :ARGUMENT+: prior (OBSOLETE!) The name of the prior distribution for the group
            ## :ARGUMENT+: correlation or precision in the internal scale
            prior = NULL,

            ## :ARGUMENT: param (OBSOLETE!) Prior parameters
            param = NULL
        )

        ## :SEEALSO: inla
    }

`inla.set.control.scopy.default` <-
    function(...) {
        ## :EXTRA:
        ## :NAME: control.scopy
        list(
            
            ## :ARGUMENT: covariate The covariate for the scopy function
            covariate = NULL,

            ## :ARGUMENT: n Number of betas
            n = 5,

            ## :ARGUMENT: model scopy model (one of 'rw1' and 'rw2')
            model = "rw2",

            ## :ARGUMENT: mean The prior mean for mean(betas)
            mean = 1.0,

            ## :ARGUMENT: prec.mean The prior precision for mean(betas)
            prec.mean =  1.0, 

            ## :ARGUMENT: prec.betas The prior precision prec(betas-mean(betas))
            prec.betas =  10.0
        )

        ## :SEEALSO: inla
    }

`inla.set.control.mix.default` <-
    function(...) {
        ## :EXTRA: The \code{control.mix} -list is set within the corresponding \code{control.family}-list a the mixture of the likelihood is likelihood spesific. (This option is EXPERIMENTAL.)
        ## :NAME: control.mix
        list(
            ## :ARGUMENT+: model The model for the random effect. Currently, only
            ## :ARGUMENT+: \code{model='gaussian'} is implemented
            model = NULL,

            ## :ARGUMENT+: hyper Definition of the hyperparameter(s) for the random effect model
            ## :ARGUMENT+: chosen
            hyper = NULL,

            ## :ARGUMENT: initial (OBSOLETE!) The initial value(s) for the hyperparameter(s)
            initial = NULL,

            ## :ARGUMENT+: fixed (OBSOLETE!) A boolean variable if hyperparmater(s) is/are fixed
            ## :ARGUMENT+: or random
            fixed = NULL,

            ## :ARGUMENT+: prior (OBSOLETE!) The name of the prior distribution(s) for the
            ## :ARGUMENT+: hyperparmater(s)
            prior = NULL,

            ## :ARGUMENT+: param (OBSOLETE!) The parameters for the prior distribution(s) for the
            ## :ARGUMENT+: hyperparmater(s)
            param = NULL,

            ## :ARGUMENT+: npoints Number of points used to do the numerical integration (default
            ## :ARGUMENT+: 101)
            npoints = 101,

            ## :ARGUMENT+: integrator The integration scheme to use (\code{default},
            ## :ARGUMENT+: \code{quadrature}, \code{simpson})
            integrator = "default"
        )

        ## :SEEALSO: inla
    }

`inla.set.control.pom.default` <-
    function(...) {
        ## :EXTRA: The \code{control.pom} -list is set within the corresponding \code{control.family}
        ## :NAME: control.pom
        list(
            ## :ARGUMENT: cdf character The cdf to use, "logit" (default) or "probit"
            cdf  = "logit",

            ## :ARGUMENT+: fast Logical Use a faster but approximate form for the probit cdf
            ## :ARGUMENT+: (default \code{FALSE})?
            fast = FALSE
        )

        ## :SEEALSO: inla
    }

`inla.set.control.link.default` <-
    function(...) {
        ## :EXTRA: The \code{control.link}-list is set within the corresponding \code{control.family}-list as the link is likelihood-familiy spesific.
        ## :NAME: control.link
        list(
            ## :ARGUMENT: model The name of the link function/model
            model = "default",

            ## :ARGUMENT+: order The \code{order} of the link function, where the interpretation
            ## :ARGUMENT+: of \code{order} is model-dependent.
            order = NULL,

            ## :ARGUMENT+: variant The \code{variant} of the link function, where the
            ## :ARGUMENT+: interpretation of \code{variant} is model-dependent.
            variant = NULL,

            ## :ARGUMENT: hyper Definition of the hyperparameter(s) for the link model chosen
            hyper = NULL,

            ## :ARGUMENT: quantile The quantile for quantile link function
            quantile = NULL,

            ## :ARGUMENT: a The parameter \code{a} in the LOGa link
            a = 1.0,

            ## :ARGUMENT: initial (OBSOLETE!) The initial value(s) for the hyperparameter(s)
            initial = NULL,

            ## :ARGUMENT+: fixed (OBSOLETE!) A boolean variable if hyperparmater(s) is/are fixed
            ## :ARGUMENT+: or random
            fixed = NULL,

            ## :ARGUMENT+: prior (OBSOLETE!) The name of the prior distribution(s) for the
            ## :ARGUMENT+: hyperparmater(s)
            prior = NULL,

            ## :ARGUMENT+: param (OBSOLETE!) The parameters for the prior distribution(s) for
            ## :ARGUMENT+: the hyperparmater(s)
            param = NULL
        )

        ## :SEEALSO: inla
    }

`inla.set.f.default` <-
    function(...) {
        list(diagonal = 1.0E-4)
    }

`inla.set.control.expert.default` <-
    function(...) {
        ## :EXTRA:
        ## :NAME: control.expert
        list(
            ## :ARGUMENT+: cpo.manual A boolean variable to decide if the inla-program is to be
            ## :ARGUMENT+: runned in a manual-cpo-mode. (EXPERT OPTION: DO NOT USE)
            cpo.manual = FALSE,

            ## :ARGUMENT+: cpo.idx  The index/indices of the data point(s) to remove. (EXPERT
            ## :ARGUMENT+: OPTION: DO NOT USE)
            cpo.idx = -1,

            ## :ARGUMENT+: disable.gaussian.check Disable the check for fast computations with a
            ## :ARGUMENT+: Gaussian likelihood and identity link (default \code{FALSE})
            disable.gaussian.check = FALSE,

            ## :ARGUMENT: jp An object of class \code{inla.jp} defining a joint prior
            jp = NULL, 

            ## :ARGUMENT+: dot.product.gain Output the gain in
            ## :ARGUMENT+: optimizing dot-products? (Default \code{FALSE})
            dot.product.gain = FALSE,

            ## :ARGUMENT+: global.constr Add a global constraint (see \code{?f} and argument
            ## :ARGUMENT+: \code{extraconstr}). Note that a global constraint does NOT
            ## :ARGUMENT+: correct the normalisation constant. 
            globalconstr = list(A = NULL, e = NULL)
        )

        ## :SEEALSO: inla
    }

`inla.set.control.compute.default` <-
    function(...) {
        ## :EXTRA:
        ## :NAME: control.compute
        list(
            ## :ARGUMENT+: openmp.strategy The computational strategy to use: 'small', 'medium',
            ## :ARGUMENT+: 'large', 'huge', 'default' and 'pardiso'.
            openmp.strategy = "default",

            ## :ARGUMENT+: hyperpar A boolean variable if the marginal for the hyperparameters
            ## :ARGUMENT+: should be computed. Default TRUE.
            hyperpar = TRUE,

            ## :ARGUMENT+: return.marginals A boolean variable if the marginals for the latent
            ## :ARGUMENT+: field should be returned (although it is computed). Default TRUE
            return.marginals = TRUE,

            ## :ARGUMENT+: return.marginals.predictor A boolean variable if the marginals for
            ## :ARGUMENT+: the linear predictor should be returned (although it is computed).
            ## :ARGUMENT+: Default FALSE
            return.marginals.predictor = FALSE,

            ## :ARGUMENT+: dic A boolean variable if the DIC-value should be computed. Default
            ## :ARGUMENT+: FALSE.
            dic = FALSE,

            ## :ARGUMENT+: mlik A boolean variable if the marginal likelihood should be computed.
            ## :ARGUMENT+: Default \code{TRUE}.
            mlik = TRUE,

            ## :ARGUMENT+: cpo A boolean variable if the cross-validated predictive measures
            ## :ARGUMENT+: (cpo, pit) should be computed (default \code{FALSE})
            cpo = FALSE,

            ## :ARGUMENT+: po A boolean variable if the predictive ordinate should be computed
            ## :ARGUMENT+: (default \code{FALSE})
            po = FALSE,

            ## :ARGUMENT+: waic A boolean variable if the Watanabe-Akaike information criteria
            ## :ARGUMENT+: should be computed (default \code{FALSE})
            waic = FALSE,

            ## :ARGUMENT+: residuals Provide estimates of `residuals` (whatever we mean by that).
            ## :ARGUMENT+: (default \code{FALSE})
            ## :ARGUMENT+: Currently only residuals base on expected (saturated) deviance are available.
            ## :ARGUMENT+: The sign of the residuals are only `very likely` correct.
            ## :ARGUMENT+: These residuals are not properly justified from a Bayesian point of
            ## :ARGUMENT+: view, hence must be used with caution. It is provided in the hope
            ## :ARGUMENT+: they would be useful. This feature is EXPERIMENTAL for the moment,
            ## :ARGUMENT+: so changes can happen at any time.
            residuals = FALSE,

            ## :ARGUMENT+: q A boolean variable if binary images of the precision matrix, the
            ## :ARGUMENT+: reordered precision matrix and the Cholesky triangle should be generated.
            ## :ARGUMENT+: (Default FALSE.)
            q = FALSE,

            ## :ARGUMENT+: config A boolean variable if the internal GMRF approximations be
            ## :ARGUMENT+: stored. (Default FALSE. EXPERIMENTAL)
            config = FALSE,

            ## :ARGUMENT+: likelihood.info A boolean variable to store likelihood-information or not.
            ## :ARGUMENT+: This option requires \code{config=TRUE} (Default FALSE. EXPERIMENTAL)
            likelihood.info = FALSE,

            ## :ARGUMENT+: smtp The sparse-matrix solver, one of 'default', 'taucs', 'band' or
            ## :ARGUMENT+: 'pardiso' (default \code{inla.getoption("smtp")}). \code{smtp='pardiso'} implies
            ## :ARGUMENT+: \code{openmp.strategy='pardiso'}.
            smtp = NULL,

            ## :ARGUMENT+: graph A boolean variable if the graph itself should be returned.
            ## :ARGUMENT+: (Default FALSE.)
            graph = FALSE,

            ## :ARGUMENT+: internal.opt A boolean variable, if to do internal online
            ## :ARGUMENT+: optimisations or not. (Default TRUE.)
            internal.opt = TRUE,

            ## :ARGUMENT+: save.memory A boolean variable, make choices which
            ## :ARGUMENT+: saves memory over accuracy. (Default FALSE.)
            save.memory = FALSE,

            ## :ARGUMENT+: control.gcpo (For experts only!) Set control variables for the gcpo.
            ## :ARGUMENT+: The intended use is to use \code{inla.group.cv}.
            ## :ARGUMENT+: Refer to \code{?inla.group.cv} and the vignette for details.
            control.gcpo = list(enable = FALSE,
                                num.level.sets = -1,
                                size.max = 32,
                                strategy = c("posterior", "prior"),
                                groups = NULL,
                                selection = NULL,
                                friends = NULL, 
                                verbose = FALSE,
                                epsilon = 0.005,
                                prior.diagonal = 1e-4, 
                                correct.hyperpar = TRUE,
                                keep = NULL,
                                remove = NULL,
                                remove.fixed = TRUE)
        )

        ## :SEEALSO: inla
    }

`inla.set.control.lp.scale.default` <-
    function(...) {
        ## :EXTRA:
        ## :NAME: control.lp.scale
        list(
            ## :ARGUMENT: hyper Definition of the hyperparameter(s)
            hyper = NULL
        )

        ## :SEEALSO: inla
    }

`inla.set.control.pardiso.default` <-
    function(...) {
        ## :EXTRA: Extra options controlling the PARDISO library
        ## :NAME: control.pardiso.default
        list(
            ## :ARGUMENT: verbose Show detailed output (default FALSE)
            verbose = FALSE,

            ## :ARGUMENT: debug Show internal debug output (default FALSE)
            debug = FALSE,

            ## :ARGUMENT: parallel.reordering Do reordering in parallel (default TRUE)
            parallel.reordering = TRUE,

            ## :ARGUMENT+: nrhs Number of right-hand sides to solve for in parallel (\code{-1}
            ## :ARGUMENT+: will determine this adapative)
            nrhs = -1
        )

        ## :SEEALSO: inla
    }

`inla.set.control.bgev.default` <-
    function(...) {
        ## :EXTRA: The \code{control.bgev}-list is set within the corresponding \code{control.family}-list as control parameters to the \code{family="bgev"}
        ## :NAME: control.bgev.default
        list(
            ## :ARGUMENT: q.location The quantile level for the location parameter
            q.location = 0.5,

            ## :ARGUMENT: q.spread The quantile level for the spread parameter (must be < 0.5)
            q.spread = 0.25,

            ## :ARGUMENT: q.mix The lower and upper quantile level for the mixing function
            q.mix = c(0.10, 0.20),

            ## :ARGUMENT: beta.ab The parameters a and b in the Beta mixing function
            beta.ab = 5L
        )

        ## :SEEALSO: inla
    }

`inla.set.control.family.default` <-
    function(...) {
        ## :EXTRA:
        ## :NAME: control.family
        list(
            ## :ARGUMENT: dummy A dummy argument that can be used as a workaround
            dummy = 0,

            ## :ARGUMENT: hyper Definition of the hyperparameters
            hyper = NULL,

            ## :ARGUMENT+: initial (OBSOLETE!) Initial value for the hyperparameter(s) of the
            ## :ARGUMENT+: likelihood in the internal scale.
            initial = NULL,

            ## :ARGUMENT+: prior (OBSOLETE!) The name of the prior distribution(s) for othe
            ## :ARGUMENT+: hyperparameter(s).
            prior = NULL,

            ## :ARGUMENT: param (OBSOLETE!) The parameters for the prior distribution
            param = NULL,

            ## :ARGUMENT+: fixed (OBSOLETE!) Boolean variable(s) to say if the hyperparameter(s)
            ## :ARGUMENT+: is fixed or random.
            fixed = NULL,

            ## :ARGUMENT+: link (OBSOLETE! Use \code{control.link=list(model=)} instead.) The
            ## :ARGUMENT+: link function to use.
            link = "default",

            ## :ARGUMENT+: sn.shape.max Maximum value for the shape-parameter for Skew Normal
            ## :ARGUMENT+: observations (default 5.0)
            sn.shape.max = 5.0,

            ## :ARGUMENT+: gev.scale.xi (Expert option, do not use unless you know what you are
            ## :ARGUMENT+: doing.) The internal scaling of the shape-parameter for the GEV distribution.
            ## :ARGUMENT+: (default 0.1)
            gev.scale.xi = 0.1,

            ## :ARGUMENT: control.bgev  See \code{?control.bgev}
            control.bgev = NULL,

            ## :ARGUMENT: cenpoisson.I The censoring interval for the censored Poisson
            cenpoisson.I = c(-1L, -1L),

            ## :ARGUMENT+: beta.censor.value The censor value for the Beta-likelihood \code{(0
            ## :ARGUMENT+: <= beta.censor.value < 1/2)}
            beta.censor.value = 0.0,

            ## :ARGUMENT+: variant This variable is used to give options for various variants of
            ## :ARGUMENT+: the likelihood,  like chosing different parameterisations for example. See the
            ## :ARGUMENT+: relevant likelihood documentations for options (does only apply to some
            ## :ARGUMENT+: likelihoods).
            variant = 0L,

            ## :ARGUMENT: control.mix See \code{?control.mix}
            control.mix = NULL,

            ## :ARGUMENT: control.pom See \code{?control.pom}
            control.pom = NULL,

            ## :ARGUMENT: control.link See \code{?control.link}
            control.link = NULL,

            ## :ARGUMENT: link.simple See \code{inla.doc("0inflated")}
            link.simple = "default"
        )

        ## :SEEALSO: inla
    }

`inla.set.control.fixed.default` <-
    function(...) {
        ## :EXTRA:
        ## :NAME: control.fixed
        list(
            ## :ARGUMENT: cdf  A list of values to compute the CDF for, for all fixed effects
            cdf = NULL,

            ## :ARGUMENT: quantiles  A list of quantiles to compute for all fixed effects
            quantiles = NULL,

            ## :ARGUMENT+: expand.factor.strategy The strategy used to expand factors into fixed
            ## :ARGUMENT+: effects based on their levels. The default strategy is us use the
            ## :ARGUMENT+: \code{model.matrix}-function for which NA's are not allowed
            ## :ARGUMENT+: (\code{expand.factor.strategy="model.matrix"}) and levels are possible removed.
            ## :ARGUMENT+: The alternative option (\code{expand.factor.strategy="inla"}) use an
            ## :ARGUMENT+: \code{inla}-spesific expansion which expand a factor into one fixed effects for
            ## :ARGUMENT+: each level, do allow for NA's and all levels are present in the model. In this
            ## :ARGUMENT+: case, factors MUST BE factors in the data.frame/list and NOT added as
            ## :ARGUMENT+: \code{.+factor(x1)+.} in the formula only.
            expand.factor.strategy = "model.matrix",

            ## :ARGUMENT+: mean Prior mean for all fixed effects except the intercept.
            ## :ARGUMENT+:  Alternatively, a named list with specific means where name=default applies to
            ## :ARGUMENT+:  unmatched names. For example \code{control.fixed=list(mean=list(a=1, b=2,
            ## :ARGUMENT+:  default=0))} assign 'mean=1' to fixed effect 'a' , 'mean=2' to effect 'b' and
            ## :ARGUMENT+: 'mean=0' to all others. (default 0.0)
            mean = 0.0,

            ## :ARGUMENT: mean.intercept Prior mean for the intercept (default 0.0)
            mean.intercept = 0.0,

            ## :ARGUMENT+: prec  Default precision for all fixed effects except the intercept.
            ## :ARGUMENT+: Alternatively, a named list with specific means where name=default applies to
            ## :ARGUMENT+: unmatched names.  For example \code{control.fixed=list(prec=list(a=1, b=2,
            ## :ARGUMENT+: default=0.01))} assign 'prec=1' to fixed effect 'a' , 'prec=2' to effect 'b' and
            ## :ARGUMENT+: 'prec=0.01' to all others. (default 0.001)
            prec = 0.001,

            ## :ARGUMENT: prec.intercept  Default precision the intercept (default 0.0)
            prec.intercept = 0.0,

            ## :ARGUMENT: compute Compute marginals for the fixed effects ? (default TRUE)
            compute = TRUE,

            ## :ARGUMENT+: correlation.matrix Compute the posterior correlation matrix for all
            ## :ARGUMENT+:  fixed effects? (default FALSE) OOPS: This option will set up appropriate linear
            ## :ARGUMENT+:  combinations and the results are shown as the posterior correlation matrix of the
            ## :ARGUMENT+:  linear combinations. This option will imply
            ## :ARGUMENT+: \code{control.inla=list(lincomb.derived.correlation.matrix=TRUE)}.
            correlation.matrix = FALSE,

            ## :ARGUMENT+: remove.names A vector of names of expanded fixed effects to remove
            ## :ARGUMENT+: from the model-matrix. This is an expert option, and should only be used if you
            ## :ARGUMENT+: know what you are doing.
            remove.names = NULL
        )

        ## :SEEALSO: inla
    }

`inla.set.control.inla.default` <-
    function(...) {
        family <- "gaussian"
        xx <- list(...)[1]
        if (!is.null(xx$family)) {
            family <- xx$family
        }

        ## :EXTRA:
        ## :NAME: control.inla
        ans <- list(
            ## :ARGUMENT+: strategy  Character The strategy to use for the approximations; one
            ## :ARGUMENT+: of 'auto' (default), 'gaussian', 'simplified.laplace', 'laplace' or 'adaptive'.
            strategy = "auto",

            ## :ARGUMENT+: int.strategy  Character The integration strategy to use; one of
            ## :ARGUMENT+: 'auto' (default),  'ccd', 'grid', 'eb' (empirical bayes),  'user' or 'user.std'.
            ## :ARGUMENT+: For the experimental mode,  then 'grid' equal 'ccd' for more than two
            ## :ARGUMENT+: hyperparameters.
            int.strategy = "auto",

            ## :ARGUMENT+: int.design  Matrix Matrix of user-defined integration points and
            ## :ARGUMENT+: weights. Each row consists theta values and the integration weight.
            ## :ARGUMENT+: (EXPERIMENTAL!).
            int.design = NULL,

            ## :ARGUMENT+: interpolator  Character The interpolator used to compute the
            ## :ARGUMENT+: marginals for the hyperparameters. One of 'auto', 'nearest', 'quadratic',
            ## :ARGUMENT+: 'weighted.distance', 'ccd', 'ccdintegrate', 'gridsum', 'gaussian'. Default is
            ## :ARGUMENT+: 'auto'.
            interpolator = "auto",

            ## :ARGUMENT+: fast Logical If TRUE, then replace conditional modes in the Laplace
            ## :ARGUMENT+: approximation with conditional expectation (default TRUE).
            fast = TRUE,

            ## :ARGUMENT+: linear.correction Logical Default TRUE for the 'strategy = laplace'
            ## :ARGUMENT+: option.
            linear.correction = NULL,

            ## :ARGUMENT+: h Numerical The step-length for the gradient calculations for the
            ## :ARGUMENT+: hyperparameters. Default 0.005.
            h = 0.005,

            ## :ARGUMENT+: dz Numerical The step-length in the standarised scale for the
            ## :ARGUMENT+: integration of the hyperparameters. Default 0.75.
            dz = 0.75,

            ## :ARGUMENT+: diff.logdens Numerical The difference of the log.density for the
            ## :ARGUMENT+: hyperpameters to stop numerical integration using int.strategy='grid'. Default 6.
            diff.logdens = 6,

            ## :ARGUMENT+: print.joint.hyper Logical If TRUE, the store also the joint
            ## :ARGUMENT+: distribution of the hyperparameters (without any costs). Default TRUE.
            print.joint.hyper = TRUE,

            ## :ARGUMENT+: force.diagonal Logical If TRUE, then force the Hessian to be
            ## :ARGUMENT+: diagonal. (Default \code{FALSE})
            force.diagonal = FALSE,

            ## :ARGUMENT+: skip.configurations Logical Skip configurations if the values at the
            ## :ARGUMENT+: main axis are to small. (Default \code{TRUE})
            skip.configurations = TRUE,

            ## :ARGUMENT+: mode.known Logical If TRUE then no optimisation is done. (Default
            ## :ARGUMENT+: FALSE.)
            mode.known = FALSE,

            ## :ARGUMENT+: adjust.weights Logical If TRUE then just more accurate integration
            ## :ARGUMENT+: weights. (Default TRUE.)
            adjust.weights = TRUE,

            ## :ARGUMENT+: tolerance Numerical The tolerance for the optimisation of the
            ## :ARGUMENT+: hyperparameters. If set, this is the default value for for 'tolerance.f/2',
            ## :ARGUMENT+: 'tolerance.g' and  'tolerance.x'; see below.
            tolerance = 0.005,

            ## :ARGUMENT+: tolerance.f Numerical The tolerance for the absolute change in the
            ## :ARGUMENT+: log posterior in the optimisation of the hyperparameters.
            tolerance.f = NULL,

            ## :ARGUMENT+: tolerance.g Numerical The tolerance for the absolute change in the
            ## :ARGUMENT+: gradient of the log posterior in the optimisation of the hyperparameters.
            tolerance.g = NULL,

            ## :ARGUMENT+: tolerance.x Numerical The tolerance for the change in the
            ## :ARGUMENT+: hyperparameters (root-mean-square) in the optimisation of the hyperparameters.
            tolerance.x = NULL,

            ## :ARGUMENT+: tolerance.step Numerical The tolerance for the change in
            ## :ARGUMENT+: root-mean_squre in the inner Newton-like optimisation of the latent field.
            tolerance.step = 0.001,

            ## :ARGUMENT+: restart Numerical To improve the optimisation, the optimiser is
            ## :ARGUMENT+: restarted at the found optimum 'restart' number of times.
            restart = 0L,

            ## :ARGUMENT+: optimiser Character The optimiser to use; one of 'gsl' or 'default'.
            optimiser = "default",

            ## :ARGUMENT: verbose Logical Run in verbose mode? (Default FALSE)
            verbose = NULL,

            ## :ARGUMENT+: reordering Character Type of reordering to use. (EXPERT OPTION; one
            ## :ARGUMENT+: of "AUTO", "DEFAULT", "IDENTITY", "REVERSEIDENTITY",  "BAND", "METIS", "GENMMD",
            ## :ARGUMENT+: "AMD", "MD", "MMD", "AMDBAR", "AMDC", "AMDBARC",  or the output from
            ## :ARGUMENT+: \code{inla.qreordering}. Default is 'auto'.)
            reordering = "auto",

            ## :ARGUMENT+: cpo.diff Numerical Threshold to define when the cpo-calculations are
            ## :ARGUMENT+: inaccurate. (EXPERT OPTION.)
            cpo.diff = NULL,

            ## :ARGUMENT+: npoints Numerical Number of points to use in the 'stratey=laplace'
            ## :ARGUMENT+: approximation (default 9)
            npoints = 9,

            ## :ARGUMENT+: cutoff Numerical The cutoff used in the 'stratey=laplace'
            ## :ARGUMENT+: approximation. (Smaller value is more accurate and more slow.) (default 1e-4)
            cutoff = 1e-4,

            ## :ARGUMENT+: adapt.hessian.mode Logical Should optimisation be continued if the
            ## :ARGUMENT+: Hessian estimate is void? (Default TRUE)
            adapt.hessian.mode = NULL,

            ## :ARGUMENT+: adapt.hessian.max.trials Numerical Number of steps in the adaptive
            ## :ARGUMENT+: Hessian optimisation
            adapt.hessian.max.trials = NULL,

            ## :ARGUMENT+: adapt.hessian.scale Numerical The scaling of the 'h' after each
            ## :ARGUMENT+: trial.
            adapt.hessian.scale = NULL,

            ## :ARGUMENT+: adaptive.max Selecting \code{strategy="adaptive"} will chose the
            ## :ARGUMENT+: default strategy for all fixed effects and model components with length less or
            ## :ARGUMENT+: equal to \code{adaptive.max}, for others, the gaussian strategy will be applied.
            adaptive.max = 25L,

            ## :ARGUMENT+: huge Logical If TRUE then try to do some of the internal
            ## :ARGUMENT+: parallisations differently. Hopefully this will be of benefite for 'HUGE' models.
            ## :ARGUMENT+: (Default FALSE.) [THIS OPTION IS OBSOLETE AND NOT USED!]
            huge = FALSE,

            ## :ARGUMENT+: step.len Numerical The step-length used to compute numerical
            ## :ARGUMENT+: derivaties of the log-likelihood (0 means \code{default} which
            ## :ARGUMENT+: depends on \code{stencil})
            step.len = 0.0, 

            ## :ARGUMENT+: stencil Numerical Number of points in the stencil used to compute the
            ## :ARGUMENT+: numerical derivaties of the log-likelihood (5, 7 or 9). (default 5)
            stencil = 5L,

            ## :ARGUMENT+: lincomb.derived.correlation.matrix Logical If TRUE compute also the
            ## :ARGUMENT+: correlations for the derived linear combinations, if FALSE do not (Default FALSE)
            lincomb.derived.correlation.matrix = FALSE,

            ## :ARGUMENT+: diagonal Numerical Expert use only! Add a this value on the diagonal
            ## :ARGUMENT+: of the joint precision matrix. (default 0.0)
            diagonal = 0.0,

            ## :ARGUMENT+: numint.maxfeval Numerical Maximum number of function evaluations in
            ## :ARGUMENT+: the the numerical integration for the hyperparameters. (Default 100000.)
            numint.maxfeval = 100000,

            ## :ARGUMENT+: numint.relerr Numerical Relative error requirement in the the
            ## :ARGUMENT+: numerical integration for the hyperparameters. (Default 1e-5)
            numint.relerr = 1e-5,

            ## :ARGUMENT+: numint.abserr Numerical Absolute error requirement in the the
            ## :ARGUMENT+: numerical integration for the hyperparameters. (Default 1e-6)
            numint.abserr = 1e-6,

            ## :ARGUMENT+: cmin Numerical The minimum value for the negative Hessian from the
            ## :ARGUMENT+: likelihood. Increasing this value will stabalise the optimisation but can
            ## :ARGUMENT+: introduce bias.  (Default -Inf)
            cmin = -Inf,

            ## :ARGUMENT+: b.strategy Character If \code{cmin} is used, either keep the linear
            ## :ARGUMENT+:     term (with \code{b.strategy="keep"}) or skip the contribution by
            ## :ARGUMENT+:     setting the linear term to zero (\code{b.strategy="skip"}). The
            ## :ARGUMENT+:     default value is \code{"skip"}
            b.strategy = "skip",

            ## :ARGUMENT+: step.factor Numerical The step factor in the Newton-Raphson algorithm
            ## :ARGUMENT+: saying how large step to take (Default 1.0) YES! setting this to a
            ## :ARGUMENT+: negative values means = 1,  EXCEPT the first time (for each thread)
            ## :ARGUMENT+: where |step.factor| is used. 
            step.factor = -0.1,

            ## :ARGUMENT+: global.node.factor Numerical The factor which defines the degree
            ## :ARGUMENT+:     required (how many neighbors), as a fraction of \code{n-1},
            ## :ARGUMENT+:     that is required to be classified as a global node and numbered
            ## :ARGUMENT+:     last (whatever the reordering routine says). Here, \code{n},
            ## :ARGUMENT+:     is the size of the graph. (Disabled if larger than 1,  default 2)
            global.node.factor = 2.0,

            ## :ARGUMENT+: global.node.degree Numerical The degree required (number of
            ## :ARGUMENT+:     neighbors) to be classified as a global node and numbered
            ## :ARGUMENT+:     last (whatever the reordering routine says).
            ## :ARGUMENT+:     (default \code{.Machine$integer.max})
            global.node.degree = .Machine$integer.max,

            ## :ARGUMENT+: stupid.search Logical Enable or disable the stupid-search-algorithm,
            ## :ARGUMENT+: if the Hessian calculations reveals that the mode is not found.
            ## :ARGUMENT+: (Default \code{TRUE}.)
            stupid.search = TRUE,

            ## :ARGUMENT+: stupid.search.max.iter Numerical Maximum number of iterations allowed
            ## :ARGUMENT+: for the stupid-search-algorithm. (default 1000)
            stupid.search.max.iter = 1000L,

            ## :ARGUMENT+: stupid.search.factor Numerical Factor (>=1) to increase the
            ## :ARGUMENT+: step-length with after each new interation. (default 1.05)
            stupid.search.factor = 1.05,

            ## :ARGUMENT+: control.vb List of arguments for various VB corrections.
            ## :ARGUMENT+: \code{enable}: Logical/Character Use this feature? If \code{"auto"} this will be selected automatically. 
            ## :ARGUMENT+: \code{strategy}: Charactor What to correct, either "mean" or "variance".
            ## :ARGUMENT+: \code{verbose}: Logical Be verbose or not. 
            ## :ARGUMENT+: \code{iter.max}: Integer Maximum number of interations. 
            ## :ARGUMENT+: \code{emergency}: Numeric If the standarized correction for the mean is larger
            ## :ARGUMENT+:           than this value, then call the vb.correction off and issue a warning
            ## :ARGUMENT+: \code{f.enable.limit}: Vector of length 2. The size limit to correct for a \code{f()}.
            ## :ARGUMENT+:         First element is for \code{strategy="mean"}.
            ## :ARGUMENT+:         Second element is for \code{strategy="mean"}.
            ## :ARGUMENT+: \code{hesssian.update}: How many times the Hessian is updated for each
            ## :ARGUMENT+:         correction (\code{strategy="variance"} only).
            ## :ARGUMENT+: \code{hesssian.strategy}: Select strategy for computing the Hessian
            ## :ARGUMENT+:         matrix for \code{strategy="variance"}, one of \code{"full"}, 
            ## :ARGUMENT+:         \code{"diagonal"}, \code{"partial"} and \code{"default"}.
            control.vb = list(enable = "auto", strategy = c("mean", "variance"), verbose = TRUE,
                              iter.max = 25, emergency = 25, f.enable.limit = c(30, 25), hessian.update = 2,
                              hessian.strategy = c("default", "full", "partial", "diagonal")), 

            ## :ARGUMENT+: num.gradient Character Set the numerical scheme to compute the
            ## :ARGUMENT+: gradient,  one of \code{"forward"} or \code{"central"} (default).
            num.gradient = "central",

            ## :ARGUMENT+: num.hessian Character Set the numerical scheme to compute the
            ## :ARGUMENT+: Hessian,  one of \code{"forward"} or \code{"central"} (default).
            num.hessian = "central",

            ## :ARGUMENT+: optimise.strategy Character THIS OPTION IS EXPERIMENTAL. Chose the
            ## :ARGUMENT+: optimiser strategy,  one of \code{"plain"} or \code{"smart"} (default)
            optimise.strategy = "smart",

            ## :ARGUMENT+: use.directions THIS OPTION IS EXPERIMENTAL. Unless \code{FALSE} or
            ## :ARGUMENT+: \code{NULL},  use directions for computing gradient and Hessian, initialised with
            ## :ARGUMENT+: \code{use.directions} if a matrix.
            use.directions = TRUE,

            ## :ARGUMENT+: constr.marginal.diagonal Add stability to \code{AQ^-1A^T} by adding a
            ## :ARGUMENT+: small diagonal term. (default \code{epsilon^0.5})
            constr.marginal.diagonal = sqrt(.Machine$double.eps),

            ## :ARGUMENT+: improved.simplified.laplace If \code{TRUE} use an experimental
            ## :ARGUMENT+: improved variant, otherwise, use the standard one.
            improved.simplified.laplace = FALSE,

            ## :ARGUMENT+: parallel.linesearch Use serial (default) or parallel line-search
            ## :ARGUMENT+: (highly experimental for the moment)
            parallel.linesearch = FALSE,
            
            ## :ARGUMENT+: compute.initial.values Compute initial values for the latent field or not.
            ## :ARGUMENT+: (experimental-mode only)
            compute.initial.values = TRUE,

            ## :ARGUMENT+: hessian.correct.skewness.only If TRUE then correct only
            ## :ARGUMENT+: skewness in the Hessian, for the hyperparameters. If FALSE (default),
            ## :ARGUMENT+: correct also variance (experimental-mode only)
            hessian.correct.skewness.only = FALSE
        )

        ## :SEEALSO: inla

        return(ans)
    }


`inla.set.control.predictor.default` <-
    function(...) {
        ## :EXTRA:
        ## :NAME: control.predictor
        list(
            ## :ARGUMENT: hyper Definition of the hyperparameters.
            hyper = NULL,

            ## :ARGUMENT+: fixed (OBSOLETE!) If the precision for the artificial noise is fixed
            ## :ARGUMENT+: or not (defualt TRUE)
            fixed = NULL,

            ## :ARGUMENT: prior (OBSOLETE!) The prior for the artificial noise
            prior = NULL,

            ## :ARGUMENT: param (OBSOLETE!) Prior parameters for the artificial noise
            param = NULL,

            ## :ARGUMENT+: initial (OBSOLETE!) The value of the log precision of the artificial
            ## :ARGUMENT+: noise
            initial = NULL,

            ## :ARGUMENT+: compute A boolean variable; should the marginals for the linear
            ## :ARGUMENT+: predictor be computed? (Default FALSE.)
            compute = FALSE,

            ## :ARGUMENT: cdf A list of values to compute the CDF for the linear predictor
            cdf = NULL,

            ## :ARGUMENT: quantiles A list of quantiles to compute for the linear predictor
            quantiles = NULL,

            ## :ARGUMENT+: cross Cross-sum-to-zero constraints with the linear predictor. All
            ## :ARGUMENT+: linear predictors with the same level of 'cross' are constrained to have sum
            ## :ARGUMENT+: zero. Use 'NA' for no contribution. 'Cross' has the same length as the linear
            ## :ARGUMENT+: predictor (including the 'A' matrix extention). (THIS IS AN EXPERIMENTAL OPTION,
            ## :ARGUMENT+: CHANGES MAY APPEAR.)
            cross = NULL,

            ## :ARGUMENT: A The observation matrix (matrix or Matrix::sparseMatrix).
            A = NULL,

            ## :ARGUMENT: precision The precision for eta* - A*eta, (default \code{exp(15)})
            precision = exp(15),

            ## :ARGUMENT+: link Define the family-connection for unobserved observations
            ## :ARGUMENT+: (\code{NA}). \code{link} is integer values which defines the family connection;
            ## :ARGUMENT+: \code{family[link[idx]]} unless \code{is.na(link[idx])} for which the
            ## :ARGUMENT+: identity-link is used. The \code{link}-argument only influence the
            ## :ARGUMENT+: \code{fitted.values} in the \code{result}-object. If \code{is.null(link)}
            ## :ARGUMENT+: (default) then the identity-link is used for all missing observations. If the
            ## :ARGUMENT+: length of \code{link} is 1, then this value is replicated with the length of the
            ## :ARGUMENT+: responce vector. If an element of the responce vector is \code{!NA} then the
            ## :ARGUMENT+: corresponding entry in \code{link} is not used (but must still be a legal value).
            ## :ARGUMENT+: Setting this variable implies \code{compute=TRUE}.
            link = NULL
        )

        ## :SEEALSO: inla
    }

`inla.set.control.mode.default` <-
    function(...) {
        ## this is internal use only...
        ## :EXTRA:
        ## :NAME: control.mode
        list(
            ## :ARGUMENT+: result Prevous result from inla(). Use the theta- and x-mode from
            ## :ARGUMENT+: this run.
            result = NULL,

            ## :ARGUMENT+: theta The theta-mode/initial values for theta. This option has
            ## :ARGUMENT+: preference over result$mode$theta.
            theta = NULL,

            ## :ARGUMENT+: x The x-mode/intitial values for x. This option has preference over
            ## :ARGUMENT+: result$mode$x.
            x = NULL,

            ## :ARGUMENT+: restart A boolean variable; should we restart the optimisation from
            ## :ARGUMENT+: this configuration or fix the mode at this configuration? (Default FALSE.)
            restart = FALSE,

            ## :ARGUMENT+: fixed A boolean variable. If TRUE then treat all thetas as known and
            ## :ARGUMENT+: fixed, and if FALSE then treat all thetas as unknown and random (default).
            fixed = FALSE
        )

        ## :SEEALSO: inla
    }

`inla.set.control.hazard.default` <-
    function(...) {
        ## :EXTRA:
        ## :NAME: control.hazard
        list(
            ## :ARGUMENT+: model The model for the baseline hazard model. One of 'rw1', 'rw2' or
            ## :ARGUMENT+: 'iid'. (Default 'rw1'.)
            model = "rw1",

            ## :ARGUMENT: hyper The definition of the hyperparameters.
            hyper = NULL,

            ## :ARGUMENT+: fixed (OBSOLETE!) A boolean variable; is the precision for 'model'
            ## :ARGUMENT+: fixed? (Default FALSE.)
            fixed = FALSE,

            ## :ARGUMENT: initial (OBSOLETE!) The initial value for the precision.
            initial = NULL,

            ## :ARGUMENT+: prior (OBSOLETE!) The prior distribution for the precision for
            ## :ARGUMENT+: 'model'
            prior = NULL,

            ## :ARGUMENT: param (OBSOLETE!) The parameters in the prior distribution
            param = NULL,

            ## :ARGUMENT+: constr A boolean variable; shall the  'model' be constrained to sum
            ## :ARGUMENT+: to zero?
            constr = TRUE,

            ## :ARGUMENT+: diagonal An extra constant added to the diagonal of the precision
            ## :ARGUMENT+: matrix
            diagonal = NULL,

            ## :ARGUMENT: n.intervals Number of intervals in the baseline hazard. (Default 15)
            n.intervals = 15,

            ## :ARGUMENT+: cutpoints The cutpoints to use. If not specified the they are compute
            ## :ARGUMENT+: from 'n.intervals' and the maximum length of the interval. (Default NULL)
            cutpoints = NULL,

            ## :ARGUMENT+: strata.name The name of the stratefication variable for the baseline
            ## :ARGUMENT+: hazard in the data.frame
            strata.name = NULL,

            ## :ARGUMENT+: scale.model Scale the baseline hazard model (RW1, RW2) so the
            ## :ARGUMENT+: generalized variance is 1. (Default
            ## :ARGUMENT+: \code{inla.getOption("scale.model.default")}.)
            scale.model = NULL
        )

        ## :SEEALSO: inla
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

inla.make.completion.function <- function(...) {
    my.eval <- function(command, envir = parent.frame(),
                        enclos = if (is.list(envir) || is.pairlist(envir)) parent.frame() else baseenv()) {
        return(eval(parse(text = command), envir, enclos))
    }

    xx <- sort(list(...)[[1L]])
    my.eval(paste("function(", paste(xx, sep = "", collapse = ", "), ") {
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

if (TRUE) {
    control.update <- inla.make.completion.function(names(inla.set.control.update.default()))
    control.lincomb <- inla.make.completion.function(names(inla.set.control.lincomb.default()))
    control.group <- inla.make.completion.function(names(inla.set.control.group.default()))
    control.mix <- inla.make.completion.function(names(inla.set.control.mix.default()))
    control.pom <- inla.make.completion.function(names(inla.set.control.pom.default()))
    control.link <- inla.make.completion.function(names(inla.set.control.link.default()))
    control.expert <- inla.make.completion.function(names(inla.set.control.expert.default()))
    control.compute <- inla.make.completion.function(names(inla.set.control.compute.default()))
    control.family <- inla.make.completion.function(names(inla.set.control.family.default()))
    control.fixed <- inla.make.completion.function(names(inla.set.control.fixed.default()))
    control.inla <- inla.make.completion.function(names(inla.set.control.inla.default()))
    control.vb <- inla.make.completion.function(names(inla.set.control.inla.default()$control.vb))
    control.gcpo <- inla.make.completion.function(names(inla.set.control.compute.default()$control.gcpo))
    control.predictor <- inla.make.completion.function(names(inla.set.control.predictor.default()))
    control.mode <- inla.make.completion.function(names(inla.set.control.mode.default()))
    control.hazard <- inla.make.completion.function(names(inla.set.control.hazard.default()))
    control.bgev <- inla.make.completion.function(names(inla.set.control.bgev.default()))
    control.lp.scale <- inla.make.completion.function(names(inla.set.control.lp.scale.default()))
    control.pardiso <- inla.make.completion.function(names(inla.set.control.pardiso.default()))
    control.scopy <- inla.make.completion.function(names(inla.set.control.scopy.default()))
} else {
    control.update <- NULL
    control.lincomb <-NULL
    control.group <- NULL
    control.mix <- NULL
    control.pom <- NULL
    control.link <- NULL
    control.expert <- NULL
    control.compute <-NULL
    control.family <- NULL
    control.fixed <- NULL
    control.inla <- NULL
    control.vb <- NULL
    control.gcpo <- NULL
    control.predictor <- NULL
    control.mode <- NULL
    control.hazard <- NULL
    control.bgev <- NULL
    control.lp.scale <- NULL
    control.pardiso <- NULL
    control.scopy <- NULL
}
