### Defines default arguments

`inla.set.control.lincomb.default` =
    function()
{
    ##:NAME: control.lincomb
    list(
         ##:ARGUMENT: precision The precision for the artificial tiny noise. Default 1e09.
         precision = NULL, 

         ##:ARGUMENT: usermap One of the predefined functions to transform the linear combinations with. Default NULL.
         usermap = NULL,

         ##:ARGUMENT: verbose Use verbose mode for linear combinations if verbose model is set globally. (Default TRUE)
         verbose = TRUE)

    ##:SEEALSO: inla
}

`inla.set.control.group.default` =
    function()
{
    ##:NAME: control.group
    list(
         ##:ARGUMENT: initial The initial value for the group correlation in the internal scale.
         initial = NULL, 

         ##:ARGUMENT: fixed A boolean variable if the group correction is assumed to be fixed or random.
         fixed = NULL,

         ##:ARGUMENT: prior The name of the prior distribution for the group correlation in the internal scale
         prior = NULL, 

         ##:ARGUMENT: param Prior parameters
         param = NULL,

         ##:ARGUMENT: model Group model (one of 'exchangable' or 'ar1')
         model = NULL)
    ##:SEEALSO: inla
}

`inla.set.f.default` =
    function()
{
    list(diagonal = 1e-5)
}

`inla.set.control.expert.default` =
    function()
{
    ##:NAME: control.expert
    list(
         ##:ARGUMENT: cpo.manual A boolean variable to decide if the inla-program is to be runned in a manual-cpo-mode. (EXPERT OPTION.)
         cpo.manual = FALSE, 

         ##:ARGUMENT: cpo.idx  The index of the data point to remove. (EXPERT OPTION.)
         cpo.idx = -1)

    ##:SEEALSO: inla
}


`inla.set.control.compute.default`=
    function()
{
    ##:NAME: control.compute
    list(
         ##:ARGUMENT: hyperpar A boolean variable if the marginal for the hyperparameters should be computed. Default TRUE.
         hyperpar=TRUE,

         ##:ARGUMENT: dic A boolean variable if the DIC-value should be computed. Default FALSE.
         dic=FALSE,

         ##:ARGUMENT: mlik A boolean variable if the marginal likelihood should be computed. Default FALSE.
         mlik=TRUE,

         ##:ARGUMENT: cpo A boolean variable if the cross-validated predictive measures (cpo, pit) should be computed
         cpo=FALSE,

         ##:ARGUMENT: q A boolean variable if binary images of the precision matrix, the reordered precision matrix and the Cholesky triangle should be generated. Default FALSE.
         q=FALSE,

         ##:ARGUMENT: smtp The sparse-matrix solver, one of 'smtp' (default) or 'band'
         smtp = NULL)

    ##:SEEALSO: inla
}

`inla.set.control.data.default`=
    function()
{
    ##:NAME: control.data
    list(
         ##:ARGUMENT: initial Initial value for the hyperparameter(s) of the likelihood in the internal scale.
         initial=NULL,

         ##:ARGUMENT: prior The name of the prior distribution(s) for othe hyperparameter(s).
         prior=NULL,

         ##:ARGUMENT: param The parameters for the prior distribution
         param=NULL,

         ##:ARGUMENT: fixed Boolean variable(s) to say if the hyperparameter(s) is fixed or random.
         fixed=NULL,

         ##:ARGUMENT: dof.max Maximum degrees of freedom for the student-t distribution before its is considered as Normal (default 100)
         dof.max=NULL,

         ##:ARGUMENT: alpha The parameter 'alpha' for the asymmetric Laplace likelihood  (default 0.5)
         alpha=0.5,

         ##:ARGUMENT: epsilon The parameter 'epsilon' for the asymmetric Laplace likelihood (default 0.01)
         epsilon = 0.01,

         ##:ARGUMENT: gamma The parameter 'gamma' for the asymmetric Laplace likelihood (default 1.0)
         gamma = 1.0,

         ##:ARGUMENT: sh.shape.max Maximum value for the shape-parameter for Skew Normal observations
         sn.shape.max = 5.0
         )

    ##:SEEALSO: inla
}

`inla.set.control.fixed.default`=
    function()
{   
    ##:NAME: control.fixed
    list(
         ##:ARGUMENT: cdf  A list of values to compute the CDF for, for all fixed effects
         cdf=NULL,

         ##:ARGUMENT: quantiles  A list of quantiles to compute for all fixed effects
         quantiles = NULL,

         ##:ARGUMENT: mean Prior mean for all fixed effects except the intercept. Alternatively, a named list with spesific means where name=default applies to unmatched names. For example \code{control.fixed=list(mean=list(a=1,b=2,default=0))} assign 'mean=1' to fixed effect 'a' , 'mean=2' to effect 'b' and 'mean=0' to all others.  
         mean = 0.0,
         
         ##:ARGUMENT: mean.intercept Prior mean for the intercept
         mean.intercept = 0.0,
         
         ##:ARGUMENT: prec  Default precision for all fixed effects except the intercept. Alternatively, a named list with spesific means where name=default applies to unmatched names.  For example \code{control.fixed=list(prec=list(a=1,b=2,default=0.01))} assign 'prec=1' to fixed effect 'a' , 'prec=2' to effect 'b' and 'prec=0.01' to all others.  
         prec=NULL, 

         ##:ARGUMENT: prec.intercept  Default precision the intercept (default 0.0)
         prec.intercept = 0.0)

    ##:SEEALSO: inla
}

`inla.set.control.inla.default`=
    function(family = "gaussian")
{
    ##:NAME: control.inla
    ans = list(
            ##:ARGUMENT: strategy  The strategy to use for the approximations; one of 'gaussian', 'simplified.laplace' (default) or 'laplace'
            strategy="simplified.laplace",

            ##:ARGUMENT: int.strategy  The integration strategy to use; one of 'ccd' (default), 'grid' or 'eb' (empirical bayes)
            int.strategy="ccd",

            ##:ARGUMENT: interpolator  The interpolator used to compute the marginals for the hyperparameters. One of 'auto', 'nearest', 'quadratic', 'weighted.distance', 'ccd', 'gridsum', 'gaussian'. Default is 'auto'.
            interpolator="auto",

            ##:ARGUMENT: linear.correction  Default TRUE for the 'strategy = laplace' option.
            linear.correction=NULL,

            ##:ARGUMENT: h The step-length for the gradient calculations for the hyperparameters. Default 0.01.
            h=0.01,

            ##:ARGUMENT: dz The step-length in the standarised scale for the integration of the hyperparameters. Default 1.0.
            dz=1.0,

            ##:ARGUMENT: diff.logdens The difference of the log.density for the hyperpameters to stop numerical integration using int.strategy='grid'. Default 2.5.
            diff.logdens=2.5,

            ##:ARGUMENT: print.joint.hyper If TRUE, the store also the joint distribution of the hyperparameters (without any costs). Default TRUE.
            print.joint.hyper=TRUE,

            ##:ARGUMENT: force.diagonal A boolean variable, if TRUE, then force the Hessian to be diagonal. (Default FALSE.)
            force.diagonal=FALSE,

            ##:ARGUMENT: skip.configurations A boolean variable; skip configurations if the values at the main axis are to small. (Default TRUE.)
            skip.configurations=TRUE,

            ##:ARGUMENT: mode.known A boolean variable: If TRUE then no optimisation is done. (Default FALSE.)
            mode.known=FALSE,

            ##:ARGUMENT: adjust.weights A boolean variable; If TRUE then just more accurate integration weights. (Default TRUE.)
            adjust.weights=TRUE,

            ##:ARGUMENT: tolerance The tolerance for the optimisation.
            tolerance = NULL,

            ##:ARGUMENT: restart To improve the optimisation, the optimiser is restarted at the found optimum 'restart' number of times.
            restart = NULL,

            ##:ARGUMENT: optimiser The optimiser to use; one of 'gsl' or 'domin'
            optimiser = NULL,

            ##:ARGUMENT: verbose A boolean variable; run in verbose mode? (Default FALSE)
            verbose = NULL,

            ##:ARGUMENT: reordering Type of reordering to use. (EXPERT OPTION; one of "AUTO", "DEFAULT", "IDENTITY", "BAND", "METIS", "GENMMD", "AMD", "MD")
            reordering = NULL,

            ##:ARGUMENT: cpo.diff Threshold to define when the cpo-calculations are inaccurate. (EXPERT OPTION.)
            cpo.diff = NULL,

            ##:ARGUMENT: npoints Number of points to use in the 'stratey=laplace' approximation
            npoints = NULL,

            ##:ARGUMENT: adapt.hessian.mode A boolean variable; should optimisation be continued if the Hessian estimate is void? (Default TRUE)
            adapt.hessian.mode = NULL,

            ##:ARGUMENT: adapt.hessian.max.trails Number of steps in the adaptive Hessian optimisation
            adapt.hessian.max.trails = NULL,

            ##:ARGUMENT: adapt.hessian.scale The scaling of the 'h' after each trial.
            adapt.hessian.scale = NULL,

            ##:ARGUMENT: huge A boolean variable; if TRUE then try to do some of the internal parallisations differently. Hopefully this will be of benefite for 'HUGE' models. (Default FALSE.)
            huge = FALSE,

            ##:ARGUMENT: step.len The step-length used to compute numerical derivaties of the log-likelihood
            step.len = NULL,

            ##:ARGUMENT: derived.only A boolen variable: if TRUE the only compute the marginals for the derived linear combinations and if FALSE, the and also the linear combinations to the graph (Default TRUE)
            derived.only = TRUE,

            ## NOT DOCUMENTED ONCES (expert only)
            diagonal = NULL,

            ##:ARGUMENT: numint.maxfeval Maximum number of function evaluations in the the numerical integration for the hyperparameters. (Default 10000.)
            numint.maxfeval = 10000,
            
            ##:ARGUMENT: numint.relerr Relative error requirement in the the numerical integration for the hyperparameters. (Default 1e-3)
            numint.relerr = 1e-3,

            ##:ARGUMENT: numint.abserr Absolute error requirement in the the numerical integration for the hyperparameters. (Default 1e-4)
            numint.abserr = 1e-4)

    ## use default Gaussian strategy if the observations are gaussian    
    if (all(inla.strcasecmp(family, "gaussian")))
        ans$strategy = "gaussian"

    ##:SEEALSO: inla

    return (ans)
}

`inla.set.control.predictor.default`=
    function()
{
    ##:NAME: control.predictor
    list(
         ##:ARGUMENT: compute A boolean variable; should the marginals for the linear predictor be computed? (Default FALSE.)
         compute=FALSE,

         ##:ARGUMENT: cdf A list of values to compute the CDF for the linear predictor
         cdf=NULL,

         ##:ARGUMENT: quantiles A list of quantiles to compute for the linear predictor
         quantiles = NULL,

         ##:ARGUMENT: fixed If the precision for the artificial noise is fixed or not (defualt TRUE)
         fixed=TRUE,

         ##:ARGUMENT: prior The prior for the artificial noise
         prior=NULL,

         ##:ARGUMENT: param Prior parameters for the artificial noise
         param=NULL,

         ##:ARGUMENT: initial The value of the log precision of the artificial noise
         initial=NULL,

         ##:ARGUMENT: cross Cross-sum-to-zero constraints
         cross=NULL,
         
         ##:ARGUMENT: predictor.usermap One of the predfined mappings to compute new marginals for
         predictor.usermap=NULL,

         ##:ARGUMENT: A The observation matrix A = list(i=..., j=..., values=...) or a filename.
         A = NULL,

         ##:ARGUMENT precision The precision for eta - A*eta,
         precision = exp(8.0))
         
    ##:SEEALSO: inla
}

`inla.set.control.results.default`=
    function()
{
    ##:NAME: control.results
    list(
         ##:ARGUMENT: return.marginals.random A boolean variable; read the marginals for the fterms? (Default TRUE)
         return.marginals.random=TRUE,

         ##:ARGUMENT: return.marginals.predictor A boolean variable; read the marginals for the linear predictor? (Default TRUE)
         return.marginals.predictor=TRUE,

         ##:ARGUMENT: image.dim The dimension to reduce the matrix-images returned by the argument 'control.compute(q=TRUE)' when displaying the images
         image.dim = 256)
    ##:SEEALSO: inla
}

`inla.set.control.mode.default`=
    function()
{
    ## this is internal use only...
    ##:NAME: control.mode
    list(
         ##:ARGUMENT: result Prevous result from inla(). Use the theta- and x-mode from this run.
         result = NULL,
         
         ##:ARGUMENT: theta The theta-mode/initial values for theta.
         theta = NULL,

         ##:ARGUMENT: x The x-mode/intitial values for x.
         x = NULL,

         ##:ARGUMENT: restart A boolean variable; should we restart the optimisation from this configuration or fix the mode at this configuration? (Default FALSE.)
         restart = FALSE)
    ##:SEEALSO: inla
}

`inla.set.control.hazard.default` = 
    function()
{
    ##:NAME: control.hazard
    list(
         ##:ARGUMENT: model The model for the baseline hazard model. One of 'rw1' or 'rw2'. (Default 'rw1'.)
         model = "rw1",

         ##:ARGUMENT: fixed A boolean variable; is the precision for 'model' fixed? (Default FALSE.)
         fixed = FALSE,
         
         ##:ARGUMENT: initial The initial value for the precision.
         initial = NULL,

         ##:ARGUMENT: constr A boolean variable; shall the  'model' be constrained to sum to zero?
         constr = TRUE,

         ##:ARGUMENT: prior The prior distribution for the precision for 'model'
         prior = NULL,

         ##:ARGUMENT: param The parameters in the prior distribution
         param = NULL,

         ##:ARGUMENT: n.intervals Number of intervals in the baseline hazard. (Default 15)
         n.intervals = 15,

         ##:ARGUMENT: cutpoints The cutpoints to use. If not specified the they are compute from 'n.intervals' and the maximum length of the interval. (Default NULL)
         cutpoints = NULL,

         ##:ARGUMENT: strata.name The name of the stratefication variable for the baseline hazard in the data.frame
         strata.name = NULL,
         
         ##:ARGUMENT: si A boolean variable; should all the Gaussian approximations configurations be written to files?
         si = FALSE)
    ##:SEEALSO: inla
}


## check control-arguments

`inla.check.control` = function(contr)
{
    ## This function will signal an error if the arguments in CONTR
    ## does not match the ones in the corresponding
    ## `inla.set.XX.default()' routine.

    ## EX: contr is `control.inla' and default arguments is found in
    ## `inla.set.control.inla.default()'

    if (!is.list(contr))
        stop("argument is not an list")
    if (length(contr) == 0)
        return(invisible())
    
    nm = paste(sys.call()[2])
    f = paste("inla.set.", nm, ".default()", sep="")
    elms = names(inla.eval(f))

    if (is.null(names(contr)))
        stop(inla.paste(c("Named elements in in control-argument `", nm, "', is required: ", contr,
                          "\n\n  Valid ones are:\n\t",
                          inla.paste(sort(elms),sep="\n\t")), sep=""))
    
    for(elm in names(contr))
        if (!is.element(elm, elms))
            stop(inla.paste(c("Name `", elm,"' in control-argument `", nm, "', is void.\n\n  Valid ones are:\n\t",
                              inla.paste(sort(elms),sep="\n\t")), sep=""))

    return(invisible())
}
