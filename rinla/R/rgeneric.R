## Export: inla.rgeneric.ar1.model 
## Export: inla.rgeneric.iid.model 
## Export: inla.rgeneric.define
## Export: inla.rgeneric.wrapper
## Export: inla.rgeneric.q

##!\name{rgeneric.define}
##!\alias{rgeneric}
##!\alias{rgeneric.define}
##!\alias{inla.rgeneric.define}
##!\alias{rgeneric.ar1.model}
##!\alias{inla.rgeneric.ar1.model}
##!\alias{rgeneric.iid.model}
##!\alias{inla.rgeneric.iid.model}
##!\alias{rgeneric.wrapper}
##!\alias{inla.rgeneric.wrapper}
##!\alias{rgeneric.q}
##!\alias{inla.rgeneric.q}
##!
##!\title{rgeneric models}
##!
##!\description{A framework for defining latent models in R}
##!
##!\usage{
##!inla.rgeneric.define(model = NULL, debug = FALSE, ...)
##!inla.rgeneric.iid.model(
##!        cmd = c("graph", "Q", "mu", "initial", "log.norm.const", "log.prior", "quit"),
##!        theta = NULL)
##!inla.rgeneric.ar1.model(
##!        cmd = c("graph", "Q", "mu", "initial", "log.norm.const", "log.prior", "quit"),
##!        theta = NULL)
##!inla.rgeneric.wrapper(
##!        cmd = c("graph", "Q", "mu", "initial", "log.norm.const", "log.prior", "quit"),
##!        model, theta = NULL)
##!inla.rgeneric.q(
##!        rmodel, 
##!        cmd = c("graph", "Q", "mu", "initial", "log.norm.const", "log.prior", "quit"),
##!        theta = NULL)
##!}
##!
##!\arguments{
##!  \item{model}{The definition of the model; see \code{inla.rgeneric.ar1.model}}
##!  \item{rmodel}{The rgeneric model-object, the output of \code{inla.rgeneric.define}}
##!  \item{debug}{Logical. Turn on/off debugging}
##!  \item{cmd}{An allowed request}
##!  \item{theta}{Values of theta}
##!  \item{...}{Named list of variables that defines the environment of \code{model}}
##!  \item{debug}{Logical. Enable debug output}
##!}
##!
##!\value{%%
##!  This allows a latent model to be 
##!  defined in \code{R}.
##!  See \code{inla.rgeneric.ar1.model} and
##!  \code{inla.rgeneric.iid.model} and the documentation for 
##!  worked out examples of how to define latent models in  this way.
##!  This will be somewhat slow and is intended for special cases and
##!  protyping. The function \code{inla.rgeneric.wrapper} is for
##!  internal use only.}
##!\author{Havard Rue \email{hrue@r-inla.org}}


`inla.rgeneric.ar1.model` = function(
    cmd = c("graph", "Q", "mu", "initial", "log.norm.const", "log.prior", "quit"),
    theta = NULL)
{
    ## this is an example of the 'rgeneric' model. here we implement
    ## the AR-1 model as described in inla.doc("ar1"), where 'rho' is
    ## the lag-1 correlation and 'prec' is the *marginal* (not
    ## conditional) precision.
    
    ## variables defined the in the define-call, are stored here
    ## (which is in the path)
    envir = parent.env(environment())

    interpret.theta = function()
    {
        ## internal helper-function to map the parameters from the internal-scale to the
        ## user-scale
        return (list(prec = exp(theta[1L]),
                     rho = 2*exp(theta[2L])/(1+exp(theta[2L])) - 1.0))
    }

    graph = function()
    {
        if (TRUE) {
            ## we can also use this easy solution, since we know that Q[i, j] is not 0 by
            ## accident... this require that 'theta' is set; see 'theta = initial()' below.
            G = Q()
        } else {
            ## return the graph of the model. the values of Q is only interpreted as zero or
            ## non-zero. return a sparse.matrix
            if (FALSE) {
                ## slow and easy: dense-matrices
                G = toeplitz(c(1, 1, rep(0, n-2L)))
                G = inla.as.sparse(G)
            } else {
                ## faster. we only need to define the upper-triangular of G
                i = c(
                    ## diagonal
                    1L, n, 2L:(n-1L),
                    ## off-diagonal
                    1L:(n-1L))
                j = c(
                    ## diagonal
                    1L, n, 2L:(n-1L),
                    ## off-diagonal
                    2L:n)
                x = 1 ## meaning that all are 1
                G = sparseMatrix(i=i, j=j, x=x, giveCsparse = FALSE)
            }
        }
        return (G)
    }

    Q = function()
    {
        ## returns the precision matrix for given parameters
        param = interpret.theta()
        if (FALSE) {
            ## slow and easy: dense-matrices
            Q = param$prec/(1-param$rho^2) * toeplitz(c(1+param$rho^2, -param$rho, rep(0, n-2L)))
            Q[1, 1] = Q[n, n] = param$prec/(1-param$rho^2)
            Q = inla.as.sparse(Q)
        } else {
            ## faster. we only need to define the upper-triangular Q!
            i = c(
                ## diagonal
                1L, n, 2L:(n-1L),
                ## off-diagonal
                1L:(n-1L))
            j = c(
                ## diagonal
                1L, n, 2L:(n-1L),
                ## off-diagonal
                2L:n)
            x = param$prec/(1-param$rho^2) *
                c(  ## diagonal
                    1L, 1L, rep(1+param$rho^2, n-2L),
                    ## off-diagonal
                    rep(-param$rho, n-1L))
            Q = sparseMatrix(i=i, j=j, x=x, giveCsparse=FALSE)
        }            
        return (Q)
    }

    mu = function()
    {
        return (numeric(0))
    }
        
    log.norm.const = function()
    {
        ## return the log(normalising constant) for the model
        param = interpret.theta()
        prec.innovation  = param$prec / (1.0 - param$rho^2)
        val = n * (- 0.5 * log(2*pi) + 0.5 * log(prec.innovation)) + 0.5 * log(1.0 - param$rho^2)
        return (val)
    }

    log.prior = function()
    {
        ## return the log-prior for the hyperparameters. the '+theta[1L]' is the log(Jacobian)
        ## for having a gamma prior on the precision and convert it into the prior for the
        ## log(precision).
        param = interpret.theta()
        val = (dgamma(param$prec, shape = 1, rate = 1, log=TRUE) + theta[1L] + 
                   dnorm(theta[2L], mean = 0, sd = 1, log=TRUE))
        return (val)
    }

    initial = function()
    {
        ## return initial values
        return (rep(1, 2))
    }

    quit = function()
    {
        return (invisible())
    }

    ## if theta is not required, it is not set. we set it here, for convenience.
    ## (see the graph() function)
    if (is.null(theta))
        theta = initial()
    
    val = do.call(match.arg(cmd), args = list())
    return (val)
}

`inla.rgeneric.iid.model` = function(
    cmd = c("graph", "Q", "mu", "initial", "log.norm.const", "log.prior", "quit"),
    theta = NULL)
{
    ## this is an example of the 'rgeneric' model. here we implement the iid model as described
    ## in inla.doc("iid"), without the scaling-option

    ## variables defined the in the define-call, are stored here
    ## (which is in the path)
    envir = parent.env(environment())
    
    interpret.theta = function()
    {
        return (list(prec = exp(theta[1L])))
    }

    graph = function()
    {
        G = Diagonal(n, x= rep(1, n))
        return (G)
    }

    Q = function()
    {
        prec = interpret.theta()$prec
        Q = Diagonal(n, x= rep(prec, n))
        return (Q)
    }

    mu = function()
    {
        return (numeric(0))
    }
    
    log.norm.const = function()
    {
        prec = interpret.theta()$prec
        val = sum(dnorm(rep(0, n), sd = 1/sqrt(prec), log=TRUE))
        return (val)
    }

    log.prior = function()
    {
        prec = interpret.theta()$prec
        val = dgamma(prec, shape = 1, rate = 1, log=TRUE) + theta[1L]
        return (val)
    }

    initial = function()
    {
        ntheta = 1
        return (rep(1, ntheta))
    }

    quit = function()
    {
        return (invisible())
    }

    ## if theta is not required, it is not set. we set it here, for convenience.
    if (is.null(theta))
        theta = initial()
    
    val = do.call(match.arg(cmd), args = list())
    return (val)
}

`inla.rgeneric.define` = function(model = NULL, debug = FALSE, ...)
{
    stopifnot(!missing(model))
    args = list(...)
    if (any(names(args) == "")) {
        stop("The '...' argument in 'inla.rgeneric.define()' needs *named* arguments.")
    }
    env = if (length(args) > 0) as.environment(args) else new.env()
    parent.env(env) = .GlobalEnv
    environment(model) = env

    rmodel = list(
        f = list(
            model = "rgeneric", 
            n = dim(model(cmd="graph", theta = NULL))[1], 
            rgeneric = list(
                definition = if (TRUE) {
                                 model
                             } else {
                                 ## did not see any speedup. maybe revisit this issue later...
                                 inla.require("compiler")
                                 compiler::cmpfun(model,
                                                  options = list(
                                                      optimize=3L,
                                                      suppressUndefined=TRUE))
                             }, 
                debug = debug
                )
            )
        )
    class(rmodel) = "inla.rgeneric"
    class(rmodel$f$rgeneric) = "inla.rgeneric"
    return (rmodel)
}

`inla.rgeneric.wrapper` = function(
    cmd = c("graph", "Q", "mu", "initial", "log.norm.const", "log.prior", "quit"),
    model, theta = NULL)
{
    debug.cat = function(...) {
        if (debug)
            cat("Rgeneric: ", ..., "\n", file = stderr())
    }

    model.orig = model
    if (is.character(model)) {
        model = get(model, envir = parent.frame())
    }
    stopifnot(inherits(model, "inla.rgeneric"))

    debug = ifelse(is.null(model$debug) || !model$debug, FALSE, TRUE)
    if (is.character(model.orig)) {
        debug.cat("Enter with cmd=", cmd, ", model=", model.orig, "theta=", theta)
    } else {
        debug.cat("Enter with cmd=", cmd, ", theta=", theta)
    }

    result = NULL
    cmd = match.arg(cmd)
    res = do.call(model$definition, args = list(cmd = cmd, theta = theta))
    time.ref = proc.time()[3]
    
    if (cmd %in% "Q") {
        Q = inla.as.sparse(res)
        debug.cat("dim(Q)", dim(Q))
        n = dim(Q)[1L]
        stopifnot(dim(Q)[1L] == dim(Q)[2L])
        idx = which(Q@i <= Q@j)
        len = length(Q@i[idx])
        result = c(n, len, Q@i[idx], Q@j[idx], Q@x[idx])
        Q = NULL
    } else if (cmd %in% "graph") {
        diag(res) = 1
        G = inla.as.sparse(res)
        stopifnot(dim(G)[1L] == dim(G)[2L])
        n = dim(G)[1L]
        idx = which(G@i <= G@j)
        len = length(G@i[idx])
        debug.cat("n", n, "len", len)
        result = c(n, len, G@i[idx], G@j[idx])
        G = NULL
    } else if (cmd %in% "mu") {
        mu = res
        debug.cat("length(mu)", length(mu))
        result = c(length(mu), mu)
    } else if (cmd %in% "initial") {
        init = res
        debug.cat("initial", init)
        result = c(length(init), init)
    } else if (cmd %in% "log.norm.const") {
        lnc = res
        debug.cat("log.norm.const", lnc)
        result = c(lnc)
    } else if (cmd %in% "log.prior") {
        lp = res
        debug.cat("log.prior", lp)
        result = c(lp)
    } else if (cmd %in% c("quit", "exit")) {
        ## nothing for the moment
        result = NULL
    } else {
        stop(paste("Unknown command", cmd))
    }
    res = NULL

    if (FALSE) {
        nm = "...cpu.time"
        envir = environment(model$definition)
        cpu.time = if (exists(nm, envir, envir)) get(nm, envir = envir) else list()
        if (is.null(cpu.time[[cmd]])) cpu.time[[cmd]] = list(total.time = 0, n.times = 0)
        cpu.time[[cmd]] =
            list(total.time = cpu.time[[cmd]]$total.time + proc.time()[3] - time.ref,
                 n.times = cpu.time[[cmd]]$n.times + 1)
        assign(nm, cpu.time, envir = envir)
    }
    
    return (as.numeric(result))
}

`inla.rgeneric.q` = function(rmodel,
                             cmd = c("graph", "Q", "mu", "initial", "log.norm.const",
                                     "log.prior", "quit"),
                             theta = NULL) 
{
    if (missing(cmd)) {
        stop("A value for argument 'cmd' is required.")
    }
    cmd = match.arg(cmd)
    rmodel.orig = rmodel
    if (is.character(rmodel)) {
        rmodel = get(rmodel, envir = parent.frame())
    }
    if (!inherits(rmodel, "inla.rgeneric")) {
        stop("Argument 'rmodel' is not of class 'inla.rgeneric' (usually the output of 'inla.rgeneric.define')")
    }
    func = rmodel$f$rgeneric$definition

    if (cmd %in% c("Q", "mu", "log.norm.const", "log.prior")) {
        ## for these we need values of 'theta': check that the length is correct
        initial = do.call(what = func, args = list(cmd = "initial", theta = NULL))
        if (length(initial) != length(theta)) {
            stop(paste0("Length of argument theta: ",  length(theta),
                        ", does not match the length of the initial values in 'rmodel': ", length(initial)))
        }
        ## just to make sure nothing else of length zero is passed
        if (length(initial) == 0)
            theta = NULL
    } else {
        theta = NULL
    }

    res = do.call(what = func, args = list(cmd = cmd, theta = theta))
    if (cmd %in% c("Q", "graph")) {
        ## since only the upper triangular matrix (diagonal included) is required return from
        ## 'do.call', then make sure its symmetric and that diag(Graph) = 1
        if (cmd %in% "Q") {
            Q = inla.as.sparse(res)
        } else {
            diag(res) = 1
            Q = inla.as.sparse(res, na.rm = TRUE, zeros.rm = TRUE)
            Q[Q != 0] = 1
        }
        n = dim(Q)[1]
        idx.eq = which(Q@i == Q@j)
        idx.gt = which(Q@i < Q@j)
        Q = sparseMatrix(i = c(Q@i[idx.eq], Q@i[idx.gt], Q@j[idx.gt]),
                         j = c(Q@j[idx.eq], Q@j[idx.gt], Q@i[idx.gt]),
                         x = c(Q@x[idx.eq], Q@x[idx.gt], Q@x[idx.gt]),
                         index1 = FALSE, 
                         dims = c(n, n),
                         giveCsparse = FALSE)
        return (Q)
    } else if (cmd %in% c("mu", "initial", "log.norm.const", "log.prior")) {
        return (c(as.numeric(res)))
    } else if (cmd %in% "quit") {
        return (NULL)
    } else {
        stop(paste("Unknown command", cmd))
    }
}
