#' cgeneric models
#' 
#' A framework for defining latent models in C
#' 
#' 
#' @aliases cgeneric cgeneric.define cgeneric.q inla.cgeneric.define
#' inla.cgeneric.q
#' @param model The name of the model function
#' @param cmodel The name of a cgeneric model-object (output from
#' \code{inla.cgeneric.define}
#' @param shlib Name of the compiled object-file with \code{model}
#' @param n The size of the model
#' @param debug Logical. Turn on/off debugging
#' @author Havard Rue \email{hrue@@r-inla.org}
#' 
#' @name cgeneric
#' @rdname cgeneric
#' @export
`inla.cgeneric.define` <- function(model = NULL, shlib = NULL, n = 0L, debug = FALSE, ...) {
    stopifnot(!missing(n))
    stopifnot(as.integer(n) > 0)
    stopifnot(!missing(model))
    stopifnot(!missing(shlib))
    shlib <- normalizePath(shlib)
    if (!file.exists(shlib)) {
        shlib <- paste0(getwd(), "/", shlib)
        stopifnot(file.exists(shlib))
    }

    args <- list(...)
    if (any(names(args) == "")) {
        stop("The '...' argument in 'inla.cgeneric.define()' needs *named* arguments.")
    }

    cmodel <- list(
        f = list(
            model = "cgeneric",
            n = as.integer(n), 
            cgeneric = list(
                model = model, 
                shlib = shlib, 
                n = as.integer(n), 
                debug = debug,
                data = inla.cgeneric.convert.arguments(n = as.integer(n), model = model, shlib = shlib, debug = as.integer(debug), ...)
            )
        )
    )
    class(cmodel) <- "inla.cgeneric"
    class(cmodel$f$cgeneric) <- "inla.cgeneric"
    return(cmodel)
}

`inla.cgeneric.convert.arguments` <-  function(n = 0L, model = NULL, shlib = NULL, debug = 0L, ...)
{
    a <- c(n = n, model = model, shlib = shlib, debug = debug, list(...))
    if (any(names(a) == "")) {
        stop("Unnamed arguments are not allowed")
    }

    res <- list(ints = list(), doubles = list(), characters = list(), matrices = list(), smatrices = list())
    for(idx in seq_along(a)) {
        xx <- a[idx][[1]]
        if (is.matrix(xx)) {
            ## store columnwise and not rowwise, so we need to revert the internal storage
            xx <- list(c(nrow(xx), ncol(xx), as.numeric(matrix(t(xx), nrow = ncol(xx), ncol = nrow(xx)))))
            names(xx) <- names(a)[idx]
            res$matrices <- c(res$matrices, xx)
        } else if (is(xx, "Matrix")) {
            if (!inherits(xx,"dgTMatrix")) {
                xx <- inla.as.sparse(xx)
            }
            xx <- list(c(nrow(xx), ncol(xx), length(xx@x), xx@i, xx@j, xx@x))
            names(xx) <- names(a)[idx]
            res$smatrices <- c(res$smatrices, xx)
        } else if (is.integer(xx)) {
            res$ints <- c(res$ints, a[idx])
        } else if (is.character(xx)) {
            res$characters <- c(res$characters, a[idx])
        } else if (is.double(xx)) {
            res$doubles <- c(res$doubles, a[idx])
        } else {
            stop(paste0("Unknown type: ", names(a)[idx], " ", idx))
        }
    }
    return (res)
}

#' @rdname cgeneric
#' @export
`inla.cgeneric.q` <- function(cmodel = NULL)
{
    stopifnot(!is.null(cmodel))
    stopifnot(inherits(cmodel$f$cgeneric, "inla.cgeneric"))

    ## need to turn off warnings for the numerics-conversions in 'split()'
    opt <- options()
    options(warn = -1)
    on.exit(options(opt))

    result <- list()
    ## add hidden options, which will make the inla-program print the contents of 'cmodel' and
    ## then exit
    tfile <- tempfile()
    cmodel$f$cgeneric$debug <- FALSE
    cmodel$f$cgeneric$.q <- TRUE
    cmodel$f$cgeneric$.q.file <- tfile

    try(inla(y ~ -1 + f(one, model = cmodel), data = data.frame(y = NA, one = 1),
             silent = 2L, verbose = FALSE), silent = TRUE)
    res <- readLines(tfile)
    unlink(tfile)
    
    ## cleanup the output
    for(i in 1:length(res)) {
        res[i] <- gsub("\t", "", res[i])
        res[i] <- gsub("[ ]+", " ", res[i])
    }

    split.char <- function(i) strsplit(res[i], "[ ]+")[[1]]
    split <- function(i) as.numeric(split.char(i))

    ## control...
    line <- 1
    stopifnot(split.char(line)[1] == "CGENERIC_BEGIN")
    line <- line + 2

    ntheta <- split(line)[3]
    line <- line+1

    theta <- numeric(ntheta)
    for(i in seq_len(ntheta)) {
        theta[i] <- split(line)[6]
        line <- line + 1
    }
    line <- line + 1
    result$theta <- theta

    nelm <- split(line)[6]
    line <- line + 1
    ii <- numeric(nelm)
    jj <- numeric(nelm)

    for(i in seq_len(nelm)) {
        ij <- split(line)[c(6, 9)]
        ii[i] <- ij[1]
        jj[i] <- ij[2]
        line <- line + 1
    }
    G <- sparseMatrix(i = ii +1, j = jj + 1)
    G <- G + t(G)
    diag(G) <- 1
    result$graph <- G

    line <- line+2
    stopifnot(split.char(line)[1] == "optimized")
    line <- line+1

    qq <- numeric(nelm)
    for(i in seq_len(nelm)) {
        qq[i] <- split(line)[6]
        line <- line + 1
    }
    Q <- sparseMatrix(i = ii + 1, j =  jj + 1, x = qq)
    dQ <- diag(Q)
    Q <- Q + t(Q)
    diag(Q) <- dQ
    result$Q <- Q

    line <- line + 1
    nelm <- split(line)[3]
    line <- line + 1
    mu <- numeric(nelm)
    for(i in seq_len(nelm)) {
        mu[i] <- split(line)[6]
        line <- line + 1
    }
    result$mu <- mu

    line <- line + 1
    result$log.prior <- split(line)[3]
    line <- line + 2
    result$log.norm.const <- split(line)[3]
    line <- line + 1
    ## control...
    stopifnot(split.char(line)[1] == "CGENERIC_END")

    return(result)
}
