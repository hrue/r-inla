#' Do a dryrun to extract some internal information upfront
#' 
#' Do a dryrun to get information about the internal storage and the list (and
#' ordering) of the hyperparameters
#' 
#' 
#' @aliases inla.dryrun dryrun
#' @param ... Same arguments as `inla`()
#' @return A list of start-index and length for each latent component and a
#' list of the hyperparameters in the model
#' @author Havard Rue \email{hrue@@r-inla.org}
#' 
#' @name dryrun
#' @rdname dryrun
#' @export

`inla.dryrun` <- function(...) {

    split.string <- function(x) {
        return (unlist(strsplit(x, "[ \t]+")))
    }
    to.ints <- function(x) {
        return (as.integer(split.string(x)))
    }

    args <- list(...)
    args$verbose <- FALSE
    args$debug <- FALSE
    args$safe <- FALSE
    args$inla.call <- inla.call.builtin()
    args$inla.arg <- "-t1 -s -m dryrun"
    r <- unclass(do.call("inla", args = args))

    off <- 1
    x <- to.ints(r[off])
    off <- off+1
    nmc <- x[1]
    totlen <- x[2]

    tag <- c()
    start <- c()
    len <- c()
    for(i in seq_len(nmc)) {
        x <- split.string(r[off])
        tag <- c(tag, x[1])
        start <- c(start, as.integer(x[2]) + 1)
        len <- c(len, as.integer(x[3]))
        off <- off + 1
    }

    if (tag[1] != "APredictor") {
        tag <- c("APredictor", tag)
        start <- c(1, start)
        len <- c(0, len)
    }

    ll <- 1:2
    pred.tag <- tag[ll]
    pred.start <- start[ll]
    pred.len <- len[ll]

    nhyper <- to.ints(r[off])
    off <- off + 1
    hyperpar <- c()
    for (i in seq_len(nhyper)) {
        hyperpar <- c(hyperpar, r[off])
        off <- off + 1
    }
    
    res <- list()
    res$classic = list(latent = list(tag = tag, start = start, length = len),
                       hyperpar = list(nhyper = nhyper, hyperpar = hyperpar))

    for(nm in c("APredictor", "Predictor")) {
        if (tag[1] == nm) {
            tag <- tag[-1]
            start <- start[-1] - start[2] + 1
            len <- len[-1]
        }
    }

    res$compact <- list(predictor = list(tag = pred.tag, start = pred.start, length = pred.len), 
                        latent = list(tag = tag, start = start, length = len),
                        hyperpar = list(nhyper = nhyper, hyperpar = hyperpar))

    return(res)
}
