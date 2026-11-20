### #' Compute log-determinants for a GMRF 
### #' 
### #' Computes the log-determinant of (Q + d * I) with optional constraints, for a possible sequence of d's
### #' 
### #' @aliases inla.qsample qsample
### #' @param Q The precision matrix or a filename containing it.
### #' @param d A possible vector of values (default d=0)
### #' @param constr Optional linear constraints; see `?INLA::f` and argument `extraconstr`
### #' @param num.threads Number of threads to use. If in format "A:B:C", then "A" is used.
### #' @return The log-determinant of (Q + d * I) 
### #' 
### #' @author Havard Rue \email{hrue@@r-inla.org}
### #' @examples
### #' 
### #' @name qlogdet
### #' @rdname qlogdet
### #' @export

`inla.qlogdet` <- function(Q, d = 0, constr, num.threads = NULL)
{
    t.dir <- inla.tempdir()
    smtp <- match.arg(inla.getOption("smtp"), c("stiles", "taucs", "band", "default"))
    stopifnot(!missing(Q))

    if (is.null(num.threads)) {
        num.threads <- inla.getOption("num.threads")
    }
    cores <- numeric(2)
    cores[1] <- as.integer(max(1, strsplit(inla.parse.num.threads(num.threads), ":")[[1]][1]))
    cores[2] <- as.integer(max(1, strsplit(inla.parse.num.threads(num.threads), ":")[[1]][2]))
    nt <- paste0("1:", prod(cores))

    Q <- inla.sparse.check(Q)
    if (is(Q, "dgTMatrix")) {
        Q.file <- inla.write.fmesher.file(Q, filename = inla.tempfile(tmpdir = t.dir))
    } else if (is.character(Q)) {
        Q.file <- Q
    } else {
        stop("This should not happen.")
    }

    out.file <- inla.tempfile(tmpdir = t.dir)
    d.file <- inla.tempfile(tmpdir = t.dir)
    inla.write.fmesher.file(matrix(d, length(d), 1), filename = d.file)

    constr.file <- inla.tempfile(tmpdir = t.dir)
    if (!missing(constr) && !is.null(constr)) {
        stopifnot(is.list(constr))
        A <- as.matrix(constr$A)
        e <- as.numeric(constr$e)
        stopifnot(ncol(A) == ncol(Q))
        stopifnot(nrow(A) == length(e))
        xx <- matrix(c(nrow(A), c(A), c(e)), ncol = 1)
        inla.write.fmesher.file(xx, filename = constr.file)
    }

    inla.set.environment()
    env <- inla.run.environment.set()
    if (inla.os("linux") || inla.os("mac") || inla.os("mac.arm64")) {
        s <- system(paste(
            shQuote(inla.call.no.remote()), "-s -m qlogdet",
            paste0("-t", nt), "-v", "-S", smtp, 
            Q.file, d.file, constr.file, out.file
        ), intern = TRUE)
    } else if (inla.os("windows")) {
        s <- system(paste(
            shQuote(inla.call.no.remote()), "-E -s -m qlogdet",
            paste0("-t", nt), "-S", smtp, 
            Q.file, d.file, constr.file, out.file
        ), intern = TRUE)
    } else {
        stop("\n\tNot supported architecture.")
    }

    inla.run.environment.unset(env)
    out <- as.numeric(inla.read.fmesher.file(out.file))
    unlink(t.dir, recursive = TRUE)
    return(out)
}
