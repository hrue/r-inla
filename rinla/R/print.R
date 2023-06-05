#' Print an INLA fit
#' 
#' Print an INLA fit
#' 
#' @param x An inla-object (output from an \code{\link{inla}}-call).
#' @param digits Number of digits to print
#' @param ...  other arguments.
#' @return None
#' @author Havard Rue
#' @seealso \code{\link{inla}}
#' 
#' @rdname print
#' @export
`print.inla` <- function(x, digits = 3L, ...) {
    nsmall <- 2L
    form <- strwrap(inla.formula2character(x$call))
    cat("Call:\n")
    for (i in seq_along(form)) {
        cat("  ", form[i], "\n")
    }

    if (inla.is.element("cpu.used", x)) {
        cat("Time used:\n  ")
        cat(
            sep = "",
            names(x$cpu.used)[1], " = ",
            format(x$cpu.used[1], digits = digits), ", ",
            names(x$cpu.used)[2], " = ",
            format(x$cpu.used[2], digits = digits), ", ",
            names(x$cpu.used)[3], " = ",
            format(x$cpu.used[3], digits = digits), ", ",
            names(x$cpu.used)[4], " = ",
            format(x$cpu.used[4], digits = digits), "\n"
        )
    }

    return(invisible())
}
