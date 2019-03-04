## Export: print!inla

##!\name{print.inla}
##!\alias{print.inla}
##!\title{Print a INLA fit }
##!\description{
##!Print a INLA fit
##!}
##!\usage{
##!\method{print}{inla}(x, digits = 3L, ...)
##!}
##!\arguments{
##!  \item{x}{An inla-object (output from an \code{\link{inla}}-call).}
##!  \item{digits}{Number of digits to print}
##!  \item{...}{ other arguments.}
##!}
##!\details{
##! None
##!}
##!\value{
##! None
##!}
##!\author{Havard Rue}
##!
##!\seealso{ \code{\link{inla}} }
##!\examples{
##!## None
##!}

`print.inla` = function(x, digits = 3L, ...)
{
    nsmall = 2L
    form = strwrap(inla.formula2character(x$call))
    cat("Call:\n")
    for (i in seq_along(form)) {
        cat("  ", form[i], "\n")
    }

    if (inla.is.element("cpu.used", x)) {
        cat("Time used:\n  ")
        cat(sep = "", 
            names(x$cpu.used)[1], " = ",
            format(x$cpu.used[1], digits = digits), ", ", 
            names(x$cpu.used)[2], " = ",
            format(x$cpu.used[2], digits = digits), ", ", 
            names(x$cpu.used)[3], " = ",
            format(x$cpu.used[3], digits = digits), ", ",
            names(x$cpu.used)[4], " = ",
            format(x$cpu.used[4], digits = digits), "\n")
    }

    return(invisible())
}
