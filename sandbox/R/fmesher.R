`inla.mesh` = function(s,manifold="R2")
{
##    mesh = list(s=s, tv=matrix(as.integer(c(1,2,3)),1,3))

    mesh.obj = inla.fmesher.mesh(s)
    mesh = (list(call.internal = mesh.obj$call,
                 prefix = mesh.obj$prefix,
                 s=mesh.obj$mesh$s, tv=mesh.obj$mesh$tv,
                 segm=mesh.obj$mesh$segm))

    mesh$call = match.call()
    class(mesh) <- "fmesher"

    return(mesh)
}

`summary.fmesher` = function(x, ...)
{
    ## provides a summary for a mesh object
    if (!inherits(x, "fmesher"))
        stop("'x' must inherit from class \"fmesher\"")

    ret = list()
    ret = c(ret, list(call=x$call))
    ret = c(ret, list(mesh=x))

    class(ret) <- "summary.fmesher"
    return (ret)
}

`print.summary.fmesher` = function(x, ...)
{
    if (!inherits(x, "summary.fmesher"))
        stop("'x' must inherit from class \"summary.fmesher\"")

    cat("\nCall:\n", deparse(x$call), "\n\n", sep = "")

    if(!is.null(x$mesh)) {
        cat("Vertices:\t", as.character(nrow(x$mesh$s)), "\n", sep="")
        cat("Triangles:\t", as.character(nrow(x$mesh$tv)), "\n", sep="")
        cat("\n")
    } else {
        cat("The mesh is empty\n\n")
    }

    if(!is.null(x$mesh$segm)) {
        cat("Boundary (groups):\t",
            as.character(max(0,nrow(x$mesh$segm$bnd$idx))),
            " (", as.character(length(unique(x$mesh$segm$bnd$grp))), ")",
            "\n", sep="")
        cat("Interior (groups):\t",
            as.character(max(0,nrow(x$mesh$segm$int$idx))),
            " (", as.character(length(unique(x$mesh$segm$int$grp))), ")",
            "\n", sep="")
        cat("\n")
    } else {
        cat("No segment information\n\n")
    }

}
