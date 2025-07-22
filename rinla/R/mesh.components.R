#' Compute connected mesh subsets
#'
#' `r lifecycle::badge("deprecated")` Compute subsets of vertices and triangles
#' in an inla.mesh object that are
#' connected by edges. This function is deprecated from INLA `25.4.10` when
#' fmesher version `0.3.0.9005` or later is installed, which had
#' `fm_mesh_components()`, replaced from fmesher `0.4.0.9001` by
#' [fmesher::fm_components()].
#'
#'
#' @return A list with elements `vertex` and `triangle`, vectors of
#' integer labels for which connected component they belong, and `info`, a
#' `data.frame` with columns
#' \item{component}{Connected component integer label.}
#' \item{nV}{The number of vertices in the component.}
#' \item{nT}{The number of triangles in the component.}
#' \item{area}{The surface area associated with the component.
#' Component labels are not comparable across
#' different meshes, but some ordering stability is guaranteed by initiating
#' each component from the lowest numbered triangle whenever a new component is
#' initiated.}
#'
#' @param mesh An `fm_mesh_2d` object
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @seealso [fmesher::fm_mesh_2d()], [fmesher::fm_rcdt_2d()]
#' @examples
#'
#' # Construct two simple meshes:
#' library(fmesher)
#' loc <- matrix(c(0, 1, 0, 1), 2, 2)
#' mesh1 <- fm_mesh_2d(loc = loc, max.edge = 0.1)
#' bnd <- fm_nonconvex_hull(loc, 0.3)
#' mesh2 <- fm_mesh_2d(boundary = bnd, max.edge = 0.1)
#'
#' # Compute connectivity information:
#' conn1 <- fm_components(mesh1)
#' conn2 <- fm_components(mesh2)
#' # One component, simply connected mesh
#' conn1$info
#' # Two disconnected components
#' conn2$info
#' @export
inla.mesh.components <- function(mesh) {
    if (utils::packageVersion("fmesher") >= "0.4.0.9001") {
        fmesher_deprecate(
            "warn",
            2L,
            when = "25.04.10",
            what = "inla.mesh.components()",
            with = "fmesher::fm_components()",
            details =
                paste0("Use the `fmesher::fm_components()` method instead, ",
                       "from `fmesher` version `0.4.0.9001`.")
        )
        return(eval(parse(
            text = paste0("fmesher::fm_components(mesh)")
        )))
    }
    if (utils::packageVersion("fmesher") >= "0.3.0.9005") {
        fmesher_deprecate(
            "warn",
            2L,
            when = "25.04.10",
            what = "inla.mesh.components()",
            with = "fmesher::fm_components()",
            details =
                paste0("Use the `fmesher::fm_components()` method instead, ",
                       "from `fmesher` version `0.4.0.9001`, or from ",
                       "`0.3.0.9005` to `0.4.0`, ",
                       "`fmesher::fm_mesh_components()`.")
        )
        return(eval(parse(
            text = paste0("fmesher::fm_mesh_components(mesh)")
        )))
    }
    
    vertex <- integer(mesh$n)
    Nt <- nrow(mesh$graph$tv)
    triangle <- integer(Nt)
    ok <- !is.na(mesh$graph$tt)
    tt <- Matrix::sparseMatrix(
        i = c(
            which(ok[, 1]),
            which(ok[, 2]),
            which(ok[, 3])
        ),
        j = c(
            mesh$graph$tt[ok[, 1], 1],
            mesh$graph$tt[ok[, 2], 2],
            mesh$graph$tt[ok[, 3], 3]
        ),
        x = rep(1, sum(ok)),
        dims = c(Nt, Nt)
    )
    component <- 0
    while (any(triangle == 0)) {
        component <- component + 1
        tri <- integer(Nt)
        tri[min(which(triangle == 0))] <- 1
        tri_prev <- integer(Nt)
        while (any(tri_prev != tri)) {
            tri_prev <- tri
            tri <- tri | as.vector(tt %*% tri > 0)
        }
        triangle[tri > 0] <- component
        vtx <- sort(unique(as.vector(mesh$graph$tv[tri > 0, ])))
        if (any(vertex[vtx] > 0)) {
            warning(paste0(
                "Corner-only connected triangles detected.\n",
                "  Vertices = ", paste0(vtx[vertex[vtx] > 0], collapse = ", "), "\n",
                "  Components = ", paste0(vertex[vtx[vertex[vtx] > 0]], collapse = ", "), "\n",
                "  New component = ", component, "\n",
                "  Vertex component information will be inconsistent.", "\n",
                "  Triangle component information will ignore corner-only connections."
            ))
        }
        vertex[vtx] <- component
    }

    if (component == 0) {
        info <- data.frame(component = integer(0), nV = integer(0), nT = integer(0), area = numeric(0))
    } else {
        fem <- inla.mesh.fem(mesh, order = 1)
        info <- do.call(
            rbind,
            lapply(
                sort(unique(triangle)),
                function(x) {
                    data.frame(
                        component = x,
                        nV = sum(vertex == x),
                        nT = sum(triangle == x),
                        area = sum(fem$ta[triangle == x])
                    )
                }
            )
        )
    }

    list(
        vertex = vertex,
        triangle = triangle,
        info = info
    )
}
