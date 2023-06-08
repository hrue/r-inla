#' Check which mesh triangles are inside a polygon
#'
#' Wrapper for the [sp::over()] method to find triangle centroids
#' or vertices inside `sp` polygon objects.
#' Deprecated since 23.06.06 in favour of `inlabru::fm_contains()` when `inlabru`
#' version `>= 2.7.0.9011` is installed.
#'
#' @param x geometry (typically a [sp::SpatialPolygons()] object) for the queries
#' @param y an [inla.mesh()] object
#' @param type the query type; either `'centroid'` (default, for triangle centroids),
#' or `'vertex'` (for mesh vertices)
#' @param ignore.CRS logical; whether to ignore the coordinate system information in `x` and `y`
#' (default `FALSE`)
#'
#' @return A vector of triangle indices (when `type` is `'centroid'`) or
#' vertex indices (when `type` is `'vertex'`)
#'
#' @author Haakon Bakka, \email{bakka@@r-inla.org}, and Finn Lindgren \email{finn.lindgren@@gmail.com}
#'
#' @examples
#' # Create a polygon and a mesh
#' obj <- sp::SpatialPolygons(list(Polygons(list(Polygon(rbind(
#'     c(0, 0),
#'     c(50, 0),
#'     c(50, 50),
#'     c(0, 50)
#' ))),
#' ID = 1
#' )),
#' proj4string = inla.CRS("longlat_globe")
#' )
#' mesh <- inla.mesh.create(globe = 2, crs = inla.CRS("sphere"))
#'
#' ## 3 vertices found in the polygon
#' inla.over_sp_mesh(obj, mesh, type = "vertex")
#'
#' ## 3 triangles found in the polygon
#' inla.over_sp_mesh(obj, mesh)
#'
#' ## Multiple transformations can lead to slightly different results due to edge cases
#' ## 4 triangles found in the polygon
#' inla.over_sp_mesh(
#'   obj,
#'   inla.spTransform(mesh, CRSobj = inla.CRS("mollweide_norm")),
#'   ignore.CRS = FALSE)
#'
#' ## Ignoring mismatching coordinate systems is rarely useful
#' ## 20 triangles "found in" the polygon
#' inla.over_sp_mesh(obj,
#'   inla.spTransform(mesh, CRSobj = inla.CRS("mollweide_norm")),
#'   ignore.CRS = TRUE)
#' @export
inla.over_sp_mesh <- function(x, y, type = c("centroid", "vertex"), ignore.CRS = FALSE) {
    if ((getNamespaceVersion("inlabru") >= "2.7.0.9011")) {
        if (inla.getOption("fmesher.evolution") >= 2L) {
            lifecycle::deprecate_soft("23.06.06",
                                      what = "inla.over_sp_mesh()",
                                      with = "inlabru::fm_contains()",
                                      details = c("fm_contains() is available from inlabru version 2.7.0.9011",
                                                  "For equivalent output, use 'unlist(fm_contains(...))'"))
            if (!missing(ignore.CRS) && isTRUE(ignore.CRS)) {
                lifecycle::deprecate_soft("23.06.06",
                                          what = "inla.over_sp_mesh(ignore.CRS)",
                                          with = "ignore.CRS will be treated as FALSE")
            }
        }
        return(unlist(inlabru::fm_contains(x = x, y = y, type = type)))
    }

    if (!inherits(y, "inla.mesh")) {
        stop(paste0(
            "'y' must be an 'inla.mesh' object, not '",
            paste0(class(y), collapse = ", "),
            "'."
        ))
    }
    type <- match.arg(type)
    if (identical(type, "centroid")) {
        ## Extract triangle centroids
        points <- (y$loc[y$graph$tv[, 1], , drop = FALSE] +
            y$loc[y$graph$tv[, 2], , drop = FALSE] +
            y$loc[y$graph$tv[, 3], , drop = FALSE]) / 3
        if (identical(y$manifold, "S2")) {
            points <- points / rowSums(points^2)^0.5
        }
    } else if (identical(type, "vertex")) {
        ## Extract vertices
        points <- y$loc
    }
    ## Convert to SpatialPoints
    if (ignore.CRS) {
        if (identical(y$manifold, "S2")) {
            stop("sp::over cannot operate on geocentric coordinates.")
        }
        ## Ignore any actual CRS information, and copy the one from x
        points <- sp::SpatialPoints(points, proj4string = inla.sp_get_crs(x))
    } else {
        ## Extract coordinate system information
        if (identical(y$manifold, "S2")) {
            crs <- inla.CRS("sphere")
        } else {
            crs <- y$crs
        }
        crs_x <- inla.sp_get_crs(x)
        ## Create SpatialPoints object and transform the coordinates.
        points <- sp::SpatialPoints(points, proj4string = crs)
        if (!is.null(inla.crs_get_wkt(crs)) &&
            !is.null(inla.crs_get_wkt(crs_x))) {
            ## Convert to the target object CRS
            points <- inla.spTransform(points, CRSobj = crs_x)
        }
    }
    ## Find indices:
    ids <- unlist(sp::over(x, points, returnList = TRUE))

    ids
}
