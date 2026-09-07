#' Check which mesh triangles are inside a polygon
#'
#' Wrapper for the `sp::over()` method to find triangle centroids
#' or vertices inside `sp` polygon objects.
#' `r lifecycle::badge("deprecated")` since 23.06.06 in favour of
#' `inlabru::fm_contains()` when `inlabru`
#' version `>= 2.7.0.9011` is installed, and since 23.08.02 in favour of
#' `fmesher::fm_contains()` when `fmesher`.
#'
#' @param x geometry (typically a `sp::SpatialPolygons()` object) for the queries
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
#' if (require("sp", quietly = TRUE)) {
#'     # Create a polygon and a mesh
#'     obj <- sp::SpatialPolygons(
#'         list(sp::Polygons(
#'             list(sp::Polygon(rbind(
#'                 c(0, 0),
#'                 c(50, 0),
#'                 c(50, 50),
#'                 c(0, 50)
#'             ))),
#'             ID = 1
#'         )),
#'         proj4string = fmesher::fm_CRS("longlat_globe")
#'     )
#'     mesh <- inla.mesh.create(globe = 2, crs = fmesher::fm_CRS("sphere"))
#'
#'     ## 3 vertices found in the polygon
#'     inla.over_sp_mesh(obj, mesh, type = "vertex")
#'
#'     ## 3 triangles found in the polygon
#'     inla.over_sp_mesh(obj, mesh)
#'
#'     ## Multiple transformations can lead to slightly different results due to edge cases
#'     ## 4 triangles found in the polygon
#'     inla.over_sp_mesh(
#'         obj,
#'         fmesher::fm_transform(mesh, crs = fmesher::fm_crs("mollweide_norm"))
#'     )
#' }
#' @export
inla.over_sp_mesh <- function(x, y, type = c("centroid", "vertex"), ignore.CRS = FALSE) {
    fmesher_deprecate("soft", 1L,
        "23.06.06",
        what = "inla.over_sp_mesh()",
        with = "fmesher::fm_contains()",
        details = c(
            "fm_contains() is available from fmesher version 0.0.9",
            "For equivalent output, use 'unlist(fm_contains(...))'"
        )
    )
    if (!missing(ignore.CRS) && isTRUE(ignore.CRS)) {
        fmesher_deprecate(
            "soft", 1L,
            "23.06.06",
            what = "inla.over_sp_mesh(ignore.CRS)",
            with = "ignore.CRS will be treated as FALSE"
        )
    }

    type <- match.arg(type)

    return(unlist(fmesher::fm_contains(x = x, y = y, type = type)))
}
