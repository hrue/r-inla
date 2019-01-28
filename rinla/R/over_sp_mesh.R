#' Check which mesh triangles are inside a polygon
#'
#' Wrapper for the \code{\link[sp]{over}} method to find triangle centroids
#' or vertices inside \code{sp} polygon objects
#'
#' @param x geometry (typically a \code{\link[sp]{SpatialPolygons}} object) for the queries
#' @param y an \code{\link{inla.mesh}} object
#' @param type the query type; either \code{'centroid'} (default, for triangle centroids),
#'   or \code{'vertex'} (for mesh vertices)
#' @param ignore.CRS logical; whether to ignore the coordinate system information in \code{x} and \code{y}
#'   (default \code{FALSE})
#'
#' @return A vector of triangle indices (when \code{type} is \code{'centroid'}) or
#'   vertex indices (when \code{type} is \code{'vertex'})
#'   
#' @author Haakon Bakka, \email{bakka@@r-inla.org}, and Finn Lindgren \email{finn.lindgren@@gmail.com}
#'
#' @examples
#' # Create a polygon and a mesh
#' obj <- sp::SpatialPolygons(list(Polygons(list(Polygon(rbind(c(0,0),
#'                                                             c(50,0),
#'                                                             c(50,50),
#'                                                             c(0,50)))),
#'                                          ID=1)),
#'                            proj4string = inla.CRS("longlat"))
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
#' inla.over_sp_mesh(obj, inla.spTransform(mesh, CRSobj=inla.CRS("mollweide")), ignore.CRS = FALSE)
#'
#' ## Ignoring mismatching coordinate systems is rarely useful
#' ## 20 triangles "found in" the polygon
#' inla.over_sp_mesh(obj, inla.spTransform(mesh, CRSobj=inla.CRS("mollweide")), ignore.CRS = TRUE)
#'
#' @export

inla.over_sp_mesh <- function(x, y, type = c("centroid", "vertex"), ignore.CRS=FALSE) {
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
      points <- points / rowSums(points ^ 2) ^ 0.5
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
    points <- sp::SpatialPoints(points, proj4string = CRS(proj4string(x)))
  } else {
    ## Extract coordinate system information
    if (identical(y$manifold, "S2")) {
      crs <- inla.CRSargs(inla.CRS("sphere"))
    } else {
      crs <- inla.CRSargs(y$crs)
    }
    ## Create SpatialPoints object and transform the coordinates.
    points <- sp::SpatialPoints(points, proj4string = inla.CRS(crs))
    if (!is.na(crs) & !is.na(proj4string(x))) {
      ## Convert to the target object CRS
      points <- inla.spTransform(points, CRSobj = inla.CRS(proj4string(x)))
    }
  }
  ## Find indices:
  ids <- unlist(sp::over(x, points, returnList = TRUE))

  ids
}
