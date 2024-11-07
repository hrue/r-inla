## Internal: inla.mesh.filter.locations
## Internal: inla.mesh.parse.segm.input inla.mesh.extract.segments


#' @importFrom fmesher fm_transform


#' @title Constraint segments for inla.mesh
#'
#' @description
#' `r lifecycle::badge("deprecated")` Use [fmesher::fm_segm()] instead.
#'
#' Constructs `inla.mesh.segment` objects that can be used to specify
#' boundary and interior constraint edges in calls to [inla.mesh()].
#'
#' @param loc Matrix of point locations, or `SpatialPoints`, or `sf`/`sfc` point
#' object.
#' @param idx Segment index sequence vector or index pair matrix.  The indices
#' refer to the rows of `loc`.  If `loc==NULL`, the indices will be
#' interpreted as indices into the point specification supplied to
#' [inla.mesh.create()].  If `is.bnd==TRUE`, defaults to linking
#' all the points in `loc`, as `c(1:nrow(loc),1L)`, otherwise
#' `1:nrow(loc)`.
#' @param grp Vector of group labels for each segment.  Set to `NULL` to
#' let the labels be chosen automatically in a call to
#' [inla.mesh.create()].
#' @param is.bnd `TRUE` if the segments are boundary segments, otherwise
#' `FALSE`.
#' @param grp.default When joining segments, use this group label for segments
#' that have `grp=NULL`.
#' @param x,y,z,nlevels,levels Parameters specifying a set of surface contours,
#' with syntax described in [contour()].
#' @param groups Vector of group ID:s, one for each contour level.
#' @param positive `TRUE` if the contours should encircle positive level
#' excursions in a counter clockwise direction.
#' @param eps Tolerance for [inla.simplify.curve()].
#' @param crs An optional `CRS` or `inla.CRS` object
#' @param ...  Additional parameters.  When joining segments, a list of
#' `inla.mesh.segment` objects.
#' @return An `inla.mesh.segment` object.
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @seealso [inla.mesh.create()], [inla.mesh.2d()]
#' @examples
#'
#' ## Create a square boundary and a diagonal interior segment
#' loc.bnd <- matrix(c(0, 0, 1, 0, 1, 1, 0, 1), 4, 2, byrow = TRUE)
#' loc.int <- matrix(c(0.9, 0.1, 0.1, 0.6), 2, 2, byrow = TRUE)
#' segm.bnd <- inla.mesh.segment(loc.bnd)
#' segm.int <- inla.mesh.segment(loc.int, is.bnd = FALSE)
#'
#' ## Points to be meshed
#' loc <- matrix(runif(10 * 2), 10, 2) * 0.9 + 0.05
#' mesh <- inla.mesh.create(loc,
#'   boundary = segm.bnd,
#'   interior = segm.int,
#'   refine = list()
#' )
#' plot(mesh)
#'
#' mesh <- inla.mesh.create(loc, interior = fm_segm_join(segm.bnd, segm.int))
#' plot(mesh)
#'
#' @export inla.mesh.segment
inla.mesh.segment <- function(...) {
  fmesher_deprecate("soft",
    2L,
    "23.08.18",
    "inla.mesh.segment()",
    "fmesher::fm_segm()",
    id = "inla.mesh.segment"
  )

  return(fmesher::fm_segm(...))
}



#' @title Draw `inla.mesh.segment` objects.
#'
#' @description
#' `r lifecycle::badge("deprecated")` Use [fmesher::lines.fm_segm()] or
#' [fmesher::lines_rgl()] instead.
#'
#' Draws a [inla.mesh.segment()] object with generic or `rgl`
#' graphics.
#'
#' @param x An [inla.mesh.segment()] object.
#' @param loc Point locations to be used if `x$loc` is `NULL`.
#' @param col Segment color specification.
#' @param colors Colors to cycle through if `col` is `NULL`.
#' @param add If `TRUE`, add to the current plot, otherwise start a new
#' plot.
#' @param xlim X axis limits for a new plot.
#' @param ylim Y axis limits for a new plot.
#' @param rgl If `TRUE`, use `rgl` for plotting.
#' @param \dots Additional parameters, passed on to graphics methods.
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @seealso [inla.mesh.segment()]
#' @export
lines.inla.mesh.segment <- function(x, loc = NULL, col = NULL,
                                    colors = c("black", "blue", "red", "green"),
                                    add = TRUE, xlim = NULL, ylim = NULL,
                                    rgl = FALSE, ...) {
  if (rgl) {
    fmesher_deprecate(
      "soft",
      2L,
      "23.08.18",
      "lines.inla.mesh.segment()",
      "fmesher::lines_rgl()"
    )
    return(fmesher::lines_rgl(
      x = fmesher::fm_as_fm(x),
      add = add,
      loc = loc,
      col = col,
      colors = colors,
      ...
    ))
  } else {
    fmesher_deprecate(
      "soft",
      2L,
      "23.08.18",
      "lines.inla.mesh.segment()",
      "fmesher::lines()"
    )
    return(lines(
      x = fmesher::fm_as_fm(x),
      add = add,
      loc = loc,
      col = col,
      colors = colors,
      xlim = xlim,
      ylim = ylim,
      ...
    ))
  }
}

#' @title Generate text RGB color specifications.
#'
#' @description
#' `r lifecycle::badge("deprecated")` Use [fmesher::fm_generate_colors()]
#' instead.
#'
#' Generates a tex RGB color specification matrix based on a color palette.
#'
#'
#' @param color `character`, `matrix` or `vector`
#' @param color.axis The min/max limit values for the color mapping.
#' @param color.n The number of colors to use in the color palette.
#' @param color.palette A color palette function.
#' @param color.truncate If `TRUE`, truncate the colors at the color axis
#' limits.
#' @param alpha Transparency/opaqueness values.
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @export inla.generate.colors
inla.generate.colors <- function(color,
                                 color.axis = NULL,
                                 color.n = 512,
                                 color.palette = cm.colors,
                                 color.truncate = FALSE,
                                 alpha = NULL) {
  fmesher_deprecate(
    "soft",
    2L,
    "23.08.18",
    "inla.generate.colors()",
    "fmesher::fm_generate_colors()"
  )
  return(fmesher::fm_generate_colors(
    color = color,
    color.axis = color.axis,
    color.n = color.n,
    color.palette = color.palette,
    color.truncate = color.truncate,
    alpha = alpha
  ))
}


#' @title Low level triangulation mesh plotting
#'
#' @description
#' `r lifecycle::badge("deprecated")` Use [fmesher::plot_rgl()] instead.
#'
#' Plots a triangulation mesh using `rgl`.
#'
#'
#' @param x A 3-column triangle-to-vertex index map matrix.
#' @param S A 3-column vertex coordinate matrix.
#' @param color Color specification.  A single named color, a vector of scalar
#' values, or a matrix of RGB values.
#' @param color.axis The min/max limit values for the color mapping.
#' @param color.n The number of colors to use in the color palette.
#' @param color.palette A color palette function.
#' @param color.truncate If `TRUE`, truncate the colors at the color axis
#' limits.
#' @param alpha Transparency/opaqueness values. See `rgl.material`.
#' @param lwd Line width for edges. See `rgl.material`.
#' @param specular Specular color. See `rgl.material`.
#' @param draw.vertices If `TRUE`, draw triangle vertices.
#' @param draw.edges If `TRUE`, draw triangle edges.
#' @param edge.color Edge color specification.
#' @param \dots Additional parameters passed to and from other methods.
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @seealso [plot.inla.mesh()]
#' @method plot inla.trimesh
#' @export
plot.inla.trimesh <- function(x, S, color = NULL, color.axis = NULL,
                              color.n = 512, color.palette = cm.colors,
                              color.truncate = FALSE, alpha = NULL,
                              lwd = 1, specular = "black",
                              draw.vertices = TRUE, draw.edges = TRUE,
                              edge.color = rgb(0.3, 0.3, 0.3), ...) {
  fmesher_deprecate(
    "soft",
    2L,
    "23.08.18",
    "plot.inla.trimesh()",
    "fmesher::plot_rgl()"
  )
  return(fmesher::plot_rgl(
    x = fmesher::fm_as_fm(x),
    col = col,
    color.axis = color.axis,
    color.n = color.n,
    color.palette = color.palette,
    alpha = alpha,
    lwd = lwd,
    specular = specular,
    draw.vertices = draw.vertices,
    draw.edges = draw.edges,
    draw.faces = TRUE,
    edge.color = edge.color,
    ...
  ))
}

## library(geometry)
## S = cbind(x=rnorm(30), y=rnorm(30), z=0)
## TV = delaunayn(S[, 1:2]) # NOTE: inconsistent triangle orders, only for test.
## trimesh(TV, S)
##
## colors = rgb(runif(30), runif(30), runif(30))
## rgl.viewpoint(0, 0, fov=20)
## plot.inla.trimesh(TV, S, colors)

## Ecol = col2rgb(color)/256
## Ecol = Ecol*0.5+(1-0.5)*0 # Rescale towards black
## Ecol = 1-Ecol # Invert
## Ecol = Ecol[, c(2, 3, 1)] # Permute
## Ecol = rgb(Ecol[1,], Ecol[2,], Ecol[3,], maxColorValue = 1)
## Ecol = Ecol[tETV]


#' @title Draw a triangulation mesh object
#'
#' @description
#' `r lifecycle::badge("deprecated")` Use [fmesher::plot.fm_mesh_2d()] or
#' [fmesher::plot_rgl()] instead.
#'
#' Plots an [inla.mesh()] object using either standard graphics or
#' with `rgl`.
#'
#'
#' @param x An [inla.mesh()] object.
#' @param col Color specification.  A single named color, a vector of scalar
#' values, or a matrix of RGB values.  Requires `rgl=TRUE`.
#' @param t.sub Optional triangle index subset to be drawn.
#' @param add If `TRUE`, adds to the current plot instead of starting a
#' new one.
#' @param lwd Line width for triangle edges.
#' @param xlim X-axis limits.
#' @param ylim Y-axis limits.
#' @param main The main plot title.  If not specified, a default title is
#' generated based on the mesh type.
#' @param rgl When `TRUE`, generates an `rgl` plot instead of a
#' generic graphics plot.  Allows 3D plotting and color surface plotting.
#' @param size Size of vertex points in `rgl` plotting. See
#' `rgl.material`.
#' @param draw.vertices If `TRUE`, draw triengle vertices.
#' @param vertex.color Color specification for all vertices.
#' @param draw.edges If `TRUE`, draw triangle edges.
#' @param edge.color Color specification for all edges.
#' @param draw.segments If `TRUE`, draw boundary and interior constraint
#' edges more prominently.
#' @param \dots Further graphics parameters, interpreted by the respective
#' plotting systems.
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @seealso [plot.inla.trimesh()]
#' @examples
#'
#' mesh <- inla.mesh.create(globe = 10)
#' plot(mesh)
#' \donttest{
#' if (require(rgl)) {
#'   plot(mesh, rgl = TRUE, col = mesh$loc[, 1])
#' }
#' }
#'
#' @method plot inla.mesh
#' @export
plot.inla.mesh <- function(x,
                           col = "white",
                           t.sub = 1:nrow(mesh$graph$tv),
                           add = FALSE,
                           lwd = 1,
                           xlim = range(mesh$loc[, 1]),
                           ylim = range(mesh$loc[, 2]),
                           main = NULL,
                           rgl = FALSE,
                           size = 2,
                           draw.vertices = FALSE,
                           vertex.color = "black",
                           draw.edges = TRUE,
                           edge.color = rgb(0.3, 0.3, 0.3),
                           draw.segments = draw.edges,
                           ...) {
  if (rgl) {
    fmesher_deprecate(
      "soft",
      2L,
      "23.08.18",
      "plot.inla.mesh()",
      "fmesher::plot_rgl()"
    )
    return(fmesher::plot_rgl(
      x = fmesher::fm_as_fm(x),
      col = col,
      lwd = lwd,
      size = size,
      draw.vertices = TRUE,
      draw.edges = TRUE,
      draw.faces = TRUE,
      draw.segments = draw.segments,
      edge.color = edge.color,
      ...
    ))
  } else {
    fmesher_deprecate(
      "soft",
      2L,
      "23.08.18",
      "plot.inla.mesh()",
      "fmesher::plot()"
    )
    plot(
      x = fmesher::fm_as_fm(x),
      col = col,
      t.sub = t.sub,
      add = add,
      pwd = lwd,
      xlim = xlim,
      ylim = ylim,
      size = size,
      draw.vertices = draw.vertices,
      vertex.color = vertex.color,
      draw.edges = draw.edges,
      edge.color = edge.color,
      draw.segments = draw.segments,
      ...
    )
    if (!add) {
      if (missing(main)) {
        if (mesh$meta$is.refined) {
          title("Constrained refined Delaunay triangulation")
        } else {
          title("Constrained Delaunay triangulation")
        }
      } else if (!is.null(main)) {
        title(main)
      }
    }
    return(invisible())
  }
}


#' @export
#' @describeIn inla.mesh.map
#' `r lifecycle::badge("deprecated")` Use [fmesher::fm_mesh_2d_map_lim()]
#' instead.
#'
#' Projection extent limit calculations
inla.mesh.map.lim <- function(loc = NULL,
                              projection =
                                c("default", "longlat", "longsinlat", "mollweide")) {
  fmesher_deprecate(
    "soft",
    2L,
    "23.08.18",
    "inla.mesh.map.lim()",
    "fmesher::fm_mesh_2d_map_lim()"
  )
  return(fmesher::fm_mesh_2d_map_lim(loc = loc, projection = projection))
}


#' @title Coordinate mappings for `inla.mesh` projections.
#'
#' @description
#' `r lifecycle::badge("deprecated")` Use [fmesher::fm_mesh_2d_map()] instead.
#'
#' Calculates coordinate mappings for `inla.mesh` projections.
#'
#'
#' @aliases inla.mesh.map inla.mesh.map.lim
#' @param loc Coordinates to be mapped.
#' @param projection The projection type.
#' @param inverse If `TRUE`, `loc` are map coordinates and
#' coordinates in the mesh domain are calculated.  If `FALSE`, `loc`
#' are coordinates in the mesh domain and the forward map projection is
#' calculated.
#' @return For `inla.mesh.map.lim`, a list: \item{xlim }{X axis limits in
#' the map domain} \item{ylim }{Y axis limits in the map domain} No attempt is
#' made to find minimal limits for partial spherical domains.
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @seealso [inla.mesh.project()]
#' @export inla.mesh.map
inla.mesh.map <- function(loc,
                          projection =
                            c("default", "longlat", "longsinlat", "mollweide"),
                          inverse = TRUE) {
  fmesher_deprecate(
    "soft",
    2L,
    "23.08.18",
    "inla.mesh.map()",
    "fmesher::fm_mesh_2d_map()"
  )
  return(fmesher::fm_mesh_2d_map(
    loc = loc,
    projection = projection,
    inverse = inverse
  ))
}



#' @title Lattice grids for inla.mesh
#'
#' @description
#' `r lifecycle::badge("deprecated")` Use [fmesher::fm_lattice_2d()] instead.
#'
#' Construct a lattice grid for [inla.mesh()]
#'
#'
#' @param x vector or grid matrix of x-values
#' @param y vector of grid matrix of y-values
#' @param z if x is a matrix, a grid matrix of z-values
#' @param dims the size of the grid, length 2 vector
#' @param units One of `c("default", "longlat", "longsinlat")`.
#' @param crs An optional `CRS` or `inla.CRS` object
#' @return An `inla.mesh.lattice` object.
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @seealso [inla.mesh()]
#' @examples
#'
#' lattice <- inla.mesh.lattice(seq(0, 1, length.out = 17), seq(0, 1, length.out = 10))
#'
#' ## Use the lattice "as-is", without refinement:
#' mesh <- inla.mesh.create(lattice = lattice, boundary = lattice$segm)
#' mesh <- inla.mesh.create(lattice = lattice, extend = FALSE)
#' plot(mesh)
#'
#' ## Refine the triangulation, with limits on triangle angles and edges:
#' mesh <- inla.mesh.create(
#'   lattice = lattice,
#'   refine = list(max.edge = 0.08),
#'   extend = FALSE
#' )
#' plot(mesh)
#'
#' ## Add an extension around the lattice, but maintain the lattice edges:
#' mesh <- inla.mesh.create(
#'   lattice = lattice,
#'   refine = list(max.edge = 0.08),
#'   interior = lattice$segm
#' )
#' plot(mesh)
#'
#' ## Only add extension:
#' mesh <- inla.mesh.create(lattice = lattice, refine = list(max.edge = 0.08))
#' plot(mesh)
#' @export inla.mesh.lattice
inla.mesh.lattice <- function(x = seq(0, 1, length.out = 2),
                              y = seq(0, 1, length.out = 2),
                              z = NULL,
                              dims =
                                if (is.matrix(x)) {
                                  dim(x)
                                } else {
                                  c(length(x), length(y))
                                },
                              units = NULL,
                              crs = NULL) {
  fmesher_deprecate(
    "soft",
    2L,
    "23.08.18",
    "inla.mesh.lattice()",
    "fmesher::fm_lattice_2d()"
  )
  return(fmesher::fm_lattice_2d(
    x = x,
    y = y,
    z = z,
    dims = dims,
    units = units,
    crs = fm_crs(crs)
  ))
}


#' @title Extract tagged boundary/internal segments.
#'
#' @description
#' `r lifecycle::badge("deprecated")` Use [fmesher::fm_segm()] instead.
#'
#' Extract boundary or internal segments tagged by group id:s.
#'
#' @aliases extract.groups extract.groups.inla.mesh.segment
#' @param segm An [inla.mesh.segment()] object.
#' @param groups The segment groups id:s to extract.
#' @param groups.new Optional vector of group id remapping; `groups[k]` in
#' the input will be replaced by `groups.new[k]` in the output.
#' @param \dots Additional arguments, passed on to other methods.
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @seealso [inla.mesh.segment()]
#' @export extract.groups
extract.groups <- function(segm, groups, groups.new = groups, ...) {
  fmesher_deprecate(
    "soft",
    2L,
    "23.08.18",
    I("`extract.groups(segm, groups)`"),
    I("`fmesher::fm_segm(segm, grp = groups)`"),
    "Note that the `groups.new` argument is not supported."
  )
  return(fmesher::fm_segm(fmesher::fm_as_fm(segm), grp = groups))
}







# Mesh creation ----

inla.mesh <- function(...) {
  args <- list(...)
  if (length(args) > 0) {
    if (inherits(args[[1]], "inla.mesh")) {
      warning("'inla.mesh(mesh, ...)' is deprecated.  Use 'inla.mesh.query(mesh, ...)' instead.")
      return(inla.mesh.query(...))
    }
  }
  stop("'inla.mesh(...)' is deprecated.  Use 'fmesher::fm_mesh_2d_inla(...)' instead.")
  return(inla.mesh.create(...))
}




#' @title Low level function for high-quality triangulations
#'
#' @description
#' `r lifecycle::badge("deprecated")` in favour of [fmesher::fm_rcdt_2d_inla()].
#'
#' Create a constrained refined Delaunay triangulation (CRDT) for a set of
#' spatial locations.
#'
#' `inla.mesh.create` generates triangular meshes on subsets of
#' \eqn{R^2}{R^2} and \eqn{S^2}{S^2}.  Use the higher level wrapper function
#' [inla.mesh.2d()] for greater control over mesh resolution and
#' coarser domain extensions.
#'
#' `inla.delaunay` is a wrapper function for obtaining the convex hull of
#' a point set and calling `inla.mesh.create` to generate the classical
#' Delaunay tringulation.
#'
#' @aliases inla.mesh inla.mesh.create inla.delaunay
#' @param loc Matrix of point locations.  Can alternatively be a
#' `SpatialPoints` or `SpatialPointsDataFrame` object.
#' @param tv A triangle-vertex index matrix, specifying an existing
#' triangulation.
#' @param boundary A list of `inla.mesh.segment` objects, generated by
#' [inla.mesh.segment()], specifying boundary constraint segments.
#' @param interior A list of `inla.mesh.segment` objects, generated by
#' [inla.mesh.segment()], specifying interior constraint segments.
#' @param extend `logical` or `list` specifying whether to extend the
#' data region, with parameters \describe{ \item{list("n")}{the number of edges
#' in the extended boundary (default=8)} \item{list("offset")}{the extension
#' distance.  If negative, interpreted as a factor relative to the approximate
#' data diameter (default=-0.10)} } Setting to `FALSE` is only useful in
#' combination `lattice` or `boundary`.
#' @param refine `logical` or `list` specifying whether to refine the
#' triangulation, with parameters \describe{ \item{list("min.angle")}{the
#' minimum allowed interior angle in any triangle.  The algorithm is guaranteed
#' to converge for `min.angle` at most 21 (default=`21`)}
#' \item{list("max.edge")}{the maximum allowed edge length in any triangle.  If
#' negative, interpreted as a relative factor in an ad hoc formula depending on
#' the data density (default=`Inf`)} \item{list("max.n.strict")}{the
#' maximum number of vertices allowed, overriding `min.angle` and
#' `max.edge` (default=-1, meaning no limit)} \item{list("max.n")}{the
#' maximum number of vertices allowed, overriding `max.edge` only
#' (default=-1, meaning no limit)} }
#' @param lattice An `inla.mesh.lattice` object, generated by
#' [inla.mesh.lattice()], specifying points on a regular lattice.
#' @param globe Subdivision resolution for a semi-regular spherical
#' triangulation with equidistant points along equidistant latitude bands.
#' @param cutoff The minimum allowed distance between points.  Point at most as
#' far apart as this are replaced by a single vertex prior to the mesh
#' refinement step.
#' @param plot.delay On Linux (and Mac if appropriate X11 libraries are
#' installed), specifying a numeric value activates a rudimentary plotting
#' system in the underlying `fmesher` program, showing the triangulation
#' algorithm at work.
#' @param data.dir Where to store the `fmesher` data files.  Defaults to
#' `tempdir()` if `keep` is `FALSE`, otherwise
#' `"inla.mesh.data"`.
#' @param keep `TRUE` if the data files should be kept in `data.dir`
#' or deleted afterwards.  Defaults to true if `data.dir` is specified,
#' otherwise false.  Warning: If `keep` is false, `data.dir` and its
#' contents will be deleted (unless set to `tempdir()`).
#' @param timings If `TRUE`, obtain timings for the mesh construction.
#' @param quality.spec List of vectors of per vertex `max.edge` target
#' specification for each location in `loc`, `boundary/interior`
#' (`segm`), and `lattice`.  Only used if refining the mesh.
#' @param crs An optional `CRS` or `inla.CRS` object
#' @param ...  Optional parameters passed on to `inla.mesh.create`.
#' @return An `inla.mesh` object.
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @seealso [inla.mesh.2d()], [inla.mesh.1d()],
#' [inla.mesh.segment()], [inla.mesh.lattice()],
#' [inla.mesh.query()]
#' @examples
#'
#' loc <- matrix(runif(10 * 2), 10, 2)
#'
#' mesh <- inla.delaunay(loc)
#' plot(mesh)
#'
#' mesh <- inla.mesh.create(loc,
#'   interior = inla.mesh.segment(idx = 1:2),
#'   extend = TRUE,
#'   refine = list(max.edge = 0.1)
#' )
#' plot(mesh)
#'
#' loc2 <- matrix(c(0, 1, 1, 0, 0, 0, 1, 1), 4, 2)
#' mesh2 <- inla.mesh.create(
#'   loc = loc,
#'   boundary = inla.mesh.segment(loc2),
#'   interior = inla.mesh.segment(idx = 1:2),
#'   quality.spec = list(segm = 0.2, loc = 0.05),
#'   refine = list(min.angle = 26)
#' )
#' plot(mesh2)
#' @export inla.mesh.create
#' @importFrom fmesher fm_CRS fm_crs
inla.mesh.create <- function(loc = NULL, tv = NULL,
                             boundary = NULL, interior = NULL,
                             extend = (missing(tv) || is.null(tv)),
                             refine = FALSE,
                             lattice = NULL,
                             globe = NULL,
                             cutoff = 1e-12,
                             plot.delay = NULL,
                             data.dir = NULL,
                             keep = (!missing(data.dir) && !is.null(data.dir)),
                             timings = FALSE,
                             quality.spec = NULL,
                             crs = NULL) {
  fmesher_deprecate(
    "soft",
    2L,
    "23.08.18",
    "inla.mesh.create()",
    "fmesher::fm_rcdt_2d_inla()"
  )
  return(fmesher::fm_rcdt_2d_inla(
    loc = loc,
    tv = tv,
    boundary = boundary,
    interior = interior,
    extend = extend,
    refine = refine,
    lattice = lattice,
    globe = globe,
    cutoff = cutoff,
    quality.spec = quality.spec,
    crs = crs
  ))
}



#' @title Constraint segment extraction for inla.mesh
#'
#' @description
#' `r lifecycle::badge("deprecated")` Use [fmesher::fm_segm()] instead.
#'
#' Constructs an list of `inla.mesh.segment` object from boundary or
#' interior constraint information in an [inla.mesh()] object.
#'
#'
#' @aliases inla.mesh.boundary inla.mesh.interior
#' @param mesh An `inla.mesh` object.
#' @param grp Group indices to extract.  If `NULL`, all boundary/interior
#' constrain groups are extracted.
#' @return A list of `inla.mesh.segment` objects.
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @seealso [inla.mesh.segment()], [inla.mesh.create()],
#' [inla.mesh.create.helper()]
#' @examples
#'
#' loc <- matrix(runif(100 * 2) * 1000, 100, 2)
#' mesh <- fmesher::fm_mesh_2d_inla(loc.domain = loc, max.edge = c(50, 500))
#' boundary <- inla.mesh.boundary(mesh)
#' interior <- inla.mesh.interior(mesh)
#' @export
#' @rdname inla.mesh.boundary
inla.mesh.boundary <- function(mesh, grp = NULL) {
  fmesher_deprecate(
    "soft",
    2L,
    "23.08.18",
    "inla.mesh.boundary()",
    I("`fmesher::fm_segm(mesh, boundary = TRUE)`")
  )
  return(list(fmesher::fm_segm(fmesher::fm_as_fm(mesh),
    boundary = TRUE,
    grp = grp
  )))
}

#' @export
#' @rdname inla.mesh.boundary
inla.mesh.interior <- function(mesh, grp = NULL) {
  fmesher_deprecate(
    "soft",
    2L,
    "23.08.18",
    "inla.mesh.interior()",
    I("`fmesher::fm_segm(mesh, boundary = FALSE)`")
  )
  return(list(fmesher::fm_segm(fmesher::fm_as_fm(mesh),
    boundary = FALSE,
    grp = grp
  )))
}


## Generate nice triangulation, with an inner domain strictly enclosed
## by an outer domain
## The inner and outer domains can have different quality parameters
## At least one of loc, loc.domain, boundary[[1]], boundary[[2]], interior
## must be non-NULL
## For more complicated multi-step meshings, study the code and write your own.


#' @title High-quality triangulations
#'
#' @description
#' `r lifecycle::badge("deprecated")` Use [fmesher::fm_mesh_2d_inla()] instead.
#'
#' Create a triangle mesh based on initial point locations, specified or
#' automatic boundaries, and mesh quality parameters.
#'
#'
#' @param loc Matrix of point locations to be used as initial triangulation
#' nodes.  Can alternatively be a `SpatialPoints` or
#' `SpatialPointsDataFrame` object.
#' @param loc.domain Matrix of point locations used to determine the domain
#' extent.  Can alternatively be a `SpatialPoints` or
#' `SpatialPointsDataFrame` object.
#' @param offset The automatic extension distance.  One or two values, for an
#' inner and an optional outer extension.  If negative, interpreted as a factor
#' relative to the approximate data diameter (default=-0.10???)
#' @param n The number of initial nodes in the automatic extensions
#' (default=16)
#' @param boundary A list of one or two [inla.mesh.segment()] objects
#' describing domain boundaries.
#' @param interior An [inla.mesh.segment()] object describing desired
#' interior edges.
#' @param max.edge The largest allowed triangle edge length.  One or two
#' values.
#' @param min.angle The smallest allowed triangle angle.  One or two values.
#' (Default=21)
#' @param cutoff The minimum allowed distance between points.  Point at most as
#' far apart as this are replaced by a single vertex prior to the mesh
#' refinement step.
#' @param max.n.strict The maximum number of vertices allowed, overriding
#' `min.angle` and `max.edge` (default=-1, meaning no limit).  One or
#' two values, where the second value gives the number of additional vertices
#' allowed for the extension.
#' @param max.n The maximum number of vertices allowed, overriding
#' `max.edge` only (default=-1, meaning no limit).  One or two values,
#' where the second value gives the number of additional vertices allowed for
#' the extension.
#' @param plot.delay On Linux (and Mac if appropriate X11 libraries are
#' installed), specifying a nonnegative numeric value activates a rudimentary
#' plotting system in the underlying `fmesher` program, showing the
#' triangulation algorithm at work, with waiting time factor `plot.delay`
#' between each step.
#'
#' On all systems, specifying any negative value activates displaying the
#' result after each step of the multi-step domain extension algorithm.
#' @param crs An optional `CRS` or `inla.CRS` object
#' @return An `inla.mesh` object.
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @seealso [inla.mesh.create()], [inla.delaunay()],
#' [inla.nonconvex.hull()]
#' @examples
#'
#' loc <- matrix(runif(10 * 2), 10, 2)
#'
#' if (require("splancs")) {
#'   boundary <- list(
#'     inla.nonconvex.hull(loc, 0.1, 0.15),
#'     inla.nonconvex.hull(loc, 0.2, 0.2)
#'   )
#'   offset <- NULL
#' } else {
#'   boundary <- NULL
#'   offset <- c(0.1, 0.2)
#' }
#' mesh <- inla.mesh.2d(loc, boundary = boundary, offset = offset, max.edge = c(0.05, 0.1))
#'
#' plot(mesh)
#' @importFrom fmesher fm_crs_is_geocent fm_crs_is_identical
#' @export inla.mesh.2d
inla.mesh.2d <- function(loc = NULL, ## Points to include in final triangulation
                         loc.domain = NULL, ## Points that determine the automatic domain
                         offset = NULL, ## Size of automatic extensions
                         n = NULL, ## Sides of automatic extension polygons
                         boundary = NULL, ## User-specified domains (list of length 2)
                         interior = NULL, ## User-specified constraints for the inner domain
                         max.edge = NULL,
                         min.angle = NULL, ## Angle constraint for the entire domain
                         cutoff = 1e-12, ## Only add input points further apart than this
                         max.n.strict = NULL,
                         max.n = NULL,
                         plot.delay = NULL,
                         crs = NULL) {
  fmesher_deprecate(
    "soft",
    2L,
    "23.08.18",
    "inla.mesh.2d()",
    "fmesher::fm_mesh_2d_inla()"
  )
  return(fmesher::fm_mesh_2d_inla(
    loc = loc,
    loc.domain = loc.domain,
    offset = offset,
    n = n,
    boundary = boundary,
    interior = interior,
    max.edge = max.edge,
    min.angle = min.angle,
    cutoff = cutoff,
    max.n.strict = max.n.strict,
    max.n = max.n,
    crs = crs
  ))
}

## Support for legacy code:


#' @title High-quality triangulations
#'
#' @description
#' `r lifecycle::badge("deprecated")` Use [fmesher::fm_mesh_2d_inla()] instead.
#'
#' Create a triangle mesh based on initial point locations, specified or
#' automatic boundaries, and mesh quality parameters.
#'
#'
#' @param points Matrix of point locations to be used as initial triangulation
#' nodes.
#' @param points.domain Matrix of point locations used to determine the domain
#' extent.
#' @param \dots Additional arguments passed on to [inla.mesh.2d()].
#' @note Obsolete function kept for backwards compatibility.  Use
#' [fmesher::fm_mesh_2d_inla()] instead.
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @seealso [inla.mesh.2d()]
#' @keywords internal
#' @export inla.mesh.create.helper
inla.mesh.create.helper <- function(points = NULL, points.domain = NULL, ...) {
  lifecycle::deprecate_warn(
    "a very long time ago",
    "inla.mesh.create.helper()",
    "fmesher::fm_mesh_2d_inla()"
  )
  return(invisible(inla.mesh.2d(loc = points, loc.domain = points.domain, ...)))
}

## Example:
## if (FALSE) {
## mesh =
## inla.mesh.create.helper(loc=matrix(runif(20),10,2)*200,
## n=c(8,16),
## offset=c(10,140),
## max.edge=c(25,1000),
## min.angle=26,
## cutoff=0,
## plot.delay=-1
## )
## }
##

#' @export
#' @describeIn inla.mesh.create
#' `r lifecycle::badge("deprecated")` Use [fmesher::fm_delaunay_2d()] instead.
inla.delaunay <- function(loc, ...) {
  fmesher_deprecate(
    "soft",
    2L,
    "23.08.18",
    "inla.delaunay()",
    "fmesher::fm_delaunay_2d()"
  )
  return(fmesher::fm_delaunay_2d(loc, ...))
}





# Queries ----


#' High-quality triangulations
#'
#' Query information about an inla.mesh object.
#'
#'
#' @param mesh An `inla.mesh` object.
#' @param \dots Query arguments.
#' \itemize{
#' \item tt.neighbours Compute
#' neighbour triangles for triangles; list of vectors: list(triangles, orders)
#'
#' \item vt.neighbours Compute neighbour triangles for vertices; list of
#' vectors: list(vertices, orders)
#' }
#' @return A list of query results.
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @seealso [inla.mesh.create()], [inla.mesh.segment()],
#' [inla.mesh.lattice()]
#' @examples
#'
#' loc <- matrix(c(0.1, 0.15), 1, 2)
#' lattice <- inla.mesh.lattice(
#'   seq(0, 1, length.out = 10),
#'   seq(0, 1, length.out = 10)
#' )
#' mesh <- inla.mesh.create(loc = loc, lattice = lattice, extend = FALSE)
#'
#' vt <- which(inla.mesh.query(mesh,
#'   vt.neighbours = list(
#'     mesh$idx$loc,
#'     4:6
#'   )
#' )$vt.neighbours)
#'
#' mesh2 <- inla.mesh.create(mesh$loc,
#'   tv = mesh$graph$tv[vt, , drop = FALSE],
#'   refine = FALSE, extend = FALSE
#' )
#' @export
inla.mesh.query <- function(mesh, ...) {
  fmesher_deprecate("soft",
    2L,
    "23.08.18",
    "inla.mesh.query()",
    details = "Does not yet have an `fmesher` alternative."
  )

  inla.require.inherits(mesh, "inla.mesh", "'mesh'")

  not.known <- function(mesh, queryname) {
    stop(paste("Query '", queryname,
      "' unknown.",
      sep = ""
    ))
  }
  not.implemented <- function(mesh, queryname) {
    stop(paste("Query '", queryname,
      "' not implemented for inla.mesh.",
      sep = ""
    ))
    ## stop(paste("Query '", queryname,
    ## "' not implemented for inla.mesh for mesh type '",
    ## model$type, "'.", sep=""))
  }

  result <- list()
  queries <- inla.parse.queries(...)
  if (length(queries) == 0L) {
    return(result)
  }

  for (query.idx in 1:length(queries)) {
    query <- names(queries)[query.idx]
    param <- queries[[query.idx]]
    answer <- NULL
    query <- (match.arg(query, c(
      "tt.neighbours",
      "vt.neighbours"
    )))
    if (identical(query, "tt.neighbours")) {
      if (is.null(param)) {
        param <- list(c(1))
      }
      if (length(param) < 2) {
        param <- c(as.list(param), list(c(1)))
      }
      nT <- nrow(mesh$graph$tv)
      i <- rep(1:nT, 3)
      j <- as.vector(mesh$graph$tt)
      i <- i[!is.na(j)]
      j <- j[!is.na(j)]
      tt <- sparseMatrix(i = i, j = j, x = rep(1, length(i)), dims = c(nT, nT))
      answer <- (sparseMatrix(
        i = param[[1]],
        j = rep(1, length(param[[1]])),
        x = 1,
        dims = c(nT, 1)
      ))
      tt0 <- (answer == 0.5) * 1
      tt1 <- answer
      order <- 0
      while (order < min(param[[2]])) {
        order <- order + 1
        tt0 <- tt1
        tt1 <- answer
        answer <- ((((tt %*% answer) > 0) - tt1 - tt0) > 0)
      }
      while (order < max(param[[2]])) {
        order <- order + 1
        answer <- ((((tt %*% answer) > 0) - tt1 - tt0) > 0)
      }

      ## not.implemented(mesh,query)
    } else if (identical(query, "vt.neighbours")) {
      if (is.null(param)) {
        param <- list(1)
      }
      if (length(param) < 2) {
        param <- c(as.list(param), list(c(1)))
      }
      nV <- nrow(mesh$loc)
      nT <- nrow(mesh$graph$tv)
      i <- rep(1:nT, 3)
      j <- as.vector(mesh$graph$tv)
      # i = i[!is.na(j)]
      # j = j[!is.na(j)]
      tv <- sparseMatrix(i = i, j = j, x = rep(1, length(i)), dims = c(nT, nV))
      vv <- (sparseMatrix(
        i = param[[1]],
        j = rep(1, length(param[[1]])),
        x = 1,
        dims = c(nV, 1)
      ))
      vt <- (tv %*% vv) > 0
      vv0 <- (vv == 0.5) * 1
      vv1 <- vv
      vt0 <- (tv %*% vv0) > 0
      vt1 <- vt
      order <- 0
      while (order < min(param[[2]])) {
        order <- order + 1
        vv0 <- vv1
        vv1 <- vv
        vv <- ((((mesh$graph$vv %*% vv) > 0) - vv1 - vv0) > 0)
        vt0 <- vt1
        vt1 <- vt
        vt <- (((tv %*% vv) > 0) - vt1 - vt0) > 0
      }
      while (order < max(param[[2]])) {
        order <- order + 1
        vv <- ((((mesh$graph$vv %*% vv) > 0) - vv1 - vv0) > 0)
        vt <- (((((tv %*% vv) > 0) - vt1 - vt0) > 0) + vt) > 0
      }
      answer <- vt

      ## not.implemented(mesh,query)
    } else if (!identical(query, "")) {
      not.known(mesh, query)
    }
    ## Expand the result list:
    result[query.idx] <- list(NULL)
    names(result)[query.idx] <- query
    ## Set the answer:
    if (!is.null(answer)) {
      result[[query.idx]] <- answer
    }
  }

  return(result)
}



#' Summarizing triangular mesh objects
#'
#' Construct and print `inla.mesh` object summaries
#'
#'
#' @aliases summary.inla.mesh print.summary.inla.mesh
#' @param object an object of class `"inla.mesh"`, usually a result of a
#' call to [inla.mesh.create()] or [inla.mesh.2d()].
#' @param x an object of class `"summary.inla.mesh"`, usually a result of
#' a call to [summary.inla.mesh()].
#' @param verbose If `TRUE`, produce a more detailed output.
#' @param \dots further arguments passed to or from other methods.
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @method summary inla.mesh
#' @export
summary.inla.mesh <- function(object, verbose = FALSE, ...) {
  x <- object
  ## provides a summary for a mesh object
  inla.require.inherits(x, "inla.mesh", "'x'")

  ret <- list(verbose = verbose)
  if (verbose) {
    ret <- (c(ret, list(
      call = x$meta$call,
      fmesher.args = x$meta$fmesher.args,
      prefix = x$meta$prefix,
      time = x$meta$time,
      is.refined = x$meta$is.refined
    )))
  }
  ret <- (c(ret, list(
    manifold = x$manifold,
    nV = x$n,
    nT = nrow(x$graph$tv),
    xlim = range(x$loc[, 1]),
    ylim = range(x$loc[, 2]),
    zlim = if (ncol(x$loc) >= 3) range(x$loc[, 3]) else c(NA_real_, NA_real_)
  )))
  ret <- c(ret, list(crs = as.character(fmesher::fm_wkt(x$crs))))
  ret <- c(ret, list(crs_proj4 = as.character(fmesher::fm_crs(x$crs)$proj4string)))

  my.segm <- function(x) {
    if (is.null(x)) {
      return(list(n = 0, grps = NULL))
    }
    n <- max(0, nrow(x$idx))
    if (max(0, length(unique(x$grp))) > 0) {
      grps <- unique(x$grp)
    } else {
      grps <- NULL
    }
    return(list(n = n, grps = grps))
  }
  if (!is.null(x$segm)) {
    ret <- (c(ret, list(
      segm.bnd = my.segm(x$segm$bnd),
      segm.int = my.segm(x$segm$int)
    )))
  } else {
    ret <- (c(ret, list(
      segm.bnd = my.segm(NULL),
      segm.int = my.segm(NULL)
    )))
  }

  class(ret) <- "summary.inla.mesh"
  return(ret)
}

#' @export
#' @method print summary.inla.mesh
#' @rdname summary.inla.mesh
print.summary.inla.mesh <- function(x, ...) {
  my.print.proc_time <- function(x, ...) {
    if (!is.matrix(x)) {
      y <- matrix(x, 1, 5)
    } else {
      y <- x
    }
    for (k in 1:nrow(x)) {
      if (!is.na(y[k, 4L])) {
        y[k, 1L] <- y[k, 1L] + y[k, 4L]
      }
      if (!is.na(y[k, 5L])) {
        y[k, 2L] <- y[k, 2L] + y[k, 5L]
      }
    }
    y <- y[, 1L:3L]
    colnames(y) <- c(gettext("user"), gettext("system"), gettext("elapsed"))
    print(y, ...)
    invisible(x)
  }


  inla.require.inherits(x, "summary.inla.mesh", "'x'")

  if (x$verbose) {
    cat("\nCall:\n")
    print(x$call)

    cat("\nfmesher:\t", x$fmesher.args, "\n", sep = "")
    cat("prefix:\t\t", x$prefix, "\n", sep = "")

    cat("\nTimings:\n")
    my.print.proc_time(x$time)
  }

  cat("\nManifold:\t", x$manifold, "\n", sep = "")
  cat("CRS/LegacyPROJ4:\t", x$crs_proj4, "\n", sep = "")
  if (x$verbose) {
    cat("CRS/WKT:\n", x$crs, "\n", sep = "")
  } else {
    cat("CRS/WKT: (only shown with verbose = TRUE)", "\n", sep = "")
  }
  if (x$verbose) {
    cat("Refined:\t", x$is.refined, "\n", sep = "")
  }
  cat("Vertices:\t", as.character(x$nV), "\n", sep = "")
  cat("Triangles:\t", as.character(x$nT), "\n", sep = "")

  my.print.segm <- function(x) {
    cat(as.character(x$n))
    if (!is.null(x$grps)) {
      n <- length(x$grps)
      cat(" (", n, " group", inla.ifelse(n == 1, "", "s"), sep = "")
      if (n <= 10) {
        cat(":", x$grps, sep = " ")
      } else {
        cat(":", x$grps[1:10], "...", sep = " ")
      }
      cat(")")
    }
    cat("\n", sep = "")
    return(invisible())
  }
  cat("Boundary segm.:\t")
  my.print.segm(x$segm.bnd)
  cat("Interior segm.:\t")
  my.print.segm(x$segm.int)
  cat("xlim:\t", x$xlim[1], " ", x$xlim[2], "\n", sep = "")
  cat("ylim:\t", x$ylim[1], " ", x$ylim[2], "\n", sep = "")
  cat("zlim:\t", x$zlim[1], " ", x$zlim[2], "\n", sep = "")
  cat("\n")

  invisible(x)
}



# Point/mesh connection methods ####

#' @title Methods for projecting to/from an inla.mesh
#'
#' @description
#' `r lifecycle::badge("deprecated")` Use [fmesher::fm_evaluate()] and
#' [fmesher::fm_evaluator()] instead.
#'
#' Calculate a lattice projection to/from an [inla.mesh()].
#'
#' The call `inla.mesh.project(mesh, loc, field=..., ...)`, is a shortcut
#' to `inla.mesh.project(inla.mesh.projector(mesh, loc), field)`.
#'
#' @aliases inla.mesh.project inla.mesh.projector inla.mesh.projector.inla.mesh
#' inla.mesh.project.inla.mesh inla.mesh.project.inla.mesh.projector
#' inla.mesh.project.inla.mesh.1d inla.mesh.projector.inla.mesh.1d
#' @param \dots Arguments passed on to [fmesher::fm_evaluate()] and
#' [fmesher::fm_evaluator()].
#' @return For `inla.mesh.project(mesh, ...)`, a list with projection
#' information.  For `inla.mesh.projector(mesh, ...)`, an
#' `inla.mesh.projector` object.  For `inla.mesh.project(projector,
#' field, ...)`, a field projected from the mesh onto the locations given by
#' the projector object.
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @seealso [inla.mesh()], [inla.mesh.1d()],
#' [inla.mesh.lattice()]
#' @examples
#'
#' n <- 20
#' loc <- matrix(runif(n * 2), n, 2)
#' mesh <- inla.mesh.create(loc, refine = list(max.edge = 0.05))
#' proj <- inla.mesh.projector(mesh)
#' field <- cos(mesh$loc[, 1] * 2 * pi * 3) * sin(mesh$loc[, 2] * 2 * pi * 7)
#' image(proj$x, proj$y, inla.mesh.project(proj, field))
#' \donttest{
#' if (require(rgl)) {
#'   plot(mesh, rgl = TRUE, col = field, draw.edges = FALSE, draw.vertices = FALSE)
#' }
#' }
#'
#' @export inla.mesh.project
inla.mesh.project <- function(...) {
  fmesher_deprecate(
    "warn",
    "23.06.07",
    "inla.mesh.project()",
    "fmesher::fm_evaluate()"
  )
  return(fmesher::fm_evaluate(...))
}


#' @export
#' @rdname inla.mesh.project
inla.mesh.projector <- function(...) {
  fmesher_deprecate(
    "warn",
    "23.06.07",
    "inla.mesh.projector()",
    "fmesher::fm_evaluator()"
  )
  return(fmesher::fm_evaluator(...))
}




# Misc ####

#' @title Basis functions for inla.mesh
#'
#' @description
#' `r lifecycle::badge("deprecated")` Use [fmesher::fm_raw_basis()] instead.
#'
#' Calculate basis functions on a 1d or 2d [inla.mesh()]
#'
#' @param mesh An `inla.mesh.1d` or `inla.mesh` object.
#' @param type `b.spline` (default) for B-spline basis functions,
#' `sph.harm` for spherical harmonics (available only for meshes on the
#' sphere)
#' @param n For B-splines, the number of basis functions in each direction (for
#' 1d meshes `n` must be a scalar, and for planar 2d meshes a 2-vector).
#' For spherical harmonics, `n` is the maximal harmonic order.
#' @param degree Degree of B-spline polynomials.  See
#' [inla.mesh.1d()].
#' @param knot.placement For B-splines on the sphere, controls the latitudinal
#' placements of knots. `"uniform.area"` (default) gives uniform spacing
#' in `sin(latitude)`, `"uniform.latitude"` gives uniform spacing in
#' latitudes.
#' @param rot.inv For spherical harmonics on a sphere, `rot.inv=TRUE`
#' gives the rotationally invariant subset of basis functions.
#' @param boundary Boundary specification, default is free boundaries.  See
#' [inla.mesh.1d()] for more information.
#' @param free.clamped If `TRUE` and `boundary` is `"free"`, the
#' boundary basis functions are clamped to 0/1 at the interval boundary by
#' repeating the boundary knots.
#' @param ... Unused
#'
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @seealso [inla.mesh.1d()] [inla.mesh.2d()]
#' @examples
#'
#' n <- 100
#' loc <- matrix(runif(n * 2), n, 2)
#' mesh <- inla.mesh.2d(loc, max.edge = 0.05)
#' basis <- inla.mesh.basis(mesh, n = c(4, 5))
#'
#' proj <- inla.mesh.projector(mesh)
#' image(proj$x, proj$y, inla.mesh.project(proj, basis[, 7]))
#' \donttest{
#' if (require(rgl)) {
#'   plot(mesh, rgl = TRUE, col = basis[, 7], draw.edges = FALSE, draw.vertices = FALSE)
#' }
#' }
#'
#' @export inla.mesh.basis
inla.mesh.basis <- function(mesh,
                            type = "b.spline",
                            n = 3,
                            degree = 2,
                            knot.placement = "uniform.area",
                            rot.inv = TRUE,
                            boundary = "free",
                            free.clamped = TRUE,
                            ...) {
  fmesher_deprecate(
    "soft",
    2L,
    "23.08.18",
    "inla.mesh.basis()",
    "fmesher::fm_raw_basis()"
  )
  return(
    fmesher::fm_raw_basis(
      mesh,
      type = type,
      n = n,
      degree = degree,
      knot.placement = knot.placement,
      rot.inv = rot.inv,
      boundary = boundary,
      free.clamped = free.clamped,
      ...
    )
  )
}


inla.parse.queries <- function(...) {
  queries <- list(...)
  if (length(queries) == 0) {
    return(queries)
  }

  ## Make sure that we have a list of names, empty or not:
  if (is.null(names(queries))) {
    q.names <- rep("", length(queries))
  } else {
    q.names <- names(queries)
  }

  ## All nameless entries must be strings with query names.  Replace
  ## empty names with those names, and set those entries to NULL.
  for (query.idx in 1:length(queries)) {
    if (q.names[[query.idx]] == "") {
      if (is.character(queries[[query.idx]])) {
        names(queries)[query.idx] <- queries[[query.idx]]
        queries[query.idx] <- list(NULL)
      } else {
        queries[query.idx] <- list(NULL)
        warning(paste("Unnamed query ignored.  Check query #",
          query.idx,
          sep = ""
        ))
      }
    }
  }

  return(queries)
}


#' @title Compute various mesh related quantities.
#'
#' @description
#' `r lifecycle::badge("deprecated")` Use the methods in the `fmesher` package
#' instead; see details below.
#'
#' Low level function for computing finite element matrices, spherical
#' harmonics, B-splines, and point mappings with barycentric triangle
#' coordinates.
#'
#' @param loc 3-column triangle vertex coordinate matrix.
#' @param tv 3-column triangle vertex index matrix.
#' @param fem `r lifecycle::badge("deprecated")` Use [fmesher::fm_fem()] instead.
#' Maximum finite element matrix order to be computed.
#' @param aniso `r lifecycle::badge("deprecated")` Use [fmesher::fm_fem()] instead.
#' A two-element list with \eqn{\gamma}{gamma} and \eqn{v}{v} for
#' an anisotropic operator \eqn{\nabla\cdot H \nabla}{div H grad}, where
#' \eqn{H=\gamma I + v v^\top}{H = gamma I + v v'}
#' @param gradients `r lifecycle::badge("deprecated")` Use [fmesher::fm_fem()] instead.
#' When `TRUE`, calculate derivative operator matrices
#' `dx`, dy, and dz.
#' @param sph0 `r lifecycle::badge("deprecated")` Use [fmesher::fm_raw_basis()] instead.
#' @param sph `r lifecycle::badge("deprecated")` Use [fmesher::fm_raw_basis()] instead.
#' @param bspline `r lifecycle::badge("deprecated")` Use [fmesher::fm_raw_basis()] instead.
#' Rotationally invariant B-splines on a sphere.  3-vector with
#' number of basis functions `n`, basis degree `degree`, and a
#' logical; `TRUE` uniform knot angles, `FALSE` for uniform spacing
#' in \eqn{\sin(latitude)}{sin(latitude)}.
#' @param points2mesh `r lifecycle::badge("deprecated")` Use [fmesher::fm_bary()] instead.
#' 3-column matrix with points to be located in the mesh.
#' @param splitlines `r lifecycle::badge("deprecated")` Use [fmesher::fm_split_lines()] or
#' [fmesher::fmesher_split_lines()] instead.
#' A list with elements `loc` (3-column coordinate
#' matrix) and `idx` (2-column index matrix) describing line segments that
#' are to be split into sub-segments at triangle boundaries.
#' @param output Names of objects to be included in the output, if different
#' from defaults.
#' @param keep When `TRUE`, for debugging purposes keep the `fmesher`
#' I/O files on disk.
#' @return A list of generated named quantities.
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @export inla.fmesher.smorg
`inla.fmesher.smorg` <- function(loc, tv,
                                 fem = NULL,
                                 aniso = NULL,
                                 gradients = FALSE,
                                 sph0 = deprecated(),
                                 sph = deprecated(),
                                 bspline = NULL,
                                 points2mesh = NULL,
                                 splitlines = NULL,
                                 output = NULL,
                                 keep = FALSE) {
  msg <-
    c(
      fem = "`fmesher::fm_fem(..., order = ...)`",
      aniso = "`fmesher::fm_fem(..., aniso = ...)`",
      gradients = "`fmesher::fm_basis(..., derivatives = TRUE)`",
      bspline = "fmesher::fmesher_spherical_bsplines()",
      points2mesh = "`fmesher::fm_bary()`",
      splitlines = paste0(
        "`fmesher::fm_split_lines()` for a high-level ",
        "interface, or `fmesher::fmesher_split_lines()` ",
        "if you need the low-level outputs."
      )
    )
  given <- !c(
    fem = is.null(fem),
    aniso = is.null(aniso),
    gradients = isTRUE(gradients),
    bspline = is.null(bspline),
    points2mesh = is.null(points2mesh),
    splitlines = is.null(splitlines)
  )
  msg <- paste("Use ",
    msg[names(given)[given]],
    " instead of `inla.fmesher.smorg(..., ",
    names(given)[given],
    " = ...)`.",
    sep = ""
  )
  names(msg) <- NULL
  fmesher_deprecate("warn",
    2L,
    "23.08.18",
    "inla.fmesher.smorg()",
    details = msg
  )

  output <- list()
  if (!is.null(fem)) {
    output <- c(
      output,
      fmesher::fm_fem(
        fmesher::fm_rcdt_2d_inla(loc = loc, tv = tv),
        order = fem,
        aniso = aniso
      )
    )
  }
  if (isTRUE(gradients)) {
    res <- fmesher::fm_basis(
      fmesher::fm_rcdt_2d_inla(loc = loc, tv = tv),
      loc = loc,
      derivatives = TRUE
    )
    output[["dx"]] <- res[["dx"]]
    output[["dy"]] <- res[["dy"]]
    output[["dz"]] <- res[["dz"]]
  }
  if (!is.null(bspline)) {
    output[["bspline"]] <-
      fmesher::fmesher_spherical_bsplines(
        loc,
        n = bspline[1],
        degree = bspline[2],
        uniform = bspline[3] > 0
      )
  }
  if (!is.null(points2mesh)) {
    res <- fmesher::fm_bary(
      fmesher::fm_rcdt_2d_inla(loc = loc, tv = tv),
      loc = points2mesh
    )
    res[["t"]][is.na(res[["t"]])] <- 0L
    res[["bary"]][is.na(res[["bary"]])] <- 0
    output[["p2m.t"]] <- res[["t"]]
    output[["p2m.b"]] <- res[["bary"]]
  }
  if (!is.null(splitlines)) {
    splt <- fmesher::fmesher_split_lines(
      mesh_loc = loc,
      mesh_tv = tv - 1L,
      loc = splitlines$loc,
      idx = splitlines$idx - 1L,
      options = list()
    )
    indexoutput <- list("split.idx", "split.t", "split.origin")
    for (name in intersect(indexoutput, names(splt))) {
      splt[[name]] <- splt[[name]] + 1L
    }
    output <- c(output, splt)
  }
  return(output)
}


# 1D mesh creation ----


#' @title Function space definition objects for 1D SPDE models.
#'
#' @description
#' `r lifecycle::badge("deprecated")` Use [fmesher::fm_mesh_1d()] instead.
#'
#' Create a 1D mesh specification `inla.mesh.1d` object, that defines a
#' function space for 1D SPDE models.
#'
#'
#' @param loc B-spline knot locations.
#' @param interval Interval domain endpoints.
#' @param boundary Boundary condition specification.  Valid conditions are
#' `c('neumann', 'dirichlet', 'free', 'cyclic')`.  Two separate values can
#' be specified, one applied to each endpoint.
#' @param degree The B-spline basis degree.  Supported values are 0, 1, and 2.
#' @param free.clamped If `TRUE`, for `'free'` boundaries, clamp the
#' basis functions to the interval endpoints.
#' @param \dots Additional option, currently unused.
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @export inla.mesh.1d
inla.mesh.1d <- function(loc,
                         interval = range(loc),
                         boundary = NULL,
                         degree = 1,
                         free.clamped = FALSE,
                         ...) {
  fmesher_deprecate(
    "soft",
    2L,
    "23.08.18",
    "inla.mesh.1d()",
    "fmesher::fm_mesh_1d()"
  )
  return(
    fmesher::fm_mesh_1d(
      loc = loc,
      interval = interval,
      boundary = boundary,
      degree = degree,
      free.clamped = free.clamped,
      ...
    )
  )
}


# 1D mesh queries ----

#' @export inla.mesh.1d.bary
#' @describeIn inla.mesh.1d.A `r lifecycle::badge("deprecated")`
#' Use [fmesher::fm_bary()] instead.
inla.mesh.1d.bary <- function(mesh, loc, method = c("linear", "nearest")) {
  fmesher_deprecate(
    "soft",
    2L,
    "23.08.18",
    "inla.mesh.1d.bary()",
    "fmesher::fm_bary()",
    details =
      c(
        "An `index` field his being included for backwards compatibility.",
        "The canconical element from `fm_bary()` is `t`."
      )
  )
  result <- fmesher::fm_bary(
    fmesher::fm_as_mesh_1d(mesh),
    loc,
    method = method
  )
  result$index <- result$t
  return(result)
}



#' @title Mapping matrix for 1D meshes
#'
#' @description
#' `r lifecycle::badge("deprecated")` Use [fmesher::fm_basis()] instead.
#'
#' Calculates barycentric coordinates and weight matrices for
#' [inla.mesh.1d()] objects.
#'
#'
#' @aliases inla.mesh.1d.A inla.mesh.1d.bary
#' @param mesh An [inla.mesh.1d()] object.
#' @param loc Coordinate values.
#' @param weights Weights to be applied to the `A` matrix rows.
#' @param derivatives If `TRUE`, also compute derivative weight matrices
#' `dA` and `d2A`.
#' @param method Interpolation method.  If not specified for
#' `inla.mesh.1d.A` (recommended), it is determined by the `mesh`
#' basis function properties.
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @export inla.mesh.1d.A
inla.mesh.1d.A <- function(mesh, loc,
                           weights = NULL,
                           derivatives = NULL,
                           method = NULL) {
  fmesher_deprecate(
    "soft",
    2L,
    "23.08.18",
    "inla.mesh.1d.A()",
    "fmesher::fm_basis()"
  )
  return(fmesher::fm_basis(
    fmesher::fm_as_mesh_1d(mesh),
    loc,
    weights = weights,
    derivatives = derivatives,
    method = method
  ))
}

#' @param mesh An inla.mesh.1d object
#' @export
#' @describeIn inla.mesh.1d `r lifecycle::badge("deprecated")`
#' Use [fmesher::fm_fem()] instead.
inla.mesh.1d.fem <- function(mesh) {
  fmesher_deprecate(
    "soft",
    2L,
    "23.08.18",
    "inla.mesh.fem()",
    "fmesher::fm_fem()"
  )
  return(fmesher::fm_fem(fmesher::fm_as_mesh_1d(mesh), order = 2))
}


# Queries ----

#' @title Diameter of a point set
#'
#' @description
#' `r lifecycle::badge("deprecated")` Use [fmesher::fm_diameter()] instead.
#'
#' Find an upper bound to the convex hull of a point set
#'
#'
#' @aliases inla.diameter inla.diameter.default inla.diameter.inla.mesh
#' inla.diameter.inla.mesh.segment inla.diameter.inla.mesh.lattice
#' inla.diameter.inla.mesh.1d
#' @param x A point set as an \eqn{n\times d}{n x d} matrix, or an
#' `inla.mesh` related object.
#' @param manifold Character string specifying the manifold type. Default is to
#' treat the point set with Euclidean \eqn{R^d} metrics. Use
#' `manifold="S2"` for great circle distances on the unit sphere (this is
#' set automatically for `inla.mesh` objects).
#' @param \dots Additional parameters passed on to other methods.
#' @return A scalar, upper bound for the diameter of the convex hull of the
#' point set.
#' @author Finn Lindgren <finn.lindgren@@gmail.com>
#' @examples
#'
#' inla.diameter(matrix(c(0, 1, 1, 0, 0, 0, 1, 1), 4, 2))
#' @export inla.diameter
inla.diameter <- function(x, ...) {
  fmesher_deprecate(
    "soft",
    2L,
    "23.08.18",
    "inla.diameter()",
    "fmesher::fm_diameter()"
  )
  return(fmesher::fm_diameter(x, ...))
}




#' @title Finite element matrices
#'
#' @description
#' `r lifecycle::badge("deprecated")` Use [fmesher::fm_fem()] instead.
#'
#' Constructs finite element matrices for [inla.mesh()] and
#' [inla.mesh.1d()] objects.
#'
#'
#' @aliases inla.mesh.fem
#' @param mesh An [inla.mesh()] or [inla.mesh.1d()] object.
#' @param order The model order.
#' @return A list of sparse matrices based on basis functions `psi_i`:
#' \item{c0 }{`c0[i,j] = < psi_i, 1 >`} \item{c1 }{`c1[i,j] = <
#' psi_i, psi_j >`} \item{g1 }{`g1[i,j] = < grad psi_i, grad psi_j >`}
#' \item{g2 }{`g2 = g1 * c0^-1 * g1`} \item{gk }{`gk = g1 * (c0^-1 *
#' g1)^(k-1)`, up to and including `k=order`}
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @seealso [inla.mesh.1d.fem()]
#' @export inla.mesh.fem
inla.mesh.fem <- function(mesh, order = 2) {
  fmesher_deprecate(
    "soft",
    2L,
    "23.08.18",
    "inla.mesh.fem()",
    "fmesher::fm_fem()"
  )
  return(fmesher::fm_fem(mesh, order = order))
}





#' @title Directional derivative matrices for functions on meshes.
#'
#' @description
#' `r lifecycle::badge("deprecated")` Use [fmesher::fm_basis()] instead.
#'
#' Calculates directional derivative matrices for functions on
#' [inla.mesh()] objects.
#'
#' @param mesh An [inla.mesh()] object.
#' @param loc Coordinates where the derivatives should be evaluated.
#' @return \item{A }{The projection matrix, `u(loc_i)=sum_j A_ij w_i`}
#' \item{dx, dy, dz }{Derivative weight matrices, `du/dx(loc_i)=sum_j
#' dx_ij w_i`, etc.}
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @export inla.mesh.deriv
inla.mesh.deriv <- function(mesh, loc) {
  fmesher_deprecate(
    "soft",
    2L,
    "23.08.18",
    "inla.mesh.deriv()",
    I("`fmesher::fm_basis(..., derivatives = TRUE)`")
  )
  return(fmesher::fm_basis(mesh, loc, derivatives = TRUE))
}



#' @title Recursive curve simplification.
#'
#' @description
#' `r lifecycle::badge("deprecated")` Use [fmesher::fm_simplify_helper()]
#' instead.
#'
#' Attempts to simplify a polygonal curve by joining nearly colinear segments.
#'
#' Uses a variation of the binary splitting Ramer-Douglas-Peucker algorithm,
#' with a width `eps` ellipse instead of a rectangle, motivated by
#' prediction ellipse for Brownian bridge.
#'
#' @param loc Coordinate matrix.
#' @param idx Index vector into `loc` specifying a polygonal curve.
#' @param eps Straightness tolerance.
#' @return An index vector into `loc` specifying the simplified polygonal
#' curve.
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @examples
#'
#' theta <- seq(0, 2 * pi, length.out = 1000)
#' loc <- cbind(cos(theta), sin(theta))
#' idx <- inla.simplify.curve(loc = loc, idx = 1:nrow(loc), eps = 0.01)
#' print(c(nrow(loc), length(idx)))
#' plot(loc, type = "l")
#' lines(loc[idx, ], col = "red")
#' @export inla.simplify.curve
inla.simplify.curve <- function(loc, idx, eps) {
  fmesher_deprecate(
    "soft",
    2L,
    "23.08.18",
    "inla.simplify.curve()",
    "fmesher::fm_simplify_helper()"
  )
  return(fmesher::fm_simplify_helper(loc = loc, idx = idx, eps = eps))
}


#' @export
#' @describeIn inla.mesh.segment `r lifecycle::badge("deprecated")` Use
#'   [fmesher::fm_segm_contour_helper()] instead.
inla.contour.segment <- function(x = seq(0, 1, length.out = nrow(z)),
                                 y = seq(0, 1, length.out = ncol(z)),
                                 z,
                                 nlevels = 10,
                                 levels = pretty(range(z, na.rm = TRUE), nlevels),
                                 groups = seq_len(length(levels)),
                                 positive = TRUE,
                                 eps = NULL,
                                 crs = NULL) {
  fmesher_deprecate(
    "soft",
    2L,
    "23.08.18",
    "inla.contour.segment()",
    "fmesher::fm_segm_contour_helper()"
  )
  return(fmesher::fm_segm_contour_helper(
    x = x,
    y = y,
    z = z,
    nlevels = nlevels,
    levels = levels,
    groups = groups,
    positive = positive,
    eps = eps,
    crs = crs
  ))
}




## Morphological dilation by "convex",
## followed by closing by "concave", with
## minimum concave curvature radius "concave".
## If the dilated set has no gaps of width between
## 2*convex*(sqrt(1+2*concave/convex) - 1)
## and
## 2*concave,
## then the minimum convex curvature radius is "convex".
## Default is concave=convex
## Special case concave=0 delegates to inla.nonconvex.hull.basic()
##
## The implementation is based on the identity
## dilation(a) & closing(b) = dilation(a+b) & erosion(b)
## where all operations are with respect to disks with the specified radii.


#' @title Nonconvex set extensions.
#'
#' @description
#' `r lifecycle::badge("deprecated")` Use [fmesher::fm_nonconvex_hull_inla()] or
#' [fmesher::fm_nonconvex_hull()] instead.
#'
#' Constructs a nonconvex boundary for a point set using morphological
#' operations.
#'
#' @details
#' Morphological dilation by `convex`, followed by closing by
#' `concave`, with minimum concave curvature radius `concave`.  If
#' the dilated set has no gaps of width between \deqn{2 convex (\sqrt{1+2
#' concave/convex} - 1)}{2*convex*(sqrt(1+2*concave/convex) - 1)} and \eqn{2
#' concave}{2*concave}, then the minimum convex curvature radius is
#' `convex`.  Special case `concave=0` delegates to
#' `inla.nonconvex.hull.basic`
#'
#' The implementation is based on the identity \deqn{dilation(a) \& closing(b)
#' = dilation(a+b) \& erosion(b)}{dilation(a) & closing(b) = dilation(a+b) &
#' erosion(b)} where all operations are with respect to disks with the
#' specified radii.
#'
#' @aliases inla.nonconvex.hull inla.nonconvex.hull.basic
#' @param points 2D point coordinates (2-column matrix).  Can alternatively be
#' a `SpatialPoints` or `SpatialPointsDataFrame` object.
#' @param convex The desired extension radius.  Also determines the smallest
#' allowed convex curvature radius.  Negative values are interpreted as
#' fractions of the approximate initial set diameter.
#' @param concave The desired minimal concave curvature radius. Default is
#' `concave=convex`.
#' @param resolution The internal computation resolution.  A warning will be
#' issued when this needs to be increased for higher accuracy, with the
#' required resolution stated.
#' @param eps The polygonal curve simplification tolerance used for simplifying
#' the resulting boundary curve.  See [inla.simplify.curve()] for
#' details.
#' @param crs An optional `CRS` or `inla.CRS` object
#' @return An [inla.mesh.segment()] object.
#' @note Requires `nndistF` from the `splancs` package.
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @examples
#'
#' if (require(splancs)) {
#'   loc <- matrix(runif(20), 10, 2)
#'   boundary <- inla.nonconvex.hull(loc, convex = 0.2)
#'   lines(boundary, add = FALSE)
#'   points(loc)
#' }
#' @export inla.nonconvex.hull
inla.nonconvex.hull <- function(points,
                                convex = -0.15,
                                concave = convex,
                                resolution = 40,
                                eps = NULL,
                                crs = NULL) {
  fmesher_deprecate(
    "soft",
    2L,
    "23.08.18",
    "inla.nonconvex.hull()",
    I("`fmesher::fm_nonconvex_hull()` or `fmesher::fm_nonconvex_hull_inla()`")
  )
  return(fmesher::fm_nonconvex_hull_inla(
    x = points,
    convex = convex,
    concave = concave,
    resolution = resolution,
    eps = eps,
    crs = crs
  ))
}

#' @export
#' @describeIn inla.nonconvex.hull
#' `r lifecycle::badge("deprecated")`
#' Use [fmesher::fm_nonconvex_hull_inla_basic()] instead.
## Based on an idea from Elias Teixeira Krainski
## Requires  splancs::nndistF
inla.nonconvex.hull.basic <- function(points,
                                      convex = -0.15,
                                      resolution = 40,
                                      eps = NULL,
                                      crs = NULL) {
  fmesher_deprecate(
    "soft",
    2L,
    "23.08.18",
    "inla.nonconvex.hull.basic()",
    I(paste0(
      "`fmesher::fm_nonconvex_hull()` or ",
      "`fmesher::fm_nonconvex_hull_inla_basic()`"
    ))
  )
  return(fmesher::fm_nonconvex_hull_inla_basic(
    x = points,
    convex = convex,
    resolution = resolution,
    eps = eps,
    crs = crs
  ))
}
