`inla.generate.colors` = function(color, color.axis = NULL, color.n=512,
  color.palette = cm.colors, color.truncate=FALSE, alpha=NULL)
  {
    if (is.character(color)) {
      colors = color
    } else if (is.vector(color) || (is.matrix(color) && (ncol(color)==1))) {
      if (is.null(color.axis))
        color.axis = c(min(color, na.rm=TRUE), max(color, na.rm=TRUE))
      if (color.truncate) {
        not.ok = ((color<color.axis[1]) |
                  (color>color.axis[2]))
      } else {
        not.ok = rep(FALSE, length(color))
      }
      cs = (pmax(color.axis[1],
                 pmin(color.axis[2], color, na.rm=TRUE), na.rm=TRUE))
      cs = (cs-color.axis[1])/(color.axis[2]-color.axis[1])
      not.ok = not.ok | is.na(cs)
      cs[not.ok] = 0.5
      if (is.null(alpha)) {
        alpha = as.numeric(!not.ok)
      } else {
        alpha[not.ok] = 0
      }

      ics = (as.numeric(cut(cs, seq(0, 1, length.out=color.n+1),
                            include.lowest=TRUE)))
      colors = color.palette(color.n)[ics]

      ## Todo: handle alpha, combining "input alpha" with "not.ok-alpha"
    } else if (is.matrix(color) && (ncol(color)==3)) {
      if (is.null(color.axis))
        color.axis = c(min(color, na.rm=TRUE), max(color, na.rm=TRUE))
      if (color.truncate) {
        not.ok = ((color[, 1]<color.axis[1]) |
                  (color[, 2]<color.axis[1]) |
                  (color[, 3]<color.axis[1]) |
                  (color[, 1]>color.axis[2]) |
                  (color[, 2]>color.axis[2]) |
                  (color[, 3]>color.axis[2]))
      } else {
        not.ok = rep(FALSE, nrow(color))
      }
      cs = matrix(
        pmax(color.axis[1],
             pmin(color.axis[2], color, na.rm=TRUE), na.rm=TRUE), dim(color))
      cs = (cs-color.axis[1])/(color.axis[2]-color.axis[1])
      not.ok = not.ok | is.na(cs[, 1]) | is.na(cs[, 2]) | is.na(cs[, 3])
      cs[not.ok,] = c(0.5, 0.5, 0.5)
      if (is.null(alpha)) {
        alpha = as.numeric(!not.ok)
      } else {
        alpha[not.ok] = 0
      }
      colors = rgb(cs[, 1], cs[, 2], cs[, 3])
    } else {
      stop("color specification must be character, matrix, or vector.")
    }

    return (list(colors=colors, alpha=alpha))
  }


`plot.inla.trimesh` = function(TV, S, color = NULL, color.axis = NULL,
  color.n=512, color.palette = cm.colors, color.truncate=FALSE, alpha=NULL,
  lwd = 1, specular = "black",
  draw.vertices=TRUE, draw.edges=TRUE, edge.color=rgb(0.3, 0.3, 0.3), ...)
{
    ## Make indices 1 based.  Deprecated and will be deactivated.
    if (min(TV) == 0) {
      warning("Zero-based indices in TV are deprecated and will be deactivated in a future version.")
        TV <- TV + 1L
    }

    colors = (inla.generate.colors(color, color.axis, color.n,
                                   color.palette, color.truncate, alpha))

    tTV = t(TV);
    tETV = t(TV[, c(1, 2, 3, 1, NA)]);
    Tx = S[tTV, 1]
    Ty = S[tTV, 2]
    Tz = S[tTV, 3]
    if (length(colors$colors) == 1) {
      ## One color
      Tcol = colors$colors
      Talpha = colors$alpha
    } else if (length(colors$colors) == nrow(S)) {
      ## One color per vertex
      Tcol = colors$colors[tTV]
      Talpha = colors$alpha[tTV]
    } else {
      ## One color per triangle
      stopifnot(length(colors$colors) == nrow(TV))
      Tcol = colors$colors[t(matrix(rep(1:nrow(TV), 3), dim(TV)))]
      Talpha = colors$alpha[t(matrix(rep(1:nrow(TV), 3), dim(TV)))]
    }
    Ex = S[tETV, 1]
    Ey = S[tETV, 2]
    Ez = S[tETV, 3]
    Ecol = edge.color
    if (draw.vertices) {
      points3d(S, color="black", ...)
    }
    if (draw.edges) {
      lines3d(Ex, Ey, Ez, color=Ecol, lwd=lwd, ...)
    }
    triangles3d(Tx, Ty, Tz, color=Tcol, specular=specular, alpha=Talpha, ...)

    return (invisible())
}

## library(geometry)
## S = cbind(x=rnorm(30), y=rnorm(30), z=0)
## TV = delaunayn(S[, 1:2]) # NOTE: inconsistent triangle orders, only for test.
## trimesh(TV, S)
##
## colors = rgb(runif(30), runif(30), runif(30))
## rgl.viewpoint(0, 0, fov=20)
## plot.inla.trimesh(TV, S, colors)

##    Ecol = col2rgb(color)/256
##    Ecol = Ecol*0.5+(1-0.5)*0 # Rescale towards black
##    Ecol = 1-Ecol # Invert
##    Ecol = Ecol[, c(2, 3, 1)] # Permute
##    Ecol = rgb(Ecol[1,], Ecol[2,], Ecol[3,], maxColorValue = 1)
##    Ecol = Ecol[tETV]




`plot.inla.fmesher.mesh` = function(m, color = "green", size = 2, lwd=2, add=FALSE, draw.vertices=TRUE, ...)
{
    ## a simple function that plots the mesh from inla.fmesher.mesh()

    if (length(color) == 1)
        color = rep(color, dim(m$mesh$s)[1])

    require(rgl)
    if (!add) {
        dev=open3d()
        view3d(0, 0, fov=0)
    } else {
        dev = NULL
    }
    if (draw.vertices)
        rgl.points(m$mesh$s[m$locations.idx, ],
                   size=2*size, lwd=lwd, color = "blue", ...)
    plot.inla.trimesh(m$mesh$tv, m$mesh$s, color = color,
                      size=size, lwd=lwd,
                      draw.vertices=draw.vertices, add=add, ...)

    return (invisible(dev))
}
