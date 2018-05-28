## Export: inla.contour.segment inla.delaunay inla.fmesher.smorg
## Export: inla.generate.colors inla.mesh inla.mesh.1d inla.mesh.1d.A
## Export: inla.mesh.1d.bary inla.mesh.1d.fem inla.mesh.2d inla.mesh.basis
## Export: inla.mesh.boundary inla.mesh.create inla.mesh.create.helper
## Export: inla.mesh.deriv inla.mesh.fem
## Export: inla.mesh.interior inla.mesh.lattice inla.mesh.map
## Export: inla.mesh.map.lim inla.mesh.query
## Export: inla.nonconvex.hull inla.nonconvex.hull.basic
## Export: inla.simplify.curve plot.inla.trimesh plot!inla.trimesh
## Internal: inla.mesh.filter.locations
## Internal: inla.mesh.parse.segm.input inla.mesh.extract.segments
##
## S3methods; also export some methods explicitly
## Export: extract.groups inla.mesh.project inla.mesh.projector
## Export: extract.groups!inla.mesh.segment
## Export: inla.mesh.segment
## Export: inla.mesh.segment.default
## Export: inla.mesh.segment.inla.mesh.segment
## Export: inla.mesh.segment!default
## Export: inla.mesh.segment!inla.mesh.segment
## Export: inla.mesh.project.inla.mesh inla.mesh.project.inla.mesh.1d
## Export: inla.mesh.project.inla.mesh.projector
## Export: inla.mesh.projector.inla.mesh inla.mesh.projector.inla.mesh.1d
## Export: inla.mesh.project!inla.mesh inla.mesh.project!inla.mesh.1d
## Export: inla.mesh.project!inla.mesh.projector
## Export: inla.mesh.projector!inla.mesh inla.mesh.projector!inla.mesh.1d
## Export: lines.inla.mesh.segment plot.inla.mesh
## Export: print.summary.inla.mesh summary.inla.mesh
## Export: lines!inla.mesh.segment plot!inla.mesh
## Export: print!summary.inla.mesh summary!inla.mesh
## Export: inla.diameter
## Export: inla.diameter!default
## Export: inla.diameter!inla.mesh.1d
## Export: inla.diameter!inla.mesh
## Export: inla.diameter!inla.mesh.segment
## Export: inla.diameter!inla.mesh.lattice



inla.mesh.segment <- function(...) {
    UseMethod("inla.mesh.segment")
}

inla.mesh.segment.default <-
  function(loc = NULL, idx = NULL, grp = NULL, is.bnd = TRUE,
           crs=NULL, ...)
{
    if ((missing(loc) || is.null(loc)) &&
        (missing(idx) || is.null(idx))) {
        stop("At most one of 'loc' and 'idx' may be missing or null.")
    }
    if (!missing(loc) && !is.null(loc)) {
      ## Handle loc given as SpatialPoints or SpatialPointsDataFrame object
      if (inherits(loc, "SpatialPoints") ||
          inherits(loc, "SpatialPointsDataFrame")) {
        loc <- inla.spTransform(coordinates(loc),
                                CRS(proj4string(loc)),
                                crs,
                                passthrough=TRUE)
      }

        if (!is.matrix(loc)) {
            loc = as.matrix(loc)
        }
        if (!is.double(loc)) {
            storage.mode(loc) = "double"
        }
        if (missing(idx) || is.null(idx))
            idx = (inla.ifelse(is.bnd,
                               c(1:nrow(loc),1),
                               c(1:nrow(loc))))
    } else {
        loc = NULL
    }

    if (!missing(idx) && !is.null(idx)) {
        if (!is.vector(idx) && !is.matrix(idx))
            stop("'idx' must be a vector or a matrix")
        if (is.vector(idx))
            idx = as.matrix(idx, nrow=length(idx), ncol=1)
        if (ncol(idx) == 1) {
            if (nrow(idx) < 2) {
                if (nrow(idx) == 1) {
                    warning("Segment specification must have at least 2, or 0, indices.")
                }
                idx <- matrix(0L, 0, 2)
            } else {
                idx = matrix(c(idx[-nrow(idx)],idx[-1]), nrow(idx)-1, 2)
            }
        }
        storage.mode(idx) <- "integer"
        if (!is.null(loc) &&
            (nrow(idx) > 0) &&
            (max(idx, na.rm=TRUE) > nrow(loc))) {
            warning("Segment indices (max=", max(idx, na.rm=TRUE),
                    ") exceed specified location list length (",
                    nrow(loc), ").")
        }
    }

    if (!missing(grp) && !is.null(grp)) {
        if (!is.vector(grp) && !is.matrix(grp))
            stop("'grp' must be a vector or a matrix")
        grp = matrix(grp, min(length(grp), nrow(idx)), 1)
        if (nrow(grp)<nrow(idx))
            grp = (matrix(c(as.vector(grp),
                            rep(grp[nrow(grp)], nrow(idx)-length(grp))),
                            nrow(idx), 1))
        storage.mode(grp) <- "integer"
    } else
        grp = NULL

    ## Filter away NAs in loc and idx
    if (!is.null(loc)) {
        idx[is.na(idx)] = 0L ## Avoid R annoyances with logical+NA indexing
        while (sum(is.na(loc))>0) {
            i = min(which(rowSums(is.na(loc))>0))
            loc = loc[-i,,drop=FALSE]
            idx[idx==i] = 0L
            idx[idx>i] = idx[idx>i]-1L
        }
        idx[idx==0L] = NA
    }
    while (sum(is.na(idx))>0) {
        i = min(which(rowSums(is.na(idx))>0))
        idx = idx[-i,,drop=FALSE]
        if (!is.null(grp))
            grp = grp[-i,,drop=FALSE]
    }

    if (!is.null(loc)) {
        ## Identify unused locations and remap indices accordingly.
        idx.new = rep(0L, nrow(loc))
        idx.new[as.vector(idx)] = 1L
        loc = loc[idx.new==1L,, drop=FALSE]
        idx.new[idx.new==1L] = seq_len(sum(idx.new))
        idx = (matrix(idx.new[as.vector(idx)],
                      nrow=nrow(idx),
                      ncol=ncol(idx)))
    }

    ret = list(loc=loc, idx=idx, grp=grp, is.bnd=is.bnd, crs=crs)
    class(ret) <- "inla.mesh.segment"
    return(invisible(ret))
}

inla.mesh.segment.inla.mesh.segment <- function(..., grp.default=0) {
    segm <- list(...)
    if (!all(unlist(lapply(segm,
                           function(x) inherits(x,"inla.mesh.segment"))))) {
        stop("All objects must be of class 'inla.mesh.segment'.")
    }

    Nloc <- unlist(lapply(segm, function(x) nrow(x$loc)))
    cumNloc <- c(0, cumsum(Nloc))
    Nidx <- unlist(lapply(segm, function(x) nrow(x$idx)))

    loc <- do.call(rbind, lapply(segm, function(x) x$loc))
    idx <- do.call(rbind, lapply(seq_along(segm),
                                function(x) segm[[x]]$idx+cumNloc[x]))
    grp <- unlist(lapply(seq_along(segm),
                         function(x) {
                             if (is.null(segm[[x]]$grp)) {
                                 rep(grp.default, Nidx[x])
                             } else {
                                 segm[[x]]$grp
                             }
                         }))
    is.bnd <- unlist(lapply(segm, function(x) x$is.bnd))
    if (!all(is.bnd) || (any(!is.bnd) && !all(!is.bnd))) {
        warning("Inconsistent 'is.bnd' attributes.  Setting 'is.bnd=FALSE'.")
        is.bnd <- FALSE
    } else {
        is.bnd <- all(is.bnd)
    }

    crs <- lapply(segm, function(x) x$crs)
    if (!is.null(crs)) {
      crs <- crs[unlist(lapply(crs,
                               function(x) !is.null(x)))]
      if (length(crs) > 0) {
        if (!all(unlist(lapply(crs,
                               function(x) identical(crs[[1]], x))))) {
          lapply(crs, function(x) show(x))
          stop("Inconsistent 'crs' attributes.")
        } else {
          crs <- crs[[1]]
        }
      } else {
        crs <- NULL
      }
    }

    invisible(inla.mesh.segment(loc=loc, idx=idx, grp=grp, is.bnd=is.bnd,
                                crs=crs))
}




lines.inla.mesh.segment <- function(x, loc=NULL, col=NULL,
                                    colors=c("black", "blue", "red", "green"),
                                    add=TRUE, xlim=NULL, ylim=NULL,
                                    rgl=FALSE, ...)
{
    segm = x
    if (!is.null(segm$loc))
        loc = segm$loc
    stopifnot(!is.null(loc), ncol(loc)>=2)
    if (ncol(loc) < 3) {
        loc = cbind(loc, 0.0)
    }
    color = col
    if (rgl) {
        stopifnot(inla.require("rgl"))
        if (!add) {
            dev = rgl::open3d()
            rgl::view3d(0, 0, fov=0)
        } else {
            dev = NULL
        }
    } else if (!add) {
        idx = unique(as.vector(segm$idx))
        if (is.null(xlim))
            xlim=range(loc[idx,1])
        if (is.null(ylim))
            ylim=range(loc[idx,2])
        plot.new()
        plot.window(xlim=xlim, ylim=ylim, ...)
    }

    grps = inla.ifelse(is.null(segm$grp), rep(0L,nrow(segm$idx)), segm$grp)
    for (grp in unique(grps)) {
        idx = which(grps==grp)
        if (is.null(col)) {
            color=colors[1+(grp%%length(colors))]
        }
        if (rgl) {
            rgl::segments3d(loc[as.vector(t(segm$idx[idx,, drop=FALSE])),,drop=FALSE],
                            color=color,
                            ...)
        } else {
            lines(loc[t(cbind(segm$idx[idx,, drop=FALSE], NA)), 1],
                  loc[t(cbind(segm$idx[idx,, drop=FALSE], NA)), 2],
                  col=color,
                  ...)
        }
    }
}



`inla.generate.colors` <- function(color,
                                   color.axis = NULL,
                                   color.n=512,
                                   color.palette = cm.colors,
                                   color.truncate=FALSE,
                                   alpha=NULL)
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


## Export as plot!inla.trimesh as well as plot.inla.trimesh,
## even though it's not actually an S3 method.
`plot.inla.trimesh` <- function(x, S, color = NULL, color.axis = NULL,
                                color.n=512, color.palette = cm.colors,
                                color.truncate=FALSE, alpha=NULL,
                                lwd = 1, specular = "black",
                                draw.vertices=TRUE, draw.edges=TRUE,
                                edge.color=rgb(0.3, 0.3, 0.3), ...)
{
    TV = x
    stopifnot(inla.require("rgl"))

    ## Make indices 1 based.  Deprecated and is deactivated.
    if (min(TV) == 0) {
        stop("Zero-based indices in TV are not supported.")
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
        rgl::points3d(S, color="black", ...)
    }
    if (draw.edges) {
        rgl::lines3d(Ex, Ey, Ez, color=Ecol, lwd=lwd, ...)
    }
    rgl::triangles3d(Tx, Ty, Tz, color=Tcol, specular=specular, alpha=Talpha, ...)

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





plot.inla.mesh <- function(x,
                           col="white",
                           t.sub=1:nrow(mesh$graph$tv),
                           add=FALSE,
                           lwd=1,
                           xlim = range(mesh$loc[,1]),
                           ylim = range(mesh$loc[,2]),
                           main = NULL,
                           rgl = FALSE,
                           size = 2,
                           draw.vertices=FALSE,
                           vertex.color="black",
                           draw.edges=TRUE,
                           edge.color=rgb(0.3, 0.3, 0.3),
                           draw.segments=draw.edges,
                           ...)
{
    inla.require.inherits(x, "inla.mesh", "'mesh'")
    mesh = x

    if (rgl) {
        stopifnot(inla.require("rgl"))
        if (!add) {
            dev = rgl::open3d()
            rgl::view3d(0, 0, fov=0)
        } else {
            dev = NULL
        }
        tv = mesh$graph$tv[t.sub,,drop=FALSE]
        if (draw.vertices) {
            idx = intersect(unique(as.vector(tv)), mesh$idx$loc)
            rgl::points3d(mesh$loc[idx,,drop=FALSE],
                          size=2*size, lwd=lwd, color = "blue", ...)
        }
        if (draw.segments) {
            if (!is.null(mesh$segm$bnd))
                lines(mesh$segm$bnd, mesh$loc, lwd=lwd+1,
                      rgl=TRUE, add=TRUE, ...)
            if (!is.null(mesh$segm$int))
                lines(mesh$segm$int, mesh$loc, lwd=lwd+1,
                      rgl=TRUE, add=TRUE, ...)
        }
        plot.inla.trimesh(tv, mesh$loc, color = col,
                          size=size, lwd=lwd,
                          draw.vertices=draw.vertices,
                          vertex.color=vertex.color,
                          draw.edges=draw.edges,
                          edge.color=edge.color,
                          ...)
        return(invisible(dev))
    } else {
        idx = cbind(mesh$graph$tv[t.sub,c(1:3,1), drop=FALSE], NA)
        x = mesh$loc[t(idx), 1]
        y = mesh$loc[t(idx), 2]

        if (!add) {
            plot.new()
            plot.window(xlim=xlim, ylim=ylim, ...)
        }
        if (draw.edges) {
            lines(x, y, type="l", col=edge.color, lwd=lwd)
        }
        tv = mesh$graph$tv[t.sub,,drop=FALSE]
        if (draw.vertices) {
            idx = unique(as.vector(tv))
            points(mesh$loc[idx,,drop=FALSE],
                   pch=20, col=vertex.color, ...)
            idx = intersect(idx, mesh$idx$loc)
            points(mesh$loc[idx,,drop=FALSE],
                   pch=20, col="blue", ...)
        }
        if (draw.segments) {
            if (!is.null(mesh$segm$bnd))
                lines(mesh$segm$bnd, mesh$loc, lwd=lwd+1, ...)
            if (!is.null(mesh$segm$int))
                lines(mesh$segm$int, mesh$loc, lwd=lwd+1, ...)
        }
        if (!add && missing(main)) {
            if (mesh$meta$is.refined) {
                title("Constrained refined Delaunay triangulation", ...)
            } else {
                title("Constrained Delaunay triangulation", ...)
            }
        }
        return(invisible())
    }
}


inla.mesh.map.lim <-
    function(loc=NULL,
             projection=
             c("default", "longlat", "longsinlat", "mollweide"))
{
    projection = match.arg(projection)
    if (identical(projection, "default")) {
        if (is.null(loc)) {
            lim = list(xlim=c(0, 1), ylim=c(0, 1))
        } else {
            lim =
                list(xlim=range(loc[,1], na.rm=TRUE),
                     ylim=range(loc[,2], na.rm=TRUE))
        }
    } else if (identical(projection, "longlat")) {
        lim = list(xlim=c(-180,180), ylim=c(-90,90))
    } else if (identical(projection, "longsinlat")) {
        lim = list(xlim=c(-180,180), ylim=c(-1,1))
    } else if (identical(projection, "mollweide")) {
        lim = list(xlim=c(-2,2), ylim=c(-1,1))
    } else {
        stop(paste("Unknown projection '", projection, "'.", sep=""))
    }
    return(lim)
}

inla.mesh.map <-
    function(loc,
             projection=
             c("default", "longlat", "longsinlat", "mollweide"),
             inverse=TRUE)
{
    projection = match.arg(projection)
    if (identical(projection, "default")) {
        return(loc)
    } else if (identical(projection, "longlat")) {
        if (inverse) {
            proj =
                cbind(cos(loc[,1]*pi/180)*cos(loc[,2]*pi/180),
                      sin(loc[,1]*pi/180)*cos(loc[,2]*pi/180),
                      sin(loc[,2]*pi/180))
        } else {
            proj =
                cbind(atan2(loc[,2], loc[,1])*180/pi,
                      asin(pmax(-1, pmin(+1, loc[,3])))*180/pi)
        }
    } else if (identical(projection, "longsinlat")) {
        if (inverse) {
            coslat = sqrt(pmax(0, 1-loc[,2]^2))
            proj =
                cbind(cos(loc[,1]*pi/180)*coslat,
                      sin(loc[,1]*pi/180)*coslat,
                      loc[,2]
                      )
        } else {
            proj =
                cbind(atan2(loc[,2], loc[,1])*180/pi,
                      loc[,3])
        }
    } else if (identical(projection, "mollweide")) {
        if (inverse) {
            ok = ((loc[,1]^2+4*loc[,2]^2) <= 4)
            cos.theta = sqrt(pmax(0, 1-loc[ok,2]^2))
            theta = atan2(loc[ok,2], cos.theta)
            sin.lat = (2*theta+sin(2*theta))/pi
            cos.lat = sqrt(pmax(0, 1-sin.lat^2))
            lon = loc[ok,1]*pi/2/(cos.theta+(cos.theta==0))
            lon[cos.theta==0] = pi/2*sign(theta[cos.theta==0])
            proj = matrix(NA, nrow(loc), 3)
            proj[ok,] = cbind(cos(lon)*cos.lat, sin(lon)*cos.lat, sin.lat)
        } else {
            lon = atan2(loc[,2],loc[,1])
            z = pmin(1, pmax(-1, loc[,3]))
            sin.theta = z
            cos.theta = sqrt(pmax(0, 1-sin.theta^2))
            ## NR-solver for sin.theta.
            ## Typically finishes after at most 7 iterations.
            ## When cos.theta=0, sin.theta is already correct, +/- 1.
            nook = (cos.theta > 0)
            for (k in 1:20) {
              if (any(nook)) {
                delta =
                  (atan2(sin.theta[nook], cos.theta[nook]) +
                     sin.theta[nook]*cos.theta[nook] - pi/2*z[nook])/
                    (2*cos.theta[nook])
                sin.theta[nook] = sin.theta[nook] - delta
                cos.theta[nook] = sqrt(1-sin.theta[nook]^2)
                nook[nook] = (abs(delta) > 1e-14)
              }
            }
            proj = cbind(2*lon/pi*cos.theta, sin.theta)
        }
    } else {
        stop(paste("Unknown projection '", projection, "'.", sep=""))
    }
    return(proj)
}


inla.mesh.lattice <- function(x=seq(0, 1, length.out=2),
                              y=seq(0, 1, length.out=2),
                              z=NULL,
                              dims =
                              if (is.matrix(x)) {
                                  dim(x)
                              } else {
                                  c(length(x), length(y))
                              },
                              units = NULL,
                              crs = NULL)
{
  if (is.null(crs)) {
    units = match.arg(units, c("default", "longlat", "longsinlat", "mollweide"))

    lim = inla.mesh.map.lim(projection=units)
    xlim = lim$xlim
    ylim = lim$ylim
  } else { ## !is.null(crs)
    if (!is.null(units)) {
      stop("Only one of 'units' and 'crs' can be non-null.")
    }
    bounds <- inla.crs.bounds(crs)
    xlim = bounds$xlim
    ylim = bounds$ylim
  }

    if (missing(x) && !missing(dims)) {
      x = seq(xlim[1], xlim[2], length.out=dims[1])
    }
    if (missing(y) && !missing(dims)) {
      y = seq(ylim[1], ylim[2], length.out=dims[2])
    }
    dims = as.integer(dims)

    if (is.matrix(x)) {
        if (!identical(dims, dim(x)) ||
            !identical(dims, dim(y)) ||
            (is.matrix(z) && !identical(dims, dim(z))))
            stop("The size of matrices 'x', 'y', and 'z' must match 'dims'.")
        loc = cbind(as.vector(x), as.vector(y), as.vector(z))
        x = NULL
        y = NULL
    } else {
        if (!identical(dims[1], length(x)) ||
            !identical(dims[2], length(y)))
            stop(paste("The lengths of vectors 'x' and 'y' (",
                       length(x),",",length(y),
                       ") must match 'dims' (",dims[1],",",dims[2],").",
                       sep=""))
        loc = (cbind(rep(x, times = dims[2]),
                     rep(y, each = dims[1])))
    }
    if (!is.double(loc))
        storage.mode(loc) = "double"

  if (is.null(crs)) {
    loc = inla.mesh.map(loc=loc, projection=units, inverse=TRUE)
  }

    ## Construct lattice boundary
    segm.idx = (c(1:(dims[1]-1),
                  dims[1]*(1:(dims[2]-1)),
                  dims[1]*dims[2]-(0:(dims[1]-2)),
                  dims[1]*((dims[2]-1):1)+1))
    segm.grp = (c(rep(1L, dims[1]-1),
                  rep(2L, dims[2]-1),
                  rep(3L, dims[1]-1),
                  rep(4L, dims[2]-1)))

    segm = (inla.mesh.segment(loc=loc[segm.idx,, drop=FALSE],
                              grp=segm.grp,
                              is.bnd=TRUE,
                              crs=crs))

    lattice = list(dims=dims, x=x, y=y, loc=loc, segm=segm, crs=crs)
    class(lattice) = "inla.mesh.lattice"
    return(lattice)
}


extract.groups <- function(...)
{
    UseMethod("extract.groups")
}

extract.groups.inla.mesh.segment <- function(segm,
                                             groups,
                                             groups.new=groups,
                                             ...)
{
    inla.require.inherits(segm, "inla.mesh.segment", "'segm'")

    if (length(groups.new)==1L) {
        groups.new = rep(groups.new, length(groups))
    }
    if (length(groups.new)!=length(groups)) {
        stop("Length of 'groups.new' (", length(groups.new),
             ") does not match length of 'groups' (",length(groups),")")
    }

    idx = c()
    segm.grp = c()
    for (k in 1:length(groups)) {
        extract.idx = which(segm$grp==groups[k])
        idx = c(idx, extract.idx)
        segm.grp = c(segm.grp, rep(groups.new[k], length(extract.idx)))
    }
    segm.idx = segm$idx[idx,, drop=FALSE]

    return(inla.mesh.segment(loc=segm$loc,
                             idx=segm.idx,
                             grp=segm.grp,
                             segm$is.bnd,
                             crs=segm$crs))
}



inla.mesh.parse.segm.input <- function(boundary=NULL,
                                       interior=NULL,
                                       segm.offset=0L,
                                       loc.offset=0L,
                                       crs=NULL)
{
###########################################
    homogenise.segm.input <- function(x, is.bnd, crs=NULL)
    {
        if (is.matrix(x) || is.vector(x)) { ## Coordinates or indices
            x = (inla.ifelse(is.matrix(x),x,
                             as.matrix(x,nrow=length(x),ncol=1)))
            if (is.integer(x)) { ## Indices
                ret = inla.mesh.segment(NULL, x, NULL, is.bnd, crs=crs)
            } else if (is.numeric(x)) { ## Coordinates
                ret = inla.mesh.segment(x, NULL, NULL, is.bnd, crs=crs)
            } else {
                stop("Segment info matrix must be numeric or integer.")
            }
        } else if (inherits(x, "inla.mesh.segment")) {
            ## Override x$is.bnd:
            ret = inla.mesh.segment(x$loc, x$idx, x$grp, is.bnd, x$crs)
        } else if (inherits(x, c("SpatialPoints", "SpatialPointsDataFrame",
                                 "Line", "Lines",
                                 "SpatialLines", "SpatialLinesDataFrame",
                                 "Polygon", "Polygons",
                                 "SpatialPolygons",
                                 "SpatialPolygonsDataFrame"))) {
            x <- as.inla.mesh.segment(x)
            ## Override x$is.bnd:
            ret <- inla.mesh.segment(x$loc, x$idx, x$grp, is.bnd, x$crs)
        } else if (!is.null(x)) {
            inla.require.inherits(NULL,
                                  c("matrix", "inla.mesh.segment"),
                                  "Segment info")
        } else {
            ret = NULL
        }
        ret
    }
##################################################
    homogenise.segm.grp <- function(input) {
        grp.idx = 0L
        for (k in 1:length(input)) if (!is.null(input[[k]])) {
            inla.require.inherits(input[[k]],
                                  "inla.mesh.segment",
                                  "Segment info list members ")
            if (is.null(input[[k]]$grp)) {
                grp.idx = grp.idx+1L
                input[[k]] = (inla.mesh.segment(input[[k]]$loc,
                                                input[[k]]$idx,
                                                grp.idx,
                                                input[[k]]$is.bnd,
                                                input[[k]]$crs))
            } else {
                grp.idx = max(grp.idx, input[[k]]$grp, na.rm=TRUE)
            }
        }
        return(input)
    }
##################################################
  join.segm.input <- function(segm, segm.offset=0L, loc.offset=0L,
                              crs=NULL)
    {
        loc = NULL
        bnd = list(loc=NULL, idx = matrix(,0,2), grp = matrix(,0,1))
        int = list(loc=NULL, idx = matrix(,0,2), grp = matrix(,0,1))
        storage.mode(bnd$idx) <- "integer"
        storage.mode(bnd$grp) <- "integer"
        storage.mode(int$idx) <- "integer"
        storage.mode(int$grp) <- "integer"
        for (k in seq_along(segm)) {
          if (!is.null(segm[[k]])) {
            inla.require.inherits(segm[[k]],
                                  "inla.mesh.segment",
                                  "Segment info list members ")
            if (!is.null(segm[[k]]$loc)) {
              local.loc <- inla.spTransform(segm[[k]], crs,
                                            passthrough=TRUE)$loc
              extra.loc.n = nrow(local.loc)
              idx.offset = segm.offset
              segm.offset = segm.offset + extra.loc.n
              loc = (inla.ifelse(is.null(loc),
                                 local.loc,
                                 rbind(loc, local.loc)))
            } else {
              idx.offset = loc.offset
            }
            if (segm[[k]]$is.bnd) {
                bnd$idx = rbind(bnd$idx, segm[[k]]$idx+idx.offset)
                bnd$grp = rbind(bnd$grp, segm[[k]]$grp)
            } else {
                int$idx = rbind(int$idx, segm[[k]]$idx+idx.offset)
                int$grp = rbind(int$grp, segm[[k]]$grp)
            }
          }
        }
        if (nrow(bnd$idx)==0)
            bnd = NULL
        if (nrow(int$idx)==0)
            int = NULL
        return(list(loc = loc,
                    bnd = (inla.ifelse(is.null(bnd),
                                       NULL,
                                       inla.mesh.segment(bnd$loc,
                                                         bnd$idx,
                                                         bnd$grp,
                                                         TRUE,
                                                         crs))),
                    int = (inla.ifelse(is.null(int),
                                       NULL,
                                       inla.mesh.segment(int$loc,
                                                         int$idx,
                                                         int$grp,
                                                         FALSE,
                                                         crs)))))
    }
###########################

    segm = (c(lapply(inla.ifelse(inherits(boundary, "list"),
                                 boundary, list(boundary)),
                     function(x){homogenise.segm.input(x, TRUE, crs=crs)}),
              lapply(inla.ifelse(inherits(interior, "list"),
                                 interior, list(interior)),
                     function(x){homogenise.segm.input(x, FALSE, crs=crs)})))
    segm = homogenise.segm.grp(segm)

  join.segm.input(segm=segm, segm.offset=segm.offset, loc.offset=loc.offset,
                  crs=crs)
}








################################
##
## Old code.  Filtering is now done in fmesher itself.
## Retained for now so that we can check if the results are the same.
##
inla.mesh.filter.locations <- function(loc, cutoff)
    {
        ## Map locations to nodes, avoiding near-duplicates.
        loc.n = nrow(loc)
        loc.dim = ncol(loc)

        loc.is.na = (rowSums(is.na(loc))>0)
        if (sum(loc.is.na)>0)
            stop("NAs in locations not yet supported.")

        node.coord = matrix(nrow=loc.n, ncol=loc.dim)
        map.loc.to.node = rep(0L, nrow=loc.n)
        excluded = c()
        loc.i = 1L
        node.i.max = 1L
        node.coord[node.i.max,] = loc[loc.i,]
        map.loc.to.node[[loc.i]] = node.i.max
        for (loc.i in 2:loc.n) {
            loc.to.node.dist =
                sqrt(rowSums((as.matrix(rep(1, node.i.max)) %*%
                              loc[loc.i,, drop=FALSE] -
                              node.coord[1:node.i.max,, drop=FALSE])^2))
            if (min(loc.to.node.dist) > cutoff) {
                node.i.max = node.i.max+1L
                node.coord[node.i.max,] = loc[loc.i,, drop=FALSE]
                map.loc.to.node[[loc.i]] = node.i.max
            } else {
                excluded = c(excluded, loc.i)
            }
        }
        ## Remove excess nodes.
        node.coord = node.coord[1:node.i.max,]

        ## Identify nearest nodes for excluded locations.
        for (loc.i in excluded) {
            loc.to.node.dist =
                sqrt(rowSums((as.matrix(rep(1, node.i.max)) %*%
                              loc[loc.i,, drop=FALSE] -
                              node.coord)^2))
            node.i = which.min(loc.to.node.dist);
            map.loc.to.node[[loc.i]] = node.i
        }

        return(list(loc = node.coord, node.idx = map.loc.to.node))
    }
############################







inla.mesh <- function(...)
{
    args = list(...)
    if (length(args)>0) {
        if (inherits(args[[1]],"inla.mesh")) {
            warning("'inla.mesh(mesh, ...)' is deprecated.  Use 'inla.mesh.query(mesh, ...)' instead.")
            return(inla.mesh.query(...))
        }
    }
    warning("'inla.mesh(...)' is deprecated.  Use 'inla.mesh.create(...)' instead.")
    return(inla.mesh.create(...))
}


inla.mesh.create <- function(loc=NULL, tv=NULL,
                             boundary=NULL, interior=NULL,
                             extend = (missing(tv) || is.null(tv)),
                             refine=FALSE,
                             lattice=NULL,
                             globe=NULL,
                             cutoff = 1e-12,
                             plot.delay = NULL,
                             data.dir,
                             keep = (!missing(data.dir) && !is.null(data.dir)),
                             timings = FALSE,
                             quality.spec=NULL,
                             crs=NULL)
{
    if (!timings) {
        system.time <- function(expr) {
            expr
            structure(c(0,0,0,0,0), class = "proc_time")
        }
    }

    time.pre = system.time({ ## Pre-processing timing start

    if (!is.null(crs) &&
        identical(inla.as.list.CRS(crs)$proj, "geocent")) {
      ## Build all geocentric meshes on a sphere, and transform afterwards,
      ## to allow general geoids.
      crs.target <- crs
      crs <- inla.CRS("sphere")
    }

    if (!(missing(loc) || is.null(loc))) {
      ## Handle loc given as SpatialPoints or SpatialPointsDataFrame object
      if (inherits(loc, "SpatialPoints") ||
          inherits(loc, "SpatialPointsDataFrame")) {
        loc <- inla.spTransform(coordinates(loc),
                                CRS(proj4string(loc)),
                                crs,
                                passthrough=TRUE)
      }

      if (!is.matrix(loc)) {
        loc = as.matrix(loc)
      }
      if (!is.double(loc)) {
        storage.mode(loc) = "double"
      }
    }
    loc.n = max(0L,nrow(loc))

    if (is.logical(extend) && extend) extend = list()
    if (is.logical(refine) && refine) refine = list()

    if (missing(lattice) || is.null(lattice)) {
        lattice = list(loc=NULL, segm=NULL)
        lattice.n = 0L
    } else {
        inla.require.inherits(lattice, "inla.mesh.lattice", "'lattice'")
        if (!is.null(tv)) {
            warning("Both 'lattice' and 'tv' specified.  Ignoring 'tv'.")
            tv = NULL
        }
        if (!inherits(extend, "list")) {
            boundary = (c(inla.ifelse(inherits(boundary, "list"),
                                      boundary, list(boundary)),
                          list(lattice$segm)))
        }

        if (!is.null(lattice[["crs"]])) {
          lattice$loc <- inla.spTransform(lattice$loc,
                                          lattice$crs,
                                          crs,
                                          passthrough=TRUE)
        }

        if (!is.matrix(lattice$loc)) {
          lattice$loc = as.matrix(lattice$loc)
        }
        if (!is.double(lattice$loc)) {
          storage.mode(lattice$loc) = "double"
        }
    }
      lattice.n = max(0L,nrow(lattice$loc))

    segm = (inla.mesh.parse.segm.input(boundary,
                                       interior,
                                       loc.n,
                                       0L,
                                       crs=crs))
    segm.n = max(0,nrow(segm$loc))
    ## Run parse again now that we know where the indices should point:
    segm = (inla.mesh.parse.segm.input(boundary,
                                       interior,
                                       0L,
                                       segm.n+lattice.n,
                                       crs=crs))

    loc0 = rbind(segm$loc, lattice$loc, loc)
    if ((!is.null(loc0)) && (nrow(loc0)>0))
        idx0 = 1:nrow(loc0)
    else
        idx0 = c()

    if (is.null(quality.spec)) {
        quality <- NULL
    } else {
        quality <- rep(NA, segm.n+lattice.n+loc.n)
        ## Order must be same as for loc0 above.
        if (!is.null(quality.spec$segm)) {
            quality[seq_len(segm.n)] <- quality.spec$segm
        }
        if (!is.null(quality.spec$lattice)) {
            quality[segm.n+seq_len(lattice.n)] <- quality.spec$lattice
        }
        if (!is.null(quality.spec$loc)) {
            quality[segm.n+lattice.n+seq_len(loc.n)] <- quality.spec$loc
        }
        ## NA:s will be replaced with max.edge settings below.
    }

    ## Where to put the files?
    if (keep) {
        if (missing(data.dir)) {
            data.dir = inla.fmesher.make.dir("inla.mesh.data")
        }
        prefix = paste(data.dir, "/mesh.", sep="")
        keep.dir = TRUE
    } else {
        if (missing(data.dir)) {
            prefix = inla.tempfile(pattern="fmesher", tmpdir=tempdir())
            prefix = paste(prefix, ".", sep="")
            keep.dir = TRUE
        } else {
            data.dir = inla.fmesher.make.dir(data.dir)
            prefix = paste(data.dir, "/mesh.", sep="")
            ## We should not try to delete the session tempdir...
            keep.dir = identical(tempdir(), dirname(prefix))
        }
    }
    prefix = inla.fmesher.make.prefix(NULL, prefix)

    if (is.null(loc0)) {
        all.args = ""
    } else {
        loc.file = fmesher.write(inla.affirm.double(loc0), prefix, "input.s")
        all.args = "--input=input.s"

        if (!is.null(tv)) {
            fmesher.write(inla.affirm.integer(tv)-1L, prefix, "input.tv")
            all.args = paste(all.args, ",input.tv", sep="")
        }
    }
    if (!missing(cutoff)) {
        all.args = paste(all.args, " --cutoff=", cutoff, sep="")
    }
    if (!is.null(segm$bnd)) {
        fmesher.write(inla.affirm.integer(segm$bnd$idx)-1L, prefix, "input.segm.bnd.idx")
        fmesher.write(inla.affirm.integer(segm$bnd$grp), prefix, "input.segm.bnd.grp")
        all.args = paste(all.args," --boundary=input.segm.bnd.idx")
        all.args = paste(all.args," --boundarygrp=input.segm.bnd.grp")
    }
    if (!is.null(segm$int)) {
        fmesher.write(inla.affirm.integer(segm$int$idx)-1L, prefix, "input.segm.int.idx")
        fmesher.write(inla.affirm.integer(segm$int$grp), prefix, "input.segm.int.grp")
        all.args = paste(all.args," --interior=input.segm.int.idx")
        all.args = paste(all.args," --interiorgrp=input.segm.int.grp")
    }
    if (!missing(globe) && !is.null(globe)) {
        all.args = paste(all.args, " --globe=", globe, sep="")
    }


    if (inherits(extend,"list")) {
        cet = c(0,0)
        cet[1] = inla.ifelse(is.null(extend$n), 16, extend$n)
        cet[2] = inla.ifelse(is.null(extend$offset), -0.1, extend$offset)
        all.args = (paste(all.args," --cet=",
                          cet[1],",", cet[2], sep=""))
    }
    if (inherits(refine,"list")) {
        rcdt = c(0,0,0)
        if (!missing(globe) && !is.null(globe)) {
          max.edge.default <- pi
        } else {
          ## Multiply by to ensure cover S2 domains
          max.edge.default <- inla.diameter(loc0) * 2
        }
        if ((inherits(extend,"list")) && (!is.null(extend$offset))) {
            max.edge.default = (max.edge.default +
                                max(0,2*extend$offset))
            max.edge.default = (max.edge.default *
                                (1+max(0,-2*extend$offset)))
        }
        rcdt[1] = inla.ifelse(is.null(refine$min.angle), 21, refine$min.angle)
        rcdt[2] = (inla.ifelse(is.null(refine$max.edge) ||
                               is.na(refine$max.edge),
                               max.edge.default,
                               refine$max.edge))
        rcdt[3] = (inla.ifelse(is.null(refine$max.edge) ||
                               is.na(refine$max.edge),
                               max.edge.default,
                               refine$max.edge))
        rcdt[2] = (inla.ifelse(is.null(refine$max.edge.extra) ||
                               is.na(refine$max.edge.extra),
                               rcdt[2], refine$max.edge.extra))
        rcdt[3] = (inla.ifelse(is.null(refine$max.edge.data) ||
                               is.na(refine$max.edge.date),
                               rcdt[3], refine$max.edge.data))
        all.args = (paste(all.args," --rcdt=",
                          rcdt[1],",", rcdt[2],",", rcdt[3], sep=""))

        if (!is.null(refine[["max.n.strict"]]) &&
            !is.na(refine$max.n.strict)) {
          rcdt_max_n0 <- refine$max.n.strict
        } else {
          rcdt_max_n0 <- -1
        }
        if (!is.null(refine[["max.n"]]) &&
            !is.na(refine$max.n)) {
          rcdt_max_n1 <- refine$max.n
        } else {
          rcdt_max_n1 <- -1
        }
        all.args = (paste(all.args,
                          " --max_n0=", rcdt_max_n0,
                          " --max_n1=", rcdt_max_n1,
                          sep=""))
        is.refined = TRUE

        if (!is.null(quality)) {
            quality[is.na(quality)] <- rcdt[3]
            quality.file <-
                fmesher.write(inla.affirm.double(as.matrix(quality)),
                              prefix, "input.quality")
            all.args <- (paste(all.args,
                               " --input=input.quality",
                               " --quality=input.quality",
                               sep=""))
        }
    } else {
        is.refined = FALSE
    }

    if (!is.null(plot.delay)) {
        all.args = paste(all.args," --x11=", plot.delay, sep="")
    }
    all.args = paste(all.args, inla.getOption("fmesher.arg"))

    }) ## Pre-processing timing end

    ## Call fmesher:
    time.fmesher = system.time({
        echoc = (inla.fmesher.call(all.args=all.args,
                                   prefix=prefix))
    })

    time.post = system.time({ ## Post-processing timing start

    ## Read the mesh:
    manifold = 1L+fmesher.read(prefix, "manifold")
    manifold = list("M", "R2", "S2")[[manifold]]

    loc = fmesher.read(prefix, "s")
    graph = (list(tv = 1L+fmesher.read(prefix, "tv"),
                  vt = 1L+fmesher.read(prefix, "vt"),
                  tt = 1L+fmesher.read(prefix, "tt"),
                  tti = 1L+fmesher.read(prefix, "tti"),
                  vv = fmesher.read(prefix, "vv")))
    graph$tv[graph$tv==0L] = NA
    graph$vt[graph$vt==0L] = NA
    graph$tt[graph$tt==0L] = NA
    graph$tti[graph$tti==0L] = NA

    ## Read the vertex input/output mapping:
    idx.all = 1L+fmesher.read(prefix, "idx")
    idx = (list(loc = (inla.ifelse(loc.n>0,
                                   idx.all[idx0[segm.n+lattice.n+(1:loc.n)]],
                                   NULL)),
                lattice = (inla.ifelse(lattice.n>0,
                                    idx.all[idx0[segm.n+(1:lattice.n)]],
                                    NULL)),
                segm = (inla.ifelse(segm.n>0,
                                    idx.all[idx0[(1:segm.n)]],
                                    NULL))))
    if (!is.null(idx$loc)) idx$loc[idx$loc == 0L] = NA
    if (!is.null(idx$lattice)) idx$lattice[idx$lattice == 0L] = NA
    if (!is.null(idx$segm)) idx$segm[idx$segm == 0L] = NA

    ## Read constraint segment information:
    segm.bnd = (inla.mesh.segment(NULL,
                                  1L+fmesher.read(prefix, "segm.bnd.idx"),
                                  fmesher.read(prefix, "segm.bnd.grp"),
                                  TRUE))
    segm.int = (inla.mesh.segment(NULL,
                                  1L+fmesher.read(prefix, "segm.int.idx"),
                                  fmesher.read(prefix, "segm.int.grp"),
                                  FALSE))

    ## Remap indices to remove unused vertices
    used = !is.na(graph$vt)
    if (!all(used)) {
        used = which(used)
        idx.map = rep(NA, nrow(loc))
        idx.map[used] = seq_len(length(used))
        loc = loc[used,,drop=FALSE]
        graph$tv = matrix(idx.map[as.vector(graph$tv)], nrow(graph$tv), 3)
        graph$vt = graph$vt[used,,drop=FALSE]
        ## graph$tt  ## No change needed
        ## graph$tti ## No change needed
        graph$vv = graph$vv[used, used, drop=FALSE]
        if (!is.null(idx$loc)) idx$loc = idx.map[idx$loc]
        if (!is.null(idx$lattice)) idx$lattice = idx.map[idx$lattice]
        if (!is.null(idx$segm)) idx$segm = idx.map[idx$segm]
        segm.bnd$idx = matrix(idx.map[segm.bnd$idx], nrow(segm.bnd$idx), 2)
        segm.int$idx = matrix(idx.map[segm.int$idx], nrow(segm.int$idx), 2)
        segm.bnd$idx[segm.bnd$idx == 0L] = NA
        segm.int$idx[segm.int$idx == 0L] = NA
    }


    if (!keep)
      unlink(paste(prefix, "*", sep=""), recursive=FALSE)
    if (!keep.dir) {
        unlink(dirname(prefix), recursive=TRUE)
    }

    if (!is.null(crs) &&
        identical(inla.as.list.CRS(crs)$proj, "geocent")) {
      if (!inla.identical.CRS(crs, crs.target)) {
        ## Target is a non-spherical geoid
        loc <- inla.spTransform(loc, crs, crs.target)
        crs <- crs.target
      }
    }

    }) ## Post-processing timing end

    time.object = system.time({ ## Object construction timing start
    mesh = (list(meta = (list(call=match.call(),
                              fmesher.args = all.args,
                              time = (rbind(pre = time.pre,
                                            fmesher = time.fmesher,
                                            post = time.post)),
                              prefix = prefix,
                              is.refined = is.refined)),
                 manifold = manifold,
                 n = nrow(loc),
                 loc = loc,
                 graph = graph,
                 segm = list(bnd=segm.bnd, int=segm.int),
                 idx = idx,
                 crs = crs))
    class(mesh) <- "inla.mesh"

    }) ## Object construction timing end

    ## Entire function timing
    time.total = time.pre+time.fmesher+time.post+time.object

    mesh$meta$time = (rbind(mesh$meta$time,
                            object=time.object,
                            total=time.total))
    return(invisible(mesh))
}











inla.mesh.extract.segments <- function(mesh.loc,
                                       mesh.idx,
                                       mesh.grp,
                                       grp=NULL,
                                       is.bnd,
                                       crs=NULL)
{
    segments = list()
    if (nrow(mesh.idx)>0) {
        if (is.null(grp)) {
            grp = unique(sort(mesh.grp))
        }
        for (g in grp) {
            extract = (mesh.grp==g)
            segments =
                c(segments,
                  list(inla.mesh.segment(mesh.loc,
                                         idx=mesh.idx[extract,,drop=FALSE],
                                         grp=mesh.grp[extract,drop=FALSE],
                                         is.bnd=is.bnd,
                                         crs=crs)) )
        }
    }
    if (length(segments)>0)
        return(segments)
    else
        return(NULL)
}

inla.mesh.boundary <- function(mesh, grp=NULL)
{
    inla.require.inherits(mesh, "inla.mesh", "'mesh'")

    return(inla.mesh.extract.segments(mesh$loc,
                                      mesh$segm$bnd$idx,
                                      mesh$segm$bnd$grp,
                                      grp,
                                      TRUE,
                                      mesh$crs))
}

inla.mesh.interior <- function(mesh, grp=NULL)
{
    inla.require.inherits(mesh, "inla.mesh", "'mesh'")

    return(inla.mesh.extract.segments(mesh$loc,
                                      mesh$segm$int$idx,
                                      mesh$segm$int$grp,
                                      grp,
                                      FALSE,
                                      mesh$crs))
}


## Generate nice triangulation, with an inner domain strictly enclosed
## by an outer domain
## The inner and outer domains can have different quality parameters
## At least one of loc, loc.domain, boundary[[1]], boundary[[2]], interior
## must be non-NULL
## For more complicated multi-step meshings, study the code and write your own.
inla.mesh.2d <-
    function(loc=NULL, ## Points to include in final triangulation
             loc.domain=NULL, ## Points that determine the automatic domain
             offset=NULL, ## Size of automatic extensions
             n=NULL, ## Sides of automatic extension polygons
             boundary=NULL, ## User-specified domains (list of length 2)
             interior=NULL, ## User-specified constraints for the inner domain
             max.edge=NULL,
             min.angle=NULL, ## Angle constraint for the entire domain
             cutoff=1e-12, ## Only add input points further apart than this
             max.n.strict=NULL,
             max.n=NULL,
             plot.delay=NULL,
             crs=NULL) ## Coordinate Reference System
    ## plot.delay: Do plotting.
    ## NULL --> No plotting
    ## <0  --> Intermediate meshes displayed at the end
    ## >0   --> Dynamical fmesher plotting
{
###########################
  unify.one.segm <- function(segm, crs=NULL) {
    if (inherits(segm, "inla.mesh.segment")) {
      segm <- inla.spTransform(segm, crs, passthrough=TRUE)
    } else if (inherits(segm, "matrix")) {
      segm <- inla.mesh.segment(loc=segm, crs=crs)
    } else {
      segm <- inla.spTransform(as.inla.mesh.segment(segm), crs,
                               passthrough=TRUE)
    }
    if (ncol(segm$loc)==2) {
      segm$loc = cbind(segm$loc, 0.0)
    }
    segm
  }
  unify.segm.coords <- function(segm, crs=NULL) {
    if (!is.null(segm)) {
      if (inherits(segm, "list")) {
        for (j in seq_along(segm)) {
          segm[[j]] <- unify.one.segm(segm[[j]], crs=crs)
        }
      } else {
        segm <- unify.one.segm(segm, crs=crs)
      }
    }
    segm
  }
###########################

  if ((missing(max.edge) || is.null(max.edge)) &&
      (missing(max.n.strict) || is.null(max.n.strict)) &&
      (missing(max.n) || is.null(max.n))) {
    stop("At least one of max.edge, max.n.strict, and max.n must be specified")
  }

  if (!is.null(crs)) {
    issphere <- inla.identical.CRS(crs, inla.CRS("sphere"))
    isgeocentric <- identical(inla.as.list.CRS(crs)[["proj"]], "geocent")
    if (isgeocentric) {
      crs.target <- crs
      crs <- inla.CRS("sphere")
    }
  }

  ## Handle loc given as SpatialPoints or SpatialPointsDataFrame object
  if (!(missing(loc) || is.null(loc)) &&
      (inherits(loc, "SpatialPoints") ||
       inherits(loc, "SpatialPointsDataFrame"))) {
    loc = inla.spTransform(coordinates(loc),
                           CRS(proj4string(loc)),
                           crs,
                           passthrough=TRUE)
  }
  if (!(missing(loc.domain) || is.null(loc.domain)) &&
      (inherits(loc.domain, "SpatialPoints") ||
       inherits(loc.domain, "SpatialPointsDataFrame"))) {
    loc.domain = inla.spTransform(coordinates(loc.domain),
                                  CRS(proj4string(loc.domain)),
                                  crs,
                                  passthrough=TRUE)
  }

    if (missing(loc) || is.null(loc)) {
        loc = matrix(c(0.0), 0, 3)
    } else if (!is.matrix(loc)) {
        loc = as.matrix(loc)
    }
    if (missing(loc.domain) || is.null(loc.domain)) {
        loc.domain = loc
    } else if (!is.matrix(loc.domain)) {
        loc.domain = as.matrix(loc.domain)
    }
    if (missing(boundary)) {
        boundary = list(NULL)
    } else {
        if (!inherits(boundary, "list"))
            boundary = list(boundary)
    }
    if (missing(interior))
        interior = NULL
    if (missing(offset) || is.null(offset)) {
        if (length(boundary) < 2)
            offset = -0.05
        else
            offset = c(-0.05, -0.15)
    }
    if (missing(n) || is.null(n))
        n = c(8)
    if (missing(max.edge) || is.null(max.edge))
        max.edge = c(NA)
    if (missing(min.angle) || is.null(min.angle))
        min.angle = c(21)
    if (missing(max.n.strict) || is.null(max.n.strict))
        max.n.strict = c(NA)
    if (missing(max.n) || is.null(max.n))
        max.n = c(NA)
    if (missing(cutoff) || is.null(cutoff))
        cutoff = 1e-12
    if (missing(plot.delay) || is.null(plot.delay))
        plot.delay = NULL

    num.layers =
        max(c(length(boundary), length(offset), length(n),
              length(min.angle), length(max.edge),
              length(max.n.strict), length(max.n)))
    if (num.layers > 2) {
        warning(paste("num.layers=", num.layers, " > 2 detected.  ",
                      "Excess information ignored.", sep=""))
        num.layers = 2
    }

    if (length(boundary) < num.layers)
        boundary = c(boundary, list(NULL))
    if (length(min.angle) < num.layers)
        min.angle = c(min.angle, min.angle)
    if (length(max.n.strict) < num.layers)
        max.n0 = c(max.n.strict, max.n.strict)
    if (length(max.n) < num.layers)
        max.n = c(max.n, max.n)
    if (length(max.edge) < num.layers)
        max.edge = c(max.edge, max.edge)
    if (length(offset) < num.layers)
        offset = c(offset, -0.15)
    if (length(n) < num.layers)
        n = c(n, 16)
    if (length(n) < num.layers)
        n = c(n, 16)

    ## Unify the dimensionality of the point input.
    if (!is.null(loc) && (ncol(loc)==2))
        loc = cbind(loc, 0.0)
    if (!is.null(loc.domain) && (ncol(loc.domain)==2))
        loc.domain = cbind(loc.domain, 0.0)
    ## Unify the dimensionality of the boundary&interior segments input
    ## and optionally transform coordinates.
    for (k in seq_len(num.layers)) {
      if (!is.null(boundary[[k]])) {
        boundary[[k]] <- unify.segm.coords(boundary[[k]], crs=crs)
      }
    }
    interior <- unify.segm.coords(interior, crs=crs)

    ## Triangulate to get inner domain boundary
    ## Constraints included only to get proper domain extent
    ## First, attach the loc points to the domain definition set
    if (!is.null(loc) && !is.null(loc.domain)) {
        loc.domain = rbind(loc.domain, loc)
    }
    mesh1 =
        inla.mesh.create(loc=loc.domain,
                         boundary=boundary[[1]],
                         interior=interior,
                         cutoff=cutoff,
                         extend=list(n=n[1], offset=offset[1]),
                         refine=FALSE,
                         plot.delay=plot.delay,
                         crs=crs)

    ## Save the resulting boundary
    boundary1 = inla.mesh.boundary(mesh1)
    interior1 = inla.mesh.interior(mesh1)

    if (!is.null(plot.delay) && (plot.delay<0)) {
        inla.dev.new()
        plot(mesh1)
    }

    ## Triangulate inner domain
    mesh2 =
        inla.mesh.create(loc=loc,
                         boundary=boundary1,
                         interior=interior1,
                         cutoff=cutoff,
                         extend=FALSE, ## Should have no effect
                         refine=
                             list(min.angle=min.angle[1],
                                  max.edge=max.edge[1],
                                  max.edge.extra=max.edge[1],
                                  max.n.strict=max.n.strict[1],
                                  max.n=max.n[1]),
                         plot.delay=plot.delay,
                         crs=crs)

    boundary2 = inla.mesh.boundary(mesh2)
    interior2 = inla.mesh.interior(mesh2)

    if (!is.null(plot.delay) && (plot.delay<0)) {
        inla.dev.new()
        plot(mesh2)
    }

    if (num.layers == 1) {
      if (!is.null(crs) && isgeocentric && !issphere) {
        mesh2$loc <- inla.spTransform(mesh2$loc, mesh2$crs, crs.target)
        mesh2$crs <- crs.target
      }
      return(invisible(mesh2))
    }

    ## Triangulate inner+outer domain
    mesh3 =
        inla.mesh.create(loc=rbind(loc, mesh2$loc),
                         boundary=boundary[[2]],
                         interior=c(boundary2, interior2),
                         cutoff=cutoff,
                         extend=list(n=n[2], offset=offset[2]),
                         refine=
                             list(min.angle=min.angle[2],
                                  max.edge=max.edge[2],
                                  max.edge.extra=max.edge[2],
                                  max.n.strict=mesh2$n + max.n.strict[2],
                                  max.n=mesh2$n + max.n[2]),
                         plot.delay=plot.delay,
                         crs=crs)

    ## Hide generated points, to match regular inla.mesh.create output
    mesh3$idx$loc = mesh3$idx$loc[seq_len(nrow(loc))]

    ## Obtain the corresponding segm indices.
    segm.loc = matrix(0.0, 0, 3)
    for (k in seq_along(boundary)) {
        if (!is.null(boundary[[k]])) {
            segm.loc = rbind(segm.loc, boundary[[k]]$loc)
        }
    }
    for (k in seq_along(interior)) {
        if (!is.null(interior[[k]])) {
            segm.loc = rbind(segm.loc, interior[[k]]$loc)
        }
    }
    if (nrow(segm.loc) > 0) {
        proj = inla.mesh.project(mesh3, loc=segm.loc)
        mesh3$idx$segm = rep(NA, nrow(segm.loc))
        if (any(proj$ok)) {
            t.idx <- proj$t[proj$ok]
            tv.idx <- max.col(proj$bary[proj$ok,,drop=FALSE],
                              ties.method="first")
            mesh3$idx$segm[proj$ok] <-
                mesh3$graph$tv[t.idx + nrow(mesh3$graph$tv)*(tv.idx-1)]
        }
    } else {
        mesh3$idx$segm = NULL
    }

    if (!is.null(crs) && isgeocentric && !issphere) {
      mesh3$loc <- inla.spTransform(mesh3$loc, mesh3$crs, crs.target)
      mesh3$crs <- crs.target
    }

    if (!is.null(plot.delay) && (plot.delay<0)) {
        inla.dev.new()
        plot(mesh3)
    }

    return(invisible(mesh3))
}

## Support for legacy code:
inla.mesh.create.helper <- function(points=NULL, points.domain=NULL, ...)
{
    return(invisible(inla.mesh.2d(loc=points, loc.domain=points.domain, ...)))
}

## Example:
##if (FALSE) {
##    mesh =
##        inla.mesh.create.helper(loc=matrix(runif(20),10,2)*200,
##                                n=c(8,16),
##                                offset=c(10,140),
##                                max.edge=c(25,1000),
##                                min.angle=26,
##                                cutoff=0,
##                                plot.delay=-1
##                                )
##}
##


inla.delaunay <- function(loc, ...)
{
    ## Handle loc given as SpatialPoints or SpatialPointsDataFrame object
    if (!(missing(loc) || is.null(loc)) &&
        (inherits(loc, "SpatialPoints") ||
         inherits(loc, "SpatialPointsDataFrame"))) {
      crs = CRS(proj4string(loc))
      loc = coordinates(loc)
    } else {
      crs = NULL
    }

    hull = chull(loc[,1],loc[,2])
    bnd = inla.mesh.segment(loc=loc[hull[length(hull):1],],is.bnd=TRUE)
    mesh =
        inla.mesh.create(loc=loc,
                         boundary=bnd,
                         extend=list(n=3),
                         refine=FALSE,
                         crs=crs,
                         ...)

    return(invisible(mesh))
}






inla.mesh.query <- function(mesh, ...)
{
    inla.require.inherits(mesh, "inla.mesh", "'mesh'")

    not.known <- function(mesh, queryname)
    {
        stop(paste("Query '", queryname,
                   "' unknown.", sep=""))
    }
    not.implemented <- function(mesh, queryname)
    {
        stop(paste("Query '", queryname,
                   "' not implemented for inla.mesh.", sep=""))
##        stop(paste("Query '", queryname,
##                   "' not implemented for inla.mesh for mesh type '",
##                   model$type, "'.", sep=""))
    }

    result = list()
    queries = inla.parse.queries(...)
    if (length(queries)==0L)
        return(result)

    for (query.idx in 1:length(queries)) {
        query = names(queries)[query.idx]
        param = queries[[query.idx]]
        answer = NULL
        query = (match.arg(query, c("tt.neighbours",
                                    "vt.neighbours"
                                    )))
        if (identical(query, "tt.neighbours")) {
            if (is.null(param))
                param = list(c(1))
            if (length(param)<2)
                param = c(as.list(param), list(c(1)))
            nT = nrow(mesh$graph$tv)
            i = rep(1:nT,3)
            j = as.vector(mesh$graph$tt)
            i = i[!is.na(j)]
            j = j[!is.na(j)]
            tt = sparseMatrix(i=i, j=j, x=rep(1,length(i)), dims=c(nT,nT))
            answer = (sparseMatrix(i=param[[1]],
                                   j=rep(1,length(param[[1]])),
                                   x=1,
                                   dims=c(nT,1)))
            tt0 = (answer == 0.5)*1
            tt1 = answer
            order = 0
            while (order<min(param[[2]])) {
                order = order+1
                tt0 = tt1
                tt1 = answer
                answer = ((((tt %*% answer) > 0) - tt1 - tt0) > 0)
            }
            while (order<max(param[[2]])) {
                order = order+1
                answer = ((((tt %*% answer) > 0) - tt1 - tt0) > 0)
            }

            ##                not.implemented(mesh,query)
        } else if (identical(query, "vt.neighbours")) {
            if (is.null(param))
                param = list(1)
            if (length(param)<2)
                param = c(as.list(param), list(c(1)))
            nV = nrow(mesh$loc)
            nT = nrow(mesh$graph$tv)
            i = rep(1:nT,3)
            j = as.vector(mesh$graph$tv)
#            i = i[!is.na(j)]
#            j = j[!is.na(j)]
            tv = sparseMatrix(i=i, j=j, x=rep(1,length(i)), dims=c(nT,nV))
            vv = (sparseMatrix(i=param[[1]],
                               j=rep(1,length(param[[1]])),
                               x=1,
                               dims=c(nV,1)))
            vt = (tv %*% vv ) > 0
            vv0 = (vv == 0.5)*1
            vv1 = vv
            vt0 = (tv %*% vv0 ) > 0
            vt1 = vt
            order = 0
            while (order<min(param[[2]])) {
                order = order+1
                vv0 = vv1
                vv1 = vv
                vv = ((((mesh$graph$vv %*% vv) > 0) - vv1 - vv0) > 0)
                vt0 = vt1
                vt1 = vt
                vt = (((tv %*% vv) > 0) - vt1 - vt0 ) > 0
            }
            while (order<max(param[[2]])) {
                order = order+1
                vv = ((((mesh$graph$vv %*% vv) > 0) - vv1 - vv0) > 0)
                vt = (((((tv %*% vv) > 0) - vt1 - vt0 ) > 0 ) + vt) > 0
            }
            answer = vt

            ##                not.implemented(mesh,query)
        } else if (!identical(query, "")) {
            not.known(mesh,query)
        }
        ## Expand the result list:
        result[query.idx] = list(NULL)
        names(result)[query.idx] = query
        ## Set the answer:
        if (!is.null(answer))
            result[[query.idx]] = answer
    }

    return(result)
}

summary.inla.mesh <- function(object, verbose=FALSE, ...)
{
    x=object
    ## provides a summary for a mesh object
    inla.require.inherits(x, "inla.mesh", "'x'")

    ret = list(verbose=verbose)
    if (verbose) {
        ret = (c(ret, list(call=x$meta$call,
                           fmesher.args=x$meta$fmesher.args,
                           prefix=x$meta$prefix,
                           time = x$meta$time,
                           is.refined = x$meta$is.refined)))
    }
    ret = (c(ret, list(manifold=x$manifold,
                       nV=x$n,
                       nT=nrow(x$graph$tv),
                       xlim=range(x$loc[,1]),
                       ylim=range(x$loc[,2]),
                       zlim=range(x$loc[,3]))))
    if (is.na(inla.CRSargs(x$crs))) {
      ret <- c(ret, list(crs="N/A"))
    } else {
      ret <- c(ret, list(crs=inla.CRSargs(x$crs)))
    }

    my.segm <- function(x) {
        if (is.null(x))
            return(list(n=0, grps=NULL))
        n = max(0, nrow(x$idx))
        if (max(0, length(unique(x$grp)))>0) {
            grps = unique(x$grp)
        } else {
            grps = NULL
        }
        return(list(n=n, grps=grps))
    }
    if(!is.null(x$segm)) {
        ret = (c(ret, list(segm.bnd=my.segm(x$segm$bnd),
                           segm.int=my.segm(x$segm$int))))
    } else {
        ret = (c(ret, list(segm.bnd=my.segm(NULL),
                           segm.int=my.segm(NULL))))
    }

    class(ret) <- "summary.inla.mesh"
    return (ret)
}

print.summary.inla.mesh <- function(x, ...)
{
    my.print.proc_time <- function (x, ...)
    {
        if (!is.matrix(x)) {
            y = matrix(x,1,5)
        } else {
            y = x
        }
        for (k in 1:nrow(x)) {
            if (!is.na(y[k,4L]))
                y[k,1L] <- y[k,1L] + y[k,4L]
            if (!is.na(y[k,5L]))
                y[k,2L] <- y[k,2L] + y[k,5L]
        }
        y <- y[,1L:3L]
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

    cat("\nManifold:\t", x$manifold, "\n", sep="")
    cat("CRS:\t", x$crs, "\n", sep="")
    if (x$verbose) {
        cat("Refined:\t", x$is.refined, "\n", sep="")
    }
    cat("Vertices:\t", as.character(x$nV), "\n", sep="")
    cat("Triangles:\t", as.character(x$nT), "\n", sep="")

    my.print.segm <- function(x) {
        cat(as.character(x$n))
        if (!is.null(x$grps)) {
            n = length(x$grps)
            cat(" (", n, " group", inla.ifelse(n==1, "", "s"), sep="")
            if (n <= 10) {
                cat(":", x$grps, sep=" ")
            } else {
                cat(":", x$grps[1:10], "...", sep=" ")
            }
            cat(")")
        }
        cat("\n", sep="")
        return(invisible())
    }
    cat("Boundary segm.:\t")
    my.print.segm(x$segm.bnd)
    cat("Interior segm.:\t")
    my.print.segm(x$segm.int)
    cat("xlim:\t", x$xlim[1], " ", x$xlim[2], "\n", sep="")
    cat("ylim:\t", x$ylim[1], " ", x$ylim[2], "\n", sep="")
    cat("zlim:\t", x$zlim[1], " ", x$zlim[2], "\n", sep="")
    cat("\n")

    invisible(x)
}



inla.mesh.project <- function(...)
{
    UseMethod("inla.mesh.project")
}

inla.mesh.project.inla.mesh <- function(mesh, loc=NULL, field=NULL,
                                        crs=NULL, ...)
{
    inla.require.inherits(mesh, "inla.mesh", "'mesh'")

    if (!is.null(mesh$crs) &&
        identical(inla.as.list.CRS(mesh$crs)[["proj"]], "geocent")) {
      crs.sphere <- inla.CRS("sphere")
      if (!inla.identical.CRS(mesh$crs, crs.sphere)) {
        ## Convert the mesh to a perfect sphere.
        mesh$loc <- inla.spTransform(mesh$loc, mesh$crs, crs.sphere)
        mesh$manifold <- "S2"
        mesh$crs <- crs.sphere
      }
    }

    ## Handle loc given as SpatialPoints or SpatialPointsDataFrame object
    ## or an explicit crs
    if (!(missing(loc) || is.null(loc))) {
      if (inherits(loc, "SpatialPoints") ||
          inherits(loc, "SpatialPointsDataFrame")) {
        if (is.null(mesh$crs)) {
          loc <- coordinates(loc)
        } else {
          loc = inla.spTransform(coordinates(loc),
                                 CRS(proj4string(loc)),
                                 mesh$crs,
                                 passthrough=FALSE)
        }
      } else if (!is.null(crs)) {
        if (!is.null(mesh$crs)) {
          loc = inla.spTransform(loc,
                                 crs,
                                 mesh$crs,
                                 passthrough=FALSE)
        }
      }
    }

    if (!missing(field) && !is.null(field)) {
      proj <- inla.mesh.projector(mesh, loc=loc, ...)
      return(inla.mesh.project(proj, field=field))
    }

    jj =
        which(rowSums(matrix(is.na(as.vector(loc)),
                             nrow=nrow(loc),
                             ncol=ncol(loc))) == 0)
    smorg = (inla.fmesher.smorg(mesh$loc,
                                mesh$graph$tv,
                                points2mesh=loc[jj,,drop=FALSE]))
    ti = matrix(0L, nrow(loc), 1)
    b = matrix(0, nrow(loc), 3)
    ti[jj,1L] = smorg$p2m.t
    b[jj,] = smorg$p2m.b

    ok = (ti[,1L] > 0L)
    ti[ti[,1L] == 0L,1L] = NA

    ii = which(ok)
    A = (sparseMatrix(dims=c(nrow(loc),mesh$n),
                      i = rep(ii, 3),
                      j = as.vector(mesh$graph$tv[ti[ii,1L],]),
                      x = as.vector(b[ii,]) ))

    list(t=ti, bary=b, A=A, ok=ok)
}

inla.mesh.project.inla.mesh.1d <- function(mesh, loc, field=NULL, ...)
{
    inla.require.inherits(mesh, "inla.mesh.1d", "'mesh'")

    if (!missing(field) && !is.null(field)) {
        proj = inla.mesh.projector(mesh, loc, ...)
        return(inla.mesh.project(proj, field))
    }

    A = inla.mesh.1d.A(mesh, loc=loc)

    return (list(A=A,
                 ok=(loc >= mesh$interval[1]) & (loc <= mesh$interval[2])))
}


inla.mesh.project.inla.mesh.projector <-
    function(projector, field, ...)
{
    inla.require.inherits(projector, "inla.mesh.projector", "'projector'")

    if (is.data.frame(field)) {
        field = as.matrix(field)
    }

    if (is.null(dim(field))) {
        if (is.null(projector$lattice)) {
            data = as.vector(projector$proj$A %*% as.vector(field))
            data[!projector$proj$ok] = NA
            return(data)
        } else {
            data = as.vector(projector$proj$A %*% as.vector(field))
            data[!projector$proj$ok] = NA
            return(matrix(data,
                          projector$lattice$dims[1],
                          projector$lattice$dims[2]))
        }
    } else {
        data = projector$proj$A %*% field
        data[!projector$proj$ok,] = NA
        return(data)
    }
}


inla.mesh.projector <- function(...)
{
    UseMethod("inla.mesh.projector")
}

inla.mesh.projector.inla.mesh <-
    function(mesh,
             loc=NULL,
             lattice=NULL,
             xlim=NULL,
             ylim=NULL,
             dims=c(100,100),
             projection=NULL,
             crs=NULL,
             ...)
{
  inla.require.inherits(mesh, "inla.mesh", "'mesh'")

  if (missing(loc) || is.null(loc)) {
    if (missing(lattice) || is.null(lattice)) {
      if (identical(mesh$manifold, "R2") &&
          (is.null(mesh$crs) || is.null(crs))) {
        units = "default"
        lim <- list(xlim=range(mesh$loc[,1]),
                    ylim=range(mesh$loc[,2]))
      } else if (identical(mesh$manifold, "S2") &&
                 (is.null(mesh$crs) || is.null(crs))) {
        projection =
          match.arg(projection, c("longlat", "longsinlat",
                                  "mollweide"))
        units = projection
        lim = inla.mesh.map.lim(loc=mesh$loc, projection=projection)
      } else {
        lim <- inla.crs.bounds(crs)
        if (identical(mesh$manifold, "R2")) {
          lim0 <- list(xlim=range(mesh$loc[,1]),
                       ylim=range(mesh$loc[,2]))
          lim$xlim[1] <- max(lim$xlim[1], lim0$xlim[1])
          lim$xlim[2] <- min(lim$xlim[2], lim0$xlim[2])
          lim$ylim[1] <- max(lim$ylim[1], lim0$ylim[1])
          lim$ylim[2] <- min(lim$ylim[2], lim0$ylim[2])
        }
      }
      if (missing(xlim) && is.null(xlim)) {
        xlim = lim$xlim
      }
      if (missing(ylim) && is.null(ylim)) {
        ylim = lim$ylim
      }
      x = seq(xlim[1], xlim[2], length.out=dims[1])
      y = seq(ylim[1], ylim[2], length.out=dims[2])
      if (is.null(mesh$crs) || is.null(crs)) {
        lattice <- inla.mesh.lattice(x=x, y=y, units = units)
      } else {
        lattice <- inla.mesh.lattice(x=x, y=y, crs=crs)
      }
    } else {
      dims = lattice$dims
      x = lattice$x
      y = lattice$y
    }

    if (is.null(mesh$crs) || is.null(lattice$crs)) {
      proj = inla.mesh.project(mesh, lattice$loc)
    } else {
      proj <- inla.mesh.project(mesh,
                                loc=lattice$loc,
                                crs=lattice$crs)
    }
    projector = list(x=x, y=y, lattice=lattice, loc=NULL, proj=proj, crs=crs)
    class(projector) = "inla.mesh.projector"
  } else {
    proj = inla.mesh.project(mesh, loc=loc, crs=crs)
    projector = list(x=NULL, y=NULL, lattice=NULL, loc=loc, proj=proj, crs=crs)
    class(projector) = "inla.mesh.projector"
  }

  return (projector)
}


inla.mesh.projector.inla.mesh.1d <-
    function(mesh,
             loc=NULL,
             xlim=mesh$interval,
             dims=100,
             ...)
  {
    inla.require.inherits(mesh, "inla.mesh.1d", "'mesh'")

    if (missing(loc) || is.null(loc)) {
        loc = seq(xlim[1], xlim[2], length.out=dims[1])
    }

    proj = inla.mesh.project(mesh, loc)
    projector = list(x=loc, lattice=NULL, loc=loc, proj=proj)
    class(projector) = "inla.mesh.projector"

    return (projector)
  }




inla.internal.make.spline.mesh <-
    function(interval, m, degree, boundary, free.clamped)
{
    boundary =
        match.arg(boundary,
                  c("neumann", "dirichlet", "free", "cyclic"))
    if (degree <= 1) {
        n = (switch(boundary,
                    neumann = m,
                    dirichlet = m+2,
                    free = m,
                    cyclic = m+1))
        if (n<2) {
            n = 2
            degree = 0
            boundary = "c"
        }
    } else {
        stopifnot(degree==2)
        n = (switch(boundary,
                    neumann = m+1,
                    dirichlet = m+1,
                    free = m-1,
                    cyclic = m))
        if (boundary=="free") {
            if (m <= 1) {
                n = 2
                degree = 0
                boundary = "c"
            } else if (m == 2) {
                n = 2
                degree = 1
            }
        } else if (boundary=="cyclic") {
            if (m <= 1) {
                n = 2
                degree = 0
            }
        }
    }
    return(inla.mesh.1d(seq(interval[1], interval[2], length=n),
                        degree=degree,
                        boundary=boundary,
                        free.clamped=free.clamped))
}


inla.mesh.basis <- function(mesh,
                            type="b.spline",
                            n=3,
                            degree=2,
                            knot.placement="uniform.area",
                            rot.inv=TRUE,
                            boundary="free",
                            free.clamped=TRUE,
                            ...)
{
    inla.require.inherits(mesh, c("inla.mesh", "inla.mesh.1d"), "'mesh'")

    type = match.arg(type, c("b.spline", "sph.harm"))
    knot.placement = (match.arg(knot.placement,
                                c("uniform.area",
                                  "uniform.latitude")))

    if (identical(type, "b.spline")) {


        if (identical(mesh$manifold, "R1") || identical(mesh$manifold, "S1")) {
            mesh1 =
                inla.internal.make.spline.mesh(mesh$interval, n, degree,
                                               boundary, free.clamped)
            basis = inla.mesh.1d.A(mesh1, mesh$loc)
        } else if (identical(mesh$manifold, "R2")) {
            if (length(n) == 2) {
                if (length(degree)==1) {
                    degree = rep(degree, 2)
                }
                if (length(boundary)==1) {
                    boundary = rep(boundary, 2)
                }
                mesh1x =
                    inla.internal.make.spline.mesh(range(mesh$loc[,1]),
                                                   n[1], degree[1],
                                                   boundary[1], free.clamped)
                mesh1y =
                    inla.internal.make.spline.mesh(range(mesh$loc[,2]),
                                                   n[2], degree[2],
                                                   boundary[2], free.clamped)
                basis =
                    inla.row.kron(inla.mesh.1d.A(mesh1y, mesh$loc[,2]),
                                  inla.mesh.1d.A(mesh1x, mesh$loc[,1]))
            } else {
                warning("Old style call for an R2 mesh detected.\n  Please supply a 2-element n-vector instead.")
                long = ((mesh$loc[,1]-min(mesh$loc[,1]))/
                        diff(range(mesh$loc[,1])))*90*pi/180
                sinlat =  (((mesh$loc[,2]-min(mesh$loc[,2]))/
                            diff(range(mesh$loc[,2])))*2-1)
                coslat = sapply(sinlat, function(x) sqrt(max(0.0,1.0-x^2)))
                loc = (matrix(c(cos(long)*coslat,
                                sin(long)*coslat,
                                sinlat), mesh$n, 3))
                knots = 0
                degree = max(0L, min(n-1L, degree))
                basis = (inla.fmesher.smorg(loc,
                                            mesh$graph$tv,
                                            bspline = c(n, degree, knots),
                                            fem=-1)$bspline)
            }
        } else if (identical(mesh$manifold, "S2")) {
            loc = mesh$loc
            knots = identical(knot.placement, "uniform.latitude")
            degree = max(0L, min(n-1L, degree))
            basis = (inla.fmesher.smorg(loc,
                                        mesh$graph$tv,
                                        bspline = c(n, degree, knots),
                                        fem=-1)$bspline)
            if (!rot.inv) {
                warning("Currently only 'rot.inv=TRUE' is supported for B-splines.")
            }
        } else {
            stop("Only know how to make B-splines on R2 and S2.")
        }
    } else if (identical(type, "sph.harm")) {
        if (!identical(mesh$manifold, "S2")) {
            stop("Only know how to make spherical harmonics on S2.")
        }
        if (rot.inv) {
            basis = (inla.fmesher.smorg(mesh$loc,
                                        mesh$graph$tv,
                                        sph0=n)$sph0)
        } else {
            basis = (inla.fmesher.smorg(mesh$loc,
                                        mesh$graph$tv,
                                        sph=n)$sph)
        }
    }

    return(basis)
}
















inla.parse.queries <-function(...)
{
    queries = list(...)
    if (length(queries)==0)
        return(queries)

    ## Make sure that we have a list of names, empty or not:
    if (is.null(names(queries))) {
        q.names = rep("", length(queries))
    } else {
        q.names = names(queries)
    }

    ## All nameless entries must be strings with query names.  Replace
    ## empty names with those names, and set those entries to NULL.
    for (query.idx in 1:length(queries)) {
        if (q.names[[query.idx]]=="") {
            if (is.character(queries[[query.idx]])) {
                names(queries)[query.idx] = queries[[query.idx]]
                queries[query.idx] = list(NULL)
            } else {
                queries[query.idx] = list(NULL)
                warning(paste("Unnamed query ignored.  Check query #",
                              query.idx, sep=""))
            }
        }
    }

    return(queries)
}










`inla.fmesher.smorg` <- function(loc, tv,
                                 fem=NULL,
                                 aniso=NULL,
                                 gradients=FALSE,
                                 sph0=NULL,
                                 sph=NULL,
                                 bspline=NULL,
                                 points2mesh=NULL,
                                 splitlines=NULL,
                                 output=NULL,
                                 keep=FALSE)
{
    prefix = inla.fmesher.make.prefix(NULL, NULL)

    n = nrow(loc)
    s.dim = ncol(loc)

    if (s.dim==1)
        stop("1-D models not implemented yet.")
    stopifnot(s.dim>=2, s.dim<=3)

    if (missing(output) || is.null(output)) {
        output.given = FALSE
        output = NULL
    } else {
        output.given = TRUE
    }
    output.fem = list("c0", "g1", "g2")
    output.aniso = list("g1aniso", "g2aniso")
    output.gradients = list("dx", "dy", "dz")
    output.sph0 = list("sph0")
    output.sph = list("sph")
    output.bspline = list("bspline")
    output.p2m = list("p2m.t", "p2m.b")
    output.splitlines <- list("split.loc", "split.idx", "split.origin",
                              "split.t", "split.b1", "split.b2")

    ## Outputs that need +1L index adjustment
    indexoutput <- list("split.idx", "split.t", "split.origin")

    fmesher.write(inla.affirm.double(loc), prefix, "s")
    fmesher.write(inla.affirm.integer(tv)-1L, prefix, "tv")
    all.args = "--smorg --input=s,tv"

    ## additional arguments
    if (!is.null(fem)) {
        all.args = paste(all.args," --fem=", max(2,fem), sep="")
        if (!output.given) output = c(output, output.fem)
    }
    if (!is.null(aniso)) {
        fmesher.write(inla.affirm.double(aniso[[1]]), prefix, "aniso.gamma")
        fmesher.write(inla.affirm.double(aniso[[2]]), prefix, "aniso.vec")

        all.args = paste(all.args," --aniso=aniso.gamma,aniso.vec", sep="")
        if (!output.given) output = c(output, output.aniso)
    }
    if (!is.null(gradients) && gradients) {
        all.args = paste(all.args," --grad", sep="")
        if (!output.given) output = c(output, output.gradients)
    }
    if (!is.null(sph0)) {
        all.args = paste(all.args," --sph0=", sph0, sep="")
        if (!output.given) output = c(output, output.sph0)
    }
    if (!is.null(sph)) {
        all.args = paste(all.args," --sph=", sph, sep="")
        if (!output.given) output = c(output, output.sph)
    }
    if (!is.null(bspline)) {
        all.args = (paste(all.args, " --bspline=",
                          bspline[1], ",", bspline[2], ",", bspline[3],
                          sep=""))
        if (!output.given) output = c(output, output.bspline)
    }
    if (!is.null(points2mesh)) {
        fmesher.write(inla.affirm.double(points2mesh), prefix, "p2m")

        all.args = paste(all.args," --points2mesh=p2m", sep="")
        if (!output.given) output = c(output, output.p2m)
    }
    if (!is.null(splitlines)) {
        fmesher.write(inla.affirm.double(splitlines$loc), prefix, "splitloc")
        fmesher.write(inla.affirm.integer(splitlines$idx)-1L, prefix, "splitidx")

        all.args = paste(all.args," --splitlines=splitloc,splitidx", sep="")
        if (!output.given) output = c(output, output.splitlines)
    }
    all.args = paste(all.args, inla.getOption("fmesher.arg"))

    echoc = inla.fmesher.call(all.args=all.args, prefix=prefix)

    result = list()
    for (name in output) {
        if (identical(name, "p2m.t"))
            if (!file.exists(paste(prefix, name, sep="")))
                result[[name]] = fmesher.read(prefix, "points2mesh.t")+1L
            else
                result[[name]] = fmesher.read(prefix, name) + 1L
        else if (identical(name, "p2m.b"))
            if (!file.exists(paste(prefix, name, sep="")))
                result[[name]] = fmesher.read(prefix, "points2mesh.b")
            else
                result[[name]] = fmesher.read(prefix, name)
        else if (name %in% indexoutput)
          result[[name]] = fmesher.read(prefix, name) + 1L
        else
          result[[name]] = fmesher.read(prefix, name)
    }

    if (!keep)
        unlink(paste(prefix, "*", sep=""), recursive=FALSE)

    return (result)
}







## Deprecated: cyclic
inla.mesh.1d <- function(loc,
                         interval=range(loc),
                         boundary=NULL,
                         degree=1,
                         free.clamped=FALSE,
                         ...)
{
    ## Note: do not change the order of these options without also
    ## changing 'basis.reduction' below.
    boundary.options = c("neumann", "dirichlet", "free", "cyclic")

    if ("cyclic" %in% names(list(...))) {
        warning(paste("Option 'cyclic' is deprecated.",
                      "  Use 'boundary=\"cyclic\"' instead.", sep=""))
        if (cyclic) {
            if (!missing(boundary)) {
                boundary = match.arg.vector(boundary, boundary.options, length=2)
                warning(paste("'cyclic=TRUE' overrides 'boundary=c(\"",
                              paste(boundary, collapse="\",\""), "\")'.", sep=""))
            }
            boundary = "cyclic"
        }
    }
    boundary = match.arg.vector(boundary, boundary.options, length=2)
    cyclic = !is.na(pmatch(boundary[1], "cyclic"))
    if (cyclic && is.na(pmatch(boundary[2], "cyclic"))) {
        stop("Inconsistent boundary specification 'boundary=c(",
             paste(boundary, collapse=","), ")'.", sep="")
    }

    loc.orig = loc
    if (cyclic)
        loc =
            (sort(unique(c(0, loc-interval[1]) %% diff(interval))) +
             interval[1])
    else
        loc =
            (sort(unique(c(interval,
                           pmax(interval[1], pmin(interval[2], loc)) ))))

    n = length(loc)

    if (loc[1]<interval[1])
        stop("All 'loc' must be >= interval[1].")
    if (loc[n]>interval[2])
        stop("All 'loc' must be <= interval[2].")

    if ( (degree<0) || (degree>2)) {
        stop(paste("'degree' must be 0, 1, or 2.  'degree=",
                   degree,
                   "' is not supported.", sep=""))
    }

    if (length(free.clamped) == 1L) {
        free.clamped = rep(free.clamped, 2)
    }


    ## Number of basis functions
    if (degree==0) {
        basis.reduction = c(0, 1, 0, 1/2) ## neu, dir, free, cyclic
    } else if (degree==1) {
        basis.reduction = c(0, 1, 0, 1/2) ## neu, dir, free, cyclic
    } else {
        basis.reduction = c(1, 1, 0, 1)
    }
    m = (n+cyclic+(degree==2)*1
         -basis.reduction[pmatch(boundary[1], boundary.options)]
         -basis.reduction[pmatch(boundary[2], boundary.options)])
##    if (m < 1+max(1,degree)) {
    if (m < 1L) {
        stop("Degree ", degree,
             " meshes must have at least ", 1L,
             " basis functions, not 'm=", m, "'.", sep="")
    }

    ## Compute representative basis midpoints.
    if ((degree==0) || (degree==1)) {
        mid = loc
        if (boundary[1] == "dirichlet") {
            mid = mid[-1]
        }
        if (boundary[2] == "dirichlet") {
            mid = mid[-(m+1)]
        }
    } else { ## degree==2
        if (cyclic) {
            mid = (loc + c(loc[-1], interval[2]))/2
        } else {
            mid = c(loc[1], (loc[-n]+loc[-1])/2, loc[n])
            mid =
                switch(boundary[1],
                       neumann = mid[-1],
                       dirichlet = mid[-1],
                       free = mid)
            mid =
                switch(boundary[2],
                       neumann = mid[-(m+1)],
                       dirichlet = mid[-(m+1)],
                       free = mid)
        }
    }

    mesh =
        list(n=n,
             m=m,
             loc=loc,
             mid=mid,
             interval=interval,
             boundary=boundary,
             cyclic=cyclic,
             manifold = inla.ifelse(cyclic, "S1", "R1"),
             degree=degree,
             free.clamped=free.clamped,
             idx = list(loc=NULL))
    class(mesh) = "inla.mesh.1d"

    if (degree<2) {
        mesh$idx$loc =
            inla.mesh.1d.bary(mesh, loc.orig, method="nearest")$index[,1]
    } else {
        if (length(mid) >= 2) {
            mesh$idx$loc =
                inla.mesh.1d.bary(inla.mesh.1d(mid, degree=0),
                                  loc.orig,
                                  method="nearest")$index[,1]
        } else {
            mesh$idx$loc = rep(1, length(loc.orig))
        }
    }

    return(invisible(mesh))
}

inla.mesh.1d.bary <- function(mesh, loc, method=c("linear", "nearest"))
{
    inla.require.inherits(mesh, "inla.mesh.1d", "'mesh'")
    method = match.arg(method)

    if (method=="linear") {
        if (mesh$cyclic) {
            mloc = c(mesh$loc-mesh$loc[1], diff(mesh$interval))
            loc = (loc-mesh$loc[1]) %% diff(mesh$interval)
        } else {
            mloc = c(mesh$loc-mesh$loc[1], diff(mesh$interval))
            loc = pmax(0, pmin(diff(mesh$interval), loc-mesh$loc[1]))
        }
    } else {
        if (mesh$cyclic) {
            mloc =
                c(mesh$loc[mesh$n]-diff(mesh$interval),
                  mesh$loc,
                  diff(mesh$interval))
            mloc = (mloc[-(mesh$n+2)] + mloc[-1])/2
            loc = (loc-mloc[1]) %% diff(mesh$interval)
            mloc = mloc-mloc[1]
        } else {
            mloc =
                c(0,
                  (mesh$loc[1:(mesh$n-1L)] +
                   mesh$loc[2:mesh$n])/2-mesh$loc[1],
                  diff(mesh$interval))
            loc = pmax(0, pmin(diff(mesh$interval), loc-mesh$loc[1]))
        }
    }

    ## Binary split method:
    do.the.split <- function(knots, loc) {
        n = length(knots)
        if (n <= 2L) {
            return(rep(1L, length(loc)))
        }
        split = 1L + (n-1L) %/% 2L ## Split point
        upper = (loc >= knots[split])
        idx = rep(0, length(loc))
        idx[!upper] = do.the.split(knots[1:split], loc[!upper])
        idx[upper] = split - 1L + do.the.split(knots[split:n], loc[upper])
        return(idx)
    }

    idx = do.the.split(mloc, loc)

    if (method=="nearest") {
        u = rep(0, length(loc))
        if (mesh$cyclic) {
            found = which(idx==(mesh$n+1L))
            idx[found] = 1L
        }
    } else { ## (method=="linear") {
        u = pmax(0, pmin(1, (loc-mloc[idx])/(mloc[idx+1L]-mloc[idx]) ))
        if (!mesh$cyclic) {
            found = which(idx==mesh$n)
            idx[found] = mesh$n-1L
            u[found] = 1
        }
    }

    if (mesh$cyclic) {
        index = matrix(c(idx, (idx %% mesh$n)+1L), length(idx), 2)
        bary = matrix(c(1-u, u), length(idx), 2)
    } else {
        index = matrix(c(idx, idx+1L), length(idx), 2)
        bary = matrix(c(1-u, u), length(idx), 2)
    }

    return(list(index=index, bary=bary))
}

inla.mesh.1d.A <- function(mesh, loc,
                           weights=NULL,
                           derivatives=NULL,
                           method=c("linear", "nearest", "quadratic")
                           )
{
    inla.require.inherits(mesh, "inla.mesh.1d", "'mesh'")
    if (missing(method)) {
        ## Compute basis based on mesh$degree and mesh$boundary
        if (mesh$degree==0) {
            info = (inla.mesh.1d.A(mesh, loc, method="nearest",
                                   weights=weights,
                                   derivatives=
                                   inla.ifelse(is.null(derivatives),
                                               FALSE, TRUE)))
            if (mesh$boundary[1] == "dirichlet") {
                info$A = info$A[,-1,drop=FALSE]
                if (!is.null(derivatives) && derivatives)
                    info$dA = info$dA[,-1,drop=FALSE]
            }
            if (mesh$boundary[2] == "dirichlet") {
                info$A = info$A[,-(mesh$m+1),drop=FALSE]
                if (!is.null(derivatives) && derivatives)
                    info$dA = info$dA[,-(mesh$m+1),drop=FALSE]
            }
            if (!is.null(derivatives) && derivatives)
                return(info)
            else
                return(info$A)
        } else if (mesh$degree==1) {
            info = (inla.mesh.1d.A(mesh, loc, method="linear",
                                   weights=weights,
                                   derivatives=
                                   inla.ifelse(is.null(derivatives),
                                               FALSE, TRUE)))
            if (mesh$boundary[1] == "dirichlet") {
                info$A = info$A[,-1,drop=FALSE]
                if (!is.null(derivatives) && derivatives)
                    info$dA = info$dA[,-1,drop=FALSE]
            }
            if (mesh$boundary[2] == "dirichlet") {
                info$A = info$A[,-(mesh$m+1),drop=FALSE]
                if (!is.null(derivatives) && derivatives)
                    info$dA = info$dA[,-(mesh$m+1),drop=FALSE]
            }
            if (!is.null(derivatives) && derivatives)
                return(info)
            else
                return(info$A)
        } else if (mesh$degree==2) {
            info =
                inla.mesh.1d.A(mesh, loc,
                               method="quadratic",
                               weights=weights,
                               derivatives=
                               inla.ifelse(
                                   is.null(derivatives),
                                   FALSE, TRUE))

            adjust.matrix <- function(mesh, A) {
                A = inla.as.dgTMatrix(A)
                i = A@i+1L
                j = A@j+1L
                x = A@x
                if (mesh$cyclic) {
                    j[j==1L] = mesh$m+1L
                    j[j==mesh$m+2L] = 2L
                    j = j-1L
                } else {
                    if (mesh$boundary[1] == "neumann") {
                        j[j==1L] = 2L
                        j = j-1L
                    } else if (mesh$boundary[1] == "dirichlet") {
                        x[j==1L] = -x[j==1L]
                        j[j==1L] = 2L
                        j = j-1L
                    } else if ((mesh$boundary[1] == "free") &&
                               mesh$free.clamped[1]) {
                        j1 = which(j==1L)
                        i = c(i, i[j1])
                        j = c(j, j[j1]+1L)
                        x = c(x, -x[j1])
                        x[j1] = 2*x[j1]
                    }
                    if (mesh$boundary[2] == "neumann") {
                        j[j==(mesh$m+1L)] = mesh$m
                    } else if (mesh$boundary[2] == "dirichlet") {
                        x[j==(mesh$m+1L)] = -x[j==(mesh$m+1L)]
                        j[j==(mesh$m+1L)] = mesh$m
                    } else if ((mesh$boundary[2] == "free") &&
                               mesh$free.clamped[2]) {
                        j1 = which(j==(mesh$m))
                        i = c(i, i[j1])
                        j = c(j, j[j1]-1L)
                        x = c(x, -x[j1])
                        x[j1] = 2*x[j1]
                    }
                }
                return(sparseMatrix(i=i, j=j, x=x,
                                    dims=c(length(loc), mesh$m)))
            }

            if (!is.null(derivatives) && derivatives) {
                return(list(A=adjust.matrix(mesh, info$A),
                            dA=adjust.matrix(mesh, info$dA),
                            d2A=adjust.matrix(mesh, info$d2A)))
            } else {
                return(adjust.matrix(mesh, info$A))
            }
        } else {
            stop(paste("'degree' must be 0, 1, or 2.  'degree=",
                       mesh$degree,
                       "' is not supported.", sep=""))
        }
    } else {
        method = match.arg(method)
        if (missing(weights) || is.null(weights)) {
            weights = rep(1, length(loc))
        }

        if (!is.na(pmatch(method, c("linear", "nearest")))) {
            idx = inla.mesh.1d.bary(mesh, loc, method)
            dA = NULL

            if (method=="linear") {
                ## Compute the n 1st order B-splines.
                i = rep(1:length(loc), times=2)
                weights.i = weights[i]
                A = (sparseMatrix(i=i,
                                  j=as.vector(idx$index),
                                  x=weights.i*as.vector(idx$bary),
                                  dims=c(length(loc), mesh$n)))
                if (!is.null(derivatives) && derivatives) {
                    if (mesh$cyclic) {
                        d = (c(mesh$loc[-1], mesh$interval[2]) -
                             mesh$loc)
                    } else {
                        d = (mesh$loc[-1] - mesh$loc[-mesh$n])
                    }
                    dA = (sparseMatrix(i=i,
                                       j=as.vector(idx$index),
                                       x=(weights.i*c(-1/d[idx$index[,1]],
                                                      1/d[idx$index[,1]])),
                                       dims=c(length(loc), mesh$n)))
                }
            } else {
                ## Nearest neighbours.
                A = (sparseMatrix(i=1:length(loc),
                                  j=idx$index[,1],
                                  x=weights*idx$bary[,1],
                                  dims=c(length(loc), mesh$n)))
            }

            if (missing(derivatives) || is.null(derivatives)) {
                return(A)
            } else {
                return(list(A=A, dA=dA))
            }
        } else { ## Quadratic
            ## Compute the n+1 2nd order B-splines.

            if (mesh$cyclic) {
                Boundary.knots =
                    (c(mesh$loc, mesh$interval[2])[c(mesh$n,2)] +
                     diff(mesh$interval)*c(-1,1))
                knots = c(mesh$loc, mesh$interval[2])
            } else {
                Boundary.knots =
                    c(mesh$loc[1] - (mesh$loc[2]-mesh$loc[1]),
                      mesh$loc[mesh$n] + (mesh$loc[mesh$n]-mesh$loc[mesh$n-1]))
                knots = mesh$loc
            }

            if (FALSE) { ## Only for debugging.
                ## Using bs():
                ## Note: Intermediate step constructs dense matrix.
                bsobj =
                    bs(x=loc,
                       knots=knots,
                       degree=2,
                       Boundary.knots=Boundary.knots)
                bsobj =
                    Matrix(as.vector(bsobj[,1:(mesh$n+mesh$cyclic+1)]),
                           nrow(bsobj), mesh$n+mesh$cyclic+1)
            }

            ## Direct calculation:
            knots = c(Boundary.knots[1], knots)
            if (mesh$cyclic) {
                idx =
                    inla.mesh.1d.bary(inla.mesh.1d(knots, boundary="free"),
                                      (loc-mesh$interval[1]) %%
                                      diff(mesh$interval) + mesh$interval[1],
                                      method="linear")
            } else {
                idx =
                    inla.mesh.1d.bary(inla.mesh.1d(knots, boundary="free"),
                                      loc, method="linear")
            }
            knots = c(knots, Boundary.knots[2])
            idx$index = idx$index[,1]-1L ## Indices into mesh intervals.

            d = knots[2:length(knots)]-knots[1:(length(knots)-1)]
            d2 = knots[3:length(knots)]-knots[1:(length(knots)-2)]

            ## Left intervals for each basis function:
            i.l = 1:length(idx$index)
            j.l = idx$index+2L
            x.l = (idx$bary[,2]*d[idx$index+1]/d2[idx$index+1] * idx$bary[,2])
            ## Right intervals for each basis function:
            i.r = 1:length(idx$index)
            j.r = idx$index
            x.r = (idx$bary[,1]*d[idx$index+1]/d2[idx$index] * idx$bary[,1])
            ## Middle intervals for each basis function:
            i.m = 1:length(idx$index)
            j.m = idx$index+1L
            x.m = (1-(idx$bary[,1]*d[idx$index+1]/d2[idx$index]* idx$bary[,1] +
                      idx$bary[,2]*d[idx$index+1]/d2[idx$index+1] * idx$bary[,2]
                      ))

            i = c(i.l, i.r, i.m)
            j = c(j.l, j.r, j.m)
            weights.i = weights[i]

            A = (sparseMatrix(i=i, j=j,
                              x=weights.i*c(x.l, x.r, x.m),
                              dims=c(length(idx$index), mesh$n+mesh$cyclic+1L)
                              ))

            dA = NULL
            d2A = NULL
            if (!is.null(derivatives) && derivatives) {
                ## dA:
                ## Left, right, middle intervals for each basis function:
                x.l = (2 / d2[idx$index+1] * idx$bary[,2])
                x.r = (-2 / d2[idx$index] * idx$bary[,1])
                x.m = (-(-2/ d2[idx$index] * idx$bary[,1] +
                          2/ d2[idx$index+1] * idx$bary[,2]
                          ))
                dA = (sparseMatrix(i=i, j=j,
                                  x=weights.i*c(x.l, x.r, x.m),
                                  dims=(c(length(idx$index),
                                          mesh$n+mesh$cyclic+1L))
                                  ))

                ## d2A:
                ## Left, right, middle intervals for each basis function:
                x.l = (2 / d[idx$index+1] / d2[idx$index+1])
                x.r = (2 / d[idx$index+1] / d2[idx$index])
                x.m = (-(2/d[idx$index+1] / d2[idx$index] +
                         2/d[idx$index+1] / d2[idx$index+1] ))
                d2A = (sparseMatrix(i=i, j=j,
                                    x=weights.i*c(x.l, x.r, x.m),
                                    dims=(c(length(idx$index),
                                            mesh$n+mesh$cyclic+1L))
                                  ))
            }

            if (missing(derivatives) || is.null(derivatives)) {
                return(A)
            } else {
                return(list(A=A, dA=dA, d2A=d2A))
            }
        }
    }
}

inla.mesh.1d.fem <- function(mesh)
{
    inla.require.inherits(mesh, "inla.mesh.1d", "'mesh'")

    ## Use the same matrices for degree 0 as for degree 1
    if ((mesh$degree==0) || (mesh$degree==1)) {
        if (mesh$cyclic) {
            loc =
                c(mesh$loc[mesh$n]-diff(mesh$interval),
                  mesh$loc,
                  mesh$loc[1]+diff(mesh$interval))
            c0 = (loc[3:length(loc)]-loc[1:(length(loc)-2)])/2
            c1.l = (loc[2:(length(loc)-1)]-loc[1:(length(loc)-2)])/6
            c1.r = (loc[3:length(loc)]-loc[2:(length(loc)-1)])/6
            c1.0 = (c1.l+c1.r)*2
            g1.l = -1/(loc[2:(length(loc)-1)]-loc[1:(length(loc)-2)])
            g1.r = -1/(loc[3:length(loc)]-loc[2:(length(loc)-1)])
            g1.0 = -g1.l-g1.r
            i.l = 1:mesh$n
            i.r = 1:mesh$n
            i.0 = 1:mesh$n
            if (mesh$n > 1) {
                j.l = c(mesh$n, 1:(mesh$n-1))
                j.r = c(2:mesh$n, 1)
                j.0 = 1:mesh$n
            } else {
                j.l = 1L
                j.r = 1L
                j.0 = 1L
            }
        } else {
            c0 =
                c((mesh$loc[2]-mesh$loc[1])/2,
                  (mesh$loc[mesh$n]-mesh$loc[mesh$n-1])/2
                  )
            if (mesh$n>2) {
                c0 =
                    c(c0[1],
                      (mesh$loc[3:mesh$n]-mesh$loc[1:(mesh$n-2)])/2,
                      c0[2]
                      )
            }
            c1.l = (mesh$loc[2:mesh$n]-mesh$loc[1:(mesh$n-1)])/6
            c1.r = c1.l
            c1.0 = (c(0, c1.l)+c(c1.r,0))*2
            g1.l = -1/(mesh$loc[2:mesh$n]-mesh$loc[1:(mesh$n-1)])
            g1.r = g1.l
            g1.0 = -c(0, g1.l)-c(g1.r,0)
            i.l = 2:mesh$n
            i.r = 1:(mesh$n-1)
            i.0 = 1:mesh$n
            j.l = 1:(mesh$n-1)
            j.r = 2:mesh$n
            j.0 = 1:mesh$n

            if (mesh$boundary[1] == "dirichlet") {
                g1.0 = g1.0[-1]; g1.l = g1.l[-1]; g1.r = g1.r[-1]
                c1.0 = c1.0[-1]; c1.l = c1.l[-1]; c1.r = c1.r[-1]
                c0 = c0[-1]
                i.l = i.l[-1]-1; i.r = i.r[-1]-1; i.0 = i.0[-1]-1
                j.l = j.l[-1]-1; j.r = j.r[-1]-1; j.0 = j.0[-1]-1
            } else if (mesh$boundary[1] == "free") {
                g1.0[1] = 0
                g1.r[1] = 0
            }
            if (mesh$boundary[2] == "dirichlet") {
                m = mesh$m
                g1.0 = g1.0[-(m+1)]; g1.l = g1.l[-m]; g1.r = g1.r[-m]
                c1.0 = c1.0[-(m+1)]; c1.l = c1.l[-m]; c1.r = c1.r[-m]
                c0 = c0[-(m+1)]
                i.l = i.l[-m]; i.r = i.r[-m]; i.0 = i.0[-(m+1)]
                j.l = j.l[-m]; j.r = j.r[-m]; j.0 = j.0[-(m+1)]
            } else if (mesh$boundary[2] == "free") {
                g1.0[mesh$m] = 0
                g1.l[mesh$m-1] = 0
            }
        }

        g1 =
            sparseMatrix(i=c(i.l, i.r, i.0),
                         j=c(j.l, j.r, j.0),
                         x=c(g1.l, g1.r, g1.0),
                         dims=c(mesh$m, mesh$m))
        c1 =
            sparseMatrix(i=c(i.l, i.r, i.0),
                         j=c(j.l, j.r, j.0),
                         x=c(c1.l, c1.r, c1.0),
                         dims=c(mesh$m, mesh$m))
        g2 = t(g1) %*% Diagonal(mesh$m, 1/c0) %*% g1
        c0 = Diagonal(mesh$m, c0)

    } else if (mesh$degree==2) {
        if (mesh$cyclic) {
            knots1 = mesh$loc
            knots2 = c(mesh$loc[-1], mesh$interval[2])
        } else {
            knots1 = mesh$loc[-mesh$n]
            knots2 = mesh$loc[-1]
        }
        knots.m = (knots1+knots2)/2
        knots.d = (knots2-knots1)/2
        ## 3-point Gaussian quadrature
        info =
            inla.mesh.1d.A(mesh,
                           loc=(c(knots.m,
                                  knots.m - knots.d*sqrt(3/5),
                                  knots.m + knots.d*sqrt(3/5))),
                           weights =
                           c(knots.d*8/9, knots.d*5/9, knots.d*5/9)^0.5,
                           derivatives = TRUE
                           )
        c1 = t(info$A) %*% info$A
        g1 = t(info$dA) %*% info$dA
        g2 = t(info$d2A) %*% info$d2A

        g01 = t(info$A) %*% info$dA
        g02 = t(info$A) %*% info$d2A
        g12 = t(info$dA) %*% info$d2A

        c0 = Diagonal(nrow(c1), rowSums(c1))

    return(list(c0=c0, c1=c1, g1=g1, g2=g2, g01=g01, g02=g02, g12=g12))

    } else {
        stop(paste("Mesh basis degree=", mesh$degree,
                   " is not supported.", sep=""))
    }

    return(list(c0=c0, c1=c1, g1=g1, g2=g2))
}









inla.diameter <- function(x, ...) {
  UseMethod("inla.diameter")
}

## Calculate upper bound for the diameter of a point set,
## by encapsulating in a circular domain.
inla.diameter.default <- function(x, manifold="", ...) {
  if (nrow(x) <= 1) {
    0
  } else {
    if (identical(manifold, "S2")) {
      distance <- function(u,v) {
        2*asin(pmin(1,
          ((u[1]-v[,1])^2 + (u[2]-v[,2])^2 + (u[3]-v[,3])^2)^0.5 / 2))
      }
      center <- colMeans(x)
      tmp <- sqrt(sum(center^2))
      if (tmp < 1e-6) {
        pi
      } else {
        center <- center/tmp
        min(pi, 2 * max(distance(center, x)))
      }
    } else {
      distance <- function(u,v) {
        d <- 0
        for (k in seq_len(ncol(v))) {
          d <- d + (u[k]-v[,k])^2
        }
        d^0.5
      }
      center <- rep(0, ncol(x))
      for (k in seq_len(ncol(x))) {
        center[k] <- mean(range(x[,k]))
      }
      2 * max(distance(center, x))
    }
  }
}

inla.diameter.inla.mesh.1d <- function(x, ...) {
  diff(x$interval)
}

inla.diameter.inla.mesh <- function(x, ...) {
  inla.diameter(x$loc, manifold=x$manifold, ...)
}

inla.diameter.inla.mesh.segment <- function(x, ...) {
  inla.diameter(x$loc, ...)
}

inla.diameter.inla.mesh.lattice <- function(x, ...) {
  inla.diameter(x$loc, ...)
}


inla.mesh.fem <- function(mesh, order=2)
{
    inla.require.inherits(mesh, c("inla.mesh", "inla.mesh.1d"), "'mesh'")
    if (inherits(mesh, "inla.mesh.1d")) {
        if (order > 2) {
            stop("order>2 not supported for 1d meshes.")
            ## TODO: Compute higher order matrices based on order 2 output.
        }
        return(inla.mesh.1d.fem(mesh))
    } else {
        ## output name list:
        if (order > 0) {
          output <- paste("g", seq_len(order), sep="")
        } else {
          output <- NULL
        }
        output = c("c0", "c1", output,
                   "va", "ta")
        return(inla.fmesher.smorg(mesh$loc, mesh$graph$tv,
                                  fem=order, output=output))
    }
}



inla.mesh.deriv <- function(mesh, loc)
{
    inla.require.inherits(mesh, c("inla.mesh", "inla.mesh.1d"), "'mesh'")
    if (inherits(mesh, "inla.mesh.1d")) {
        stop("1d derivatives are not implemented.")
    }
    n.mesh = mesh$n

    info = inla.mesh.project(mesh, loc=loc)

    ii = which(info$ok)
    n.ok = sum(info$ok)
    tv = mesh$graph$tv[info$t[ii,1L],]
    e1 = mesh$loc[tv[,3],]-mesh$loc[tv[,2],]
    e2 = mesh$loc[tv[,1],]-mesh$loc[tv[,3],]
    e3 = mesh$loc[tv[,2],]-mesh$loc[tv[,1],]
    n1 = e2-e1*matrix(rowSums(e1*e2)/rowSums(e1*e1),n.ok,3)
    n2 = e3-e2*matrix(rowSums(e2*e3)/rowSums(e2*e2),n.ok,3)
    n3 = e1-e3*matrix(rowSums(e3*e1)/rowSums(e3*e3),n.ok,3)
    g1 = n1/matrix(rowSums(n1*n1),n.ok,3)
    g2 = n2/matrix(rowSums(n2*n2),n.ok,3)
    g3 = n3/matrix(rowSums(n3*n3),n.ok,3)
    x = cbind(g1[,1],g2[,1],g3[,1])
    y = cbind(g1[,2],g2[,2],g3[,2])
    z = cbind(g1[,3],g2[,3],g3[,3])
    dx = (sparseMatrix(dims=c(nrow(loc),n.mesh),
                       i = rep(ii, 3),
                       j = as.vector(tv),
                       x = as.vector(x) ))
    dy = (sparseMatrix(dims=c(nrow(loc),n.mesh),
                       i = rep(ii, 3),
                       j = as.vector(tv),
                       x = as.vector(y) ))
    dz = (sparseMatrix(dims=c(nrow(loc),n.mesh),
                       i = rep(ii, 3),
                       j = as.vector(tv),
                       x = as.vector(z) ))

    return (list(A=info$A, dx=dx, dy=dy, dz=dz))
}









inla.simplify.curve <- function(loc, idx, eps) {
    ## Variation of Ramer-Douglas-Peucker
    ## Uses width epsilon ellipse instead of rectangle,
    ## motivated by prediction ellipse for Brownian bridge
    n = length(idx)
    if ((n==2) || (eps==0)) {
        return(idx)
    }
    segm = loc[idx[n],]-loc[idx[1],]
    segm.len = sum(segm^2)^0.5
    if (segm.len <= 1e-12) {
        ## End point same as start; closed curve.  Split.
        len2 = ((loc[idx[2:(n-1)], 1] - loc[idx[1], 1])^2 +
                (loc[idx[2:(n-1)], 2] - loc[idx[1], 2])^2 )
        split = which.max(len2)+1L
    } else {
        segm.mid = (loc[idx[n],]+loc[idx[1],])/2
        segm = segm/segm.len
        segm.perp = c(-segm[2], segm[1])
        vec = (cbind(loc[idx[2:(n-1)], 1] - segm.mid[1],
                     loc[idx[2:(n-1)], 2] - segm.mid[2]))
        ## Always split if any point is outside the circle
        epsi = min(c(eps, segm.len/2))
        dist1 = abs(vec[,1]*segm[1] + vec[,2]*segm[2])/(segm.len/2)*epsi
        dist2 = abs(vec[,1]*segm.perp[1] + vec[,2]*segm.perp[2])
        dist = (dist1^2 + dist2^2)^0.5

        ## Find the furthest point, in the ellipse metric, and
        ## check if it inside the radius (radius=segm.len/2)
        split = which.max(dist)+1L
        if (dist[split-1L] < epsi) {
            ## Flat segment, eliminate.
            return(idx[c(1, n)])
        }
        ## Split at the furthest point.
        split = which.max(dist)+1L
    }

    ## Do the split recursively:
    return(c(inla.simplify.curve(loc, idx[1L:split], eps),
             inla.simplify.curve(loc, idx[split:n], eps)[-1L]
             ))
}


inla.contour.segment <-
    function(x = seq(0, 1, length.out = nrow(z)),
             y = seq(0, 1, length.out = ncol(z)),
             z, nlevels = 10,
             levels = pretty(range(z, na.rm=TRUE), nlevels),
             groups = seq_len(length(levels)),
             positive = TRUE,
             eps = NULL,
             crs = NULL)
{
    ## Input checking from contourLines:
    if (missing(z)) {
        if (!missing(x)) {
            if (is.list(x)) {
                z <- x$z
                y <- x$y
                x <- x$x
            }
            else {
                z <- x
                x <- seq.int(0, 1, length.out = nrow(z))
            }
        }
        else stop("no 'z' matrix specified")
    }
    else if (is.list(x)) {
        y <- x$y
        x <- x$x
    }

    if (is.null(eps)) {
        eps = min(c(min(diff(x)), min(diff(y))))/8
    }
    ## End of input checking.

    ## Get contour pieces
    curves = contourLines(x,y,z,levels=levels)

    ## Make a mesh for easy gradient interpolation:
    latt = inla.mesh.lattice(x,y)
    mesh =
        inla.mesh.create(lattice=latt,
                         boundary=latt$segm,
                         extend=(list(n=3,
                                      offset=(max(diff(range(x)),
                                                  diff(range(y)))*0.1)
                                      ))
                         )
    ## Map function values to mesh indexing:
    zz = rep(0, prod(dim(z)))
    zz[mesh$idx$lattice] = as.vector(z)

    ## Mapping from level to group value:
    level2grp = function(level) {
        if (length(groups)==1)
            return(groups)
        for (k in seq_along(groups)) {
            if (levels[k] == level) {
                return(groups[k])
            }
        }
        return(0)
    }

    ## Join all contour pieces into a single mesh.segment
    ## Different levels can later be identified via the grp indices.
    loc = matrix(0,0,2)
    idx = matrix(0,0,2)
    grp = c()
    for (k in seq_len(length(curves))) {
        curve.loc = cbind(curves[[k]]$x, curves[[k]]$y)
        curve.n = nrow(curve.loc)

        ## Extract the rotated gradients along the curve
        curve.mid = (curve.loc[1:(curve.n-1),]+curve.loc[2:curve.n,])/2
        A = inla.mesh.deriv(mesh, loc=curve.mid)
        ## Gradients rotated 90 degrees CW, i.e. to the direction
        ## of CCW curves around positive excursions:
        grid.diff = cbind(A$dy %*% zz, -A$dx %*% zz)

        ## Determine the CCW/CW orientation
        curve.diff = diff(curve.loc)
        ccw = (sum(curve.diff * grid.diff) >= 0) ## True if in CCW direction
        if ((ccw && positive) | (!ccw && !positive)) {
            curve.idx = 1:curve.n
        } else {
            curve.idx = curve.n:1
        }

        ## Filter short line segments:
        curve.idx =
            inla.simplify.curve(curve.loc,
                                curve.idx,
                                eps = eps)

        ## Reorder, making sure any unused points are removed:
        curve.loc = curve.loc[curve.idx,,drop=FALSE]
        curve.n = nrow(curve.loc)
        curve.idx = cbind(1L:(curve.n-1L), 2L:curve.n)

        ## Check if the curve is closed, and adjust if it is:
        if (max(abs(curve.loc[1,] - curve.loc[curve.n,])) < 1e-12) {
            curve.loc = curve.loc[-curve.n,,drop=FALSE]
            curve.n = nrow(curve.loc)
            curve.idx = cbind(1L:curve.n, c(2L:curve.n, 1L))
        }

        ## Add the curve:
        offset = nrow(loc)
        loc = rbind(loc, curve.loc)
        idx = rbind(idx, curve.idx+offset)
        grp = c(grp, rep(level2grp(curves[[k]]$level), curve.n-1L))
    }

  inla.mesh.segment(loc=loc, idx=idx, grp=grp, is.bnd=FALSE, crs=crs)
}


## Based on an idea from Elias Teixeira Krainski
## Requires  splancs::nndistF
inla.nonconvex.hull.basic <-
    function(points, convex=-0.15, resolution=40, eps=NULL, crs=NULL)
{
  if (!(missing(points) || is.null(points)) &&
      (inherits(points, "SpatialPoints") ||
       inherits(points, "SpatialPointsDataFrame"))) {
    points <- inla.spTransform(coordinates(points),
                               CRS(proj4string(points)),
                               crs,
                               passthrough=TRUE)
  }

  if (length(convex)==1)
    convex = rep(convex,2)
  if (length(resolution)==1)
    resolution = rep(resolution,2)

    lim = rbind(range(points[,1]), range(points[,2]))
    ex = convex
    if (convex[1]<0) {ex[1] = -convex[1]*diff(lim[1,])}
    if (convex[2]<0) {ex[2] = -convex[2]*diff(lim[2,])}

    domain = c(diff(lim[1,]), diff(lim[2,])) + 2*ex
    dif = domain/(resolution-1)
    if (any(dif > min(convex))) {
        req.res = ceiling(domain/convex+1)
        warning(paste("Resolution (",
                      paste(resolution,collapse=","),
                      ") too small for convex (",
                      paste(convex,collapse=","),
                      ").\n",
                      "Resolution >=(",
                      paste(req.res,collapse=","),
                      ") required for more accurate results.",
                      sep=""))
    }

    ax = list(seq(lim[1,1] - ex[1], lim[1,2] + ex[1], length=resolution[1]),
        seq(lim[2,1] - ex[2], lim[2,2] + ex[2], length=resolution[2]))
    xy = as.matrix(expand.grid(ax[[1]], ax[[2]]))
    tr = diag(c(1/ex[1],1/ex[2]))

    stopifnot(inla.require("splancs"))
    z = (matrix(splancs::nndistF(points%*%tr, xy%*%tr),
                resolution[1], resolution[2]))
    segm = inla.contour.segment(ax[[1]], ax[[2]], z,
        levels=c(1), positive=FALSE, eps=eps, crs=crs)
    return(segm)
}


## Morphological dilation by "convex",
## followed by closing by "concave", with
## minimum concave curvature radius "concave".
## If the dilated set has no gaps of width between
##   2*convex*(sqrt(1+2*concave/convex) - 1)
## and
##   2*concave,
## then the minimum convex curvature radius is "convex".
## Default is concave=convex
## Special case concave=0 delegates to inla.nonconvex.hull.basic()
##
## The implementation is based on the identity
##   dilation(a) & closing(b) = dilation(a+b) & erosion(b)
## where all operations are with respect to disks with the specified radii.
inla.nonconvex.hull <-
  function(points, convex=-0.15, concave=convex, resolution=40, eps=NULL,
           crs=NULL)
{
  if (!(missing(points) || is.null(points)) &&
      (inherits(points, "SpatialPoints") ||
       inherits(points, "SpatialPointsDataFrame"))) {
    points <- inla.spTransform(coordinates(points),
                               CRS(proj4string(points)),
                               crs,
                               passthrough=TRUE)
  }

    if (length(resolution)==1)
        resolution = rep(resolution,2)
    lim = rbind(range(points[,1]), range(points[,2]))

    approx.diam = max(diff(lim[1,]),diff(lim[2,]))
    if (convex<0) {convex = -convex*approx.diam}
    if (concave<0) {concave = -concave*approx.diam}
    if (concave==0) {
      return(inla.nonconvex.hull.basic(points, convex, resolution, eps,
                                       crs=crs))
    }

    ex = convex+concave
    domain = c(diff(lim[1,]), diff(lim[2,])) + 2*ex
    dif = domain/(resolution-1)
    if (max(dif) > min(convex,concave)) {
        req.res = ceiling(domain/min(convex,concave)+1)
        warning(paste("Resolution (",
                      paste(resolution,collapse=","),
                      ") too small for convex/concave radius (",
                      convex,",",concave,
                      ").\n",
                      "Resolution >=(",
                      paste(req.res,collapse=","),
                      ") required for more accurate results.",
                      sep=""))
    }
    ax =
        list(
            seq(lim[1,1] - ex, lim[1,2] + ex, length=resolution[1]),
            seq(lim[2,1] - ex, lim[2,2] + ex, length=resolution[2])
            )
    xy = as.matrix(expand.grid(ax[[1]], ax[[2]]))

    stopifnot(inla.require("splancs"))
    z = (matrix(splancs::nndistF(points, xy),
                resolution[1],resolution[2]))
    segm.dilation =
        inla.contour.segment(ax[[1]], ax[[2]], z,
                             levels=c(convex+concave),
                             positive=TRUE,
                             eps=0) ## Don't simplify curve at this stage
    mesh.dilation =
        inla.mesh.create(loc=xy,
                         boundary=segm.dilation,
                         extend=(list(n=3,
                                      offset=(max(diff(ax[[1]]),
                                                  diff(ax[[2]]))*0.1)
                                      ))
                         )

    ## This filtering is not necessary; the inla.mesh.create() should
    ## have removed all unused points. 2013-02-24 /FL
    ## points.dilation =
    ##    mesh.dilation$loc[unique(as.vector(mesh.dilation$graph$tv)),]

    z = (matrix(splancs::nndistF(mesh.dilation$loc, xy),
                resolution[1],resolution[2]))
    segm.closing =
        inla.contour.segment(ax[[1]],ax[[2]],z,
                             levels=c(concave),
                             positive=TRUE,
                             eps=eps)

    segm.closing$crs <- crs

    segm.closing
}
