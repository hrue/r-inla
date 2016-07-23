## Internal: update.crs
## Internal: inla.internal.sp2segment.join
## Internal: inla.identical.CRS
##
## S3methods; also export some methods explicitly
## Export: inla.sp2segment
## Export: as.inla.mesh.segment
## Export: as.inla.mesh.segment!Polygon
## Export: as.inla.mesh.segment!Polygons
## Export: as.inla.mesh.segment!SpatialPolygons
## Export: as.inla.mesh.segment!SpatialPolygonsDataFrame
## Export: as.inla.mesh.segment!Line
## Export: as.inla.mesh.segment!Lines
## Export: as.inla.mesh.segment!SpatialLines
## Export: as.inla.mesh.segment!SpatialLinesDataFrame
## Export: as.inla.mesh.segment!SpatialPoints
## Export: as.inla.mesh.segment!SpatialPointsDataFrame
## Export: inla.CRS
## Export: inla.spTransform
## Export: inla.spTransform!SpatialPoints
## Export: inla.spTransform!inla.mesh
## Export: inla.spTransform!inla.mesh.segment
## Export: inla.spTransform!default


inla.CRS <- function(projargs = NA_character_, doCheckCRSArgs = TRUE,
                     params=NULL, orient=NULL) {
  predef <- list(
    longlat = "+proj=longlat +ellps=sphere +a=1 +b=1",
    sphere = "+proj=geocent +ellps=sphere +a=1 +b=1 +units=m",
    mollweide = "+proj=moll +ellps=sphere +units=m +a=0.7071067811865476 +b=0.7071067811865476",
    lambert = "+proj=cea +ellps=sphere +lat_ts=0 +units=m +a=1 +b=1")
  if (projargs %in% names(predef)) {
    x <- CRS(predef[[projargs]], doCheckCRSArgs=doCheckCRSArgs)
  } else {
    x <- CRS(projargs, doCheckCRSArgs=doCheckCRSArgs)
  }
  if (!is.null(orient)) {
    x <- list(crs=x, orient=orient)
    class(x) <- "inla.CRS"
  }
  x
}


inla.crs.transform.orient <- function(x, orient, inverse) {
  warning("inla.crs.transform.orient NOT IMPLEMENTED")
  x
}


## CRS proj4 string for name=value pair list
inla.list2proj4string <- function(x) {
  paste(lapply(names(x),
               function(xx) {
    if (is.na(x[[xx]])) {
      paste("+", xx, sep="")
    } else {
      paste("+", xx, "=", x[[xx]], sep="")
    }
  }),
  collapse=" ")
}

## List of name=value pairs from CRS proj4 string
inla.proj4string2list <- function(x) {
  if (!is.character(x)) {
    stop("proj4string must be of class character")
  }
  do.call(c, lapply(strsplit(x=strsplit(x=paste(" ", x, sep=""),
                                        split=" \\+")[[1]][-1],
                             split="="),
                    function(x) {
               xx <- list(x[2])
               names(xx) <- x[1]
               xx
             }))
}


## +proj=longlat in (-180,360)x(-90,90)
## +proj=moll in (-2,2)x(-1,1) scaled by +a and +b, and +units
## +proj=lambert in (?,?)x(?,?) scaled by +a and +b, and +units
inla.spTransformBounds <- function(crs) {
  if (inherits(crs, "inla.CRS")) {
    args <- inla.proj4string2list(rgdal::CRSargs(crs$crs))
  } else {
    args <- inla.proj4string2list(rgdal::CRSargs(crs))
  }
  if (args[["proj"]] == "longlat") {
    bounds <- list(type="rectangle", xlim=c(-180,360), ylim=c(-90,90))
  } else if (args[["proj"]] == "moll") {
    axis <- c(2, 1)
    center <- c(0,0)
    if (!is.null(args[["a"]])) {
      axis[1] <- axis[1] * as.numeric(args$a) / sqrt(1/2)
    }
    if (!is.null(args[["b"]])) {
      axis[2] <- axis[2] * as.numeric(args$b) / sqrt(1/2)
    }
    bounds <- list(type="ellipse", axis=axis, center=c(0,0))
  } else {
    warning("Could not identify transformation shape.")
    bounds <- list(type="rectangle", xlim=c(-Inf, Inf), ylim=c(-Inf, Inf))
  }
  if (bounds$type == "rectangle") {
    bounds$polygon <- cbind(bounds$xlim[c(1,2,2,1,1)],
                            bounds$ylim[c(1,1,2,2,1)])
  } else if (bounds$type == "ellipse") {
    theta <- seq(0, 2*pi, length=1000)
    bounds$polygon <- cbind(bounds$center[1] + bounds$axis[1]*cos(theta),
                            bounds$center[2] + bounds$axis[2]*sin(theta))
  } else {
    stop("Unknown transformation type. This should not happen.")
  }
  bounds
}
## TRUE/FALSE for points inside/outside projection domain.
inla.spTransformBounds.ok <- function(x, bounds) {
  inla.require.inherits(x, "matrix")
  (sp::point.in.polygon(x[,1], x[,2],
                        bounds$polygon[,1], bounds$polygon[,2])
    > 0)
}



update.crs <- function(crs, newcrs, mismatch.allowed) {
  if (is.null(crs)) {
    newcrs
  } else {
    if (!mismatch.allowed && !identical(crs, newcrs)) {
      show(crs)
      show(newcrs)
      stop("CRS information mismatch.")
    }
    crs
  }
}


inla.identical.CRS <- function(crs0, crs1, crsonly=FALSE) {
  if (!crsonly) {
    stop("'crsonly=FALSE' not implemented.")
  } else {
    if (inherits(crs0, "inla.CRS")) {
      crs0 <- crs0$crs
    }
    if (inherits(crs1, "inla.CRS")) {
      crs1 <- crs1$crs
    }
    identical(crs0, crs1)
  }
}

## Low level transformation of raw coordinates.
inla.spTransform.default <- function(x, crs0, crs1, passthrough=FALSE, ...) {
  ok0 <- (!is.null(crs0) &&
          ((inherits(crs0, "CRS") && !is.na(rgdal::CRSargs(crs0))) ||
           (inherits(crs0, "inla.CRS"))))
  ok1 <- (!is.null(crs1) &&
          ((inherits(crs1, "CRS") && !is.na(rgdal::CRSargs(crs1))) ||
           (inherits(crs1, "inla.CRS"))))
  if (ok0 && ok1) {
    if (ncol(x) == 2) {
      x <- cbind(x, 0)
    }
    onsphere <- inla.identical.CRS(crs0, inla.CRS("sphere"), crsonly=TRUE)
    if (onsphere) {
      ok <- TRUE
    } else {
      bounds <- inla.spTransformBounds(crs0)
      ok <- inla.spTransformBounds.ok(x, bounds)
      if (!all(ok)) {
        xx <- x
      }
    }
    if (inherits(crs0, "inla.CRS")) {
      if (!onsphere) {
        x <- spTransform(SpatialPoints(x[ok,], proj4string=crs0$crs),
                         inla.CRS("sphere"))
      }
      if (!is.null(crs0$orient)) {
        x <- SpatialPoints(inla.crs.transform.orient(coordinates(x),
                                                     crs0$orient,
                                                     inverse=TRUE),
                           proj4string=inla.CRS("sphere"))
      }
      onshpere <- TRUE
    } else {
      x <- SpatialPoints(x[ok,], proj4string=crs0)
    }
    if (inherits(crs1, "inla.CRS")) {
      if (!onsphere) {
        x <- spTransform(x, inla.CRS("sphere"))
      }
      if (!is.null(crs1$orient)) {
        x <- SpatialPoints(inla.crs.transform.orient(coordinates(x),
                                                     crs1$orient,
                                                     inverse=FALSE),
                           proj4string=inla.CRS("sphere"))
      }
      x <- spTransform(x, crs1$crs)
    } else {
      x <- spTransform(x, crs1)
    }
    if (!all(ok)) {
      xx[ok,] <- coordinates(x)
      xx[!ok,] <- NA
      x <- xx
    }
  } else if (!passthrough) {
    if (!ok0) {
      stop("'crs0' is an invalid coordinate reference object.")
    }
    if (!ok1) {
      stop("'crs1' is an invalid coordinate reference object.")
    }
  }
  if (is.matrix(x)) {
    invisible(x)
  } else {
    invisible(coordinates(x))
  }
}

inla.spTransform.SpatialPoints <- function(x, CRSobj, passthrough=FALSE, ...) {
  ok0 <- !is.na(proj4string(x))
  ok1 <- (!missing(CRSobj) && is.null(CRSobj) &&
          (inherits(CRSobj, "CRS") && !is.na(rgdal::CRSargs(CRSobj))))
  if (ok0 && ok1) {
    invisible(spTransform(x, CRSobj=CRSobj))
  } else if (ok1) { ## Know: !ok0 && ok1
    if (!passthrough) {
      stop("Invalid origin CRS for SpatialPoints")
    }
    invisible(SpatialPoints(coordinates(x), proj4string=CRSobj))
  } else { ## Know: (ok0 || !ok0) && !ok1
    if (!passthrough) {
      stop("Invalid target CRS for SpatialPoints")
    }
    invisible(SpatialPoints(coordinates(x), proj4string=inla.CRS()))
  }
}

inla.spTransform.inla.mesh.segment <- function(x, CRSobj, passthrough=FALSE, ...) {
  x$loc <- inla.spTransform(x$loc, x$crs, CRSobj, passthrough=passthrough)
  x$crs <- CRSobj
  invisible(x)
}

inla.spTransform.inla.mesh <- function(x, CRSobj, passthrough=FALSE, ...) {
  x$loc <- inla.spTransform(x$loc, x$crs, CRSobj, passthrough=passthrough)
  x$crs <- CRSobj
  invisible(x)
}

inla.spTransform <- function(x, ...) {
  UseMethod("inla.spTransform")
}


## Input: list of segments, all closed polygons.
inla.internal.sp2segment.join <- function(inp, grp=NULL, closed=TRUE) {
    crs <- NULL
    if (length(inp) > 0) {
      out.loc = matrix(0,0,ncol(inp[[1]]$loc))
      for (k in seq_along(inp)) {
        crs <- update.crs(crs, inp[[k]]$loc, mismatch.allowed=FALSE)
      }
    } else {
        out.loc = matrix(0,0,2)
    }
    out.idx = matrix(0L,0,2)
    if (is.null(grp)) {
        out.grp = NULL
    } else {
        out.grp = integer(0)
    }
    for (k in seq_along(inp)) {
        inp.loc = inp[[k]]$loc
        inp.idx = inp[[k]]$idx
        inp.grp = inp[[k]]$grp
        offset = nrow(out.loc)
        n = nrow(as.matrix(inp.idx))
        if (closed) {
            if (!is.null(grp) && is.null(inp.grp)) {
                inp.grp = rep(grp[k], n)
            }
            if (ncol(as.matrix(inp.idx))==1) {
                inp.idx = cbind(inp.idx, inp.idx[c(2:n,1)])
            }
        } else {
            if (!is.null(grp) && is.null(inp.grp)) {
                inp.grp = rep(grp[k], n-1)
            }
            if (ncol(as.matrix(inp.idx))==1) {
                inp.idx = cbind(inp.idx[-n], inp.idx[-1])
            }
        }
        out.loc = rbind(out.loc, inp.loc)
        out.idx = rbind(out.idx, inp.idx+offset)
        if (!is.null(grp)) {
            out.grp = c(out.grp, inp.grp)
        }
    }
    inla.mesh.segment(loc=out.loc, idx=out.idx, grp=out.grp, is.bnd=FALSE,
                      crs=crs)
}


as.inla.mesh.segment <-
    function(sp, ...)
{
    UseMethod("as.inla.mesh.segment")
}

inla.sp2segment <-
    function(sp, ...)
{
    UseMethod("as.inla.mesh.segment")
}



as.inla.mesh.segment.SpatialPoints <-
    function (sp, reverse=FALSE, grp = NULL, is.bnd=TRUE, ...)
  {
    crs <- CRS(proj4string(sp))
    loc <- coordinates(sp)

    n = dim(loc)[1L]
    if (reverse) {
        idx <- seq(n, 1L, length=n)
    } else {
        idx <- seq_len(n)
    }
    inla.mesh.segment(loc = loc, idx = idx, grp = grp, is.bnd = is.bnd,
                      crs=crs)
  }

as.inla.mesh.segment.SpatialPointsDataFrame <-
    function (sp, ...)
{
    as.inla.mesh.segment.SpatialLines(sp, ...)
}



as.inla.mesh.segment.Line <-
    function(sp, reverse=FALSE, crs=NULL, ...)
{
    loc = sp@coords
    n = dim(loc)[1L]
    if (reverse) {
        idx <- seq(n, 1L, length=n)
    } else {
        idx <- seq_len(n)
    }
    inla.mesh.segment(loc = loc, idx = idx, is.bnd = FALSE, crs=crs)
}

as.inla.mesh.segment.Lines <-
    function (sp, join = TRUE, crs=NULL, ...)
{
    segm <- as.list(lapply(sp@Lines,
                           function(x) as.inla.mesh.segment(x, crs=crs, ...)))
    if (join)
        segm = inla.internal.sp2segment.join(segm, grp = NULL, closed=FALSE)
    segm
}

as.inla.mesh.segment.SpatialLines <-
    function (sp, join = TRUE, grp = NULL, ...)
  {
    crs <- CRS(proj4string(sp))
    segm = list()
    for (k in 1:length(sp@lines))
      segm[[k]] = as.inla.mesh.segment(sp@lines[[k]], join = TRUE,
                                       crs = crs, ...)
    if (join) {
        if (missing(grp)) {
            grp = 1:length(segm)
        }
        segm = inla.internal.sp2segment.join(segm, grp = grp, closed=FALSE)
    }
    segm
}

as.inla.mesh.segment.SpatialLinesDataFrame <-
    function (sp, ...)
{
    as.inla.mesh.segment.SpatialLines(sp, ...)
}

as.inla.mesh.segment.SpatialPolygons <-
    function(sp, join=TRUE, grp=NULL, ...)
{
    crs <- CRS(proj4string(sp))
    segm = list()
    for (k in 1:length(sp@polygons))
        segm[[k]] = as.inla.mesh.segment(sp@polygons[[k]], join=TRUE, crs=crs)
    if (join) {
        if (missing(grp)) {
            grp = 1:length(segm)
        }
        segm = inla.internal.sp2segment.join(segm, grp=grp)
    }
    segm
}

as.inla.mesh.segment.SpatialPolygonsDataFrame <-
    function(sp, ...)
{
    as.inla.mesh.segment.SpatialPolygons(sp, ...)
}

as.inla.mesh.segment.Polygons <-
    function(sp, join=TRUE, crs=NULL, ...)
{
    crs <- CRS(proj4string(sp))
    segm = as.list(lapply(sp@Polygons,
                          function (x) as.inla.mesh.segment(x, crs=crs)))
    if (join)
        segm = inla.internal.sp2segment.join(segm, grp=NULL)
    segm
}

as.inla.mesh.segment.Polygon <-
    function(sp, crs=NULL, ...)
{
    loc = sp@coords[-dim(sp@coords)[1L],,drop=FALSE]
    n = dim(loc)[1L]
    if (sp@hole)
        if (sp@ringDir==1)
            idx = c(1L:n,1L)
        else
            idx = c(1L,seq(n,1L,length.out=n))
    else
        if (sp@ringDir==1)
            idx = c(1L,seq(n,1L,length.out=n))
        else
            idx = c(1L:n,1L)
    inla.mesh.segment(loc=loc, idx=idx, is.bnd=TRUE, crs=crs)
}
