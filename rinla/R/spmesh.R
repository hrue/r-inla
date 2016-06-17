## Internal: update.crs
## Internal: inla.internal.sp2segment.join
## Internal: safe.spTransform
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
## Internal: spTransform.inla.mesh
## Internal: spTransform.inla.mesh.segment


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


safe.spTransform <- function(x, CRSobj, ...) {
  if (!is.null(CRSobj) &&
      !is.na(CRSargs(CRSobj)) &&
      !is.na(proj4string(x))) {
    spTransform(x, CRSobj)
  } else {
    x
  }
}

spTransform.inla.mesh.segment <- function(x, CRSobj, ...) {
  if (!is.null(CRSobj) &&
      !is.na(CRSargs(CRSobj)) &&
      !is.null(x$crs) &&
      !is.na(CRSargs(x$crs))) {
    x$loc <-
      coordinates(spTransform(
        SpatialPoints(x$loc, proj4string=x$crs),
        CRSobj))
    x$crs <- CRSobj
  }
  invisible(x)
}

spTransform.inla.mesh <- function(x, CRSobj, ...) {
  if (!is.null(CRSobj) &&
      !is.na(CRSargs(CRSobj)) &&
      !is.null(x$crs) &&
      !is.na(CRSargs(x$crs))) {
    x$loc <-
      coordinates(spTransform(
        SpatialPoints(x$loc, proj4string=x$crs),
        CRSobj))
    x$crs <- CRSobj
  }
  invisible(x)
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
