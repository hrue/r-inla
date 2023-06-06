## Internal: internal.update.crs
## Internal: inla.internal.sp2segment.join
## Internal: inla.crs.transform.oblique
## Internal: inla.rotmat3213
## Internal: inla.rotmat3123
## Internal: inla.crs.graticule
## Internal: inla.crs.tissot
## Internal: internal.clip
## Internal: inla.crs.bounds
## Internal: inla.crs.bounds.check




#' PROJ6 detection
#'
#' Detect whether PROJ6 is available for INLA. Deprecated and always returns `TRUE`..
#'
#' @aliases inla.has_PROJ6 inla.not_for_PROJ6 inla.not_for_PROJ4
#' inla.fallback_PROJ6 inla.requires_PROJ6
#' @details `inla.has_PROJ6` is called to check if PROJ6&GDAL3 are available.
#' @return For `inla.has_PROJ6`, always returns `TRUE`. Previously: logical; `TRUE` if PROJ6 is available,
#' `FALSE` otherwise
#' @examples
#'
#' inla.has_PROJ6()
#' @export inla.has_PROJ6
#' @rdname inla.has_PROJ6
inla.has_PROJ6 <- function() {
    lifecycle::deprecate_warn(
        "23.01.01",
        "inla.has_PROJ6()",
        details = "PROJ4 is no longer supported, and inla.has_PROJ6() always returns TRUE"
    )
    TRUE
}

#' @describeIn inla.has_PROJ6 Called to warn about using old PROJ4
#' features even though PROJ6 is available
#' @param fun The name of the calling function

inla.not_for_PROJ6 <- function(fun) {
    if (inla.has_PROJ6()) {
        stack <- sys.calls()
        stack <- lapply(as.list(stack), function(x) as.character(deparse(x)))
        warning(paste0(
            "'",
            fun,
            "' should not be used with PROJ6\n",
            "Call stack for developer debugging:\n",
            paste0(lapply(stack, function(x) substr(x, 1, 70)),
                collapse = "\n"
            )
        ))
    }
}

#' @describeIn inla.has_PROJ6 Called to give an error when
#' calling methods that are only available for PROJ6

inla.not_for_PROJ4 <- function(fun) {
    if (!inla.has_PROJ6()) {
        stop(paste0(
            "'",
            fun,
            "' is not supported for PROJ4"
        ))
    }
}

#' @describeIn inla.has_PROJ6 Called to warn about falling back
#' to using old PROJ4 methods when a RPOJ6 method hasn't been implemented

inla.fallback_PROJ6 <- function(fun) {
    if (inla.has_PROJ6()) {
        stop(paste0(
            "'",
            fun,
            "' method for PROJ6 not implemented. PROJ4 is no longer supported."
        ))
    }
}

#' @describeIn inla.has_PROJ6 Called to give an error when PROJ6
#' is required but not available

inla.requires_PROJ6 <- function(fun) {
    if (!inla.has_PROJ6()) {
        stop(paste0(
            "'",
            fun,
            "' requires PROJ6"
        ))
    }
}






#' @title Extract CRS information
#'
#' @description
#' Wrapper for CRS(projargs) (PROJ4) and CRS(wkt) for `sp::Spatial`
#' objects.
#'
#' This function is a convenience method to workaround PROJ4/PROJ6 differences,
#' and the lack of a crs extraction method for Spatial objects.
#'
#' @param x A `sp::Spatial` object
#' @return A `CRS` object, or NULL if no valid CRS identified
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @examples
#'
#' \dontrun{
#' if (interactive()) {
#'     s <- sp::SpatialPoints(matrix(1:6, 3, 2), proj4string = inla.CRS("sphere"))
#'     inla.sp_get_crs(s)
#' }
#' }
#'
#' @export inla.sp_get_crs
inla.sp_get_crs <- function(x) {
    if (is.null(x)) {
        return(NULL)
    }
    suppressWarnings(crs <- sp::CRS(SRS_string = sp::wkt(x)))
    crs
}





internal.clip <- function(bounds, coords, eps = 0.05) {
    ## Clip 2D coordinate matrix of polylines and generate a list of Line objects
    ## bounds is from inla.crs.bounds
    ## This implementation only removes "long" line segments.
    thelines <- list()
    ## Rudimentary cutting:
    toolong <-
        which(c(
            TRUE,
            (diff(coords[, 1]) / diff(bounds$xlim))^2
                + (diff(coords[, 2]) / diff(bounds$ylim))^2
            > eps^2,
            TRUE
        ))
    start <- toolong[-length(toolong)]
    ending <- toolong[-1] - 1
    for (i in seq_along(start)) {
        if (start[i] < ending[i]) {
            thelines <-
                c(
                    thelines,
                    list(sp::Line(coords[start[i]:ending[i], 1:2]))
                )
        }
    }
    thelines
}




#' Handling CRS/WKT
#'
#' Get and set CRS object or WKT string properties.
#'
#'
#' @aliases inla.wkt_is_geocent inla.crs_is_geocent
#' inla.wkt_get_ellipsoid_radius inla.crs_get_ellipsoid_radius
#' inla.wkt_set_ellipsoid_radius inla.crs_set_ellipsoid_radius
#' inla.wkt_unit_params crs_wkt inla.wkt_get_lengthunit inla.wkt_set_lengthunit
#' inla.crs_get_wkt inla.crs_get_lengthunit inla.crs_set_lengthunit
#' @param wkt A WKT2 character string
#' @param crs A `sp::CRS` or `inla.CRS` object
#' @param unit character, name of a unit. Supported names are "metre",
#' "kilometre", and the aliases "meter", "m", International metre",
#' "kilometer", and "km", as defined by `inla.wkt_unit_params` or the
#' `params` argument. (For legacy PROJ4 use, only "m" and "km" are
#' supported)
#' @param params Length unit definitions, in the list format produced by
#' `inla.wkt_unit_params()`, Default: NULL, which invokes
#' `inla.wkt_unit_params()`
#' @return For `inla.wkt_unit_params`, a list of named unit definitions
#'
#' For `inla.wkt_get_lengthunit`, a list of length units used in the wkt
#' string, excluding the ellipsoid radius unit.
#'
#' For `inla.wkt_set_lengthunit`, a WKT2 string with altered length units.
#' Note that the length unit for the ellipsoid radius is unchanged.
#'
#' For `inla.crs_get_wkt`, WKT2 string.
#'
#' For `inla.crs_get_lengthunit`, a list of length units used in the wkt
#' string, excluding the ellipsoid radius unit. (For legacy PROJ4 code, the raw
#' units from the proj4string are returned, if present.)
#'
#' For `inla.crs_set_lengthunit`, a `sp::CRS` object with altered
#' length units. Note that the length unit for the ellipsoid radius is
#' unchanged.
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @seealso [inla.sp_get_crs()]
#' @examples
#'
#' \dontrun{
#'   c1 <- inla.CRS("globe")
#'   inla.crs_get_lengthunit(c1)
#'   c2 <- inla.crs_set_lengthunit(c1, "km")
#'   inla.crs_get_lengthunit(c2)
#' }
#'
#' @rdname crs_wkt
#' @export inla.wkt_is_geocent
inla.wkt_is_geocent <- function(wkt) {
    if (is.null(wkt) || identical(wkt, "")) {
        return(FALSE)
    }
    # See https://proceedings.esri.com/library/userconf/proc17/tech-workshops/tw_2588-212.pdf
    geo_crs_items <- c(
        "GEODCRS", "GEOGCRS",
        "BASEGEODCRS", "BASEGEOGCRS"
    )
    wt <- inla.as.wkt_tree.wkt(wkt)
    if (identical(wt[["label"]], "BOUNDCRS")) {
        wt <- inla.wkt_tree_get_item(wt, "SOURCECRS")
        wt <- inla.wkt_tree_get_item(wt, c("PROJCRS", geo_crs_items))
    }
    if (identical(wt[["label"]], "PROJCRS")) {
        wt <- inla.wkt_tree_get_item(wt, geo_crs_items)
    }
    if (!(wt[["label"]] %in% geo_crs_items)) {
        return(FALSE)
    }
    cs <- inla.wkt_tree_get_item(wt, "CS")
    if (is.null(cs)) {
        return(FALSE)
    }
    cart <- ((cs[["params"]][[1]] == "Cartesian") &&
        (cs[["params"]][[2]] == "3"))
    if (!cart) {
        return(FALSE)
    }
    axis_names <- c('"(X)"', '"(Y)"', '"(Z)"')
    axis_types <- c("geocentricX", "geocentricY", "geocentricZ")
    for (k in seq_len(3)) {
        axis <- inla.wkt_tree_get_item(wt, "AXIS", k)
        if (!((axis[["params"]][[1]] == axis_names[k]) &&
            (axis[["params"]][[2]] == axis_types[k]))) {
            return(FALSE)
        }
    }
    TRUE
}

#' @export
#' @rdname crs_wkt

inla.crs_is_geocent <- function(crs) {
    wkt <- inla.crs_get_wkt(crs)
    result <- inla.wkt_is_geocent(wkt)
    result
}


#' @rdname crs_wkt
#' @export

inla.wkt_get_ellipsoid_radius <- function(wkt) {
    geo_crs_items <- c(
        "GEODCRS", "GEOGCRS",
        "BASEGEODCRS", "BASEGEOGCRS"
    )
    wt <- inla.as.wkt_tree.wkt(wkt)

    if (identical(wt[["label"]], "BOUNDCRS")) {
        wt <- inla.wkt_tree_get_item(wt, "SOURCECRS")
        wt <- inla.wkt_tree_get_item(wt, c("PROJCRS", geo_crs_items))
    }
    if (identical(wt[["label"]], "PROJCRS")) {
        wt <- inla.wkt_tree_get_item(wt, geo_crs_items)
    }
    if (is.null(wt) || !(wt[["label"]] %in% geo_crs_items)) {
        stop("Ellipsoid settings not found")
    }

    datum <- inla.wkt_tree_get_item(wt, c("DATUM", "ENSEMBLE"))
    if (is.null(datum)) {
        stop("Ellipsoid settings not found")
    }
    ellipsoid <- inla.wkt_tree_get_item(datum, "ELLIPSOID")
    if (is.null(ellipsoid)) {
        stop("Ellipsoid settings not found")
    }
    as.numeric(ellipsoid[["params"]][[2]])
}

#' @rdname crs_wkt
#' @export

inla.crs_get_ellipsoid_radius <- function(crs) {
    inla.wkt_get_ellipsoid_radius(inla.crs_get_wkt(crs))
}


#' @rdname crs_wkt
#' @param radius numeric
#' @export

inla.wkt_set_ellipsoid_radius <- function(wkt, radius) {
    geo_crs_items <- c(
        "GEODCRS", "GEOGCRS",
        "BASEGEODCRS", "BASEGEOGCRS"
    )

    set_radius <- function(wt) {
        if (is.null(wt)) {
            stop("Ellipsoid settings not found")
        } else if (wt[["label"]] %in% geo_crs_items) {
            datum <- inla.wkt_tree_get_item(wt, c("DATUM", "ENSEMBLE"))
            null_datum <- is.null(datum)
            if (null_datum) {
                stop("Ellipsoid settings not found")
            }
            ellipsoid <- inla.wkt_tree_get_item(datum, "ELLIPSOID")
            if (is.null(ellipsoid)) {
                stop("Ellipsoid settings not found")
            }
            ellipsoid[["params"]][[2]] <- as.character(radius)
            datum <- inla.wkt_tree_set_item(datum, ellipsoid)
            wt <- inla.wkt_tree_set_item(wt, datum)
        } else if (wt[["label"]] %in% c("BOUNDCRS", "SOURCECRS", "PROJCRS")) {
            wt_sub <- inla.wkt_tree_get_item(
                wt,
                c(
                    "BOUNDCRS", "SOURCECRS", "PROJCRS",
                    geo_crs_items
                )
            )
            if (is.null(wt_sub)) {
                stop("Ellipsoid settings not found")
            }
            wt_sub_new <- set_radius(wt_sub)
            wt <- inla.wkt_tree_set_item(wt, wt_sub_new)
        } else {
            stop("Ellipsoid settings not found")
        }
        wt
    }

    wt <- inla.as.wkt_tree.wkt(wkt)
    wt <- set_radius(wt)
    inla.as.wkt.wkt_tree(wt)
}

#' @rdname crs_wkt
#' @export

inla.crs_set_ellipsoid_radius <- function(crs, radius) {
    wkt <- inla.crs_get_wkt(crs)
    wkt <- inla.wkt_set_ellipsoid_radius(wkt, radius)
    new_crs <- inla.CRS(SRS_string = wkt)
    if (inherits(crs, "inla.CRS")) {
        crs$crs <- new_crs
        crs
    } else {
        new_crs
    }
}





inla.crs.graticule <- function(x, by = c(15, 15, 45), add = FALSE, do.plot = TRUE,
                               eps = 0.05, ...) {
    inla.fallback_PROJ6("inla.crs.graticule")

    ## Graticule
    if (is.null(by)) {
        return(invisible(list()))
    }
    if (length(by) < 2) {
        by <- by[c(1, 1, 1)]
    } else if (length(by) < 3) {
        by <- by[c(1, 2, 1)]
    }
    n <- c(floor(180 / by[1]), ceiling(90 / by[2]) - 1, floor(180 / by[3]))
    bounds <- inla.crs.bounds(x)
    if (by[1] > 0) {
        special.poles <- (by[1] != by[3]) && (by[2] > 0)
        lon <- ((1 - n[1]):n[1]) * by[1]
        if (special.poles) {
            lat <- seq(-n[2] * by[2], n[2] * by[2], length = 91)
        } else {
            lat <- seq(-90 + 1e-6, 90 - 1e-6, length = 91)
        }
        meridians <- as.matrix(expand.grid(lat, lon)[, 2:1])
        proj.mer.coords <- inla.spTransform(meridians, inla.CRS("longlat"), x)
        proj.mer.coords1 <- matrix(
            proj.mer.coords[, 1], length(lat),
            length(lon)
        )
        proj.mer.coords2 <- matrix(
            proj.mer.coords[, 2], length(lat),
            length(lon)
        )

        mer.coords <-
            unlist(lapply(
                seq_along(lon),
                function(k) {
                    internal.clip(bounds, cbind(
                        proj.mer.coords1[, k, drop = FALSE],
                        proj.mer.coords2[, k, drop = FALSE]
                    ),
                    eps = eps
                    )
                }
            ),
            recursive = FALSE
            )

        if (special.poles) {
            if (by[3] > 0) {
                lon <- ((1 - n[3]):n[3]) * by[3]
                lat <- seq(-90 + 1e-6, -n[2] * by[2], length = ceiling((90 - n[2] * by[2]) / 2) + 1)
                meridians <- as.matrix(expand.grid(lat, lon)[, 2:1])
                proj.mer.coords <- inla.spTransform(meridians, inla.CRS("longlat"), x)
                proj.mer.coords1 <- matrix(
                    proj.mer.coords[, 1], length(lat),
                    length(lon)
                )
                proj.mer.coords2 <- matrix(
                    proj.mer.coords[, 2], length(lat),
                    length(lon)
                )
                mer.coords <-
                    c(
                        mer.coords,
                        unlist(lapply(
                            seq_along(lon),
                            function(k) {
                                internal.clip(bounds, cbind(
                                    proj.mer.coords1[, k, drop = FALSE],
                                    proj.mer.coords2[, k, drop = FALSE]
                                ),
                                eps = eps
                                )
                            }
                        ),
                        recursive = FALSE
                        )
                    )

                lat <- seq(n[2] * by[2], 90 - 1e-6, length = ceiling((90 - n[2] * by[2]) / 2) + 1)
                meridians <- as.matrix(expand.grid(lat, lon)[, 2:1])
                proj.mer.coords <- inla.spTransform(meridians, inla.CRS("longlat"), x)
                proj.mer.coords1 <- matrix(
                    proj.mer.coords[, 1], length(lat),
                    length(lon)
                )
                proj.mer.coords2 <- matrix(
                    proj.mer.coords[, 2], length(lat),
                    length(lon)
                )
                mer.coords <-
                    c(
                        mer.coords,
                        unlist(lapply(
                            seq_along(lon),
                            function(k) {
                                internal.clip(bounds, cbind(
                                    proj.mer.coords1[, k, drop = FALSE],
                                    proj.mer.coords2[, k, drop = FALSE]
                                ),
                                eps = eps
                                )
                            }
                        ),
                        recursive = FALSE
                        )
                    )
            }
        }

        proj.mer <-
            sp::SpatialLines(list(sp::Lines(mer.coords, ID = "meridians")),
                proj4string = CRS(inla.CRSargs(x))
            )
        if (do.plot) {
            args <- list(x = proj.mer, ...)
            args <- args[intersect(
                names(args),
                union(
                    names(formals(plot.default)),
                    union(
                        names(formals(sp:::plot.Spatial)),
                        names(formals(sp:::plotSpatialLines))
                    )
                )
            )]
            if (add) {
                do.call(plot, c(list(add = TRUE), args))
            } else {
                do.call(plot, args)
                add <- TRUE
            }
        }
    } else {
        proj.mer <- NULL
    }
    if (by[2] > 0) {
        lon <- seq(-180 + 1e-6, 180 - 1e-6, length = 181)
        lat <- ((-n[2]):n[2]) * by[2]
        parallels <- as.matrix(expand.grid(lon, lat))
        proj.par.coords <- inla.spTransform(parallels, inla.CRS("longlat"), x)
        proj.par.coords1 <- matrix(
            proj.par.coords[, 1], length(lon),
            length(lat)
        )
        proj.par.coords2 <- matrix(
            proj.par.coords[, 2], length(lon),
            length(lat)
        )
        proj.par <-
            sp::SpatialLines(list(sp::Lines(
                unlist(lapply(
                    seq_along(lat),
                    function(k) {
                        internal.clip(bounds, cbind(
                            proj.par.coords1[, k],
                            proj.par.coords2[, k]
                        ), eps = eps)
                    }
                ),
                recursive = FALSE
                ),
                ID = "parallels"
            )),
            proj4string = CRS(inla.CRSargs(x))
            )
        if (do.plot) {
            args <- list(x = proj.par, ...)
            args <- args[intersect(
                names(args),
                union(
                    names(formals(plot.default)),
                    union(
                        names(formals(sp:::plot.Spatial)),
                        names(formals(sp:::plotSpatialLines))
                    )
                )
            )]
            if (add) {
                do.call(plot, c(list(add = TRUE), args))
            } else {
                do.call(plot, args)
            }
        }
    } else {
        proj.par <- NULL
    }
    invisible(list(meridians = proj.mer, parallels = proj.par))
}

inla.crs.tissot <- function(x, by = c(30, 30, 30), add = FALSE, do.plot = TRUE,
                            eps = 0.05, diff.eps = 1e-2, ...) {
    inla.fallback_PROJ6("inla.crs.graticule")

    if (is.null(by)) {
        return(invisible(list()))
    }
    if (length(by) < 2) {
        by <- c(by[1], by[1], 30)
    } else if (length(by) < 3) {
        by <- c(by[1:2], 30)
    }
    bounds <- inla.crs.bounds(x)
    n <- c(floor(180 / by[1]), ceiling(90 / by[2]) - 1)
    lon <- ((1 - n[1]):n[1]) * by[1]
    lat <- ((-n[2]):n[2]) * by[2]
    loc0.lon <- loc0.lat <- loc0 <-
        cbind(as.matrix(expand.grid(lat, lon)[, 2:1]), 0)
    loc0.lon[, 1] <- loc0.lon[, 1] + diff.eps / cos(loc0.lat[, 2] * pi / 180)
    loc0.lat[, 2] <- loc0.lat[, 2] + diff.eps
    crs.longlat <- inla.CRS("longlat")

    loc1 <- inla.spTransform(loc0, crs.longlat, x)
    loc1.lon <- inla.spTransform(loc0.lon, crs.longlat, x)
    loc1.lat <- inla.spTransform(loc0.lat, crs.longlat, x)
    ok <- (rowSums(is.na(loc1)) +
        rowSums(is.na(loc1.lon)) +
        rowSums(is.na(loc1.lat)) == 0)
    loc1 <- loc1[ok, , drop = FALSE]
    loc1.lon <- loc1.lon[ok, , drop = FALSE]
    loc1.lat <- loc1.lat[ok, , drop = FALSE]

    diff.lon <- (loc1.lon - loc1) / eps
    diff.lat <- (loc1.lat - loc1) / eps

    scale <- by[3]
    theta <- seq(0, 2 * pi, length = 181)
    ct <- cos(theta) * scale
    st <- sin(theta) * scale

    collection <-
        sp::SpatialLines(list(sp::Lines(
            unlist(lapply(
                seq_len(nrow(loc1)),
                function(k) {
                    loc1.ellipse <-
                        cbind(
                            loc1[k, 1] + diff.lon[k, 1] * ct + diff.lat[k, 1] * st,
                            loc1[k, 2] + diff.lon[k, 2] * ct + diff.lat[k, 2] * st
                        )
                    internal.clip(bounds, loc1.ellipse, eps = eps)
                }
            ),
            recursive = FALSE
            ),
            ID = "parallels"
        )),
        proj4string = CRS(inla.CRSargs(x))
        )
    if (do.plot) {
        args <- list(x = collection, ...)
        args <- args[intersect(
            names(args),
            union(
                names(formals(plot.default)),
                union(
                    names(formals(sp:::plot.Spatial)),
                    names(formals(sp:::plotSpatialLines))
                )
            )
        )]
        if (add) {
            do.call(plot, c(list(add = TRUE), args))
        } else {
            do.call(plot, args)
        }
    }

    invisible(list(tissot = collection))
}



#' Plot CRS and inla.CRS objects
#'
#' Plot the outline of a CRS or inla.CRS projection, with optional graticules
#' (transformed parallels and meridians) and Tissot indicatrices.
#'
#'
#' @aliases plot.inla.CRS plot.CRS
#' @param x A `CRS` or [inla.CRS()] object.
#' @param xlim Optional x-axis limits.
#' @param ylim Optional y-axis limits.
#' @param outline Logical, if `TRUE`, draw the outline of the projection.
#' @param graticule Vector of length at most 3, to plot meridians with spacing
#' `graticule[1]` degrees and parallels with spacing `graticule[2]`
#' degrees. `graticule[3]` optionally specifies the spacing above and
#' below the first and last parallel.  When `graticule[1]==0` no meridians
#' are drawn, and when `graticule[2]==0` no parallels are drawn. Use
#' `graticule=NULL` to skip drawing a graticule.
#' @param tissot Vector of length at most 3, to plot Tissot's indicatrices with
#' spacing `tissot[1]` degrees and parallels with spacing `tissot[2]`
#' degrees. `tissot[3]` specifices a scaling factor.  Use
#' `tissot=NULL` to skip drawing a Tissot's indicatrices.
#' @param asp The aspect ratio for the plot, default 1.
#' @param add If `TRUE`, add the projecton plot to an existing plot.
#' @param eps Clipping tolerance for rudimentary boundary clipping
#' @param \dots Additional arguments passed on to the internal calls to
#' `plot` and `lines`.
#' @author Finn Lindgren <finn.lindgren@@gmail.com>
#' @seealso [inla.CRS()]
#' @examples
#'
#' if (require("sf") && require("sp")) {
#'     oblique <- c(0, 45, 45, 0)
#'     for (projtype in c("longlat", "lambert", "mollweide", "hammer")) {
#'         plot(inla.CRS(projtype), main = projtype)
#'         plot(inla.CRS(projtype, oblique = oblique), main = paste("oblique", projtype))
#'     }
#' }
#' @method plot inla.CRS
#' @export
plot.inla.CRS <- function(x, xlim = NULL, ylim = NULL,
                          outline = TRUE,
                          graticule = c(15, 15, 45),
                          tissot = c(30, 30, 30),
                          asp = 1,
                          add = FALSE,
                          eps = 0.05,
                          ...) {
    inla.fallback_PROJ6("plot.inla.CRS")

    bounds <- inla.crs.bounds(x)
    if (is.null(xlim)) xlim <- bounds$xlim
    if (is.null(ylim)) ylim <- bounds$ylim
    if (!add) {
        args <- list(x = NA, type = "n", xlim = xlim, ylim = ylim, asp = asp, ...)
        args <- args[intersect(names(args), names(formals(plot.default)))]
        do.call(plot, args)
    }
    ## Outline
    if (outline) {
        args <- list(x = bounds$polygon, ...)
        args <- args[intersect(names(args), union(
            names(formals(lines.default)),
            names(formals(plot.xy))
        ))]
        do.call(lines, args)
    }
    ## Graticule
    inla.crs.graticule(x, by = graticule, add = TRUE, do.plot = TRUE, eps = eps, ...)
    ## Tissot
    inla.crs.tissot(x, by = tissot, add = TRUE, do.plot = TRUE, eps = eps, ...)
    invisible(NULL)
}

#' @export
#' @method plot CRS
#' @rdname plot.inla.CRS
plot.CRS <- function(x, xlim = NULL, ylim = NULL,
                     outline = TRUE,
                     graticule = c(15, 15, 45),
                     tissot = c(30, 30, 30),
                     asp = 1,
                     add = FALSE,
                     eps = 0.05,
                     ...) {
    invisible(plot.inla.CRS(x,
        xlim = xlim, ylim = ylim,
        outline = outline, graticule = graticule, tissot = tissot,
        asp = asp, add = add, ...
    ))
}




#' Create a coordinate reference system object
#'
#' Creates either a CRS object or an inla.CRS object, describing a coordinate
#' reference system
#'
#' The first two elements of the `oblique` vector are the (longitude,
#' latitude) coordinates for the oblique centre point. The third value
#' (orientation) is a counterclockwise rotation angle for an observer looking
#' at the centre point from outside the sphere. The fourth value is the
#' quasi-longitude (orbit angle) for a rotation along the oblique observers
#' equator.
#'
#' Simple oblique: `oblique=c(0, 45)`
#'
#' Polar: `oblique=c(0, 90)`
#'
#' Quasi-transversal: `oblique=c(0, 0, 90)`
#'
#' Satellite orbit viewpoint: `oblique=c(lon0-time*v1, 0, orbitangle,
#' orbit0+time*v2)`, where `lon0` is the longitude at which a satellite
#' orbit crosses the equator at `time=0`, when the satellite is at an
#' angle `orbit0` further along in its orbit.  The orbital angle relative
#' to the equatorial plane is `orbitangle`, and `v1` and `v2`
#' are the angular velocities of the planet and the satellite, respectively.
#' Note that "forward" from the satellite's point of view is "to the right" in
#' the projection.
#'
#' When `oblique[2]` or `oblique[3]` are non-zero, the resulting
#' projection is only correct for perfect spheres.
#'
#' @aliases inla.CRS inla.wkt_predef
#' @param projargs Either 1) a projection argument string suitable as input to
#' `sp::CRS`, or 2) an existing `CRS` object, or 3) a shortcut
#' reference string to a predefined projection; run
#' `names(inla.wkt_predef())` for valid predefined projections.
#' @param doCheckCRSArgs ignored.
#' @param args An optional list of name/value pairs to add to and/or override
#' the PROJ4 arguments in `projargs`.  `name=value` is converted to
#' `"+name=value"`, and `name=NA` is converted to `"+name"`.
#' @param oblique Vector of length at most 4 of rotation angles (in degrees)
#' for an oblique projection, all values defaulting to zero. The values
#' indicate (longitude, latitude, orientation, orbit), as explained in the
#' Details section below.
#' @param SRS_string a WKT2 string defining the coordinate system; see
#' `sp::CRS`. This takes precedence over `projargs`.
#' @param \dots Additional parameters. Not currently in use.
#' @return Either an `sp::CRS` object or an `inla.CRS` object,
#' depending on if the coordinate reference system described by the parameters
#' can be expressed with a pure `sp::CRS` object or not.
#'
#' An S3 `inla.CRS` object is a list, usually (but not necessarily)
#' containing at least one element: \item{crs }{The basic `sp::CRS`
#' object}
#'
#' `inla.wkt_predef` returns a WKT2 string defining a projection
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @seealso [sp::CRS()], [crs_wkt()],
#' [inla.sp_get_crs()] [plot.CRS()],
#' [inla.identical.CRS()]
#' @examples
#'
#'
#' if (require("sf") && require("sp")) {
#'     crs1 <- inla.CRS("longlat_globe")
#'     crs2 <- inla.CRS("lambert_globe")
#'     crs3 <- inla.CRS("mollweide_norm")
#'     crs4 <- inla.CRS("hammer_globe")
#'     crs5 <- inla.CRS("sphere")
#'     crs6 <- inla.CRS("globe")
#' }
#' \dontrun{
#' names(inla.wkt_predef())
#' }
#'
#' @export inla.CRS
inla.CRS <- function(projargs = NULL, doCheckCRSArgs = FALSE,
                     args = NULL, oblique = NULL,
                     SRS_string = NULL,
                     ...) {
    if (identical(projargs, "")) {
        projargs <- NULL
    }
    if (is.null(SRS_string) &&
        !is.null(projargs) &&
        !is.na(projargs) &&
        is.character(projargs)) {
        if (projargs %in% c("hammer", "lambert", "longlat", "mollweide")) {
            warning(
                paste0(
                    "Use of old predefined projection '",
                    projargs,
                    "' is deprecated. Converting to '",
                    projargs,
                    "_norm'"
                )
            )
            projargs <- paste0(projargs, "_norm")
        }
        predef <- inla.wkt_predef()
        if (projargs %in% names(predef)) {
            SRS_string <- predef[[projargs]]
            projargs <- NULL
        } else {
            warning(
                "'inla.CRS' should be given a SRS_string for PROJ6 or a known keyword for a predefined string given in projargs. Using fallback PROJ4 method."
            )
            if (!is.null(args)) {
                x <- CRS(projargs)
                if (typeof(args) != "list") {
                    stop("'args' must be NULL or a list of name=value pairs.")
                }
                xargs <- inla.as.list.CRS(x)
                for (name in names(args)) {
                    xargs[[name]] <- args[[name]]
                }
                projargs <- inla.as.CRSargs.list(xargs)
            }
            SRS_string <- sf::st_crs(projargs)$wkt
            projargs <- NULL
        }
    }
    
    if (!is.null(SRS_string)) {
        if (!is.null(projargs)) {
            warning("SRS_string specified. Ignoring non-null projargs.")
        }
        x <-
            CRS(SRS_string = SRS_string)
    } else if (inherits(projargs, "CRS")) {
        x <- projargs
    } else {
        x <- CRS(NA_character_)
    }
    
    if (!is.null(oblique)) {
        stopifnot(is.vector(oblique))
        if (length(oblique) < 4) {
            oblique <- c(oblique, rep(0, 4 - length(oblique)))
        }
        x <- list(crs = x, oblique = oblique)
        class(x) <- "inla.CRS"
    }
    x
}

#' @return `inla.wkt_predef` returns a WKT2 string defining a projection
#' @examples
#' \dontrun{
#' names(inla.wkt_predef())
#' }
#' @export
#' @rdname inla.CRS

inla.wkt_predef <- function() {
    list(
        hammer_norm = 'PROJCRS["unknown",BASEGEOGCRS["unknown",DATUM["unknown",ELLIPSOID["unknown",0.707106781186548,0,LENGTHUNIT["metre",1,ID["EPSG",9001]]]],PRIMEM["Reference meridian",0,ANGLEUNIT["degree",0.0174532925199433,ID["EPSG",9122]]]],CONVERSION["unknown",METHOD["PROJ hammer"]],CS[Cartesian,2],AXIS["(E)",east,ORDER[1],LENGTHUNIT["metre",1,ID["EPSG",9001]]],AXIS["(N)",north,ORDER[2],LENGTHUNIT["metre",1,ID["EPSG",9001]]]]',
        lambert_norm = 'PROJCRS["unknown",BASEGEOGCRS["unknown",DATUM["unknown",ELLIPSOID["unknown",1,0,LENGTHUNIT["metre",1,ID["EPSG",9001]]]],PRIMEM["Reference meridian",0,ANGLEUNIT["degree",0.0174532925199433,ID["EPSG",9122]]]],CONVERSION["unknown",METHOD["Lambert Cylindrical Equal Area (Spherical)",ID["EPSG",9834]],PARAMETER["Latitude of 1st standard parallel",0,ANGLEUNIT["degree",0.0174532925199433],ID["EPSG",8823]],PARAMETER["Longitude of natural origin",0,ANGLEUNIT["degree",0.0174532925199433],ID["EPSG",8802]],PARAMETER["False easting",0,LENGTHUNIT["metre",1],ID["EPSG",8806]],PARAMETER["False northing",0,LENGTHUNIT["metre",1],ID["EPSG",8807]]],CS[Cartesian,2],AXIS["(E)",east,ORDER[1],LENGTHUNIT["metre",1,ID["EPSG",9001]]],AXIS["(N)",north,ORDER[2],LENGTHUNIT["metre",1,ID["EPSG",9001]]]]',
        longlat_norm = 'GEOGCRS["unknown",DATUM["unknown",ELLIPSOID["unknown",1,0,LENGTHUNIT["metre",1,ID["EPSG",9001]]]],PRIMEM["Reference meridian",0,ANGLEUNIT["degree",0.0174532925199433,ID["EPSG",9122]]],CS[ellipsoidal,2],AXIS["longitude",east,ORDER[1],ANGLEUNIT["degree",0.0174532925199433,ID["EPSG",9122]]],AXIS["latitude",north,ORDER[2],ANGLEUNIT["degree",0.0174532925199433,ID["EPSG",9122]]]]',
        mollweide_norm = 'PROJCRS["unknown",BASEGEOGCRS["unknown",DATUM["unknown",ELLIPSOID["unknown",0.707106781186548,0,LENGTHUNIT["metre",1,ID["EPSG",9001]]]],PRIMEM["Reference meridian",0,ANGLEUNIT["degree",0.0174532925199433,ID["EPSG",9122]]]],CONVERSION["unknown",METHOD["Mollweide"],PARAMETER["Longitude of natural origin",0,ANGLEUNIT["degree",0.0174532925199433],ID["EPSG",8802]],PARAMETER["False easting",0,LENGTHUNIT["metre",1],ID["EPSG",8806]],PARAMETER["False northing",0,LENGTHUNIT["metre",1],ID["EPSG",8807]]],CS[Cartesian,2],AXIS["(E)",east,ORDER[1],LENGTHUNIT["metre",1,ID["EPSG",9001]]],AXIS["(N)",north,ORDER[2],LENGTHUNIT["metre",1,ID["EPSG",9001]]]]',
        hammer_globe = 'PROJCRS["unknown",BASEGEOGCRS["unknown",DATUM["Unknown based on Normal Sphere (r=6370997) ellipsoid",ELLIPSOID["Normal Sphere (r=6370997)",6370997,0,LENGTHUNIT["metre",1,ID["EPSG",9001]]]],PRIMEM["Greenwich",0,ANGLEUNIT["degree",0.0174532925199433],ID["EPSG",8901]]],CONVERSION["unknown",METHOD["PROJ hammer"]],CS[Cartesian,2],AXIS["(E)",east,ORDER[1],LENGTHUNIT["kilometre",1000,ID["EPSG",9036]]],AXIS["(N)",north,ORDER[2],LENGTHUNIT["kilometre",1000,ID["EPSG",9036]]]]',
        lambert_globe = 'PROJCRS["unknown",BASEGEOGCRS["unknown",DATUM["Unknown based on Normal Sphere (r=6370997) ellipsoid",ELLIPSOID["Normal Sphere (r=6370997)",6370997,0,LENGTHUNIT["metre",1,ID["EPSG",9001]]]],PRIMEM["Greenwich",0,ANGLEUNIT["degree",0.0174532925199433],ID["EPSG",8901]]],CONVERSION["unknown",METHOD["Lambert Cylindrical Equal Area (Spherical)",ID["EPSG",9834]],PARAMETER["Latitude of 1st standard parallel",0,ANGLEUNIT["degree",0.0174532925199433],ID["EPSG",8823]],PARAMETER["Longitude of natural origin",0,ANGLEUNIT["degree",0.0174532925199433],ID["EPSG",8802]],PARAMETER["False easting",0,LENGTHUNIT["kilometre",1000],ID["EPSG",8806]],PARAMETER["False northing",0,LENGTHUNIT["kilometre",1000],ID["EPSG",8807]]],CS[Cartesian,2],AXIS["(E)",east,ORDER[1],LENGTHUNIT["kilometre",1000,ID["EPSG",9036]]],AXIS["(N)",north,ORDER[2],LENGTHUNIT["kilometre",1000,ID["EPSG",9036]]]]',
        longlat_globe = 'GEOGCRS["unknown",DATUM["Unknown based on Normal Sphere (r=6370997) ellipsoid",ELLIPSOID["Normal Sphere (r=6370997)",6370997,0,LENGTHUNIT["metre",1,ID["EPSG",9001]]]],PRIMEM["Greenwich",0,ANGLEUNIT["degree",0.0174532925199433],ID["EPSG",8901]],CS[ellipsoidal,2],AXIS["longitude",east,ORDER[1],ANGLEUNIT["degree",0.0174532925199433,ID["EPSG",9122]]],AXIS["latitude",north,ORDER[2],ANGLEUNIT["degree",0.0174532925199433,ID["EPSG",9122]]]]',
        mollweide_globe = 'PROJCRS["unknown",BASEGEOGCRS["unknown",DATUM["Unknown based on Normal Sphere (r=6370997) ellipsoid",ELLIPSOID["Normal Sphere (r=6370997)",6370997,0,LENGTHUNIT["metre",1,ID["EPSG",9001]]]],PRIMEM["Greenwich",0,ANGLEUNIT["degree",0.0174532925199433],ID["EPSG",8901]]],CONVERSION["unknown",METHOD["Mollweide"],PARAMETER["Longitude of natural origin",0,ANGLEUNIT["degree",0.0174532925199433],ID["EPSG",8802]],PARAMETER["False easting",0,LENGTHUNIT["kilometre",1000],ID["EPSG",8806]],PARAMETER["False northing",0,LENGTHUNIT["kilometre",1000],ID["EPSG",8807]]],CS[Cartesian,2],AXIS["(E)",east,ORDER[1],LENGTHUNIT["kilometre",1000,ID["EPSG",9036]]],AXIS["(N)",north,ORDER[2],LENGTHUNIT["kilometre",1000,ID["EPSG",9036]]]]',
        sphere = 'GEODCRS["unknown",DATUM["unknown",ELLIPSOID["unknown",1,0,LENGTHUNIT["metre",1,ID["EPSG",9001]]]],PRIMEM["Reference meridian",0,ANGLEUNIT["degree",0.0174532925199433,ID["EPSG",9122]]],CS[Cartesian,3],AXIS["(X)",geocentricX,ORDER[1],LENGTHUNIT["metre",1,ID["EPSG",9001]]],AXIS["(Y)",geocentricY,ORDER[2],LENGTHUNIT["metre",1,ID["EPSG",9001]]],AXIS["(Z)",geocentricZ,ORDER[3],LENGTHUNIT["metre",1,ID["EPSG",9001]]]]',
        globe = 'GEODCRS["unknown",DATUM["Unknown based on Normal Sphere (r=6370997) ellipsoid",ELLIPSOID["Normal Sphere (r=6370997)",6370997,0,LENGTHUNIT["metre",1,ID["EPSG",9001]]]],PRIMEM["Greenwich",0,ANGLEUNIT["degree",0.0174532925199433],ID["EPSG",8901]],CS[Cartesian,3],AXIS["(X)",geocentricX,ORDER[1],LENGTHUNIT["kilometre",1000,ID["EPSG",9036]]],AXIS["(Y)",geocentricY,ORDER[2],LENGTHUNIT["kilometre",1000,ID["EPSG",9036]]],AXIS["(Z)",geocentricZ,ORDER[3],LENGTHUNIT["kilometre",1000,ID["EPSG",9036]]]]'
    )
}





#' Internal WKT handling
#'
#' Conversion between WKT and a tree representation
#'
#'
#' @aliases inla.as.wkt_tree.wkt inla.as.wkt.wkt_tree inla.wkt_tree_get_item
#' inla.wkt_tree_set_item
#' @param x A WKT2 string, or a `wkt_tree` list structure
#' @param \dots Unused
#' @param item character vector with item labels identifying a parameter item
#' entry.
#' @param duplicate For items that have more than one match, `duplicate`
#' indicates the index number of the desired version. Default: 1
#' @param item_tree An item tree identifying a parameter item entry
#' @export inla.as.wkt_tree.wkt
#' @rdname wkt_tree
inla.as.wkt_tree.wkt <- function(x, ...) {
    # Basic parsing of WKT string
    # ITEM[Param1, Param2, ...]
    # Param can be a constant or an ITEM[...]
    parse_item <- function(x) {
        # Parse item label
        n <- regexpr("\\[", x)
        item_label <- substr(x, 1, n - 1)
        x <- substr(x, n + 1, nchar(x))
        params <- list()
        # Parse parameters
        done <- FALSE
        while (!done) {
            # If [ comes before , or ], it's an Item, otherwise a constant
            n <- regexpr("[],[]", x)
            if (n < 1) {
                # Nothing found
                done <- TRUE
                break
            }
            if (substr(x, n, n) == "[") {
                n <- regexpr("[^ ]", x)
                x <- substr(x, n, nchar(x))
                item_info <- parse_item(x)
                params[[length(params) + 1]] <- item_info$item
                x <- item_info$x
            } else {
                const <- substr(x, 1, n - 1)
                params[[length(params) + 1]] <- const
                x <- substr(x, n, nchar(x))
            }
            n <- regexpr("[],]", x)
            done <- (substr(x, n, n) == "]")
            x <- substr(x, n + 1, nchar(x))
        }
        list(item = list(label = item_label, params = params), x = x)
    }

    x <- gsub("\n", "", x)
    item_info <- parse_item(x)
    item <- item_info$item
    item
}

#' @rdname wkt_tree
#' @param pretty logical
#' @export

inla.as.wkt.wkt_tree <- function(x, pretty = FALSE, ...) {
    construct_item <- function(x, level) {
        paste0(
            if (pretty) {
                paste0(rep("    ", level), collapse = "")
            } else {
                ""
            },
            x[["label"]],
            "[",
            paste0(vapply(
                x[["params"]],
                function(param) {
                    if (!is.list(param)) {
                        paste0(param)
                    } else {
                        paste0(
                            if (pretty) {
                                "\n"
                            } else {
                                ""
                            },
                            construct_item(param,
                                level = level + 1
                            )
                        )
                    }
                },
                ""
            ),
            collapse = ","
            ),
            "]"
        )
    }
    construct_item(x, 0)
}

#' @param item character vector with item labels identifying a parameter item
#' entry.
#' @param duplicate For items that have more than one match, `duplicate`
#' indicates the index number of the desired version. Default: 1
#' @rdname wkt_tree
#' @export

inla.wkt_tree_get_item <- function(x, item, duplicate = 1) {
    for (k in seq_along(x[["params"]])) {
        if (is.list(x[["params"]][[k]]) &&
            (!is.null(x[["params"]][[k]][["label"]])) &&
            (x[["params"]][[k]][["label"]] %in% item)) {
            if (duplicate == 1) {
                return(x[["params"]][[k]])
            }
            duplicate <- duplicate - 1
        }
    }
    NULL
}

#' @param item_tree An item tree identifying a parameter item entry
#' @rdname wkt_tree
#' @export

inla.wkt_tree_set_item <- function(x, item_tree, duplicate = 1) {
    success <- FALSE
    for (k in seq_along(x[["params"]])) {
        if (is.list(x[["params"]][[k]]) && (x[["params"]][[k]][["label"]] == item_tree[["label"]])) {
            if (duplicate == 1) {
                x[["params"]][[k]] <- item_tree
            }
            duplicate <- duplicate - 1
            success <- TRUE
            break
        }
    }
    if (!success) {
        x[["params"]] <- c(x[["params"]], list(item_tree))
    }
    x
}



#' Show expanded CRS arguments
#'
#' Wrapper for `sp::CRS` and `inla.CRS` objects to extract the
#' coordinate reference system argument string.
#' Should no longer be used with PROJ6;
#' see [inla.crs_get_wkt()]
#'
#' @aliases inla.CRSargs inla.as.list.CRS inla.as.list.CRSargs inla.as.CRS.list
#' inla.as.CRSargs.list
#' @param x An `sp::CRS` or `inla.CRS` object (for
#' `inla.CRSargs` and `inla.as.list.CRS`), a character string (for
#' `inla.as.list.CRSargs`), or a list (for `inla.as.CRS.list` and
#' `inla.as.CRSargs.list`).
#' @param \dots Additional arguments passed on to other methods.
#' @return For `inla.CRSargs` and `inla.as.CRSargs.list`, a character
#' string with PROJ.4 arguments.
#'
#' For `inla.as.list.CRS` and `inla.as.list.CRSargs`, a list of
#' name/value pairs.
#'
#' For `inla.as.CRS.list`, a `CRS` or `inla.CRS` object.
#' @author Finn Lindgren <finn.lindgren@@gmail.com>
#' @seealso [inla.CRS()]
#' @examples
#'
#' if (require("sf") && require("sp")) {
#'     crs0 <- inla.CRS("longlat")
#'     p4s <- inla.CRSargs(crs0)
#'     lst <- inla.as.list.CRSargs(p4s)
#'     crs1 <- inla.as.CRS.list(lst)
#'     lst$a <- 2
#'     crs2 <- inla.CRS(p4s, args = lst)
#'     print(inla.CRSargs(crs0))
#'     print(inla.CRSargs(crs1))
#'     print(inla.CRSargs(crs2))
#' }
#' @export
#' @rdname CRSargs
inla.CRSargs <- function(x, ...) {
    inla.not_for_PROJ6("inla.CRSargs")

    if (inherits(x, "inla.CRS")) {
        x <- x[["crs"]]
    }
    if (is.null(x)) {
        as.character(NA)
    } else {
        inla.require("sf", stop.on.error = TRUE)
        sf::st_crs(x)$proj4string
    }
}


#' @details * `inla.as.CRSargs.list`: CRS proj4 string for name=value pair list
#'
#' @export
#' @rdname CRSargs
inla.as.CRSargs.list <- function(x, ...) {
    paste(lapply(
        names(x),
        function(xx) {
            if (is.na(x[[xx]])) {
                paste("+", xx, sep = "")
            } else {
                paste("+", xx, "=", x[[xx]], sep = "")
            }
        }
    ),
    collapse = " "
    )
}

#' @details * `inla.as.list.CRSargs`: List of name=value pairs from CRS proj4 string
#'
#' @export
#' @rdname CRSargs
inla.as.list.CRSargs <- function(x, ...) {
    inla.not_for_PROJ6("inla.as.list.CRSargs")

    if (is.na(x)) {
        return(list())
    }
    if (!is.character(x)) {
        stop("proj4string must be of class character")
    }
    do.call(c, lapply(
        strsplit(
            x = strsplit(
                x = paste(" ", x, sep = ""),
                split = " \\+"
            )[[1]][-1],
            split = "="
        ),
        function(x) {
            xx <- list(x[2])
            names(xx) <- x[1]
            xx
        }
    ))
}


#' @export
#' @rdname CRSargs
inla.as.list.CRS <- function(x, ...) {
    inla.not_for_PROJ6("inla.as.list.CRS")
    inla.as.list.CRSargs(inla.CRSargs(x))
}

#' @export
#' @rdname CRSargs
inla.as.CRS.list <- function(x, ...) {
    inla.not_for_PROJ6("inla.as.CRS.list")
    inla.CRS(args = x)
}





#' @param crs A `sp::CRS` or `inla.CRS` object
#' @param wkt A WKT2 character string
#' @param unit character, name of a unit. Supported names are
#' "metre", "kilometre", and the aliases "meter", "m", International metre",
#' "kilometer", and "km", as defined by `inla.wkt_unit_params` or the
#' `params` argument. (For legacy PROJ4 use, only "m" and "km" are
#' supported)
#' @param params Length unit definitions, in the list format produced by
#' `inla.wkt_unit_params()`, Default: NULL, which invokes
#' `inla.wkt_unit_params()`
#' @return For `inla.wkt_unit_params`, a
#' list of named unit definitions
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @examples
#' \dontrun{
#'   c1 <- inla.CRS("globe")
#'   inla.crs_get_lengthunit(c1)
#'   c2 <- inla.crs_set_lengthunit(c1, "km")
#'   inla.crs_get_lengthunit(c2)
#' }
#' @export
#' @seealso [inla.sp_get_crs()]
#' @aliases crs_wkt
#' @rdname crs_wkt

inla.wkt_unit_params <- function() {
    params <- list(
        "metre" =
            list(
                '"metre"',
                "1",
                list(
                    label = "ID",
                    params = list('"EPSG"', "9001")
                )
            ),
        "kilometre" =
            list(
                '"kilometre"',
                "1000",
                list(
                    label = "ID",
                    params = list('"EPSG"', "9036")
                )
            )
    )
    params[["meter"]] <- params[["metre"]]
    params[["m"]] <- params[["metre"]]
    params[["International metre"]] <- params[["metre"]]
    params[["kilometer"]] <- params[["kilometre"]]
    params[["km"]] <- params[["kilometre"]]
    params
}

#' @export
#' @rdname crs_wkt
#' @return For `inla.wkt_get_lengthunit`, a
#' list of length units used in the wkt string, excluding the ellipsoid radius
#' unit.

inla.wkt_get_lengthunit <- function(wkt) {
    extract <- function(wt) {
        # 1. Recursively find LENGTHUNIT, except within ELLIPSOID
        # 2. Return unit

        if (wt[["label"]] == "LENGTHUNIT") {
            result <- list(wt[["params"]])
        } else if (wt[["label"]] != "ELLIPSOID") {
            result <- list()
            for (k in seq_along(wt$param)) {
                if (is.list(wt[["params"]][[k]])) {
                    result <- c(result, extract(wt[["params"]][[k]]))
                }
            }
        } else {
            result <- list()
        }
        result
    }

    wt <- inla.as.wkt_tree.wkt(wkt)
    params <- unique(extract(wt))
    names(params) <-
        vapply(
            params,
            function(x) {
                gsub('"', "", x[[1]])
            },
            ""
        )
    params
}

#' @export
#' @rdname crs_wkt
#' @return For `inla.wkt_set_lengthunit`, a
#' WKT2 string with altered length units.
#' Note that the length unit for the ellipsoid radius is unchanged.

inla.wkt_set_lengthunit <- function(wkt, unit, params = NULL) {
    convert <- function(wt, unit) {
        # 1. Recursively find LENGTHUNIT, except within ELLIPSOID
        # 2. Change unit

        if (wt[["label"]] == "LENGTHUNIT") {
            wt[["params"]] <- unit
        } else if (wt[["label"]] != "ELLIPSOID") {
            for (k in seq_along(wt$param)) {
                if (is.list(wt[["params"]][[k]])) {
                    wt[["params"]][[k]] <- convert(wt[["params"]][[k]], unit)
                }
            }
        }
        wt
    }

    if (is.null(params)) {
        params <- inla.wkt_unit_params()
    }
    if (!(unit %in% names(params))) {
        warning(paste0(
            "'inla.wkt_set_lengthunit' unit conversion to '",
            unit,
            "' not supported. Unit left unchanged."
        ))
        return(wkt)
    }

    wt <- inla.as.wkt_tree.wkt(wkt)
    wt <- convert(wt, params[[unit]])
    inla.as.wkt.wkt_tree(wt)
}

#' @return For `inla.crs_get_wkt`, WKT2 string.
#' @export
#' @rdname crs_wkt

inla.crs_get_wkt <- function(crs) {
    if (inherits(crs, "inla.CRS")) {
        crs <- crs[["crs"]]
    }

    if (is.null(crs)) {
        return(NULL)
    }

    comment(crs)
}

#' @return For `inla.crs_get_lengthunit`, a
#' list of length units used in the wkt string, excluding the ellipsoid radius
#' unit. (For legacy PROJ4 code, the raw units from the proj4string are
#' returned, if present.)
#' @export
#' @rdname crs_wkt

inla.crs_get_lengthunit <- function(crs) {
    inla.wkt_get_lengthunit(inla.crs_get_wkt(crs))
}

#' @return For `inla.crs_set_lengthunit`, a `sp::CRS` object with
#' altered length units.
#' Note that the length unit for the ellipsoid radius is unchanged.
#' @export
#' @rdname crs_wkt

inla.crs_set_lengthunit <- function(crs, unit, params = NULL) {
    if (inherits(crs, "inla.CRS")) {
        crs_ <- crs
        crs <- crs[["crs"]]
    } else {
        crs_ <- NULL
    }
    x <- sp::CRS(SRS_string = inla.wkt_set_lengthunit(inla.crs_get_wkt(crs),
                                                      unit,
                                                      params = params
    ))
    if (!is.null(crs_)) {
        crs_[["crs"]] <- x
        x <- crs_
    }
    x
}


inla.rotmat3213 <- function(rot) {
    cs <- cos(rot[1])
    sn <- sin(rot[1])
    R <- matrix(c(
        cs, -sn, 0,
        sn, cs, 0,
        0, 0, 1
    ), 3, 3)
    cs <- cos(rot[2])
    sn <- sin(rot[2])
    R <- R %*% matrix(c(
        cs, 0, sn,
        0, 1, 0,
        -sn, 0, cs
    ), 3, 3)
    cs <- cos(rot[3])
    sn <- sin(rot[3])
    R <- R %*% matrix(c(
        1, 0, 0,
        0, cs, -sn,
        0, sn, cs
    ), 3, 3)
    cs <- cos(rot[4])
    sn <- sin(rot[4])
    R <- R %*% matrix(c(
        cs, -sn, 0,
        sn, cs, 0,
        0, 0, 1
    ), 3, 3)
    R
}

inla.rotmat3123 <- function(rot) {
    cs <- cos(rot[4])
    sn <- sin(rot[4])
    R <- matrix(c(
        cs, -sn, 0,
        sn, cs, 0,
        0, 0, 1
    ), 3, 3)
    cs <- cos(rot[3])
    sn <- sin(rot[3])
    R <- R %*% matrix(c(
        1, 0, 0,
        0, cs, -sn,
        0, sn, cs
    ), 3, 3)
    cs <- cos(rot[2])
    sn <- sin(rot[2])
    R <- R %*% matrix(c(
        cs, 0, sn,
        0, 1, 0,
        -sn, 0, cs
    ), 3, 3)
    cs <- cos(rot[1])
    sn <- sin(rot[1])
    R <- R %*% matrix(c(
        cs, -sn, 0,
        sn, cs, 0,
        0, 0, 1
    ), 3, 3)
    R
}


inla.crs.transform.oblique <- function(x, oblique, to.oblique = TRUE) {
    if (to.oblique) {
        ## Transform to oblique orientation
        ## 1) Rotate -oblique[1] around (0,0,1)
        ## 2) Rotate +oblique[2] around (0,1,0)
        ## 3) Rotate -oblique[3] around (1,0,0)
        ## 3) Rotate -oblique[4] around (0,0,1)
        x %*% inla.rotmat3213(c(-1, 1, -1, -1) * oblique * pi / 180)
    } else {
        ## Transform back from oblique orientation
        ## 1) Rotate +oblique[4] around (0,0,1)
        ## 2) Rotate +oblique[3] around (1,0,0)
        ## 3) Rotate -oblique[2] around (0,1,0)
        ## 4) Rotate +oblique[1] around (0,0,1)
        x %*% inla.rotmat3123(c(1, -1, 1, 1) * oblique * pi / 180)
    }
}



inla.wkt_tree_projection_type <- function(wt) {
    axis1 <- inla.wkt_tree_get_item(wt, "AXIS", 1)
    axis2 <- inla.wkt_tree_get_item(wt, "AXIS", 2)
    if (identical(axis1[["params"]][[1]], '"longitude"') &&
        identical(axis2[["params"]][[1]], '"latitude"')) {
        return("longlat")
    }
    conversion <- inla.wkt_tree_get_item(wt, "CONVERSION")
    if (!is.null(conversion)) {
        method <- inla.wkt_tree_get_item(conversion, "METHOD")
        if (identical(method[["params"]][[1]], '"Lambert Cylindrical Equal Area (Sherical)"')) {
            return("lambert")
        }
        if (identical(method[["params"]][[1]], '"Mollweide"')) {
            return("mollweide")
        }
        if (identical(method[["params"]][[1]], '"PROJ hammer"')) {
            return("hammer")
        }
        if (identical(method[["params"]][[1]], '"tmerc"')) {
            return("tmerc")
        }
    }
    NULL
}

inla.wkt_projection_type <- function(wkt) {
    wt <- inla.as.wkt_tree.wkt(wkt)
    inla.wkt_tree_projection_type(wt)
}

inla.crs_projection_type <- function(crs) {
    wkt <- inla.crs_get_wkt(crs)
    inla.wkt_projection_type(wkt)
}

## +proj=longlat in (-180,180)x(-90,90)
## +proj=moll in (-2,2)x(-1,1) scaled by +a and +b, and +units
## +proj=lambert in (-pi,pi)x(-1,1) scaled by +a and +b, and +units
inla.crs.bounds <- function(crs, warn.unknown = FALSE) {
    wkt <- inla.crs_get_wkt(crs)
    wt <- inla.as.wkt_tree.wkt(wkt)
    type <- inla.wkt_tree_projection_type(wt)
    
    if (is.null(type)) {
        if (inla.wkt_is_geocent(wkt)) {
            bounds <-
                list(
                    type = "rectangle",
                    xlim = c(-Inf, Inf),
                    ylim = c(-Inf, Inf)
                )
        } else {
            if (warn.unknown) {
                warning(
                    "Could not determine shape of transformation bounds. Using infinite rectangle."
                )
            }
            bounds <-
                list(
                    type = "rectangle",
                    xlim = c(-Inf, Inf),
                    ylim = c(-Inf, Inf)
                )
        }
    } else if (type == "longlat") {
        bounds <-
            list(
                type = "rectangle",
                xlim = c(-180, 180),
                ylim = c(-90, 90)
            )
    } else if (type == "lambert") {
        axis <- c(pi, 1)
        radius <- inla.wkt_get_ellipsoid_radius(wkt)
        axis[1] <- axis[1] * radius
        # TODO: handle eccentricity
        axis[2] <- axis[2] * sqrt(radius) * sqrt(radius)
        # TODO: Handle units"
        bounds <- list(
            type = "rectangle",
            xlim = c(-1, 1) * axis[1],
            ylim = c(-1, 1) * axis[2]
        )
    } else if (type %in% c("mollweide", "hammer")) {
        axis <- c(2, 1)
        center <- c(0, 0)
        radius <- inla.wkt_get_ellipsoid_radius(wkt)
        axis[1] <- axis[1] * radius / sqrt(1 / 2)
        axis[2] <- axis[2] * radius / sqrt(1 / 2)
        # TODO: Handle "units"
        bounds <- list(
            type = "ellipse",
            axis = axis,
            center = center,
            xlim = center[1] + c(-1, 1) * axis[1],
            ylim = center[2] + c(-1, 1) * axis[2]
        )
    } else if (type == "tmerc") {
        bounds <-
            list(
                type = "rectangle",
                xlim = c(-Inf, Inf),
                ylim = c(-Inf, Inf)
            )
    } else {
        stop("'inla.crs.bounds' internal error: transformation detected but not handled.")
    }
    
    if (bounds$type == "rectangle") {
        bounds$polygon <- cbind(
            bounds$xlim[c(1, 2, 2, 1, 1)],
            bounds$ylim[c(1, 1, 2, 2, 1)]
        )
    } else if (bounds$type == "ellipse") {
        theta <- seq(0, 2 * pi, length = 1000)
        bounds$polygon <- cbind(
            bounds$center[1] + bounds$axis[1] * cos(theta),
            bounds$center[2] + bounds$axis[2] * sin(theta)
        )
    } else {
        stop("Unknown transformation type. This should not happen.")
    }
    bounds
}

## TRUE/FALSE for points inside/outside projection domain.
inla.crs.bounds.check <- function(x, bounds) {
    inla.require.inherits(x, "matrix")
    if (all(is.finite(bounds$xlim)) && all(is.finite(bounds$ylim))) {
        (sp::point.in.polygon(
            x[, 1], x[, 2],
            bounds$polygon[, 1], bounds$polygon[, 2]
        )
        > 0)
    } else {
        ok <- rep(TRUE, nrow(x))
    }
}



internal.update.crs <- function(crs, newcrs, mismatch.allowed) {
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




#' Test CRS and inla.CRS for equality
#'
#' Wrapper for identical, optionally testing only the CRS part of two objects
#'
#'
#' @param crs0 A `CRS` or `inla.CRS` object.
#' @param crs1 A `CRS` or `inla.CRS` object.
#' @param crsonly Logical. If `TRUE`, only the `CRS` part of a
#' `inla.CRS` object is compared.
#' @author Finn Lindgren <finn.lindgren@@gmail.com>
#' @seealso [inla.CRS()]
#' @examples
#'
#' crs0 <- inla.CRS("longlat")
#' crs1 <- inla.CRS("longlat", oblique = c(0, 90))
#' print(c(
#'     inla.identical.CRS(crs0, crs0),
#'     inla.identical.CRS(crs0, crs1),
#'     inla.identical.CRS(crs0, crs1, crsonly = TRUE)
#' ))
#' @export inla.identical.CRS
inla.identical.CRS <- function(crs0, crs1, crsonly = FALSE) {
    if (!crsonly) {
        identical(crs0, crs1)
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






#' Wrapper method for `sp::spTransform`
#'
#' Handles transformation of various inla objects accorting to coordinate
#' reference systems of `sp::CRS` or `inla.CRS` class.
#'
#'
#' @aliases inla.spTransform inla.spTransform.default
#' inla.spTransform.SpatialPoints inla.spTransform.inla.mesh.lattice
#' inla.spTransform.inla.mesh.segment inla.spTransform.inla.mesh
#' @param x The object that should be transformed from it's current CRS to a
#' new CRS
#' @param crs0 The source `sp::CRS` or `inla.CRS` object
#' @param crs1 The target `sp::CRS` or `inla.CRS` object
#' @param CRSobj The target `sp::CRS` or `inla.CRS` object
#' @param passthrough default FALSE. Setting to TRUE allows objects with no CRS
#' information to be passed through without transformation.
#' @param \dots Potential additional arguments
#' @return The object is returned with its coordinates transformed
#' @author Finn Lindgren <finn.lindgren@@gmail.com>
#' @seealso [inla.CRS()]
#' @examples
#'
#' if (require("sf") && require("sp")) {
#'     latt <- inla.mesh.lattice(-10:10, 40:60)
#'     mesh1 <- inla.mesh.create(
#'         lattice = latt, extend = FALSE, refine = FALSE,
#'         crs = inla.CRS("longlat")
#'     )
#'     mesh2 <- inla.spTransform(mesh1, inla.CRS("lambert"))
#'     summary(mesh1)
#'     summary(mesh2)
#' }
#' @export inla.spTransform
inla.spTransform <- function(x, ...) {
    UseMethod("inla.spTransform")
}



#' @details `inla.spTransform.default` Low level transformation of raw coordinates.
#' @export
#' @rdname inla.spTransform
inla.spTransform.default <- function(x, crs0, crs1, passthrough = FALSE, ...) {
    # PROJ6
    ok0 <- (!is.null(crs0) &&
                ((
                    inherits(crs0, "CRS") && !is.null(inla.crs_get_wkt(crs0))
                ) ||
                    (inherits(crs0, "inla.CRS"))))
    ok1 <- (!is.null(crs1) &&
                ((
                    inherits(crs1, "CRS") && !is.null(inla.crs_get_wkt(crs1))
                ) ||
                    (inherits(crs1, "inla.CRS"))))
    if (ok0 && ok1) {
        if (ncol(x) == 2) {
            x <- cbind(x, 0)
        }
        sphere_radius_0 <- inla.crs_get_ellipsoid_radius(crs0)
        sphere_radius_1 <- inla.crs_get_ellipsoid_radius(crs1)
        different_radii <- (sphere_radius_0 != sphere_radius_1)
        longlat_norm <- inla.CRS("longlat_norm")
        longlat_0 <-
            inla.crs_set_ellipsoid_radius(longlat_norm, sphere_radius_0)
        longlat_1 <-
            inla.crs_set_ellipsoid_radius(longlat_norm, sphere_radius_1)
        
        crs_sphere <- inla.CRS("sphere")
        onsphere_0 <-
            inla.identical.CRS(crs0, crs_sphere, crsonly = TRUE)
        onsphere_1 <-
            inla.identical.CRS(crs1, crs_sphere, crsonly = TRUE)
        is_geocentric_0 <- inla.crs_is_geocent(crs0)
        is_geocentric_1 <- inla.crs_is_geocent(crs1)
        if (is_geocentric_0) {
            ok <- TRUE
        } else {
            bounds <- inla.crs.bounds(crs0)
            if (identical(inla.crs_projection_type(crs0), "longlat")) {
                ## Wrap longitudes to [-180,180]
                needswrap <- (x[, 1] < -180) | (x[, 1] > 180)
                if (any(needswrap)) {
                    x[needswrap, 1] <- ((x[needswrap, 1] + 180) %% 360) - 180
                }
            }
            ok <- inla.crs.bounds.check(x, bounds)
            if (!all(ok)) {
                xx <- x
            }
        }
        do_work_on_sphere <-
            inherits(crs0, "inla.CRS") ||
            inherits(crs1, "inla.CRS") ||
            different_radii
        if (inherits(crs0, "inla.CRS")) {
            crs0crs <- crs0$crs
            crs0oblique <- crs0$oblique
        } else {
            crs0crs <- crs0
            crs0oblique <- NULL
        }
        if (inherits(crs1, "inla.CRS")) {
            crs1crs <- crs1$crs
            crs1oblique <- crs1$oblique
        } else {
            crs1crs <- crs1
            crs1oblique <- NULL
        }
        x <-
            SpatialPoints(x[ok, , drop = FALSE], proj4string = crs0crs)
        if (do_work_on_sphere) {
            if (!onsphere_0) {
                if (sphere_radius_0 != 1) {
                    x <- spTransform(x, longlat_0)
                    proj4string(x) <-
                        CRS(NA_character_) # Reset CRS to avoid warning
                    proj4string(x) <- longlat_norm
                }
                x <- spTransform(x, crs_sphere)
            }
            if (!is.null(crs0oblique)) {
                x <- SpatialPoints(
                    inla.crs.transform.oblique(coordinates(x),
                                               crs0oblique,
                                               to.oblique = FALSE),
                    proj4string = crs_sphere
                )
            }
            
            if (!is.null(crs1oblique)) {
                x <- SpatialPoints(
                    inla.crs.transform.oblique(coordinates(x),
                                               crs1oblique,
                                               to.oblique = TRUE),
                    proj4string = crs_sphere
                )
            }
            if (sphere_radius_1 != 1) {
                x <- spTransform(x, longlat_norm)
                proj4string(x) <-
                    CRS(NA_character_) # Reset CRS to avoid warning
                proj4string(x) <- longlat_1
            }
        }
        
        x <- spTransform(x, crs1crs)
        
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


#' @export
#' @rdname inla.spTransform
inla.spTransform.SpatialPoints <- function(x, CRSobj, passthrough = FALSE, ...) {
    crs_x <- inla.sp_get_crs(x)
    ok0 <- !is.null(inla.crs_get_wkt(crs_x))
    ok1 <- (!missing(CRSobj) && !is.null(CRSobj) &&
                (inherits(CRSobj, "CRS") &&
                     !is.null(inla.crs_get_wkt(CRSobj))))
    if (ok0 && ok1) {
        invisible(SpatialPoints(
            inla.spTransform(coordinates(x),
                             crs_x,
                             CRSobj),
            proj4string = CRSobj
        ))
    } else if (ok1) {
        ## Know: !ok0 && ok1
        if (!passthrough) {
            stop("Invalid origin CRS for SpatialPoints")
        }
        invisible(SpatialPoints(coordinates(x), proj4string = CRSobj))
    } else {
        ## Know: (ok0 || !ok0) && !ok1
        if (!passthrough) {
            stop("Invalid target CRS for SpatialPoints")
        }
        invisible(SpatialPoints(coordinates(x), proj4string = inla.CRS()))
    }
}

#' @export
#' @rdname inla.spTransform
inla.spTransform.SpatialPointsDataFrame <- function(x,
                                                    CRSobj,
                                                    passthrough = FALSE,
                                                    ...) {
    ok1 <- (!missing(CRSobj) && !is.null(CRSobj) &&
        (inherits(CRSobj, "CRS") && !is.null(inla.crs_get_wkt(CRSobj))))
    if (!ok1 && !passthrough) {
        stop("Invalid target CRS for SpatialPointsDataFrame")
    }

    x_no_df <- SpatialPoints(coordinates(x),
        proj4string = inla.sp_get_crs(x)
    )
    x_no_df <- inla.spTransform(x_no_df,
        CRSobj = CRSobj,
        passthrough = passthrough
    )
    if (ok1) {
        invisible(SpatialPointsDataFrame(coordinates(x_no_df),
            proj4string = CRSobj,
            data = x@data
        ))
    } else {
        invisible(SpatialPointsDataFrame(coordinates(x_no_df),
            proj4string = inla.CRS(),
            data = x@data
        ))
    }
}

#' @export
#' @rdname inla.spTransform
inla.spTransform.inla.mesh.lattice <- function(x, CRSobj, passthrough = FALSE, ...) {
    x$segm <- inla.spTransform(x$segm, CRSobj, passthrough = passthrough)
    x$loc <- inla.spTransform(x$loc, x$crs, CRSobj, passthrough = passthrough)
    x$crs <- CRSobj
    invisible(x)
}

#' @export
#' @rdname inla.spTransform
inla.spTransform.inla.mesh.segment <- function(x, CRSobj, passthrough = FALSE, ...) {
    x$loc <- inla.spTransform(x$loc, x$crs, CRSobj, passthrough = passthrough)
    x$crs <- CRSobj
    invisible(x)
}




inla.crs_detect_manifold <- function(crs) {
    if (inla.crs_is_geocent(crs)) {
        manifold <- "S2"
    } else {
        manifold <- "R2"
    }
    manifold
}

#' @export
#' @rdname inla.spTransform
inla.spTransform.inla.mesh <- function(x, CRSobj, passthrough = FALSE, ...) {
    x$loc <- inla.spTransform(x$loc, x$crs, CRSobj, passthrough = passthrough)
    x$manifold <- inla.crs_detect_manifold(CRSobj)
    x$crs <- CRSobj
    invisible(x)
}




## Input: list of segments, all closed polygons.
inla.internal.sp2segment.join <- function(inp, grp = NULL, closed = TRUE) {
    crs <- NULL
    if (length(inp) > 0) {
        out.loc <- matrix(0, 0, ncol(inp[[1]]$loc))
        for (k in seq_along(inp)) {
            crs <- internal.update.crs(crs, inp[[k]]$crs, mismatch.allowed = FALSE)
        }
    } else {
        out.loc <- matrix(0, 0, 2)
    }
    out.idx <- matrix(0L, 0, 2)
    if (is.null(grp)) {
        out.grp <- NULL
    } else {
        out.grp <- integer(0)
    }
    for (k in seq_along(inp)) {
        inp.loc <- inp[[k]]$loc
        inp.idx <- inp[[k]]$idx
        inp.grp <- inp[[k]]$grp
        offset <- nrow(out.loc)
        n <- nrow(as.matrix(inp.idx))
        if (closed) {
            if (!is.null(grp) && is.null(inp.grp)) {
                inp.grp <- rep(grp[k], n)
            }
            if (ncol(as.matrix(inp.idx)) == 1) {
                inp.idx <- cbind(inp.idx, inp.idx[c(2:n, 1)])
            }
        } else {
            if (!is.null(grp) && is.null(inp.grp)) {
                inp.grp <- rep(grp[k], n - 1)
            }
            if (ncol(as.matrix(inp.idx)) == 1) {
                inp.idx <- cbind(inp.idx[-n], inp.idx[-1])
            }
        }
        out.loc <- rbind(out.loc, inp.loc)
        out.idx <- rbind(out.idx, inp.idx + offset)
        if (!is.null(grp)) {
            out.grp <- c(out.grp, inp.grp)
        }
    }
    inla.mesh.segment(
        loc = out.loc, idx = out.idx, grp = out.grp, is.bnd = FALSE,
        crs = crs
    )
}




#' Convert `sp` curve objects to `inla.mesh.segment` objects.
#'
#' Convert `sp` curve objects to `inla.mesh.segment` objects.
#'
#'
#' @aliases as.inla.mesh.segment as.inla.mesh.segment.Line
#' as.inla.mesh.segment.Lines as.inla.mesh.segment.SpatialLines
#' as.inla.mesh.segment.SpatialLinesDataFrame as.inla.mesh.segment.Polygon
#' as.inla.mesh.segment.Polygons as.inla.mesh.segment.SpatialPolygons
#' as.inla.mesh.segment.SpatialPolygonsDataFrame inla.sp2segment
#' @param sp An `sp` polygon object of class `Polygon`,
#' `Polygons`, `SpatialPolygons`, or `SpatialPolygonsDataFrame`.
#' @param join If `TRUE`, join multiple polygons into a single segment
#' (possibly non-simply connected).
#' @param grp Group ID specification for each polygon, as used by
#' [inla.mesh.segment()], one ID per polygon.
#' @param reverse Logical, indicating if the line sequence should be traversed
#' backwards.
#' @param crs An optional `CRS` or `inla.CRS` object
#' @param \dots Additional arguments passed on to other methods.
#' @param is.bnd logical
#' @return A [inla.mesh.segment()] object, or a list of
#' [inla.mesh.segment()] objects.
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @seealso [inla.mesh.segment()]
#' @export as.inla.mesh.segment
as.inla.mesh.segment <-
    function(sp, ...) {
        UseMethod("as.inla.mesh.segment")
    }

#' @export
#' @rdname as.inla.mesh.segment
inla.sp2segment <-
    function(sp, ...) {
        UseMethod("as.inla.mesh.segment")
    }



#' @export
#' @rdname as.inla.mesh.segment
as.inla.mesh.segment.SpatialPoints <-
    function(sp, reverse = FALSE, grp = NULL, is.bnd = TRUE, ...) {
        crs <- inla.sp_get_crs(sp)
        loc <- coordinates(sp)

        n <- dim(loc)[1L]
        if (reverse) {
            idx <- seq(n, 1L, length = n)
        } else {
            idx <- seq_len(n)
        }
        inla.mesh.segment(
            loc = loc, idx = idx, grp = grp, is.bnd = is.bnd,
            crs = crs
        )
    }

#' @export
#' @rdname as.inla.mesh.segment
as.inla.mesh.segment.SpatialPointsDataFrame <-
    function(sp, ...) {
        as.inla.mesh.segment.SpatialLines(sp, ...)
    }



#' @export
#' @rdname as.inla.mesh.segment
as.inla.mesh.segment.Line <-
    function(sp, reverse = FALSE, crs = NULL, ...) {
        loc <- sp@coords
        n <- dim(loc)[1L]
        if (reverse) {
            idx <- seq(n, 1L, length = n)
        } else {
            idx <- seq_len(n)
        }
        inla.mesh.segment(loc = loc, idx = idx, is.bnd = FALSE, crs = crs)
    }

#' @export
#' @rdname as.inla.mesh.segment
as.inla.mesh.segment.Lines <-
    function(sp, join = TRUE, crs = NULL, ...) {
        segm <- as.list(lapply(
            sp@Lines,
            function(x) as.inla.mesh.segment(x, crs = crs, ...)
        ))
        if (join) {
              segm <- inla.internal.sp2segment.join(segm, grp = NULL, closed = FALSE)
          }
        segm
    }

#' @export
#' @rdname as.inla.mesh.segment
as.inla.mesh.segment.SpatialLines <-
    function(sp, join = TRUE, grp = NULL, ...) {
        crs <- inla.sp_get_crs(sp)
        segm <- list()
        for (k in 1:length(sp@lines)) {
              segm[[k]] <- as.inla.mesh.segment(sp@lines[[k]],
                  join = TRUE,
                  crs = crs, ...
              )
          }
        if (join) {
            if (missing(grp)) {
                grp <- 1:length(segm)
            }
            segm <- inla.internal.sp2segment.join(segm, grp = grp, closed = FALSE)
        }
        segm
    }

#' @export
#' @rdname as.inla.mesh.segment
as.inla.mesh.segment.SpatialLinesDataFrame <-
    function(sp, ...) {
        as.inla.mesh.segment.SpatialLines(sp, ...)
    }

#' @export
#' @rdname as.inla.mesh.segment
as.inla.mesh.segment.SpatialPolygons <-
    function(sp, join = TRUE, grp = NULL, ...) {
        crs <- inla.sp_get_crs(sp)
        segm <- list()
        for (k in 1:length(sp@polygons)) {
              segm[[k]] <- as.inla.mesh.segment(sp@polygons[[k]], join = TRUE, crs = crs)
          }
        if (join) {
            if (missing(grp)) {
                grp <- 1:length(segm)
            }
            segm <- inla.internal.sp2segment.join(segm, grp = grp)
        }
        segm
    }

#' @export
#' @rdname as.inla.mesh.segment
as.inla.mesh.segment.SpatialPolygonsDataFrame <-
    function(sp, ...) {
        as.inla.mesh.segment.SpatialPolygons(sp, ...)
    }

#' @export
#' @rdname as.inla.mesh.segment
as.inla.mesh.segment.Polygons <-
    function(sp, join = TRUE, crs = NULL, ...) {
        segm <- as.list(lapply(
            sp@Polygons,
            function(x) as.inla.mesh.segment(x, crs = crs)
        ))
        if (join) {
              segm <- inla.internal.sp2segment.join(segm, grp = NULL)
          }
        segm
    }

#' @export
#' @rdname as.inla.mesh.segment
as.inla.mesh.segment.Polygon <-
    function(sp, crs = NULL, ...) {
        loc <- sp@coords[-dim(sp@coords)[1L], , drop = FALSE]
        n <- dim(loc)[1L]
        if (sp@hole) {
              if (sp@ringDir == 1) {
                    idx <- c(1L:n, 1L)
                } else {
                    idx <- c(1L, seq(n, 1L, length.out = n))
                }
          } else
        if (sp@ringDir == 1) {
              idx <- c(1L, seq(n, 1L, length.out = n))
          } else {
              idx <- c(1L:n, 1L)
          }
        inla.mesh.segment(loc = loc, idx = idx, is.bnd = TRUE, crs = crs)
    }
