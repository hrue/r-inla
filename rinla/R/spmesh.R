#' PROJ6 detection
#'
#' Detect whether PROJ6 is available for INLA. Deprecated and always returns `TRUE`.
#'
#' @aliases inla.has_PROJ6 inla.not_for_PROJ6 inla.not_for_PROJ4
#' inla.fallback_PROJ6 inla.requires_PROJ6
#' @details `inla.has_PROJ6` is called to check if PROJ6&GDAL3 are available.
#' @return For `inla.has_PROJ6`, always returns `TRUE`. Previously: logical; `TRUE` if PROJ6 is available,
#' `FALSE` otherwise
#' @examples
#' \dontrun{
#' inla.has_PROJ6()
#' }
#' @export inla.has_PROJ6
#' @describeIn inla.has_PROJ6 `r lifecycle::badge("deprecated")`
inla.has_PROJ6 <- function() {
    lifecycle::deprecate_warn(
        "23.01.01",
        "inla.has_PROJ6()",
        details = "inla.has_PROJ6() should no longer be used"
    )
    stack <- sys.calls()
    stack <- lapply(as.list(stack), function(x) as.character(deparse(x)))
    stop(paste0(
        "'inla.has_PROJ6' should no longer be used.\n",
        "Call stack for developer debugging:\n",
        paste0(lapply(stack, function(x) substr(x, 1, 70)),
            collapse = "\n"
        )
    ))
}

#' @describeIn inla.has_PROJ6 `r lifecycle::badge("deprecated")`
#' Called to warn about using old PROJ4
#' features even though PROJ6 is available
#' @param fun The name of the calling function

inla.not_for_PROJ6 <- function(fun) {
    lifecycle::deprecate_soft(
        when = "2023.04.19",
        what = "inla.not_for_PROJ6()",
        details = "Should no longer be used."
    )
    stack <- sys.calls()
    stack <- lapply(as.list(stack), function(x) as.character(deparse(x)))
    stop(paste0(
        "'inla.not_for_PROJ6' should no longer be used.\n",
        "Call stack for developer debugging:\n",
        paste0(lapply(stack, function(x) substr(x, 1, 70)),
            collapse = "\n"
        )
    ))
}

#' @describeIn inla.has_PROJ6 `r lifecycle::badge("deprecated")`
#' Called to give an error when
#' calling methods that are only available for PROJ6

inla.not_for_PROJ4 <- function(fun) {
    lifecycle::deprecate_soft(
        when = "2023.04.19",
        what = "inla.not_for_PROJ4()",
        details = "Should no longer be used."
    )
    stack <- sys.calls()
    stack <- lapply(as.list(stack), function(x) as.character(deparse(x)))
    stop(paste0(
        "'inla.not_for_PROJ4' should no longer be used.\n",
        "Call stack for developer debugging:\n",
        paste0(lapply(stack, function(x) substr(x, 1, 70)),
            collapse = "\n"
        )
    ))
}

#' @describeIn inla.has_PROJ6 `r lifecycle::badge("deprecated")`
#' Called to warn about falling back
#' to using old PROJ4 methods when a PROJ6 method hasn't been implemented

inla.fallback_PROJ6 <- function(fun) {
    lifecycle::deprecate_soft(
        when = "2023.04.19",
        what = "inla.fallback_PROJ6()",
        details = "Should no longer be used."
    )
    stack <- sys.calls()
    stack <- lapply(as.list(stack), function(x) as.character(deparse(x)))
    stop(paste0(
        "'inla.fallback_PROJ6' should no longer be used.\n",
        "Call stack for developer debugging:\n",
        paste0(lapply(stack, function(x) substr(x, 1, 70)),
            collapse = "\n"
        )
    ))
}

#' @describeIn inla.has_PROJ6 `r lifecycle::badge("deprecated")`
#' Called to give an error when PROJ6
#' is required but not available

inla.requires_PROJ6 <- function(fun) {
    lifecycle::deprecate_soft(
        when = "2023.04.19",
        what = "inla.requires_PROJ6()",
        details = "Should no longer be used."
    )
    stack <- sys.calls()
    stack <- lapply(as.list(stack), function(x) as.character(deparse(x)))
    stop(paste0(
        "'inla.requires_PROJ6' should no longer be used.\n",
        "Call stack for developer debugging:\n",
        paste0(lapply(stack, function(x) substr(x, 1, 70)),
            collapse = "\n"
        )
    ))
}






#' @title Extract CRS information
#'
#' @description
#' `r lifecycle::badge("deprecated")` Use `fmesher::fm_CRS()` instead.
#'
#' Wrapper for `CRS(projargs)` (PROJ4) and `CRS(wkt)` for `sp::Spatial`
#' objects.
#'
#' This function is a convenience method to workaround PROJ4/PROJ6 differences,
#' and the lack of a crs extraction method for Spatial objects.
#'
#' @param x A `sp::Spatial` object
#' @return A `CRS` object, or NULL if no valid CRS identified
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @examples
#' \dontrun{
#' if (require("sp", quietly = TRUE) && interactive()) {
#'     s <- sp::SpatialPoints(matrix(1:6, 3, 2), proj4string = fmesher::fm_CRS("sphere"))
#'     inla.sp_get_crs(s)
#' }
#' }
#'
#' @export inla.sp_get_crs
inla.sp_get_crs <- function(x) {
    fmesher_deprecate(
        "soft",
        1L,
        when = "2023.06.06",
        what = "inla.sp_get_crs()",
        with = "fmesher::fm_CRS()"
    )
    if (is.null(x)) {
        return(NULL)
    }
    fm_CRS(x)
}








#' @title Handling CRS/WKT
#'
#' @description
#' `r lifecycle::badge("deprecated")` in favour of [fmesher::fm_wkt()] and related
#' methods.
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
#' \dontrun{
#' c1 <- fmesher::fm_crs("globe")
#' inla.crs_get_lengthunit(c1)
#' c2 <- inla.crs_set_lengthunit(c1, "metre")
#' inla.crs_get_lengthunit(c2)
#' }
#'
#' @rdname crs_wkt
#' @name crs_wkt
NULL

#' @describeIn crs_wkt `r lifecycle::badge("deprecated")` in favour of
#' [fmesher::fm_wkt_is_geocent()]
#' @export inla.wkt_is_geocent
inla.wkt_is_geocent <- function(wkt) {
    fmesher_deprecate(
        "soft",
        when = "2023.06.06",
        what = "inla.wkt_is_geocent()",
        with = "fmesher::fm_wkt_is_geocent()"
    )
    fmesher::fm_wkt_is_geocent(wkt)
}

#' @export
#' @describeIn crs_wkt `r lifecycle::badge("deprecated")` in favour of
#' [fmesher::fm_crs_is_geocent()]

inla.crs_is_geocent <- function(crs) {
    fmesher_deprecate(
        "soft",
        1L,
        when = "2023.06.06",
        what = "inla.crs_is_geocent()",
        with = "fmesher::fm_crs_is_geocent()"
    )
    fmesher::fm_crs_is_geocent(crs)
}


#' @describeIn crs_wkt `r lifecycle::badge("deprecated")` in favour of
#' [fmesher::fm_ellipsoid_radius()]
#' @export

inla.wkt_get_ellipsoid_radius <- function(wkt) {
    fmesher_deprecate(
        "soft",
        1L,
        when = "2023.06.06",
        what = "inla.wkt_get_ellipsoid_radius()",
        with = "fmesher::fm_ellipsoid_radius()"
    )
    fmesher::fm_ellipsoid_radius(wkt)
}

#' @describeIn crs_wkt `r lifecycle::badge("deprecated")` in favour of
#' [fmesher::fm_ellipsoid_radius()]
#' @export

inla.crs_get_ellipsoid_radius <- function(crs) {
    fmesher_deprecate(
        "soft",
        1L,
        when = "2023.06.06",
        what = "inla.crs_get_ellipsoid_radius()",
        with = "fmesher::fm_ellipsoid_radius()"
    )
    fmesher::fm_ellipsoid_radius(crs)
}


#' @describeIn crs_wkt `r lifecycle::badge("deprecated")` in favour of
#' [fmesher::fm_wkt_set_ellipsoid_radius()]
#' @param radius numeric
#' @export

inla.wkt_set_ellipsoid_radius <- function(wkt, radius) {
    fmesher_deprecate(
        "soft",
        1L,
        when = "2023.06.06",
        what = "inla.wkt_set_ellipsoid_radius()",
        with = "`fmesher::fm_ellipsoid_radius<-`()"
    )
    fmesher::fm_ellipsoid_radius(wkt) <- radius
    wkt
}

#' @describeIn crs_wkt `r lifecycle::badge("deprecated")` in favour of
#' [fmesher::fm_ellipsoid_radius<-()]
#' @export

inla.crs_set_ellipsoid_radius <- function(crs, radius) {
    fmesher_deprecate(
        "soft",
        1L,
        when = "2023.06.06",
        what = "inla.crs_set_ellipsoid_radius()",
        with = "`fmesher::fm_ellipsoid_radius<-`()"
    )
    fmesher::fm_ellipsoid_radius(crs) <- radius
    crs
}






#' Create a coordinate reference system object
#'
#' @description
#' `r lifecycle::badge("deprecated")` in favour of
#' [fmesher::fm_CRS()]
#'
#' Creates either a CRS object or an inla.CRS object, describing a coordinate
#' reference system.
#'
#' @param \dots Arguments passed on to `fmesher::fm_CRS(...)`.
#' @param args list of named proj4 arguments.
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
#' @seealso [fmesher::fm_crs()], [fmesher::fm_wkt()],
#' [fmesher::fm_crs_is_identical()]
#' @examples
#' if (require("sf")) {
#'     crs1 <- fmesher::fm_crs("longlat_globe")
#'     crs2 <- fmesher::fm_crs("lambert_globe")
#'     crs3 <- fmesher::fm_crs("mollweide_norm")
#'     crs4 <- fmesher::fm_crs("hammer_globe")
#'     crs5 <- fmesher::fm_crs("sphere")
#'     crs6 <- fmesher::fm_crs("globe")
#' }
#' \dontrun{
#' names(inla.wkt_predef())
#' }
#'
#' @export inla.CRS
inla.CRS <- function(..., args = NULL) {
    fmesher_deprecate(
        "soft",
        when = "2023.06.06",
        what = "inla.CRS()",
        with = "fmesher::fm_CRS()"
    )

    crs <- fm_CRS(...)
    if (!is.null(args) && (inherits(crs, "CRS") || fmesher::fm_crs_is_null(crs))) {
        if (fmesher::fm_crs_is_null(crs)) {
            crs <- fm_CRS(
                inla.as.CRSargs.list(args)
            )
        } else {
            crs <- fm_CRS(
                paste0(
                    fmesher::fm_proj4string(crs),
                    " ",
                    inla.as.CRSargs.list(args)
                )
            )
        }
    }
    crs
}

#' @return `inla.wkt_predef` returns a WKT2 string defining a projection
#' @examples
#' \dontrun{
#' names(inla.wkt_predef())
#' }
#' @export
#' @describeIn inla.CRS `r lifecycle::badge("deprecated")` in favour of
#' [fmesher::fm_wkt_predef()]

inla.wkt_predef <- function() {
    fmesher_deprecate(
        "soft",
        1L,
        when = "2023.06.06",
        what = "inla.wkt_predef()",
        with = "fmesher::fm_wkt_predef()"
    )
    fmesher::fm_wkt_predef()
}





#' @title Internal WKT handling
#'
#' @description
#' `r lifecycle::badge("deprecated")` in favour of
#' [fmesher::fm_wkt_as_wkt_tree()].
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
    fmesher_deprecate(
        "soft",
        1L,
        when = "2023.06.06",
        what = "inla.as.wkt_tree.wkt()",
        with = "fmesher::fm_wkt_as_wkt_tree()"
    )
    fmesher::fm_wkt_as_wkt_tree(x, ...)
}

#' @rdname wkt_tree
#' @param pretty logical
#' @export

inla.as.wkt.wkt_tree <- function(x, pretty = FALSE, ...) {
    fmesher_deprecate(
        "soft",
        1L,
        when = "2023.06.06",
        what = "inla.as.wkt.wkt_tree()",
        with = "fmesher::fm_wkt_tree_as_wkt()"
    )
    fmesher::fm_wkt_tree_as_wkt(x, pretty = pretty, ...)
}

#' @param item character vector with item labels identifying a parameter item
#' entry.
#' @param duplicate For items that have more than one match, `duplicate`
#' indicates the index number of the desired version. Default: 1
#' @rdname wkt_tree
#' @export

inla.wkt_tree_get_item <- function(x, item, duplicate = 1) {
    fmesher_deprecate(
        "soft",
        1L,
        when = "2023.06.06",
        what = "inla.wkt_tree_get_item()",
        with = "fmesher::fm_wkt_tree_get_item()"
    )
    fmesher::fm_wkt_tree_get_item(x, item = item, duplicate = duplicate)
}

#' @param item_tree An item tree identifying a parameter item entry
#' @rdname wkt_tree
#' @export

inla.wkt_tree_set_item <- function(x, item_tree, duplicate = 1) {
    fmesher_deprecate(
        "soft",
        1L,
        when = "2023.06.06",
        what = "inla.wkt_tree_set_item()",
        with = "fmesher::fm_wkt_tree_set_item()"
    )
    fmesher::fm_wkt_tree_set_item(x, item_tree = item_tree, duplicate = duplicate)
}



#' Show expanded CRS arguments
#'
#' Wrapper for `sp::CRS` and `inla.CRS` objects to extract the
#' coordinate reference system argument string.
#' 'r lifecycle::badge("deprecated")` in favour of [fmesher::fm_proj4string()],
#' or [fmesher::fm_wkt()] for WKT2 representations.
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
#' if (require("sf") && require("sp") && require("fmesher")) {
#'     crs0 <- fm_CRS("longlat_norm")
#'     p4s <- fm_proj4string(crs0)
#'     lst <- inla.as.list.CRSargs(p4s)
#'     crs1 <- inla.as.CRS.list(lst)
#'     lst$a <- 2
#'     crs2 <- fm_CRS(p4s, args = lst)
#'     print(fm_proj4string(crs0))
#'     print(fm_proj4string(crs1))
#'     print(fm_proj4string(crs2))
#' }
#' @export
#' @rdname CRSargs
inla.CRSargs <- function(x, ...) {
    fmesher_deprecate(
        "soft",
        1L,
        when = "2023.06.06",
        what = "inla.CRSargs()",
        with = "fmesher::fm_proj4string()"
    )

    return(fmesher::fm_proj4string(x))
}


#' @details * `inla.as.CRSargs.list`: CRS proj4 string for name=value pair list
#'
#' @export
#' @rdname CRSargs
inla.as.CRSargs.list <- function(x, ...) {
    paste(
        lapply(
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
    inla.as.list.CRSargs(fmesher::fm_proj4string(x))
}

#' @export
#' @rdname CRSargs
inla.as.CRS.list <- function(x, ...) {
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
#' c1 <- inla.CRS("globe")
#' inla.crs_get_lengthunit(c1)
#' c2 <- inla.crs_set_lengthunit(c1, "km")
#' inla.crs_get_lengthunit(c2)
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
#' @importFrom fmesher fm_length_unit

inla.wkt_get_lengthunit <- function(wkt) {
    fmesher_deprecate(
        "soft",
        1L,
        when = "2023.06.06",
        what = "inla.wkt_get_lengthunit()",
        with = "fmesher::fm_length_unit()"
    )

    fmesher::fm_length_unit(wkt)
}

#' @export
#' @rdname crs_wkt
#' @return For `inla.wkt_set_lengthunit`, a
#' WKT2 string with altered length units.
#' Note that the length unit for the ellipsoid radius is unchanged.
#' @importFrom fmesher `fm_length_unit<-`

inla.wkt_set_lengthunit <- function(wkt, unit, params = NULL) {
    fmesher_deprecate(
        "soft",
        1L,
        when = "2023.06.06",
        what = "inla.wkt_set_lengthunit()",
        with = "`fmesher::fm_length_unit<-`()"
    )

    if (!is.null(params)) {
        warning("Non-NULL 'params' for inla.wkt_set_lengthunit ignored.")
    }
    fmesher::fm_length_unit(wkt) <- unit
    wkt
}

#' @return For `inla.crs_get_wkt`, WKT2 string.
#' @export
#' @rdname crs_wkt

inla.crs_get_wkt <- function(crs) {
    fmesher_deprecate(
        "soft",
        1L,
        when = "2023.06.06",
        what = "inla.crs_get_wkt()",
        with = "fmesher::fm_wkt()"
    )

    fmesher::fm_wkt(crs)
}

#' @return For `inla.crs_get_lengthunit`, a
#' list of length units used in the wkt string, excluding the ellipsoid radius
#' unit. (For legacy PROJ4 code, the raw units from the proj4string are
#' returned, if present.)
#' @export
#' @rdname crs_wkt

inla.crs_get_lengthunit <- function(crs) {
    fmesher_deprecate(
        "soft",
        1L,
        when = "2023.06.06",
        what = "inla.crs_get_lengthunit()",
        with = "fmesher::fm_length_unit()"
    )

    fmesher::fm_length_unit(crs)
}

#' @return For `inla.crs_set_lengthunit`, a `sp::CRS` object with
#' altered length units.
#' Note that the length unit for the ellipsoid radius is unchanged.
#' @export
#' @rdname crs_wkt

inla.crs_set_lengthunit <- function(crs, unit, params = NULL) {
    fmesher_deprecate(
        "soft",
        1L,
        when = "2023.06.06",
        what = "inla.crs_set_lengthunit()",
        with = "`fmesher::fm_length_unit<-`()"
    )

    if (!is.null(params)) {
        warning("Non-NULL 'params' for inla.crs_set_lengthunit ignored.")
    }
    fmesher::fm_length_unit(crs) <- unit
    crs
}





inla.wkt_tree_projection_type <- function(wt) {
    fmesher_deprecate(
        "soft",
        1L,
        when = "23.06.06",
        what = "inla.wkt_tree_projection_type()",
        with = "fmesher::fm_wkt_tree_projection_type()"
    )

    return(fmesher::fm_wkt_tree_projection_type(wt))
}

inla.wkt_projection_type <- function(wkt) {
    fmesher_deprecate(
        "soft",
        1L,
        when = "23.06.06",
        what = "inla.wkt_projection_type()",
        with = "fmesher::fm_wkt_projection_type()"
    )

    return(fmesher::fm_wkt_projection_type(wkt))
}

inla.crs_projection_type <- function(crs) {
    fmesher_deprecate(
        "soft",
        1L,
        when = "23.06.06",
        what = "inla.crs_projection_type()",
        with = "fmesher::fm_crs_projection_type()"
    )

    return(fmesher::fm_crs_projection_type(crs))
}

## +proj=longlat in (-180,180)x(-90,90)
## +proj=moll in (-2,2)x(-1,1) scaled by +a and +b, and +units
## +proj=lambert in (-pi,pi)x(-1,1) scaled by +a and +b, and +units
inla.crs.bounds <- function(crs, warn.unknown = FALSE) {
    fmesher_deprecate(
        "soft",
        2L,
        when = "23.08.18",
        what = "inla.crs.bounds()",
        with = "fmesher::fm_crs_bounds()"
    )
    return(fmesher::fm_crs_bounds(crs, warn.unknown = warn.unknown))
}







#' @title Test CRS and inla.CRS for equality
#'
#' @description
#' `r lifecycle::badge("deprecated")` Use [fmesher::fm_crs_is_identical()]
#' instead.
#'
#' Wrapper for identical, optionally testing only the CRS part of two objects
#' Deprecated in favour of [fmesher::fm_crs_is_identical()]
#'
#' @param ... Arguments passed on to [fmesher::fm_crs_is_identical()]
#'
#' @export inla.identical.CRS
inla.identical.CRS <- function(...) {
    fmesher_deprecate(
        "soft",
        1L,
        when = "2023.06.06",
        what = "inla.identical.CRS()",
        with = "fmesher::fm_crs_is_identical()"
    )

    fmesher::fm_crs_is_identical(...)
}






#' @title Wrapper method for `fmesher::fm_transform`
#'
#' @description
#' `r lifecycle::badge("deprecated")` in favour of [fmesher::fm_transform()].
#'
#' Handles transformation of various inla objects according to coordinate
#' reference systems of `sf::crs`, `sp::CRS` or `inla.CRS` class.
#'
#'
#' @param x The object that should be transformed from it's current CRS to a
#' new CRS
#' @param CRSobj passed on as the `crs` argument to [fmesher::fm_transform()].
#' @param \dots Potential other arguments for [fmesher::fm_transform()].
#' @return The object is returned with its coordinates transformed
#' @author Finn Lindgren <finn.lindgren@@gmail.com>
#' @seealso [inla.CRS()]
#' @examples
#'
#' if (require("sf") && require("sp") && require("fmesher")) {
#'     latt <- inla.mesh.lattice(-10:10, 40:60)
#'     mesh1 <- inla.mesh.create(
#'         lattice = latt, extend = FALSE, refine = FALSE,
#'         crs = fm_CRS("longlat_norm")
#'     )
#'     mesh2 <- fm_transform(mesh1, fm_crs("lambert_globe"))
#'     print(summary(mesh1))
#'     print(summary(mesh2))
#' }
#' @export inla.spTransform
inla.spTransform <- function(x, CRSobj, ...) {
    fmesher_deprecate(
        "soft",
        1L,
        when = "2023.06.06",
        what = "inla.spTransform()",
        with = "fmesher::fm_transform()"
    )

    fmesher::fm_transform(x, crs = CRSobj, ...)
}




inla.crs_detect_manifold <- function(crs) {
    fmesher_deprecate(
        "soft",
        1L,
        when = "2023.06.06",
        what = "inla.crs_detect_manifold()",
        with = "fmesher::fm_crs_detect_manifold()"
    )

    fmesher::fm_crs_detect_manifold(crs)
}





#' @title Convert `sp` objects to `inla.mesh.segment` objects.
#'
#' @description `r lifecycle::badge("superseded")` by [fmesher::fm_as_segm()]
#'
#' @param sp An `sp` polygon object of class `Polygon`,
#' `Polygons`, `SpatialPolygons`, or `SpatialPolygonsDataFrame`.
#' @param \dots Additional arguments passed on to `fmesher::fm_as_segm()`.
#' @return A [inla.mesh.segment()] object, or a list of
#' [inla.mesh.segment()] objects.
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @seealso [inla.mesh.segment()]
#' @export as.inla.mesh.segment
#' @importFrom fmesher fm_as_segm
as.inla.mesh.segment <-
    function(sp, ...) {
        fmesher_deprecate(
            "soft",
            1L,
            when = "2023.06.06",
            what = "as.inla.mesh.segment()",
            with = "fmesher::fm_as_segm()"
        )

        fmesher::fm_as_segm(sp, ...)
    }

#' @export
#' @describeIn as.inla.mesh.segment `r lifecycle::badge("superseded")`
#' by [fmesher::fm_as_segm()]
inla.sp2segment <-
    function(sp, ...) {
        fmesher_deprecate(
            "soft",
            1L,
            when = "2023.06.06",
            what = "inla.sp2segment()",
            with = "fmesher::fm_as_segm()"
        )

        fmesher::fm_as_segm(sp, ...)
    }
