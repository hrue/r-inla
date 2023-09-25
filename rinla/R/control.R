#' @title inla control argument handlers
#' @description
#' Methods for controlling [inla()] and [f()] `control` arguments
#' 
#' @name inla-control
#' @rdname inla-control
#' @keywords internal
#' @export
ctrl_class <- function(x) {
  cl <- class(x)
  cl[min(grep("^ctrl_", cl))]
}

#' @describeIn inla-control Returns the type-part of the control class name, i.e.
#' without the "ctrl_" prefix.
#' @export
ctrl_type <- function(x) {
  UseMethod("ctrl_type")
}
#' @rdname inla-control
#' @export
ctrl_type.default <- function(x) {
  sub("^ctrl_", replacement = "", x = ctrl_class(x))
}
#' @rdname inla-control
#' @export
ctrl_type.character <- function(x) {
  sub("^ctrl_|^control\\.", replacement = "", x = x)
}

#' @rdname inla-control
#' @export
ctrl_object <- function(x, the_type) {
  structure(if (is.null(x)) list() else x, class = paste0("ctrl_", the_type))
}

#' @rdname inla-control
#' @export
ctrl_defaults <- function(x) {
  UseMethod("ctrl_defaults")
}

ctrl_defaults_call <- function(x) {
  the_call <- sub(pattern = "^ctrl_", replacement = "control.", x = x)
  the_call <- gsub(pattern = "_", replacement = ".", x = the_call)
  the_call
}

#' @rdname inla-control
#' @export
ctrl_defaults.character <- function(x) {
  do.call(ctrl_defaults_call(x), list())
}

#' @rdname inla-control
#' @export
ctrl_defaults.default <- function(x) {
  do.call(ctrl_defaults_call(ctrl_class(x)), list())
}

#' @rdname inla-control
#' @export
ctrl_update <- function(x, ...) {
  UseMethod("ctrl_update")
}

#' @describeIn inla-control Merges the `x` object with the defaults for control
#' type `ctrl_type(x)`. Recursively handles `control.*` elements, with `NULL`
#' defaults also being replaced by their corresponding defaults before merging.
#' @export
ctrl_update.default <- function(x, ...) {
  def <- ctrl_defaults(x)
  def_names <- names(def)
  x_names <- names(x)
  if (any(!(x_names %in% def_names))) {
    warning(paste0("Control name ", paste0(
      "'",
      x_names[!(x_names %in% def_names)],
      "'",
      collapse = ", "
    ), " appears to be invalid for ", ctrl_class(x), "."),
    immediate. = TRUE)
  }
  x_names <- intersect(x_names, def_names)
  # Handle sub-controls
  sub.controls <- def_names[grep(pattern = "^control\\.", x = def_names)]
  for (nm in sub.controls) {
    the_type <- ctrl_type(nm)
    # If NULL in def, remove if present but NULL in x
    if (nm %in% x_names) {
      def[[nm]] <- x[[nm]]
    }
    if (!is.null(def[[nm]])) {
      def[[nm]] <- ctrl_update(ctrl_object(def[[nm]], the_type))
    }
  }
  # Copy over everything else
  copy_names <- setdiff(x_names, sub.controls)
  def[copy_names] <- x[copy_names]

  def
}


ctrl_handle_hyper <- function(x, model, section) {
  x$hyper <- inla.set.hyper(
    model,
    section,
    x[["hyper"]],
    x[["initial"]],
    x[["fixed"]],
    x[["prior"]],
    x[["param"]]
  )
  x
}

#' @rdname inla-control
#' @export
ctrl_update.ctrl_family <- function(x, model, ...) {
  y <- NextMethod()
  y <- ctrl_handle_hyper(y, model = model, section = "likelihood")
  y
}

#' @rdname inla-control
#' @export
ctrl_update.ctrl_mix <- function(x, ...) {
  y <- NextMethod()
  y <- ctrl_handle_hyper(y, y[["model"]], "mix")
  y
}

#' @rdname inla-control
#' @export
ctrl_update.ctrl_link <- function(x, ...) {
  y <- NextMethod()
  y <- ctrl_handle_hyper(y, y[["model"]], "link")
  y
}

#' @rdname inla-control
#' @export
ctrl_update.ctrl_lp_scale <- function(x, ...) {
  y <- NextMethod()
  y <- ctrl_handle_hyper(y, "lp.scale", "lp.scale")
  y
}

#' @rdname inla-control
#' @export
ctrl_update.ctrl_predictor <- function(x, ...) {
  y <- NextMethod()
  y <- ctrl_handle_hyper(y, "predictor", "predictor")
  y
}

#' @describeIn inla-control If `residuals is `TRUE`, also sets `dic = TRUE`.
#' @export
ctrl_update.ctrl_compute <- function(x, ...) {
  y <- NextMethod()
  if (isTRUE(y[["residuals"]])) {
    y$dic <- TRUE
  }
  y
}



#' @rdname inla-control
#' @export
ctrl_update.ctrl_group <- function(x, ...) {
  y <- NextMethod()
  y <- ctrl_handle_hyper(y, y[["model"]], "group")
  y
}

#' @rdname inla-control
#' @export
ctrl_update.ctrl_scopy <- function(x, ...) {
  y <- NextMethod()
  y <- ctrl_handle_hyper(y, y[["model"]], "scopy")
  y
}

#' @rdname inla-control
#' @export
ctrl_update.ctrl_hazard <- function(x, ...) {
  y <- NextMethod()
  y <- ctrl_handle_hyper(y, y[["model"]], "hazard")
  y
}
