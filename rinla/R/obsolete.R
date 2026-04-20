## Nothing to export

## obsolete functions goes here.

`inla.obsolete` <- function(old, new) {
    stop(inla.paste(c(
        "INLA-function `", as.character(old)[1], "' is obsolete. Use function `",
        as.character(new), "' instead."
    ), sep = ""))
}
