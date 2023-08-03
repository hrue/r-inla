## Nothing to Export

`inla.fmesher.call` <- function(
                                fmesher.call = inla.getOption("fmesher.call"),
                                all.args, prefix,
                                timeout = inla.getOption("fmesher.timeout")) {
    # Change to _stop once evolution status 2L has eliminated
    # all standalone fmesher binary calls
    fmesher_deprecate_warn(2L,
                           "23.08.03",
                           "inla.fmesher.call()",
                           details = "With `fmesher.evolution = 2L`, no standalone fmesher binary calls should be made.",
                           always = TRUE)
    
    if (inla.os("linux") || inla.os("mac") || inla.os("mac.arm64")) {
        echoc <- system(paste(shQuote(fmesher.call), all.args, shQuote(prefix)),
            timeout = ceiling(timeout)
        )
    }
    else if (inla.os("windows")) {
        echoc <- try(system(paste(
            shQuote(fmesher.call),
            all.args,
            shQuote(prefix)
        ),
        timeout = ceiling(timeout)
        ),
        silent = TRUE
        )
        echoc <- 0
    }
    else {
          stop("\n\tNot supported architecture.")
      }
}

`inla.fmesher.make.prefix` <- function(dir = tempdir(), prefix = NULL) {
    if (is.null(prefix)) {
        if (is.null(dir)) {
              prefix <- paste(tempfile(), ".", sep = "")
          } else {
            inla.dir.create(dir)
            prefix <- paste(tempfile(tmpdir = dir), ".", sep = "")
        }
    } else {
        inla.dir.create(inla.dirname(prefix))
    }
    return(prefix)
}
