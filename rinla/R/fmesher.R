## Nothing to Export

`inla.fmesher.call` <- function(
                                fmesher.call = inla.getOption("fmesher.call"),
                                all.args, prefix,
                                timeout = inla.getOption("fmesher.timeout")) {
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
