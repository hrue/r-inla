### Package setup for the Germany4 smoke test.
### Kept in a file rather than inline `Rscript -e`: a multi-line -e
### argument crashed R outright on the Windows runner (exit 139, no
### output at all), and a script file sidesteps that shell quoting path.

cat("R:", R.version.string, "\n")
cat("platform:", R.version$platform, "\n")

## ADD the INLA repo to whatever the runner already has; do not replace it.
## setup-r points CRAN at the Posit package manager, which serves prebuilt
## binaries for this platform. Overriding that with cloud.r-project.org
## forced SOURCE installs of every dependency on Linux, and the heavy ones
## (sf needs GDAL/GEOS/PROJ, fmesher needs compiling) failed for want of
## system libraries this test has no reason to need.
## The Windows runner leaves CRAN as the unset placeholder "@CRAN@", which
## install.packages rejects ("trying to use CRAN without setting a mirror"),
## so fill that one in rather than replacing the whole option.
repos <- getOption("repos")
if (!length(repos) || is.na(repos["CRAN"]) || !nzchar(repos["CRAN"]) ||
    identical(unname(repos["CRAN"]), "@CRAN@"))
    repos["CRAN"] <- "https://cloud.r-project.org"
repos <- c(inla = "https://inla.r-inla-download.org/R/testing", repos)
cat("repos:\n"); print(repos)

install.packages(c("remotes", "INLAtools"), repos = repos)

## INLA: binary if this platform has one in the testing repo, else source
tryCatch(install.packages("INLA", repos = repos),
         error = function(e) message("binary INLA install failed: ", conditionMessage(e)))
if (!requireNamespace("INLA", quietly = TRUE))
    install.packages("INLA", repos = repos, type = "source")

## graphpcor from GitHub, not CRAN, because this test calls into its
## cgeneric library FROM R and that is the one path with no macOS copy.
##
## Two different mechanisms are easy to conflate. FITTING a registered
## cgeneric model needs nothing on disk: the model functions are compiled
## into the inla binary (external-packages -> cgeneric-table.h) and the
## mapper ignores the shlib argument for registered names, which is why
## the mac binary carries inla_cgeneric_{graphpcor,treepcor,kronecker}
## and the board's fbesag test passes there. But graphpcor also evaluates
## its models in R (cgeneric_get, used here to compose the kronecker),
## and THAT dlopens a real shared library, which INLA ships under
## bin/<platform>/64bit/external/. Those files come and go: INLA 26.08.07
## carried five of them (fbesag, graphpcor, INLAspacetime, INLAtools,
## rSPDE) in both the source package and the macOS binary, while in
## 26.08.22 the source and macOS packages carry NONE and only the Windows
## zip still has them. CRAN's macOS binary of graphpcor ships no compiled
## library either, so on 26.08.22 there is nothing to load on macOS at
## all. Building graphpcor from its own sources gives every platform its
## own libs/graphpcor.so, so this test keeps working either way.
remotes::install_github("eliaskrainski/graphpcor", upgrade = "never")

cat("INLA:", as.character(packageVersion("INLA")), "\n")
cat("graphpcor:", as.character(packageVersion("graphpcor")), "\n")

## Populate INLA's external/ tree from the packages' own libraries.
##
## Resolving the library in R is not enough. The R side falls back to the
## package's own libs/ and warns "Changed `shlib` to ...", but the model
## object still travels to the inla PROGRAM carrying the original path,
## and inla dlopens it there: on a 26.08.22 macOS install that aborts the
## fit with "Failed to load shared library 'inla_cgeneric_generic0':
## dlopen(.../external/INLAtools/libINLAtools.so) (no such file)".
##
## So put a copy where inla looks. inla.external.lib() defines that
## location, and note it spells the file lib<pkg>.so on EVERY platform
## (macOS included: the 26.08.07 mac package shipped Mach-O objects named
## .so), so ask it rather than building the path by hand. Only fills gaps:
## where INLA already ships the library (Windows) nothing is touched.
library(INLA)
for (p in c("INLAtools", "graphpcor")) {
    own <- list.files(file.path(find.package(p), "libs"),
                      pattern = "\\.(so|dll|dylib)$", full.names = TRUE)
    target <- eval(call("inla.external.lib", as.name(p)))
    if (file.exists(target)) {
        cat("external ", p, ": already shipped by INLA\n", sep = "")
    } else if (length(own)) {
        dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
        ok <- file.copy(own[1], target, overwrite = TRUE)
        cat("external ", p, ": installed from ", own[1], " -> ", target,
            if (ok) "" else "  (COPY FAILED)", "\n", sep = "")
    } else {
        cat("external ", p, ": NO library available\n", sep = "")
    }
    if (!file.exists(target))
        stop("no cgeneric library for ", p, " on this platform")
}
