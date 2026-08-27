### Package setup for the Germany4 smoke test.
### Kept in a file rather than inline `Rscript -e`: a multi-line -e
### argument crashed R outright on the Windows runner (exit 139, no
### output at all), and a script file sidesteps that shell quoting path.

cat("R:", R.version.string, "\n")
cat("platform:", R.version$platform, "\n")

repos <- c(inla = "https://inla.r-inla-download.org/R/testing",
           cran = "https://cloud.r-project.org")

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

## Report both places the cgeneric library could come from, then require
## that at least one of them exists: without it the model fails later
## inside cgeneric_get() with a dlopen error, far from the cause.
own <- list.files(file.path(find.package("graphpcor"), "libs"))
ext <- tryCatch(
    list.files(file.path(find.package("INLA"), "bin",
                         if (.Platform$OS.type == "windows") "windows"
                         else if (Sys.info()[["sysname"]] == "Darwin")
                             if (R.version$arch == "aarch64") "mac.arm64" else "mac"
                         else "linux",
                         "64bit", "external", "graphpcor")),
    error = function(e) character(0))
cat("graphpcor own libs:", if (length(own)) own else "NONE", "\n")
cat("INLA external graphpcor:", if (length(ext)) ext else "NONE", "\n")
if (!length(own) && !length(ext))
    stop("no cgeneric library for graphpcor on this platform")
