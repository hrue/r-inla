## The model checks, shared by every platform.
##
## This used to be a heredoc inside ci/smoke-test.sh, which meant Windows
## could not run it: that lane never calls smoke-test.sh, so it had no
## cgeneric coverage at all and only an ad-hoc rgeneric check. Keeping the
## checks in one file instead of two means Windows cannot silently drift into
## testing something weaker than Linux and macOS.
##
## Usage:  Rscript --vanilla ci/smoke-models.R <inla-binary> [testmodel-dir] [smtp]
argv <- commandArgs(trailingOnly = TRUE)
library(INLA)
inla.setOption(inla.call = argv[1])
inla.setOption(num.threads = "2:1")

## INLA_SMTP names a second sparse-matrix backend to check. The models
## below then run TWICE -- once on the default backend, once on this one --
## and the two marginal likelihoods must agree. Linking against a solver
## proves nothing about whether it computes the same posterior.
smtp <- if (length(argv) >= 3 && nzchar(argv[3])) argv[3] else ""

set.seed(1)
n   <- 300
x   <- rnorm(n)
idx <- 1:n
y   <- 1 + 0.7 * x + rnorm(n, sd = 0.3) + rnorm(n, sd = 0.2)

r <- inla(y ~ x + f(idx, model = "iid"),
          data   = data.frame(y = y, x = x, idx = idx),
          family = "gaussian")

stopifnot(is.finite(r$mlik[1]))
b <- r$summary.fixed
stopifnot(abs(b["(Intercept)", "mean"] - 1.0) < 0.1)
stopifnot(abs(b["x", "mean"] - 0.7) < 0.1)
cat("smoke test OK: mlik =", r$mlik[1], "\n")

## The same model through the alternative backend. Both are exact
## factorizations of the same precision matrix, so the marginal likelihoods
## differ only by arithmetic ordering -- a mismatch beyond that means the
## backend is wrong, not merely different.
if (nzchar(smtp)) {
    cat("== backend:", smtp, "==\n")
    inla.setOption(smtp = smtp)
    rs <- inla(y ~ x + f(idx, model = "iid"),
               data   = data.frame(y = y, x = x, idx = idx),
               family = "gaussian")
    stopifnot(is.finite(rs$mlik[1]))
    ## The backends factorize the same matrices but sit under different BLAS
    ## libraries, and INLA re-optimizes per run -- tiny arithmetic
    ## differences move the optimizer's stopping point, which shows up as
    ## ~1e-2 in mlik. The bar here is "same posterior", not "same bits":
    ## the parameter estimates must recover the truth exactly as the
    ## default-backend run does, and mlik must agree to well under anything
    ## a real defect (wrong logdet, wrong solve) would produce.
    bs <- rs$summary.fixed
    stopifnot(abs(bs["(Intercept)", "mean"] - 1.0) < 0.1)
    stopifnot(abs(bs["x", "mean"] - 0.7) < 0.1)
    d <- abs(rs$mlik[1] - r$mlik[1])
    cat(sprintf("mlik default=%.6f  %s=%.6f  |diff|=%.2e\n",
                r$mlik[1], smtp, rs$mlik[1], d))
    if (d > 5e-2)
        stop(sprintf("%s disagrees with the default backend by %.3e", smtp, d))
    cat(smtp, "agrees with the default backend\n")
}

## An rgeneric model: the only part of the binary that calls back INTO R
## while it runs. It is worth its own check because nothing else here
## exercises that path, and in a build with INLA_WITH_LIBR_DLOPEN it also
## proves the runtime library loading -- the binary locates the libR of
## whatever R is driving it and evaluates the model's R functions through
## it. A plain "the model converged" result is not evidence of that.
cat("== rgeneric ==\n")
rg <- inla.rgeneric.define(inla.rgeneric.iid.model, n = 30)
set.seed(2)
yy  <- rnorm(30)
iid <- 1:30
rr  <- inla(yy ~ -1 + f(iid, model = rg),
            data = data.frame(yy = yy, iid = iid),
            family = "gaussian")
stopifnot(is.finite(rr$mlik[1]))
stopifnot(nrow(rr$summary.random$iid) == 30)
cat("rgeneric OK: mlik =", rr$mlik[1], "\n")

## Built-in cgeneric models. The model functions live in the BINARY
## (cgeneric-mapper: the shlib argument is ignored for registered names), so
## this is the only part of the suite that exercises the external-packages
## objects and the mapper. Under INLA_SMTP each model runs on both backends
## and the marginal likelihoods must agree -- the cgeneric path hands the
## backend a precision assembled per iteration, which nothing above covers.
check_both <- function(fit, label, tol = 5e-2) {
    old <- inla.getOption("smtp"); on.exit(inla.setOption(smtp = old), add = TRUE)
    r1 <- fit()
    stopifnot(is.finite(r1$mlik[1]))
    cat(sprintf("%s: mlik(default) = %.6f\n", label, r1$mlik[1]))
    if (nzchar(smtp)) {
        inla.setOption(smtp = smtp)
        r2 <- fit()
        stopifnot(is.finite(r2$mlik[1]))
        d <- abs(r2$mlik[1] - r1$mlik[1])
        cat(sprintf("%s: mlik(%s) = %.6f  |diff| = %.2e\n", label, smtp, r2$mlik[1], d))
        if (d > tol)
            stop(sprintf("%s: %s disagrees with the default backend by %.3e",
                         label, smtp, d))
    }
    invisible(r1)
}

## fbesag on a chain graph. The input is the intrinsic PRECISION (D - A):
## get_fbesag scales it via inla.scale.model, which aborts on an adjacency
## (indefinite matrix) -- the demo's all-ones matrix is equally degenerate.
cat("== cgeneric: fbesag ==\n")
library(fbesag)
nf  <- 20
Af  <- matrix(0, nf, nf)
for (i in 1:(nf - 1)) { Af[i, i + 1] <- 1; Af[i + 1, i] <- 1 }
Qf  <- diag(rowSums(Af)) - Af
## The checkout-installed INLA package ships no precompiled external libs,
## so point fbesag at its OWN installed shlib (R CMD INSTALL compiled it from
## the same src/fbesag.c that is baked into the binary). The binary resolves
## the registered model name from its cgeneric table either way; the file
## only has to EXIST for inla.cgeneric.define's R-side check.
fbso <- system.file("libs", paste0("fbesag", .Platform$dynlib.ext), package = "fbesag")
stopifnot(nzchar(fbso))
mfb <- fbesag::get_fbesag(graph = Qf, id = rep(1:2, each = nf / 2),
                          sd_gamma = 0.15, param = list(p1 = 1, p2 = 1e-5),
                          useINLAprecomp = FALSE, libpath = fbso)
set.seed(3)
yf <- rnorm(nf, 0, 0.1)
check_both(function()
    inla(y ~ 1 + f(idx, model = mfb),
         data = list(y = yf, idx = 1:nf), family = "gaussian"),
    "fbesag")

## INLAspacetime model 102 (sstspde, also compiled in): a small space-time
## field on a coarse mesh -- the heaviest precision structure in this suite.
cat("== cgeneric: spacetime (model 102) ==\n")
library(INLAspacetime)
library(fmesher)
smesh <- fm_mesh_2d(cbind(c(0, 1, 0, 1, 0.5), c(0, 0, 1, 1, 0.5)),
                    max.edge = 0.6, offset = 0.3)
tmesh <- fm_mesh_1d(1:4)
## Pass the CRAN package's own compiled shlib explicitly: the checkout-
## installed INLA ships no precompiled external libs, and the shorter
## useINLAprecomp = FALSE route trips an upstream bug (its branch never
## assigns 'hasverbose'). The binary resolves model 102 from its cgeneric
## table regardless; the file only satisfies the R-side existence check.
stlib <- system.file("libs", paste0("INLAspacetime", .Platform$dynlib.ext),
                     package = "INLAspacetime")
stopifnot(nzchar(stlib))
stm <- stModel.define(smesh, tmesh, model = "102",
                      control.priors = list(prs    = c(0.5, 0.5),
                                            prt    = c(2, 0.5),
                                            psigma = c(1, 0.5)),
                      libpath = stlib)
set.seed(4)
nst <- 60
loc <- cbind(runif(nst), runif(nst))
tim <- sample(1:4, nst, replace = TRUE)
Ast <- inla.spde.make.A(smesh, loc = loc, group = tim, group.mesh = tmesh)
yst <- rnorm(nst, 0, 0.3)
stk <- inla.stack(data = list(y = yst), A = list(Ast, 1),
                  effects = list(field = 1:(smesh$n * tmesh$n),
                                 intercept = rep(1, nst)))
check_both(function()
    inla(y ~ -1 + intercept + f(field, model = stm),
         data = inla.stack.data(stk),
         control.predictor = list(A = inla.stack.A(stk)),
         family = "gaussian"),
    "spacetime-102")
cat("cgeneric OK\n")

outdir <- argv[2]
if (!is.na(outdir) && nzchar(outdir)) {
    unlink(outdir, recursive = TRUE)
    r2 <- inla(y ~ x + f(idx, model = "iid"),
               data = data.frame(y = y, x = x, idx = idx),
               family = "gaussian",
               keep = TRUE, working.directory = outdir)
    ## the files may sit in a subdirectory of the working directory
    ini <- list.files(outdir, pattern = "^Model\\.ini$",
                      recursive = TRUE, full.names = TRUE)
    if (length(ini) == 0) {
        cat("contents of", outdir, ":\n")
        print(list.files(outdir, recursive = TRUE))
        stop("no Model.ini found under ", outdir)
    }
    ## make the export relocatable: the ini references its data files by
    ## absolute path, which only exists on this machine
    dir  <- dirname(normalizePath(ini[1]))
    txt  <- readLines(ini[1])
    txt  <- gsub(paste0(dir, "/"), "./", txt, fixed = TRUE)
    writeLines(txt, ini[1])
    cat("test model exported:", ini[1], "\n")
}
