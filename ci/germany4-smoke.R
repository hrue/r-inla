### Germany4-style smoke test with a synthetic spatial graph.
### Fits the same model classes as the Germany4 multivariate disease
### mapping example (a scaled besag per disease, then a treepcor cgeneric
### kronecker generic0 joint model) but generates a random connected
### graph and Poisson data, so it runs anywhere with no data files.
### The inla binary under test comes from the INLA_CALL environment
### variable; unset, the installed package default is used.

library(graphpcor)
library(INLA)
library(Matrix)

inla.call <- Sys.getenv("INLA_CALL", unset = NA)
if (!is.na(inla.call)) {
    stopifnot(file.exists(inla.call))
    inla.setOption(inla.call = inla.call)
}
inla.setOption(smtp = "stiles")
cat("== inla.call:", inla.getOption("inla.call"), "\n")
cat("== smtp:     ", inla.getOption("smtp"), "\n")

set.seed(1)
K <- 4
n.areas <- 60

## random connected undirected graph: a chain (keeps it connected, which
## scale.model needs) plus random extra edges for realistic degree spread
ii <- c(1:(n.areas - 1), sample(1:n.areas, 2 * n.areas, replace = TRUE))
jj <- c(2:n.areas,       sample(1:n.areas, 2 * n.areas, replace = TRUE))
keep <- ii != jj
ii2 <- pmin(ii[keep], jj[keep]); jj2 <- pmax(ii[keep], jj[keep])
graphG <- sparseMatrix(i = c(ii2, jj2), j = c(jj2, ii2), x = 1,
                       dims = c(n.areas, n.areas))
graphG@x[] <- 1   ## dedupe multi-edges back to 0/1
cat("== graph: n =", n.areas, ", edges =", sum(graphG) / 2, "\n")

vnames <- c("oral", "osph", "lary", "lung")

## synthetic observed/expected counts per disease
E.all <- runif(K * n.areas, 0.5, 3)
eta   <- rnorm(K * n.areas, 0, 0.3)
Obs.all <- rpois(K * n.areas, E.all * exp(eta))

### Model 1: scaled besag for each disease
pcprec <- list(prior = "pc.prec", param = c(1, 0.01))
m1f <- Obs ~ f(id.area, model = "besag", graph = graphG,
               scale.model = TRUE, hyper = list(theta = pcprec))
lres.m1 <- lapply(1:K, function(k) {
    idx <- (k - 1) * n.areas + 1:n.areas
    r <- inla(formula = m1f, family = "poisson",
              data = data.frame(id.area = 1:n.areas,
                                Obs = Obs.all[idx],
                                expected = E.all[idx]),
              E = expected)
    stopifnot(is.finite(r$mlik[1]))
    r$mlik[1]
})
cat("== model 1 (besag x", K, "diseases) mlik:",
    sprintf("%.2f", unlist(lres.m1)), "\n")

### Model 2: treepcor cgeneric + generic0 kronecker (fails with
### 'inla_cgeneric_kronecker not found' on binaries older than the
### cgeneric mapper that introduced it)
tree3 <- treepcor(p1 ~ p2 + c4,
                  p2 ~ p3 + c2 + c3,
                  p3 ~ c1)
tree3model <- cgeneric(model = tree3, lambda = 5,
                       sigma.prior.reference = c(1, 1, 1, 1),
                       sigma.prior.probability = c(0.1, 0.1, 0.1, 0.1))

R.spatial <- inla.as.sparse(
    Diagonal(n.areas, x = colSums(graphG)) - graphG)
cBesag <- cgeneric(model = "generic0", R = R.spatial, param = c(1, NA))
tree3Besag <- kronecker(tree3model, cBesag)

ldata <- data.frame(
    iddata = 1:(K * n.areas),
    disease = factor(rep(vnames, each = n.areas), levels = vnames),
    Obs = Obs.all,
    E = E.all)

m2f <- Obs ~ 0 + disease + f(iddata, model = tree3Besag)
res.m2 <- inla(formula = m2f, family = "poisson",
               data = ldata, E = E)

stopifnot(is.finite(res.m2$mlik[1]))
cat("== model 2 (treepcor kronecker besag) mlik:",
    sprintf("%.2f", res.m2$mlik[1]), "\n")
cat("== fixed effects:\n")
print(round(res.m2$summary.fixed[, 1:2], 3))
cat("GERMANY4 SMOKE OK\n")
