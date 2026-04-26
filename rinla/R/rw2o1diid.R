#' 1D Second-order random walk plus iid noise (BYM2-style parameterisation)
#'
#' Convenience constructor for a 1D second-order random walk mixed with iid
#' noise, using the BYM2-style `(tau, phi)` parameterisation. The latent
#' field has precision
#' \deqn{Q(\tau, \phi) = \tau \, \bigl( \phi \, R + (1-\phi)\, I \bigr),}
#' where `R` is the RW2 structure matrix scaled to generalised variance one
#' and `I` is the identity. The `tau` parameter controls the total marginal
#' precision and `phi` in `(0, 1)` controls the fraction of variance
#' attributable to the structured RW2 component.
#'
#' Both hyperparameters use PC priors. The PC prior on `tau` is the standard
#' [`inla.pc.dprec()`] prior, parameterised by `(u, alpha)` such that
#' `P(1/sqrt(tau) > u) = alpha`. The PC prior on `phi` is computed from the
#' eigenvalues of the scaled RW2 structure matrix using the same KLD-based
#' derivation as `bym2` (see [`inla.pc.bym.phi()`]), parameterised by
#' `(u, alpha)` such that `P(phi < u) = alpha`.
#'
#' Because this model is implemented as an rgeneric, INLA reports the
#' hyperparameter marginals on the internal scale: `theta_1 = log(tau)`
#' (named `"log precision for <label>"`) and `theta_2 = logit(phi)`
#' (named `"logit phi for <label>"`), where `<label>` is the is the variable
#' name passed as the first argument to f() (e.g. "time" for f(time, ...)).
#' To obtain marginals for `tau` and `phi` on the user scale, use
#' [`inla.rw2o1diid.hyperpar()`].
#'
#' @param n Integer. Length of the chain. Must be `>= 5` (the minimum size
#'     for which RW2 is well-defined).
#' @param prior.tau Named list with entries `u` and `alpha` for the PC prior
#'     on the precision `tau`. Default `list(u = 1, alpha = 0.01)`.
#' @param prior.phi Named list with entries `u` and `alpha` for the PC prior
#'     on the mixing parameter `phi`. Default `list(u = 0.5, alpha = 0.5)`,
#'     matching the `bym2` default.
#' @param debug Logical. Passed to [`inla.rgeneric.define()`].
#' @return A model-specification object that can be passed as the `model`
#'     argument to [`f()`].
#' @author Antonio R. Vargas
#' @seealso [`f()`], [`inla.rw2o1diid.hyperpar()`],
#'     [`inla.rgeneric.define()`], `inla.doc("bym2")`
#' @examples
#' \dontrun{
#'   n <- 100
#'   time <- 1:n
#'   x <- cumsum(rnorm(n, sd = 0.3)) + rnorm(n, sd = 1)
#'
#'   ## defaults: PC prior on tau with (u, alpha) = (1, 0.01)
#'   ## and on phi with (u, alpha) = (0.5, 0.5)
#'   r <- inla(
#'     y ~ f(time, model = inla.rw2o1diid(n)) - 1,
#'     data = data.frame(y, time)
#'   )
#'
#'   ## custom priors
#'   r <- inla(
#'     y ~ f(
#'       time,
#'       model = inla.rw2o1diid(
#'         n,
#'         prior.tau = list(u = 2, alpha = 0.05),
#'         prior.phi = list(u = 0.7, alpha = 0.3)
#'       )
#'     ) - 1,
#'     data = data.frame(y, time)
#'   )
#'
#'   ## user-scale tau and phi marginals
#'   hp <- inla.rw2o1diid.hyperpar(r, "time")
#'   inla.zmarginal(hp$tau)
#'   inla.zmarginal(hp$phi)
#' }
#' @rdname rw2o1diid
#' @export
`inla.rw2o1diid` <- function(
  n,
  prior.tau = list(u = 1, alpha = 0.01),
  prior.phi = list(u = 0.5, alpha = 0.5),
  debug = FALSE
) {
  stopifnot(is.numeric(n), length(n) == 1, n >= 5)
  stopifnot(is.list(prior.tau), all(c("u", "alpha") %in% names(prior.tau)))
  stopifnot(is.list(prior.phi), all(c("u", "alpha") %in% names(prior.phi)))
  n <- as.integer(n)

  ####
  # Build the phi PC log-prior

  # Calculate variance components of the random walk
  R_scaled <- inla.rw(n, order = 2, scale.model = TRUE, sparse = TRUE)
  R_dense <- as.matrix(R_scaled)
  eig_R <- pmax(0, eigen(R_dense, symmetric = TRUE, only.values = TRUE)$values)
  marg_var <- diag(inla.ginv(R_dense, rankdef = 2))

  # - inla.pc.bym.phi() with return.as.table = FALSE returns a cubic-spline
  #   closure giving log-density on the external phi scale
  # - local() trims the captured environment to just the spline fun
  prior_phi_fn <- local({
    f <- inla.pc.bym.phi(
      eigenvalues = eig_R,
      marginal.variances = marg_var,
      rankdef = 2, # 2nd order RW has nullity 2
      u = prior.phi$u,
      alpha = prior.phi$alpha,
      scale.model = TRUE,
      return.as.table = FALSE,
      adjust.for.con.comp = TRUE
    )
    function(phi) f(phi)
  })

  rmodel <- inla.rgeneric.define(
    model = inla.rw2o1diid.model,
    debug = debug,
    n = n,
    R_scaled = R_scaled,
    eig_R = eig_R,
    prior_phi_fn = prior_phi_fn,
    prior_u_tau = prior.tau$u,
    prior_alpha_tau = prior.tau$alpha
  )

  # Internal-scale names for the hyperparameter output. INLA appends
  # " for <label>" (where <label> is the f() term name) and capitalises
  # the first letter, so these render as e.g. "log precision for time".
  rmodel$f$hyper <- list(
    theta1 = list(name = "log precision"),
    theta2 = list(name = "logit phi")
  )

  rmodel
}

#' @rdname rw2o1diid
#' @export
`inla.rw2o1diid.model` <- function(
  cmd = c("graph", "Q", "mu", "initial", "log.norm.const", "log.prior", "quit"),
  theta = NULL
) {
  interpret.theta <- function() {
    list(
      tau = exp(theta[1]),
      phi = plogis(theta[2])
    )
  }

  graph <- function() {
    R_scaled + Matrix::Diagonal(n)
  }

  Q <- function() {
    param <- interpret.theta()
    param$tau * (param$phi * R_scaled + (1 - param$phi) * Matrix::Diagonal(n))
  }

  mu <- function() {
    numeric(0)
  }

  log.norm.const <- function() {
    param <- interpret.theta()
    eig_Q <- param$tau * (param$phi * eig_R + (1 - param$phi))
    0.5 * sum(log(eig_Q)) - 0.5 * n * log(2 * pi)
  }

  log.prior <- function() {
    param <- interpret.theta()

    # inla.pc.dprec returns log-density wrt tau (external scale) so we need to
    # add the Jacobian theta[1] = log(tau)
    log_prior_tau <- INLA::inla.pc.dprec(
      param$tau,
      u = prior_u_tau,
      alpha = prior_alpha_tau,
      log = TRUE
    ) +
      theta[1]

    # prior_phi_fn returns log-density wrt phi (external scale) so we need to
    # add the Jacobian log(phi) + log(1 - phi)
    log_prior_phi <- prior_phi_fn(param$phi) +
      log(param$phi) +
      log1p(-param$phi)

    log_prior_tau + log_prior_phi
  }

  initial <- function() {
    c(4, 0)
  }

  quit <- function() {
    invisible()
  }

  if (!length(theta)) {
    theta <- initial()
  }

  do.call(match.arg(cmd), args = list())
}

#' Extract user-scale hyperparameter marginals for an `inla.rw2o1diid` term
#'
#' INLA reports `inla.rw2o1diid` hyperparameter marginals on the internal
#' (rgeneric) scale: `theta_1 = log(tau)` and `theta_2 = logit(phi)`. This
#' helper transforms them back to the user scale via [`inla.tmarginal()`],
#' returning marginals for `tau` (the total marginal precision) and `phi`
#' (the structured-variance fraction).
#'
#' @param result An `"inla"` object returned by [`inla()`].
#' @param name Character. The label of the `f()` term that used
#'     `inla.rw2o1diid()`, i.e. the first argument to `f()`. For
#'     `f(time, model = inla.rw2o1diid(n))` this is `"time"`.
#' @return A named list with elements `tau` and `phi`, each an
#'     `inla.marginal` object suitable for [`inla.zmarginal()`],
#'     [`inla.qmarginal()`], etc.
#' @author Antonio R. Vargas
#' @seealso [`inla.rw2o1diid()`], [`inla.tmarginal()`], [`inla.zmarginal()`]
#' @examples
#' \dontrun{
#'   r <- inla(
#'     y ~ f(time, model = inla.rw2o1diid(n)) - 1,
#'     data = data.frame(y, time)
#'   )
#'   hp <- inla.rw2o1diid.hyperpar(r, "time")
#'   inla.zmarginal(hp$tau)
#'   inla.zmarginal(hp$phi)
#' }
#' @export
`inla.rw2o1diid.hyperpar` <- function(result, name) {
  stopifnot(inherits(result, "inla"), is.character(name), length(name) == 1)

  hp_names <- names(result$marginals.hyperpar)
  find_one <- function(prefix) {
    target <- tolower(paste0(prefix, " for ", name))
    hit <- hp_names[tolower(hp_names) == target]
    if (length(hit) != 1) NULL else hit
  }
  key_tau <- find_one("log precision")
  key_phi <- find_one("logit phi")
  if (is.null(key_tau) || is.null(key_phi)) {
    stop(sprintf(
      "Could not find rw2o1diid hyperparameter marginals for term '%s'. Available names: %s",
      name,
      paste(shQuote(hp_names), collapse = ", ")
    ))
  }
  list(
    tau = inla.tmarginal(exp, result$marginals.hyperpar[[key_tau]]),
    phi = inla.tmarginal(plogis, result$marginals.hyperpar[[key_phi]])
  )
}
