#' 1D Second-order random walk plus iid noise (BYM2-style parameterisation)
#'
#' Convenience constructor for a 1D second-order random walk mixed with iid
#' noise, using a BYM2-style `(tau, phi)` parameterisation. The latent
#' field has precision
#' \deqn{Q(\tau, \phi) = \tau \, \bigl( \phi \, R + (1-\phi)\, I \bigr),}
#' where `R` is the RW2 structure matrix scaled so that the geometric mean
#' of the diagonal of its (Moore-Penrose) pseudoinverse is one, and `I` is
#' the identity. The `tau` parameter controls the overall precision of the
#' field and `phi` in `(0, 1)` controls the relative weight of the
#' structured RW2 component vs. the iid component.
#'
#' Note that this is a *precision-mixture* parameterisation, which is not
#' exactly equivalent to BYM2's *variance-mixture* parameterisation. BYM2
#' constructs the field as
#' \deqn{b = \tau^{-1/2}\bigl(\sqrt{\phi}\, u + \sqrt{1-\phi}\, v\bigr),}
#' so that the marginal variance of `b` is exactly `(1/tau)` and `phi` is
#' exactly the fraction of marginal variance contributed by the structured
#' component `u`. The two parameterisations agree at `phi = 0` and (under
#' the scaling above) approximately at `phi = 1`, but differ in the
#' interior: under our `Q`, the marginal variance of the field is only
#' approximately `1/tau`, and `phi` is only approximately the
#' structured-variance fraction. The precision-mixture form is used here
#' because it admits a single n-dimensional latent representation suitable
#' for `rgeneric`, whereas the variance-mixture form requires a 2n latent
#' space (which is what `bym2` uses internally).
#'
#' Both hyperparameters use PC priors. The PC prior on `tau` is the standard
#' [`inla.pc.dprec()`] prior, parameterised by `(u, alpha)` such that
#' `P(1/sqrt(tau) > u) = alpha`. The PC prior on `phi` is computed from the
#' KLD between the precision-mixture flexible model and the iid base model
#' (see [`inla.pc.rw2o1diid.phi()`]), parameterised by `(u, alpha)` such
#' that `P(phi < u) = alpha`.
#'
#' Because this model is implemented as an rgeneric, INLA reports the
#' hyperparameter marginals on the internal scale under generic names:
#' `"Theta1 for <label>"` is `log(tau)` and `"Theta2 for <label>"` is
#' `logit(phi)`, where `<label>` is the variable name passed as the first
#' argument to f() (e.g. "time" for f(time, ...)). To obtain marginals
#' for `tau` and `phi` on the user scale, use [`inla.rw2o1diid.hyperpar()`].
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
#'   y <- cumsum(rnorm(n, sd = 0.3)) + rnorm(n, sd = 1)
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

  R_scaled <- inla.rw(n, order = 2, scale.model = TRUE, sparse = TRUE)
  eig_R <- pmax(
    0,
    eigen(as.matrix(R_scaled), symmetric = TRUE, only.values = TRUE)$values
  )

  # local() trims the captured environment to just the spline fun
  prior_phi_fn <- local({
    f <- inla.pc.rw2o1diid.phi(
      eigenvalues = eig_R,
      u = prior.phi$u,
      alpha = prior.phi$alpha
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

  # Sum-to-zero constraint on the latent field. This matches what bym2 and
  # the built-in rw2 do: the constant null mode of R is constrained out
  # (otherwise tau's interpretation as overall precision drifts as phi -> 1
  # because the constant null direction inflates the marginal variance).
  rmodel$f$extraconstr <- list(
    A = matrix(1, nrow = 1, ncol = n),
    e = 0
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
#' (rgeneric) scale under generic names: `"Theta1 for <label>"` is
#' `log(tau)` and `"Theta2 for <label>"` is `logit(phi)`. This helper
#' transforms them back to the user scale via [`inla.tmarginal()`],
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
  key_tau <- find_one("Theta1")
  key_phi <- find_one("Theta2")
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

#' PC prior on phi for the rw2o1diid (precision-mixture) model
#'
#' Builds the penalised-complexity prior on the mixing parameter `phi` of an
#' `inla.rw2o1diid` term. The flexible model has covariance
#' \deqn{\Sigma(\phi) = \tau^{-1}\bigl(\phi R + (1-\phi) I\bigr)^{-1}}
#' and the base model is `phi = 0` (pure iid with covariance `tau^{-1} I`).
#' With `tau` held fixed, the Gaussian KLD between flexible and base reduces
#' to
#' \deqn{2\,\mathrm{KLD}(\phi) = \sum_i \frac{1}{\phi\lambda_i + (1-\phi)} - n + \sum_i \log\bigl(\phi\lambda_i + (1-\phi)\bigr),}
#' where the `lambda_i` are the eigenvalues of `R` (including any zeros from
#' the null space). The PC distance is `dist(phi) = sqrt(2*KLD(phi))`, and
#' the prior on distance is `Exp(lambda)` with `lambda` chosen so that
#' `P(phi < u) = alpha`.
#'
#' This is the precision-mixture analogue of [`inla.pc.bym.phi()`], which
#' derives the PC prior from BYM2's *variance-mixture* KLD. The two priors
#' agree at `phi = 0` and approximately at `phi = 1`, but differ in the
#' interior.
#'
#' @param eigenvalues Numeric vector of eigenvalues of the (scaled) RW2
#'     structure matrix `R`, including any zeros from the null space.
#' @param u,alpha Calibration parameters for the PC prior on `phi`, such
#'     that `P(phi < u) = alpha`. Both must lie strictly in `(0, 1)`.
#' @param return.as.table Currently unused; retained for parity with
#'     [`inla.pc.bym.phi()`]. The function always returns a closure.
#' @return A function of one argument `phi` returning the log-density of the
#'     PC prior on `phi`, evaluated on `(0, 1)`.
#' @seealso [`inla.rw2o1diid()`], [`inla.pc.bym.phi()`]
`inla.pc.rw2o1diid.phi` <- function(
  eigenvalues,
  u,
  alpha,
  return.as.table = FALSE
) {
  stopifnot(alpha > 0, alpha < 1, u > 0, u < 1)

  eigenvalues_minus_1 <- pmax(0, eigenvalues) - 1

  # Algebraically equivalent to sum(1/w) - n + sum(log(w)) with
  # w = phi*eig + (1-phi), but avoids the cancellation between sum(1/w) and n
  # at small phi (and the log1p form is more accurate there too).
  twice_kld <- function(phi) {
    v <- phi * eigenvalues_minus_1
    sum(-v / (1 + v)) + sum(log1p(v))
  }

  # Evaluate distance dist(phi) = sqrt(2*KLD) on a coarse logit-phi grid, then
  # build a spline of log(dist) vs logit(phi) so we can resample/differentiate.
  logit_phi_grid <- seq(-25, 25, length.out = 1000)
  phi_grid <- plogis(logit_phi_grid)
  dist_grid <- vapply(
    phi_grid,
    function(phi) {
      val <- twice_kld(phi)
      if (val >= 0) sqrt(val) else NA_real_
    },
    numeric(1)
  )

  good_idx <- !is.na(dist_grid)
  dist_grid <- dist_grid[good_idx]
  logit_phi_grid <- logit_phi_grid[good_idx]

  log_dist_spline <- splinefun(logit_phi_grid, log(dist_grid))
  dist_at_logit_phi <- function(logit_phi, deriv = 0) {
    if (deriv == 0) {
      exp(log_dist_spline(logit_phi))
    } else if (deriv == 1) {
      exp(log_dist_spline(logit_phi)) *
        log_dist_spline(logit_phi, deriv = 1)
    } else {
      stop("deriv must be 0 or 1")
    }
  }
  dist_at_phi <- function(phi) dist_at_logit_phi(qlogis(phi))

  logit_phi_dense <- seq(
    min(logit_phi_grid),
    max(logit_phi_grid),
    length.out = 10000
  )
  phi_dense <- plogis(logit_phi_dense)
  dist_dense <- dist_at_phi(phi_dense)

  # P(phi < u) = P(dist < dist(u)) = 1 - exp(-lambda * dist(u)) = alpha,
  # so lambda = -log(1 - alpha) / dist(u).
  lambda <- -log(1 - alpha) / dist_at_phi(u)

  # log p(logit_phi) = log p(dist) + log|d(dist)/d(logit_phi)|
  log_jacobian_logit <- log(abs(dist_at_logit_phi(logit_phi_dense, deriv = 1)))
  log_prior_logit <- log(lambda) - lambda * dist_dense + log_jacobian_logit

  # Renormalise so that the prior integrates to 1 over logit(phi).
  trapezoid_weights <- (c(0, diff(logit_phi_dense)) +
    c(diff(logit_phi_dense), 0)) /
    2
  norm_const <- sum(exp(log_prior_logit) * trapezoid_weights)
  log_prior_logit <- log_prior_logit - log(norm_const)

  # https://cran.r-project.org/web/packages/Rmpfr/vignettes/log1mexp-note.pdf
  log1pexp <- function(x) {
    x_bins <- .bincode(x, c(-Inf, -37, 18, 33.3, Inf))

    x_in_bin <- which(x_bins == 1)
    if (length(x_in_bin) > 0) {
      x[x_in_bin] <- exp(x[x_in_bin])
    }

    x_in_bin <- which(x_bins == 2)
    if (length(x_in_bin) > 0) {
      x[x_in_bin] <- log1p(exp(x[x_in_bin]))
    }

    x_in_bin <- which(x_bins == 3)
    if (length(x_in_bin) > 0) {
      x[x_in_bin] <- x[x_in_bin] + exp(-x[x_in_bin])
    }

    x
  }

  # Convert log prior to phi scale via Jacobian:
  #   log p(phi) = log p(logit_phi) - log(phi*(1-phi))
  #              = log p(logit_phi) + logit_phi + 2*log(1 + exp(-logit_phi))
  log_prior_phi_spline <- splinefun(
    logit_phi_dense,
    log_prior_logit + logit_phi_dense + 2 * log1pexp(-logit_phi_dense)
  )
  log_prior_phi <- function(phi) log_prior_phi_spline(qlogis(phi))

  log_prior_phi
}
