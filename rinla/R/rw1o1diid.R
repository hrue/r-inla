#' 1D First-order random walk plus iid noise (BYM2-style parameterisation)
#'
#' Convenience constructor for a 1D first-order random walk mixed with iid
#' noise, using the BYM2-style `(tau, phi)` parameterisation. The latent
#' field is
#' \deqn{x = \tau^{-1/2} \left(\sqrt{\phi}\, u + \sqrt{1-\phi}\, v\right),}
#' where `u` is a (scaled) RW1 on a linear chain and `v` is iid Gaussian.
#' The `tau` parameter controls the total marginal precision and `phi` in
#' `(0, 1)` controls the fraction of variance attributable to the structured
#' RW1 component.
#'
#' Both hyperparameters use PC priors. The PC prior on `tau` is the standard
#' [`inla.pc.dprec()`] prior, parameterised by `(u, alpha)` such that
#' `P(1/sqrt(tau) > u) = alpha`. The PC prior on `phi` is the adaptive prior
#' computed by `bym2` from the chain-graph Laplacian (see
#' [`inla.pc.bym.phi()`]), parameterised by `(u, alpha)` such that
#' `P(phi < u) = alpha`.
#'
#' @param n Integer. Length of the chain. Must be `>= 3`.
#' @param prior.tau Named list with entries `u` and `alpha` for the PC prior
#'     on the precision `tau`. Default `list(u = 1, alpha = 0.01)`.
#' @param prior.phi Named list with entries `u` and `alpha` for the PC prior
#'     on the mixing parameter `phi`. Default `list(u = 0.5, alpha = 0.5)`.
#' @return A model-specification object that can be passed as the `model`
#'     argument to [`f()`].
#' @author Antonio R. Vargas
#' @seealso [`f()`], `inla.doc("bym2")`
#' @examples
#' \dontrun{
#'   n <- 100
#'   time <- 1:n
#'   y <- cumsum(rnorm(n, sd = 0.3)) + rnorm(n, sd = 1)
#'
#'   ## defaults: PC prior on tau with (u, alpha) = (1, 0.01)
#'   ## and on phi with (u, alpha) = (0.5, 0.5)
#'   r <- inla(
#'     y ~ f(time, model = inla.rw1o1diid(n)) - 1,
#'     data = data.frame(y, time)
#'   )
#'
#'   ## custom priors
#'   r <- inla(
#'     y ~ f(
#'       time,
#'       model = inla.rw1o1diid(
#'         n,
#'         prior.tau = list(u = 2, alpha = 0.05),
#'         prior.phi = list(u = 0.7, alpha = 0.3)
#'       )
#'     ) - 1,
#'     data = data.frame(y, time)
#'   )
#' }
#' @rdname rw1o1diid
#' @export
`inla.rw1o1diid` <- function(
  n,
  prior.tau = list(u = 1, alpha = 0.01),
  prior.phi = list(u = 0.5, alpha = 0.5)
) {
  stopifnot(is.numeric(n), length(n) == 1, n >= 3)
  stopifnot(is.list(prior.tau), all(c("u", "alpha") %in% names(prior.tau)))
  stopifnot(is.list(prior.phi), all(c("u", "alpha") %in% names(prior.phi)))
  n <- as.integer(n)

  # Construct linear chain graph
  i_idx <- c(seq_len(n - 1), seq.int(2, n))
  j_idx <- c(seq.int(2, n), seq_len(n - 1))
  graph <- Matrix::sparseMatrix(i = i_idx, j = j_idx, x = 1, dims = c(n, n))

  hyper <- list(
    prec = list(prior = "pc.prec", param = c(prior.tau$u, prior.tau$alpha)),
    phi = list(prior = "pc", param = c(prior.phi$u, prior.phi$alpha))
  )

  ret <- list(
    f = list(
      model = "bym2",
      graph = graph,
      n = n,
      scale.model = TRUE,
      hyper = hyper
    )
  )
  class(ret) <- "inla.model.class"

  ret
}
