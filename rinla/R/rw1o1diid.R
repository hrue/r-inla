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
#' Internally this is implemented as the existing `bym2` model on the linear
#' chain graph `1 - 2 - 3 - ... - n`, whose graph Laplacian coincides with
#' the RW1 structure matrix. All `bym2` machinery (PC prior on `phi`, hyper
#' overrides, posterior tooling) applies unchanged.
#'
#' @param n Integer. Length of the chain. Must be `>= 3`.
#' @return An object of class `"inla.model.class"` suitable for passing as
#'     the `model` argument to [`f()`].
#' @author Antonio Vargas
#' @seealso inla.doc("bym2")
#' @examples
#' \dontrun{
#'   n <- 100
#'   time <- 1:n
#'   y <- cumsum(rnorm(n, sd = 0.3)) + rnorm(n, sd = 1)
#'   r <- inla(y ~ f(time, model = inla.rw1o1diid(n)) - 1,
#'             data = data.frame(y, time))
#' }
#' @rdname rw1o1diid
#' @export
`inla.rw1o1diid` <- function(n) {
  if (!is.numeric(n) || length(n) != 1) {
    stop("'n' must be a single numeric value.")
  }
  if (n < 3) {
    stop("'n' must be at least 3.")
  }
  n <- as.integer(n)

  i_idx <- c(seq_len(n - 1), seq.int(2, n))
  j_idx <- c(seq.int(2, n), seq_len(n - 1))
  graph <- Matrix::sparseMatrix(i = i_idx, j = j_idx, x = 1, dims = c(n, n))

  ret <- list(
    f = list(
      model = "bym2",
      graph = graph,
      n = n,
      scale.model = TRUE
    )
  )
  class(ret) <- "inla.model.class"
  return(ret)
}
