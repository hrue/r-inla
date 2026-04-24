## =============================================================
## Rate-Based Mean-Parametrized Gamma-Count Regression
## Models: log(lambda_i) = log(E) + x_i' gamma,  where
## lambda_i = E[N(T_i)] / T_i
## =============================================================

## -------------------------------------------------------------
## Exact renewal-process simulator
## -------------------------------------------------------------
rGC_renewal <- function(n, alpha, beta, T) {
    y <- integer(n)
    for (i in seq_len(n)) {
        t <- 0
        k <- 0L
        repeat {
            t <- t + rgamma(1, shape = alpha, rate = beta[i])
            if (t > T[i]) break
            k <- k + 1L
        }
        y[i] <- k
    }
    y
}

## -------------------------------------------------------------
## Exact Gamma-count PMF (vectorized over T)
## -------------------------------------------------------------
dGC_renewal <- function(y, alpha, beta, T, log = FALSE) {
    Fy1 <- pgamma(T, shape = (y + 1) * alpha, rate = beta)
    Fy0 <- ifelse(y == 0, 1, pgamma(T, shape = y * alpha, rate = beta))
    p <- Fy0 - Fy1
    p <- pmax(p, 1e-14)
    if (log) log(p) else p
}

## -------------------------------------------------------------
## Negative log-likelihood
## -------------------------------------------------------------
negloglik_gc_rate <- function(par, y, X, T, E) {
    if (missing(E))
        E <- rep(1, length(T))
    
    log_alpha <- par[1]
    gamma <- par[-1]

    alpha <- exp(log_alpha)
    if (!is.finite(alpha) || alpha < 0.05) return(1e12)

    ## Model log(rate): lambda = exp(X %*% gamma)
    eta <- X %*% gamma
    lambda <- as.numeric(E * exp(eta))
    if (any(!is.finite(lambda)) || any(lambda < 1e-6)) return(1e12)

    ## Mean count: mu = lambda * T
    mu <- lambda * T

    ## Link: beta = (alpha * mu + 0.5) / T = alpha * lambda + 0.5 / T
    beta <- alpha * lambda + 0.5 / T
    if (any(beta <= 0)) return(1e12)

    ll <- sum(dGC_renewal(y, alpha, beta, T, log = TRUE))
    if (!is.finite(ll)) return(1e12)

    -ll
}

## -------------------------------------------------------------
## Simulate data
## -------------------------------------------------------------
n <- 300
x <- rnorm(n, sd = 0.3)
X <- cbind(1, x)
E <- runif(n)

## Varying observation windows
T_obs <- runif(n, 0.5, 3)  # e.g., different durations per unit

alpha_true <- 2.0
gamma_true <- c(1.0, -0.6)  # log(rate) = gamma0 + gamma1 * x

lambda_true <- as.numeric(E * exp(X %*% gamma_true))  # true rate
mu_true <- lambda_true * T_obs                    # expected count
beta_true <- alpha_true * lambda_true + 0.5 / T_obs

y <- rGC_renewal(n, alpha_true, beta_true, T_obs)

## -------------------------------------------------------------
## ML estimation
## -------------------------------------------------------------
## Initialize with Poisson log-rate
glm_init <- glm(y ~ x, family = poisson(), offset = log(E) + log(T_obs))
gamma_init <- coef(glm_init)
par_init <- c(log(1), gamma_init)

fit <- optim(
    par = par_init,
    fn = negloglik_gc_rate,
    y = y, X = X, T = T_obs, E = E, 
    method = "BFGS", control = list(maxit = 1000, reltol = 1e-10))
alpha_est <- exp(fit$par[1])
gamma_est <- fit$par[-1]

r <- inla(Y ~ 1 + x,
          family = "gammacountmean",
          data = list(Y = inla.mdata(y, T_obs, E), x = x))
print(c(alpha_est,  gamma_est))
summary(r)
