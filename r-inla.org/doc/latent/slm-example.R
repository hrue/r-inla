## Example using the Boston dataset from package spdep

require(INLA)
require(spdep)
data(boston)

## Index for the latent model
n <- nrow(boston.c)
boston.c$idx <- 1:n

## Define adjacency using a row-standardised matrix
lw <- nb2listw(boston.soi)
W <- as(as_dgRMatrix_listw(lw), "CsparseMatrix")

## Model definition
f1 <- log(CMEDV) ~ CRIM + ZN + INDUS + CHAS + I(NOX^2)+ I(RM^2) +  AGE + 
  log(DIS) + log(RAD) + TAX + PTRATIO + B + log(LSTAT)
mmatrix <- model.matrix(f1, boston.c)

## Zero-variance for error term
zero.variance = list(prec=list(initial = 25, fixed=TRUE))

## Compute eigenvalues for SLM model, used to obtain rho.min and
## rho.max
e = eigenw(lw)
re.idx = which(abs(Im(e)) < 1e-6)
rho.max = 1/max(Re(e[re.idx]))
rho.min = 1/min(Re(e[re.idx]))
rho = mean(c(rho.min, rho.max))

## Precision matrix for beta coeffients' prior
betaprec <- .0001
Q.beta = Diagonal(n=ncol(mmatrix), betaprec)


## Priors on the hyperparameters
hyper = list(
        prec = list(
                prior = "loggamma",
                param = c(0.01, 0.01)), 
        rho = list(
                initial=0,
                prior = "logitbeta",
                param = c(1,1)))

## Fit model
slmm1 <- inla( log(CMEDV) ~ -1 +
              f(idx, model="slm",
                args.slm=list(
                        rho.min = rho.min,
                        rho.max = rho.max,
                        W=W,
                        X=mmatrix,
                        Q.beta=Q.beta),
                hyper=hyper),
              data=boston.c, family="gaussian",
              control.family = list(hyper=zero.variance),
              control.compute=list(dic=TRUE, cpo=TRUE)
              )
summary(slmm1)

## Summary of the coefficients (at the end of the vector of random effects)
slmm1$summary.random$idx[n+1:ncol(mmatrix),]

## Re-scale rho to real scale
rhomarg <- inla.tmarginal(function(x){rho.min+x*(rho.max-rho.min)},
                          slmm1$marginals.hyperpar[[2]])
inla.zmarginal(rhomarg)

## Maximum likelihood estimate of model (used for comparison)
summary(m2 <- lagsarlm(f1, boston.c, lw))
