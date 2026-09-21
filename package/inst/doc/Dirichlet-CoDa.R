## ----setup, include=FALSE-----------------------------------------------------
set.seed(123)
library(INLA)
inla.setOption(num.threads="1:1")
inla.setOption(smtp="taucs")
if (file.exists("myinit.R")) source("myinit.R")
library(knitr)
library(rmarkdown)
knitr::opts_chunk$set(echo=TRUE, cache=FALSE, message=FALSE,warning=FALSE)
knitr::opts_chunk$set(fig.path="figures/conditional-logit/")

## ----warning=FALSE------------------------------------------------------------
  ### --- 0. Loading libraries --- ####
library(INLA)
library(dplyr)
library(ggplot2)

## -------------------------------------------------------------------
set.seed(201803)
inla.seed = sample.int(n=1E6, size=1)
options(width=70, digits=3)

## -------------------------------------------------------------------
### --- 1. Simulation --- ####
# Parameters for the simulation
D <- 3
N <- 1000
sigma2 <- c(0.5, 0.4)
cov_param <- 0.1
sigma_diag <- sqrt(sigma2 + cov_param)
hypers_lik <- data.frame(hypers = c(sigma2, cov_param),
                         name1 = c("sigma2.1", "sigma2.2", "gamma"))
# We create the correlation parameters based on the previous idea
# We are going to have ((D-1)^2 - (D-1))/2 rhos
rho <- diag(1/sigma_diag) %*% matrix(cov_param, D-1, D-1) %*% diag(1/sigma_diag)
diag(rho) <- 1
rho

## -------------------------------------------------------------------
x = runif(N)-0.5
# - mean 0 to not affect intercept
betas = matrix(c(-1, 3, -1, 5), nrow = D-1, byrow = TRUE)
X <- data.frame(1, x) %>% as.matrix(.)
lin.pred <- X %*% t(betas) 

## -------------------------------------------------------------------
### ------- 1.2.3. Constructing the likelihood --- ####
Sigma <- matrix(sigma_diag, ncol = 1) %*% matrix(sigma_diag, nrow = 1)
Sigma <- Sigma*rho

lin.pred %>%
  apply(., 1, function(z)
    MASS::mvrnorm( n  = 1,
             mu = z,
             Sigma = Sigma)) %>%
  t(.)-> alry

## -------------------------------------------------------------------
  y.simplex <- compositions::alrInv(alry)
  y.simplex <- as.numeric(t(y.simplex)) %>% matrix(., ncol = D, byrow = TRUE)
  colnames(y.simplex) <- paste0("y", 1:D)  
  data <- data.frame(alry, y.simplex, x)
colnames(data)[1:(D-1)] <- c(paste0("alry.", 1:(D-1)))
data %>% head(.)

## ----fig.cap = "Simulated data using alr-coordinates in terms of x"----
### Alr coordinates
data %>% 
  tidyr::pivot_longer(., cols = ,starts_with("alr"), 
                      names_to = "y.names", values_to = "y.resp") %>%
  ggplot(data = .) +
  geom_point(aes(x = x, y = y.resp, fill = x), shape = 21, size = 2) +
  ylab("alr") +
  facet_wrap(~y.names) +
  theme_bw() +
  theme(legend.position = "bottom") -> p_alr

#pdf("simulated_data.pdf", width = 8, height = 6)
p_alr
#dev.off()

## -------------------------------------------------------------------
  data$id.z <- 1:dim(data)[1]

## -------------------------------------------------------------------
data_ext <- data %>%
  tidyr::pivot_longer(., cols = all_of(paste0("alry.", 1:(D-1))),
                      names_to  = "y.names",
                      values_to = "y.resp") %>%
  .[order(ordered(.$y.names)),]
data_ext$y.names <- ordered(data_ext$y.names)
head(data_ext)

## -------------------------------------------------------------------
names_y <- paste0("alry.", 1:(D-1))
1:length(names_y) %>%
  lapply(., function(i){
    data_ext %>%
      dplyr::filter(y.names == names_y[i]) -> data_comp_i
    #Response
    y_alr <- matrix(ncol = names_y %>% length(.), nrow = dim(data_comp_i)[1])
    y_alr[, i] <- data_comp_i$y.resp
  }) -> y.resp

1:length(names_y) %>%
  lapply(., function(i){
    y_aux <- data_ext %>%
      dplyr::select(y.resp, y.names) %>%
      dplyr::filter(y.names == names_y[i]) %>%
      dplyr::select(y.resp) %>%
      as.matrix(.)
    aux_vec <- rep(NA, (D-1))
    aux_vec[i] <- 1
    kronecker(aux_vec, y_aux)
  }) -> y_list

y_tot <- do.call(cbind, y_list)
y_tot %>% head(.)

## -------------------------------------------------------------------
variables <- c("intercept", data %>%
                 dplyr::select(starts_with("x")) %>%
                 colnames(.))
id.names <- paste0("id.", variables)
id.variables <- rep(data_ext$y.names %>% as.factor(.) %>% as.numeric(.), 
                    length(variables)) %>%
  matrix(., ncol = length(variables), byrow = FALSE)
colnames(id.variables) <- id.names

variables
id.variables %>% head(.)

## -------------------------------------------------------------------
stk.est <- inla.stack(data    = list(resp = y_tot),
                      A       = list(1),
                      effects = list(cbind(data_ext %>%
                                             dplyr::select(starts_with("x")),
                                           data_ext %>%
                                             dplyr::select(starts_with("id.z")),
                                           id.variables,
                                           intercept = 1)),
                      tag     = 'est')

## -------------------------------------------------------------------
  # Have different parameters for fixed effects, and do not include spatial random effects.
list_prior <- rep(list(list(prior = "pc.prec", param = c(1, 0.01))), D-1)

### Fitting the model
formula.typeII <- resp ~ -1 +
  f(id.intercept, intercept,
    model   = "iid",
    initial = log(1/1000),
    fixed   = TRUE) +
  f(id.x, x,
    model   = "iid",
    initial = log(1/1000),
    fixed   = TRUE) +
  f(id.z,
    model = "iid",
    hyper = list(prec = list(prior = "pc.prec",
                             param = c(1, 0.01))), constr = TRUE)
model.typeII <- inla(formula.typeII,
                     family         = rep("gaussian", D-1),
                     data           = inla.stack.data(stk.est),
                     control.compute = list(config = TRUE),
                     control.predictor = list(A = inla.stack.A(stk.est),
                                              compute = TRUE),
                     control.family = list_prior,
                     inla.mode = "experimental" ,
                     verbose = FALSE)

## -------------------------------------------------------------------
### Posterior distribution of the fixed effects
data_fixed <- rbind(data.frame(inla.smarginal(model.typeII$marginals.random$id.x$index.1),
                               alr = "alr(y1/y3)",
                               var = "x",
                               param = "beta1",
                               real  = betas[1,2]),
                    data.frame(inla.smarginal(model.typeII$marginals.random$id.x$index.2),
                               alr = "alr(y2/y3)",
                               var = "x",
                               param = "beta1",
                               real  = betas[2,2]),
                    data.frame(inla.smarginal(model.typeII$marginals.random$id.intercept$index.1),
                               alr = "alr(y1/y3)",
                               var = "intercept",
                               param = "beta0",
                               real = betas[1,1]),
                    data.frame(inla.smarginal(model.typeII$marginals.random$id.intercept$index.2),
                               alr = "alr(y2/y3)",
                               var = "intercept",
                               param = "beta0",
                               real = betas[2,1]))

p_fixed <- ggplot() +
  geom_line(data = data_fixed, aes(x = x, y = y), size = 0.9) +
  #ggtitle("Effect of the covariate bio12") +
  theme_bw() +
  geom_vline(data = data_fixed, aes(xintercept = real), col = "red4") +
 # scale_color_manual(values=c("#E75F00", "#56B4E9"))+
  theme(legend.position = "bottom") +
  facet_wrap(~param + alr, ncol = D-1, scales = "free") +
  xlab(expression(beta^(d))) +
  ylab(expression(p(beta^(d) *'|'* theta))) +
  theme(legend.title = element_blank())

#pdf("posterior_fixed.pdf", width = 6, height = 5)
p_fixed
#dev.off()

## ----fig.asp = 0.5--------------------------------------------------
### Posterior distribution of the hyperparameters
prec <- list(sigma2.1 = model.typeII$marginals.hyperpar$`Precision for the Gaussian observations`,
             sigma2.2 = model.typeII$marginals.hyperpar$`Precision for the Gaussian observations[2]`,
             gamma = model.typeII$marginals.hyper$`Precision for id.z`)

hyper <- lapply(1:length(prec),
                function(x){
                  inla.smarginal(inla.tmarginal(prec[[x]], fun = function(y)(1/y))) %>%
                    data.frame(.)
                })
names(hyper) <- names(prec)

hyper.df <- lapply(1:length(hyper),
                   function(x){
                     cbind(data.frame(hyper[[x]]), name1 = names(hyper)[x])
                   })  %>%
  do.call(rbind.data.frame, .)

hyper.df$name1 <- ordered(hyper.df$name1,
                          levels = c("sigma2.1", "sigma2.2",
                                     "gamma"))
p.hyper <- ggplot(hyper.df) +
  geom_line(aes(x = x, y = y)) +
  geom_vline(data = hypers_lik, aes(xintercept = hypers), col = "red4") +
  facet_wrap(~ name1, scales = "free") +
  theme_bw() +
  xlab(expression(theta)) +
  ylab(expression(p(theta*'|'*y)))

#pdf("marginals_hyperpar.pdf", width = 6, height = 3)
print(p.hyper)
#dev.off()

## -------------------------------------------------------------------
sim <- 1000
x.pred <- seq(-0.5, 0.5, 0.3)
n.pred <- length(x.pred)
cat("\n ----------------------------------------------- \n")
cat("Creating the data.frame for predictions \n")

data_pred <- data.frame(intercept = 1,
                        x = rep(x.pred, D-1))
id.z.pred  <- rep((N + 1):(N + n.pred), D - 1) #random effect z to model the correlation

# Category
id.cat_pred <- rep(1:(D - 1), rep(n.pred, D - 1))
#Index for covariates
variables_pred <- c("intercept", data_pred %>% 
                      dplyr::select(starts_with("x")) %>% 
                      colnames(.))
id.names_pred <- paste0("id.", variables_pred)
id.variables_pred <- rep(id.cat_pred, length(variables_pred)) %>% 
  matrix(., ncol = length(variables_pred), byrow = FALSE)
colnames(id.variables_pred) <- id.names_pred

## -------------------------------------------------------------------
stk.pred <- inla.stack(data    = list(resp = matrix(NA, ncol = D - 1, 
                                                    nrow = n.pred*(D - 1))),
                       A       = list(1),
                       effects = list(cbind(data_pred,
                                            id.z = id.z.pred,
                                            id.variables_pred)),
                       tag     = 'pred')
### --- Total stack
stk <- inla.stack(stk.est, stk.pred)

## -------------------------------------------------------------------
mod.pred <- inla(formula.typeII, 
                 family         = rep("gaussian", D - 1),
                 data              = inla.stack.data(stk), 
                 control.compute   = list(config = TRUE),
                 control.predictor = list(A = inla.stack.A(stk), compute = TRUE, link = 1),
                 control.mode      = list(theta = model.typeII$mode$theta, restart = TRUE), 
                 control.family = list_prior,   
                 num.threads       = 2,
                 inla.mode = "experimental" , 
                 verbose           = FALSE)

## -------------------------------------------------------------------
pred.values.mean <- mod.pred$summary.fitted.values$mean[inla.stack.index(stk, 'pred')$data] %>% 
  matrix(., ncol = D - 1, byrow = FALSE)

post_sim_pred <- inla.posterior.sample(n = sim, result = mod.pred)
post_sim_predictor <- inla.posterior.sample.eval(fun = function(...){
  APredictor}, post_sim_pred, return.matrix = TRUE)
post_sim_idz <- inla.posterior.sample.eval(fun = function(...){
  id.z}, post_sim_pred, return.matrix = TRUE)

ind.pred <- inla.stack.index(stk, 'pred')$data
ind.idz <- inla.stack.index(stk, 'est')$data #This is the shared random effect
ind.idz <- ind.idz[1:(length(ind.idz)/(D - 1))]

post_sim_predictor[ind.pred, ] <- post_sim_predictor[ind.pred, ]- 
  kronecker(rep(1, D-1), post_sim_idz[-ind.idz,])

post_sim_pred_alr <- post_sim_predictor[ind.pred,]

#Computing mean and sd
pred_alr_summary <- t(apply(post_sim_pred_alr, 1, function(x){c(mean(x), sd(x))}))
pred_alr_summary <- data.frame(pred_alr_summary, 
                               y.names = rep(names_y, rep(n.pred, D-1)),
                               x.pred = rep(x.pred, D-1))
colnames(pred_alr_summary)[1:2] <- c("mean", "sd")

pred_alr_summary

## -------------------------------------------------------------------
###  Prediction in the simplex --- #####
apply(post_sim_predictor[ind.pred,], 2, function(x){
  alr_pred <- matrix(x, ncol = D - 1)
  pred_simplex <- compositions::alrInv(alr_pred)
  as.numeric(t(pred_simplex)) #Byrows
}) -> post_sim_pred_simplex

#Computing credible intervals
pred_simplex_summary <- t(apply(post_sim_pred_simplex, 1, function(x){c(mean(x), sd(x))}))
pred_simplex_summary <- data.frame(pred_simplex_summary, 
                                   y.names = rep(c("y1", "y2", "y3"), n.pred),
                                   x.pred  = rep(x.pred, rep(D, n.pred)))
colnames(pred_simplex_summary)[1:2] <- c("mean", "sd")

pred_simplex_summary

