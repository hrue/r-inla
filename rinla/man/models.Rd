
%% DO NOT EDIT!
%% This file is generated automatically from models.R
\name{inla.models}
\alias{inla.models}
\title{Valid models in INLA}
\description{
This page describe the models implemented in \code{inla}, divided into sections: latent, group, mix, link, predictor, hazard, likelihood, prior, wrapper .
}

\usage{
inla.models()
}
\value{
 	 Valid sections are: latent, group, mix, link, predictor, hazard, likelihood, prior, wrapper 

 	 \describe{
		\item{Section `latent'.}{
			 Valid models in this section are:
			 \describe{
				\item{Model `linear'.}{Number of hyperparmeters are 0.
				 }
				\item{Model `iid'.}{Number of hyperparmeters are 1.
					 \describe{
						\item{Hyperparameter `theta'}{
							 \describe{
							\item{hyperid = }{`1001'}
							\item{name = }{`log precision'}
							\item{short.name = }{`prec'}
							\item{prior = }{`loggamma'}
							\item{param = }{`1 5e-05'}
							\item{initial = }{`4'}
							\item{fixed = }{`FALSE'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{doc = }{`Gaussian random effects in dim=1'}
								\item{constr = }{`FALSE'}
								\item{nrow.ncol = }{`FALSE'}
								\item{augmented = }{`FALSE'}
								\item{aug.factor = }{`1'}
								\item{aug.constr = }{`NULL'}
								\item{n.div.by = }{`NULL'}
								\item{n.required = }{`FALSE'}
								\item{set.default.values = }{`FALSE'}
								\item{pdf = }{`indep'}
							 }
						 }
					 }
				 }
				\item{Model `mec'.}{Number of hyperparmeters are 4.
					 \describe{
						\item{Hyperparameter `theta1'}{
							 \describe{
							\item{hyperid = }{`2001'}
							\item{name = }{`beta'}
							\item{short.name = }{`b'}
							\item{prior = }{`gaussian'}
							\item{param = }{`1 0.001'}
							\item{initial = }{`1'}
							\item{fixed = }{`FALSE'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta2'}{
							 \describe{
							\item{hyperid = }{`2002'}
							\item{name = }{`prec.u'}
							\item{short.name = }{`prec'}
							\item{prior = }{`loggamma'}
							\item{param = }{`1 1e-04'}
							\item{initial = }{`9.21034037197618'}
							\item{fixed = }{`TRUE'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
						\item{Hyperparameter `theta3'}{
							 \describe{
							\item{hyperid = }{`2003'}
							\item{name = }{`mean.x'}
							\item{short.name = }{`mu.x'}
							\item{prior = }{`gaussian'}
							\item{param = }{`0 1e-04'}
							\item{initial = }{`0'}
							\item{fixed = }{`TRUE'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta4'}{
							 \describe{
							\item{hyperid = }{`2004'}
							\item{name = }{`prec.x'}
							\item{short.name = }{`prec.x'}
							\item{prior = }{`loggamma'}
							\item{param = }{`1 10000'}
							\item{initial = }{`-9.21034037197618'}
							\item{fixed = }{`TRUE'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{doc = }{`Classical measurement error model'}
								\item{constr = }{`FALSE'}
								\item{nrow.ncol = }{`FALSE'}
								\item{augmented = }{`FALSE'}
								\item{aug.factor = }{`1'}
								\item{aug.constr = }{`NULL'}
								\item{n.div.by = }{`NULL'}
								\item{n.required = }{`FALSE'}
								\item{set.default.values = }{`FALSE'}
								\item{pdf = }{`mec'}
							 }
						 }
					 }
				 }
				\item{Model `meb'.}{Number of hyperparmeters are 2.
					 \describe{
						\item{Hyperparameter `theta1'}{
							 \describe{
							\item{hyperid = }{`3001'}
							\item{name = }{`beta'}
							\item{short.name = }{`b'}
							\item{prior = }{`gaussian'}
							\item{param = }{`1 0.001'}
							\item{initial = }{`1'}
							\item{fixed = }{`FALSE'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta2'}{
							 \describe{
							\item{hyperid = }{`3002'}
							\item{name = }{`prec.u'}
							\item{short.name = }{`prec'}
							\item{prior = }{`loggamma'}
							\item{param = }{`1 1e-04'}
							\item{initial = }{`6.90775527898214'}
							\item{fixed = }{`FALSE'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{doc = }{`Berkson measurement error model'}
								\item{constr = }{`FALSE'}
								\item{nrow.ncol = }{`FALSE'}
								\item{augmented = }{`FALSE'}
								\item{aug.factor = }{`1'}
								\item{aug.constr = }{`NULL'}
								\item{n.div.by = }{`NULL'}
								\item{n.required = }{`FALSE'}
								\item{set.default.values = }{`FALSE'}
								\item{pdf = }{`meb'}
							 }
						 }
					 }
				 }
				\item{Model `rgeneric'.}{Number of hyperparmeters are 0.
				 }
				\item{Model `rw1'.}{Number of hyperparmeters are 1.
					 \describe{
						\item{Hyperparameter `theta'}{
							 \describe{
							\item{hyperid = }{`4001'}
							\item{name = }{`log precision'}
							\item{short.name = }{`prec'}
							\item{prior = }{`loggamma'}
							\item{param = }{`1 5e-05'}
							\item{initial = }{`4'}
							\item{fixed = }{`FALSE'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{doc = }{`Random walk of order 1'}
								\item{constr = }{`TRUE'}
								\item{nrow.ncol = }{`FALSE'}
								\item{augmented = }{`FALSE'}
								\item{aug.factor = }{`1'}
								\item{aug.constr = }{`NULL'}
								\item{n.div.by = }{`NULL'}
								\item{n.required = }{`FALSE'}
								\item{set.default.values = }{`FALSE'}
								\item{min.diff = }{`1e-05'}
								\item{pdf = }{`rw1'}
							 }
						 }
					 }
				 }
				\item{Model `rw2'.}{Number of hyperparmeters are 1.
					 \describe{
						\item{Hyperparameter `theta'}{
							 \describe{
							\item{hyperid = }{`5001'}
							\item{name = }{`log precision'}
							\item{short.name = }{`prec'}
							\item{prior = }{`loggamma'}
							\item{param = }{`1 5e-05'}
							\item{initial = }{`4'}
							\item{fixed = }{`FALSE'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{doc = }{`Random walk of order 2'}
								\item{constr = }{`TRUE'}
								\item{nrow.ncol = }{`FALSE'}
								\item{augmented = }{`FALSE'}
								\item{aug.factor = }{`1'}
								\item{aug.constr = }{`NULL'}
								\item{n.div.by = }{`NULL'}
								\item{n.required = }{`FALSE'}
								\item{set.default.values = }{`FALSE'}
								\item{min.diff = }{`0.001'}
								\item{pdf = }{`rw2'}
							 }
						 }
					 }
				 }
				\item{Model `crw2'.}{Number of hyperparmeters are 1.
					 \describe{
						\item{Hyperparameter `theta'}{
							 \describe{
							\item{hyperid = }{`6001'}
							\item{name = }{`log precision'}
							\item{short.name = }{`prec'}
							\item{prior = }{`loggamma'}
							\item{param = }{`1 5e-05'}
							\item{initial = }{`4'}
							\item{fixed = }{`FALSE'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{doc = }{`Exact solution to the random walk of order 2'}
								\item{constr = }{`TRUE'}
								\item{nrow.ncol = }{`FALSE'}
								\item{augmented = }{`FALSE'}
								\item{aug.factor = }{`2'}
								\item{aug.constr = }{`1'}
								\item{n.div.by = }{`NULL'}
								\item{n.required = }{`FALSE'}
								\item{set.default.values = }{`FALSE'}
								\item{min.diff = }{`0.001'}
								\item{pdf = }{`crw2'}
							 }
						 }
					 }
				 }
				\item{Model `seasonal'.}{Number of hyperparmeters are 1.
					 \describe{
						\item{Hyperparameter `theta'}{
							 \describe{
							\item{hyperid = }{`7001'}
							\item{name = }{`log precision'}
							\item{short.name = }{`prec'}
							\item{prior = }{`loggamma'}
							\item{param = }{`1 5e-05'}
							\item{initial = }{`4'}
							\item{fixed = }{`FALSE'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{doc = }{`Seasonal model for time series'}
								\item{constr = }{`FALSE'}
								\item{nrow.ncol = }{`FALSE'}
								\item{augmented = }{`FALSE'}
								\item{aug.factor = }{`1'}
								\item{aug.constr = }{`NULL'}
								\item{n.div.by = }{`NULL'}
								\item{n.required = }{`FALSE'}
								\item{set.default.values = }{`FALSE'}
								\item{pdf = }{`seasonal'}
							 }
						 }
					 }
				 }
				\item{Model `besag'.}{Number of hyperparmeters are 1.
					 \describe{
						\item{Hyperparameter `theta'}{
							 \describe{
							\item{hyperid = }{`8001'}
							\item{name = }{`log precision'}
							\item{short.name = }{`prec'}
							\item{prior = }{`loggamma'}
							\item{param = }{`1 5e-05'}
							\item{initial = }{`4'}
							\item{fixed = }{`FALSE'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{doc = }{`The Besag area model (CAR-model)'}
								\item{constr = }{`TRUE'}
								\item{nrow.ncol = }{`FALSE'}
								\item{augmented = }{`FALSE'}
								\item{aug.factor = }{`1'}
								\item{aug.constr = }{`NULL'}
								\item{n.div.by = }{`NULL'}
								\item{n.required = }{`TRUE'}
								\item{set.default.values = }{`TRUE'}
								\item{pdf = }{`besag'}
							 }
						 }
					 }
				 }
				\item{Model `besag2'.}{Number of hyperparmeters are 2.
					 \describe{
						\item{Hyperparameter `theta1'}{
							 \describe{
							\item{hyperid = }{`9001'}
							\item{name = }{`log precision'}
							\item{short.name = }{`prec'}
							\item{prior = }{`loggamma'}
							\item{param = }{`1 5e-05'}
							\item{initial = }{`4'}
							\item{fixed = }{`FALSE'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
						\item{Hyperparameter `theta2'}{
							 \describe{
							\item{hyperid = }{`9002'}
							\item{name = }{`scaling parameter'}
							\item{short.name = }{`a'}
							\item{prior = }{`loggamma'}
							\item{param = }{`10 10'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{doc = }{`The shared Besag model'}
								\item{constr = }{`FALSE'}
								\item{nrow.ncol = }{`FALSE'}
								\item{augmented = }{`FALSE'}
								\item{aug.factor = }{`1'}
								\item{aug.constr = }{`1 2'}
								\item{n.div.by = }{`2'}
								\item{n.required = }{`TRUE'}
								\item{set.default.values = }{`TRUE'}
								\item{pdf = }{`besag2'}
							 }
						 }
					 }
				 }
				\item{Model `bym'.}{Number of hyperparmeters are 2.
					 \describe{
						\item{Hyperparameter `theta1'}{
							 \describe{
							\item{hyperid = }{`10001'}
							\item{name = }{`log unstructured precision'}
							\item{short.name = }{`prec.unstruct'}
							\item{prior = }{`loggamma'}
							\item{param = }{`1 5e-04'}
							\item{initial = }{`4'}
							\item{fixed = }{`FALSE'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
						\item{Hyperparameter `theta2'}{
							 \describe{
							\item{hyperid = }{`10002'}
							\item{name = }{`log spatial precision'}
							\item{short.name = }{`prec.spatial'}
							\item{prior = }{`loggamma'}
							\item{param = }{`1 5e-04'}
							\item{initial = }{`4'}
							\item{fixed = }{`FALSE'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{doc = }{`The BYM-model (Besag-York-Mollier model)'}
								\item{constr = }{`TRUE'}
								\item{nrow.ncol = }{`FALSE'}
								\item{augmented = }{`TRUE'}
								\item{aug.factor = }{`2'}
								\item{aug.constr = }{`2'}
								\item{n.div.by = }{`NULL'}
								\item{n.required = }{`TRUE'}
								\item{set.default.values = }{`TRUE'}
								\item{pdf = }{`bym'}
							 }
						 }
					 }
				 }
				\item{Model `bym2'.}{Number of hyperparmeters are 2.
					 \describe{
						\item{Hyperparameter `theta1'}{
							 \describe{
							\item{hyperid = }{`11001'}
							\item{name = }{`log precision'}
							\item{short.name = }{`prec'}
							\item{prior = }{`pc.prec'}
							\item{param = }{`1 0.01'}
							\item{initial = }{`4'}
							\item{fixed = }{`FALSE'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
						\item{Hyperparameter `theta2'}{
							 \describe{
							\item{hyperid = }{`11002'}
							\item{name = }{`logit phi'}
							\item{short.name = }{`phi'}
							\item{prior = }{`pc'}
							\item{param = }{`0.5 0.5'}
							\item{initial = }{`-3'}
							\item{fixed = }{`FALSE'}
							\item{to.theta = }{`\code{function(x) log(x/(1-x))}'}
							\item{from.theta = }{`\code{function(x) exp(x)/(1+exp(x))}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{doc = }{`The BYM-model with the PC priors'}
								\item{constr = }{`TRUE'}
								\item{nrow.ncol = }{`FALSE'}
								\item{augmented = }{`TRUE'}
								\item{aug.factor = }{`2'}
								\item{aug.constr = }{`2'}
								\item{n.div.by = }{`NULL'}
								\item{n.required = }{`TRUE'}
								\item{set.default.values = }{`TRUE'}
								\item{status = }{`experimental'}
								\item{pdf = }{`bym2'}
							 }
						 }
					 }
				 }
				\item{Model `besagproper'.}{Number of hyperparmeters are 2.
					 \describe{
						\item{Hyperparameter `theta1'}{
							 \describe{
							\item{hyperid = }{`12001'}
							\item{name = }{`log precision'}
							\item{short.name = }{`prec'}
							\item{prior = }{`loggamma'}
							\item{param = }{`1 5e-04'}
							\item{initial = }{`2'}
							\item{fixed = }{`FALSE'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
						\item{Hyperparameter `theta2'}{
							 \describe{
							\item{hyperid = }{`12002'}
							\item{name = }{`log diagonal'}
							\item{short.name = }{`diag'}
							\item{prior = }{`loggamma'}
							\item{param = }{`1 1'}
							\item{initial = }{`1'}
							\item{fixed = }{`FALSE'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{doc = }{`A proper version of the Besag model'}
								\item{constr = }{`FALSE'}
								\item{nrow.ncol = }{`FALSE'}
								\item{augmented = }{`FALSE'}
								\item{aug.factor = }{`1'}
								\item{aug.constr = }{`NULL'}
								\item{n.div.by = }{`NULL'}
								\item{n.required = }{`TRUE'}
								\item{set.default.values = }{`TRUE'}
								\item{status = }{`experimental'}
								\item{pdf = }{`besagproper'}
							 }
						 }
					 }
				 }
				\item{Model `besagproper2'.}{Number of hyperparmeters are 2.
					 \describe{
						\item{Hyperparameter `theta1'}{
							 \describe{
							\item{hyperid = }{`13001'}
							\item{name = }{`log precision'}
							\item{short.name = }{`prec'}
							\item{prior = }{`loggamma'}
							\item{param = }{`1 5e-04'}
							\item{initial = }{`2'}
							\item{fixed = }{`FALSE'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
						\item{Hyperparameter `theta2'}{
							 \describe{
							\item{hyperid = }{`13002'}
							\item{name = }{`logit lambda'}
							\item{short.name = }{`lambda'}
							\item{prior = }{`gaussian'}
							\item{param = }{`0 0.45'}
							\item{initial = }{`3'}
							\item{fixed = }{`FALSE'}
							\item{to.theta = }{`\code{function(x) log(x/(1-x))}'}
							\item{from.theta = }{`\code{function(x) exp(x)/(1+exp(x))}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{doc = }{`An alternative proper version of the Besag model'}
								\item{constr = }{`FALSE'}
								\item{nrow.ncol = }{`FALSE'}
								\item{augmented = }{`FALSE'}
								\item{aug.factor = }{`1'}
								\item{aug.constr = }{`NULL'}
								\item{n.div.by = }{`NULL'}
								\item{n.required = }{`TRUE'}
								\item{set.default.values = }{`TRUE'}
								\item{status = }{`experimental'}
								\item{pdf = }{`besagproper2'}
							 }
						 }
					 }
				 }
				\item{Model `fgn'.}{Number of hyperparmeters are 2.
					 \describe{
						\item{Hyperparameter `theta1'}{
							 \describe{
							\item{hyperid = }{`13101'}
							\item{name = }{`log precision'}
							\item{short.name = }{`prec'}
							\item{prior = }{`pc.prec'}
							\item{param = }{`3 0.01'}
							\item{initial = }{`1'}
							\item{fixed = }{`FALSE'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
						\item{Hyperparameter `theta2'}{
							 \describe{
							\item{hyperid = }{`13102'}
							\item{name = }{`logit H'}
							\item{short.name = }{`H'}
							\item{prior = }{`pcfgnh'}
							\item{param = }{`0.9 0.1'}
							\item{initial = }{`2'}
							\item{fixed = }{`FALSE'}
							\item{to.theta = }{`\code{function(x) log((2*x-1)/(2*(1-x)))}'}
							\item{from.theta = }{`\code{function(x) 0.5 + 0.5*exp(x)/(1+exp(x))}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{doc = }{`Fractional Gaussian noise model'}
								\item{constr = }{`FALSE'}
								\item{nrow.ncol = }{`FALSE'}
								\item{augmented = }{`TRUE'}
								\item{aug.factor = }{`5'}
								\item{aug.constr = }{`1'}
								\item{n.div.by = }{`NULL'}
								\item{n.required = }{`FALSE'}
								\item{set.default.values = }{`TRUE'}
								\item{order.default = }{`4'}
								\item{order.defined = }{`3 4'}
								\item{pdf = }{`fgn'}
							 }
						 }
					 }
				 }
				\item{Model `fgn2'.}{Number of hyperparmeters are 2.
					 \describe{
						\item{Hyperparameter `theta1'}{
							 \describe{
							\item{hyperid = }{`13101'}
							\item{name = }{`log precision'}
							\item{short.name = }{`prec'}
							\item{prior = }{`pc.prec'}
							\item{param = }{`3 0.01'}
							\item{initial = }{`1'}
							\item{fixed = }{`FALSE'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
						\item{Hyperparameter `theta2'}{
							 \describe{
							\item{hyperid = }{`13102'}
							\item{name = }{`logit H'}
							\item{short.name = }{`H'}
							\item{prior = }{`pcfgnh'}
							\item{param = }{`0.9 0.1'}
							\item{initial = }{`2'}
							\item{fixed = }{`FALSE'}
							\item{to.theta = }{`\code{function(x) log((2*x-1)/(2*(1-x)))}'}
							\item{from.theta = }{`\code{function(x) 0.5 + 0.5*exp(x)/(1+exp(x))}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{doc = }{`Fractional Gaussian noise model (alt 2)'}
								\item{constr = }{`FALSE'}
								\item{nrow.ncol = }{`FALSE'}
								\item{augmented = }{`TRUE'}
								\item{aug.factor = }{`4'}
								\item{aug.constr = }{`1'}
								\item{n.div.by = }{`NULL'}
								\item{n.required = }{`FALSE'}
								\item{set.default.values = }{`TRUE'}
								\item{order.default = }{`4'}
								\item{order.defined = }{`3 4'}
								\item{pdf = }{`fgn'}
							 }
						 }
					 }
				 }
				\item{Model `ar1'.}{Number of hyperparmeters are 3.
					 \describe{
						\item{Hyperparameter `theta1'}{
							 \describe{
							\item{hyperid = }{`14001'}
							\item{name = }{`log precision'}
							\item{short.name = }{`prec'}
							\item{prior = }{`loggamma'}
							\item{param = }{`1 5e-05'}
							\item{initial = }{`4'}
							\item{fixed = }{`FALSE'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
						\item{Hyperparameter `theta2'}{
							 \describe{
							\item{hyperid = }{`14002'}
							\item{name = }{`logit lag one correlation'}
							\item{short.name = }{`rho'}
							\item{prior = }{`normal'}
							\item{param = }{`0 0.15'}
							\item{initial = }{`2'}
							\item{fixed = }{`FALSE'}
							\item{to.theta = }{`\code{function(x) log((1+x)/(1-x))}'}
							\item{from.theta = }{`\code{function(x) 2*exp(x)/(1+exp(x))-1}'}
							}
						}
						\item{Hyperparameter `theta3'}{
							 \describe{
							\item{hyperid = }{`14003'}
							\item{name = }{`mean'}
							\item{short.name = }{`mean'}
							\item{prior = }{`normal'}
							\item{param = }{`0 1'}
							\item{initial = }{`0'}
							\item{fixed = }{`TRUE'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{doc = }{`Auto-regressive model of order 1 (AR(1))'}
								\item{constr = }{`FALSE'}
								\item{nrow.ncol = }{`FALSE'}
								\item{augmented = }{`FALSE'}
								\item{aug.factor = }{`1'}
								\item{aug.constr = }{`NULL'}
								\item{n.div.by = }{`NULL'}
								\item{n.required = }{`FALSE'}
								\item{set.default.values = }{`FALSE'}
								\item{pdf = }{`ar1'}
							 }
						 }
					 }
				 }
				\item{Model `ar1c'.}{Number of hyperparmeters are 2.
					 \describe{
						\item{Hyperparameter `theta1'}{
							 \describe{
							\item{hyperid = }{`14101'}
							\item{name = }{`log precision'}
							\item{short.name = }{`prec'}
							\item{prior = }{`pc.prec'}
							\item{param = }{`1 0.01'}
							\item{initial = }{`4'}
							\item{fixed = }{`FALSE'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
						\item{Hyperparameter `theta2'}{
							 \describe{
							\item{hyperid = }{`14102'}
							\item{name = }{`logit lag one correlation'}
							\item{short.name = }{`rho'}
							\item{prior = }{`pc.cor0'}
							\item{param = }{`0.5 0.5'}
							\item{initial = }{`2'}
							\item{fixed = }{`FALSE'}
							\item{to.theta = }{`\code{function(x) log((1+x)/(1-x))}'}
							\item{from.theta = }{`\code{function(x) 2*exp(x)/(1+exp(x))-1}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{doc = }{`Auto-regressive model of order 1 w/covariates'}
								\item{constr = }{`FALSE'}
								\item{nrow.ncol = }{`FALSE'}
								\item{augmented = }{`FALSE'}
								\item{aug.factor = }{`1'}
								\item{aug.constr = }{`NULL'}
								\item{n.div.by = }{`NULL'}
								\item{n.required = }{`FALSE'}
								\item{set.default.values = }{`TRUE'}
								\item{status = }{`experimental'}
								\item{pdf = }{`ar1c'}
							 }
						 }
					 }
				 }
				\item{Model `ar'.}{Number of hyperparmeters are 11.
					 \describe{
						\item{Hyperparameter `theta1'}{
							 \describe{
							\item{hyperid = }{`15001'}
							\item{name = }{`log precision'}
							\item{short.name = }{`prec'}
							\item{initial = }{`4'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`pc.prec'}
							\item{param = }{`3 0.01'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
						\item{Hyperparameter `theta2'}{
							 \describe{
							\item{hyperid = }{`15002'}
							\item{name = }{`pacf1'}
							\item{short.name = }{`pacf1'}
							\item{initial = }{`1'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`pc.cor0'}
							\item{param = }{`0.5 0.5'}
							\item{to.theta = }{`\code{function(x) log((1+x)/(1-x))}'}
							\item{from.theta = }{`\code{function(x) 2*exp(x)/(1+exp(x))-1}'}
							}
						}
						\item{Hyperparameter `theta3'}{
							 \describe{
							\item{hyperid = }{`15003'}
							\item{name = }{`pacf2'}
							\item{short.name = }{`pacf2'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`pc.cor0'}
							\item{param = }{`0.5 0.4'}
							\item{to.theta = }{`\code{function(x) log((1+x)/(1-x))}'}
							\item{from.theta = }{`\code{function(x) 2*exp(x)/(1+exp(x))-1}'}
							}
						}
						\item{Hyperparameter `theta4'}{
							 \describe{
							\item{hyperid = }{`15004'}
							\item{name = }{`pacf3'}
							\item{short.name = }{`pacf3'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`pc.cor0'}
							\item{param = }{`0.5 0.3'}
							\item{to.theta = }{`\code{function(x) log((1+x)/(1-x))}'}
							\item{from.theta = }{`\code{function(x) 2*exp(x)/(1+exp(x))-1}'}
							}
						}
						\item{Hyperparameter `theta5'}{
							 \describe{
							\item{hyperid = }{`15005'}
							\item{name = }{`pacf4'}
							\item{short.name = }{`pacf4'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`pc.cor0'}
							\item{param = }{`0.5 0.2'}
							\item{to.theta = }{`\code{function(x) log((1+x)/(1-x))}'}
							\item{from.theta = }{`\code{function(x) 2*exp(x)/(1+exp(x))-1}'}
							}
						}
						\item{Hyperparameter `theta6'}{
							 \describe{
							\item{hyperid = }{`15006'}
							\item{name = }{`pacf5'}
							\item{short.name = }{`pacf5'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`pc.cor0'}
							\item{param = }{`0.5 0.1'}
							\item{to.theta = }{`\code{function(x) log((1+x)/(1-x))}'}
							\item{from.theta = }{`\code{function(x) 2*exp(x)/(1+exp(x))-1}'}
							}
						}
						\item{Hyperparameter `theta7'}{
							 \describe{
							\item{hyperid = }{`15007'}
							\item{name = }{`pacf6'}
							\item{short.name = }{`pacf6'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`pc.cor0'}
							\item{param = }{`0.5 0.1'}
							\item{to.theta = }{`\code{function(x) log((1+x)/(1-x))}'}
							\item{from.theta = }{`\code{function(x) 2*exp(x)/(1+exp(x))-1}'}
							}
						}
						\item{Hyperparameter `theta8'}{
							 \describe{
							\item{hyperid = }{`15008'}
							\item{name = }{`pacf7'}
							\item{short.name = }{`pacf7'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`pc.cor0'}
							\item{param = }{`0.5 0.1'}
							\item{to.theta = }{`\code{function(x) log((1+x)/(1-x))}'}
							\item{from.theta = }{`\code{function(x) 2*exp(x)/(1+exp(x))-1}'}
							}
						}
						\item{Hyperparameter `theta9'}{
							 \describe{
							\item{hyperid = }{`15009'}
							\item{name = }{`pacf8'}
							\item{short.name = }{`pacf8'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`pc.cor0'}
							\item{param = }{`0.5 0.1'}
							\item{to.theta = }{`\code{function(x) log((1+x)/(1-x))}'}
							\item{from.theta = }{`\code{function(x) 2*exp(x)/(1+exp(x))-1}'}
							}
						}
						\item{Hyperparameter `theta10'}{
							 \describe{
							\item{hyperid = }{`15010'}
							\item{name = }{`pacf9'}
							\item{short.name = }{`pacf9'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`pc.cor0'}
							\item{param = }{`0.5 0.1'}
							\item{to.theta = }{`\code{function(x) log((1+x)/(1-x))}'}
							\item{from.theta = }{`\code{function(x) 2*exp(x)/(1+exp(x))-1}'}
							}
						}
						\item{Hyperparameter `theta11'}{
							 \describe{
							\item{hyperid = }{`15011'}
							\item{name = }{`pacf10'}
							\item{short.name = }{`pacf10'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`pc.cor0'}
							\item{param = }{`0.5 0.1'}
							\item{to.theta = }{`\code{function(x) log((1+x)/(1-x))}'}
							\item{from.theta = }{`\code{function(x) 2*exp(x)/(1+exp(x))-1}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{doc = }{`Auto-regressive model of order p (AR(p))'}
								\item{constr = }{`FALSE'}
								\item{nrow.ncol = }{`FALSE'}
								\item{augmented = }{`FALSE'}
								\item{aug.factor = }{`1'}
								\item{aug.constr = }{`NULL'}
								\item{n.div.by = }{`NULL'}
								\item{n.required = }{`FALSE'}
								\item{set.default.values = }{`FALSE'}
								\item{pdf = }{`ar'}
							 }
						 }
					 }
				 }
				\item{Model `ou'.}{Number of hyperparmeters are 2.
					 \describe{
						\item{Hyperparameter `theta1'}{
							 \describe{
							\item{hyperid = }{`16001'}
							\item{name = }{`log precision'}
							\item{short.name = }{`prec'}
							\item{prior = }{`loggamma'}
							\item{param = }{`1 5e-05'}
							\item{initial = }{`4'}
							\item{fixed = }{`FALSE'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
						\item{Hyperparameter `theta2'}{
							 \describe{
							\item{hyperid = }{`16002'}
							\item{name = }{`log phi'}
							\item{short.name = }{`phi'}
							\item{prior = }{`normal'}
							\item{param = }{`0 0.2'}
							\item{initial = }{`-1'}
							\item{fixed = }{`FALSE'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{doc = }{`The Ornstein-Uhlenbeck process'}
								\item{constr = }{`FALSE'}
								\item{nrow.ncol = }{`FALSE'}
								\item{augmented = }{`FALSE'}
								\item{aug.factor = }{`1'}
								\item{aug.constr = }{`NULL'}
								\item{n.div.by = }{`NULL'}
								\item{n.required = }{`FALSE'}
								\item{set.default.values = }{`FALSE'}
								\item{pdf = }{`ou'}
							 }
						 }
					 }
				 }
				\item{Model `intslope'.}{Number of hyperparmeters are 13.
					 \describe{
						\item{Hyperparameter `theta1'}{
							 \describe{
							\item{hyperid = }{`16101'}
							\item{name = }{`log precision1'}
							\item{short.name = }{`prec1'}
							\item{initial = }{`4'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`wishart2d'}
							\item{param = }{`4 1 1 0'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
						\item{Hyperparameter `theta2'}{
							 \describe{
							\item{hyperid = }{`16102'}
							\item{name = }{`log precision2'}
							\item{short.name = }{`prec2'}
							\item{initial = }{`4'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
						\item{Hyperparameter `theta3'}{
							 \describe{
							\item{hyperid = }{`16103'}
							\item{name = }{`logit correlation'}
							\item{short.name = }{`cor'}
							\item{initial = }{`4'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) log((1+x)/(1-x))}'}
							\item{from.theta = }{`\code{function(x) 2*exp(x)/(1+exp(x))-1}'}
							}
						}
						\item{Hyperparameter `theta4'}{
							 \describe{
							\item{hyperid = }{`16104'}
							\item{name = }{`gamma1'}
							\item{short.name = }{`g1'}
							\item{initial = }{`1'}
							\item{fixed = }{`TRUE'}
							\item{prior = }{`normal'}
							\item{param = }{`1 36'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta5'}{
							 \describe{
							\item{hyperid = }{`16105'}
							\item{name = }{`gamma2'}
							\item{short.name = }{`g2'}
							\item{initial = }{`1'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`normal'}
							\item{param = }{`1 36'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta6'}{
							 \describe{
							\item{hyperid = }{`16106'}
							\item{name = }{`gamma3'}
							\item{short.name = }{`g3'}
							\item{initial = }{`1'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`normal'}
							\item{param = }{`1 36'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta7'}{
							 \describe{
							\item{hyperid = }{`16107'}
							\item{name = }{`gamma4'}
							\item{short.name = }{`g4'}
							\item{initial = }{`1'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`normal'}
							\item{param = }{`1 36'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta8'}{
							 \describe{
							\item{hyperid = }{`16108'}
							\item{name = }{`gamma5'}
							\item{short.name = }{`g5'}
							\item{initial = }{`1'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`normal'}
							\item{param = }{`1 36'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta9'}{
							 \describe{
							\item{hyperid = }{`16109'}
							\item{name = }{`gamma6'}
							\item{short.name = }{`g6'}
							\item{initial = }{`1'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`normal'}
							\item{param = }{`1 36'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta10'}{
							 \describe{
							\item{hyperid = }{`16110'}
							\item{name = }{`gamma7'}
							\item{short.name = }{`g7'}
							\item{initial = }{`1'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`normal'}
							\item{param = }{`1 36'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta11'}{
							 \describe{
							\item{hyperid = }{`16111'}
							\item{name = }{`gamma8'}
							\item{short.name = }{`g8'}
							\item{initial = }{`1'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`normal'}
							\item{param = }{`1 36'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta12'}{
							 \describe{
							\item{hyperid = }{`16112'}
							\item{name = }{`gamma9'}
							\item{short.name = }{`g9'}
							\item{initial = }{`1'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`normal'}
							\item{param = }{`1 36'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta13'}{
							 \describe{
							\item{hyperid = }{`16113'}
							\item{name = }{`gamma10'}
							\item{short.name = }{`g10'}
							\item{initial = }{`1'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`normal'}
							\item{param = }{`1 36'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{doc = }{`Intecept-slope model with Wishart-prior'}
								\item{constr = }{`FALSE'}
								\item{nrow.ncol = }{`FALSE'}
								\item{augmented = }{`FALSE'}
								\item{aug.factor = }{`1'}
								\item{aug.constr = }{`NULL'}
								\item{n.div.by = }{`NULL'}
								\item{n.required = }{`FALSE'}
								\item{set.default.values = }{`TRUE'}
								\item{status = }{`experimental'}
								\item{pdf = }{`intslope'}
							 }
						 }
					 }
				 }
				\item{Model `generic'.}{Number of hyperparmeters are 1.
					 \describe{
						\item{Hyperparameter `theta'}{
							 \describe{
							\item{hyperid = }{`17001'}
							\item{name = }{`log precision'}
							\item{short.name = }{`prec'}
							\item{prior = }{`loggamma'}
							\item{param = }{`1 5e-05'}
							\item{initial = }{`4'}
							\item{fixed = }{`FALSE'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{doc = }{`A generic model'}
								\item{constr = }{`FALSE'}
								\item{nrow.ncol = }{`FALSE'}
								\item{augmented = }{`FALSE'}
								\item{aug.factor = }{`1'}
								\item{aug.constr = }{`NULL'}
								\item{n.div.by = }{`NULL'}
								\item{n.required = }{`TRUE'}
								\item{set.default.values = }{`TRUE'}
								\item{pdf = }{`generic0'}
							 }
						 }
					 }
				 }
				\item{Model `generic0'.}{Number of hyperparmeters are 1.
					 \describe{
						\item{Hyperparameter `theta'}{
							 \describe{
							\item{hyperid = }{`18001'}
							\item{name = }{`log precision'}
							\item{short.name = }{`prec'}
							\item{prior = }{`loggamma'}
							\item{param = }{`1 5e-05'}
							\item{initial = }{`4'}
							\item{fixed = }{`FALSE'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{doc = }{`A generic model (type 0)'}
								\item{constr = }{`FALSE'}
								\item{nrow.ncol = }{`FALSE'}
								\item{augmented = }{`FALSE'}
								\item{aug.factor = }{`1'}
								\item{aug.constr = }{`NULL'}
								\item{n.div.by = }{`NULL'}
								\item{n.required = }{`TRUE'}
								\item{set.default.values = }{`TRUE'}
								\item{pdf = }{`generic0'}
							 }
						 }
					 }
				 }
				\item{Model `generic1'.}{Number of hyperparmeters are 2.
					 \describe{
						\item{Hyperparameter `theta1'}{
							 \describe{
							\item{hyperid = }{`19001'}
							\item{name = }{`log precision'}
							\item{short.name = }{`prec'}
							\item{prior = }{`loggamma'}
							\item{param = }{`1 5e-05'}
							\item{initial = }{`4'}
							\item{fixed = }{`FALSE'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
						\item{Hyperparameter `theta2'}{
							 \describe{
							\item{hyperid = }{`19002'}
							\item{name = }{`beta'}
							\item{short.name = }{`beta'}
							\item{initial = }{`2'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`gaussian'}
							\item{param = }{`0 0.1'}
							\item{to.theta = }{`\code{function(x) log(x/(1-x))}'}
							\item{from.theta = }{`\code{function(x) exp(x)/(1+exp(x))}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{doc = }{`A generic model (type 1)'}
								\item{constr = }{`FALSE'}
								\item{nrow.ncol = }{`FALSE'}
								\item{augmented = }{`FALSE'}
								\item{aug.factor = }{`1'}
								\item{aug.constr = }{`NULL'}
								\item{n.div.by = }{`NULL'}
								\item{n.required = }{`TRUE'}
								\item{set.default.values = }{`TRUE'}
								\item{pdf = }{`generic1'}
							 }
						 }
					 }
				 }
				\item{Model `generic2'.}{Number of hyperparmeters are 2.
					 \describe{
						\item{Hyperparameter `theta1'}{
							 \describe{
							\item{hyperid = }{`20001'}
							\item{name = }{`log precision cmatrix'}
							\item{short.name = }{`prec'}
							\item{initial = }{`4'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`loggamma'}
							\item{param = }{`1 5e-05'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
						\item{Hyperparameter `theta2'}{
							 \describe{
							\item{hyperid = }{`20002'}
							\item{name = }{`log precision random'}
							\item{short.name = }{`prec.random'}
							\item{initial = }{`4'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`loggamma'}
							\item{param = }{`1 0.001'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{doc = }{`A generic model (type 2)'}
								\item{constr = }{`FALSE'}
								\item{nrow.ncol = }{`FALSE'}
								\item{augmented = }{`FALSE'}
								\item{aug.factor = }{`2'}
								\item{aug.constr = }{`2'}
								\item{n.div.by = }{`NULL'}
								\item{n.required = }{`TRUE'}
								\item{set.default.values = }{`TRUE'}
								\item{pdf = }{`generic2'}
							 }
						 }
					 }
				 }
				\item{Model `generic3'.}{Number of hyperparmeters are 11.
					 \describe{
						\item{Hyperparameter `theta1'}{
							 \describe{
							\item{hyperid = }{`21001'}
							\item{name = }{`log precision1'}
							\item{short.name = }{`prec1'}
							\item{initial = }{`4'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`loggamma'}
							\item{param = }{`1 5e-05'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
						\item{Hyperparameter `theta2'}{
							 \describe{
							\item{hyperid = }{`21002'}
							\item{name = }{`log precision2'}
							\item{short.name = }{`prec2'}
							\item{initial = }{`4'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`loggamma'}
							\item{param = }{`1 5e-05'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
						\item{Hyperparameter `theta3'}{
							 \describe{
							\item{hyperid = }{`21003'}
							\item{name = }{`log precision3'}
							\item{short.name = }{`prec3'}
							\item{initial = }{`4'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`loggamma'}
							\item{param = }{`1 5e-05'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
						\item{Hyperparameter `theta4'}{
							 \describe{
							\item{hyperid = }{`21004'}
							\item{name = }{`log precision4'}
							\item{short.name = }{`prec4'}
							\item{initial = }{`4'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`loggamma'}
							\item{param = }{`1 5e-05'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
						\item{Hyperparameter `theta5'}{
							 \describe{
							\item{hyperid = }{`21005'}
							\item{name = }{`log precision5'}
							\item{short.name = }{`prec5'}
							\item{initial = }{`4'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`loggamma'}
							\item{param = }{`1 5e-05'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
						\item{Hyperparameter `theta6'}{
							 \describe{
							\item{hyperid = }{`21006'}
							\item{name = }{`log precision6'}
							\item{short.name = }{`prec6'}
							\item{initial = }{`4'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`loggamma'}
							\item{param = }{`1 5e-05'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
						\item{Hyperparameter `theta7'}{
							 \describe{
							\item{hyperid = }{`21007'}
							\item{name = }{`log precision7'}
							\item{short.name = }{`prec7'}
							\item{initial = }{`4'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`loggamma'}
							\item{param = }{`1 5e-05'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
						\item{Hyperparameter `theta8'}{
							 \describe{
							\item{hyperid = }{`21008'}
							\item{name = }{`log precision8'}
							\item{short.name = }{`prec8'}
							\item{initial = }{`4'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`loggamma'}
							\item{param = }{`1 5e-05'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
						\item{Hyperparameter `theta9'}{
							 \describe{
							\item{hyperid = }{`21009'}
							\item{name = }{`log precision9'}
							\item{short.name = }{`prec9'}
							\item{initial = }{`4'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`loggamma'}
							\item{param = }{`1 5e-05'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
						\item{Hyperparameter `theta10'}{
							 \describe{
							\item{hyperid = }{`21010'}
							\item{name = }{`log precision10'}
							\item{short.name = }{`prec10'}
							\item{initial = }{`4'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`loggamma'}
							\item{param = }{`1 5e-05'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
						\item{Hyperparameter `theta11'}{
							 \describe{
							\item{hyperid = }{`21011'}
							\item{name = }{`log precision common'}
							\item{short.name = }{`prec.common'}
							\item{initial = }{`0'}
							\item{fixed = }{`TRUE'}
							\item{prior = }{`loggamma'}
							\item{param = }{`1 5e-05'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{doc = }{`A generic model (type 3)'}
								\item{constr = }{`FALSE'}
								\item{nrow.ncol = }{`FALSE'}
								\item{augmented = }{`FALSE'}
								\item{aug.factor = }{`1'}
								\item{aug.constr = }{`NULL'}
								\item{n.div.by = }{`NULL'}
								\item{n.required = }{`TRUE'}
								\item{set.default.values = }{`TRUE'}
								\item{status = }{`experimental'}
								\item{pdf = }{`generic3'}
							 }
						 }
					 }
				 }
				\item{Model `spde'.}{Number of hyperparmeters are 4.
					 \describe{
						\item{Hyperparameter `theta1'}{
							 \describe{
							\item{hyperid = }{`22001'}
							\item{name = }{`theta.T'}
							\item{short.name = }{`T'}
							\item{initial = }{`2'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`normal'}
							\item{param = }{`0 1'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta2'}{
							 \describe{
							\item{hyperid = }{`22002'}
							\item{name = }{`theta.K'}
							\item{short.name = }{`K'}
							\item{initial = }{`-2'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`normal'}
							\item{param = }{`0 1'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta3'}{
							 \describe{
							\item{hyperid = }{`22003'}
							\item{name = }{`theta.KT'}
							\item{short.name = }{`KT'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`normal'}
							\item{param = }{`0 1'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta4'}{
							 \describe{
							\item{hyperid = }{`22004'}
							\item{name = }{`theta.OC'}
							\item{short.name = }{`OC'}
							\item{initial = }{`-20'}
							\item{fixed = }{`TRUE'}
							\item{prior = }{`normal'}
							\item{param = }{`0 0.2'}
							\item{to.theta = }{`\code{function(x) log(x/(1-x))}'}
							\item{from.theta = }{`\code{function(x) exp(x)/(1+exp(x))}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{doc = }{`A SPDE model'}
								\item{constr = }{`FALSE'}
								\item{nrow.ncol = }{`FALSE'}
								\item{augmented = }{`FALSE'}
								\item{aug.factor = }{`1'}
								\item{aug.constr = }{`NULL'}
								\item{n.div.by = }{`NULL'}
								\item{n.required = }{`TRUE'}
								\item{set.default.values = }{`TRUE'}
								\item{pdf = }{`spde'}
							 }
						 }
					 }
				 }
				\item{Model `spde2'.}{Number of hyperparmeters are 100.
					 \describe{
						\item{Hyperparameter `theta1'}{
							 \describe{
							\item{hyperid = }{`23001'}
							\item{name = }{`theta1'}
							\item{short.name = }{`t1'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`mvnorm'}
							\item{param = }{`1 1'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta2'}{
							 \describe{
							\item{hyperid = }{`23002'}
							\item{name = }{`theta2'}
							\item{short.name = }{`t2'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta3'}{
							 \describe{
							\item{hyperid = }{`23003'}
							\item{name = }{`theta3'}
							\item{short.name = }{`t3'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta4'}{
							 \describe{
							\item{hyperid = }{`23004'}
							\item{name = }{`theta4'}
							\item{short.name = }{`t4'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta5'}{
							 \describe{
							\item{hyperid = }{`23005'}
							\item{name = }{`theta5'}
							\item{short.name = }{`t5'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta6'}{
							 \describe{
							\item{hyperid = }{`23006'}
							\item{name = }{`theta6'}
							\item{short.name = }{`t6'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta7'}{
							 \describe{
							\item{hyperid = }{`23007'}
							\item{name = }{`theta7'}
							\item{short.name = }{`t7'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta8'}{
							 \describe{
							\item{hyperid = }{`23008'}
							\item{name = }{`theta8'}
							\item{short.name = }{`t8'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta9'}{
							 \describe{
							\item{hyperid = }{`23009'}
							\item{name = }{`theta9'}
							\item{short.name = }{`t9'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta10'}{
							 \describe{
							\item{hyperid = }{`23010'}
							\item{name = }{`theta10'}
							\item{short.name = }{`t10'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta11'}{
							 \describe{
							\item{hyperid = }{`23011'}
							\item{name = }{`theta11'}
							\item{short.name = }{`t11'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta12'}{
							 \describe{
							\item{hyperid = }{`23012'}
							\item{name = }{`theta12'}
							\item{short.name = }{`t12'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta13'}{
							 \describe{
							\item{hyperid = }{`23013'}
							\item{name = }{`theta13'}
							\item{short.name = }{`t13'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta14'}{
							 \describe{
							\item{hyperid = }{`23014'}
							\item{name = }{`theta14'}
							\item{short.name = }{`t14'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta15'}{
							 \describe{
							\item{hyperid = }{`23015'}
							\item{name = }{`theta15'}
							\item{short.name = }{`t15'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta16'}{
							 \describe{
							\item{hyperid = }{`23016'}
							\item{name = }{`theta16'}
							\item{short.name = }{`t16'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta17'}{
							 \describe{
							\item{hyperid = }{`23017'}
							\item{name = }{`theta17'}
							\item{short.name = }{`t17'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta18'}{
							 \describe{
							\item{hyperid = }{`23018'}
							\item{name = }{`theta18'}
							\item{short.name = }{`t18'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta19'}{
							 \describe{
							\item{hyperid = }{`23019'}
							\item{name = }{`theta19'}
							\item{short.name = }{`t19'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta20'}{
							 \describe{
							\item{hyperid = }{`23020'}
							\item{name = }{`theta20'}
							\item{short.name = }{`t20'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta21'}{
							 \describe{
							\item{hyperid = }{`23021'}
							\item{name = }{`theta21'}
							\item{short.name = }{`t21'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta22'}{
							 \describe{
							\item{hyperid = }{`23022'}
							\item{name = }{`theta22'}
							\item{short.name = }{`t22'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta23'}{
							 \describe{
							\item{hyperid = }{`23023'}
							\item{name = }{`theta23'}
							\item{short.name = }{`t23'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta24'}{
							 \describe{
							\item{hyperid = }{`23024'}
							\item{name = }{`theta24'}
							\item{short.name = }{`t24'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta25'}{
							 \describe{
							\item{hyperid = }{`23025'}
							\item{name = }{`theta25'}
							\item{short.name = }{`t25'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta26'}{
							 \describe{
							\item{hyperid = }{`23026'}
							\item{name = }{`theta26'}
							\item{short.name = }{`t26'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta27'}{
							 \describe{
							\item{hyperid = }{`23027'}
							\item{name = }{`theta27'}
							\item{short.name = }{`t27'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta28'}{
							 \describe{
							\item{hyperid = }{`23028'}
							\item{name = }{`theta28'}
							\item{short.name = }{`t28'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta29'}{
							 \describe{
							\item{hyperid = }{`23029'}
							\item{name = }{`theta29'}
							\item{short.name = }{`t29'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta30'}{
							 \describe{
							\item{hyperid = }{`23030'}
							\item{name = }{`theta30'}
							\item{short.name = }{`t30'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta31'}{
							 \describe{
							\item{hyperid = }{`23031'}
							\item{name = }{`theta31'}
							\item{short.name = }{`t31'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta32'}{
							 \describe{
							\item{hyperid = }{`23032'}
							\item{name = }{`theta32'}
							\item{short.name = }{`t32'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta33'}{
							 \describe{
							\item{hyperid = }{`23033'}
							\item{name = }{`theta33'}
							\item{short.name = }{`t33'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta34'}{
							 \describe{
							\item{hyperid = }{`23034'}
							\item{name = }{`theta34'}
							\item{short.name = }{`t34'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta35'}{
							 \describe{
							\item{hyperid = }{`23035'}
							\item{name = }{`theta35'}
							\item{short.name = }{`t35'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta36'}{
							 \describe{
							\item{hyperid = }{`23036'}
							\item{name = }{`theta36'}
							\item{short.name = }{`t36'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta37'}{
							 \describe{
							\item{hyperid = }{`23037'}
							\item{name = }{`theta37'}
							\item{short.name = }{`t37'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta38'}{
							 \describe{
							\item{hyperid = }{`23038'}
							\item{name = }{`theta38'}
							\item{short.name = }{`t38'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta39'}{
							 \describe{
							\item{hyperid = }{`23039'}
							\item{name = }{`theta39'}
							\item{short.name = }{`t39'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta40'}{
							 \describe{
							\item{hyperid = }{`23040'}
							\item{name = }{`theta40'}
							\item{short.name = }{`t40'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta41'}{
							 \describe{
							\item{hyperid = }{`23041'}
							\item{name = }{`theta41'}
							\item{short.name = }{`t41'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta42'}{
							 \describe{
							\item{hyperid = }{`23042'}
							\item{name = }{`theta42'}
							\item{short.name = }{`t42'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta43'}{
							 \describe{
							\item{hyperid = }{`23043'}
							\item{name = }{`theta43'}
							\item{short.name = }{`t43'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta44'}{
							 \describe{
							\item{hyperid = }{`23044'}
							\item{name = }{`theta44'}
							\item{short.name = }{`t44'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta45'}{
							 \describe{
							\item{hyperid = }{`23045'}
							\item{name = }{`theta45'}
							\item{short.name = }{`t45'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta46'}{
							 \describe{
							\item{hyperid = }{`23046'}
							\item{name = }{`theta46'}
							\item{short.name = }{`t46'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta47'}{
							 \describe{
							\item{hyperid = }{`23047'}
							\item{name = }{`theta47'}
							\item{short.name = }{`t47'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta48'}{
							 \describe{
							\item{hyperid = }{`23048'}
							\item{name = }{`theta48'}
							\item{short.name = }{`t48'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta49'}{
							 \describe{
							\item{hyperid = }{`23049'}
							\item{name = }{`theta49'}
							\item{short.name = }{`t49'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta50'}{
							 \describe{
							\item{hyperid = }{`23050'}
							\item{name = }{`theta50'}
							\item{short.name = }{`t50'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta51'}{
							 \describe{
							\item{hyperid = }{`23051'}
							\item{name = }{`theta51'}
							\item{short.name = }{`t51'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta52'}{
							 \describe{
							\item{hyperid = }{`23052'}
							\item{name = }{`theta52'}
							\item{short.name = }{`t52'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta53'}{
							 \describe{
							\item{hyperid = }{`23053'}
							\item{name = }{`theta53'}
							\item{short.name = }{`t53'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta54'}{
							 \describe{
							\item{hyperid = }{`23054'}
							\item{name = }{`theta54'}
							\item{short.name = }{`t54'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta55'}{
							 \describe{
							\item{hyperid = }{`23055'}
							\item{name = }{`theta55'}
							\item{short.name = }{`t55'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta56'}{
							 \describe{
							\item{hyperid = }{`23056'}
							\item{name = }{`theta56'}
							\item{short.name = }{`t56'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta57'}{
							 \describe{
							\item{hyperid = }{`23057'}
							\item{name = }{`theta57'}
							\item{short.name = }{`t57'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta58'}{
							 \describe{
							\item{hyperid = }{`23058'}
							\item{name = }{`theta58'}
							\item{short.name = }{`t58'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta59'}{
							 \describe{
							\item{hyperid = }{`23059'}
							\item{name = }{`theta59'}
							\item{short.name = }{`t59'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta60'}{
							 \describe{
							\item{hyperid = }{`23060'}
							\item{name = }{`theta60'}
							\item{short.name = }{`t60'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta61'}{
							 \describe{
							\item{hyperid = }{`23061'}
							\item{name = }{`theta61'}
							\item{short.name = }{`t61'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta62'}{
							 \describe{
							\item{hyperid = }{`23062'}
							\item{name = }{`theta62'}
							\item{short.name = }{`t62'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta63'}{
							 \describe{
							\item{hyperid = }{`23063'}
							\item{name = }{`theta63'}
							\item{short.name = }{`t63'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta64'}{
							 \describe{
							\item{hyperid = }{`23064'}
							\item{name = }{`theta64'}
							\item{short.name = }{`t64'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta65'}{
							 \describe{
							\item{hyperid = }{`23065'}
							\item{name = }{`theta65'}
							\item{short.name = }{`t65'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta66'}{
							 \describe{
							\item{hyperid = }{`23066'}
							\item{name = }{`theta66'}
							\item{short.name = }{`t66'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta67'}{
							 \describe{
							\item{hyperid = }{`23067'}
							\item{name = }{`theta67'}
							\item{short.name = }{`t67'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta68'}{
							 \describe{
							\item{hyperid = }{`23068'}
							\item{name = }{`theta68'}
							\item{short.name = }{`t68'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta69'}{
							 \describe{
							\item{hyperid = }{`23069'}
							\item{name = }{`theta69'}
							\item{short.name = }{`t69'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta70'}{
							 \describe{
							\item{hyperid = }{`23070'}
							\item{name = }{`theta70'}
							\item{short.name = }{`t70'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta71'}{
							 \describe{
							\item{hyperid = }{`23071'}
							\item{name = }{`theta71'}
							\item{short.name = }{`t71'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta72'}{
							 \describe{
							\item{hyperid = }{`23072'}
							\item{name = }{`theta72'}
							\item{short.name = }{`t72'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta73'}{
							 \describe{
							\item{hyperid = }{`23073'}
							\item{name = }{`theta73'}
							\item{short.name = }{`t73'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta74'}{
							 \describe{
							\item{hyperid = }{`23074'}
							\item{name = }{`theta74'}
							\item{short.name = }{`t74'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta75'}{
							 \describe{
							\item{hyperid = }{`23075'}
							\item{name = }{`theta75'}
							\item{short.name = }{`t75'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta76'}{
							 \describe{
							\item{hyperid = }{`23076'}
							\item{name = }{`theta76'}
							\item{short.name = }{`t76'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta77'}{
							 \describe{
							\item{hyperid = }{`23077'}
							\item{name = }{`theta77'}
							\item{short.name = }{`t77'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta78'}{
							 \describe{
							\item{hyperid = }{`23078'}
							\item{name = }{`theta78'}
							\item{short.name = }{`t78'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta79'}{
							 \describe{
							\item{hyperid = }{`23079'}
							\item{name = }{`theta79'}
							\item{short.name = }{`t79'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta80'}{
							 \describe{
							\item{hyperid = }{`23080'}
							\item{name = }{`theta80'}
							\item{short.name = }{`t80'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta81'}{
							 \describe{
							\item{hyperid = }{`23081'}
							\item{name = }{`theta81'}
							\item{short.name = }{`t81'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta82'}{
							 \describe{
							\item{hyperid = }{`23082'}
							\item{name = }{`theta82'}
							\item{short.name = }{`t82'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta83'}{
							 \describe{
							\item{hyperid = }{`23083'}
							\item{name = }{`theta83'}
							\item{short.name = }{`t83'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta84'}{
							 \describe{
							\item{hyperid = }{`23084'}
							\item{name = }{`theta84'}
							\item{short.name = }{`t84'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta85'}{
							 \describe{
							\item{hyperid = }{`23085'}
							\item{name = }{`theta85'}
							\item{short.name = }{`t85'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta86'}{
							 \describe{
							\item{hyperid = }{`23086'}
							\item{name = }{`theta86'}
							\item{short.name = }{`t86'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta87'}{
							 \describe{
							\item{hyperid = }{`23087'}
							\item{name = }{`theta87'}
							\item{short.name = }{`t87'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta88'}{
							 \describe{
							\item{hyperid = }{`23088'}
							\item{name = }{`theta88'}
							\item{short.name = }{`t88'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta89'}{
							 \describe{
							\item{hyperid = }{`23089'}
							\item{name = }{`theta89'}
							\item{short.name = }{`t89'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta90'}{
							 \describe{
							\item{hyperid = }{`23090'}
							\item{name = }{`theta90'}
							\item{short.name = }{`t90'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta91'}{
							 \describe{
							\item{hyperid = }{`23091'}
							\item{name = }{`theta91'}
							\item{short.name = }{`t91'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta92'}{
							 \describe{
							\item{hyperid = }{`23092'}
							\item{name = }{`theta92'}
							\item{short.name = }{`t92'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta93'}{
							 \describe{
							\item{hyperid = }{`23093'}
							\item{name = }{`theta93'}
							\item{short.name = }{`t93'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta94'}{
							 \describe{
							\item{hyperid = }{`23094'}
							\item{name = }{`theta94'}
							\item{short.name = }{`t94'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta95'}{
							 \describe{
							\item{hyperid = }{`23095'}
							\item{name = }{`theta95'}
							\item{short.name = }{`t95'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta96'}{
							 \describe{
							\item{hyperid = }{`23096'}
							\item{name = }{`theta96'}
							\item{short.name = }{`t96'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta97'}{
							 \describe{
							\item{hyperid = }{`23097'}
							\item{name = }{`theta97'}
							\item{short.name = }{`t97'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta98'}{
							 \describe{
							\item{hyperid = }{`23098'}
							\item{name = }{`theta98'}
							\item{short.name = }{`t98'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta99'}{
							 \describe{
							\item{hyperid = }{`23099'}
							\item{name = }{`theta99'}
							\item{short.name = }{`t99'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta100'}{
							 \describe{
							\item{hyperid = }{`23100'}
							\item{name = }{`theta100'}
							\item{short.name = }{`t100'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{doc = }{`A SPDE2 model'}
								\item{constr = }{`FALSE'}
								\item{nrow.ncol = }{`FALSE'}
								\item{augmented = }{`FALSE'}
								\item{aug.factor = }{`1'}
								\item{aug.constr = }{`NULL'}
								\item{n.div.by = }{`NULL'}
								\item{n.required = }{`TRUE'}
								\item{set.default.values = }{`TRUE'}
								\item{pdf = }{`spde2'}
							 }
						 }
					 }
				 }
				\item{Model `spde3'.}{Number of hyperparmeters are 100.
					 \describe{
						\item{Hyperparameter `theta1'}{
							 \describe{
							\item{hyperid = }{`24001'}
							\item{name = }{`theta1'}
							\item{short.name = }{`t1'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`mvnorm'}
							\item{param = }{`1 1'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta2'}{
							 \describe{
							\item{hyperid = }{`24002'}
							\item{name = }{`theta2'}
							\item{short.name = }{`t2'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta3'}{
							 \describe{
							\item{hyperid = }{`24003'}
							\item{name = }{`theta3'}
							\item{short.name = }{`t3'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta4'}{
							 \describe{
							\item{hyperid = }{`24004'}
							\item{name = }{`theta4'}
							\item{short.name = }{`t4'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta5'}{
							 \describe{
							\item{hyperid = }{`24005'}
							\item{name = }{`theta5'}
							\item{short.name = }{`t5'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta6'}{
							 \describe{
							\item{hyperid = }{`24006'}
							\item{name = }{`theta6'}
							\item{short.name = }{`t6'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta7'}{
							 \describe{
							\item{hyperid = }{`24007'}
							\item{name = }{`theta7'}
							\item{short.name = }{`t7'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta8'}{
							 \describe{
							\item{hyperid = }{`24008'}
							\item{name = }{`theta8'}
							\item{short.name = }{`t8'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta9'}{
							 \describe{
							\item{hyperid = }{`24009'}
							\item{name = }{`theta9'}
							\item{short.name = }{`t9'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta10'}{
							 \describe{
							\item{hyperid = }{`24010'}
							\item{name = }{`theta10'}
							\item{short.name = }{`t10'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta11'}{
							 \describe{
							\item{hyperid = }{`24011'}
							\item{name = }{`theta11'}
							\item{short.name = }{`t11'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta12'}{
							 \describe{
							\item{hyperid = }{`24012'}
							\item{name = }{`theta12'}
							\item{short.name = }{`t12'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta13'}{
							 \describe{
							\item{hyperid = }{`24013'}
							\item{name = }{`theta13'}
							\item{short.name = }{`t13'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta14'}{
							 \describe{
							\item{hyperid = }{`24014'}
							\item{name = }{`theta14'}
							\item{short.name = }{`t14'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta15'}{
							 \describe{
							\item{hyperid = }{`24015'}
							\item{name = }{`theta15'}
							\item{short.name = }{`t15'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta16'}{
							 \describe{
							\item{hyperid = }{`24016'}
							\item{name = }{`theta16'}
							\item{short.name = }{`t16'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta17'}{
							 \describe{
							\item{hyperid = }{`24017'}
							\item{name = }{`theta17'}
							\item{short.name = }{`t17'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta18'}{
							 \describe{
							\item{hyperid = }{`24018'}
							\item{name = }{`theta18'}
							\item{short.name = }{`t18'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta19'}{
							 \describe{
							\item{hyperid = }{`24019'}
							\item{name = }{`theta19'}
							\item{short.name = }{`t19'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta20'}{
							 \describe{
							\item{hyperid = }{`24020'}
							\item{name = }{`theta20'}
							\item{short.name = }{`t20'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta21'}{
							 \describe{
							\item{hyperid = }{`24021'}
							\item{name = }{`theta21'}
							\item{short.name = }{`t21'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta22'}{
							 \describe{
							\item{hyperid = }{`24022'}
							\item{name = }{`theta22'}
							\item{short.name = }{`t22'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta23'}{
							 \describe{
							\item{hyperid = }{`24023'}
							\item{name = }{`theta23'}
							\item{short.name = }{`t23'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta24'}{
							 \describe{
							\item{hyperid = }{`24024'}
							\item{name = }{`theta24'}
							\item{short.name = }{`t24'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta25'}{
							 \describe{
							\item{hyperid = }{`24025'}
							\item{name = }{`theta25'}
							\item{short.name = }{`t25'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta26'}{
							 \describe{
							\item{hyperid = }{`24026'}
							\item{name = }{`theta26'}
							\item{short.name = }{`t26'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta27'}{
							 \describe{
							\item{hyperid = }{`24027'}
							\item{name = }{`theta27'}
							\item{short.name = }{`t27'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta28'}{
							 \describe{
							\item{hyperid = }{`24028'}
							\item{name = }{`theta28'}
							\item{short.name = }{`t28'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta29'}{
							 \describe{
							\item{hyperid = }{`24029'}
							\item{name = }{`theta29'}
							\item{short.name = }{`t29'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta30'}{
							 \describe{
							\item{hyperid = }{`24030'}
							\item{name = }{`theta30'}
							\item{short.name = }{`t30'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta31'}{
							 \describe{
							\item{hyperid = }{`24031'}
							\item{name = }{`theta31'}
							\item{short.name = }{`t31'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta32'}{
							 \describe{
							\item{hyperid = }{`24032'}
							\item{name = }{`theta32'}
							\item{short.name = }{`t32'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta33'}{
							 \describe{
							\item{hyperid = }{`24033'}
							\item{name = }{`theta33'}
							\item{short.name = }{`t33'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta34'}{
							 \describe{
							\item{hyperid = }{`24034'}
							\item{name = }{`theta34'}
							\item{short.name = }{`t34'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta35'}{
							 \describe{
							\item{hyperid = }{`24035'}
							\item{name = }{`theta35'}
							\item{short.name = }{`t35'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta36'}{
							 \describe{
							\item{hyperid = }{`24036'}
							\item{name = }{`theta36'}
							\item{short.name = }{`t36'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta37'}{
							 \describe{
							\item{hyperid = }{`24037'}
							\item{name = }{`theta37'}
							\item{short.name = }{`t37'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta38'}{
							 \describe{
							\item{hyperid = }{`24038'}
							\item{name = }{`theta38'}
							\item{short.name = }{`t38'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta39'}{
							 \describe{
							\item{hyperid = }{`24039'}
							\item{name = }{`theta39'}
							\item{short.name = }{`t39'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta40'}{
							 \describe{
							\item{hyperid = }{`24040'}
							\item{name = }{`theta40'}
							\item{short.name = }{`t40'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta41'}{
							 \describe{
							\item{hyperid = }{`24041'}
							\item{name = }{`theta41'}
							\item{short.name = }{`t41'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta42'}{
							 \describe{
							\item{hyperid = }{`24042'}
							\item{name = }{`theta42'}
							\item{short.name = }{`t42'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta43'}{
							 \describe{
							\item{hyperid = }{`24043'}
							\item{name = }{`theta43'}
							\item{short.name = }{`t43'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta44'}{
							 \describe{
							\item{hyperid = }{`24044'}
							\item{name = }{`theta44'}
							\item{short.name = }{`t44'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta45'}{
							 \describe{
							\item{hyperid = }{`24045'}
							\item{name = }{`theta45'}
							\item{short.name = }{`t45'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta46'}{
							 \describe{
							\item{hyperid = }{`24046'}
							\item{name = }{`theta46'}
							\item{short.name = }{`t46'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta47'}{
							 \describe{
							\item{hyperid = }{`24047'}
							\item{name = }{`theta47'}
							\item{short.name = }{`t47'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta48'}{
							 \describe{
							\item{hyperid = }{`24048'}
							\item{name = }{`theta48'}
							\item{short.name = }{`t48'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta49'}{
							 \describe{
							\item{hyperid = }{`24049'}
							\item{name = }{`theta49'}
							\item{short.name = }{`t49'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta50'}{
							 \describe{
							\item{hyperid = }{`24050'}
							\item{name = }{`theta50'}
							\item{short.name = }{`t50'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta51'}{
							 \describe{
							\item{hyperid = }{`24051'}
							\item{name = }{`theta51'}
							\item{short.name = }{`t51'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta52'}{
							 \describe{
							\item{hyperid = }{`24052'}
							\item{name = }{`theta52'}
							\item{short.name = }{`t52'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta53'}{
							 \describe{
							\item{hyperid = }{`24053'}
							\item{name = }{`theta53'}
							\item{short.name = }{`t53'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta54'}{
							 \describe{
							\item{hyperid = }{`24054'}
							\item{name = }{`theta54'}
							\item{short.name = }{`t54'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta55'}{
							 \describe{
							\item{hyperid = }{`24055'}
							\item{name = }{`theta55'}
							\item{short.name = }{`t55'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta56'}{
							 \describe{
							\item{hyperid = }{`24056'}
							\item{name = }{`theta56'}
							\item{short.name = }{`t56'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta57'}{
							 \describe{
							\item{hyperid = }{`24057'}
							\item{name = }{`theta57'}
							\item{short.name = }{`t57'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta58'}{
							 \describe{
							\item{hyperid = }{`24058'}
							\item{name = }{`theta58'}
							\item{short.name = }{`t58'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta59'}{
							 \describe{
							\item{hyperid = }{`24059'}
							\item{name = }{`theta59'}
							\item{short.name = }{`t59'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta60'}{
							 \describe{
							\item{hyperid = }{`24060'}
							\item{name = }{`theta60'}
							\item{short.name = }{`t60'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta61'}{
							 \describe{
							\item{hyperid = }{`24061'}
							\item{name = }{`theta61'}
							\item{short.name = }{`t61'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta62'}{
							 \describe{
							\item{hyperid = }{`24062'}
							\item{name = }{`theta62'}
							\item{short.name = }{`t62'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta63'}{
							 \describe{
							\item{hyperid = }{`24063'}
							\item{name = }{`theta63'}
							\item{short.name = }{`t63'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta64'}{
							 \describe{
							\item{hyperid = }{`24064'}
							\item{name = }{`theta64'}
							\item{short.name = }{`t64'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta65'}{
							 \describe{
							\item{hyperid = }{`24065'}
							\item{name = }{`theta65'}
							\item{short.name = }{`t65'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta66'}{
							 \describe{
							\item{hyperid = }{`24066'}
							\item{name = }{`theta66'}
							\item{short.name = }{`t66'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta67'}{
							 \describe{
							\item{hyperid = }{`24067'}
							\item{name = }{`theta67'}
							\item{short.name = }{`t67'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta68'}{
							 \describe{
							\item{hyperid = }{`24068'}
							\item{name = }{`theta68'}
							\item{short.name = }{`t68'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta69'}{
							 \describe{
							\item{hyperid = }{`24069'}
							\item{name = }{`theta69'}
							\item{short.name = }{`t69'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta70'}{
							 \describe{
							\item{hyperid = }{`24070'}
							\item{name = }{`theta70'}
							\item{short.name = }{`t70'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta71'}{
							 \describe{
							\item{hyperid = }{`24071'}
							\item{name = }{`theta71'}
							\item{short.name = }{`t71'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta72'}{
							 \describe{
							\item{hyperid = }{`24072'}
							\item{name = }{`theta72'}
							\item{short.name = }{`t72'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta73'}{
							 \describe{
							\item{hyperid = }{`24073'}
							\item{name = }{`theta73'}
							\item{short.name = }{`t73'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta74'}{
							 \describe{
							\item{hyperid = }{`24074'}
							\item{name = }{`theta74'}
							\item{short.name = }{`t74'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta75'}{
							 \describe{
							\item{hyperid = }{`24075'}
							\item{name = }{`theta75'}
							\item{short.name = }{`t75'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta76'}{
							 \describe{
							\item{hyperid = }{`24076'}
							\item{name = }{`theta76'}
							\item{short.name = }{`t76'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta77'}{
							 \describe{
							\item{hyperid = }{`24077'}
							\item{name = }{`theta77'}
							\item{short.name = }{`t77'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta78'}{
							 \describe{
							\item{hyperid = }{`24078'}
							\item{name = }{`theta78'}
							\item{short.name = }{`t78'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta79'}{
							 \describe{
							\item{hyperid = }{`24079'}
							\item{name = }{`theta79'}
							\item{short.name = }{`t79'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta80'}{
							 \describe{
							\item{hyperid = }{`24080'}
							\item{name = }{`theta80'}
							\item{short.name = }{`t80'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta81'}{
							 \describe{
							\item{hyperid = }{`24081'}
							\item{name = }{`theta81'}
							\item{short.name = }{`t81'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta82'}{
							 \describe{
							\item{hyperid = }{`24082'}
							\item{name = }{`theta82'}
							\item{short.name = }{`t82'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta83'}{
							 \describe{
							\item{hyperid = }{`24083'}
							\item{name = }{`theta83'}
							\item{short.name = }{`t83'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta84'}{
							 \describe{
							\item{hyperid = }{`24084'}
							\item{name = }{`theta84'}
							\item{short.name = }{`t84'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta85'}{
							 \describe{
							\item{hyperid = }{`24085'}
							\item{name = }{`theta85'}
							\item{short.name = }{`t85'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta86'}{
							 \describe{
							\item{hyperid = }{`24086'}
							\item{name = }{`theta86'}
							\item{short.name = }{`t86'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta87'}{
							 \describe{
							\item{hyperid = }{`24087'}
							\item{name = }{`theta87'}
							\item{short.name = }{`t87'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta88'}{
							 \describe{
							\item{hyperid = }{`24088'}
							\item{name = }{`theta88'}
							\item{short.name = }{`t88'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta89'}{
							 \describe{
							\item{hyperid = }{`24089'}
							\item{name = }{`theta89'}
							\item{short.name = }{`t89'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta90'}{
							 \describe{
							\item{hyperid = }{`24090'}
							\item{name = }{`theta90'}
							\item{short.name = }{`t90'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta91'}{
							 \describe{
							\item{hyperid = }{`24091'}
							\item{name = }{`theta91'}
							\item{short.name = }{`t91'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta92'}{
							 \describe{
							\item{hyperid = }{`24092'}
							\item{name = }{`theta92'}
							\item{short.name = }{`t92'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta93'}{
							 \describe{
							\item{hyperid = }{`24093'}
							\item{name = }{`theta93'}
							\item{short.name = }{`t93'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta94'}{
							 \describe{
							\item{hyperid = }{`24094'}
							\item{name = }{`theta94'}
							\item{short.name = }{`t94'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta95'}{
							 \describe{
							\item{hyperid = }{`24095'}
							\item{name = }{`theta95'}
							\item{short.name = }{`t95'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta96'}{
							 \describe{
							\item{hyperid = }{`24096'}
							\item{name = }{`theta96'}
							\item{short.name = }{`t96'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta97'}{
							 \describe{
							\item{hyperid = }{`24097'}
							\item{name = }{`theta97'}
							\item{short.name = }{`t97'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta98'}{
							 \describe{
							\item{hyperid = }{`24098'}
							\item{name = }{`theta98'}
							\item{short.name = }{`t98'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta99'}{
							 \describe{
							\item{hyperid = }{`24099'}
							\item{name = }{`theta99'}
							\item{short.name = }{`t99'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta100'}{
							 \describe{
							\item{hyperid = }{`24100'}
							\item{name = }{`theta100'}
							\item{short.name = }{`t100'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{doc = }{`A SPDE3 model'}
								\item{constr = }{`FALSE'}
								\item{nrow.ncol = }{`FALSE'}
								\item{augmented = }{`FALSE'}
								\item{aug.factor = }{`1'}
								\item{aug.constr = }{`NULL'}
								\item{n.div.by = }{`NULL'}
								\item{n.required = }{`TRUE'}
								\item{set.default.values = }{`TRUE'}
								\item{pdf = }{`spde3'}
							 }
						 }
					 }
				 }
				\item{Model `iid1d'.}{Number of hyperparmeters are 1.
					 \describe{
						\item{Hyperparameter `theta'}{
							 \describe{
							\item{hyperid = }{`25001'}
							\item{name = }{`precision'}
							\item{short.name = }{`prec'}
							\item{initial = }{`4'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`wishart1d'}
							\item{param = }{`2 1e-04'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{doc = }{`Gaussian random effect in dim=1 with Wishart prior'}
								\item{constr = }{`FALSE'}
								\item{nrow.ncol = }{`FALSE'}
								\item{augmented = }{`FALSE'}
								\item{aug.factor = }{`1'}
								\item{aug.constr = }{`NULL'}
								\item{n.div.by = }{`NULL'}
								\item{n.required = }{`FALSE'}
								\item{set.default.values = }{`TRUE'}
								\item{pdf = }{`iid123d'}
							 }
						 }
					 }
				 }
				\item{Model `iid2d'.}{Number of hyperparmeters are 3.
					 \describe{
						\item{Hyperparameter `theta1'}{
							 \describe{
							\item{hyperid = }{`26001'}
							\item{name = }{`log precision1'}
							\item{short.name = }{`prec1'}
							\item{initial = }{`4'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`wishart2d'}
							\item{param = }{`4 1 1 0'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
						\item{Hyperparameter `theta2'}{
							 \describe{
							\item{hyperid = }{`26002'}
							\item{name = }{`log precision2'}
							\item{short.name = }{`prec2'}
							\item{initial = }{`4'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
						\item{Hyperparameter `theta3'}{
							 \describe{
							\item{hyperid = }{`26003'}
							\item{name = }{`logit correlation'}
							\item{short.name = }{`cor'}
							\item{initial = }{`4'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) log((1+x)/(1-x))}'}
							\item{from.theta = }{`\code{function(x) 2*exp(x)/(1+exp(x))-1}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{doc = }{`Gaussian random effect in dim=2 with Wishart prior'}
								\item{constr = }{`FALSE'}
								\item{nrow.ncol = }{`FALSE'}
								\item{augmented = }{`TRUE'}
								\item{aug.factor = }{`1'}
								\item{aug.constr = }{`1 2'}
								\item{n.div.by = }{`2'}
								\item{n.required = }{`TRUE'}
								\item{set.default.values = }{`TRUE'}
								\item{pdf = }{`iid123d'}
							 }
						 }
					 }
				 }
				\item{Model `iid3d'.}{Number of hyperparmeters are 6.
					 \describe{
						\item{Hyperparameter `theta1'}{
							 \describe{
							\item{hyperid = }{`27001'}
							\item{name = }{`log precision1'}
							\item{short.name = }{`prec1'}
							\item{initial = }{`4'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`wishart3d'}
							\item{param = }{`7 1 1 1 0 0 0'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
						\item{Hyperparameter `theta2'}{
							 \describe{
							\item{hyperid = }{`27002'}
							\item{name = }{`log precision2'}
							\item{short.name = }{`prec2'}
							\item{initial = }{`4'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
						\item{Hyperparameter `theta3'}{
							 \describe{
							\item{hyperid = }{`27003'}
							\item{name = }{`log precision3'}
							\item{short.name = }{`prec3'}
							\item{initial = }{`4'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
						\item{Hyperparameter `theta4'}{
							 \describe{
							\item{hyperid = }{`27004'}
							\item{name = }{`logit correlation12'}
							\item{short.name = }{`cor12'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) log((1+x)/(1-x))}'}
							\item{from.theta = }{`\code{function(x) 2*exp(x)/(1+exp(x))-1}'}
							}
						}
						\item{Hyperparameter `theta5'}{
							 \describe{
							\item{hyperid = }{`27005'}
							\item{name = }{`logit correlation13'}
							\item{short.name = }{`cor13'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) log((1+x)/(1-x))}'}
							\item{from.theta = }{`\code{function(x) 2*exp(x)/(1+exp(x))-1}'}
							}
						}
						\item{Hyperparameter `theta6'}{
							 \describe{
							\item{hyperid = }{`27006'}
							\item{name = }{`logit correlation23'}
							\item{short.name = }{`cor23'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) log((1+x)/(1-x))}'}
							\item{from.theta = }{`\code{function(x) 2*exp(x)/(1+exp(x))-1}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{doc = }{`Gaussian random effect in dim=3 with Wishart prior'}
								\item{constr = }{`FALSE'}
								\item{nrow.ncol = }{`FALSE'}
								\item{augmented = }{`TRUE'}
								\item{aug.factor = }{`1'}
								\item{aug.constr = }{`1 2 3'}
								\item{n.div.by = }{`3'}
								\item{n.required = }{`TRUE'}
								\item{set.default.values = }{`TRUE'}
								\item{pdf = }{`iid123d'}
							 }
						 }
					 }
				 }
				\item{Model `iid4d'.}{Number of hyperparmeters are 10.
					 \describe{
						\item{Hyperparameter `theta1'}{
							 \describe{
							\item{hyperid = }{`28001'}
							\item{name = }{`log precision1'}
							\item{short.name = }{`prec1'}
							\item{initial = }{`4'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`wishart4d'}
							\item{param = }{`11 1 1 1 1 0 0 0 0 0 0'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
						\item{Hyperparameter `theta2'}{
							 \describe{
							\item{hyperid = }{`28002'}
							\item{name = }{`log precision2'}
							\item{short.name = }{`prec2'}
							\item{initial = }{`4'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
						\item{Hyperparameter `theta3'}{
							 \describe{
							\item{hyperid = }{`28003'}
							\item{name = }{`log precision3'}
							\item{short.name = }{`prec3'}
							\item{initial = }{`4'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
						\item{Hyperparameter `theta4'}{
							 \describe{
							\item{hyperid = }{`28004'}
							\item{name = }{`log precision4'}
							\item{short.name = }{`prec4'}
							\item{initial = }{`4'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
						\item{Hyperparameter `theta5'}{
							 \describe{
							\item{hyperid = }{`28005'}
							\item{name = }{`logit correlation12'}
							\item{short.name = }{`cor12'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) log((1+x)/(1-x))}'}
							\item{from.theta = }{`\code{function(x) 2*exp(x)/(1+exp(x))-1}'}
							}
						}
						\item{Hyperparameter `theta6'}{
							 \describe{
							\item{hyperid = }{`28006'}
							\item{name = }{`logit correlation13'}
							\item{short.name = }{`cor13'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) log((1+x)/(1-x))}'}
							\item{from.theta = }{`\code{function(x) 2*exp(x)/(1+exp(x))-1}'}
							}
						}
						\item{Hyperparameter `theta7'}{
							 \describe{
							\item{hyperid = }{`28007'}
							\item{name = }{`logit correlation14'}
							\item{short.name = }{`cor14'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) log((1+x)/(1-x))}'}
							\item{from.theta = }{`\code{function(x) 2*exp(x)/(1+exp(x))-1}'}
							}
						}
						\item{Hyperparameter `theta8'}{
							 \describe{
							\item{hyperid = }{`28008'}
							\item{name = }{`logit correlation23'}
							\item{short.name = }{`cor23'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) log((1+x)/(1-x))}'}
							\item{from.theta = }{`\code{function(x) 2*exp(x)/(1+exp(x))-1}'}
							}
						}
						\item{Hyperparameter `theta9'}{
							 \describe{
							\item{hyperid = }{`28009'}
							\item{name = }{`logit correlation24'}
							\item{short.name = }{`cor24'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) log((1+x)/(1-x))}'}
							\item{from.theta = }{`\code{function(x) 2*exp(x)/(1+exp(x))-1}'}
							}
						}
						\item{Hyperparameter `theta10'}{
							 \describe{
							\item{hyperid = }{`28010'}
							\item{name = }{`logit correlation34'}
							\item{short.name = }{`cor34'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) log((1+x)/(1-x))}'}
							\item{from.theta = }{`\code{function(x) 2*exp(x)/(1+exp(x))-1}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{doc = }{`Gaussian random effect in dim=4 with Wishart prior'}
								\item{constr = }{`FALSE'}
								\item{nrow.ncol = }{`FALSE'}
								\item{augmented = }{`TRUE'}
								\item{aug.factor = }{`1'}
								\item{aug.constr = }{`1 2 3 4'}
								\item{n.div.by = }{`4'}
								\item{n.required = }{`TRUE'}
								\item{set.default.values = }{`TRUE'}
								\item{pdf = }{`iid123d'}
							 }
						 }
					 }
				 }
				\item{Model `iid5d'.}{Number of hyperparmeters are 15.
					 \describe{
						\item{Hyperparameter `theta1'}{
							 \describe{
							\item{hyperid = }{`29001'}
							\item{name = }{`log precision1'}
							\item{short.name = }{`prec1'}
							\item{initial = }{`4'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`wishart5d'}
							\item{param = }{`16 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
						\item{Hyperparameter `theta2'}{
							 \describe{
							\item{hyperid = }{`29002'}
							\item{name = }{`log precision2'}
							\item{short.name = }{`prec2'}
							\item{initial = }{`4'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
						\item{Hyperparameter `theta3'}{
							 \describe{
							\item{hyperid = }{`29003'}
							\item{name = }{`log precision3'}
							\item{short.name = }{`prec3'}
							\item{initial = }{`4'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
						\item{Hyperparameter `theta4'}{
							 \describe{
							\item{hyperid = }{`29004'}
							\item{name = }{`log precision4'}
							\item{short.name = }{`prec4'}
							\item{initial = }{`4'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
						\item{Hyperparameter `theta5'}{
							 \describe{
							\item{hyperid = }{`29005'}
							\item{name = }{`log precision5'}
							\item{short.name = }{`prec5'}
							\item{initial = }{`4'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
						\item{Hyperparameter `theta6'}{
							 \describe{
							\item{hyperid = }{`29006'}
							\item{name = }{`logit correlation12'}
							\item{short.name = }{`cor12'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) log((1+x)/(1-x))}'}
							\item{from.theta = }{`\code{function(x) 2*exp(x)/(1+exp(x))-1}'}
							}
						}
						\item{Hyperparameter `theta7'}{
							 \describe{
							\item{hyperid = }{`29007'}
							\item{name = }{`logit correlation13'}
							\item{short.name = }{`cor13'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) log((1+x)/(1-x))}'}
							\item{from.theta = }{`\code{function(x) 2*exp(x)/(1+exp(x))-1}'}
							}
						}
						\item{Hyperparameter `theta8'}{
							 \describe{
							\item{hyperid = }{`29008'}
							\item{name = }{`logit correlation14'}
							\item{short.name = }{`cor14'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) log((1+x)/(1-x))}'}
							\item{from.theta = }{`\code{function(x) 2*exp(x)/(1+exp(x))-1}'}
							}
						}
						\item{Hyperparameter `theta9'}{
							 \describe{
							\item{hyperid = }{`29009'}
							\item{name = }{`logit correlation15'}
							\item{short.name = }{`cor15'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) log((1+x)/(1-x))}'}
							\item{from.theta = }{`\code{function(x) 2*exp(x)/(1+exp(x))-1}'}
							}
						}
						\item{Hyperparameter `theta10'}{
							 \describe{
							\item{hyperid = }{`29010'}
							\item{name = }{`logit correlation23'}
							\item{short.name = }{`cor23'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) log((1+x)/(1-x))}'}
							\item{from.theta = }{`\code{function(x) 2*exp(x)/(1+exp(x))-1}'}
							}
						}
						\item{Hyperparameter `theta11'}{
							 \describe{
							\item{hyperid = }{`29011'}
							\item{name = }{`logit correlation24'}
							\item{short.name = }{`cor24'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) log((1+x)/(1-x))}'}
							\item{from.theta = }{`\code{function(x) 2*exp(x)/(1+exp(x))-1}'}
							}
						}
						\item{Hyperparameter `theta12'}{
							 \describe{
							\item{hyperid = }{`29012'}
							\item{name = }{`logit correlation25'}
							\item{short.name = }{`cor25'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) log((1+x)/(1-x))}'}
							\item{from.theta = }{`\code{function(x) 2*exp(x)/(1+exp(x))-1}'}
							}
						}
						\item{Hyperparameter `theta13'}{
							 \describe{
							\item{hyperid = }{`29013'}
							\item{name = }{`logit correlation34'}
							\item{short.name = }{`cor34'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) log((1+x)/(1-x))}'}
							\item{from.theta = }{`\code{function(x) 2*exp(x)/(1+exp(x))-1}'}
							}
						}
						\item{Hyperparameter `theta14'}{
							 \describe{
							\item{hyperid = }{`29014'}
							\item{name = }{`logit correlation35'}
							\item{short.name = }{`cor35'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) log((1+x)/(1-x))}'}
							\item{from.theta = }{`\code{function(x) 2*exp(x)/(1+exp(x))-1}'}
							}
						}
						\item{Hyperparameter `theta15'}{
							 \describe{
							\item{hyperid = }{`29015'}
							\item{name = }{`logit correlation45'}
							\item{short.name = }{`cor45'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) log((1+x)/(1-x))}'}
							\item{from.theta = }{`\code{function(x) 2*exp(x)/(1+exp(x))-1}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{doc = }{`Gaussian random effect in dim=5 with Wishart prior'}
								\item{constr = }{`FALSE'}
								\item{nrow.ncol = }{`FALSE'}
								\item{augmented = }{`TRUE'}
								\item{aug.factor = }{`1'}
								\item{aug.constr = }{`1 2 3 4 5'}
								\item{n.div.by = }{`5'}
								\item{n.required = }{`TRUE'}
								\item{set.default.values = }{`TRUE'}
								\item{pdf = }{`iid123d'}
							 }
						 }
					 }
				 }
				\item{Model `2diid'.}{Number of hyperparmeters are 3.
					 \describe{
						\item{Hyperparameter `theta1'}{
							 \describe{
							\item{hyperid = }{`30001'}
							\item{name = }{`log precision1'}
							\item{short.name = }{`prec1'}
							\item{initial = }{`4'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`loggamma'}
							\item{param = }{`1 5e-05'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
						\item{Hyperparameter `theta2'}{
							 \describe{
							\item{hyperid = }{`30002'}
							\item{name = }{`log precision2'}
							\item{short.name = }{`prec2'}
							\item{initial = }{`4'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`loggamma'}
							\item{param = }{`1 5e-05'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
						\item{Hyperparameter `theta3'}{
							 \describe{
							\item{hyperid = }{`30003'}
							\item{name = }{`correlation'}
							\item{short.name = }{`cor'}
							\item{initial = }{`4'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`normal'}
							\item{param = }{`0 0.15'}
							\item{to.theta = }{`\code{function(x) log((1+x)/(1-x))}'}
							\item{from.theta = }{`\code{function(x) 2*exp(x)/(1+exp(x))-1}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{doc = }{`(This model is obsolute)'}
								\item{constr = }{`FALSE'}
								\item{nrow.ncol = }{`FALSE'}
								\item{augmented = }{`FALSE'}
								\item{aug.factor = }{`1'}
								\item{aug.constr = }{`1 2'}
								\item{n.div.by = }{`2'}
								\item{n.required = }{`TRUE'}
								\item{set.default.values = }{`TRUE'}
								\item{pdf = }{`iid123d'}
							 }
						 }
					 }
				 }
				\item{Model `z'.}{Number of hyperparmeters are 1.
					 \describe{
						\item{Hyperparameter `theta'}{
							 \describe{
							\item{hyperid = }{`31001'}
							\item{name = }{`log precision'}
							\item{short.name = }{`prec'}
							\item{initial = }{`4'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`loggamma'}
							\item{param = }{`1 5e-05'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{doc = }{`The z-model in a classical mixed model formulation'}
								\item{constr = }{`FALSE'}
								\item{nrow.ncol = }{`FALSE'}
								\item{augmented = }{`FALSE'}
								\item{aug.factor = }{`1'}
								\item{aug.constr = }{`NULL'}
								\item{n.div.by = }{`NULL'}
								\item{n.required = }{`TRUE'}
								\item{set.default.values = }{`TRUE'}
								\item{pdf = }{`z'}
								\item{status = }{`experimental'}
							 }
						 }
					 }
				 }
				\item{Model `rw2d'.}{Number of hyperparmeters are 1.
					 \describe{
						\item{Hyperparameter `theta'}{
							 \describe{
							\item{hyperid = }{`32001'}
							\item{name = }{`log precision'}
							\item{short.name = }{`prec'}
							\item{initial = }{`4'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`loggamma'}
							\item{param = }{`1 5e-05'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{doc = }{`Thin-plate spline model'}
								\item{constr = }{`TRUE'}
								\item{nrow.ncol = }{`TRUE'}
								\item{augmented = }{`FALSE'}
								\item{aug.factor = }{`1'}
								\item{aug.constr = }{`NULL'}
								\item{n.div.by = }{`NULL'}
								\item{n.required = }{`FALSE'}
								\item{set.default.values = }{`TRUE'}
								\item{pdf = }{`rw2d'}
							 }
						 }
					 }
				 }
				\item{Model `rw2diid'.}{Number of hyperparmeters are 2.
					 \describe{
						\item{Hyperparameter `theta1'}{
							 \describe{
							\item{hyperid = }{`33001'}
							\item{name = }{`log precision'}
							\item{short.name = }{`prec'}
							\item{prior = }{`pc.prec'}
							\item{param = }{`1 0.01'}
							\item{initial = }{`4'}
							\item{fixed = }{`FALSE'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
						\item{Hyperparameter `theta2'}{
							 \describe{
							\item{hyperid = }{`33002'}
							\item{name = }{`logit phi'}
							\item{short.name = }{`phi'}
							\item{prior = }{`pc'}
							\item{param = }{`0.5 0.5'}
							\item{initial = }{`3'}
							\item{fixed = }{`FALSE'}
							\item{to.theta = }{`\code{function(x) log(x/(1-x))}'}
							\item{from.theta = }{`\code{function(x) exp(x)/(1+exp(x))}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{doc = }{`Thin-plate spline with iid noise'}
								\item{constr = }{`TRUE'}
								\item{nrow.ncol = }{`TRUE'}
								\item{augmented = }{`TRUE'}
								\item{aug.factor = }{`2'}
								\item{aug.constr = }{`2'}
								\item{n.div.by = }{`NULL'}
								\item{n.required = }{`FALSE'}
								\item{set.default.values = }{`TRUE'}
								\item{status = }{`experimental'}
								\item{pdf = }{`rw2diid'}
							 }
						 }
					 }
				 }
				\item{Model `slm'.}{Number of hyperparmeters are 2.
					 \describe{
						\item{Hyperparameter `theta1'}{
							 \describe{
							\item{hyperid = }{`34001'}
							\item{name = }{`log precision'}
							\item{short.name = }{`prec'}
							\item{initial = }{`4'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`loggamma'}
							\item{param = }{`1 5e-05'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
						\item{Hyperparameter `theta2'}{
							 \describe{
							\item{hyperid = }{`34002'}
							\item{name = }{`rho'}
							\item{short.name = }{`rho'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`normal'}
							\item{param = }{`0 10'}
							\item{to.theta = }{`\code{function(x) log(x/(1-x))}'}
							\item{from.theta = }{`\code{function(x) 1/(1+exp(-x))}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{doc = }{`Spatial lag model'}
								\item{constr = }{`FALSE'}
								\item{nrow.ncol = }{`FALSE'}
								\item{augmented = }{`FALSE'}
								\item{aug.factor = }{`1'}
								\item{aug.constr = }{`NULL'}
								\item{n.div.by = }{`NULL'}
								\item{n.required = }{`TRUE'}
								\item{set.default.values = }{`TRUE'}
								\item{pdf = }{`slm'}
								\item{status = }{`experimental'}
							 }
						 }
					 }
				 }
				\item{Model `matern2d'.}{Number of hyperparmeters are 2.
					 \describe{
						\item{Hyperparameter `theta1'}{
							 \describe{
							\item{hyperid = }{`35001'}
							\item{name = }{`log precision'}
							\item{short.name = }{`prec'}
							\item{initial = }{`4'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`loggamma'}
							\item{param = }{`1 5e-05'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
						\item{Hyperparameter `theta2'}{
							 \describe{
							\item{hyperid = }{`35002'}
							\item{name = }{`log range'}
							\item{short.name = }{`range'}
							\item{initial = }{`2'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`loggamma'}
							\item{param = }{`1 0.01'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{doc = }{`Matern covariance function on a regular grid'}
								\item{constr = }{`FALSE'}
								\item{nrow.ncol = }{`TRUE'}
								\item{augmented = }{`FALSE'}
								\item{aug.factor = }{`1'}
								\item{aug.constr = }{`NULL'}
								\item{n.div.by = }{`NULL'}
								\item{n.required = }{`FALSE'}
								\item{set.default.values = }{`TRUE'}
								\item{pdf = }{`matern2d'}
							 }
						 }
					 }
				 }
				\item{Model `dmatern'.}{Number of hyperparmeters are 3.
					 \describe{
						\item{Hyperparameter `theta1'}{
							 \describe{
							\item{hyperid = }{`35101'}
							\item{name = }{`log precision'}
							\item{short.name = }{`prec'}
							\item{initial = }{`3'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`pc.prec'}
							\item{param = }{`1 0.01'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
						\item{Hyperparameter `theta2'}{
							 \describe{
							\item{hyperid = }{`35102'}
							\item{name = }{`log range'}
							\item{short.name = }{`range'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`pc.range'}
							\item{param = }{`1 0.5'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
						\item{Hyperparameter `theta3'}{
							 \describe{
							\item{hyperid = }{`35103'}
							\item{name = }{`log nu'}
							\item{short.name = }{`nu'}
							\item{initial = }{`-0.693147180559945'}
							\item{fixed = }{`TRUE'}
							\item{prior = }{`loggamma'}
							\item{param = }{`0.5 1'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{doc = }{`Dense Matern field'}
								\item{constr = }{`FALSE'}
								\item{nrow.ncol = }{`FALSE'}
								\item{augmented = }{`FALSE'}
								\item{aug.factor = }{`1'}
								\item{aug.constr = }{`NULL'}
								\item{n.div.by = }{`NULL'}
								\item{n.required = }{`TRUE'}
								\item{set.default.values = }{`TRUE'}
								\item{status = }{`experimental'}
								\item{pdf = }{`dmatern'}
							 }
						 }
					 }
				 }
				\item{Model `copy'.}{Number of hyperparmeters are 1.
					 \describe{
						\item{Hyperparameter `theta'}{
							 \describe{
							\item{hyperid = }{`36001'}
							\item{name = }{`beta'}
							\item{short.name = }{`b'}
							\item{initial = }{`1'}
							\item{fixed = }{`TRUE'}
							\item{prior = }{`normal'}
							\item{param = }{`1 10'}
							\item{to.theta = }{`\code{function(x, REPLACE.ME.low, REPLACE.ME.high) {} \code{                             if (all(is.infinite(c(low, high))) || low == high) {} \code{                                 return (x)} \code{                             } else if (all(is.finite(c(low, high)))) {} \code{                                 stopifnot(low < high)} \code{                                 return (log( - (low - x)/(high -x)))} \code{                             } else if (is.finite(low) && is.infinite(high) && high > low) {} \code{                                 return (log(x-low))} \code{                             } else {} \code{                                 stop("Condition not yet implemented")} \code{                             }} \code{                         }}'}
							\item{from.theta = }{`\code{function(x, REPLACE.ME.low, REPLACE.ME.high) {} \code{                             if (all(is.infinite(c(low, high))) || low == high) {} \code{                                 return (x)} \code{                             } else if (all(is.finite(c(low, high)))) {} \code{                                 stopifnot(low < high)} \code{                                 return (low + exp(x)/(1+exp(x)) * (high - low))} \code{                             } else if (is.finite(low) && is.infinite(high) && high > low) {} \code{                                 return (low + exp(x))} \code{                             } else {} \code{                                 stop("Condition not yet implemented")} \code{                             }} \code{                         }}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{doc = }{`Create a copy of a model component'}
								\item{constr = }{`FALSE'}
								\item{nrow.ncol = }{`FALSE'}
								\item{augmented = }{`FALSE'}
								\item{aug.factor = }{`1'}
								\item{aug.constr = }{`NULL'}
								\item{n.div.by = }{`NULL'}
								\item{n.required = }{`FALSE'}
								\item{set.default.values = }{`FALSE'}
								\item{pdf = }{`NA'}
							 }
						 }
					 }
				 }
				\item{Model `clinear'.}{Number of hyperparmeters are 1.
					 \describe{
						\item{Hyperparameter `theta'}{
							 \describe{
							\item{hyperid = }{`37001'}
							\item{name = }{`beta'}
							\item{short.name = }{`b'}
							\item{initial = }{`1'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`normal'}
							\item{param = }{`1 10'}
							\item{to.theta = }{`\code{function(x, REPLACE.ME.low, REPLACE.ME.high) {} \code{                             if (all(is.infinite(c(low, high))) || low == high) {} \code{                                 stopifnot(low < high)} \code{                                 return (x)} \code{                             } else if (all(is.finite(c(low, high)))) {} \code{                                 stopifnot(low < high)} \code{                                 return (log( - (low - x)/(high -x)))} \code{                             } else if (is.finite(low) && is.infinite(high) && high > low) {} \code{                                 return (log(x-low))} \code{                             } else {} \code{                                 stop("Condition not yet implemented")} \code{                             }} \code{                         }}'}
							\item{from.theta = }{`\code{function(x, REPLACE.ME.low, REPLACE.ME.high) {} \code{                             if (all(is.infinite(c(low, high))) || low == high) {} \code{                                 stopifnot(low < high)} \code{                                 return (x)} \code{                             } else if (all(is.finite(c(low, high)))) {} \code{                                 stopifnot(low < high)} \code{                                 return (low + exp(x)/(1+exp(x)) * (high - low))} \code{                             } else if (is.finite(low) && is.infinite(high) && high > low) {} \code{                                 return (low + exp(x))} \code{                             } else {} \code{                                 stop("Condition not yet implemented")} \code{                             }} \code{                         }}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{doc = }{`Constrained linear effect'}
								\item{constr = }{`FALSE'}
								\item{nrow.ncol = }{`FALSE'}
								\item{augmented = }{`FALSE'}
								\item{aug.factor = }{`1'}
								\item{aug.constr = }{`NULL'}
								\item{n.div.by = }{`NULL'}
								\item{n.required = }{`FALSE'}
								\item{set.default.values = }{`FALSE'}
								\item{pdf = }{`clinear'}
							 }
						 }
					 }
				 }
				\item{Model `sigm'.}{Number of hyperparmeters are 3.
					 \describe{
						\item{Hyperparameter `theta1'}{
							 \describe{
							\item{hyperid = }{`38001'}
							\item{name = }{`beta'}
							\item{short.name = }{`b'}
							\item{initial = }{`1'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`normal'}
							\item{param = }{`1 10'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta2'}{
							 \describe{
							\item{hyperid = }{`38002'}
							\item{name = }{`loghalflife'}
							\item{short.name = }{`halflife'}
							\item{initial = }{`3'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`loggamma'}
							\item{param = }{`3 1'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
						\item{Hyperparameter `theta3'}{
							 \describe{
							\item{hyperid = }{`38003'}
							\item{name = }{`logshape'}
							\item{short.name = }{`shape'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`loggamma'}
							\item{param = }{`10 10'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{doc = }{`Sigmoidal effect of a covariate'}
								\item{constr = }{`FALSE'}
								\item{nrow.ncol = }{`FALSE'}
								\item{augmented = }{`FALSE'}
								\item{aug.factor = }{`1'}
								\item{aug.constr = }{`NULL'}
								\item{n.div.by = }{`NULL'}
								\item{n.required = }{`FALSE'}
								\item{set.default.values = }{`FALSE'}
								\item{status = }{`experimental'}
								\item{pdf = }{`sigm'}
							 }
						 }
					 }
				 }
				\item{Model `revsigm'.}{Number of hyperparmeters are 3.
					 \describe{
						\item{Hyperparameter `theta1'}{
							 \describe{
							\item{hyperid = }{`39001'}
							\item{name = }{`beta'}
							\item{short.name = }{`b'}
							\item{initial = }{`1'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`normal'}
							\item{param = }{`1 10'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta2'}{
							 \describe{
							\item{hyperid = }{`39002'}
							\item{name = }{`loghalflife'}
							\item{short.name = }{`halflife'}
							\item{initial = }{`3'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`loggamma'}
							\item{param = }{`3 1'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
						\item{Hyperparameter `theta3'}{
							 \describe{
							\item{hyperid = }{`39003'}
							\item{name = }{`logshape'}
							\item{short.name = }{`shape'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`loggamma'}
							\item{param = }{`10 10'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{doc = }{`Reverse sigmoidal effect of a covariate'}
								\item{constr = }{`FALSE'}
								\item{nrow.ncol = }{`FALSE'}
								\item{augmented = }{`FALSE'}
								\item{aug.factor = }{`1'}
								\item{aug.constr = }{`NULL'}
								\item{n.div.by = }{`NULL'}
								\item{n.required = }{`FALSE'}
								\item{set.default.values = }{`FALSE'}
								\item{status = }{`experimental'}
								\item{pdf = }{`sigm'}
							 }
						 }
					 }
				 }
				\item{Model `log1exp'.}{Number of hyperparmeters are 3.
					 \describe{
						\item{Hyperparameter `theta1'}{
							 \describe{
							\item{hyperid = }{`39011'}
							\item{name = }{`beta'}
							\item{short.name = }{`b'}
							\item{initial = }{`1'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`normal'}
							\item{param = }{`0 1'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta2'}{
							 \describe{
							\item{hyperid = }{`39012'}
							\item{name = }{`alpha'}
							\item{short.name = }{`a'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`normal'}
							\item{param = }{`0 1'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta3'}{
							 \describe{
							\item{hyperid = }{`39013'}
							\item{name = }{`gamma'}
							\item{short.name = }{`g'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`normal'}
							\item{param = }{`0 1'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{doc = }{`A nonlinear model of a covariate'}
								\item{constr = }{`FALSE'}
								\item{nrow.ncol = }{`FALSE'}
								\item{augmented = }{`FALSE'}
								\item{aug.factor = }{`1'}
								\item{aug.constr = }{`NULL'}
								\item{n.div.by = }{`NULL'}
								\item{n.required = }{`FALSE'}
								\item{set.default.values = }{`FALSE'}
								\item{status = }{`experimental'}
								\item{pdf = }{`log1exp'}
							 }
						 }
					 }
				 }
				\item{Model `logdist'.}{Number of hyperparmeters are 3.
					 \describe{
						\item{Hyperparameter `theta1'}{
							 \describe{
							\item{hyperid = }{`39021'}
							\item{name = }{`beta'}
							\item{short.name = }{`b'}
							\item{initial = }{`1'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`normal'}
							\item{param = }{`0 1'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta2'}{
							 \describe{
							\item{hyperid = }{`39022'}
							\item{name = }{`alpha1'}
							\item{short.name = }{`a1'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`loggamma'}
							\item{param = }{`0.1 1'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
						\item{Hyperparameter `theta3'}{
							 \describe{
							\item{hyperid = }{`39023'}
							\item{name = }{`alpha2'}
							\item{short.name = }{`a2'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`loggamma'}
							\item{param = }{`0.1 1'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{doc = }{`A nonlinear model of a covariate'}
								\item{constr = }{`FALSE'}
								\item{nrow.ncol = }{`FALSE'}
								\item{augmented = }{`FALSE'}
								\item{aug.factor = }{`1'}
								\item{aug.constr = }{`NULL'}
								\item{n.div.by = }{`NULL'}
								\item{n.required = }{`FALSE'}
								\item{set.default.values = }{`FALSE'}
								\item{status = }{`experimental'}
								\item{pdf = }{`logdist'}
							 }
						 }
					 }
				 }
			 }
		 }
		\item{Section `group'.}{
			 Valid models in this section are:
			 \describe{
				\item{Model `exchangeable'.}{Number of hyperparmeters are 1.
					 \describe{
						\item{Hyperparameter `theta'}{
							 \describe{
							\item{hyperid = }{`40001'}
							\item{name = }{`logit correlation'}
							\item{short.name = }{`rho'}
							\item{initial = }{`1'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`normal'}
							\item{param = }{`0 0.2'}
							\item{to.theta = }{`\code{function(x, REPLACE.ME.ngroup) log((1+x*(ngroup-1))/(1-x))}'}
							\item{from.theta = }{`\code{function(x, REPLACE.ME.ngroup) (exp(x)-1)/(exp(x) + ngroup -1)}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{doc = }{`Exchangeable correlations'}
							 }
						 }
					 }
				 }
				\item{Model `exchangeablepos'.}{Number of hyperparmeters are 1.
					 \describe{
						\item{Hyperparameter `theta'}{
							 \describe{
							\item{hyperid = }{`40101'}
							\item{name = }{`logit correlation'}
							\item{short.name = }{`rho'}
							\item{initial = }{`1'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`pc.cor0'}
							\item{param = }{`0.5 0.5'}
							\item{to.theta = }{`\code{function(x) log(x/(1-x))}'}
							\item{from.theta = }{`\code{function(x) exp(x)/(1+exp(x))}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{doc = }{`Exchangeable positive correlations'}
							 }
						 }
					 }
				 }
				\item{Model `ar1'.}{Number of hyperparmeters are 1.
					 \describe{
						\item{Hyperparameter `theta'}{
							 \describe{
							\item{hyperid = }{`41001'}
							\item{name = }{`logit correlation'}
							\item{short.name = }{`rho'}
							\item{initial = }{`2'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`normal'}
							\item{param = }{`0 0.15'}
							\item{to.theta = }{`\code{function(x) log((1+x)/(1-x))}'}
							\item{from.theta = }{`\code{function(x) 2*exp(x)/(1+exp(x))-1}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{doc = }{`AR(1) correlations'}
							 }
						 }
					 }
				 }
				\item{Model `ar'.}{Number of hyperparmeters are 11.
					 \describe{
						\item{Hyperparameter `theta1'}{
							 \describe{
							\item{hyperid = }{`42001'}
							\item{name = }{`log precision'}
							\item{short.name = }{`prec'}
							\item{initial = }{`0'}
							\item{fixed = }{`TRUE'}
							\item{prior = }{`pc.prec'}
							\item{param = }{`3 0.01'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
						\item{Hyperparameter `theta2'}{
							 \describe{
							\item{hyperid = }{`42002'}
							\item{name = }{`pacf1'}
							\item{short.name = }{`pacf1'}
							\item{initial = }{`2'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`pc.cor0'}
							\item{param = }{`0.5 0.5'}
							\item{to.theta = }{`\code{function(x) log((1+x)/(1-x))}'}
							\item{from.theta = }{`\code{function(x) 2*exp(x)/(1+exp(x))-1}'}
							}
						}
						\item{Hyperparameter `theta3'}{
							 \describe{
							\item{hyperid = }{`42003'}
							\item{name = }{`pacf2'}
							\item{short.name = }{`pacf2'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`pc.cor0'}
							\item{param = }{`0.5 0.4'}
							\item{to.theta = }{`\code{function(x) log((1+x)/(1-x))}'}
							\item{from.theta = }{`\code{function(x) 2*exp(x)/(1+exp(x))-1}'}
							}
						}
						\item{Hyperparameter `theta4'}{
							 \describe{
							\item{hyperid = }{`42004'}
							\item{name = }{`pacf3'}
							\item{short.name = }{`pacf3'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`pc.cor0'}
							\item{param = }{`0.5 0.3'}
							\item{to.theta = }{`\code{function(x) log((1+x)/(1-x))}'}
							\item{from.theta = }{`\code{function(x) 2*exp(x)/(1+exp(x))-1}'}
							}
						}
						\item{Hyperparameter `theta5'}{
							 \describe{
							\item{hyperid = }{`42005'}
							\item{name = }{`pacf4'}
							\item{short.name = }{`pacf4'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`pc.cor0'}
							\item{param = }{`0.5 0.2'}
							\item{to.theta = }{`\code{function(x) log((1+x)/(1-x))}'}
							\item{from.theta = }{`\code{function(x) 2*exp(x)/(1+exp(x))-1}'}
							}
						}
						\item{Hyperparameter `theta6'}{
							 \describe{
							\item{hyperid = }{`42006'}
							\item{name = }{`pacf5'}
							\item{short.name = }{`pacf5'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`pc.cor0'}
							\item{param = }{`0.5 0.1'}
							\item{to.theta = }{`\code{function(x) log((1+x)/(1-x))}'}
							\item{from.theta = }{`\code{function(x) 2*exp(x)/(1+exp(x))-1}'}
							}
						}
						\item{Hyperparameter `theta7'}{
							 \describe{
							\item{hyperid = }{`42007'}
							\item{name = }{`pacf6'}
							\item{short.name = }{`pacf6'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`pc.cor0'}
							\item{param = }{`0.5 0.1'}
							\item{to.theta = }{`\code{function(x) log((1+x)/(1-x))}'}
							\item{from.theta = }{`\code{function(x) 2*exp(x)/(1+exp(x))-1}'}
							}
						}
						\item{Hyperparameter `theta8'}{
							 \describe{
							\item{hyperid = }{`42008'}
							\item{name = }{`pacf7'}
							\item{short.name = }{`pacf7'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`pc.cor0'}
							\item{param = }{`0.5 0.1'}
							\item{to.theta = }{`\code{function(x) log((1+x)/(1-x))}'}
							\item{from.theta = }{`\code{function(x) 2*exp(x)/(1+exp(x))-1}'}
							}
						}
						\item{Hyperparameter `theta9'}{
							 \describe{
							\item{hyperid = }{`42009'}
							\item{name = }{`pacf8'}
							\item{short.name = }{`pacf8'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`pc.cor0'}
							\item{param = }{`0.5 0.1'}
							\item{to.theta = }{`\code{function(x) log((1+x)/(1-x))}'}
							\item{from.theta = }{`\code{function(x) 2*exp(x)/(1+exp(x))-1}'}
							}
						}
						\item{Hyperparameter `theta10'}{
							 \describe{
							\item{hyperid = }{`42010'}
							\item{name = }{`pacf9'}
							\item{short.name = }{`pacf9'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`pc.cor0'}
							\item{param = }{`0.5 0.1'}
							\item{to.theta = }{`\code{function(x) log((1+x)/(1-x))}'}
							\item{from.theta = }{`\code{function(x) 2*exp(x)/(1+exp(x))-1}'}
							}
						}
						\item{Hyperparameter `theta11'}{
							 \describe{
							\item{hyperid = }{`42011'}
							\item{name = }{`pacf10'}
							\item{short.name = }{`pacf10'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`pc.cor0'}
							\item{param = }{`0.5 0.1'}
							\item{to.theta = }{`\code{function(x) log((1+x)/(1-x))}'}
							\item{from.theta = }{`\code{function(x) 2*exp(x)/(1+exp(x))-1}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{doc = }{`AR(p) correlations'}
							 }
						 }
					 }
				 }
				\item{Model `rw1'.}{Number of hyperparmeters are 1.
					 \describe{
						\item{Hyperparameter `theta'}{
							 \describe{
							\item{hyperid = }{`43001'}
							\item{name = }{`log precision'}
							\item{short.name = }{`prec'}
							\item{prior = }{`loggamma'}
							\item{param = }{`1 5e-05'}
							\item{initial = }{`0'}
							\item{fixed = }{`TRUE'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{doc = }{`Random walk of order 1'}
							 }
						 }
					 }
				 }
				\item{Model `rw2'.}{Number of hyperparmeters are 1.
					 \describe{
						\item{Hyperparameter `theta'}{
							 \describe{
							\item{hyperid = }{`44001'}
							\item{name = }{`log precision'}
							\item{short.name = }{`prec'}
							\item{prior = }{`loggamma'}
							\item{param = }{`1 5e-05'}
							\item{initial = }{`0'}
							\item{fixed = }{`TRUE'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{doc = }{`Random walk of order 2'}
							 }
						 }
					 }
				 }
				\item{Model `besag'.}{Number of hyperparmeters are 1.
					 \describe{
						\item{Hyperparameter `theta'}{
							 \describe{
							\item{hyperid = }{`45001'}
							\item{name = }{`log precision'}
							\item{short.name = }{`prec'}
							\item{prior = }{`loggamma'}
							\item{param = }{`1 5e-05'}
							\item{initial = }{`0'}
							\item{fixed = }{`TRUE'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{doc = }{`Besag model'}
							 }
						 }
					 }
				 }
				\item{Model `iid'.}{Number of hyperparmeters are 1.
					 \describe{
						\item{Hyperparameter `theta'}{
							 \describe{
							\item{hyperid = }{`46001'}
							\item{name = }{`log precision'}
							\item{short.name = }{`prec'}
							\item{prior = }{`loggamma'}
							\item{param = }{`1 5e-05'}
							\item{initial = }{`0'}
							\item{fixed = }{`TRUE'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{doc = }{`Independent model'}
							 }
						 }
					 }
				 }
			 }
		 }
		\item{Section `mix'.}{
			 Valid models in this section are:
			 \describe{
				\item{Model `gaussian'.}{Number of hyperparmeters are 1.
					 \describe{
						\item{Hyperparameter `theta'}{
							 \describe{
							\item{hyperid = }{`47001'}
							\item{name = }{`log precision'}
							\item{short.name = }{`prec'}
							\item{prior = }{`pc.prec'}
							\item{param = }{`1 0.01'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{doc = }{`Gaussian mixture'}
							 }
						 }
					 }
				 }
				\item{Model `loggamma'.}{Number of hyperparmeters are 1.
					 \describe{
						\item{Hyperparameter `theta'}{
							 \describe{
							\item{hyperid = }{`47101'}
							\item{name = }{`log precision'}
							\item{short.name = }{`prec'}
							\item{prior = }{`pc.mgamma'}
							\item{param = }{`4.8'}
							\item{initial = }{`4'}
							\item{fixed = }{`FALSE'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{doc = }{`LogGamma mixture'}
							 }
						 }
					 }
				 }
				\item{Model `mloggamma'.}{Number of hyperparmeters are 1.
					 \describe{
						\item{Hyperparameter `theta'}{
							 \describe{
							\item{hyperid = }{`47201'}
							\item{name = }{`log precision'}
							\item{short.name = }{`prec'}
							\item{prior = }{`pc.mgamma'}
							\item{param = }{`4.8'}
							\item{initial = }{`4'}
							\item{fixed = }{`FALSE'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{doc = }{`Minus-LogGamma mixture'}
							 }
						 }
					 }
				 }
			 }
		 }
		\item{Section `link'.}{
			 Valid models in this section are:
			 \describe{
				\item{Model `default'.}{Number of hyperparmeters are 0.
				 }
				\item{Model `cloglog'.}{Number of hyperparmeters are 0.
				 }
				\item{Model `loglog'.}{Number of hyperparmeters are 0.
				 }
				\item{Model `identity'.}{Number of hyperparmeters are 0.
				 }
				\item{Model `inverse'.}{Number of hyperparmeters are 0.
				 }
				\item{Model `log'.}{Number of hyperparmeters are 0.
				 }
				\item{Model `neglog'.}{Number of hyperparmeters are 0.
				 }
				\item{Model `logit'.}{Number of hyperparmeters are 0.
				 }
				\item{Model `probit'.}{Number of hyperparmeters are 0.
				 }
				\item{Model `cauchit'.}{Number of hyperparmeters are 0.
				 }
				\item{Model `tan'.}{Number of hyperparmeters are 0.
				 }
				\item{Model `quantile'.}{Number of hyperparmeters are 0.
				 }
				\item{Model `pquantile'.}{Number of hyperparmeters are 0.
				 }
				\item{Model `sslogit'.}{Number of hyperparmeters are 2.
					 \describe{
						\item{Hyperparameter `theta1'}{
							 \describe{
							\item{hyperid = }{`48001'}
							\item{name = }{`sensitivity'}
							\item{short.name = }{`sens'}
							\item{prior = }{`logitbeta'}
							\item{param = }{`10 5'}
							\item{initial = }{`1'}
							\item{fixed = }{`FALSE'}
							\item{to.theta = }{`\code{function(x) log(x/(1-x))}'}
							\item{from.theta = }{`\code{function(x) exp(x)/(1+exp(x))}'}
							}
						}
						\item{Hyperparameter `theta2'}{
							 \describe{
							\item{hyperid = }{`48002'}
							\item{name = }{`specificity'}
							\item{short.name = }{`spec'}
							\item{prior = }{`logitbeta'}
							\item{param = }{`10 5'}
							\item{initial = }{`1'}
							\item{fixed = }{`FALSE'}
							\item{to.theta = }{`\code{function(x) log(x/(1-x))}'}
							\item{from.theta = }{`\code{function(x) exp(x)/(1+exp(x))}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{doc = }{`Logit link with sensitivity and specificity'}
								\item{pdf = }{`NA'}
							 }
						 }
					 }
				 }
				\item{Model `logoffset'.}{Number of hyperparmeters are 1.
					 \describe{
						\item{Hyperparameter `theta'}{
							 \describe{
							\item{hyperid = }{`49001'}
							\item{name = }{`beta'}
							\item{short.name = }{`b'}
							\item{prior = }{`normal'}
							\item{param = }{`0 100'}
							\item{initial = }{`0'}
							\item{fixed = }{`TRUE'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{doc = }{`Log-link with an offset'}
								\item{pdf = }{`logoffset'}
							 }
						 }
					 }
				 }
				\item{Model `logitoffset'.}{Number of hyperparmeters are 1.
					 \describe{
						\item{Hyperparameter `theta'}{
							 \describe{
							\item{hyperid = }{`49011'}
							\item{name = }{`prob'}
							\item{short.name = }{`p'}
							\item{prior = }{`normal'}
							\item{param = }{`-1 100'}
							\item{initial = }{`-1'}
							\item{fixed = }{`FALSE'}
							\item{to.theta = }{`\code{function(x) log(x/(1-x))}'}
							\item{from.theta = }{`\code{function(x) exp(x)/(1+exp(x))}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{doc = }{`Logit-link with an offset'}
								\item{status = }{`experimental'}
								\item{pdf = }{`logitoffset'}
							 }
						 }
					 }
				 }
				\item{Model `test1'.}{Number of hyperparmeters are 1.
					 \describe{
						\item{Hyperparameter `theta'}{
							 \describe{
							\item{hyperid = }{`50001'}
							\item{name = }{`beta'}
							\item{short.name = }{`b'}
							\item{prior = }{`normal'}
							\item{param = }{`0 100'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{doc = }{`A test1-link function (experimental)'}
								\item{pdf = }{`NA'}
							 }
						 }
					 }
				 }
				\item{Model `special1'.}{Number of hyperparmeters are 11.
					 \describe{
						\item{Hyperparameter `theta1'}{
							 \describe{
							\item{hyperid = }{`51001'}
							\item{name = }{`log precision'}
							\item{short.name = }{`prec'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`loggamma'}
							\item{param = }{`1 1'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta2'}{
							 \describe{
							\item{hyperid = }{`51002'}
							\item{name = }{`beta1'}
							\item{short.name = }{`beta1'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`mvnorm'}
							\item{param = }{`0 100'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta3'}{
							 \describe{
							\item{hyperid = }{`51003'}
							\item{name = }{`beta2'}
							\item{short.name = }{`beta2'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta4'}{
							 \describe{
							\item{hyperid = }{`51004'}
							\item{name = }{`beta3'}
							\item{short.name = }{`beta3'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta5'}{
							 \describe{
							\item{hyperid = }{`51005'}
							\item{name = }{`beta4'}
							\item{short.name = }{`beta4'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta6'}{
							 \describe{
							\item{hyperid = }{`51006'}
							\item{name = }{`beta5'}
							\item{short.name = }{`beta5'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta7'}{
							 \describe{
							\item{hyperid = }{`51007'}
							\item{name = }{`beta6'}
							\item{short.name = }{`beta6'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta8'}{
							 \describe{
							\item{hyperid = }{`51008'}
							\item{name = }{`beta7'}
							\item{short.name = }{`beta7'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta9'}{
							 \describe{
							\item{hyperid = }{`51009'}
							\item{name = }{`beta8'}
							\item{short.name = }{`beta8'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta10'}{
							 \describe{
							\item{hyperid = }{`51010'}
							\item{name = }{`beta9'}
							\item{short.name = }{`beta9'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta11'}{
							 \describe{
							\item{hyperid = }{`51011'}
							\item{name = }{`beta10'}
							\item{short.name = }{`beta10'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{doc = }{`A special1-link function (experimental)'}
								\item{pdf = }{`NA'}
							 }
						 }
					 }
				 }
				\item{Model `special2'.}{Number of hyperparmeters are 1.
					 \describe{
						\item{Hyperparameter `theta'}{
							 \describe{
							\item{hyperid = }{`52001'}
							\item{name = }{`beta'}
							\item{short.name = }{`b'}
							\item{prior = }{`normal'}
							\item{param = }{`0 10'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{doc = }{`A special2-link function (experimental)'}
								\item{pdf = }{`NA'}
							 }
						 }
					 }
				 }
			 }
		 }
		\item{Section `predictor'.}{
			 Valid models in this section are:
			 \describe{
				\item{Model `predictor'.}{Number of hyperparmeters are 1.
					 \describe{
						\item{Hyperparameter `theta'}{
							 \describe{
							\item{hyperid = }{`53001'}
							\item{name = }{`log precision'}
							\item{short.name = }{`prec'}
							\item{initial = }{`12'}
							\item{fixed = }{`TRUE'}
							\item{prior = }{`loggamma'}
							\item{param = }{`1 1e-05'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{doc = }{`(not used)'}
							 }
						 }
					 }
				 }
			 }
		 }
		\item{Section `hazard'.}{
			 Valid models in this section are:
			 \describe{
				\item{Model `rw1'.}{Number of hyperparmeters are 1.
					 \describe{
						\item{Hyperparameter `theta'}{
							 \describe{
							\item{hyperid = }{`54001'}
							\item{name = }{`log precision'}
							\item{short.name = }{`prec'}
							\item{initial = }{`4'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`loggamma'}
							\item{param = }{`1 5e-05'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{doc = }{`A random walk of order 1 for the log-hazard'}
							 }
						 }
					 }
				 }
				\item{Model `rw2'.}{Number of hyperparmeters are 1.
					 \describe{
						\item{Hyperparameter `theta'}{
							 \describe{
							\item{hyperid = }{`55001'}
							\item{name = }{`log precision'}
							\item{short.name = }{`prec'}
							\item{initial = }{`4'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`loggamma'}
							\item{param = }{`1 5e-05'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{doc = }{`A random walk of order 2 for the log-hazard'}
							 }
						 }
					 }
				 }
			 }
		 }
		\item{Section `likelihood'.}{
			 Valid models in this section are:
			 \describe{
				\item{Model `poisson'.}{Number of hyperparmeters are 0.
				 }
				\item{Model `contpoisson'.}{Number of hyperparmeters are 0.
				 }
				\item{Model `qcontpoisson'.}{Number of hyperparmeters are 0.
				 }
				\item{Model `cenpoisson'.}{Number of hyperparmeters are 0.
				 }
				\item{Model `gpoisson'.}{Number of hyperparmeters are 2.
					 \describe{
						\item{Hyperparameter `theta1'}{
							 \describe{
							\item{hyperid = }{`56001'}
							\item{name = }{`overdispersion'}
							\item{short.name = }{`phi'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`loggamma'}
							\item{param = }{`1 1'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
						\item{Hyperparameter `theta2'}{
							 \describe{
							\item{hyperid = }{`56002'}
							\item{name = }{`p'}
							\item{short.name = }{`p'}
							\item{initial = }{`1'}
							\item{fixed = }{`TRUE'}
							\item{prior = }{`normal'}
							\item{param = }{`1 100'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{doc = }{`The generalized Poisson likelihood'}
								\item{survival = }{`FALSE'}
								\item{discrete = }{`TRUE'}
								\item{link = }{`default log logoffset'}
								\item{pdf = }{`gpoisson'}
								\item{status = }{`experimental'}
							 }
						 }
					 }
				 }
				\item{Model `binomial'.}{Number of hyperparmeters are 0.
				 }
				\item{Model `pom'.}{Number of hyperparmeters are 10.
					 \describe{
						\item{Hyperparameter `theta1'}{
							 \describe{
							\item{hyperid = }{`57101'}
							\item{name = }{`theta1'}
							\item{short.name = }{`theta1'}
							\item{initial = }{`NA'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`dirichlet'}
							\item{param = }{`3'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta2'}{
							 \describe{
							\item{hyperid = }{`57102'}
							\item{name = }{`theta2'}
							\item{short.name = }{`theta2'}
							\item{initial = }{`NA'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
						\item{Hyperparameter `theta3'}{
							 \describe{
							\item{hyperid = }{`57103'}
							\item{name = }{`theta3'}
							\item{short.name = }{`theta3'}
							\item{initial = }{`NA'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
						\item{Hyperparameter `theta4'}{
							 \describe{
							\item{hyperid = }{`57104'}
							\item{name = }{`theta4'}
							\item{short.name = }{`theta4'}
							\item{initial = }{`NA'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
						\item{Hyperparameter `theta5'}{
							 \describe{
							\item{hyperid = }{`57105'}
							\item{name = }{`theta5'}
							\item{short.name = }{`theta5'}
							\item{initial = }{`NA'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
						\item{Hyperparameter `theta6'}{
							 \describe{
							\item{hyperid = }{`57106'}
							\item{name = }{`theta6'}
							\item{short.name = }{`theta6'}
							\item{initial = }{`NA'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
						\item{Hyperparameter `theta7'}{
							 \describe{
							\item{hyperid = }{`57107'}
							\item{name = }{`theta7'}
							\item{short.name = }{`theta7'}
							\item{initial = }{`NA'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
						\item{Hyperparameter `theta8'}{
							 \describe{
							\item{hyperid = }{`57108'}
							\item{name = }{`theta8'}
							\item{short.name = }{`theta8'}
							\item{initial = }{`NA'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
						\item{Hyperparameter `theta9'}{
							 \describe{
							\item{hyperid = }{`57109'}
							\item{name = }{`theta9'}
							\item{short.name = }{`theta9'}
							\item{initial = }{`NA'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
						\item{Hyperparameter `theta10'}{
							 \describe{
							\item{hyperid = }{`57110'}
							\item{name = }{`theta10'}
							\item{short.name = }{`theta10'}
							\item{initial = }{`NA'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{doc = }{`Likelihood for the proportional odds model'}
								\item{status = }{`experimental'}
								\item{survival = }{`FALSE'}
								\item{discrete = }{`TRUE'}
								\item{link = }{`default identity'}
								\item{pdf = }{`pom'}
							 }
						 }
					 }
				 }
				\item{Model `gamma'.}{Number of hyperparmeters are 1.
					 \describe{
						\item{Hyperparameter `theta'}{
							 \describe{
							\item{hyperid = }{`58001'}
							\item{name = }{`precision parameter'}
							\item{short.name = }{`prec'}
							\item{initial = }{`4.60517018598809'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`loggamma'}
							\item{param = }{`1 0.01'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{doc = }{`The Gamma likelihood'}
								\item{survival = }{`FALSE'}
								\item{discrete = }{`FALSE'}
								\item{link = }{`default log quantile'}
								\item{pdf = }{`gamma'}
							 }
						 }
					 }
				 }
				\item{Model `gammasurv'.}{Number of hyperparmeters are 1.
					 \describe{
						\item{Hyperparameter `theta'}{
							 \describe{
							\item{hyperid = }{`58101'}
							\item{name = }{`precision parameter'}
							\item{short.name = }{`prec'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`loggamma'}
							\item{param = }{`1 0.01'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{doc = }{`The Gamma likelihood (survival)'}
								\item{survival = }{`TRUE'}
								\item{discrete = }{`FALSE'}
								\item{status = }{`experimental'}
								\item{link = }{`default log quantile'}
								\item{pdf = }{`gammasurv'}
							 }
						 }
					 }
				 }
				\item{Model `gammacount'.}{Number of hyperparmeters are 1.
					 \describe{
						\item{Hyperparameter `theta'}{
							 \describe{
							\item{hyperid = }{`59001'}
							\item{name = }{`log alpha'}
							\item{short.name = }{`alpha'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`pc.gammacount'}
							\item{param = }{`3'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{doc = }{`A Gamma generalisation of the Poisson likelihood'}
								\item{survival = }{`FALSE'}
								\item{discrete = }{`FALSE'}
								\item{link = }{`default log'}
								\item{status = }{`experimental'}
								\item{pdf = }{`gammacount'}
							 }
						 }
					 }
				 }
				\item{Model `qkumar'.}{Number of hyperparmeters are 1.
					 \describe{
						\item{Hyperparameter `theta'}{
							 \describe{
							\item{hyperid = }{`60001'}
							\item{name = }{`precision parameter'}
							\item{short.name = }{`prec'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`loggamma'}
							\item{param = }{`1 0.001'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{doc = }{`A quantile version of the Kumar likelihood'}
								\item{survival = }{`FALSE'}
								\item{discrete = }{`FALSE'}
								\item{link = }{`default logit cauchit'}
								\item{pdf = }{`qkumar'}
							 }
						 }
					 }
				 }
				\item{Model `qloglogistic'.}{Number of hyperparmeters are 1.
					 \describe{
						\item{Hyperparameter `theta'}{
							 \describe{
							\item{hyperid = }{`60011'}
							\item{name = }{`log alpha'}
							\item{short.name = }{`alpha'}
							\item{initial = }{`1'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`loggamma'}
							\item{param = }{`25 25'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{doc = }{`A quantile loglogistic likelihood'}
								\item{status = }{`changed:Oct.25.2017'}
								\item{survival = }{`FALSE'}
								\item{discrete = }{`FALSE'}
								\item{link = }{`default log neglog'}
								\item{pdf = }{`qloglogistic'}
							 }
						 }
					 }
				 }
				\item{Model `qloglogisticsurv'.}{Number of hyperparmeters are 1.
					 \describe{
						\item{Hyperparameter `theta'}{
							 \describe{
							\item{hyperid = }{`60021'}
							\item{name = }{`log alpha'}
							\item{short.name = }{`alpha'}
							\item{initial = }{`1'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`loggamma'}
							\item{param = }{`25 25'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{doc = }{`A quantile loglogistic likelihood (survival)'}
								\item{status = }{`changed:Oct.25.2017'}
								\item{survival = }{`TRUE'}
								\item{discrete = }{`FALSE'}
								\item{link = }{`default log neglog'}
								\item{pdf = }{`qloglogistic'}
							 }
						 }
					 }
				 }
				\item{Model `beta'.}{Number of hyperparmeters are 1.
					 \describe{
						\item{Hyperparameter `theta'}{
							 \describe{
							\item{hyperid = }{`61001'}
							\item{name = }{`precision parameter'}
							\item{short.name = }{`phi'}
							\item{initial = }{`2.30258509299405'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`loggamma'}
							\item{param = }{`1 0.1'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{doc = }{`The Beta likelihood'}
								\item{survival = }{`FALSE'}
								\item{discrete = }{`FALSE'}
								\item{link = }{`default logit cauchit probit cloglog loglog'}
								\item{pdf = }{`beta'}
							 }
						 }
					 }
				 }
				\item{Model `betabinomial'.}{Number of hyperparmeters are 1.
					 \describe{
						\item{Hyperparameter `theta'}{
							 \describe{
							\item{hyperid = }{`62001'}
							\item{name = }{`overdispersion'}
							\item{short.name = }{`rho'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`gaussian'}
							\item{param = }{`0 0.4'}
							\item{to.theta = }{`\code{function(x) log(x/(1-x))}'}
							\item{from.theta = }{`\code{function(x) exp(x)/(1+exp(x))}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{doc = }{`The Beta-Binomial likelihood'}
								\item{survival = }{`FALSE'}
								\item{discrete = }{`TRUE'}
								\item{link = }{`default logit cauchit probit cloglog loglog'}
								\item{pdf = }{`betabinomial'}
							 }
						 }
					 }
				 }
				\item{Model `cbinomial'.}{Number of hyperparmeters are 0.
				 }
				\item{Model `nbinomial'.}{Number of hyperparmeters are 1.
					 \describe{
						\item{Hyperparameter `theta'}{
							 \describe{
							\item{hyperid = }{`63001'}
							\item{name = }{`size'}
							\item{short.name = }{`size'}
							\item{initial = }{`2.30258509299405'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`pc.mgamma'}
							\item{param = }{`7'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{doc = }{`The negBinomial likelihood'}
								\item{survival = }{`FALSE'}
								\item{discrete = }{`TRUE'}
								\item{link = }{`default log logoffset quantile'}
								\item{pdf = }{`nbinomial'}
							 }
						 }
					 }
				 }
				\item{Model `nbinomial2'.}{Number of hyperparmeters are 0.
				 }
				\item{Model `simplex'.}{Number of hyperparmeters are 1.
					 \describe{
						\item{Hyperparameter `theta'}{
							 \describe{
							\item{hyperid = }{`64001'}
							\item{name = }{`log precision'}
							\item{short.name = }{`prec'}
							\item{initial = }{`4'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`loggamma'}
							\item{param = }{`1 5e-05'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{doc = }{`The simplex likelihood'}
								\item{survival = }{`FALSE'}
								\item{discrete = }{`FALSE'}
								\item{link = }{`default logit cauchit probit cloglog loglog'}
								\item{pdf = }{`simplex'}
							 }
						 }
					 }
				 }
				\item{Model `gaussian'.}{Number of hyperparmeters are 2.
					 \describe{
						\item{Hyperparameter `theta1'}{
							 \describe{
							\item{hyperid = }{`65001'}
							\item{name = }{`log precision'}
							\item{short.name = }{`prec'}
							\item{initial = }{`4'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`loggamma'}
							\item{param = }{`1 5e-05'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
						\item{Hyperparameter `theta2'}{
							 \describe{
							\item{hyperid = }{`65002'}
							\item{name = }{`log precision offset'}
							\item{short.name = }{`precoffset'}
							\item{initial = }{`72.0873067782343'}
							\item{fixed = }{`TRUE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{doc = }{`The Gaussian likelihoood'}
								\item{survival = }{`FALSE'}
								\item{discrete = }{`FALSE'}
								\item{link = }{`default identity logit cauchit log logoffset'}
								\item{pdf = }{`gaussian'}
							 }
						 }
					 }
				 }
				\item{Model `circularnormal'.}{Number of hyperparmeters are 1.
					 \describe{
						\item{Hyperparameter `theta'}{
							 \describe{
							\item{hyperid = }{`67001'}
							\item{name = }{`log precision parameter'}
							\item{short.name = }{`prec'}
							\item{initial = }{`2'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`loggamma'}
							\item{param = }{`1 0.01'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{doc = }{`The circular Gaussian likelihoood'}
								\item{survival = }{`FALSE'}
								\item{discrete = }{`FALSE'}
								\item{link = }{`default tan'}
								\item{pdf = }{`circular-normal'}
								\item{status = }{`experimental'}
							 }
						 }
					 }
				 }
				\item{Model `wrappedcauchy'.}{Number of hyperparmeters are 1.
					 \describe{
						\item{Hyperparameter `theta'}{
							 \describe{
							\item{hyperid = }{`68001'}
							\item{name = }{`log precision parameter'}
							\item{short.name = }{`prec'}
							\item{initial = }{`2'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`loggamma'}
							\item{param = }{`1 0.005'}
							\item{to.theta = }{`\code{function(x) log(x/(1-x))}'}
							\item{from.theta = }{`\code{function(x) exp(x)/(1+exp(x))}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{doc = }{`The wrapped Cauchy likelihoood'}
								\item{survival = }{`FALSE'}
								\item{discrete = }{`FALSE'}
								\item{link = }{`default tan'}
								\item{pdf = }{`wrapped-cauchy'}
								\item{status = }{`disabled'}
							 }
						 }
					 }
				 }
				\item{Model `iidgamma'.}{Number of hyperparmeters are 2.
					 \describe{
						\item{Hyperparameter `theta1'}{
							 \describe{
							\item{hyperid = }{`69001'}
							\item{name = }{`logshape'}
							\item{short.name = }{`shape'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`loggamma'}
							\item{param = }{`100 100'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
						\item{Hyperparameter `theta2'}{
							 \describe{
							\item{hyperid = }{`69002'}
							\item{name = }{`lograte'}
							\item{short.name = }{`rate'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`loggamma'}
							\item{param = }{`100 100'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{doc = }{`(experimental)'}
								\item{survival = }{`FALSE'}
								\item{discrete = }{`FALSE'}
								\item{link = }{`default identity'}
								\item{pdf = }{`iidgamma'}
								\item{status = }{`experimental'}
							 }
						 }
					 }
				 }
				\item{Model `iidlogitbeta'.}{Number of hyperparmeters are 2.
					 \describe{
						\item{Hyperparameter `theta1'}{
							 \describe{
							\item{hyperid = }{`70001'}
							\item{name = }{`log.a'}
							\item{short.name = }{`a'}
							\item{initial = }{`1'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`loggamma'}
							\item{param = }{`1 1'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
						\item{Hyperparameter `theta2'}{
							 \describe{
							\item{hyperid = }{`70002'}
							\item{name = }{`log.b'}
							\item{short.name = }{`b'}
							\item{initial = }{`1'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`loggamma'}
							\item{param = }{`1 1'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{doc = }{`(experimental)'}
								\item{survival = }{`FALSE'}
								\item{discrete = }{`FALSE'}
								\item{link = }{`default logit'}
								\item{pdf = }{`iidlogitbeta'}
								\item{status = }{`experimental'}
							 }
						 }
					 }
				 }
				\item{Model `loggammafrailty'.}{Number of hyperparmeters are 1.
					 \describe{
						\item{Hyperparameter `theta'}{
							 \describe{
							\item{hyperid = }{`71001'}
							\item{name = }{`log precision'}
							\item{short.name = }{`prec'}
							\item{initial = }{`4'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`loggamma'}
							\item{param = }{`1 5e-05'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{doc = }{`(experimental)'}
								\item{survival = }{`FALSE'}
								\item{discrete = }{`FALSE'}
								\item{link = }{`default identity'}
								\item{pdf = }{`loggammafrailty'}
								\item{status = }{`experimental'}
							 }
						 }
					 }
				 }
				\item{Model `logistic'.}{Number of hyperparmeters are 1.
					 \describe{
						\item{Hyperparameter `theta'}{
							 \describe{
							\item{hyperid = }{`72001'}
							\item{name = }{`log precision'}
							\item{short.name = }{`prec'}
							\item{initial = }{`1'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`loggamma'}
							\item{param = }{`1 5e-05'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{doc = }{`The Logistic likelihoood'}
								\item{survival = }{`FALSE'}
								\item{discrete = }{`FALSE'}
								\item{link = }{`default identity'}
								\item{pdf = }{`logistic'}
							 }
						 }
					 }
				 }
				\item{Model `skewnormal'.}{Number of hyperparmeters are 2.
					 \describe{
						\item{Hyperparameter `theta1'}{
							 \describe{
							\item{hyperid = }{`73001'}
							\item{name = }{`inverse.scale'}
							\item{short.name = }{`iscale'}
							\item{initial = }{`4'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`loggamma'}
							\item{param = }{`1 5e-05'}
							}
						}
						\item{Hyperparameter `theta2'}{
							 \describe{
							\item{hyperid = }{`73002'}
							\item{name = }{`skewness'}
							\item{short.name = }{`skew'}
							\item{initial = }{`4'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`gaussian'}
							\item{param = }{`0 10'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{doc = }{`The Skew-Normal likelihoood'}
								\item{survival = }{`FALSE'}
								\item{discrete = }{`FALSE'}
								\item{link = }{`default identity'}
								\item{pdf = }{`sn'}
							 }
						 }
					 }
				 }
				\item{Model `sn'.}{Number of hyperparmeters are 2.
					 \describe{
						\item{Hyperparameter `theta1'}{
							 \describe{
							\item{hyperid = }{`74001'}
							\item{name = }{`log inverse scale'}
							\item{short.name = }{`iscale'}
							\item{initial = }{`4'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`loggamma'}
							\item{param = }{`1 5e-05'}
							}
						}
						\item{Hyperparameter `theta2'}{
							 \describe{
							\item{hyperid = }{`74002'}
							\item{name = }{`logit skewness'}
							\item{short.name = }{`skew'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`gaussian'}
							\item{param = }{`0 10'}
							\item{to.theta = }{`\code{function(x, shape.max = 1) log((1+x/shape.max)/(1-x/shape.max))}'}
							\item{from.theta = }{`\code{function(x, shape.max = 1) shape.max*(2*exp(x)/(1+exp(x))-1)}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{doc = }{`The Skew-Normal likelihoood'}
								\item{survival = }{`FALSE'}
								\item{discrete = }{`FALSE'}
								\item{link = }{`default identity'}
								\item{pdf = }{`sn'}
							 }
						 }
					 }
				 }
				\item{Model `sn2'.}{Number of hyperparmeters are 2.
					 \describe{
						\item{Hyperparameter `theta1'}{
							 \describe{
							\item{hyperid = }{`75001'}
							\item{name = }{`log precision'}
							\item{short.name = }{`prec'}
							\item{initial = }{`1'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`loggamma'}
							\item{param = }{`1 5e-05'}
							}
						}
						\item{Hyperparameter `theta2'}{
							 \describe{
							\item{hyperid = }{`75002'}
							\item{name = }{`logit skewness'}
							\item{short.name = }{`skew'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`gaussian'}
							\item{param = }{`0 10'}
							\item{to.theta = }{`\code{function(x) log((1+x)/(1-x))}'}
							\item{from.theta = }{`\code{function(x) (2*exp(x)/(1+exp(x))-1)}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{doc = }{`The Skew-Normal likelihoood (alt param)'}
								\item{survival = }{`FALSE'}
								\item{discrete = }{`FALSE'}
								\item{link = }{`default identity'}
								\item{status = }{`experimental'}
								\item{pdf = }{`sn2'}
							 }
						 }
					 }
				 }
				\item{Model `gev'.}{Number of hyperparmeters are 2.
					 \describe{
						\item{Hyperparameter `theta1'}{
							 \describe{
							\item{hyperid = }{`76001'}
							\item{name = }{`log precision'}
							\item{short.name = }{`prec'}
							\item{initial = }{`4'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`loggamma'}
							\item{param = }{`1 5e-05'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
						\item{Hyperparameter `theta2'}{
							 \describe{
							\item{hyperid = }{`76002'}
							\item{name = }{`gev parameter'}
							\item{short.name = }{`gev'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`gaussian'}
							\item{param = }{`0 25'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{doc = }{`The Generalized Extreme Value likelihood'}
								\item{survival = }{`FALSE'}
								\item{discrete = }{`FALSE'}
								\item{link = }{`default identity'}
								\item{status = }{`experimental'}
								\item{pdf = }{`gev'}
							 }
						 }
					 }
				 }
				\item{Model `lognormal'.}{Number of hyperparmeters are 1.
					 \describe{
						\item{Hyperparameter `theta'}{
							 \describe{
							\item{hyperid = }{`77101'}
							\item{name = }{`log precision'}
							\item{short.name = }{`prec'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`loggamma'}
							\item{param = }{`1 5e-05'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{doc = }{`The log-Normal likelihood'}
								\item{survival = }{`FALSE'}
								\item{discrete = }{`FALSE'}
								\item{link = }{`default identity'}
								\item{pdf = }{`lognormal'}
							 }
						 }
					 }
				 }
				\item{Model `lognormalsurv'.}{Number of hyperparmeters are 1.
					 \describe{
						\item{Hyperparameter `theta'}{
							 \describe{
							\item{hyperid = }{`78001'}
							\item{name = }{`log precision'}
							\item{short.name = }{`prec'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`loggamma'}
							\item{param = }{`1 5e-05'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{doc = }{`The log-Normal likelihood (survival)'}
								\item{survival = }{`TRUE'}
								\item{discrete = }{`FALSE'}
								\item{link = }{`default identity'}
								\item{pdf = }{`lognormal'}
							 }
						 }
					 }
				 }
				\item{Model `exponential'.}{Number of hyperparmeters are 0.
				 }
				\item{Model `exponentialsurv'.}{Number of hyperparmeters are 0.
				 }
				\item{Model `coxph'.}{Number of hyperparmeters are 0.
				 }
				\item{Model `weibull'.}{Number of hyperparmeters are 1.
					 \describe{
						\item{Hyperparameter `theta'}{
							 \describe{
							\item{hyperid = }{`79001'}
							\item{name = }{`log alpha'}
							\item{short.name = }{`alpha'}
							\item{initial = }{`0.1'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`loggamma'}
							\item{param = }{`25 25'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{doc = }{`The Weibull likelihood'}
								\item{survival = }{`FALSE'}
								\item{discrete = }{`FALSE'}
								\item{link = }{`default log neglog quantile'}
								\item{pdf = }{`weibull'}
							 }
						 }
					 }
				 }
				\item{Model `weibullsurv'.}{Number of hyperparmeters are 1.
					 \describe{
						\item{Hyperparameter `theta'}{
							 \describe{
							\item{hyperid = }{`79101'}
							\item{name = }{`log alpha'}
							\item{short.name = }{`alpha'}
							\item{initial = }{`0.1'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`loggamma'}
							\item{param = }{`25 25'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{doc = }{`The Weibull likelihood (survival)'}
								\item{survival = }{`TRUE'}
								\item{discrete = }{`FALSE'}
								\item{link = }{`default log neglog quantile'}
								\item{pdf = }{`weibull'}
							 }
						 }
					 }
				 }
				\item{Model `loglogistic'.}{Number of hyperparmeters are 1.
					 \describe{
						\item{Hyperparameter `theta'}{
							 \describe{
							\item{hyperid = }{`80001'}
							\item{name = }{`log alpha'}
							\item{short.name = }{`alpha'}
							\item{initial = }{`1'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`loggamma'}
							\item{param = }{`25 25'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{doc = }{`The loglogistic likelihood'}
								\item{status = }{`changed:Oct.25.2017'}
								\item{survival = }{`FALSE'}
								\item{discrete = }{`FALSE'}
								\item{link = }{`default log neglog'}
								\item{pdf = }{`loglogistic'}
							 }
						 }
					 }
				 }
				\item{Model `loglogisticsurv'.}{Number of hyperparmeters are 1.
					 \describe{
						\item{Hyperparameter `theta'}{
							 \describe{
							\item{hyperid = }{`80011'}
							\item{name = }{`log alpha'}
							\item{short.name = }{`alpha'}
							\item{initial = }{`1'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`loggamma'}
							\item{param = }{`25 25'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{doc = }{`The loglogistic likelihood (survival)'}
								\item{status = }{`changed:Oct.25.2017'}
								\item{survival = }{`TRUE'}
								\item{discrete = }{`FALSE'}
								\item{link = }{`default log neglog'}
								\item{pdf = }{`loglogistic'}
							 }
						 }
					 }
				 }
				\item{Model `weibullcure'.}{Number of hyperparmeters are 2.
					 \describe{
						\item{Hyperparameter `theta1'}{
							 \describe{
							\item{hyperid = }{`81001'}
							\item{name = }{`log alpha'}
							\item{short.name = }{`a'}
							\item{initial = }{`0.1'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`loggamma'}
							\item{param = }{`25 25'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
						\item{Hyperparameter `theta2'}{
							 \describe{
							\item{hyperid = }{`81002'}
							\item{name = }{`logit probability'}
							\item{short.name = }{`prob'}
							\item{initial = }{`-1'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`gaussian'}
							\item{param = }{`-1 0.2'}
							\item{to.theta = }{`\code{function(x) log(x/(1-x))}'}
							\item{from.theta = }{`\code{function(x) exp(x)/(1+exp(x))}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{doc = }{`The Weibull-cure likelihood (survival)'}
								\item{survival = }{`TRUE'}
								\item{discrete = }{`FALSE'}
								\item{link = }{`default log neglog'}
								\item{pdf = }{`weibullcure'}
							 }
						 }
					 }
				 }
				\item{Model `stochvol'.}{Number of hyperparmeters are 1.
					 \describe{
						\item{Hyperparameter `theta'}{
							 \describe{
							\item{hyperid = }{`82001'}
							\item{name = }{`log precision'}
							\item{short.name = }{`prec'}
							\item{initial = }{`500'}
							\item{fixed = }{`TRUE'}
							\item{prior = }{`loggamma'}
							\item{param = }{`1 0.005'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{doc = }{`The Gaussian stochvol likelihood'}
								\item{survival = }{`FALSE'}
								\item{discrete = }{`FALSE'}
								\item{link = }{`default log'}
								\item{pdf = }{`stochvolgaussian'}
							 }
						 }
					 }
				 }
				\item{Model `stochvolt'.}{Number of hyperparmeters are 1.
					 \describe{
						\item{Hyperparameter `theta'}{
							 \describe{
							\item{hyperid = }{`83001'}
							\item{name = }{`log degrees of freedom'}
							\item{short.name = }{`dof'}
							\item{initial = }{`4'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`pc.dof'}
							\item{param = }{`15 0.5'}
							\item{to.theta = }{`\code{function(x) log(x-2)}'}
							\item{from.theta = }{`\code{function(x) 2+exp(x)}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{doc = }{`The Student-t stochvol likelihood'}
								\item{survival = }{`FALSE'}
								\item{discrete = }{`FALSE'}
								\item{link = }{`default log'}
								\item{pdf = }{`stochvolt'}
							 }
						 }
					 }
				 }
				\item{Model `stochvolnig'.}{Number of hyperparmeters are 2.
					 \describe{
						\item{Hyperparameter `theta1'}{
							 \describe{
							\item{hyperid = }{`84001'}
							\item{name = }{`skewness'}
							\item{short.name = }{`skew'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`gaussian'}
							\item{param = }{`0 10'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta2'}{
							 \describe{
							\item{hyperid = }{`84002'}
							\item{name = }{`shape'}
							\item{short.name = }{`shape'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`loggamma'}
							\item{param = }{`1 0.5'}
							\item{to.theta = }{`\code{function(x) log(x-1)}'}
							\item{from.theta = }{`\code{function(x) 1+exp(x)}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{doc = }{`The Normal inverse Gaussian stochvol likelihood'}
								\item{survival = }{`FALSE'}
								\item{discrete = }{`FALSE'}
								\item{link = }{`default log'}
								\item{pdf = }{`stochvolnig'}
							 }
						 }
					 }
				 }
				\item{Model `zeroinflatedpoisson0'.}{Number of hyperparmeters are 1.
					 \describe{
						\item{Hyperparameter `theta'}{
							 \describe{
							\item{hyperid = }{`85001'}
							\item{name = }{`logit probability'}
							\item{short.name = }{`prob'}
							\item{initial = }{`-1'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`gaussian'}
							\item{param = }{`-1 0.2'}
							\item{to.theta = }{`\code{function(x) log(x/(1-x))}'}
							\item{from.theta = }{`\code{function(x) exp(x)/(1+exp(x))}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{doc = }{`Zero-inflated Poisson, type 0'}
								\item{survival = }{`FALSE'}
								\item{discrete = }{`FALSE'}
								\item{link = }{`default log'}
								\item{pdf = }{`zeroinflated'}
							 }
						 }
					 }
				 }
				\item{Model `zeroinflatedpoisson1'.}{Number of hyperparmeters are 1.
					 \describe{
						\item{Hyperparameter `theta'}{
							 \describe{
							\item{hyperid = }{`86001'}
							\item{name = }{`logit probability'}
							\item{short.name = }{`prob'}
							\item{initial = }{`-1'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`gaussian'}
							\item{param = }{`-1 0.2'}
							\item{to.theta = }{`\code{function(x) log(x/(1-x))}'}
							\item{from.theta = }{`\code{function(x) exp(x)/(1+exp(x))}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{doc = }{`Zero-inflated Poisson, type 1'}
								\item{survival = }{`FALSE'}
								\item{discrete = }{`FALSE'}
								\item{link = }{`default log'}
								\item{pdf = }{`zeroinflated'}
							 }
						 }
					 }
				 }
				\item{Model `zeroinflatedpoisson2'.}{Number of hyperparmeters are 1.
					 \describe{
						\item{Hyperparameter `theta'}{
							 \describe{
							\item{hyperid = }{`87001'}
							\item{name = }{`log alpha'}
							\item{short.name = }{`a'}
							\item{initial = }{`0.693147180559945'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`gaussian'}
							\item{param = }{`0.693147180559945 1'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{doc = }{`Zero-inflated Poisson, type 2'}
								\item{survival = }{`FALSE'}
								\item{discrete = }{`FALSE'}
								\item{link = }{`default log'}
								\item{pdf = }{`zeroinflated'}
							 }
						 }
					 }
				 }
				\item{Model `zeroinflatedbetabinomial0'.}{Number of hyperparmeters are 2.
					 \describe{
						\item{Hyperparameter `theta1'}{
							 \describe{
							\item{hyperid = }{`88001'}
							\item{name = }{`overdispersion'}
							\item{short.name = }{`rho'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`gaussian'}
							\item{param = }{`0 0.4'}
							\item{to.theta = }{`\code{function(x) log(x/(1-x))}'}
							\item{from.theta = }{`\code{function(x) exp(x)/(1+exp(x))}'}
							}
						}
						\item{Hyperparameter `theta2'}{
							 \describe{
							\item{hyperid = }{`88002'}
							\item{name = }{`logit probability'}
							\item{short.name = }{`prob'}
							\item{initial = }{`-1'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`gaussian'}
							\item{param = }{`-1 0.2'}
							\item{to.theta = }{`\code{function(x) log(x/(1-x))}'}
							\item{from.theta = }{`\code{function(x) exp(x)/(1+exp(x))}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{doc = }{`Zero-inflated Beta-Binomial, type 0'}
								\item{survival = }{`FALSE'}
								\item{discrete = }{`TRUE'}
								\item{link = }{`default logit cauchit probit cloglog loglog'}
								\item{pdf = }{`zeroinflated'}
							 }
						 }
					 }
				 }
				\item{Model `zeroinflatedbetabinomial1'.}{Number of hyperparmeters are 2.
					 \describe{
						\item{Hyperparameter `theta1'}{
							 \describe{
							\item{hyperid = }{`89001'}
							\item{name = }{`overdispersion'}
							\item{short.name = }{`rho'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`gaussian'}
							\item{param = }{`0 0.4'}
							\item{to.theta = }{`\code{function(x) log(x/(1-x))}'}
							\item{from.theta = }{`\code{function(x) exp(x)/(1+exp(x))}'}
							}
						}
						\item{Hyperparameter `theta2'}{
							 \describe{
							\item{hyperid = }{`89002'}
							\item{name = }{`logit probability'}
							\item{short.name = }{`prob'}
							\item{initial = }{`-1'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`gaussian'}
							\item{param = }{`-1 0.2'}
							\item{to.theta = }{`\code{function(x) log(x/(1-x))}'}
							\item{from.theta = }{`\code{function(x) exp(x)/(1+exp(x))}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{doc = }{`Zero-inflated Beta-Binomial, type 1'}
								\item{survival = }{`FALSE'}
								\item{discrete = }{`TRUE'}
								\item{link = }{`default logit cauchit probit cloglog loglog'}
								\item{pdf = }{`zeroinflated'}
							 }
						 }
					 }
				 }
				\item{Model `zeroinflatedbinomial0'.}{Number of hyperparmeters are 1.
					 \describe{
						\item{Hyperparameter `theta'}{
							 \describe{
							\item{hyperid = }{`90001'}
							\item{name = }{`logit probability'}
							\item{short.name = }{`prob'}
							\item{initial = }{`-1'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`gaussian'}
							\item{param = }{`-1 0.2'}
							\item{to.theta = }{`\code{function(x) log(x/(1-x))}'}
							\item{from.theta = }{`\code{function(x) exp(x)/(1+exp(x))}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{doc = }{`Zero-inflated Binomial, type 0'}
								\item{survival = }{`FALSE'}
								\item{discrete = }{`FALSE'}
								\item{link = }{`default logit cauchit probit cloglog loglog'}
								\item{pdf = }{`zeroinflated'}
							 }
						 }
					 }
				 }
				\item{Model `zeroinflatedbinomial1'.}{Number of hyperparmeters are 1.
					 \describe{
						\item{Hyperparameter `theta'}{
							 \describe{
							\item{hyperid = }{`91001'}
							\item{name = }{`logit probability'}
							\item{short.name = }{`prob'}
							\item{initial = }{`-1'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`gaussian'}
							\item{param = }{`-1 0.2'}
							\item{to.theta = }{`\code{function(x) log(x/(1-x))}'}
							\item{from.theta = }{`\code{function(x) exp(x)/(1+exp(x))}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{doc = }{`Zero-inflated Binomial, type 1'}
								\item{survival = }{`FALSE'}
								\item{discrete = }{`FALSE'}
								\item{link = }{`default logit cauchit probit cloglog loglog'}
								\item{pdf = }{`zeroinflated'}
							 }
						 }
					 }
				 }
				\item{Model `zeroinflatedbinomial2'.}{Number of hyperparmeters are 1.
					 \describe{
						\item{Hyperparameter `theta'}{
							 \describe{
							\item{hyperid = }{`92001'}
							\item{name = }{`alpha'}
							\item{short.name = }{`alpha'}
							\item{initial = }{`-1'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`gaussian'}
							\item{param = }{`-1 0.2'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{doc = }{`Zero-inflated Binomial, type 2'}
								\item{survival = }{`FALSE'}
								\item{discrete = }{`FALSE'}
								\item{link = }{`default logit cauchit probit cloglog loglog'}
								\item{pdf = }{`zeroinflated'}
							 }
						 }
					 }
				 }
				\item{Model `zeroninflatedbinomial2'.}{Number of hyperparmeters are 2.
					 \describe{
						\item{Hyperparameter `theta1'}{
							 \describe{
							\item{hyperid = }{`93001'}
							\item{name = }{`alpha1'}
							\item{short.name = }{`alpha1'}
							\item{initial = }{`-1'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`gaussian'}
							\item{param = }{`-1 0.2'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
						\item{Hyperparameter `theta2'}{
							 \describe{
							\item{hyperid = }{`93002'}
							\item{name = }{`alpha2'}
							\item{short.name = }{`alpha2'}
							\item{initial = }{`-1'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`gaussian'}
							\item{param = }{`-1 0.2'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{doc = }{`Zero and N inflated binomial, type 2'}
								\item{survival = }{`FALSE'}
								\item{discrete = }{`FALSE'}
								\item{link = }{`default logit cauchit probit cloglog loglog'}
								\item{pdf = }{`NA'}
							 }
						 }
					 }
				 }
				\item{Model `zeroninflatedbinomial3'.}{Number of hyperparmeters are 2.
					 \describe{
						\item{Hyperparameter `theta1'}{
							 \describe{
							\item{hyperid = }{`93101'}
							\item{name = }{`alpha0'}
							\item{short.name = }{`alpha0'}
							\item{initial = }{`1'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`loggamma'}
							\item{param = }{`1 1'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
						\item{Hyperparameter `theta2'}{
							 \describe{
							\item{hyperid = }{`93102'}
							\item{name = }{`alphaN'}
							\item{short.name = }{`alphaN'}
							\item{initial = }{`1'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`loggamma'}
							\item{param = }{`1 1'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{doc = }{`Zero and N inflated binomial, type 3'}
								\item{status = }{`experimental'}
								\item{survival = }{`FALSE'}
								\item{discrete = }{`FALSE'}
								\item{link = }{`default logit cauchit probit cloglog loglog'}
								\item{pdf = }{`zeroinflated'}
							 }
						 }
					 }
				 }
				\item{Model `zeroinflatedbetabinomial2'.}{Number of hyperparmeters are 2.
					 \describe{
						\item{Hyperparameter `theta1'}{
							 \describe{
							\item{hyperid = }{`94001'}
							\item{name = }{`log alpha'}
							\item{short.name = }{`a'}
							\item{initial = }{`0.693147180559945'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`gaussian'}
							\item{param = }{`0.693147180559945 1'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
						\item{Hyperparameter `theta2'}{
							 \describe{
							\item{hyperid = }{`94002'}
							\item{name = }{`beta'}
							\item{short.name = }{`b'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`gaussian'}
							\item{param = }{`0 1'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{doc = }{`Zero inflated Beta-Binomial, type 2'}
								\item{survival = }{`FALSE'}
								\item{discrete = }{`FALSE'}
								\item{link = }{`default logit cauchit probit cloglog loglog'}
								\item{pdf = }{`zeroinflated'}
							 }
						 }
					 }
				 }
				\item{Model `zeroinflatednbinomial0'.}{Number of hyperparmeters are 2.
					 \describe{
						\item{Hyperparameter `theta1'}{
							 \describe{
							\item{hyperid = }{`95001'}
							\item{name = }{`log size'}
							\item{short.name = }{`size'}
							\item{initial = }{`2.30258509299405'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`pc.mgamma'}
							\item{param = }{`7'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
						\item{Hyperparameter `theta2'}{
							 \describe{
							\item{hyperid = }{`95002'}
							\item{name = }{`logit probability'}
							\item{short.name = }{`prob'}
							\item{initial = }{`-1'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`gaussian'}
							\item{param = }{`-1 0.2'}
							\item{to.theta = }{`\code{function(x) log(x/(1-x))}'}
							\item{from.theta = }{`\code{function(x) exp(x)/(1+exp(x))}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{doc = }{`Zero inflated negBinomial, type 0'}
								\item{survival = }{`FALSE'}
								\item{discrete = }{`FALSE'}
								\item{link = }{`default log'}
								\item{pdf = }{`zeroinflated'}
							 }
						 }
					 }
				 }
				\item{Model `zeroinflatednbinomial1'.}{Number of hyperparmeters are 2.
					 \describe{
						\item{Hyperparameter `theta1'}{
							 \describe{
							\item{hyperid = }{`96001'}
							\item{name = }{`log size'}
							\item{short.name = }{`size'}
							\item{initial = }{`2.30258509299405'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`pc.mgamma'}
							\item{param = }{`7'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
						\item{Hyperparameter `theta2'}{
							 \describe{
							\item{hyperid = }{`96002'}
							\item{name = }{`logit probability'}
							\item{short.name = }{`prob'}
							\item{initial = }{`-1'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`gaussian'}
							\item{param = }{`-1 0.2'}
							\item{to.theta = }{`\code{function(x) log(x/(1-x))}'}
							\item{from.theta = }{`\code{function(x) exp(x)/(1+exp(x))}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{doc = }{`Zero inflated negBinomial, type 1'}
								\item{survival = }{`FALSE'}
								\item{discrete = }{`FALSE'}
								\item{link = }{`default log'}
								\item{pdf = }{`zeroinflated'}
							 }
						 }
					 }
				 }
				\item{Model `zeroinflatednbinomial1strata2'.}{Number of hyperparmeters are 11.
					 \describe{
						\item{Hyperparameter `theta1'}{
							 \describe{
							\item{hyperid = }{`97001'}
							\item{name = }{`log size'}
							\item{short.name = }{`size'}
							\item{initial = }{`2.30258509299405'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`pc.mgamma'}
							\item{param = }{`7'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
						\item{Hyperparameter `theta2'}{
							 \describe{
							\item{hyperid = }{`97002'}
							\item{name = }{`logit probability 1'}
							\item{short.name = }{`prob1'}
							\item{initial = }{`-1'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`gaussian'}
							\item{param = }{`-1 0.2'}
							\item{to.theta = }{`\code{function(x) log(x/(1-x))}'}
							\item{from.theta = }{`\code{function(x) exp(x)/(1+exp(x))}'}
							}
						}
						\item{Hyperparameter `theta3'}{
							 \describe{
							\item{hyperid = }{`97003'}
							\item{name = }{`logit probability 2'}
							\item{short.name = }{`prob2'}
							\item{initial = }{`-1'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`gaussian'}
							\item{param = }{`-1 0.2'}
							\item{to.theta = }{`\code{function(x) log(x/(1-x))}'}
							\item{from.theta = }{`\code{function(x) exp(x)/(1+exp(x))}'}
							}
						}
						\item{Hyperparameter `theta4'}{
							 \describe{
							\item{hyperid = }{`97004'}
							\item{name = }{`logit probability 3'}
							\item{short.name = }{`prob3'}
							\item{initial = }{`-1'}
							\item{fixed = }{`TRUE'}
							\item{prior = }{`gaussian'}
							\item{param = }{`-1 0.2'}
							\item{to.theta = }{`\code{function(x) log(x/(1-x))}'}
							\item{from.theta = }{`\code{function(x) exp(x)/(1+exp(x))}'}
							}
						}
						\item{Hyperparameter `theta5'}{
							 \describe{
							\item{hyperid = }{`97005'}
							\item{name = }{`logit probability 4'}
							\item{short.name = }{`prob4'}
							\item{initial = }{`-1'}
							\item{fixed = }{`TRUE'}
							\item{prior = }{`gaussian'}
							\item{param = }{`-1 0.2'}
							\item{to.theta = }{`\code{function(x) log(x/(1-x))}'}
							\item{from.theta = }{`\code{function(x) exp(x)/(1+exp(x))}'}
							}
						}
						\item{Hyperparameter `theta6'}{
							 \describe{
							\item{hyperid = }{`97006'}
							\item{name = }{`logit probability 5'}
							\item{short.name = }{`prob5'}
							\item{initial = }{`-1'}
							\item{fixed = }{`TRUE'}
							\item{prior = }{`gaussian'}
							\item{param = }{`-1 0.2'}
							\item{to.theta = }{`\code{function(x) log(x/(1-x))}'}
							\item{from.theta = }{`\code{function(x) exp(x)/(1+exp(x))}'}
							}
						}
						\item{Hyperparameter `theta7'}{
							 \describe{
							\item{hyperid = }{`97007'}
							\item{name = }{`logit probability 6'}
							\item{short.name = }{`prob6'}
							\item{initial = }{`-1'}
							\item{fixed = }{`TRUE'}
							\item{prior = }{`gaussian'}
							\item{param = }{`-1 0.2'}
							\item{to.theta = }{`\code{function(x) log(x/(1-x))}'}
							\item{from.theta = }{`\code{function(x) exp(x)/(1+exp(x))}'}
							}
						}
						\item{Hyperparameter `theta8'}{
							 \describe{
							\item{hyperid = }{`97008'}
							\item{name = }{`logit probability 7'}
							\item{short.name = }{`prob7'}
							\item{initial = }{`-1'}
							\item{fixed = }{`TRUE'}
							\item{prior = }{`gaussian'}
							\item{param = }{`-1 0.2'}
							\item{to.theta = }{`\code{function(x) log(x/(1-x))}'}
							\item{from.theta = }{`\code{function(x) exp(x)/(1+exp(x))}'}
							}
						}
						\item{Hyperparameter `theta9'}{
							 \describe{
							\item{hyperid = }{`97009'}
							\item{name = }{`logit probability 8'}
							\item{short.name = }{`prob8'}
							\item{initial = }{`-1'}
							\item{fixed = }{`TRUE'}
							\item{prior = }{`gaussian'}
							\item{param = }{`-1 0.2'}
							\item{to.theta = }{`\code{function(x) log(x/(1-x))}'}
							\item{from.theta = }{`\code{function(x) exp(x)/(1+exp(x))}'}
							}
						}
						\item{Hyperparameter `theta10'}{
							 \describe{
							\item{hyperid = }{`97010'}
							\item{name = }{`logit probability 9'}
							\item{short.name = }{`prob9'}
							\item{initial = }{`-1'}
							\item{fixed = }{`TRUE'}
							\item{prior = }{`gaussian'}
							\item{param = }{`-1 0.2'}
							\item{to.theta = }{`\code{function(x) log(x/(1-x))}'}
							\item{from.theta = }{`\code{function(x) exp(x)/(1+exp(x))}'}
							}
						}
						\item{Hyperparameter `theta11'}{
							 \describe{
							\item{hyperid = }{`97011'}
							\item{name = }{`logit probability 10'}
							\item{short.name = }{`prob10'}
							\item{initial = }{`-1'}
							\item{fixed = }{`TRUE'}
							\item{prior = }{`gaussian'}
							\item{param = }{`-1 0.2'}
							\item{to.theta = }{`\code{function(x) log(x/(1-x))}'}
							\item{from.theta = }{`\code{function(x) exp(x)/(1+exp(x))}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{doc = }{`Zero inflated negBinomial, type 1, strata 2'}
								\item{status = }{`experimental'}
								\item{survival = }{`FALSE'}
								\item{discrete = }{`FALSE'}
								\item{link = }{`default log'}
								\item{pdf = }{`zeroinflated'}
							 }
						 }
					 }
				 }
				\item{Model `zeroinflatednbinomial1strata3'.}{Number of hyperparmeters are 11.
					 \describe{
						\item{Hyperparameter `theta1'}{
							 \describe{
							\item{hyperid = }{`98001'}
							\item{name = }{`logit probability'}
							\item{short.name = }{`prob'}
							\item{initial = }{`-1'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`gaussian'}
							\item{param = }{`-1 0.2'}
							\item{to.theta = }{`\code{function(x) log(x/(1-x))}'}
							\item{from.theta = }{`\code{function(x) exp(x)/(1+exp(x))}'}
							}
						}
						\item{Hyperparameter `theta2'}{
							 \describe{
							\item{hyperid = }{`98002'}
							\item{name = }{`log size 1'}
							\item{short.name = }{`size1'}
							\item{initial = }{`2.30258509299405'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`pc.mgamma'}
							\item{param = }{`7'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
						\item{Hyperparameter `theta3'}{
							 \describe{
							\item{hyperid = }{`98003'}
							\item{name = }{`log size 2'}
							\item{short.name = }{`size2'}
							\item{initial = }{`2.30258509299405'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`pc.mgamma'}
							\item{param = }{`7'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
						\item{Hyperparameter `theta4'}{
							 \describe{
							\item{hyperid = }{`98004'}
							\item{name = }{`log size 3'}
							\item{short.name = }{`size3'}
							\item{initial = }{`2.30258509299405'}
							\item{fixed = }{`TRUE'}
							\item{prior = }{`pc.mgamma'}
							\item{param = }{`7'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
						\item{Hyperparameter `theta5'}{
							 \describe{
							\item{hyperid = }{`98005'}
							\item{name = }{`log size 4'}
							\item{short.name = }{`size4'}
							\item{initial = }{`2.30258509299405'}
							\item{fixed = }{`TRUE'}
							\item{prior = }{`pc.mgamma'}
							\item{param = }{`7'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
						\item{Hyperparameter `theta6'}{
							 \describe{
							\item{hyperid = }{`98006'}
							\item{name = }{`log size 5'}
							\item{short.name = }{`size5'}
							\item{initial = }{`2.30258509299405'}
							\item{fixed = }{`TRUE'}
							\item{prior = }{`pc.mgamma'}
							\item{param = }{`7'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
						\item{Hyperparameter `theta7'}{
							 \describe{
							\item{hyperid = }{`98007'}
							\item{name = }{`log size 6'}
							\item{short.name = }{`size6'}
							\item{initial = }{`2.30258509299405'}
							\item{fixed = }{`TRUE'}
							\item{prior = }{`pc.mgamma'}
							\item{param = }{`7'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
						\item{Hyperparameter `theta8'}{
							 \describe{
							\item{hyperid = }{`98008'}
							\item{name = }{`log size 7'}
							\item{short.name = }{`size7'}
							\item{initial = }{`2.30258509299405'}
							\item{fixed = }{`TRUE'}
							\item{prior = }{`pc.mgamma'}
							\item{param = }{`7'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
						\item{Hyperparameter `theta9'}{
							 \describe{
							\item{hyperid = }{`98009'}
							\item{name = }{`log size 8'}
							\item{short.name = }{`size8'}
							\item{initial = }{`2.30258509299405'}
							\item{fixed = }{`TRUE'}
							\item{prior = }{`pc.mgamma'}
							\item{param = }{`7'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
						\item{Hyperparameter `theta10'}{
							 \describe{
							\item{hyperid = }{`98010'}
							\item{name = }{`log size 9'}
							\item{short.name = }{`size9'}
							\item{initial = }{`2.30258509299405'}
							\item{fixed = }{`TRUE'}
							\item{prior = }{`pc.mgamma'}
							\item{param = }{`7'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
						\item{Hyperparameter `theta11'}{
							 \describe{
							\item{hyperid = }{`98011'}
							\item{name = }{`log size 10'}
							\item{short.name = }{`size10'}
							\item{initial = }{`2.30258509299405'}
							\item{fixed = }{`TRUE'}
							\item{prior = }{`pc.mgamma'}
							\item{param = }{`7'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{doc = }{`Zero inflated negBinomial, type 1, strata 3'}
								\item{status = }{`experimental'}
								\item{survival = }{`FALSE'}
								\item{discrete = }{`FALSE'}
								\item{link = }{`default log'}
								\item{pdf = }{`zeroinflated'}
							 }
						 }
					 }
				 }
				\item{Model `zeroinflatednbinomial2'.}{Number of hyperparmeters are 2.
					 \describe{
						\item{Hyperparameter `theta1'}{
							 \describe{
							\item{hyperid = }{`99001'}
							\item{name = }{`log size'}
							\item{short.name = }{`size'}
							\item{initial = }{`2.30258509299405'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`pc.mgamma'}
							\item{param = }{`7'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
						\item{Hyperparameter `theta2'}{
							 \describe{
							\item{hyperid = }{`99002'}
							\item{name = }{`log alpha'}
							\item{short.name = }{`a'}
							\item{initial = }{`0.693147180559945'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`gaussian'}
							\item{param = }{`2 1'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{doc = }{`Zero inflated negBinomial, type 2'}
								\item{survival = }{`FALSE'}
								\item{discrete = }{`FALSE'}
								\item{link = }{`default log'}
								\item{pdf = }{`zeroinflated'}
							 }
						 }
					 }
				 }
				\item{Model `t'.}{Number of hyperparmeters are 2.
					 \describe{
						\item{Hyperparameter `theta1'}{
							 \describe{
							\item{hyperid = }{`100001'}
							\item{name = }{`log precision'}
							\item{short.name = }{`prec'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`loggamma'}
							\item{param = }{`1 5e-05'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
						\item{Hyperparameter `theta2'}{
							 \describe{
							\item{hyperid = }{`100002'}
							\item{name = }{`log degrees of freedom'}
							\item{short.name = }{`dof'}
							\item{initial = }{`5'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`pc.dof'}
							\item{param = }{`15 0.5'}
							\item{to.theta = }{`\code{function(x) log(x-2)}'}
							\item{from.theta = }{`\code{function(x) 2+exp(x)}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{doc = }{`Student-t likelihood'}
								\item{survival = }{`FALSE'}
								\item{discrete = }{`FALSE'}
								\item{link = }{`default identity'}
								\item{pdf = }{`student-t'}
							 }
						 }
					 }
				 }
				\item{Model `tstrata'.}{Number of hyperparmeters are 11.
					 \describe{
						\item{Hyperparameter `theta1'}{
							 \describe{
							\item{hyperid = }{`101001'}
							\item{name = }{`log degrees of freedom'}
							\item{short.name = }{`dof'}
							\item{initial = }{`4'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`pc.dof'}
							\item{param = }{`15 0.5'}
							\item{to.theta = }{`\code{function(x) log(x-5)}'}
							\item{from.theta = }{`\code{function(x) 5+exp(x)}'}
							}
						}
						\item{Hyperparameter `theta2'}{
							 \describe{
							\item{hyperid = }{`101002'}
							\item{name = }{`log precision1'}
							\item{short.name = }{`prec1'}
							\item{initial = }{`2'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`loggamma'}
							\item{param = }{`1 5e-05'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
						\item{Hyperparameter `theta3'}{
							 \describe{
							\item{hyperid = }{`101003'}
							\item{name = }{`log precision2'}
							\item{short.name = }{`prec2'}
							\item{initial = }{`2'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`loggamma'}
							\item{param = }{`1 5e-05'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
						\item{Hyperparameter `theta4'}{
							 \describe{
							\item{hyperid = }{`101004'}
							\item{name = }{`log precision3'}
							\item{short.name = }{`prec3'}
							\item{initial = }{`2'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`loggamma'}
							\item{param = }{`1 5e-05'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
						\item{Hyperparameter `theta5'}{
							 \describe{
							\item{hyperid = }{`101005'}
							\item{name = }{`log precision4'}
							\item{short.name = }{`prec4'}
							\item{initial = }{`2'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`loggamma'}
							\item{param = }{`1 5e-05'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
						\item{Hyperparameter `theta6'}{
							 \describe{
							\item{hyperid = }{`101006'}
							\item{name = }{`log precision5'}
							\item{short.name = }{`prec5'}
							\item{initial = }{`2'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`loggamma'}
							\item{param = }{`1 5e-05'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
						\item{Hyperparameter `theta7'}{
							 \describe{
							\item{hyperid = }{`101007'}
							\item{name = }{`log precision6'}
							\item{short.name = }{`prec6'}
							\item{initial = }{`2'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`loggamma'}
							\item{param = }{`1 5e-05'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
						\item{Hyperparameter `theta8'}{
							 \describe{
							\item{hyperid = }{`101008'}
							\item{name = }{`log precision7'}
							\item{short.name = }{`prec7'}
							\item{initial = }{`2'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`loggamma'}
							\item{param = }{`1 5e-05'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
						\item{Hyperparameter `theta9'}{
							 \describe{
							\item{hyperid = }{`101009'}
							\item{name = }{`log precision8'}
							\item{short.name = }{`prec8'}
							\item{initial = }{`2'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`loggamma'}
							\item{param = }{`1 5e-05'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
						\item{Hyperparameter `theta10'}{
							 \describe{
							\item{hyperid = }{`101010'}
							\item{name = }{`log precision9'}
							\item{short.name = }{`prec9'}
							\item{initial = }{`2'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`loggamma'}
							\item{param = }{`1 5e-05'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
						\item{Hyperparameter `theta11'}{
							 \describe{
							\item{hyperid = }{`101011'}
							\item{name = }{`log precision10'}
							\item{short.name = }{`prec10'}
							\item{initial = }{`2'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`loggamma'}
							\item{param = }{`1 5e-05'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{doc = }{`A stratified version of the Student-t likelihood'}
								\item{survival = }{`FALSE'}
								\item{discrete = }{`FALSE'}
								\item{link = }{`default identity'}
								\item{pdf = }{`tstrata'}
							 }
						 }
					 }
				 }
				\item{Model `nmix'.}{Number of hyperparmeters are 15.
					 \describe{
						\item{Hyperparameter `theta1'}{
							 \describe{
							\item{hyperid = }{`101101'}
							\item{name = }{`beta1'}
							\item{short.name = }{`beta1'}
							\item{initial = }{`2.30258509299405'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`normal'}
							\item{param = }{`0 0.5'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta2'}{
							 \describe{
							\item{hyperid = }{`101102'}
							\item{name = }{`beta2'}
							\item{short.name = }{`beta2'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`normal'}
							\item{param = }{`0 1'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta3'}{
							 \describe{
							\item{hyperid = }{`101103'}
							\item{name = }{`beta3'}
							\item{short.name = }{`beta3'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`normal'}
							\item{param = }{`0 1'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta4'}{
							 \describe{
							\item{hyperid = }{`101104'}
							\item{name = }{`beta4'}
							\item{short.name = }{`beta4'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`normal'}
							\item{param = }{`0 1'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta5'}{
							 \describe{
							\item{hyperid = }{`101105'}
							\item{name = }{`beta5'}
							\item{short.name = }{`beta5'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`normal'}
							\item{param = }{`0 1'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta6'}{
							 \describe{
							\item{hyperid = }{`101106'}
							\item{name = }{`beta6'}
							\item{short.name = }{`beta6'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`normal'}
							\item{param = }{`0 1'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta7'}{
							 \describe{
							\item{hyperid = }{`101107'}
							\item{name = }{`beta7'}
							\item{short.name = }{`beta7'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`normal'}
							\item{param = }{`0 1'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta8'}{
							 \describe{
							\item{hyperid = }{`101108'}
							\item{name = }{`beta8'}
							\item{short.name = }{`beta8'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`normal'}
							\item{param = }{`0 1'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta9'}{
							 \describe{
							\item{hyperid = }{`101109'}
							\item{name = }{`beta9'}
							\item{short.name = }{`beta9'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`normal'}
							\item{param = }{`0 1'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta10'}{
							 \describe{
							\item{hyperid = }{`101110'}
							\item{name = }{`beta10'}
							\item{short.name = }{`beta10'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`normal'}
							\item{param = }{`0 1'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta11'}{
							 \describe{
							\item{hyperid = }{`101111'}
							\item{name = }{`beta11'}
							\item{short.name = }{`beta11'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`normal'}
							\item{param = }{`0 1'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta12'}{
							 \describe{
							\item{hyperid = }{`101112'}
							\item{name = }{`beta12'}
							\item{short.name = }{`beta12'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`normal'}
							\item{param = }{`0 1'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta13'}{
							 \describe{
							\item{hyperid = }{`101113'}
							\item{name = }{`beta13'}
							\item{short.name = }{`beta13'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`normal'}
							\item{param = }{`0 1'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta14'}{
							 \describe{
							\item{hyperid = }{`101114'}
							\item{name = }{`beta14'}
							\item{short.name = }{`beta14'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`normal'}
							\item{param = }{`0 1'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta15'}{
							 \describe{
							\item{hyperid = }{`101115'}
							\item{name = }{`beta15'}
							\item{short.name = }{`beta15'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`normal'}
							\item{param = }{`0 1'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{doc = }{`Binomial-Poisson mixture'}
								\item{status = }{`experimental'}
								\item{survival = }{`FALSE'}
								\item{discrete = }{`TRUE'}
								\item{link = }{`default logit probit'}
								\item{pdf = }{`nmix'}
							 }
						 }
					 }
				 }
				\item{Model `nmixnb'.}{Number of hyperparmeters are 16.
					 \describe{
						\item{Hyperparameter `theta1'}{
							 \describe{
							\item{hyperid = }{`101121'}
							\item{name = }{`beta1'}
							\item{short.name = }{`beta1'}
							\item{initial = }{`2.30258509299405'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`normal'}
							\item{param = }{`0 0.5'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta2'}{
							 \describe{
							\item{hyperid = }{`101122'}
							\item{name = }{`beta2'}
							\item{short.name = }{`beta2'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`normal'}
							\item{param = }{`0 1'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta3'}{
							 \describe{
							\item{hyperid = }{`101123'}
							\item{name = }{`beta3'}
							\item{short.name = }{`beta3'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`normal'}
							\item{param = }{`0 1'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta4'}{
							 \describe{
							\item{hyperid = }{`101124'}
							\item{name = }{`beta4'}
							\item{short.name = }{`beta4'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`normal'}
							\item{param = }{`0 1'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta5'}{
							 \describe{
							\item{hyperid = }{`101125'}
							\item{name = }{`beta5'}
							\item{short.name = }{`beta5'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`normal'}
							\item{param = }{`0 1'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta6'}{
							 \describe{
							\item{hyperid = }{`101126'}
							\item{name = }{`beta6'}
							\item{short.name = }{`beta6'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`normal'}
							\item{param = }{`0 1'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta7'}{
							 \describe{
							\item{hyperid = }{`101127'}
							\item{name = }{`beta7'}
							\item{short.name = }{`beta7'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`normal'}
							\item{param = }{`0 1'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta8'}{
							 \describe{
							\item{hyperid = }{`101128'}
							\item{name = }{`beta8'}
							\item{short.name = }{`beta8'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`normal'}
							\item{param = }{`0 1'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta9'}{
							 \describe{
							\item{hyperid = }{`101129'}
							\item{name = }{`beta9'}
							\item{short.name = }{`beta9'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`normal'}
							\item{param = }{`0 1'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta10'}{
							 \describe{
							\item{hyperid = }{`101130'}
							\item{name = }{`beta10'}
							\item{short.name = }{`beta10'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`normal'}
							\item{param = }{`0 1'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta11'}{
							 \describe{
							\item{hyperid = }{`101131'}
							\item{name = }{`beta11'}
							\item{short.name = }{`beta11'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`normal'}
							\item{param = }{`0 1'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta12'}{
							 \describe{
							\item{hyperid = }{`101132'}
							\item{name = }{`beta12'}
							\item{short.name = }{`beta12'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`normal'}
							\item{param = }{`0 1'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta13'}{
							 \describe{
							\item{hyperid = }{`101133'}
							\item{name = }{`beta13'}
							\item{short.name = }{`beta13'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`normal'}
							\item{param = }{`0 1'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta14'}{
							 \describe{
							\item{hyperid = }{`101134'}
							\item{name = }{`beta14'}
							\item{short.name = }{`beta14'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`normal'}
							\item{param = }{`0 1'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta15'}{
							 \describe{
							\item{hyperid = }{`101135'}
							\item{name = }{`beta15'}
							\item{short.name = }{`beta15'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`normal'}
							\item{param = }{`0 1'}
							\item{to.theta = }{`\code{function(x) x}'}
							\item{from.theta = }{`\code{function(x) x}'}
							}
						}
						\item{Hyperparameter `theta16'}{
							 \describe{
							\item{hyperid = }{`101136'}
							\item{name = }{`overdispersion'}
							\item{short.name = }{`overdispersion'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`pc.gamma'}
							\item{param = }{`7'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{doc = }{`NegBinomial-Poisson mixture'}
								\item{status = }{`experimental'}
								\item{survival = }{`FALSE'}
								\item{discrete = }{`TRUE'}
								\item{link = }{`default logit probit'}
								\item{pdf = }{`nmixnb'}
							 }
						 }
					 }
				 }
				\item{Model `gp'.}{Number of hyperparmeters are 1.
					 \describe{
						\item{Hyperparameter `theta'}{
							 \describe{
							\item{hyperid = }{`101201'}
							\item{name = }{`shape'}
							\item{short.name = }{`xi'}
							\item{initial = }{`-2.30258509299405'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`loggamma'}
							\item{param = }{`1 15'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{doc = }{`Generalized Pareto likelihood'}
								\item{status = }{`experimental'}
								\item{survival = }{`FALSE'}
								\item{discrete = }{`TRUE'}
								\item{link = }{`default quantile'}
								\item{pdf = }{`genPareto'}
							 }
						 }
					 }
				 }
				\item{Model `logperiodogram'.}{Number of hyperparmeters are 0.
				 }
			 }
		 }
		\item{Section `prior'.}{
			 Valid models in this section are:
			 \describe{
				\item{Model `normal'.}{					 Number of parameters in the prior = 2 
				 }
				\item{Model `gaussian'.}{					 Number of parameters in the prior = 2 
				 }
				\item{Model `wishart1d'.}{					 Number of parameters in the prior = 2 
				 }
				\item{Model `wishart2d'.}{					 Number of parameters in the prior = 4 
				 }
				\item{Model `wishart3d'.}{					 Number of parameters in the prior = 7 
				 }
				\item{Model `wishart4d'.}{					 Number of parameters in the prior = 11 
				 }
				\item{Model `wishart5d'.}{					 Number of parameters in the prior = 16 
				 }
				\item{Model `loggamma'.}{					 Number of parameters in the prior = 2 
				 }
				\item{Model `minuslogsqrtruncnormal'.}{					 Number of parameters in the prior = 2 
				 }
				\item{Model `logtnormal'.}{					 Number of parameters in the prior = 2 
				 }
				\item{Model `logtgaussian'.}{					 Number of parameters in the prior = 2 
				 }
				\item{Model `flat'.}{					 Number of parameters in the prior = 0 
				 }
				\item{Model `logflat'.}{					 Number of parameters in the prior = 0 
				 }
				\item{Model `logiflat'.}{					 Number of parameters in the prior = 0 
				 }
				\item{Model `mvnorm'.}{					 Number of parameters in the prior = -1 
				 }
				\item{Model `pc.ar'.}{					 Number of parameters in the prior = 1 
				 }
				\item{Model `dirichlet'.}{					 Number of parameters in the prior = 1 
				 }
				\item{Model `none'.}{					 Number of parameters in the prior = 0 
				 }
				\item{Model `invalid'.}{					 Number of parameters in the prior = 0 
				 }
				\item{Model `betacorrelation'.}{					 Number of parameters in the prior = 2 
				 }
				\item{Model `logitbeta'.}{					 Number of parameters in the prior = 2 
				 }
				\item{Model `pc.prec'.}{					 Number of parameters in the prior = 2 
				 }
				\item{Model `pc.dof'.}{					 Number of parameters in the prior = 2 
				 }
				\item{Model `pc.cor0'.}{					 Number of parameters in the prior = 2 
				 }
				\item{Model `pc.cor1'.}{					 Number of parameters in the prior = 2 
				 }
				\item{Model `pc.fgnh'.}{					 Number of parameters in the prior = 2 
				 }
				\item{Model `pc.spde.GA'.}{					 Number of parameters in the prior = 4 
				 }
				\item{Model `pc.matern'.}{					 Number of parameters in the prior = 3 
				 }
				\item{Model `pc.range'.}{					 Number of parameters in the prior = 2 
				 }
				\item{Model `pc.gamma'.}{					 Number of parameters in the prior = 1 
				 }
				\item{Model `pc.mgamma'.}{					 Number of parameters in the prior = 1 
				 }
				\item{Model `pc.gammacount'.}{					 Number of parameters in the prior = 1 
				 }
				\item{Model `pc'.}{					 Number of parameters in the prior = 2 
				 }
				\item{Model `ref.ar'.}{					 Number of parameters in the prior = 0 
				 }
				\item{Model `pom'.}{					 Number of parameters in the prior = 0 
				 }
				\item{Model `jeffreystdf'.}{					 Number of parameters in the prior = 0 
				 }
				\item{Model `expression:'.}{					 Number of parameters in the prior = -1 
				 }
				\item{Model `table:'.}{					 Number of parameters in the prior = -1 
				 }
			 }
		 }
		\item{Section `wrapper'.}{
			 Valid models in this section are:
			 \describe{
				\item{Model `joint'.}{Number of hyperparmeters are 1.
					 \describe{
						\item{Hyperparameter `theta'}{
							 \describe{
							\item{hyperid = }{`102001'}
							\item{name = }{`log precision'}
							\item{short.name = }{`prec'}
							\item{initial = }{`0'}
							\item{fixed = }{`TRUE'}
							\item{prior = }{`loggamma'}
							\item{param = }{`1 5e-05'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{doc = }{`(experimental)'}
								\item{constr = }{`FALSE'}
								\item{nrow.ncol = }{`FALSE'}
								\item{augmented = }{`FALSE'}
								\item{aug.factor = }{`1'}
								\item{aug.constr = }{`NULL'}
								\item{n.div.by = }{`NULL'}
								\item{n.required = }{`FALSE'}
								\item{set.default.values = }{`FALSE'}
								\item{pdf = }{`NA'}
							 }
						 }
					 }
				 }
			 }
		 }
	 } 
} 
\examples{
## How to set hyperparameters to pass as the argument 'hyper'. This
## format is compatible with the old style (using 'initial', 'fixed',
## 'prior', 'param'), but the new style using 'hyper' take preceedence
## over the old style. The two styles can also be mixed. The old style
## might be removed from the code in the future...

## Only a subset need to be given
   hyper = list(theta = list(initial = 2))
## The `name' can be used instead of 'theta', or 'theta1', 'theta2',...
   hyper = list(precision = list(initial = 2))
   hyper = list(precision = list(prior = "flat", param = numeric(0)))
   hyper = list(theta2 = list(initial=3), theta1 = list(prior = "gaussian"))
## The 'short.name' can be used instead of 'name'
   hyper = list(rho = list(param = c(0,1)))
}
