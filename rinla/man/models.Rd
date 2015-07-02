
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
								\item{constr = }{`FALSE'}
								\item{nrow.ncol = }{`FALSE'}
								\item{augmented = }{`FALSE'}
								\item{aug.factor = }{`1'}
								\item{aug.constr = }{`NULL'}
								\item{n.div.by = }{`NULL'}
								\item{n.required = }{`FALSE'}
								\item{set.default.values = }{`FALSE'}
								\item{status = }{`experimental'}
								\item{pdf = }{`mec'}
							 }
						 }
					 }
				 }
				\item{Model `meb'.}{Number of hyperparmeters are 2.
					 \describe{
						\item{Hyperparameter `theta1'}{
							 \describe{
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
								\item{constr = }{`FALSE'}
								\item{nrow.ncol = }{`FALSE'}
								\item{augmented = }{`FALSE'}
								\item{aug.factor = }{`1'}
								\item{aug.constr = }{`NULL'}
								\item{n.div.by = }{`NULL'}
								\item{n.required = }{`FALSE'}
								\item{set.default.values = }{`FALSE'}
								\item{status = }{`experimental'}
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
								\item{constr = }{`TRUE'}
								\item{nrow.ncol = }{`FALSE'}
								\item{augmented = }{`FALSE'}
								\item{aug.factor = }{`1'}
								\item{aug.constr = }{`NULL'}
								\item{n.div.by = }{`NULL'}
								\item{n.required = }{`FALSE'}
								\item{set.default.values = }{`FALSE'}
								\item{pdf = }{`rw1'}
							 }
						 }
					 }
				 }
				\item{Model `rw2'.}{Number of hyperparmeters are 1.
					 \describe{
						\item{Hyperparameter `theta'}{
							 \describe{
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
								\item{constr = }{`TRUE'}
								\item{nrow.ncol = }{`FALSE'}
								\item{augmented = }{`FALSE'}
								\item{aug.factor = }{`1'}
								\item{aug.constr = }{`NULL'}
								\item{n.div.by = }{`NULL'}
								\item{n.required = }{`FALSE'}
								\item{set.default.values = }{`FALSE'}
								\item{pdf = }{`rw2'}
							 }
						 }
					 }
				 }
				\item{Model `crw2'.}{Number of hyperparmeters are 1.
					 \describe{
						\item{Hyperparameter `theta'}{
							 \describe{
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
								\item{constr = }{`TRUE'}
								\item{nrow.ncol = }{`FALSE'}
								\item{augmented = }{`FALSE'}
								\item{aug.factor = }{`2'}
								\item{aug.constr = }{`1'}
								\item{n.div.by = }{`NULL'}
								\item{n.required = }{`FALSE'}
								\item{set.default.values = }{`FALSE'}
								\item{pdf = }{`crw2'}
							 }
						 }
					 }
				 }
				\item{Model `seasonal'.}{Number of hyperparmeters are 1.
					 \describe{
						\item{Hyperparameter `theta'}{
							 \describe{
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
							\item{name = }{`logit phi'}
							\item{short.name = }{`phi'}
							\item{prior = }{`pc'}
							\item{param = }{`0.5 -1'}
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
				\item{Model `ar1'.}{Number of hyperparmeters are 2.
					 \describe{
						\item{Hyperparameter `theta1'}{
							 \describe{
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
					 }
					 \describe{
						\item{Properties:}{
							\describe{
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
				\item{Model `ar'.}{Number of hyperparmeters are 11.
					 \describe{
						\item{Hyperparameter `theta1'}{
							 \describe{
							\item{name = }{`log precision'}
							\item{short.name = }{`prec'}
							\item{initial = }{`4'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`pc.prec'}
							\item{param = }{`1 0.01'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
						\item{Hyperparameter `theta2'}{
							 \describe{
							\item{name = }{`pacf1'}
							\item{short.name = }{`pacf1'}
							\item{initial = }{`1'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`pc.ar'}
							\item{param = }{`1'}
							\item{to.theta = }{`\code{function(x) log((1+x)/(1-x))}'}
							\item{from.theta = }{`\code{function(x) 2*exp(x)/(1+exp(x))-1}'}
							}
						}
						\item{Hyperparameter `theta3'}{
							 \describe{
							\item{name = }{`pacf2'}
							\item{short.name = }{`pacf2'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) log((1+x)/(1-x))}'}
							\item{from.theta = }{`\code{function(x) 2*exp(x)/(1+exp(x))-1}'}
							}
						}
						\item{Hyperparameter `theta4'}{
							 \describe{
							\item{name = }{`pacf3'}
							\item{short.name = }{`pacf3'}
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
							\item{name = }{`pacf4'}
							\item{short.name = }{`pacf4'}
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
							\item{name = }{`pacf5'}
							\item{short.name = }{`pacf5'}
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
							\item{name = }{`pacf6'}
							\item{short.name = }{`pacf6'}
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
							\item{name = }{`pacf7'}
							\item{short.name = }{`pacf7'}
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
							\item{name = }{`pacf8'}
							\item{short.name = }{`pacf8'}
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
							\item{name = }{`pacf9'}
							\item{short.name = }{`pacf9'}
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
							\item{name = }{`pacf10'}
							\item{short.name = }{`pacf10'}
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
								\item{constr = }{`FALSE'}
								\item{nrow.ncol = }{`FALSE'}
								\item{augmented = }{`FALSE'}
								\item{aug.factor = }{`1'}
								\item{aug.constr = }{`NULL'}
								\item{n.div.by = }{`NULL'}
								\item{n.required = }{`FALSE'}
								\item{set.default.values = }{`FALSE'}
								\item{status = }{`experimental'}
								\item{pdf = }{`ar'}
							 }
						 }
					 }
				 }
				\item{Model `ou'.}{Number of hyperparmeters are 2.
					 \describe{
						\item{Hyperparameter `theta1'}{
							 \describe{
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
				\item{Model `generic'.}{Number of hyperparmeters are 1.
					 \describe{
						\item{Hyperparameter `theta'}{
							 \describe{
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
							\item{name = }{`logit phi'}
							\item{short.name = }{`phi'}
							\item{prior = }{`pc'}
							\item{param = }{`0.5 -1'}
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
				\item{Model `copy'.}{Number of hyperparmeters are 1.
					 \describe{
						\item{Hyperparameter `theta'}{
							 \describe{
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
			 }
		 }
		\item{Section `group'.}{
			 Valid models in this section are:
			 \describe{
				\item{Model `exchangeable'.}{Number of hyperparmeters are 1.
					 \describe{
						\item{Hyperparameter `theta'}{
							 \describe{
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
							 }
						 }
					 }
				 }
				\item{Model `ar1'.}{Number of hyperparmeters are 1.
					 \describe{
						\item{Hyperparameter `theta'}{
							 \describe{
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
							 }
						 }
					 }
				 }
				\item{Model `ar'.}{Number of hyperparmeters are 11.
					 \describe{
						\item{Hyperparameter `theta1'}{
							 \describe{
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
						\item{Hyperparameter `theta2'}{
							 \describe{
							\item{name = }{`pacf1'}
							\item{short.name = }{`pacf1'}
							\item{initial = }{`2'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`mvnorm'}
							\item{param = }{`0 0.15'}
							\item{to.theta = }{`\code{function(x) log((1+x)/(1-x))}'}
							\item{from.theta = }{`\code{function(x) 2*exp(x)/(1+exp(x))-1}'}
							}
						}
						\item{Hyperparameter `theta3'}{
							 \describe{
							\item{name = }{`pacf2'}
							\item{short.name = }{`pacf2'}
							\item{initial = }{`0'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`none'}
							\item{param = }{`'}
							\item{to.theta = }{`\code{function(x) log((1+x)/(1-x))}'}
							\item{from.theta = }{`\code{function(x) 2*exp(x)/(1+exp(x))-1}'}
							}
						}
						\item{Hyperparameter `theta4'}{
							 \describe{
							\item{name = }{`pacf3'}
							\item{short.name = }{`pacf3'}
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
							\item{name = }{`pacf4'}
							\item{short.name = }{`pacf4'}
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
							\item{name = }{`pacf5'}
							\item{short.name = }{`pacf5'}
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
							\item{name = }{`pacf6'}
							\item{short.name = }{`pacf6'}
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
							\item{name = }{`pacf7'}
							\item{short.name = }{`pacf7'}
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
							\item{name = }{`pacf8'}
							\item{short.name = }{`pacf8'}
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
							\item{name = }{`pacf9'}
							\item{short.name = }{`pacf9'}
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
							\item{name = }{`pacf10'}
							\item{short.name = }{`pacf10'}
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
							 }
						 }
					 }
				 }
				\item{Model `rw1'.}{Number of hyperparmeters are 1.
					 \describe{
						\item{Hyperparameter `theta'}{
							 \describe{
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
							 }
						 }
					 }
				 }
				\item{Model `rw2'.}{Number of hyperparmeters are 1.
					 \describe{
						\item{Hyperparameter `theta'}{
							 \describe{
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
							 }
						 }
					 }
				 }
				\item{Model `besag'.}{Number of hyperparmeters are 1.
					 \describe{
						\item{Hyperparameter `theta'}{
							 \describe{
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
							 }
						 }
					 }
				 }
				\item{Model `I'.}{Number of hyperparmeters are 1.
					 \describe{
						\item{Hyperparameter `theta'}{
							 \describe{
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
							\item{name = }{`log precision'}
							\item{short.name = }{`prec'}
							\item{prior = }{`loggamma'}
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
				\item{Model `identity'.}{Number of hyperparmeters are 0.
				 }
				\item{Model `log'.}{Number of hyperparmeters are 0.
				 }
				\item{Model `logit'.}{Number of hyperparmeters are 0.
				 }
				\item{Model `probit'.}{Number of hyperparmeters are 0.
				 }
				\item{Model `tan'.}{Number of hyperparmeters are 0.
				 }
				\item{Model `sslogit'.}{Number of hyperparmeters are 2.
					 \describe{
						\item{Hyperparameter `theta1'}{
							 \describe{
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
								\item{pdf = }{`NA'}
							 }
						 }
					 }
				 }
				\item{Model `logoffset'.}{Number of hyperparmeters are 1.
					 \describe{
						\item{Hyperparameter `theta'}{
							 \describe{
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
								\item{pdf = }{`logoffset'}
							 }
						 }
					 }
				 }
				\item{Model `test1'.}{Number of hyperparmeters are 1.
					 \describe{
						\item{Hyperparameter `theta'}{
							 \describe{
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
								\item{pdf = }{`NA'}
							 }
						 }
					 }
				 }
				\item{Model `special2'.}{Number of hyperparmeters are 1.
					 \describe{
						\item{Hyperparameter `theta'}{
							 \describe{
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
								\item{pdf = }{`NA'}
							 }
						 }
					 }
				 }
				\item{Model `special1'.}{Number of hyperparmeters are 11.
					 \describe{
						\item{Hyperparameter `theta1'}{
							 \describe{
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
							\item{name = }{`log precision'}
							\item{short.name = }{`prec'}
							\item{initial = }{`11'}
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
							 }
						 }
					 }
				 }
				\item{Model `rw2'.}{Number of hyperparmeters are 1.
					 \describe{
						\item{Hyperparameter `theta'}{
							 \describe{
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
				\item{Model `gpoisson'.}{Number of hyperparmeters are 2.
					 \describe{
						\item{Hyperparameter `theta1'}{
							 \describe{
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
				\item{Model `testbinomial1'.}{Number of hyperparmeters are 2.
					 \describe{
						\item{Hyperparameter `theta1'}{
							 \describe{
							\item{name = }{`sensitivity'}
							\item{short.name = }{`s'}
							\item{initial = }{`3'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`logitbeta'}
							\item{param = }{`2 1'}
							\item{to.theta = }{`\code{function(x) log(x/(1-x))}'}
							\item{from.theta = }{`\code{function(x) exp(x)/(1+exp(x))}'}
							}
						}
						\item{Hyperparameter `theta2'}{
							 \describe{
							\item{name = }{`specificity'}
							\item{short.name = }{`e'}
							\item{initial = }{`3'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`logitbeta'}
							\item{param = }{`2 1'}
							\item{to.theta = }{`\code{function(x) log(x/(1-x))}'}
							\item{from.theta = }{`\code{function(x) exp(x)/(1+exp(x))}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{status = }{`experimental'}
								\item{survival = }{`FALSE'}
								\item{discrete = }{`TRUE'}
								\item{link = }{`default logit probit cloglog log'}
								\item{pdf = }{`testbinomial1'}
							 }
						 }
					 }
				 }
				\item{Model `gamma'.}{Number of hyperparmeters are 1.
					 \describe{
						\item{Hyperparameter `theta'}{
							 \describe{
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
								\item{survival = }{`FALSE'}
								\item{discrete = }{`FALSE'}
								\item{link = }{`default log'}
								\item{pdf = }{`gamma'}
							 }
						 }
					 }
				 }
				\item{Model `gammacount'.}{Number of hyperparmeters are 1.
					 \describe{
						\item{Hyperparameter `theta'}{
							 \describe{
							\item{name = }{`log alpha'}
							\item{short.name = }{`alpha'}
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
								\item{survival = }{`FALSE'}
								\item{discrete = }{`FALSE'}
								\item{link = }{`default log'}
								\item{status = }{`experimental'}
								\item{pdf = }{`gammacount'}
							 }
						 }
					 }
				 }
				\item{Model `beta'.}{Number of hyperparmeters are 1.
					 \describe{
						\item{Hyperparameter `theta'}{
							 \describe{
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
								\item{survival = }{`FALSE'}
								\item{discrete = }{`FALSE'}
								\item{link = }{`default logit probit cloglog'}
								\item{pdf = }{`beta'}
							 }
						 }
					 }
				 }
				\item{Model `betabinomial'.}{Number of hyperparmeters are 1.
					 \describe{
						\item{Hyperparameter `theta'}{
							 \describe{
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
								\item{survival = }{`FALSE'}
								\item{discrete = }{`TRUE'}
								\item{link = }{`default logit probit cloglog'}
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
							\item{name = }{`size'}
							\item{short.name = }{`size'}
							\item{initial = }{`2.30258509299405'}
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
								\item{survival = }{`FALSE'}
								\item{discrete = }{`TRUE'}
								\item{link = }{`default log logoffset'}
								\item{pdf = }{`nbinomial'}
							 }
						 }
					 }
				 }
				\item{Model `simplex'.}{Number of hyperparmeters are 1.
					 \describe{
						\item{Hyperparameter `theta'}{
							 \describe{
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
								\item{survival = }{`FALSE'}
								\item{discrete = }{`FALSE'}
								\item{link = }{`default logit probit cloglog'}
								\item{pdf = }{`simplex'}
							 }
						 }
					 }
				 }
				\item{Model `gaussian'.}{Number of hyperparmeters are 1.
					 \describe{
						\item{Hyperparameter `theta'}{
							 \describe{
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
								\item{survival = }{`FALSE'}
								\item{discrete = }{`FALSE'}
								\item{link = }{`default identity logit log'}
								\item{pdf = }{`gaussian'}
							 }
						 }
					 }
				 }
				\item{Model `normal'.}{Number of hyperparmeters are 1.
					 \describe{
						\item{Hyperparameter `theta'}{
							 \describe{
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
								\item{survival = }{`FALSE'}
								\item{discrete = }{`FALSE'}
								\item{link = }{`default identity'}
								\item{pdf = }{`gaussian'}
							 }
						 }
					 }
				 }
				\item{Model `circularnormal'.}{Number of hyperparmeters are 1.
					 \describe{
						\item{Hyperparameter `theta'}{
							 \describe{
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
								\item{survival = }{`FALSE'}
								\item{discrete = }{`FALSE'}
								\item{link = }{`default identity'}
								\item{pdf = }{`loggammafrailty'}
							 }
						 }
					 }
				 }
				\item{Model `logistic'.}{Number of hyperparmeters are 1.
					 \describe{
						\item{Hyperparameter `theta'}{
							 \describe{
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
								\item{survival = }{`FALSE'}
								\item{discrete = }{`FALSE'}
								\item{link = }{`default identity'}
								\item{status = }{`experimental'}
								\item{pdf = }{`gev'}
							 }
						 }
					 }
				 }
				\item{Model `laplace'.}{Number of hyperparmeters are 1.
					 \describe{
						\item{Hyperparameter `theta'}{
							 \describe{
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
								\item{survival = }{`FALSE'}
								\item{discrete = }{`FALSE'}
								\item{link = }{`default identity'}
								\item{status = }{`disabled'}
								\item{pdf = }{`laplace'}
							 }
						 }
					 }
				 }
				\item{Model `lognormal'.}{Number of hyperparmeters are 1.
					 \describe{
						\item{Hyperparameter `theta'}{
							 \describe{
							\item{name = }{`log precision'}
							\item{short.name = }{`prec'}
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
				\item{Model `coxph'.}{Number of hyperparmeters are 0.
				 }
				\item{Model `weibull'.}{Number of hyperparmeters are 1.
					 \describe{
						\item{Hyperparameter `theta'}{
							 \describe{
							\item{name = }{`log alpha'}
							\item{short.name = }{`a'}
							\item{initial = }{`0'}
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
								\item{survival = }{`TRUE'}
								\item{discrete = }{`FALSE'}
								\item{link = }{`default log'}
								\item{pdf = }{`weibull'}
							 }
						 }
					 }
				 }
				\item{Model `loglogistic'.}{Number of hyperparmeters are 1.
					 \describe{
						\item{Hyperparameter `theta'}{
							 \describe{
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
								\item{survival = }{`TRUE'}
								\item{discrete = }{`FALSE'}
								\item{link = }{`default log'}
								\item{pdf = }{`loglogistic'}
							 }
						 }
					 }
				 }
				\item{Model `weibullcure'.}{Number of hyperparmeters are 2.
					 \describe{
						\item{Hyperparameter `theta1'}{
							 \describe{
							\item{name = }{`log alpha'}
							\item{short.name = }{`a'}
							\item{initial = }{`4'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`loggamma'}
							\item{param = }{`25 25'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
						\item{Hyperparameter `theta2'}{
							 \describe{
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
								\item{survival = }{`TRUE'}
								\item{discrete = }{`FALSE'}
								\item{link = }{`default log'}
								\item{pdf = }{`NA'}
							 }
						 }
					 }
				 }
				\item{Model `stochvol'.}{Number of hyperparmeters are 1.
					 \describe{
						\item{Hyperparameter `theta'}{
							 \describe{
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
							\item{name = }{`log degrees of freedom'}
							\item{short.name = }{`dof'}
							\item{initial = }{`4'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`loggamma'}
							\item{param = }{`1 0.5'}
							\item{to.theta = }{`\code{function(x) log(x-2)}'}
							\item{from.theta = }{`\code{function(x) 2+exp(x)}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
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
								\item{survival = }{`FALSE'}
								\item{discrete = }{`TRUE'}
								\item{link = }{`default logit probit cloglog'}
								\item{pdf = }{`zeroinflated'}
							 }
						 }
					 }
				 }
				\item{Model `zeroinflatedbetabinomial1'.}{Number of hyperparmeters are 2.
					 \describe{
						\item{Hyperparameter `theta1'}{
							 \describe{
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
								\item{survival = }{`FALSE'}
								\item{discrete = }{`TRUE'}
								\item{link = }{`default logit probit cloglog'}
								\item{pdf = }{`zeroinflated'}
							 }
						 }
					 }
				 }
				\item{Model `zeroinflatedbinomial0'.}{Number of hyperparmeters are 1.
					 \describe{
						\item{Hyperparameter `theta'}{
							 \describe{
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
								\item{survival = }{`FALSE'}
								\item{discrete = }{`FALSE'}
								\item{link = }{`default logit probit cloglog'}
								\item{pdf = }{`zeroinflated'}
							 }
						 }
					 }
				 }
				\item{Model `zeroinflatedbinomial1'.}{Number of hyperparmeters are 1.
					 \describe{
						\item{Hyperparameter `theta'}{
							 \describe{
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
								\item{survival = }{`FALSE'}
								\item{discrete = }{`FALSE'}
								\item{link = }{`default logit probit cloglog'}
								\item{pdf = }{`zeroinflated'}
							 }
						 }
					 }
				 }
				\item{Model `zeroinflatedbinomial2'.}{Number of hyperparmeters are 1.
					 \describe{
						\item{Hyperparameter `theta'}{
							 \describe{
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
								\item{survival = }{`FALSE'}
								\item{discrete = }{`FALSE'}
								\item{link = }{`default logit probit cloglog'}
								\item{pdf = }{`zeroinflated'}
							 }
						 }
					 }
				 }
				\item{Model `zeroninflatedbinomial2'.}{Number of hyperparmeters are 2.
					 \describe{
						\item{Hyperparameter `theta1'}{
							 \describe{
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
								\item{survival = }{`FALSE'}
								\item{discrete = }{`FALSE'}
								\item{link = }{`default logit probit cloglog'}
								\item{pdf = }{`NA'}
							 }
						 }
					 }
				 }
				\item{Model `zeroinflatedbetabinomial2'.}{Number of hyperparmeters are 2.
					 \describe{
						\item{Hyperparameter `theta1'}{
							 \describe{
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
								\item{survival = }{`FALSE'}
								\item{discrete = }{`FALSE'}
								\item{link = }{`default logit probit cloglog'}
								\item{pdf = }{`zeroinflated'}
							 }
						 }
					 }
				 }
				\item{Model `zeroinflatednbinomial0'.}{Number of hyperparmeters are 2.
					 \describe{
						\item{Hyperparameter `theta1'}{
							 \describe{
							\item{name = }{`log size'}
							\item{short.name = }{`size'}
							\item{initial = }{`2.30258509299405'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`loggamma'}
							\item{param = }{`1 1'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
						\item{Hyperparameter `theta2'}{
							 \describe{
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
							\item{name = }{`log size'}
							\item{short.name = }{`size'}
							\item{initial = }{`2.30258509299405'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`loggamma'}
							\item{param = }{`1 1'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
						\item{Hyperparameter `theta2'}{
							 \describe{
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
								\item{survival = }{`FALSE'}
								\item{discrete = }{`FALSE'}
								\item{link = }{`default log'}
								\item{pdf = }{`zeroinflated'}
							 }
						 }
					 }
				 }
				\item{Model `zeroinflatednbinomial1strata2'.}{Number of hyperparmeters are 3.
					 \describe{
						\item{Hyperparameter `theta1'}{
							 \describe{
							\item{name = }{`log size'}
							\item{short.name = }{`size'}
							\item{initial = }{`2.30258509299405'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`loggamma'}
							\item{param = }{`1 1'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
						\item{Hyperparameter `theta2'}{
							 \describe{
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
					 }
					 \describe{
						\item{Properties:}{
							\describe{
								\item{status = }{`experimental'}
								\item{survival = }{`FALSE'}
								\item{discrete = }{`FALSE'}
								\item{link = }{`default log'}
								\item{pdf = }{`zeroinflated'}
							 }
						 }
					 }
				 }
				\item{Model `zeroinflatednbinomial1strata3'.}{Number of hyperparmeters are 3.
					 \describe{
						\item{Hyperparameter `theta1'}{
							 \describe{
							\item{name = }{`log size 1'}
							\item{short.name = }{`size1'}
							\item{initial = }{`2.30258509299405'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`loggamma'}
							\item{param = }{`1 1'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
						\item{Hyperparameter `theta2'}{
							 \describe{
							\item{name = }{`log size 2'}
							\item{short.name = }{`size2'}
							\item{initial = }{`2.30258509299405'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`loggamma'}
							\item{param = }{`1 1'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
						\item{Hyperparameter `theta3'}{
							 \describe{
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
							\item{name = }{`log size'}
							\item{short.name = }{`size'}
							\item{initial = }{`2.30258509299405'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`loggamma'}
							\item{param = }{`1 1'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x)}'}
							}
						}
						\item{Hyperparameter `theta2'}{
							 \describe{
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
							\item{name = }{`log degrees of freedom'}
							\item{short.name = }{`dof'}
							\item{initial = }{`5'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`loggamma'}
							\item{param = }{`1 0.5'}
							\item{to.theta = }{`\code{function(x) log(x-2)}'}
							\item{from.theta = }{`\code{function(x) 2+exp(x)}'}
							}
						}
					 }
					 \describe{
						\item{Properties:}{
							\describe{
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
							\item{name = }{`log degrees of freedom'}
							\item{short.name = }{`dof'}
							\item{initial = }{`4'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`loggamma'}
							\item{param = }{`1 0.01'}
							\item{to.theta = }{`\code{function(x) log(x-5)}'}
							\item{from.theta = }{`\code{function(x) 5+exp(x)}'}
							}
						}
						\item{Hyperparameter `theta2'}{
							 \describe{
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
								\item{survival = }{`FALSE'}
								\item{discrete = }{`FALSE'}
								\item{link = }{`default identity'}
								\item{pdf = }{`tstrata'}
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
				\item{Model `none'.}{					 Number of parameters in the prior = 0 
				 }
				\item{Model `betacorrelation'.}{					 Number of parameters in the prior = 2 
				 }
				\item{Model `logitbeta'.}{					 Number of parameters in the prior = 2 
				 }
				\item{Model `pc.prec'.}{					 Number of parameters in the prior = 2 
				 }
				\item{Model `pc.dof'.}{					 Number of parameters in the prior = 2 
				 }
				\item{Model `pc.rho0'.}{					 Number of parameters in the prior = 2 
				 }
				\item{Model `pc.rho1'.}{					 Number of parameters in the prior = 2 
				 }
				\item{Model `pc.spde.GA'.}{					 Number of parameters in the prior = 4 
				 }
				\item{Model `pc'.}{					 Number of parameters in the prior = 2 
				 }
				\item{Model `ref.ar'.}{					 Number of parameters in the prior = 0 
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
