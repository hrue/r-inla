
\name{inla.models}
\alias{inla.models}
\title{Valid models in INLA}
\description{
This page describe the models implemented in \code{inla}, divided into sections: latent, group, predictor, hazard, likelihood, prior, wrapper .
}

\usage{
inla.models()
}
\value{
 	 Valid sections are: latent, group, predictor, hazard, likelihood, prior, wrapper 

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
							\item{param = }{`1 0.001'}
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
							\item{prior = }{`normal'}
							\item{param = }{`0 5e-05'}
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
							\item{param = }{`0 0.001'}
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
								\item{n.required = }{`TRUE'}
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
							\item{from.theta = }{`\code{function(x) exp(x)                                            }'}
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
							\item{from.theta = }{`\code{function(x) exp(x)                                            }'}
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
							\item{from.theta = }{`\code{function(x) exp(x)                                            }'}
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
								\item{n.required = }{`FALSE'}
								\item{set.default.values = }{`FALSE'}
								\item{pdf = }{`NA'}
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
								\item{constr = }{`TRUE'}
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
							\item{to.theta = }{`\code{function(x, low = -Inf, high = Inf) {} \code{                                                    if (low == -Inf && high == Inf) {} \code{                                                        return (x)} \code{                                                    } else {} \code{                                                        stopifnot((low != -Inf) &&  (high != Inf) && (low < high))} \code{                                                        return (log( - (low - x)/(high -x))) } \code{                                                    }} \code{                                                }}'}
							\item{from.theta = }{`\code{function(x, low = -Inf, high = Inf) {} \code{                                                    if (low == -Inf && high == Inf) {} \code{                                                        return (x)} \code{                                                    } else {} \code{                                                        stopifnot((low != -Inf) &&  (high != Inf) && (low < high))} \code{                                                        return (low + exp(x)/(1+exp(x)) * (high - low))} \code{                                                    }} \code{                                                }}'}
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
							\item{to.theta = }{`\code{function(x, ngroup=2) log((1+x*(ngroup-1))/(1-x))}'}
							\item{from.theta = }{`\code{function(x, ngroup=2) (exp(x)-1)/(exp(x) + ngroup -1)}'}
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
				\item{Model `binomial'.}{Number of hyperparmeters are 0.
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
							\item{param = }{`1 100'}
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
								\item{link = }{`default log'}
								\item{pdf = }{`nbinomial'}
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
								\item{link = }{`default identity'}
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
							\item{initial = }{`4'}
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
							\item{param = }{`0 6.25'}
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
								\item{pdf = }{`laplace'}
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
				\item{Model `stochvol'.}{Number of hyperparmeters are 0.
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
							\item{from.theta = }{`\code{function(x) exp(x) }'}
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
							\item{initial = }{`3'}
							\item{fixed = }{`FALSE'}
							\item{prior = }{`loggamma'}
							\item{param = }{`1 5e-05'}
							\item{to.theta = }{`\code{function(x) log(x)}'}
							\item{from.theta = }{`\code{function(x) exp(x) }'}
							}
						}
						\item{Hyperparameter `theta2'}{
							 \describe{
							\item{name = }{`log degrees of freedom'}
							\item{short.name = }{`dof'}
							\item{initial = }{`0'}
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
								\item{pdf = }{`Student-t'}
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
				\item{Model `loggamma'.}{					 Number of parameters in the prior = 2 
				 }
				\item{Model `minuslogsqrtruncnormal'.}{					 Number of parameters in the prior = 1 
				 }
				\item{Model `logtnormal'.}{					 Number of parameters in the prior = 1 
				 }
				\item{Model `logtgaussian'.}{					 Number of parameters in the prior = 1 
				 }
				\item{Model `flat'.}{					 Number of parameters in the prior = 0 
				 }
				\item{Model `logflat'.}{					 Number of parameters in the prior = 0 
				 }
				\item{Model `logiflat'.}{					 Number of parameters in the prior = 0 
				 }
				\item{Model `none'.}{					 Number of parameters in the prior = 0 
				 }
				\item{Model `betacorrelation'.}{					 Number of parameters in the prior = 2 
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
