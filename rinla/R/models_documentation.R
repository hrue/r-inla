## DO NOT EDIT!
## This file is generated automatically from models.R
#' @title Valid models in INLA
#'
#' @name inla.models
#' @rdname models
#' @aliases inla.models
#'
#' @description
#' This page describe the models implemented in `inla`, divided into sections:
#' latent, group, scopy, mix, link, predictor, hazard, likelihood, prior, wrapper, lp.scale.
#'
#' @usage
#' inla.models()
#'
#' @return
#' Valid sections are:
#' latent, group, scopy, mix, link, predictor, hazard, likelihood, prior, wrapper, lp.scale.
#' @section 'latent':
#'
#'    Valid models in this section are:
#'      \describe{
#'       \item{Model 'linear'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'Alternative interface to an fixed effect'}
#'               \item{constr = }{'FALSE'}
#'               \item{nrow.ncol = }{'FALSE'}
#'               \item{augmented = }{'FALSE'}
#'               \item{aug.factor = }{'1'}
#'               \item{aug.constr = }{'NULL'}
#'               \item{n.div.by = }{'NULL'}
#'               \item{n.required = }{'FALSE'}
#'               \item{set.default.values = }{'FALSE'}
#'               \item{pdf = }{'linear'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 0.
#'        }
#'       \item{Model 'iid'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'Gaussian random effects in dim=1'}
#'               \item{constr = }{'FALSE'}
#'               \item{nrow.ncol = }{'FALSE'}
#'               \item{augmented = }{'FALSE'}
#'               \item{aug.factor = }{'1'}
#'               \item{aug.constr = }{'NULL'}
#'               \item{n.div.by = }{'NULL'}
#'               \item{n.required = }{'FALSE'}
#'               \item{set.default.values = }{'FALSE'}
#'               \item{pdf = }{'indep'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 1.
#'          \describe{
#'           \item{Hyperparameter 'theta'}{
#'              \describe{
#'             \item{hyperid = }{1001}
#'             \item{name = }{log precision}
#'             \item{short.name = }{prec}
#'             \item{prior = }{loggamma}
#'             \item{param = }{1 5e-05}
#'             \item{initial = }{4}
#'             \item{fixed = }{FALSE}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'mec'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'Classical measurement error model'}
#'               \item{constr = }{'FALSE'}
#'               \item{nrow.ncol = }{'FALSE'}
#'               \item{augmented = }{'FALSE'}
#'               \item{aug.factor = }{'1'}
#'               \item{aug.constr = }{'NULL'}
#'               \item{n.div.by = }{'NULL'}
#'               \item{n.required = }{'FALSE'}
#'               \item{set.default.values = }{'FALSE'}
#'               \item{pdf = }{'mec'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 4.
#'          \describe{
#'           \item{Hyperparameter 'theta1'}{
#'              \describe{
#'             \item{hyperid = }{2001}
#'             \item{name = }{beta}
#'             \item{short.name = }{b}
#'             \item{prior = }{gaussian}
#'             \item{param = }{1 0.001}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta2'}{
#'              \describe{
#'             \item{hyperid = }{2002}
#'             \item{name = }{prec.u}
#'             \item{short.name = }{prec}
#'             \item{prior = }{loggamma}
#'             \item{param = }{1 1e-04}
#'             \item{initial = }{9.21034037197618}
#'             \item{fixed = }{TRUE}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta3'}{
#'              \describe{
#'             \item{hyperid = }{2003}
#'             \item{name = }{mean.x}
#'             \item{short.name = }{mu.x}
#'             \item{prior = }{gaussian}
#'             \item{param = }{0 1e-04}
#'             \item{initial = }{0}
#'             \item{fixed = }{TRUE}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta4'}{
#'              \describe{
#'             \item{hyperid = }{2004}
#'             \item{name = }{prec.x}
#'             \item{short.name = }{prec.x}
#'             \item{prior = }{loggamma}
#'             \item{param = }{1 10000}
#'             \item{initial = }{-9.21034037197618}
#'             \item{fixed = }{TRUE}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'meb'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'Berkson measurement error model'}
#'               \item{constr = }{'FALSE'}
#'               \item{nrow.ncol = }{'FALSE'}
#'               \item{augmented = }{'FALSE'}
#'               \item{aug.factor = }{'1'}
#'               \item{aug.constr = }{'NULL'}
#'               \item{n.div.by = }{'NULL'}
#'               \item{n.required = }{'FALSE'}
#'               \item{set.default.values = }{'FALSE'}
#'               \item{pdf = }{'meb'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 2.
#'          \describe{
#'           \item{Hyperparameter 'theta1'}{
#'              \describe{
#'             \item{hyperid = }{3001}
#'             \item{name = }{beta}
#'             \item{short.name = }{b}
#'             \item{prior = }{gaussian}
#'             \item{param = }{1 0.001}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta2'}{
#'              \describe{
#'             \item{hyperid = }{3002}
#'             \item{name = }{prec.u}
#'             \item{short.name = }{prec}
#'             \item{prior = }{loggamma}
#'             \item{param = }{1 1e-04}
#'             \item{initial = }{6.90775527898214}
#'             \item{fixed = }{FALSE}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'rgeneric'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'Generic latent model specified using R'}
#'               \item{constr = }{'FALSE'}
#'               \item{nrow.ncol = }{'FALSE'}
#'               \item{augmented = }{'FALSE'}
#'               \item{aug.factor = }{'1'}
#'               \item{aug.constr = }{'NULL'}
#'               \item{n.div.by = }{'NULL'}
#'               \item{n.required = }{'TRUE'}
#'               \item{set.default.values = }{'TRUE'}
#'               \item{status = }{'experimental'}
#'               \item{pdf = }{'rgeneric'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 0.
#'        }
#'       \item{Model 'cgeneric'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'Generic latent model specified using C'}
#'               \item{constr = }{'FALSE'}
#'               \item{nrow.ncol = }{'FALSE'}
#'               \item{augmented = }{'FALSE'}
#'               \item{aug.factor = }{'1'}
#'               \item{aug.constr = }{'NULL'}
#'               \item{n.div.by = }{'NULL'}
#'               \item{n.required = }{'TRUE'}
#'               \item{set.default.values = }{'TRUE'}
#'               \item{status = }{'experimental'}
#'               \item{pdf = }{'rgeneric'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 0.
#'        }
#'       \item{Model 'rw1'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'Random walk of order 1'}
#'               \item{constr = }{'TRUE'}
#'               \item{nrow.ncol = }{'FALSE'}
#'               \item{augmented = }{'FALSE'}
#'               \item{aug.factor = }{'1'}
#'               \item{aug.constr = }{'NULL'}
#'               \item{n.div.by = }{'NULL'}
#'               \item{n.required = }{'FALSE'}
#'               \item{set.default.values = }{'FALSE'}
#'               \item{min.diff = }{'1e-05'}
#'               \item{pdf = }{'rw1'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 1.
#'          \describe{
#'           \item{Hyperparameter 'theta'}{
#'              \describe{
#'             \item{hyperid = }{4001}
#'             \item{name = }{log precision}
#'             \item{short.name = }{prec}
#'             \item{prior = }{loggamma}
#'             \item{param = }{1 5e-05}
#'             \item{initial = }{4}
#'             \item{fixed = }{FALSE}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'rw2'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'Random walk of order 2'}
#'               \item{constr = }{'TRUE'}
#'               \item{nrow.ncol = }{'FALSE'}
#'               \item{augmented = }{'FALSE'}
#'               \item{aug.factor = }{'1'}
#'               \item{aug.constr = }{'NULL'}
#'               \item{n.div.by = }{'NULL'}
#'               \item{n.required = }{'FALSE'}
#'               \item{set.default.values = }{'FALSE'}
#'               \item{min.diff = }{'0.001'}
#'               \item{pdf = }{'rw2'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 1.
#'          \describe{
#'           \item{Hyperparameter 'theta'}{
#'              \describe{
#'             \item{hyperid = }{5001}
#'             \item{name = }{log precision}
#'             \item{short.name = }{prec}
#'             \item{prior = }{loggamma}
#'             \item{param = }{1 5e-05}
#'             \item{initial = }{4}
#'             \item{fixed = }{FALSE}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'crw2'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'Exact solution to the random walk of order 2'}
#'               \item{constr = }{'TRUE'}
#'               \item{nrow.ncol = }{'FALSE'}
#'               \item{augmented = }{'FALSE'}
#'               \item{aug.factor = }{'2'}
#'               \item{aug.constr = }{'1'}
#'               \item{n.div.by = }{'NULL'}
#'               \item{n.required = }{'FALSE'}
#'               \item{set.default.values = }{'FALSE'}
#'               \item{min.diff = }{'0.001'}
#'               \item{pdf = }{'crw2'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 1.
#'          \describe{
#'           \item{Hyperparameter 'theta'}{
#'              \describe{
#'             \item{hyperid = }{6001}
#'             \item{name = }{log precision}
#'             \item{short.name = }{prec}
#'             \item{prior = }{loggamma}
#'             \item{param = }{1 5e-05}
#'             \item{initial = }{4}
#'             \item{fixed = }{FALSE}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'seasonal'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'Seasonal model for time series'}
#'               \item{constr = }{'FALSE'}
#'               \item{nrow.ncol = }{'FALSE'}
#'               \item{augmented = }{'FALSE'}
#'               \item{aug.factor = }{'1'}
#'               \item{aug.constr = }{'NULL'}
#'               \item{n.div.by = }{'NULL'}
#'               \item{n.required = }{'FALSE'}
#'               \item{set.default.values = }{'FALSE'}
#'               \item{pdf = }{'seasonal'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 1.
#'          \describe{
#'           \item{Hyperparameter 'theta'}{
#'              \describe{
#'             \item{hyperid = }{7001}
#'             \item{name = }{log precision}
#'             \item{short.name = }{prec}
#'             \item{prior = }{loggamma}
#'             \item{param = }{1 5e-05}
#'             \item{initial = }{4}
#'             \item{fixed = }{FALSE}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'besag'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'The Besag area model (CAR-model)'}
#'               \item{constr = }{'TRUE'}
#'               \item{nrow.ncol = }{'FALSE'}
#'               \item{augmented = }{'FALSE'}
#'               \item{aug.factor = }{'1'}
#'               \item{aug.constr = }{'NULL'}
#'               \item{n.div.by = }{'NULL'}
#'               \item{n.required = }{'TRUE'}
#'               \item{set.default.values = }{'TRUE'}
#'               \item{pdf = }{'besag'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 1.
#'          \describe{
#'           \item{Hyperparameter 'theta'}{
#'              \describe{
#'             \item{hyperid = }{8001}
#'             \item{name = }{log precision}
#'             \item{short.name = }{prec}
#'             \item{prior = }{loggamma}
#'             \item{param = }{1 5e-05}
#'             \item{initial = }{4}
#'             \item{fixed = }{FALSE}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'besag2'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'The shared Besag model'}
#'               \item{constr = }{'TRUE'}
#'               \item{nrow.ncol = }{'FALSE'}
#'               \item{augmented = }{'FALSE'}
#'               \item{aug.factor = }{'1'}
#'               \item{aug.constr = }{'1 2'}
#'               \item{n.div.by = }{'2'}
#'               \item{n.required = }{'TRUE'}
#'               \item{set.default.values = }{'TRUE'}
#'               \item{pdf = }{'besag2'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 2.
#'          \describe{
#'           \item{Hyperparameter 'theta1'}{
#'              \describe{
#'             \item{hyperid = }{9001}
#'             \item{name = }{log precision}
#'             \item{short.name = }{prec}
#'             \item{prior = }{loggamma}
#'             \item{param = }{1 5e-05}
#'             \item{initial = }{4}
#'             \item{fixed = }{FALSE}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta2'}{
#'              \describe{
#'             \item{hyperid = }{9002}
#'             \item{name = }{scaling parameter}
#'             \item{short.name = }{a}
#'             \item{prior = }{loggamma}
#'             \item{param = }{10 10}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'bym'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'The BYM-model (Besag-York-Mollier model)'}
#'               \item{constr = }{'TRUE'}
#'               \item{nrow.ncol = }{'FALSE'}
#'               \item{augmented = }{'TRUE'}
#'               \item{aug.factor = }{'2'}
#'               \item{aug.constr = }{'2'}
#'               \item{n.div.by = }{'NULL'}
#'               \item{n.required = }{'TRUE'}
#'               \item{set.default.values = }{'TRUE'}
#'               \item{pdf = }{'bym'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 2.
#'          \describe{
#'           \item{Hyperparameter 'theta1'}{
#'              \describe{
#'             \item{hyperid = }{10001}
#'             \item{name = }{log unstructured precision}
#'             \item{short.name = }{prec.unstruct}
#'             \item{prior = }{loggamma}
#'             \item{param = }{1 5e-04}
#'             \item{initial = }{4}
#'             \item{fixed = }{FALSE}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta2'}{
#'              \describe{
#'             \item{hyperid = }{10002}
#'             \item{name = }{log spatial precision}
#'             \item{short.name = }{prec.spatial}
#'             \item{prior = }{loggamma}
#'             \item{param = }{1 5e-04}
#'             \item{initial = }{4}
#'             \item{fixed = }{FALSE}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'bym2'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'The BYM-model with the PC priors'}
#'               \item{constr = }{'TRUE'}
#'               \item{nrow.ncol = }{'FALSE'}
#'               \item{augmented = }{'TRUE'}
#'               \item{aug.factor = }{'2'}
#'               \item{aug.constr = }{'2'}
#'               \item{n.div.by = }{'NULL'}
#'               \item{n.required = }{'TRUE'}
#'               \item{set.default.values = }{'TRUE'}
#'               \item{status = }{'experimental'}
#'               \item{pdf = }{'bym2'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 2.
#'          \describe{
#'           \item{Hyperparameter 'theta1'}{
#'              \describe{
#'             \item{hyperid = }{11001}
#'             \item{name = }{log precision}
#'             \item{short.name = }{prec}
#'             \item{prior = }{pc.prec}
#'             \item{param = }{1 0.01}
#'             \item{initial = }{4}
#'             \item{fixed = }{FALSE}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta2'}{
#'              \describe{
#'             \item{hyperid = }{11002}
#'             \item{name = }{logit phi}
#'             \item{short.name = }{phi}
#'             \item{prior = }{pc}
#'             \item{param = }{0.5 0.5}
#'             \item{initial = }{-3}
#'             \item{fixed = }{FALSE}
#'             \item{to.theta = }{`function(x) log(x / (1 - x))`}
#'             \item{from.theta = }{`function(x) exp(x) / (1 + exp(x))`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'besagproper'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'A proper version of the Besag model'}
#'               \item{constr = }{'FALSE'}
#'               \item{nrow.ncol = }{'FALSE'}
#'               \item{augmented = }{'FALSE'}
#'               \item{aug.factor = }{'1'}
#'               \item{aug.constr = }{'NULL'}
#'               \item{n.div.by = }{'NULL'}
#'               \item{n.required = }{'TRUE'}
#'               \item{set.default.values = }{'TRUE'}
#'               \item{status = }{'experimental'}
#'               \item{pdf = }{'besagproper'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 2.
#'          \describe{
#'           \item{Hyperparameter 'theta1'}{
#'              \describe{
#'             \item{hyperid = }{12001}
#'             \item{name = }{log precision}
#'             \item{short.name = }{prec}
#'             \item{prior = }{loggamma}
#'             \item{param = }{1 5e-04}
#'             \item{initial = }{2}
#'             \item{fixed = }{FALSE}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta2'}{
#'              \describe{
#'             \item{hyperid = }{12002}
#'             \item{name = }{log diagonal}
#'             \item{short.name = }{diag}
#'             \item{prior = }{loggamma}
#'             \item{param = }{1 1}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'besagproper2'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'An alternative proper version of the Besag model'}
#'               \item{constr = }{'FALSE'}
#'               \item{nrow.ncol = }{'FALSE'}
#'               \item{augmented = }{'FALSE'}
#'               \item{aug.factor = }{'1'}
#'               \item{aug.constr = }{'NULL'}
#'               \item{n.div.by = }{'NULL'}
#'               \item{n.required = }{'TRUE'}
#'               \item{set.default.values = }{'TRUE'}
#'               \item{status = }{'experimental'}
#'               \item{pdf = }{'besagproper2'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 2.
#'          \describe{
#'           \item{Hyperparameter 'theta1'}{
#'              \describe{
#'             \item{hyperid = }{13001}
#'             \item{name = }{log precision}
#'             \item{short.name = }{prec}
#'             \item{prior = }{loggamma}
#'             \item{param = }{1 5e-04}
#'             \item{initial = }{2}
#'             \item{fixed = }{FALSE}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta2'}{
#'              \describe{
#'             \item{hyperid = }{13002}
#'             \item{name = }{logit lambda}
#'             \item{short.name = }{lambda}
#'             \item{prior = }{gaussian}
#'             \item{param = }{0 0.45}
#'             \item{initial = }{3}
#'             \item{fixed = }{FALSE}
#'             \item{to.theta = }{`function(x) log(x / (1 - x))`}
#'             \item{from.theta = }{`function(x) exp(x) / (1 + exp(x))`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'fgn'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'Fractional Gaussian noise model'}
#'               \item{constr = }{'FALSE'}
#'               \item{nrow.ncol = }{'FALSE'}
#'               \item{augmented = }{'TRUE'}
#'               \item{aug.factor = }{'5'}
#'               \item{aug.constr = }{'1'}
#'               \item{n.div.by = }{'NULL'}
#'               \item{n.required = }{'FALSE'}
#'               \item{set.default.values = }{'TRUE'}
#'               \item{order.default = }{'4'}
#'               \item{order.defined = }{'3 4'}
#'               \item{pdf = }{'fgn'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 2.
#'          \describe{
#'           \item{Hyperparameter 'theta1'}{
#'              \describe{
#'             \item{hyperid = }{13101}
#'             \item{name = }{log precision}
#'             \item{short.name = }{prec}
#'             \item{prior = }{pc.prec}
#'             \item{param = }{3 0.01}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta2'}{
#'              \describe{
#'             \item{hyperid = }{13102}
#'             \item{name = }{logit H}
#'             \item{short.name = }{H}
#'             \item{prior = }{pcfgnh}
#'             \item{param = }{0.9 0.1}
#'             \item{initial = }{2}
#'             \item{fixed = }{FALSE}
#'             \item{to.theta = }{`function(x) log((2 * x - 1) / (2 * (1 - x)))`}
#'             \item{from.theta = }{`function(x) 0.5 + 0.5 * exp(x) / (1 + exp(x))`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'fgn2'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'Fractional Gaussian noise model (alt 2)'}
#'               \item{constr = }{'FALSE'}
#'               \item{nrow.ncol = }{'FALSE'}
#'               \item{augmented = }{'TRUE'}
#'               \item{aug.factor = }{'4'}
#'               \item{aug.constr = }{'1'}
#'               \item{n.div.by = }{'NULL'}
#'               \item{n.required = }{'FALSE'}
#'               \item{set.default.values = }{'TRUE'}
#'               \item{order.default = }{'4'}
#'               \item{order.defined = }{'3 4'}
#'               \item{pdf = }{'fgn'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 2.
#'          \describe{
#'           \item{Hyperparameter 'theta1'}{
#'              \describe{
#'             \item{hyperid = }{13111}
#'             \item{name = }{log precision}
#'             \item{short.name = }{prec}
#'             \item{prior = }{pc.prec}
#'             \item{param = }{3 0.01}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta2'}{
#'              \describe{
#'             \item{hyperid = }{13112}
#'             \item{name = }{logit H}
#'             \item{short.name = }{H}
#'             \item{prior = }{pcfgnh}
#'             \item{param = }{0.9 0.1}
#'             \item{initial = }{2}
#'             \item{fixed = }{FALSE}
#'             \item{to.theta = }{`function(x) log((2 * x - 1) / (2 * (1 - x)))`}
#'             \item{from.theta = }{`function(x) 0.5 + 0.5 * exp(x) / (1 + exp(x))`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'ar1'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'Auto-regressive model of order 1 (AR(1))'}
#'               \item{constr = }{'FALSE'}
#'               \item{nrow.ncol = }{'FALSE'}
#'               \item{augmented = }{'FALSE'}
#'               \item{aug.factor = }{'1'}
#'               \item{aug.constr = }{'NULL'}
#'               \item{n.div.by = }{'NULL'}
#'               \item{n.required = }{'FALSE'}
#'               \item{set.default.values = }{'FALSE'}
#'               \item{pdf = }{'ar1'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 3.
#'          \describe{
#'           \item{Hyperparameter 'theta1'}{
#'              \describe{
#'             \item{hyperid = }{14001}
#'             \item{name = }{log precision}
#'             \item{short.name = }{prec}
#'             \item{prior = }{loggamma}
#'             \item{param = }{1 5e-05}
#'             \item{initial = }{4}
#'             \item{fixed = }{FALSE}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta2'}{
#'              \describe{
#'             \item{hyperid = }{14002}
#'             \item{name = }{logit lag one correlation}
#'             \item{short.name = }{rho}
#'             \item{prior = }{normal}
#'             \item{param = }{0 0.15}
#'             \item{initial = }{2}
#'             \item{fixed = }{FALSE}
#'             \item{to.theta = }{`function(x) log((1 + x) / (1 - x))`}
#'             \item{from.theta = }{`function(x) 2 * exp(x) / (1 + exp(x)) - 1`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta3'}{
#'              \describe{
#'             \item{hyperid = }{14003}
#'             \item{name = }{mean}
#'             \item{short.name = }{mean}
#'             \item{prior = }{normal}
#'             \item{param = }{0 1}
#'             \item{initial = }{0}
#'             \item{fixed = }{TRUE}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'ar1c'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'Auto-regressive model of order 1 w/covariates'}
#'               \item{constr = }{'FALSE'}
#'               \item{nrow.ncol = }{'FALSE'}
#'               \item{augmented = }{'FALSE'}
#'               \item{aug.factor = }{'1'}
#'               \item{aug.constr = }{'NULL'}
#'               \item{n.div.by = }{'NULL'}
#'               \item{n.required = }{'FALSE'}
#'               \item{set.default.values = }{'TRUE'}
#'               \item{status = }{'experimental'}
#'               \item{pdf = }{'ar1c'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 2.
#'          \describe{
#'           \item{Hyperparameter 'theta1'}{
#'              \describe{
#'             \item{hyperid = }{14101}
#'             \item{name = }{log precision}
#'             \item{short.name = }{prec}
#'             \item{prior = }{pc.prec}
#'             \item{param = }{1 0.01}
#'             \item{initial = }{4}
#'             \item{fixed = }{FALSE}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta2'}{
#'              \describe{
#'             \item{hyperid = }{14102}
#'             \item{name = }{logit lag one correlation}
#'             \item{short.name = }{rho}
#'             \item{prior = }{pc.cor0}
#'             \item{param = }{0.5 0.5}
#'             \item{initial = }{2}
#'             \item{fixed = }{FALSE}
#'             \item{to.theta = }{`function(x) log((1 + x) / (1 - x))`}
#'             \item{from.theta = }{`function(x) 2 * exp(x) / (1 + exp(x)) - 1`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'ar'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'Auto-regressive model of order p (AR(p))'}
#'               \item{constr = }{'FALSE'}
#'               \item{nrow.ncol = }{'FALSE'}
#'               \item{augmented = }{'FALSE'}
#'               \item{aug.factor = }{'1'}
#'               \item{aug.constr = }{'NULL'}
#'               \item{n.div.by = }{'NULL'}
#'               \item{n.required = }{'FALSE'}
#'               \item{set.default.values = }{'FALSE'}
#'               \item{pdf = }{'ar'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 11.
#'          \describe{
#'           \item{Hyperparameter 'theta1'}{
#'              \describe{
#'             \item{hyperid = }{15001}
#'             \item{name = }{log precision}
#'             \item{short.name = }{prec}
#'             \item{initial = }{4}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{pc.prec}
#'             \item{param = }{3 0.01}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta2'}{
#'              \describe{
#'             \item{hyperid = }{15002}
#'             \item{name = }{pacf1}
#'             \item{short.name = }{pacf1}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{pc.cor0}
#'             \item{param = }{0.5 0.5}
#'             \item{to.theta = }{`function(x) log((1 + x) / (1 - x))`}
#'             \item{from.theta = }{`function(x) 2 * exp(x) / (1 + exp(x)) - 1`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta3'}{
#'              \describe{
#'             \item{hyperid = }{15003}
#'             \item{name = }{pacf2}
#'             \item{short.name = }{pacf2}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{pc.cor0}
#'             \item{param = }{0.5 0.4}
#'             \item{to.theta = }{`function(x) log((1 + x) / (1 - x))`}
#'             \item{from.theta = }{`function(x) 2 * exp(x) / (1 + exp(x)) - 1`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta4'}{
#'              \describe{
#'             \item{hyperid = }{15004}
#'             \item{name = }{pacf3}
#'             \item{short.name = }{pacf3}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{pc.cor0}
#'             \item{param = }{0.5 0.3}
#'             \item{to.theta = }{`function(x) log((1 + x) / (1 - x))`}
#'             \item{from.theta = }{`function(x) 2 * exp(x) / (1 + exp(x)) - 1`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta5'}{
#'              \describe{
#'             \item{hyperid = }{15005}
#'             \item{name = }{pacf4}
#'             \item{short.name = }{pacf4}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{pc.cor0}
#'             \item{param = }{0.5 0.2}
#'             \item{to.theta = }{`function(x) log((1 + x) / (1 - x))`}
#'             \item{from.theta = }{`function(x) 2 * exp(x) / (1 + exp(x)) - 1`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta6'}{
#'              \describe{
#'             \item{hyperid = }{15006}
#'             \item{name = }{pacf5}
#'             \item{short.name = }{pacf5}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{pc.cor0}
#'             \item{param = }{0.5 0.1}
#'             \item{to.theta = }{`function(x) log((1 + x) / (1 - x))`}
#'             \item{from.theta = }{`function(x) 2 * exp(x) / (1 + exp(x)) - 1`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta7'}{
#'              \describe{
#'             \item{hyperid = }{15007}
#'             \item{name = }{pacf6}
#'             \item{short.name = }{pacf6}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{pc.cor0}
#'             \item{param = }{0.5 0.1}
#'             \item{to.theta = }{`function(x) log((1 + x) / (1 - x))`}
#'             \item{from.theta = }{`function(x) 2 * exp(x) / (1 + exp(x)) - 1`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta8'}{
#'              \describe{
#'             \item{hyperid = }{15008}
#'             \item{name = }{pacf7}
#'             \item{short.name = }{pacf7}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{pc.cor0}
#'             \item{param = }{0.5 0.1}
#'             \item{to.theta = }{`function(x) log((1 + x) / (1 - x))`}
#'             \item{from.theta = }{`function(x) 2 * exp(x) / (1 + exp(x)) - 1`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta9'}{
#'              \describe{
#'             \item{hyperid = }{15009}
#'             \item{name = }{pacf8}
#'             \item{short.name = }{pacf8}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{pc.cor0}
#'             \item{param = }{0.5 0.1}
#'             \item{to.theta = }{`function(x) log((1 + x) / (1 - x))`}
#'             \item{from.theta = }{`function(x) 2 * exp(x) / (1 + exp(x)) - 1`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta10'}{
#'              \describe{
#'             \item{hyperid = }{15010}
#'             \item{name = }{pacf9}
#'             \item{short.name = }{pacf9}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{pc.cor0}
#'             \item{param = }{0.5 0.1}
#'             \item{to.theta = }{`function(x) log((1 + x) / (1 - x))`}
#'             \item{from.theta = }{`function(x) 2 * exp(x) / (1 + exp(x)) - 1`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta11'}{
#'              \describe{
#'             \item{hyperid = }{15011}
#'             \item{name = }{pacf10}
#'             \item{short.name = }{pacf10}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{pc.cor0}
#'             \item{param = }{0.5 0.1}
#'             \item{to.theta = }{`function(x) log((1 + x) / (1 - x))`}
#'             \item{from.theta = }{`function(x) 2 * exp(x) / (1 + exp(x)) - 1`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'ou'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'The Ornstein-Uhlenbeck process'}
#'               \item{constr = }{'FALSE'}
#'               \item{nrow.ncol = }{'FALSE'}
#'               \item{augmented = }{'FALSE'}
#'               \item{aug.factor = }{'1'}
#'               \item{aug.constr = }{'NULL'}
#'               \item{n.div.by = }{'NULL'}
#'               \item{n.required = }{'FALSE'}
#'               \item{set.default.values = }{'FALSE'}
#'               \item{pdf = }{'ou'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 2.
#'          \describe{
#'           \item{Hyperparameter 'theta1'}{
#'              \describe{
#'             \item{hyperid = }{16001}
#'             \item{name = }{log precision}
#'             \item{short.name = }{prec}
#'             \item{prior = }{loggamma}
#'             \item{param = }{1 5e-05}
#'             \item{initial = }{4}
#'             \item{fixed = }{FALSE}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta2'}{
#'              \describe{
#'             \item{hyperid = }{16002}
#'             \item{name = }{log phi}
#'             \item{short.name = }{phi}
#'             \item{prior = }{normal}
#'             \item{param = }{0 0.2}
#'             \item{initial = }{-1}
#'             \item{fixed = }{FALSE}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'intslope'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'Intecept-slope model with Wishart-prior'}
#'               \item{constr = }{'FALSE'}
#'               \item{nrow.ncol = }{'FALSE'}
#'               \item{augmented = }{'FALSE'}
#'               \item{aug.factor = }{'1'}
#'               \item{aug.constr = }{'NULL'}
#'               \item{n.div.by = }{'NULL'}
#'               \item{n.required = }{'FALSE'}
#'               \item{set.default.values = }{'TRUE'}
#'               \item{status = }{'experimental'}
#'               \item{pdf = }{'intslope'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 13.
#'          \describe{
#'           \item{Hyperparameter 'theta1'}{
#'              \describe{
#'             \item{hyperid = }{16101}
#'             \item{name = }{log precision1}
#'             \item{short.name = }{prec1}
#'             \item{initial = }{4}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{wishart2d}
#'             \item{param = }{4 1 1 0}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta2'}{
#'              \describe{
#'             \item{hyperid = }{16102}
#'             \item{name = }{log precision2}
#'             \item{short.name = }{prec2}
#'             \item{initial = }{4}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta3'}{
#'              \describe{
#'             \item{hyperid = }{16103}
#'             \item{name = }{logit correlation}
#'             \item{short.name = }{cor}
#'             \item{initial = }{4}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) log((1 + x) / (1 - x))`}
#'             \item{from.theta = }{`function(x) 2 * exp(x) / (1 + exp(x)) - 1`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta4'}{
#'              \describe{
#'             \item{hyperid = }{16104}
#'             \item{name = }{gamma1}
#'             \item{short.name = }{g1}
#'             \item{initial = }{1}
#'             \item{fixed = }{TRUE}
#'             \item{prior = }{normal}
#'             \item{param = }{1 36}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta5'}{
#'              \describe{
#'             \item{hyperid = }{16105}
#'             \item{name = }{gamma2}
#'             \item{short.name = }{g2}
#'             \item{initial = }{1}
#'             \item{fixed = }{TRUE}
#'             \item{prior = }{normal}
#'             \item{param = }{1 36}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta6'}{
#'              \describe{
#'             \item{hyperid = }{16106}
#'             \item{name = }{gamma3}
#'             \item{short.name = }{g3}
#'             \item{initial = }{1}
#'             \item{fixed = }{TRUE}
#'             \item{prior = }{normal}
#'             \item{param = }{1 36}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta7'}{
#'              \describe{
#'             \item{hyperid = }{16107}
#'             \item{name = }{gamma4}
#'             \item{short.name = }{g4}
#'             \item{initial = }{1}
#'             \item{fixed = }{TRUE}
#'             \item{prior = }{normal}
#'             \item{param = }{1 36}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta8'}{
#'              \describe{
#'             \item{hyperid = }{16108}
#'             \item{name = }{gamma5}
#'             \item{short.name = }{g5}
#'             \item{initial = }{1}
#'             \item{fixed = }{TRUE}
#'             \item{prior = }{normal}
#'             \item{param = }{1 36}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta9'}{
#'              \describe{
#'             \item{hyperid = }{16109}
#'             \item{name = }{gamma6}
#'             \item{short.name = }{g6}
#'             \item{initial = }{1}
#'             \item{fixed = }{TRUE}
#'             \item{prior = }{normal}
#'             \item{param = }{1 36}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta10'}{
#'              \describe{
#'             \item{hyperid = }{16110}
#'             \item{name = }{gamma7}
#'             \item{short.name = }{g7}
#'             \item{initial = }{1}
#'             \item{fixed = }{TRUE}
#'             \item{prior = }{normal}
#'             \item{param = }{1 36}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta11'}{
#'              \describe{
#'             \item{hyperid = }{16111}
#'             \item{name = }{gamma8}
#'             \item{short.name = }{g8}
#'             \item{initial = }{1}
#'             \item{fixed = }{TRUE}
#'             \item{prior = }{normal}
#'             \item{param = }{1 36}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta12'}{
#'              \describe{
#'             \item{hyperid = }{16112}
#'             \item{name = }{gamma9}
#'             \item{short.name = }{g9}
#'             \item{initial = }{1}
#'             \item{fixed = }{TRUE}
#'             \item{prior = }{normal}
#'             \item{param = }{1 36}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta13'}{
#'              \describe{
#'             \item{hyperid = }{16113}
#'             \item{name = }{gamma10}
#'             \item{short.name = }{g10}
#'             \item{initial = }{1}
#'             \item{fixed = }{TRUE}
#'             \item{prior = }{normal}
#'             \item{param = }{1 36}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'generic'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'A generic model'}
#'               \item{constr = }{'FALSE'}
#'               \item{nrow.ncol = }{'FALSE'}
#'               \item{augmented = }{'FALSE'}
#'               \item{aug.factor = }{'1'}
#'               \item{aug.constr = }{'NULL'}
#'               \item{n.div.by = }{'NULL'}
#'               \item{n.required = }{'TRUE'}
#'               \item{set.default.values = }{'TRUE'}
#'               \item{pdf = }{'generic0'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 1.
#'          \describe{
#'           \item{Hyperparameter 'theta'}{
#'              \describe{
#'             \item{hyperid = }{17001}
#'             \item{name = }{log precision}
#'             \item{short.name = }{prec}
#'             \item{prior = }{loggamma}
#'             \item{param = }{1 5e-05}
#'             \item{initial = }{4}
#'             \item{fixed = }{FALSE}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'generic0'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'A generic model (type 0)'}
#'               \item{constr = }{'FALSE'}
#'               \item{nrow.ncol = }{'FALSE'}
#'               \item{augmented = }{'FALSE'}
#'               \item{aug.factor = }{'1'}
#'               \item{aug.constr = }{'NULL'}
#'               \item{n.div.by = }{'NULL'}
#'               \item{n.required = }{'TRUE'}
#'               \item{set.default.values = }{'TRUE'}
#'               \item{pdf = }{'generic0'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 1.
#'          \describe{
#'           \item{Hyperparameter 'theta'}{
#'              \describe{
#'             \item{hyperid = }{18001}
#'             \item{name = }{log precision}
#'             \item{short.name = }{prec}
#'             \item{prior = }{loggamma}
#'             \item{param = }{1 5e-05}
#'             \item{initial = }{4}
#'             \item{fixed = }{FALSE}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'generic1'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'A generic model (type 1)'}
#'               \item{constr = }{'FALSE'}
#'               \item{nrow.ncol = }{'FALSE'}
#'               \item{augmented = }{'FALSE'}
#'               \item{aug.factor = }{'1'}
#'               \item{aug.constr = }{'NULL'}
#'               \item{n.div.by = }{'NULL'}
#'               \item{n.required = }{'TRUE'}
#'               \item{set.default.values = }{'TRUE'}
#'               \item{pdf = }{'generic1'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 2.
#'          \describe{
#'           \item{Hyperparameter 'theta1'}{
#'              \describe{
#'             \item{hyperid = }{19001}
#'             \item{name = }{log precision}
#'             \item{short.name = }{prec}
#'             \item{prior = }{loggamma}
#'             \item{param = }{1 5e-05}
#'             \item{initial = }{4}
#'             \item{fixed = }{FALSE}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta2'}{
#'              \describe{
#'             \item{hyperid = }{19002}
#'             \item{name = }{beta}
#'             \item{short.name = }{beta}
#'             \item{initial = }{2}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{gaussian}
#'             \item{param = }{0 0.1}
#'             \item{to.theta = }{`function(x) log(x / (1 - x))`}
#'             \item{from.theta = }{`function(x) exp(x) / (1 + exp(x))`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'generic2'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'A generic model (type 2)'}
#'               \item{constr = }{'FALSE'}
#'               \item{nrow.ncol = }{'FALSE'}
#'               \item{augmented = }{'FALSE'}
#'               \item{aug.factor = }{'2'}
#'               \item{aug.constr = }{'2'}
#'               \item{n.div.by = }{'NULL'}
#'               \item{n.required = }{'TRUE'}
#'               \item{set.default.values = }{'TRUE'}
#'               \item{pdf = }{'generic2'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 2.
#'          \describe{
#'           \item{Hyperparameter 'theta1'}{
#'              \describe{
#'             \item{hyperid = }{20001}
#'             \item{name = }{log precision cmatrix}
#'             \item{short.name = }{prec}
#'             \item{initial = }{4}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{loggamma}
#'             \item{param = }{1 5e-05}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta2'}{
#'              \describe{
#'             \item{hyperid = }{20002}
#'             \item{name = }{log precision random}
#'             \item{short.name = }{prec.random}
#'             \item{initial = }{4}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{loggamma}
#'             \item{param = }{1 0.001}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'generic3'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'A generic model (type 3)'}
#'               \item{constr = }{'FALSE'}
#'               \item{nrow.ncol = }{'FALSE'}
#'               \item{augmented = }{'FALSE'}
#'               \item{aug.factor = }{'1'}
#'               \item{aug.constr = }{'NULL'}
#'               \item{n.div.by = }{'NULL'}
#'               \item{n.required = }{'TRUE'}
#'               \item{set.default.values = }{'TRUE'}
#'               \item{status = }{'experimental'}
#'               \item{pdf = }{'generic3'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 11.
#'          \describe{
#'           \item{Hyperparameter 'theta1'}{
#'              \describe{
#'             \item{hyperid = }{21001}
#'             \item{name = }{log precision1}
#'             \item{short.name = }{prec1}
#'             \item{initial = }{4}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{loggamma}
#'             \item{param = }{1 5e-05}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta2'}{
#'              \describe{
#'             \item{hyperid = }{21002}
#'             \item{name = }{log precision2}
#'             \item{short.name = }{prec2}
#'             \item{initial = }{4}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{loggamma}
#'             \item{param = }{1 5e-05}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta3'}{
#'              \describe{
#'             \item{hyperid = }{21003}
#'             \item{name = }{log precision3}
#'             \item{short.name = }{prec3}
#'             \item{initial = }{4}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{loggamma}
#'             \item{param = }{1 5e-05}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta4'}{
#'              \describe{
#'             \item{hyperid = }{21004}
#'             \item{name = }{log precision4}
#'             \item{short.name = }{prec4}
#'             \item{initial = }{4}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{loggamma}
#'             \item{param = }{1 5e-05}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta5'}{
#'              \describe{
#'             \item{hyperid = }{21005}
#'             \item{name = }{log precision5}
#'             \item{short.name = }{prec5}
#'             \item{initial = }{4}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{loggamma}
#'             \item{param = }{1 5e-05}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta6'}{
#'              \describe{
#'             \item{hyperid = }{21006}
#'             \item{name = }{log precision6}
#'             \item{short.name = }{prec6}
#'             \item{initial = }{4}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{loggamma}
#'             \item{param = }{1 5e-05}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta7'}{
#'              \describe{
#'             \item{hyperid = }{21007}
#'             \item{name = }{log precision7}
#'             \item{short.name = }{prec7}
#'             \item{initial = }{4}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{loggamma}
#'             \item{param = }{1 5e-05}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta8'}{
#'              \describe{
#'             \item{hyperid = }{21008}
#'             \item{name = }{log precision8}
#'             \item{short.name = }{prec8}
#'             \item{initial = }{4}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{loggamma}
#'             \item{param = }{1 5e-05}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta9'}{
#'              \describe{
#'             \item{hyperid = }{21009}
#'             \item{name = }{log precision9}
#'             \item{short.name = }{prec9}
#'             \item{initial = }{4}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{loggamma}
#'             \item{param = }{1 5e-05}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta10'}{
#'              \describe{
#'             \item{hyperid = }{21010}
#'             \item{name = }{log precision10}
#'             \item{short.name = }{prec10}
#'             \item{initial = }{4}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{loggamma}
#'             \item{param = }{1 5e-05}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta11'}{
#'              \describe{
#'             \item{hyperid = }{21011}
#'             \item{name = }{log precision common}
#'             \item{short.name = }{prec.common}
#'             \item{initial = }{0}
#'             \item{fixed = }{TRUE}
#'             \item{prior = }{loggamma}
#'             \item{param = }{1 5e-05}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'spde'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'A SPDE model'}
#'               \item{constr = }{'FALSE'}
#'               \item{nrow.ncol = }{'FALSE'}
#'               \item{augmented = }{'FALSE'}
#'               \item{aug.factor = }{'1'}
#'               \item{aug.constr = }{'NULL'}
#'               \item{n.div.by = }{'NULL'}
#'               \item{n.required = }{'TRUE'}
#'               \item{set.default.values = }{'TRUE'}
#'               \item{pdf = }{'spde'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 4.
#'          \describe{
#'           \item{Hyperparameter 'theta1'}{
#'              \describe{
#'             \item{hyperid = }{22001}
#'             \item{name = }{theta.T}
#'             \item{short.name = }{T}
#'             \item{initial = }{2}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 1}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta2'}{
#'              \describe{
#'             \item{hyperid = }{22002}
#'             \item{name = }{theta.K}
#'             \item{short.name = }{K}
#'             \item{initial = }{-2}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 1}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta3'}{
#'              \describe{
#'             \item{hyperid = }{22003}
#'             \item{name = }{theta.KT}
#'             \item{short.name = }{KT}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 1}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta4'}{
#'              \describe{
#'             \item{hyperid = }{22004}
#'             \item{name = }{theta.OC}
#'             \item{short.name = }{OC}
#'             \item{initial = }{-20}
#'             \item{fixed = }{TRUE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 0.2}
#'             \item{to.theta = }{`function(x) log(x / (1 - x))`}
#'             \item{from.theta = }{`function(x) exp(x) / (1 + exp(x))`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'spde2'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'A SPDE2 model'}
#'               \item{constr = }{'FALSE'}
#'               \item{nrow.ncol = }{'FALSE'}
#'               \item{augmented = }{'FALSE'}
#'               \item{aug.factor = }{'1'}
#'               \item{aug.constr = }{'NULL'}
#'               \item{n.div.by = }{'NULL'}
#'               \item{n.required = }{'TRUE'}
#'               \item{set.default.values = }{'TRUE'}
#'               \item{pdf = }{'spde2'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 100.
#'          \describe{
#'           \item{Hyperparameter 'theta1'}{
#'              \describe{
#'             \item{hyperid = }{23001}
#'             \item{name = }{theta1}
#'             \item{short.name = }{t1}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{mvnorm}
#'             \item{param = }{1 1}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta2'}{
#'              \describe{
#'             \item{hyperid = }{23002}
#'             \item{name = }{theta2}
#'             \item{short.name = }{t2}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta3'}{
#'              \describe{
#'             \item{hyperid = }{23003}
#'             \item{name = }{theta3}
#'             \item{short.name = }{t3}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta4'}{
#'              \describe{
#'             \item{hyperid = }{23004}
#'             \item{name = }{theta4}
#'             \item{short.name = }{t4}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta5'}{
#'              \describe{
#'             \item{hyperid = }{23005}
#'             \item{name = }{theta5}
#'             \item{short.name = }{t5}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta6'}{
#'              \describe{
#'             \item{hyperid = }{23006}
#'             \item{name = }{theta6}
#'             \item{short.name = }{t6}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta7'}{
#'              \describe{
#'             \item{hyperid = }{23007}
#'             \item{name = }{theta7}
#'             \item{short.name = }{t7}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta8'}{
#'              \describe{
#'             \item{hyperid = }{23008}
#'             \item{name = }{theta8}
#'             \item{short.name = }{t8}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta9'}{
#'              \describe{
#'             \item{hyperid = }{23009}
#'             \item{name = }{theta9}
#'             \item{short.name = }{t9}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta10'}{
#'              \describe{
#'             \item{hyperid = }{23010}
#'             \item{name = }{theta10}
#'             \item{short.name = }{t10}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta11'}{
#'              \describe{
#'             \item{hyperid = }{23011}
#'             \item{name = }{theta11}
#'             \item{short.name = }{t11}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta12'}{
#'              \describe{
#'             \item{hyperid = }{23012}
#'             \item{name = }{theta12}
#'             \item{short.name = }{t12}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta13'}{
#'              \describe{
#'             \item{hyperid = }{23013}
#'             \item{name = }{theta13}
#'             \item{short.name = }{t13}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta14'}{
#'              \describe{
#'             \item{hyperid = }{23014}
#'             \item{name = }{theta14}
#'             \item{short.name = }{t14}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta15'}{
#'              \describe{
#'             \item{hyperid = }{23015}
#'             \item{name = }{theta15}
#'             \item{short.name = }{t15}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta16'}{
#'              \describe{
#'             \item{hyperid = }{23016}
#'             \item{name = }{theta16}
#'             \item{short.name = }{t16}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta17'}{
#'              \describe{
#'             \item{hyperid = }{23017}
#'             \item{name = }{theta17}
#'             \item{short.name = }{t17}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta18'}{
#'              \describe{
#'             \item{hyperid = }{23018}
#'             \item{name = }{theta18}
#'             \item{short.name = }{t18}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta19'}{
#'              \describe{
#'             \item{hyperid = }{23019}
#'             \item{name = }{theta19}
#'             \item{short.name = }{t19}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta20'}{
#'              \describe{
#'             \item{hyperid = }{23020}
#'             \item{name = }{theta20}
#'             \item{short.name = }{t20}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta21'}{
#'              \describe{
#'             \item{hyperid = }{23021}
#'             \item{name = }{theta21}
#'             \item{short.name = }{t21}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta22'}{
#'              \describe{
#'             \item{hyperid = }{23022}
#'             \item{name = }{theta22}
#'             \item{short.name = }{t22}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta23'}{
#'              \describe{
#'             \item{hyperid = }{23023}
#'             \item{name = }{theta23}
#'             \item{short.name = }{t23}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta24'}{
#'              \describe{
#'             \item{hyperid = }{23024}
#'             \item{name = }{theta24}
#'             \item{short.name = }{t24}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta25'}{
#'              \describe{
#'             \item{hyperid = }{23025}
#'             \item{name = }{theta25}
#'             \item{short.name = }{t25}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta26'}{
#'              \describe{
#'             \item{hyperid = }{23026}
#'             \item{name = }{theta26}
#'             \item{short.name = }{t26}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta27'}{
#'              \describe{
#'             \item{hyperid = }{23027}
#'             \item{name = }{theta27}
#'             \item{short.name = }{t27}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta28'}{
#'              \describe{
#'             \item{hyperid = }{23028}
#'             \item{name = }{theta28}
#'             \item{short.name = }{t28}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta29'}{
#'              \describe{
#'             \item{hyperid = }{23029}
#'             \item{name = }{theta29}
#'             \item{short.name = }{t29}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta30'}{
#'              \describe{
#'             \item{hyperid = }{23030}
#'             \item{name = }{theta30}
#'             \item{short.name = }{t30}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta31'}{
#'              \describe{
#'             \item{hyperid = }{23031}
#'             \item{name = }{theta31}
#'             \item{short.name = }{t31}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta32'}{
#'              \describe{
#'             \item{hyperid = }{23032}
#'             \item{name = }{theta32}
#'             \item{short.name = }{t32}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta33'}{
#'              \describe{
#'             \item{hyperid = }{23033}
#'             \item{name = }{theta33}
#'             \item{short.name = }{t33}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta34'}{
#'              \describe{
#'             \item{hyperid = }{23034}
#'             \item{name = }{theta34}
#'             \item{short.name = }{t34}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta35'}{
#'              \describe{
#'             \item{hyperid = }{23035}
#'             \item{name = }{theta35}
#'             \item{short.name = }{t35}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta36'}{
#'              \describe{
#'             \item{hyperid = }{23036}
#'             \item{name = }{theta36}
#'             \item{short.name = }{t36}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta37'}{
#'              \describe{
#'             \item{hyperid = }{23037}
#'             \item{name = }{theta37}
#'             \item{short.name = }{t37}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta38'}{
#'              \describe{
#'             \item{hyperid = }{23038}
#'             \item{name = }{theta38}
#'             \item{short.name = }{t38}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta39'}{
#'              \describe{
#'             \item{hyperid = }{23039}
#'             \item{name = }{theta39}
#'             \item{short.name = }{t39}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta40'}{
#'              \describe{
#'             \item{hyperid = }{23040}
#'             \item{name = }{theta40}
#'             \item{short.name = }{t40}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta41'}{
#'              \describe{
#'             \item{hyperid = }{23041}
#'             \item{name = }{theta41}
#'             \item{short.name = }{t41}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta42'}{
#'              \describe{
#'             \item{hyperid = }{23042}
#'             \item{name = }{theta42}
#'             \item{short.name = }{t42}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta43'}{
#'              \describe{
#'             \item{hyperid = }{23043}
#'             \item{name = }{theta43}
#'             \item{short.name = }{t43}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta44'}{
#'              \describe{
#'             \item{hyperid = }{23044}
#'             \item{name = }{theta44}
#'             \item{short.name = }{t44}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta45'}{
#'              \describe{
#'             \item{hyperid = }{23045}
#'             \item{name = }{theta45}
#'             \item{short.name = }{t45}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta46'}{
#'              \describe{
#'             \item{hyperid = }{23046}
#'             \item{name = }{theta46}
#'             \item{short.name = }{t46}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta47'}{
#'              \describe{
#'             \item{hyperid = }{23047}
#'             \item{name = }{theta47}
#'             \item{short.name = }{t47}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta48'}{
#'              \describe{
#'             \item{hyperid = }{23048}
#'             \item{name = }{theta48}
#'             \item{short.name = }{t48}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta49'}{
#'              \describe{
#'             \item{hyperid = }{23049}
#'             \item{name = }{theta49}
#'             \item{short.name = }{t49}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta50'}{
#'              \describe{
#'             \item{hyperid = }{23050}
#'             \item{name = }{theta50}
#'             \item{short.name = }{t50}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta51'}{
#'              \describe{
#'             \item{hyperid = }{23051}
#'             \item{name = }{theta51}
#'             \item{short.name = }{t51}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta52'}{
#'              \describe{
#'             \item{hyperid = }{23052}
#'             \item{name = }{theta52}
#'             \item{short.name = }{t52}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta53'}{
#'              \describe{
#'             \item{hyperid = }{23053}
#'             \item{name = }{theta53}
#'             \item{short.name = }{t53}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta54'}{
#'              \describe{
#'             \item{hyperid = }{23054}
#'             \item{name = }{theta54}
#'             \item{short.name = }{t54}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta55'}{
#'              \describe{
#'             \item{hyperid = }{23055}
#'             \item{name = }{theta55}
#'             \item{short.name = }{t55}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta56'}{
#'              \describe{
#'             \item{hyperid = }{23056}
#'             \item{name = }{theta56}
#'             \item{short.name = }{t56}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta57'}{
#'              \describe{
#'             \item{hyperid = }{23057}
#'             \item{name = }{theta57}
#'             \item{short.name = }{t57}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta58'}{
#'              \describe{
#'             \item{hyperid = }{23058}
#'             \item{name = }{theta58}
#'             \item{short.name = }{t58}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta59'}{
#'              \describe{
#'             \item{hyperid = }{23059}
#'             \item{name = }{theta59}
#'             \item{short.name = }{t59}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta60'}{
#'              \describe{
#'             \item{hyperid = }{23060}
#'             \item{name = }{theta60}
#'             \item{short.name = }{t60}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta61'}{
#'              \describe{
#'             \item{hyperid = }{23061}
#'             \item{name = }{theta61}
#'             \item{short.name = }{t61}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta62'}{
#'              \describe{
#'             \item{hyperid = }{23062}
#'             \item{name = }{theta62}
#'             \item{short.name = }{t62}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta63'}{
#'              \describe{
#'             \item{hyperid = }{23063}
#'             \item{name = }{theta63}
#'             \item{short.name = }{t63}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta64'}{
#'              \describe{
#'             \item{hyperid = }{23064}
#'             \item{name = }{theta64}
#'             \item{short.name = }{t64}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta65'}{
#'              \describe{
#'             \item{hyperid = }{23065}
#'             \item{name = }{theta65}
#'             \item{short.name = }{t65}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta66'}{
#'              \describe{
#'             \item{hyperid = }{23066}
#'             \item{name = }{theta66}
#'             \item{short.name = }{t66}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta67'}{
#'              \describe{
#'             \item{hyperid = }{23067}
#'             \item{name = }{theta67}
#'             \item{short.name = }{t67}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta68'}{
#'              \describe{
#'             \item{hyperid = }{23068}
#'             \item{name = }{theta68}
#'             \item{short.name = }{t68}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta69'}{
#'              \describe{
#'             \item{hyperid = }{23069}
#'             \item{name = }{theta69}
#'             \item{short.name = }{t69}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta70'}{
#'              \describe{
#'             \item{hyperid = }{23070}
#'             \item{name = }{theta70}
#'             \item{short.name = }{t70}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta71'}{
#'              \describe{
#'             \item{hyperid = }{23071}
#'             \item{name = }{theta71}
#'             \item{short.name = }{t71}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta72'}{
#'              \describe{
#'             \item{hyperid = }{23072}
#'             \item{name = }{theta72}
#'             \item{short.name = }{t72}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta73'}{
#'              \describe{
#'             \item{hyperid = }{23073}
#'             \item{name = }{theta73}
#'             \item{short.name = }{t73}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta74'}{
#'              \describe{
#'             \item{hyperid = }{23074}
#'             \item{name = }{theta74}
#'             \item{short.name = }{t74}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta75'}{
#'              \describe{
#'             \item{hyperid = }{23075}
#'             \item{name = }{theta75}
#'             \item{short.name = }{t75}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta76'}{
#'              \describe{
#'             \item{hyperid = }{23076}
#'             \item{name = }{theta76}
#'             \item{short.name = }{t76}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta77'}{
#'              \describe{
#'             \item{hyperid = }{23077}
#'             \item{name = }{theta77}
#'             \item{short.name = }{t77}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta78'}{
#'              \describe{
#'             \item{hyperid = }{23078}
#'             \item{name = }{theta78}
#'             \item{short.name = }{t78}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta79'}{
#'              \describe{
#'             \item{hyperid = }{23079}
#'             \item{name = }{theta79}
#'             \item{short.name = }{t79}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta80'}{
#'              \describe{
#'             \item{hyperid = }{23080}
#'             \item{name = }{theta80}
#'             \item{short.name = }{t80}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta81'}{
#'              \describe{
#'             \item{hyperid = }{23081}
#'             \item{name = }{theta81}
#'             \item{short.name = }{t81}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta82'}{
#'              \describe{
#'             \item{hyperid = }{23082}
#'             \item{name = }{theta82}
#'             \item{short.name = }{t82}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta83'}{
#'              \describe{
#'             \item{hyperid = }{23083}
#'             \item{name = }{theta83}
#'             \item{short.name = }{t83}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta84'}{
#'              \describe{
#'             \item{hyperid = }{23084}
#'             \item{name = }{theta84}
#'             \item{short.name = }{t84}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta85'}{
#'              \describe{
#'             \item{hyperid = }{23085}
#'             \item{name = }{theta85}
#'             \item{short.name = }{t85}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta86'}{
#'              \describe{
#'             \item{hyperid = }{23086}
#'             \item{name = }{theta86}
#'             \item{short.name = }{t86}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta87'}{
#'              \describe{
#'             \item{hyperid = }{23087}
#'             \item{name = }{theta87}
#'             \item{short.name = }{t87}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta88'}{
#'              \describe{
#'             \item{hyperid = }{23088}
#'             \item{name = }{theta88}
#'             \item{short.name = }{t88}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta89'}{
#'              \describe{
#'             \item{hyperid = }{23089}
#'             \item{name = }{theta89}
#'             \item{short.name = }{t89}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta90'}{
#'              \describe{
#'             \item{hyperid = }{23090}
#'             \item{name = }{theta90}
#'             \item{short.name = }{t90}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta91'}{
#'              \describe{
#'             \item{hyperid = }{23091}
#'             \item{name = }{theta91}
#'             \item{short.name = }{t91}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta92'}{
#'              \describe{
#'             \item{hyperid = }{23092}
#'             \item{name = }{theta92}
#'             \item{short.name = }{t92}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta93'}{
#'              \describe{
#'             \item{hyperid = }{23093}
#'             \item{name = }{theta93}
#'             \item{short.name = }{t93}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta94'}{
#'              \describe{
#'             \item{hyperid = }{23094}
#'             \item{name = }{theta94}
#'             \item{short.name = }{t94}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta95'}{
#'              \describe{
#'             \item{hyperid = }{23095}
#'             \item{name = }{theta95}
#'             \item{short.name = }{t95}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta96'}{
#'              \describe{
#'             \item{hyperid = }{23096}
#'             \item{name = }{theta96}
#'             \item{short.name = }{t96}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta97'}{
#'              \describe{
#'             \item{hyperid = }{23097}
#'             \item{name = }{theta97}
#'             \item{short.name = }{t97}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta98'}{
#'              \describe{
#'             \item{hyperid = }{23098}
#'             \item{name = }{theta98}
#'             \item{short.name = }{t98}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta99'}{
#'              \describe{
#'             \item{hyperid = }{23099}
#'             \item{name = }{theta99}
#'             \item{short.name = }{t99}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta100'}{
#'              \describe{
#'             \item{hyperid = }{23100}
#'             \item{name = }{theta100}
#'             \item{short.name = }{t100}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'spde3'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'A SPDE3 model'}
#'               \item{constr = }{'FALSE'}
#'               \item{nrow.ncol = }{'FALSE'}
#'               \item{augmented = }{'FALSE'}
#'               \item{aug.factor = }{'1'}
#'               \item{aug.constr = }{'NULL'}
#'               \item{n.div.by = }{'NULL'}
#'               \item{n.required = }{'TRUE'}
#'               \item{set.default.values = }{'TRUE'}
#'               \item{pdf = }{'spde3'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 100.
#'          \describe{
#'           \item{Hyperparameter 'theta1'}{
#'              \describe{
#'             \item{hyperid = }{24001}
#'             \item{name = }{theta1}
#'             \item{short.name = }{t1}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{mvnorm}
#'             \item{param = }{1 1}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta2'}{
#'              \describe{
#'             \item{hyperid = }{24002}
#'             \item{name = }{theta2}
#'             \item{short.name = }{t2}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta3'}{
#'              \describe{
#'             \item{hyperid = }{24003}
#'             \item{name = }{theta3}
#'             \item{short.name = }{t3}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta4'}{
#'              \describe{
#'             \item{hyperid = }{24004}
#'             \item{name = }{theta4}
#'             \item{short.name = }{t4}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta5'}{
#'              \describe{
#'             \item{hyperid = }{24005}
#'             \item{name = }{theta5}
#'             \item{short.name = }{t5}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta6'}{
#'              \describe{
#'             \item{hyperid = }{24006}
#'             \item{name = }{theta6}
#'             \item{short.name = }{t6}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta7'}{
#'              \describe{
#'             \item{hyperid = }{24007}
#'             \item{name = }{theta7}
#'             \item{short.name = }{t7}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta8'}{
#'              \describe{
#'             \item{hyperid = }{24008}
#'             \item{name = }{theta8}
#'             \item{short.name = }{t8}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta9'}{
#'              \describe{
#'             \item{hyperid = }{24009}
#'             \item{name = }{theta9}
#'             \item{short.name = }{t9}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta10'}{
#'              \describe{
#'             \item{hyperid = }{24010}
#'             \item{name = }{theta10}
#'             \item{short.name = }{t10}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta11'}{
#'              \describe{
#'             \item{hyperid = }{24011}
#'             \item{name = }{theta11}
#'             \item{short.name = }{t11}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta12'}{
#'              \describe{
#'             \item{hyperid = }{24012}
#'             \item{name = }{theta12}
#'             \item{short.name = }{t12}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta13'}{
#'              \describe{
#'             \item{hyperid = }{24013}
#'             \item{name = }{theta13}
#'             \item{short.name = }{t13}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta14'}{
#'              \describe{
#'             \item{hyperid = }{24014}
#'             \item{name = }{theta14}
#'             \item{short.name = }{t14}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta15'}{
#'              \describe{
#'             \item{hyperid = }{24015}
#'             \item{name = }{theta15}
#'             \item{short.name = }{t15}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta16'}{
#'              \describe{
#'             \item{hyperid = }{24016}
#'             \item{name = }{theta16}
#'             \item{short.name = }{t16}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta17'}{
#'              \describe{
#'             \item{hyperid = }{24017}
#'             \item{name = }{theta17}
#'             \item{short.name = }{t17}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta18'}{
#'              \describe{
#'             \item{hyperid = }{24018}
#'             \item{name = }{theta18}
#'             \item{short.name = }{t18}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta19'}{
#'              \describe{
#'             \item{hyperid = }{24019}
#'             \item{name = }{theta19}
#'             \item{short.name = }{t19}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta20'}{
#'              \describe{
#'             \item{hyperid = }{24020}
#'             \item{name = }{theta20}
#'             \item{short.name = }{t20}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta21'}{
#'              \describe{
#'             \item{hyperid = }{24021}
#'             \item{name = }{theta21}
#'             \item{short.name = }{t21}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta22'}{
#'              \describe{
#'             \item{hyperid = }{24022}
#'             \item{name = }{theta22}
#'             \item{short.name = }{t22}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta23'}{
#'              \describe{
#'             \item{hyperid = }{24023}
#'             \item{name = }{theta23}
#'             \item{short.name = }{t23}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta24'}{
#'              \describe{
#'             \item{hyperid = }{24024}
#'             \item{name = }{theta24}
#'             \item{short.name = }{t24}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta25'}{
#'              \describe{
#'             \item{hyperid = }{24025}
#'             \item{name = }{theta25}
#'             \item{short.name = }{t25}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta26'}{
#'              \describe{
#'             \item{hyperid = }{24026}
#'             \item{name = }{theta26}
#'             \item{short.name = }{t26}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta27'}{
#'              \describe{
#'             \item{hyperid = }{24027}
#'             \item{name = }{theta27}
#'             \item{short.name = }{t27}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta28'}{
#'              \describe{
#'             \item{hyperid = }{24028}
#'             \item{name = }{theta28}
#'             \item{short.name = }{t28}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta29'}{
#'              \describe{
#'             \item{hyperid = }{24029}
#'             \item{name = }{theta29}
#'             \item{short.name = }{t29}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta30'}{
#'              \describe{
#'             \item{hyperid = }{24030}
#'             \item{name = }{theta30}
#'             \item{short.name = }{t30}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta31'}{
#'              \describe{
#'             \item{hyperid = }{24031}
#'             \item{name = }{theta31}
#'             \item{short.name = }{t31}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta32'}{
#'              \describe{
#'             \item{hyperid = }{24032}
#'             \item{name = }{theta32}
#'             \item{short.name = }{t32}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta33'}{
#'              \describe{
#'             \item{hyperid = }{24033}
#'             \item{name = }{theta33}
#'             \item{short.name = }{t33}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta34'}{
#'              \describe{
#'             \item{hyperid = }{24034}
#'             \item{name = }{theta34}
#'             \item{short.name = }{t34}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta35'}{
#'              \describe{
#'             \item{hyperid = }{24035}
#'             \item{name = }{theta35}
#'             \item{short.name = }{t35}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta36'}{
#'              \describe{
#'             \item{hyperid = }{24036}
#'             \item{name = }{theta36}
#'             \item{short.name = }{t36}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta37'}{
#'              \describe{
#'             \item{hyperid = }{24037}
#'             \item{name = }{theta37}
#'             \item{short.name = }{t37}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta38'}{
#'              \describe{
#'             \item{hyperid = }{24038}
#'             \item{name = }{theta38}
#'             \item{short.name = }{t38}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta39'}{
#'              \describe{
#'             \item{hyperid = }{24039}
#'             \item{name = }{theta39}
#'             \item{short.name = }{t39}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta40'}{
#'              \describe{
#'             \item{hyperid = }{24040}
#'             \item{name = }{theta40}
#'             \item{short.name = }{t40}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta41'}{
#'              \describe{
#'             \item{hyperid = }{24041}
#'             \item{name = }{theta41}
#'             \item{short.name = }{t41}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta42'}{
#'              \describe{
#'             \item{hyperid = }{24042}
#'             \item{name = }{theta42}
#'             \item{short.name = }{t42}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta43'}{
#'              \describe{
#'             \item{hyperid = }{24043}
#'             \item{name = }{theta43}
#'             \item{short.name = }{t43}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta44'}{
#'              \describe{
#'             \item{hyperid = }{24044}
#'             \item{name = }{theta44}
#'             \item{short.name = }{t44}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta45'}{
#'              \describe{
#'             \item{hyperid = }{24045}
#'             \item{name = }{theta45}
#'             \item{short.name = }{t45}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta46'}{
#'              \describe{
#'             \item{hyperid = }{24046}
#'             \item{name = }{theta46}
#'             \item{short.name = }{t46}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta47'}{
#'              \describe{
#'             \item{hyperid = }{24047}
#'             \item{name = }{theta47}
#'             \item{short.name = }{t47}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta48'}{
#'              \describe{
#'             \item{hyperid = }{24048}
#'             \item{name = }{theta48}
#'             \item{short.name = }{t48}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta49'}{
#'              \describe{
#'             \item{hyperid = }{24049}
#'             \item{name = }{theta49}
#'             \item{short.name = }{t49}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta50'}{
#'              \describe{
#'             \item{hyperid = }{24050}
#'             \item{name = }{theta50}
#'             \item{short.name = }{t50}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta51'}{
#'              \describe{
#'             \item{hyperid = }{24051}
#'             \item{name = }{theta51}
#'             \item{short.name = }{t51}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta52'}{
#'              \describe{
#'             \item{hyperid = }{24052}
#'             \item{name = }{theta52}
#'             \item{short.name = }{t52}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta53'}{
#'              \describe{
#'             \item{hyperid = }{24053}
#'             \item{name = }{theta53}
#'             \item{short.name = }{t53}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta54'}{
#'              \describe{
#'             \item{hyperid = }{24054}
#'             \item{name = }{theta54}
#'             \item{short.name = }{t54}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta55'}{
#'              \describe{
#'             \item{hyperid = }{24055}
#'             \item{name = }{theta55}
#'             \item{short.name = }{t55}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta56'}{
#'              \describe{
#'             \item{hyperid = }{24056}
#'             \item{name = }{theta56}
#'             \item{short.name = }{t56}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta57'}{
#'              \describe{
#'             \item{hyperid = }{24057}
#'             \item{name = }{theta57}
#'             \item{short.name = }{t57}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta58'}{
#'              \describe{
#'             \item{hyperid = }{24058}
#'             \item{name = }{theta58}
#'             \item{short.name = }{t58}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta59'}{
#'              \describe{
#'             \item{hyperid = }{24059}
#'             \item{name = }{theta59}
#'             \item{short.name = }{t59}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta60'}{
#'              \describe{
#'             \item{hyperid = }{24060}
#'             \item{name = }{theta60}
#'             \item{short.name = }{t60}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta61'}{
#'              \describe{
#'             \item{hyperid = }{24061}
#'             \item{name = }{theta61}
#'             \item{short.name = }{t61}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta62'}{
#'              \describe{
#'             \item{hyperid = }{24062}
#'             \item{name = }{theta62}
#'             \item{short.name = }{t62}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta63'}{
#'              \describe{
#'             \item{hyperid = }{24063}
#'             \item{name = }{theta63}
#'             \item{short.name = }{t63}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta64'}{
#'              \describe{
#'             \item{hyperid = }{24064}
#'             \item{name = }{theta64}
#'             \item{short.name = }{t64}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta65'}{
#'              \describe{
#'             \item{hyperid = }{24065}
#'             \item{name = }{theta65}
#'             \item{short.name = }{t65}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta66'}{
#'              \describe{
#'             \item{hyperid = }{24066}
#'             \item{name = }{theta66}
#'             \item{short.name = }{t66}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta67'}{
#'              \describe{
#'             \item{hyperid = }{24067}
#'             \item{name = }{theta67}
#'             \item{short.name = }{t67}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta68'}{
#'              \describe{
#'             \item{hyperid = }{24068}
#'             \item{name = }{theta68}
#'             \item{short.name = }{t68}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta69'}{
#'              \describe{
#'             \item{hyperid = }{24069}
#'             \item{name = }{theta69}
#'             \item{short.name = }{t69}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta70'}{
#'              \describe{
#'             \item{hyperid = }{24070}
#'             \item{name = }{theta70}
#'             \item{short.name = }{t70}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta71'}{
#'              \describe{
#'             \item{hyperid = }{24071}
#'             \item{name = }{theta71}
#'             \item{short.name = }{t71}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta72'}{
#'              \describe{
#'             \item{hyperid = }{24072}
#'             \item{name = }{theta72}
#'             \item{short.name = }{t72}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta73'}{
#'              \describe{
#'             \item{hyperid = }{24073}
#'             \item{name = }{theta73}
#'             \item{short.name = }{t73}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta74'}{
#'              \describe{
#'             \item{hyperid = }{24074}
#'             \item{name = }{theta74}
#'             \item{short.name = }{t74}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta75'}{
#'              \describe{
#'             \item{hyperid = }{24075}
#'             \item{name = }{theta75}
#'             \item{short.name = }{t75}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta76'}{
#'              \describe{
#'             \item{hyperid = }{24076}
#'             \item{name = }{theta76}
#'             \item{short.name = }{t76}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta77'}{
#'              \describe{
#'             \item{hyperid = }{24077}
#'             \item{name = }{theta77}
#'             \item{short.name = }{t77}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta78'}{
#'              \describe{
#'             \item{hyperid = }{24078}
#'             \item{name = }{theta78}
#'             \item{short.name = }{t78}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta79'}{
#'              \describe{
#'             \item{hyperid = }{24079}
#'             \item{name = }{theta79}
#'             \item{short.name = }{t79}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta80'}{
#'              \describe{
#'             \item{hyperid = }{24080}
#'             \item{name = }{theta80}
#'             \item{short.name = }{t80}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta81'}{
#'              \describe{
#'             \item{hyperid = }{24081}
#'             \item{name = }{theta81}
#'             \item{short.name = }{t81}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta82'}{
#'              \describe{
#'             \item{hyperid = }{24082}
#'             \item{name = }{theta82}
#'             \item{short.name = }{t82}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta83'}{
#'              \describe{
#'             \item{hyperid = }{24083}
#'             \item{name = }{theta83}
#'             \item{short.name = }{t83}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta84'}{
#'              \describe{
#'             \item{hyperid = }{24084}
#'             \item{name = }{theta84}
#'             \item{short.name = }{t84}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta85'}{
#'              \describe{
#'             \item{hyperid = }{24085}
#'             \item{name = }{theta85}
#'             \item{short.name = }{t85}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta86'}{
#'              \describe{
#'             \item{hyperid = }{24086}
#'             \item{name = }{theta86}
#'             \item{short.name = }{t86}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta87'}{
#'              \describe{
#'             \item{hyperid = }{24087}
#'             \item{name = }{theta87}
#'             \item{short.name = }{t87}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta88'}{
#'              \describe{
#'             \item{hyperid = }{24088}
#'             \item{name = }{theta88}
#'             \item{short.name = }{t88}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta89'}{
#'              \describe{
#'             \item{hyperid = }{24089}
#'             \item{name = }{theta89}
#'             \item{short.name = }{t89}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta90'}{
#'              \describe{
#'             \item{hyperid = }{24090}
#'             \item{name = }{theta90}
#'             \item{short.name = }{t90}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta91'}{
#'              \describe{
#'             \item{hyperid = }{24091}
#'             \item{name = }{theta91}
#'             \item{short.name = }{t91}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta92'}{
#'              \describe{
#'             \item{hyperid = }{24092}
#'             \item{name = }{theta92}
#'             \item{short.name = }{t92}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta93'}{
#'              \describe{
#'             \item{hyperid = }{24093}
#'             \item{name = }{theta93}
#'             \item{short.name = }{t93}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta94'}{
#'              \describe{
#'             \item{hyperid = }{24094}
#'             \item{name = }{theta94}
#'             \item{short.name = }{t94}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta95'}{
#'              \describe{
#'             \item{hyperid = }{24095}
#'             \item{name = }{theta95}
#'             \item{short.name = }{t95}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta96'}{
#'              \describe{
#'             \item{hyperid = }{24096}
#'             \item{name = }{theta96}
#'             \item{short.name = }{t96}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta97'}{
#'              \describe{
#'             \item{hyperid = }{24097}
#'             \item{name = }{theta97}
#'             \item{short.name = }{t97}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta98'}{
#'              \describe{
#'             \item{hyperid = }{24098}
#'             \item{name = }{theta98}
#'             \item{short.name = }{t98}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta99'}{
#'              \describe{
#'             \item{hyperid = }{24099}
#'             \item{name = }{theta99}
#'             \item{short.name = }{t99}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta100'}{
#'              \describe{
#'             \item{hyperid = }{24100}
#'             \item{name = }{theta100}
#'             \item{short.name = }{t100}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'iid1d'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'Gaussian random effect in dim=1 with Wishart prior'}
#'               \item{constr = }{'FALSE'}
#'               \item{nrow.ncol = }{'FALSE'}
#'               \item{augmented = }{'FALSE'}
#'               \item{aug.factor = }{'1'}
#'               \item{aug.constr = }{'NULL'}
#'               \item{n.div.by = }{'NULL'}
#'               \item{n.required = }{'FALSE'}
#'               \item{set.default.values = }{'TRUE'}
#'               \item{pdf = }{'iid123d'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 1.
#'          \describe{
#'           \item{Hyperparameter 'theta'}{
#'              \describe{
#'             \item{hyperid = }{25001}
#'             \item{name = }{precision}
#'             \item{short.name = }{prec}
#'             \item{initial = }{4}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{wishart1d}
#'             \item{param = }{2 1e-04}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'iid2d'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'Gaussian random effect in dim=2 with Wishart prior'}
#'               \item{constr = }{'FALSE'}
#'               \item{nrow.ncol = }{'FALSE'}
#'               \item{augmented = }{'TRUE'}
#'               \item{aug.factor = }{'1'}
#'               \item{aug.constr = }{'1 2'}
#'               \item{n.div.by = }{'2'}
#'               \item{n.required = }{'TRUE'}
#'               \item{set.default.values = }{'TRUE'}
#'               \item{pdf = }{'iid123d'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 3.
#'          \describe{
#'           \item{Hyperparameter 'theta1'}{
#'              \describe{
#'             \item{hyperid = }{26001}
#'             \item{name = }{log precision1}
#'             \item{short.name = }{prec1}
#'             \item{initial = }{4}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{wishart2d}
#'             \item{param = }{4 1 1 0}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta2'}{
#'              \describe{
#'             \item{hyperid = }{26002}
#'             \item{name = }{log precision2}
#'             \item{short.name = }{prec2}
#'             \item{initial = }{4}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta3'}{
#'              \describe{
#'             \item{hyperid = }{26003}
#'             \item{name = }{logit correlation}
#'             \item{short.name = }{cor}
#'             \item{initial = }{4}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) log((1 + x) / (1 - x))`}
#'             \item{from.theta = }{`function(x) 2 * exp(x) / (1 + exp(x)) - 1`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'iid3d'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'Gaussian random effect in dim=3 with Wishart prior'}
#'               \item{constr = }{'FALSE'}
#'               \item{nrow.ncol = }{'FALSE'}
#'               \item{augmented = }{'TRUE'}
#'               \item{aug.factor = }{'1'}
#'               \item{aug.constr = }{'1 2 3'}
#'               \item{n.div.by = }{'3'}
#'               \item{n.required = }{'TRUE'}
#'               \item{set.default.values = }{'TRUE'}
#'               \item{pdf = }{'iid123d'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 6.
#'          \describe{
#'           \item{Hyperparameter 'theta1'}{
#'              \describe{
#'             \item{hyperid = }{27001}
#'             \item{name = }{log precision1}
#'             \item{short.name = }{prec1}
#'             \item{initial = }{4}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{wishart3d}
#'             \item{param = }{7 1 1 1 0 0 0}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta2'}{
#'              \describe{
#'             \item{hyperid = }{27002}
#'             \item{name = }{log precision2}
#'             \item{short.name = }{prec2}
#'             \item{initial = }{4}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta3'}{
#'              \describe{
#'             \item{hyperid = }{27003}
#'             \item{name = }{log precision3}
#'             \item{short.name = }{prec3}
#'             \item{initial = }{4}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta4'}{
#'              \describe{
#'             \item{hyperid = }{27004}
#'             \item{name = }{logit correlation12}
#'             \item{short.name = }{cor12}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) log((1 + x) / (1 - x))`}
#'             \item{from.theta = }{`function(x) 2 * exp(x) / (1 + exp(x)) - 1`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta5'}{
#'              \describe{
#'             \item{hyperid = }{27005}
#'             \item{name = }{logit correlation13}
#'             \item{short.name = }{cor13}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) log((1 + x) / (1 - x))`}
#'             \item{from.theta = }{`function(x) 2 * exp(x) / (1 + exp(x)) - 1`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta6'}{
#'              \describe{
#'             \item{hyperid = }{27006}
#'             \item{name = }{logit correlation23}
#'             \item{short.name = }{cor23}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) log((1 + x) / (1 - x))`}
#'             \item{from.theta = }{`function(x) 2 * exp(x) / (1 + exp(x)) - 1`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'iid4d'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'Gaussian random effect in dim=4 with Wishart prior'}
#'               \item{constr = }{'FALSE'}
#'               \item{nrow.ncol = }{'FALSE'}
#'               \item{augmented = }{'TRUE'}
#'               \item{aug.factor = }{'1'}
#'               \item{aug.constr = }{'1 2 3 4'}
#'               \item{n.div.by = }{'4'}
#'               \item{n.required = }{'TRUE'}
#'               \item{set.default.values = }{'TRUE'}
#'               \item{pdf = }{'iid123d'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 10.
#'          \describe{
#'           \item{Hyperparameter 'theta1'}{
#'              \describe{
#'             \item{hyperid = }{28001}
#'             \item{name = }{log precision1}
#'             \item{short.name = }{prec1}
#'             \item{initial = }{4}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{wishart4d}
#'             \item{param = }{11 1 1 1 1 0 0 0 0 0 0}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta2'}{
#'              \describe{
#'             \item{hyperid = }{28002}
#'             \item{name = }{log precision2}
#'             \item{short.name = }{prec2}
#'             \item{initial = }{4}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta3'}{
#'              \describe{
#'             \item{hyperid = }{28003}
#'             \item{name = }{log precision3}
#'             \item{short.name = }{prec3}
#'             \item{initial = }{4}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta4'}{
#'              \describe{
#'             \item{hyperid = }{28004}
#'             \item{name = }{log precision4}
#'             \item{short.name = }{prec4}
#'             \item{initial = }{4}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta5'}{
#'              \describe{
#'             \item{hyperid = }{28005}
#'             \item{name = }{logit correlation12}
#'             \item{short.name = }{cor12}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) log((1 + x) / (1 - x))`}
#'             \item{from.theta = }{`function(x) 2 * exp(x) / (1 + exp(x)) - 1`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta6'}{
#'              \describe{
#'             \item{hyperid = }{28006}
#'             \item{name = }{logit correlation13}
#'             \item{short.name = }{cor13}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) log((1 + x) / (1 - x))`}
#'             \item{from.theta = }{`function(x) 2 * exp(x) / (1 + exp(x)) - 1`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta7'}{
#'              \describe{
#'             \item{hyperid = }{28007}
#'             \item{name = }{logit correlation14}
#'             \item{short.name = }{cor14}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) log((1 + x) / (1 - x))`}
#'             \item{from.theta = }{`function(x) 2 * exp(x) / (1 + exp(x)) - 1`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta8'}{
#'              \describe{
#'             \item{hyperid = }{28008}
#'             \item{name = }{logit correlation23}
#'             \item{short.name = }{cor23}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) log((1 + x) / (1 - x))`}
#'             \item{from.theta = }{`function(x) 2 * exp(x) / (1 + exp(x)) - 1`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta9'}{
#'              \describe{
#'             \item{hyperid = }{28009}
#'             \item{name = }{logit correlation24}
#'             \item{short.name = }{cor24}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) log((1 + x) / (1 - x))`}
#'             \item{from.theta = }{`function(x) 2 * exp(x) / (1 + exp(x)) - 1`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta10'}{
#'              \describe{
#'             \item{hyperid = }{28010}
#'             \item{name = }{logit correlation34}
#'             \item{short.name = }{cor34}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) log((1 + x) / (1 - x))`}
#'             \item{from.theta = }{`function(x) 2 * exp(x) / (1 + exp(x)) - 1`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'iid5d'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'Gaussian random effect in dim=5 with Wishart prior'}
#'               \item{constr = }{'FALSE'}
#'               \item{nrow.ncol = }{'FALSE'}
#'               \item{augmented = }{'TRUE'}
#'               \item{aug.factor = }{'1'}
#'               \item{aug.constr = }{'1 2 3 4 5'}
#'               \item{n.div.by = }{'5'}
#'               \item{n.required = }{'TRUE'}
#'               \item{set.default.values = }{'TRUE'}
#'               \item{pdf = }{'iid123d'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 15.
#'          \describe{
#'           \item{Hyperparameter 'theta1'}{
#'              \describe{
#'             \item{hyperid = }{29001}
#'             \item{name = }{log precision1}
#'             \item{short.name = }{prec1}
#'             \item{initial = }{4}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{wishart5d}
#'             \item{param = }{16 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta2'}{
#'              \describe{
#'             \item{hyperid = }{29002}
#'             \item{name = }{log precision2}
#'             \item{short.name = }{prec2}
#'             \item{initial = }{4}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta3'}{
#'              \describe{
#'             \item{hyperid = }{29003}
#'             \item{name = }{log precision3}
#'             \item{short.name = }{prec3}
#'             \item{initial = }{4}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta4'}{
#'              \describe{
#'             \item{hyperid = }{29004}
#'             \item{name = }{log precision4}
#'             \item{short.name = }{prec4}
#'             \item{initial = }{4}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta5'}{
#'              \describe{
#'             \item{hyperid = }{29005}
#'             \item{name = }{log precision5}
#'             \item{short.name = }{prec5}
#'             \item{initial = }{4}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta6'}{
#'              \describe{
#'             \item{hyperid = }{29006}
#'             \item{name = }{logit correlation12}
#'             \item{short.name = }{cor12}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) log((1 + x) / (1 - x))`}
#'             \item{from.theta = }{`function(x) 2 * exp(x) / (1 + exp(x)) - 1`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta7'}{
#'              \describe{
#'             \item{hyperid = }{29007}
#'             \item{name = }{logit correlation13}
#'             \item{short.name = }{cor13}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) log((1 + x) / (1 - x))`}
#'             \item{from.theta = }{`function(x) 2 * exp(x) / (1 + exp(x)) - 1`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta8'}{
#'              \describe{
#'             \item{hyperid = }{29008}
#'             \item{name = }{logit correlation14}
#'             \item{short.name = }{cor14}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) log((1 + x) / (1 - x))`}
#'             \item{from.theta = }{`function(x) 2 * exp(x) / (1 + exp(x)) - 1`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta9'}{
#'              \describe{
#'             \item{hyperid = }{29009}
#'             \item{name = }{logit correlation15}
#'             \item{short.name = }{cor15}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) log((1 + x) / (1 - x))`}
#'             \item{from.theta = }{`function(x) 2 * exp(x) / (1 + exp(x)) - 1`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta10'}{
#'              \describe{
#'             \item{hyperid = }{29010}
#'             \item{name = }{logit correlation23}
#'             \item{short.name = }{cor23}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) log((1 + x) / (1 - x))`}
#'             \item{from.theta = }{`function(x) 2 * exp(x) / (1 + exp(x)) - 1`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta11'}{
#'              \describe{
#'             \item{hyperid = }{29011}
#'             \item{name = }{logit correlation24}
#'             \item{short.name = }{cor24}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) log((1 + x) / (1 - x))`}
#'             \item{from.theta = }{`function(x) 2 * exp(x) / (1 + exp(x)) - 1`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta12'}{
#'              \describe{
#'             \item{hyperid = }{29012}
#'             \item{name = }{logit correlation25}
#'             \item{short.name = }{cor25}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) log((1 + x) / (1 - x))`}
#'             \item{from.theta = }{`function(x) 2 * exp(x) / (1 + exp(x)) - 1`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta13'}{
#'              \describe{
#'             \item{hyperid = }{29013}
#'             \item{name = }{logit correlation34}
#'             \item{short.name = }{cor34}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) log((1 + x) / (1 - x))`}
#'             \item{from.theta = }{`function(x) 2 * exp(x) / (1 + exp(x)) - 1`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta14'}{
#'              \describe{
#'             \item{hyperid = }{29014}
#'             \item{name = }{logit correlation35}
#'             \item{short.name = }{cor35}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) log((1 + x) / (1 - x))`}
#'             \item{from.theta = }{`function(x) 2 * exp(x) / (1 + exp(x)) - 1`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta15'}{
#'              \describe{
#'             \item{hyperid = }{29015}
#'             \item{name = }{logit correlation45}
#'             \item{short.name = }{cor45}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) log((1 + x) / (1 - x))`}
#'             \item{from.theta = }{`function(x) 2 * exp(x) / (1 + exp(x)) - 1`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'iidkd'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'Gaussian random effect in dim=k with Wishart prior'}
#'               \item{constr = }{'FALSE'}
#'               \item{nrow.ncol = }{'FALSE'}
#'               \item{augmented = }{'TRUE'}
#'               \item{aug.factor = }{'1'}
#'               \item{aug.constr = }{'1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20'}
#'               \item{n.div.by = }{'-1'}
#'               \item{n.required = }{'TRUE'}
#'               \item{set.default.values = }{'TRUE'}
#'               \item{status = }{'experimental'}
#'               \item{pdf = }{'iidkd'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 210.
#'          \describe{
#'           \item{Hyperparameter 'theta1'}{
#'              \describe{
#'             \item{hyperid = }{29101}
#'             \item{name = }{theta1}
#'             \item{short.name = }{theta1}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{wishartkd}
#'             \item{param = }{21 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576 1048576}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta2'}{
#'              \describe{
#'             \item{hyperid = }{29102}
#'             \item{name = }{theta2}
#'             \item{short.name = }{theta2}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta3'}{
#'              \describe{
#'             \item{hyperid = }{29103}
#'             \item{name = }{theta3}
#'             \item{short.name = }{theta3}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta4'}{
#'              \describe{
#'             \item{hyperid = }{29104}
#'             \item{name = }{theta4}
#'             \item{short.name = }{theta4}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta5'}{
#'              \describe{
#'             \item{hyperid = }{29105}
#'             \item{name = }{theta5}
#'             \item{short.name = }{theta5}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta6'}{
#'              \describe{
#'             \item{hyperid = }{29106}
#'             \item{name = }{theta6}
#'             \item{short.name = }{theta6}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta7'}{
#'              \describe{
#'             \item{hyperid = }{29107}
#'             \item{name = }{theta7}
#'             \item{short.name = }{theta7}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta8'}{
#'              \describe{
#'             \item{hyperid = }{29108}
#'             \item{name = }{theta8}
#'             \item{short.name = }{theta8}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta9'}{
#'              \describe{
#'             \item{hyperid = }{29109}
#'             \item{name = }{theta9}
#'             \item{short.name = }{theta9}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta10'}{
#'              \describe{
#'             \item{hyperid = }{29110}
#'             \item{name = }{theta10}
#'             \item{short.name = }{theta10}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta11'}{
#'              \describe{
#'             \item{hyperid = }{29111}
#'             \item{name = }{theta11}
#'             \item{short.name = }{theta11}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta12'}{
#'              \describe{
#'             \item{hyperid = }{29112}
#'             \item{name = }{theta12}
#'             \item{short.name = }{theta12}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta13'}{
#'              \describe{
#'             \item{hyperid = }{29113}
#'             \item{name = }{theta13}
#'             \item{short.name = }{theta13}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta14'}{
#'              \describe{
#'             \item{hyperid = }{29114}
#'             \item{name = }{theta14}
#'             \item{short.name = }{theta14}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta15'}{
#'              \describe{
#'             \item{hyperid = }{29115}
#'             \item{name = }{theta15}
#'             \item{short.name = }{theta15}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta16'}{
#'              \describe{
#'             \item{hyperid = }{29116}
#'             \item{name = }{theta16}
#'             \item{short.name = }{theta16}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta17'}{
#'              \describe{
#'             \item{hyperid = }{29117}
#'             \item{name = }{theta17}
#'             \item{short.name = }{theta17}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta18'}{
#'              \describe{
#'             \item{hyperid = }{29118}
#'             \item{name = }{theta18}
#'             \item{short.name = }{theta18}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta19'}{
#'              \describe{
#'             \item{hyperid = }{29119}
#'             \item{name = }{theta19}
#'             \item{short.name = }{theta19}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta20'}{
#'              \describe{
#'             \item{hyperid = }{29120}
#'             \item{name = }{theta20}
#'             \item{short.name = }{theta20}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta21'}{
#'              \describe{
#'             \item{hyperid = }{29121}
#'             \item{name = }{theta21}
#'             \item{short.name = }{theta21}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta22'}{
#'              \describe{
#'             \item{hyperid = }{29122}
#'             \item{name = }{theta22}
#'             \item{short.name = }{theta22}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta23'}{
#'              \describe{
#'             \item{hyperid = }{29123}
#'             \item{name = }{theta23}
#'             \item{short.name = }{theta23}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta24'}{
#'              \describe{
#'             \item{hyperid = }{29124}
#'             \item{name = }{theta24}
#'             \item{short.name = }{theta24}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta25'}{
#'              \describe{
#'             \item{hyperid = }{29125}
#'             \item{name = }{theta25}
#'             \item{short.name = }{theta25}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta26'}{
#'              \describe{
#'             \item{hyperid = }{29126}
#'             \item{name = }{theta26}
#'             \item{short.name = }{theta26}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta27'}{
#'              \describe{
#'             \item{hyperid = }{29127}
#'             \item{name = }{theta27}
#'             \item{short.name = }{theta27}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta28'}{
#'              \describe{
#'             \item{hyperid = }{29128}
#'             \item{name = }{theta28}
#'             \item{short.name = }{theta28}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta29'}{
#'              \describe{
#'             \item{hyperid = }{29129}
#'             \item{name = }{theta29}
#'             \item{short.name = }{theta29}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta30'}{
#'              \describe{
#'             \item{hyperid = }{29130}
#'             \item{name = }{theta30}
#'             \item{short.name = }{theta30}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta31'}{
#'              \describe{
#'             \item{hyperid = }{29131}
#'             \item{name = }{theta31}
#'             \item{short.name = }{theta31}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta32'}{
#'              \describe{
#'             \item{hyperid = }{29132}
#'             \item{name = }{theta32}
#'             \item{short.name = }{theta32}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta33'}{
#'              \describe{
#'             \item{hyperid = }{29133}
#'             \item{name = }{theta33}
#'             \item{short.name = }{theta33}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta34'}{
#'              \describe{
#'             \item{hyperid = }{29134}
#'             \item{name = }{theta34}
#'             \item{short.name = }{theta34}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta35'}{
#'              \describe{
#'             \item{hyperid = }{29135}
#'             \item{name = }{theta35}
#'             \item{short.name = }{theta35}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta36'}{
#'              \describe{
#'             \item{hyperid = }{29136}
#'             \item{name = }{theta36}
#'             \item{short.name = }{theta36}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta37'}{
#'              \describe{
#'             \item{hyperid = }{29137}
#'             \item{name = }{theta37}
#'             \item{short.name = }{theta37}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta38'}{
#'              \describe{
#'             \item{hyperid = }{29138}
#'             \item{name = }{theta38}
#'             \item{short.name = }{theta38}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta39'}{
#'              \describe{
#'             \item{hyperid = }{29139}
#'             \item{name = }{theta39}
#'             \item{short.name = }{theta39}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta40'}{
#'              \describe{
#'             \item{hyperid = }{29140}
#'             \item{name = }{theta40}
#'             \item{short.name = }{theta40}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta41'}{
#'              \describe{
#'             \item{hyperid = }{29141}
#'             \item{name = }{theta41}
#'             \item{short.name = }{theta41}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta42'}{
#'              \describe{
#'             \item{hyperid = }{29142}
#'             \item{name = }{theta42}
#'             \item{short.name = }{theta42}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta43'}{
#'              \describe{
#'             \item{hyperid = }{29143}
#'             \item{name = }{theta43}
#'             \item{short.name = }{theta43}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta44'}{
#'              \describe{
#'             \item{hyperid = }{29144}
#'             \item{name = }{theta44}
#'             \item{short.name = }{theta44}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta45'}{
#'              \describe{
#'             \item{hyperid = }{29145}
#'             \item{name = }{theta45}
#'             \item{short.name = }{theta45}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta46'}{
#'              \describe{
#'             \item{hyperid = }{29146}
#'             \item{name = }{theta46}
#'             \item{short.name = }{theta46}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta47'}{
#'              \describe{
#'             \item{hyperid = }{29147}
#'             \item{name = }{theta47}
#'             \item{short.name = }{theta47}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta48'}{
#'              \describe{
#'             \item{hyperid = }{29148}
#'             \item{name = }{theta48}
#'             \item{short.name = }{theta48}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta49'}{
#'              \describe{
#'             \item{hyperid = }{29149}
#'             \item{name = }{theta49}
#'             \item{short.name = }{theta49}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta50'}{
#'              \describe{
#'             \item{hyperid = }{29150}
#'             \item{name = }{theta50}
#'             \item{short.name = }{theta50}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta51'}{
#'              \describe{
#'             \item{hyperid = }{29151}
#'             \item{name = }{theta51}
#'             \item{short.name = }{theta51}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta52'}{
#'              \describe{
#'             \item{hyperid = }{29152}
#'             \item{name = }{theta52}
#'             \item{short.name = }{theta52}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta53'}{
#'              \describe{
#'             \item{hyperid = }{29153}
#'             \item{name = }{theta53}
#'             \item{short.name = }{theta53}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta54'}{
#'              \describe{
#'             \item{hyperid = }{29154}
#'             \item{name = }{theta54}
#'             \item{short.name = }{theta54}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta55'}{
#'              \describe{
#'             \item{hyperid = }{29155}
#'             \item{name = }{theta55}
#'             \item{short.name = }{theta55}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta56'}{
#'              \describe{
#'             \item{hyperid = }{29156}
#'             \item{name = }{theta56}
#'             \item{short.name = }{theta56}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta57'}{
#'              \describe{
#'             \item{hyperid = }{29157}
#'             \item{name = }{theta57}
#'             \item{short.name = }{theta57}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta58'}{
#'              \describe{
#'             \item{hyperid = }{29158}
#'             \item{name = }{theta58}
#'             \item{short.name = }{theta58}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta59'}{
#'              \describe{
#'             \item{hyperid = }{29159}
#'             \item{name = }{theta59}
#'             \item{short.name = }{theta59}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta60'}{
#'              \describe{
#'             \item{hyperid = }{29160}
#'             \item{name = }{theta60}
#'             \item{short.name = }{theta60}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta61'}{
#'              \describe{
#'             \item{hyperid = }{29161}
#'             \item{name = }{theta61}
#'             \item{short.name = }{theta61}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta62'}{
#'              \describe{
#'             \item{hyperid = }{29162}
#'             \item{name = }{theta62}
#'             \item{short.name = }{theta62}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta63'}{
#'              \describe{
#'             \item{hyperid = }{29163}
#'             \item{name = }{theta63}
#'             \item{short.name = }{theta63}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta64'}{
#'              \describe{
#'             \item{hyperid = }{29164}
#'             \item{name = }{theta64}
#'             \item{short.name = }{theta64}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta65'}{
#'              \describe{
#'             \item{hyperid = }{29165}
#'             \item{name = }{theta65}
#'             \item{short.name = }{theta65}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta66'}{
#'              \describe{
#'             \item{hyperid = }{29166}
#'             \item{name = }{theta66}
#'             \item{short.name = }{theta66}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta67'}{
#'              \describe{
#'             \item{hyperid = }{29167}
#'             \item{name = }{theta67}
#'             \item{short.name = }{theta67}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta68'}{
#'              \describe{
#'             \item{hyperid = }{29168}
#'             \item{name = }{theta68}
#'             \item{short.name = }{theta68}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta69'}{
#'              \describe{
#'             \item{hyperid = }{29169}
#'             \item{name = }{theta69}
#'             \item{short.name = }{theta69}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta70'}{
#'              \describe{
#'             \item{hyperid = }{29170}
#'             \item{name = }{theta70}
#'             \item{short.name = }{theta70}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta71'}{
#'              \describe{
#'             \item{hyperid = }{29171}
#'             \item{name = }{theta71}
#'             \item{short.name = }{theta71}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta72'}{
#'              \describe{
#'             \item{hyperid = }{29172}
#'             \item{name = }{theta72}
#'             \item{short.name = }{theta72}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta73'}{
#'              \describe{
#'             \item{hyperid = }{29173}
#'             \item{name = }{theta73}
#'             \item{short.name = }{theta73}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta74'}{
#'              \describe{
#'             \item{hyperid = }{29174}
#'             \item{name = }{theta74}
#'             \item{short.name = }{theta74}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta75'}{
#'              \describe{
#'             \item{hyperid = }{29175}
#'             \item{name = }{theta75}
#'             \item{short.name = }{theta75}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta76'}{
#'              \describe{
#'             \item{hyperid = }{29176}
#'             \item{name = }{theta76}
#'             \item{short.name = }{theta76}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta77'}{
#'              \describe{
#'             \item{hyperid = }{29177}
#'             \item{name = }{theta77}
#'             \item{short.name = }{theta77}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta78'}{
#'              \describe{
#'             \item{hyperid = }{29178}
#'             \item{name = }{theta78}
#'             \item{short.name = }{theta78}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta79'}{
#'              \describe{
#'             \item{hyperid = }{29179}
#'             \item{name = }{theta79}
#'             \item{short.name = }{theta79}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta80'}{
#'              \describe{
#'             \item{hyperid = }{29180}
#'             \item{name = }{theta80}
#'             \item{short.name = }{theta80}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta81'}{
#'              \describe{
#'             \item{hyperid = }{29181}
#'             \item{name = }{theta81}
#'             \item{short.name = }{theta81}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta82'}{
#'              \describe{
#'             \item{hyperid = }{29182}
#'             \item{name = }{theta82}
#'             \item{short.name = }{theta82}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta83'}{
#'              \describe{
#'             \item{hyperid = }{29183}
#'             \item{name = }{theta83}
#'             \item{short.name = }{theta83}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta84'}{
#'              \describe{
#'             \item{hyperid = }{29184}
#'             \item{name = }{theta84}
#'             \item{short.name = }{theta84}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta85'}{
#'              \describe{
#'             \item{hyperid = }{29185}
#'             \item{name = }{theta85}
#'             \item{short.name = }{theta85}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta86'}{
#'              \describe{
#'             \item{hyperid = }{29186}
#'             \item{name = }{theta86}
#'             \item{short.name = }{theta86}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta87'}{
#'              \describe{
#'             \item{hyperid = }{29187}
#'             \item{name = }{theta87}
#'             \item{short.name = }{theta87}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta88'}{
#'              \describe{
#'             \item{hyperid = }{29188}
#'             \item{name = }{theta88}
#'             \item{short.name = }{theta88}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta89'}{
#'              \describe{
#'             \item{hyperid = }{29189}
#'             \item{name = }{theta89}
#'             \item{short.name = }{theta89}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta90'}{
#'              \describe{
#'             \item{hyperid = }{29190}
#'             \item{name = }{theta90}
#'             \item{short.name = }{theta90}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta91'}{
#'              \describe{
#'             \item{hyperid = }{29191}
#'             \item{name = }{theta91}
#'             \item{short.name = }{theta91}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta92'}{
#'              \describe{
#'             \item{hyperid = }{29192}
#'             \item{name = }{theta92}
#'             \item{short.name = }{theta92}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta93'}{
#'              \describe{
#'             \item{hyperid = }{29193}
#'             \item{name = }{theta93}
#'             \item{short.name = }{theta93}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta94'}{
#'              \describe{
#'             \item{hyperid = }{29194}
#'             \item{name = }{theta94}
#'             \item{short.name = }{theta94}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta95'}{
#'              \describe{
#'             \item{hyperid = }{29195}
#'             \item{name = }{theta95}
#'             \item{short.name = }{theta95}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta96'}{
#'              \describe{
#'             \item{hyperid = }{29196}
#'             \item{name = }{theta96}
#'             \item{short.name = }{theta96}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta97'}{
#'              \describe{
#'             \item{hyperid = }{29197}
#'             \item{name = }{theta97}
#'             \item{short.name = }{theta97}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta98'}{
#'              \describe{
#'             \item{hyperid = }{29198}
#'             \item{name = }{theta98}
#'             \item{short.name = }{theta98}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta99'}{
#'              \describe{
#'             \item{hyperid = }{29199}
#'             \item{name = }{theta99}
#'             \item{short.name = }{theta99}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta100'}{
#'              \describe{
#'             \item{hyperid = }{29200}
#'             \item{name = }{theta100}
#'             \item{short.name = }{theta100}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta101'}{
#'              \describe{
#'             \item{hyperid = }{29201}
#'             \item{name = }{theta101}
#'             \item{short.name = }{theta101}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta102'}{
#'              \describe{
#'             \item{hyperid = }{29202}
#'             \item{name = }{theta102}
#'             \item{short.name = }{theta102}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta103'}{
#'              \describe{
#'             \item{hyperid = }{29203}
#'             \item{name = }{theta103}
#'             \item{short.name = }{theta103}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta104'}{
#'              \describe{
#'             \item{hyperid = }{29204}
#'             \item{name = }{theta104}
#'             \item{short.name = }{theta104}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta105'}{
#'              \describe{
#'             \item{hyperid = }{29205}
#'             \item{name = }{theta105}
#'             \item{short.name = }{theta105}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta106'}{
#'              \describe{
#'             \item{hyperid = }{29206}
#'             \item{name = }{theta106}
#'             \item{short.name = }{theta106}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta107'}{
#'              \describe{
#'             \item{hyperid = }{29207}
#'             \item{name = }{theta107}
#'             \item{short.name = }{theta107}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta108'}{
#'              \describe{
#'             \item{hyperid = }{29208}
#'             \item{name = }{theta108}
#'             \item{short.name = }{theta108}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta109'}{
#'              \describe{
#'             \item{hyperid = }{29209}
#'             \item{name = }{theta109}
#'             \item{short.name = }{theta109}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta110'}{
#'              \describe{
#'             \item{hyperid = }{29210}
#'             \item{name = }{theta110}
#'             \item{short.name = }{theta110}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta111'}{
#'              \describe{
#'             \item{hyperid = }{29211}
#'             \item{name = }{theta111}
#'             \item{short.name = }{theta111}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta112'}{
#'              \describe{
#'             \item{hyperid = }{29212}
#'             \item{name = }{theta112}
#'             \item{short.name = }{theta112}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta113'}{
#'              \describe{
#'             \item{hyperid = }{29213}
#'             \item{name = }{theta113}
#'             \item{short.name = }{theta113}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta114'}{
#'              \describe{
#'             \item{hyperid = }{29214}
#'             \item{name = }{theta114}
#'             \item{short.name = }{theta114}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta115'}{
#'              \describe{
#'             \item{hyperid = }{29215}
#'             \item{name = }{theta115}
#'             \item{short.name = }{theta115}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta116'}{
#'              \describe{
#'             \item{hyperid = }{29216}
#'             \item{name = }{theta116}
#'             \item{short.name = }{theta116}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta117'}{
#'              \describe{
#'             \item{hyperid = }{29217}
#'             \item{name = }{theta117}
#'             \item{short.name = }{theta117}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta118'}{
#'              \describe{
#'             \item{hyperid = }{29218}
#'             \item{name = }{theta118}
#'             \item{short.name = }{theta118}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta119'}{
#'              \describe{
#'             \item{hyperid = }{29219}
#'             \item{name = }{theta119}
#'             \item{short.name = }{theta119}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta120'}{
#'              \describe{
#'             \item{hyperid = }{29220}
#'             \item{name = }{theta120}
#'             \item{short.name = }{theta120}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta121'}{
#'              \describe{
#'             \item{hyperid = }{29221}
#'             \item{name = }{theta121}
#'             \item{short.name = }{theta121}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta122'}{
#'              \describe{
#'             \item{hyperid = }{29222}
#'             \item{name = }{theta122}
#'             \item{short.name = }{theta122}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta123'}{
#'              \describe{
#'             \item{hyperid = }{29223}
#'             \item{name = }{theta123}
#'             \item{short.name = }{theta123}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta124'}{
#'              \describe{
#'             \item{hyperid = }{29224}
#'             \item{name = }{theta124}
#'             \item{short.name = }{theta124}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta125'}{
#'              \describe{
#'             \item{hyperid = }{29225}
#'             \item{name = }{theta125}
#'             \item{short.name = }{theta125}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta126'}{
#'              \describe{
#'             \item{hyperid = }{29226}
#'             \item{name = }{theta126}
#'             \item{short.name = }{theta126}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta127'}{
#'              \describe{
#'             \item{hyperid = }{29227}
#'             \item{name = }{theta127}
#'             \item{short.name = }{theta127}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta128'}{
#'              \describe{
#'             \item{hyperid = }{29228}
#'             \item{name = }{theta128}
#'             \item{short.name = }{theta128}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta129'}{
#'              \describe{
#'             \item{hyperid = }{29229}
#'             \item{name = }{theta129}
#'             \item{short.name = }{theta129}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta130'}{
#'              \describe{
#'             \item{hyperid = }{29230}
#'             \item{name = }{theta130}
#'             \item{short.name = }{theta130}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta131'}{
#'              \describe{
#'             \item{hyperid = }{29231}
#'             \item{name = }{theta131}
#'             \item{short.name = }{theta131}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta132'}{
#'              \describe{
#'             \item{hyperid = }{29232}
#'             \item{name = }{theta132}
#'             \item{short.name = }{theta132}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta133'}{
#'              \describe{
#'             \item{hyperid = }{29233}
#'             \item{name = }{theta133}
#'             \item{short.name = }{theta133}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta134'}{
#'              \describe{
#'             \item{hyperid = }{29234}
#'             \item{name = }{theta134}
#'             \item{short.name = }{theta134}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta135'}{
#'              \describe{
#'             \item{hyperid = }{29235}
#'             \item{name = }{theta135}
#'             \item{short.name = }{theta135}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta136'}{
#'              \describe{
#'             \item{hyperid = }{29236}
#'             \item{name = }{theta136}
#'             \item{short.name = }{theta136}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta137'}{
#'              \describe{
#'             \item{hyperid = }{29237}
#'             \item{name = }{theta137}
#'             \item{short.name = }{theta137}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta138'}{
#'              \describe{
#'             \item{hyperid = }{29238}
#'             \item{name = }{theta138}
#'             \item{short.name = }{theta138}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta139'}{
#'              \describe{
#'             \item{hyperid = }{29239}
#'             \item{name = }{theta139}
#'             \item{short.name = }{theta139}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta140'}{
#'              \describe{
#'             \item{hyperid = }{29240}
#'             \item{name = }{theta140}
#'             \item{short.name = }{theta140}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta141'}{
#'              \describe{
#'             \item{hyperid = }{29241}
#'             \item{name = }{theta141}
#'             \item{short.name = }{theta141}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta142'}{
#'              \describe{
#'             \item{hyperid = }{29242}
#'             \item{name = }{theta142}
#'             \item{short.name = }{theta142}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta143'}{
#'              \describe{
#'             \item{hyperid = }{29243}
#'             \item{name = }{theta143}
#'             \item{short.name = }{theta143}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta144'}{
#'              \describe{
#'             \item{hyperid = }{29244}
#'             \item{name = }{theta144}
#'             \item{short.name = }{theta144}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta145'}{
#'              \describe{
#'             \item{hyperid = }{29245}
#'             \item{name = }{theta145}
#'             \item{short.name = }{theta145}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta146'}{
#'              \describe{
#'             \item{hyperid = }{29246}
#'             \item{name = }{theta146}
#'             \item{short.name = }{theta146}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta147'}{
#'              \describe{
#'             \item{hyperid = }{29247}
#'             \item{name = }{theta147}
#'             \item{short.name = }{theta147}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta148'}{
#'              \describe{
#'             \item{hyperid = }{29248}
#'             \item{name = }{theta148}
#'             \item{short.name = }{theta148}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta149'}{
#'              \describe{
#'             \item{hyperid = }{29249}
#'             \item{name = }{theta149}
#'             \item{short.name = }{theta149}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta150'}{
#'              \describe{
#'             \item{hyperid = }{29250}
#'             \item{name = }{theta150}
#'             \item{short.name = }{theta150}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta151'}{
#'              \describe{
#'             \item{hyperid = }{29251}
#'             \item{name = }{theta151}
#'             \item{short.name = }{theta151}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta152'}{
#'              \describe{
#'             \item{hyperid = }{29252}
#'             \item{name = }{theta152}
#'             \item{short.name = }{theta152}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta153'}{
#'              \describe{
#'             \item{hyperid = }{29253}
#'             \item{name = }{theta153}
#'             \item{short.name = }{theta153}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta154'}{
#'              \describe{
#'             \item{hyperid = }{29254}
#'             \item{name = }{theta154}
#'             \item{short.name = }{theta154}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta155'}{
#'              \describe{
#'             \item{hyperid = }{29255}
#'             \item{name = }{theta155}
#'             \item{short.name = }{theta155}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta156'}{
#'              \describe{
#'             \item{hyperid = }{29256}
#'             \item{name = }{theta156}
#'             \item{short.name = }{theta156}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta157'}{
#'              \describe{
#'             \item{hyperid = }{29257}
#'             \item{name = }{theta157}
#'             \item{short.name = }{theta157}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta158'}{
#'              \describe{
#'             \item{hyperid = }{29258}
#'             \item{name = }{theta158}
#'             \item{short.name = }{theta158}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta159'}{
#'              \describe{
#'             \item{hyperid = }{29259}
#'             \item{name = }{theta159}
#'             \item{short.name = }{theta159}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta160'}{
#'              \describe{
#'             \item{hyperid = }{29260}
#'             \item{name = }{theta160}
#'             \item{short.name = }{theta160}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta161'}{
#'              \describe{
#'             \item{hyperid = }{29261}
#'             \item{name = }{theta161}
#'             \item{short.name = }{theta161}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta162'}{
#'              \describe{
#'             \item{hyperid = }{29262}
#'             \item{name = }{theta162}
#'             \item{short.name = }{theta162}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta163'}{
#'              \describe{
#'             \item{hyperid = }{29263}
#'             \item{name = }{theta163}
#'             \item{short.name = }{theta163}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta164'}{
#'              \describe{
#'             \item{hyperid = }{29264}
#'             \item{name = }{theta164}
#'             \item{short.name = }{theta164}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta165'}{
#'              \describe{
#'             \item{hyperid = }{29265}
#'             \item{name = }{theta165}
#'             \item{short.name = }{theta165}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta166'}{
#'              \describe{
#'             \item{hyperid = }{29266}
#'             \item{name = }{theta166}
#'             \item{short.name = }{theta166}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta167'}{
#'              \describe{
#'             \item{hyperid = }{29267}
#'             \item{name = }{theta167}
#'             \item{short.name = }{theta167}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta168'}{
#'              \describe{
#'             \item{hyperid = }{29268}
#'             \item{name = }{theta168}
#'             \item{short.name = }{theta168}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta169'}{
#'              \describe{
#'             \item{hyperid = }{29269}
#'             \item{name = }{theta169}
#'             \item{short.name = }{theta169}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta170'}{
#'              \describe{
#'             \item{hyperid = }{29270}
#'             \item{name = }{theta170}
#'             \item{short.name = }{theta170}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta171'}{
#'              \describe{
#'             \item{hyperid = }{29271}
#'             \item{name = }{theta171}
#'             \item{short.name = }{theta171}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta172'}{
#'              \describe{
#'             \item{hyperid = }{29272}
#'             \item{name = }{theta172}
#'             \item{short.name = }{theta172}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta173'}{
#'              \describe{
#'             \item{hyperid = }{29273}
#'             \item{name = }{theta173}
#'             \item{short.name = }{theta173}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta174'}{
#'              \describe{
#'             \item{hyperid = }{29274}
#'             \item{name = }{theta174}
#'             \item{short.name = }{theta174}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta175'}{
#'              \describe{
#'             \item{hyperid = }{29275}
#'             \item{name = }{theta175}
#'             \item{short.name = }{theta175}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta176'}{
#'              \describe{
#'             \item{hyperid = }{29276}
#'             \item{name = }{theta176}
#'             \item{short.name = }{theta176}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta177'}{
#'              \describe{
#'             \item{hyperid = }{29277}
#'             \item{name = }{theta177}
#'             \item{short.name = }{theta177}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta178'}{
#'              \describe{
#'             \item{hyperid = }{29278}
#'             \item{name = }{theta178}
#'             \item{short.name = }{theta178}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta179'}{
#'              \describe{
#'             \item{hyperid = }{29279}
#'             \item{name = }{theta179}
#'             \item{short.name = }{theta179}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta180'}{
#'              \describe{
#'             \item{hyperid = }{29280}
#'             \item{name = }{theta180}
#'             \item{short.name = }{theta180}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta181'}{
#'              \describe{
#'             \item{hyperid = }{29281}
#'             \item{name = }{theta181}
#'             \item{short.name = }{theta181}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta182'}{
#'              \describe{
#'             \item{hyperid = }{29282}
#'             \item{name = }{theta182}
#'             \item{short.name = }{theta182}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta183'}{
#'              \describe{
#'             \item{hyperid = }{29283}
#'             \item{name = }{theta183}
#'             \item{short.name = }{theta183}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta184'}{
#'              \describe{
#'             \item{hyperid = }{29284}
#'             \item{name = }{theta184}
#'             \item{short.name = }{theta184}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta185'}{
#'              \describe{
#'             \item{hyperid = }{29285}
#'             \item{name = }{theta185}
#'             \item{short.name = }{theta185}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta186'}{
#'              \describe{
#'             \item{hyperid = }{29286}
#'             \item{name = }{theta186}
#'             \item{short.name = }{theta186}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta187'}{
#'              \describe{
#'             \item{hyperid = }{29287}
#'             \item{name = }{theta187}
#'             \item{short.name = }{theta187}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta188'}{
#'              \describe{
#'             \item{hyperid = }{29288}
#'             \item{name = }{theta188}
#'             \item{short.name = }{theta188}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta189'}{
#'              \describe{
#'             \item{hyperid = }{29289}
#'             \item{name = }{theta189}
#'             \item{short.name = }{theta189}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta190'}{
#'              \describe{
#'             \item{hyperid = }{29290}
#'             \item{name = }{theta190}
#'             \item{short.name = }{theta190}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta191'}{
#'              \describe{
#'             \item{hyperid = }{29291}
#'             \item{name = }{theta191}
#'             \item{short.name = }{theta191}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta192'}{
#'              \describe{
#'             \item{hyperid = }{29292}
#'             \item{name = }{theta192}
#'             \item{short.name = }{theta192}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta193'}{
#'              \describe{
#'             \item{hyperid = }{29293}
#'             \item{name = }{theta193}
#'             \item{short.name = }{theta193}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta194'}{
#'              \describe{
#'             \item{hyperid = }{29294}
#'             \item{name = }{theta194}
#'             \item{short.name = }{theta194}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta195'}{
#'              \describe{
#'             \item{hyperid = }{29295}
#'             \item{name = }{theta195}
#'             \item{short.name = }{theta195}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta196'}{
#'              \describe{
#'             \item{hyperid = }{29296}
#'             \item{name = }{theta196}
#'             \item{short.name = }{theta196}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta197'}{
#'              \describe{
#'             \item{hyperid = }{29297}
#'             \item{name = }{theta197}
#'             \item{short.name = }{theta197}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta198'}{
#'              \describe{
#'             \item{hyperid = }{29298}
#'             \item{name = }{theta198}
#'             \item{short.name = }{theta198}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta199'}{
#'              \describe{
#'             \item{hyperid = }{29299}
#'             \item{name = }{theta199}
#'             \item{short.name = }{theta199}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta200'}{
#'              \describe{
#'             \item{hyperid = }{29300}
#'             \item{name = }{theta200}
#'             \item{short.name = }{theta200}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta201'}{
#'              \describe{
#'             \item{hyperid = }{29301}
#'             \item{name = }{theta201}
#'             \item{short.name = }{theta201}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta202'}{
#'              \describe{
#'             \item{hyperid = }{29302}
#'             \item{name = }{theta202}
#'             \item{short.name = }{theta202}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta203'}{
#'              \describe{
#'             \item{hyperid = }{29303}
#'             \item{name = }{theta203}
#'             \item{short.name = }{theta203}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta204'}{
#'              \describe{
#'             \item{hyperid = }{29304}
#'             \item{name = }{theta204}
#'             \item{short.name = }{theta204}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta205'}{
#'              \describe{
#'             \item{hyperid = }{29305}
#'             \item{name = }{theta205}
#'             \item{short.name = }{theta205}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta206'}{
#'              \describe{
#'             \item{hyperid = }{29306}
#'             \item{name = }{theta206}
#'             \item{short.name = }{theta206}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta207'}{
#'              \describe{
#'             \item{hyperid = }{29307}
#'             \item{name = }{theta207}
#'             \item{short.name = }{theta207}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta208'}{
#'              \describe{
#'             \item{hyperid = }{29308}
#'             \item{name = }{theta208}
#'             \item{short.name = }{theta208}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta209'}{
#'              \describe{
#'             \item{hyperid = }{29309}
#'             \item{name = }{theta209}
#'             \item{short.name = }{theta209}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta210'}{
#'              \describe{
#'             \item{hyperid = }{29310}
#'             \item{name = }{theta210}
#'             \item{short.name = }{theta210}
#'             \item{initial = }{1048576}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model '2diid'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'(This model is obsolute)'}
#'               \item{constr = }{'FALSE'}
#'               \item{nrow.ncol = }{'FALSE'}
#'               \item{augmented = }{'FALSE'}
#'               \item{aug.factor = }{'1'}
#'               \item{aug.constr = }{'1 2'}
#'               \item{n.div.by = }{'2'}
#'               \item{n.required = }{'TRUE'}
#'               \item{set.default.values = }{'TRUE'}
#'               \item{pdf = }{'iid123d'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 3.
#'          \describe{
#'           \item{Hyperparameter 'theta1'}{
#'              \describe{
#'             \item{hyperid = }{30001}
#'             \item{name = }{log precision1}
#'             \item{short.name = }{prec1}
#'             \item{initial = }{4}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{loggamma}
#'             \item{param = }{1 5e-05}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta2'}{
#'              \describe{
#'             \item{hyperid = }{30002}
#'             \item{name = }{log precision2}
#'             \item{short.name = }{prec2}
#'             \item{initial = }{4}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{loggamma}
#'             \item{param = }{1 5e-05}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta3'}{
#'              \describe{
#'             \item{hyperid = }{30003}
#'             \item{name = }{correlation}
#'             \item{short.name = }{cor}
#'             \item{initial = }{4}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 0.15}
#'             \item{to.theta = }{`function(x) log((1 + x) / (1 - x))`}
#'             \item{from.theta = }{`function(x) 2 * exp(x) / (1 + exp(x)) - 1`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'z'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'The z-model in a classical mixed model formulation'}
#'               \item{constr = }{'FALSE'}
#'               \item{nrow.ncol = }{'FALSE'}
#'               \item{augmented = }{'FALSE'}
#'               \item{aug.factor = }{'1'}
#'               \item{aug.constr = }{'NULL'}
#'               \item{n.div.by = }{'NULL'}
#'               \item{n.required = }{'TRUE'}
#'               \item{set.default.values = }{'TRUE'}
#'               \item{pdf = }{'z'}
#'               \item{status = }{'experimental'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 1.
#'          \describe{
#'           \item{Hyperparameter 'theta'}{
#'              \describe{
#'             \item{hyperid = }{31001}
#'             \item{name = }{log precision}
#'             \item{short.name = }{prec}
#'             \item{initial = }{4}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{loggamma}
#'             \item{param = }{1 5e-05}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'rw2d'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'Thin-plate spline model'}
#'               \item{constr = }{'TRUE'}
#'               \item{nrow.ncol = }{'TRUE'}
#'               \item{augmented = }{'FALSE'}
#'               \item{aug.factor = }{'1'}
#'               \item{aug.constr = }{'NULL'}
#'               \item{n.div.by = }{'NULL'}
#'               \item{n.required = }{'FALSE'}
#'               \item{set.default.values = }{'TRUE'}
#'               \item{pdf = }{'rw2d'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 1.
#'          \describe{
#'           \item{Hyperparameter 'theta'}{
#'              \describe{
#'             \item{hyperid = }{32001}
#'             \item{name = }{log precision}
#'             \item{short.name = }{prec}
#'             \item{initial = }{4}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{loggamma}
#'             \item{param = }{1 5e-05}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'rw2diid'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'Thin-plate spline with iid noise'}
#'               \item{constr = }{'TRUE'}
#'               \item{nrow.ncol = }{'TRUE'}
#'               \item{augmented = }{'TRUE'}
#'               \item{aug.factor = }{'2'}
#'               \item{aug.constr = }{'2'}
#'               \item{n.div.by = }{'NULL'}
#'               \item{n.required = }{'FALSE'}
#'               \item{set.default.values = }{'TRUE'}
#'               \item{status = }{'experimental'}
#'               \item{pdf = }{'rw2diid'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 2.
#'          \describe{
#'           \item{Hyperparameter 'theta1'}{
#'              \describe{
#'             \item{hyperid = }{33001}
#'             \item{name = }{log precision}
#'             \item{short.name = }{prec}
#'             \item{prior = }{pc.prec}
#'             \item{param = }{1 0.01}
#'             \item{initial = }{4}
#'             \item{fixed = }{FALSE}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta2'}{
#'              \describe{
#'             \item{hyperid = }{33002}
#'             \item{name = }{logit phi}
#'             \item{short.name = }{phi}
#'             \item{prior = }{pc}
#'             \item{param = }{0.5 0.5}
#'             \item{initial = }{3}
#'             \item{fixed = }{FALSE}
#'             \item{to.theta = }{`function(x) log(x / (1 - x))`}
#'             \item{from.theta = }{`function(x) exp(x) / (1 + exp(x))`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'slm'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'Spatial lag model'}
#'               \item{constr = }{'FALSE'}
#'               \item{nrow.ncol = }{'FALSE'}
#'               \item{augmented = }{'FALSE'}
#'               \item{aug.factor = }{'1'}
#'               \item{aug.constr = }{'NULL'}
#'               \item{n.div.by = }{'NULL'}
#'               \item{n.required = }{'TRUE'}
#'               \item{set.default.values = }{'TRUE'}
#'               \item{pdf = }{'slm'}
#'               \item{status = }{'experimental'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 2.
#'          \describe{
#'           \item{Hyperparameter 'theta1'}{
#'              \describe{
#'             \item{hyperid = }{34001}
#'             \item{name = }{log precision}
#'             \item{short.name = }{prec}
#'             \item{initial = }{4}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{loggamma}
#'             \item{param = }{1 5e-05}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta2'}{
#'              \describe{
#'             \item{hyperid = }{34002}
#'             \item{name = }{rho}
#'             \item{short.name = }{rho}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 10}
#'             \item{to.theta = }{`function(x) log(x / (1 - x))`}
#'             \item{from.theta = }{`function(x) 1 / (1 + exp(-x))`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'matern2d'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'Matern covariance function on a regular grid'}
#'               \item{constr = }{'FALSE'}
#'               \item{nrow.ncol = }{'TRUE'}
#'               \item{augmented = }{'FALSE'}
#'               \item{aug.factor = }{'1'}
#'               \item{aug.constr = }{'NULL'}
#'               \item{n.div.by = }{'NULL'}
#'               \item{n.required = }{'FALSE'}
#'               \item{set.default.values = }{'TRUE'}
#'               \item{pdf = }{'matern2d'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 2.
#'          \describe{
#'           \item{Hyperparameter 'theta1'}{
#'              \describe{
#'             \item{hyperid = }{35001}
#'             \item{name = }{log precision}
#'             \item{short.name = }{prec}
#'             \item{initial = }{4}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{loggamma}
#'             \item{param = }{1 5e-05}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta2'}{
#'              \describe{
#'             \item{hyperid = }{35002}
#'             \item{name = }{log range}
#'             \item{short.name = }{range}
#'             \item{initial = }{2}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{loggamma}
#'             \item{param = }{1 0.01}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'dmatern'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'Dense Matern field'}
#'               \item{constr = }{'FALSE'}
#'               \item{nrow.ncol = }{'FALSE'}
#'               \item{augmented = }{'FALSE'}
#'               \item{aug.factor = }{'1'}
#'               \item{aug.constr = }{'NULL'}
#'               \item{n.div.by = }{'NULL'}
#'               \item{n.required = }{'TRUE'}
#'               \item{set.default.values = }{'TRUE'}
#'               \item{status = }{'experimental'}
#'               \item{pdf = }{'dmatern'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 3.
#'          \describe{
#'           \item{Hyperparameter 'theta1'}{
#'              \describe{
#'             \item{hyperid = }{35101}
#'             \item{name = }{log precision}
#'             \item{short.name = }{prec}
#'             \item{initial = }{3}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{pc.prec}
#'             \item{param = }{1 0.01}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta2'}{
#'              \describe{
#'             \item{hyperid = }{35102}
#'             \item{name = }{log range}
#'             \item{short.name = }{range}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{pc.range}
#'             \item{param = }{1 0.5}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta3'}{
#'              \describe{
#'             \item{hyperid = }{35103}
#'             \item{name = }{log nu}
#'             \item{short.name = }{nu}
#'             \item{initial = }{-0.693147180559945}
#'             \item{fixed = }{TRUE}
#'             \item{prior = }{loggamma}
#'             \item{param = }{0.5 1}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'copy'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'Create a copy of a model component'}
#'               \item{constr = }{'FALSE'}
#'               \item{nrow.ncol = }{'FALSE'}
#'               \item{augmented = }{'FALSE'}
#'               \item{aug.factor = }{'1'}
#'               \item{aug.constr = }{'NULL'}
#'               \item{n.div.by = }{'NULL'}
#'               \item{n.required = }{'FALSE'}
#'               \item{set.default.values = }{'FALSE'}
#'               \item{pdf = }{'copy'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 1.
#'          \describe{
#'           \item{Hyperparameter 'theta'}{
#'              \describe{
#'             \item{hyperid = }{36001}
#'             \item{name = }{beta}
#'             \item{short.name = }{b}
#'             \item{initial = }{0}
#'             \item{fixed = }{TRUE}
#'             \item{prior = }{normal}
#'             \item{param = }{1 10}
#'             \item{to.theta = }{`function(x, REPLACE.ME.low, REPLACE.ME.high) {` `                                if (all(is.infinite(c(low, high))) || low == high) {` `                                    return(x)` `                                } else if (all(is.finite(c(low, high)))) {` `                                    stopifnot(low < high)` `                                    return(log(-(low - x) / (high - x)))` `                                } else if (is.finite(low) && is.infinite(high) && high > low) {` `                                    return(log(x - low))` `                                } else {` `                                    stop("Condition not yet implemented")` `                                }` `                            }`}
#'             \item{from.theta = }{`function(x, REPLACE.ME.low, REPLACE.ME.high) {` `                                if (all(is.infinite(c(low, high))) || low == high) {` `                                    return(x)` `                                } else if (all(is.finite(c(low, high)))) {` `                                    stopifnot(low < high)` `                                    return(low + exp(x) / (1 + exp(x)) * (high - low))` `                                } else if (is.finite(low) && is.infinite(high) && high > low) {` `                                    return(low + exp(x))` `                                } else {` `                                    stop("Condition not yet implemented")` `                                }` `                            }`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'scopy'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'Create a scopy of a model component'}
#'               \item{constr = }{'FALSE'}
#'               \item{nrow.ncol = }{'FALSE'}
#'               \item{augmented = }{'FALSE'}
#'               \item{aug.factor = }{'1'}
#'               \item{aug.constr = }{'NULL'}
#'               \item{n.div.by = }{'NULL'}
#'               \item{n.required = }{'FALSE'}
#'               \item{set.default.values = }{'FALSE'}
#'               \item{status = }{'experimental'}
#'               \item{pdf = }{'scopy'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 15.
#'          \describe{
#'           \item{Hyperparameter 'theta1'}{
#'              \describe{
#'             \item{hyperid = }{36101}
#'             \item{name = }{beta1}
#'             \item{short.name = }{b1}
#'             \item{initial = }{0.1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta2'}{
#'              \describe{
#'             \item{hyperid = }{36102}
#'             \item{name = }{beta2}
#'             \item{short.name = }{b2}
#'             \item{initial = }{0.1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta3'}{
#'              \describe{
#'             \item{hyperid = }{36103}
#'             \item{name = }{beta3}
#'             \item{short.name = }{b3}
#'             \item{initial = }{0.1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta4'}{
#'              \describe{
#'             \item{hyperid = }{36104}
#'             \item{name = }{beta4}
#'             \item{short.name = }{b4}
#'             \item{initial = }{0.1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta5'}{
#'              \describe{
#'             \item{hyperid = }{36105}
#'             \item{name = }{beta5}
#'             \item{short.name = }{b5}
#'             \item{initial = }{0.1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta6'}{
#'              \describe{
#'             \item{hyperid = }{36106}
#'             \item{name = }{beta6}
#'             \item{short.name = }{b6}
#'             \item{initial = }{0.1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta7'}{
#'              \describe{
#'             \item{hyperid = }{36107}
#'             \item{name = }{beta7}
#'             \item{short.name = }{b7}
#'             \item{initial = }{0.1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta8'}{
#'              \describe{
#'             \item{hyperid = }{36108}
#'             \item{name = }{beta8}
#'             \item{short.name = }{b8}
#'             \item{initial = }{0.1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta9'}{
#'              \describe{
#'             \item{hyperid = }{36109}
#'             \item{name = }{beta9}
#'             \item{short.name = }{b9}
#'             \item{initial = }{0.1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta10'}{
#'              \describe{
#'             \item{hyperid = }{36110}
#'             \item{name = }{beta10}
#'             \item{short.name = }{b10}
#'             \item{initial = }{0.1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta11'}{
#'              \describe{
#'             \item{hyperid = }{36111}
#'             \item{name = }{beta11}
#'             \item{short.name = }{b11}
#'             \item{initial = }{0.1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta12'}{
#'              \describe{
#'             \item{hyperid = }{36112}
#'             \item{name = }{beta12}
#'             \item{short.name = }{b12}
#'             \item{initial = }{0.1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta13'}{
#'              \describe{
#'             \item{hyperid = }{36113}
#'             \item{name = }{beta13}
#'             \item{short.name = }{b13}
#'             \item{initial = }{0.1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta14'}{
#'              \describe{
#'             \item{hyperid = }{36114}
#'             \item{name = }{beta14}
#'             \item{short.name = }{b14}
#'             \item{initial = }{0.1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta15'}{
#'              \describe{
#'             \item{hyperid = }{36115}
#'             \item{name = }{beta15}
#'             \item{short.name = }{b15}
#'             \item{initial = }{0.1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'clinear'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'Constrained linear effect'}
#'               \item{constr = }{'FALSE'}
#'               \item{nrow.ncol = }{'FALSE'}
#'               \item{augmented = }{'FALSE'}
#'               \item{aug.factor = }{'1'}
#'               \item{aug.constr = }{'NULL'}
#'               \item{n.div.by = }{'NULL'}
#'               \item{n.required = }{'FALSE'}
#'               \item{set.default.values = }{'FALSE'}
#'               \item{pdf = }{'clinear'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 1.
#'          \describe{
#'           \item{Hyperparameter 'theta'}{
#'              \describe{
#'             \item{hyperid = }{37001}
#'             \item{name = }{beta}
#'             \item{short.name = }{b}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{1 10}
#'             \item{to.theta = }{`function(x, REPLACE.ME.low, REPLACE.ME.high) {` `                                if (all(is.infinite(c(low, high))) || low == high) {` `                                    stopifnot(low < high)` `                                    return(x)` `                                } else if (all(is.finite(c(low, high)))) {` `                                    stopifnot(low < high)` `                                    return(log(-(low - x) / (high - x)))` `                                } else if (is.finite(low) && is.infinite(high) && high > low) {` `                                    return(log(x - low))` `                                } else {` `                                    stop("Condition not yet implemented")` `                                }` `                            }`}
#'             \item{from.theta = }{`function(x, REPLACE.ME.low, REPLACE.ME.high) {` `                                if (all(is.infinite(c(low, high))) || low == high) {` `                                    stopifnot(low < high)` `                                    return(x)` `                                } else if (all(is.finite(c(low, high)))) {` `                                    stopifnot(low < high)` `                                    return(low + exp(x) / (1 + exp(x)) * (high - low))` `                                } else if (is.finite(low) && is.infinite(high) && high > low) {` `                                    return(low + exp(x))` `                                } else {` `                                    stop("Condition not yet implemented")` `                                }` `                            }`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'sigm'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'Sigmoidal effect of a covariate'}
#'               \item{constr = }{'FALSE'}
#'               \item{nrow.ncol = }{'FALSE'}
#'               \item{augmented = }{'FALSE'}
#'               \item{aug.factor = }{'1'}
#'               \item{aug.constr = }{'NULL'}
#'               \item{n.div.by = }{'NULL'}
#'               \item{n.required = }{'FALSE'}
#'               \item{set.default.values = }{'FALSE'}
#'               \item{status = }{'experimental'}
#'               \item{pdf = }{'sigm'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 3.
#'          \describe{
#'           \item{Hyperparameter 'theta1'}{
#'              \describe{
#'             \item{hyperid = }{38001}
#'             \item{name = }{beta}
#'             \item{short.name = }{b}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{1 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta2'}{
#'              \describe{
#'             \item{hyperid = }{38002}
#'             \item{name = }{loghalflife}
#'             \item{short.name = }{halflife}
#'             \item{initial = }{3}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{loggamma}
#'             \item{param = }{3 1}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta3'}{
#'              \describe{
#'             \item{hyperid = }{38003}
#'             \item{name = }{logshape}
#'             \item{short.name = }{shape}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{loggamma}
#'             \item{param = }{10 10}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'revsigm'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'Reverse sigmoidal effect of a covariate'}
#'               \item{constr = }{'FALSE'}
#'               \item{nrow.ncol = }{'FALSE'}
#'               \item{augmented = }{'FALSE'}
#'               \item{aug.factor = }{'1'}
#'               \item{aug.constr = }{'NULL'}
#'               \item{n.div.by = }{'NULL'}
#'               \item{n.required = }{'FALSE'}
#'               \item{set.default.values = }{'FALSE'}
#'               \item{status = }{'experimental'}
#'               \item{pdf = }{'sigm'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 3.
#'          \describe{
#'           \item{Hyperparameter 'theta1'}{
#'              \describe{
#'             \item{hyperid = }{39001}
#'             \item{name = }{beta}
#'             \item{short.name = }{b}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{1 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta2'}{
#'              \describe{
#'             \item{hyperid = }{39002}
#'             \item{name = }{loghalflife}
#'             \item{short.name = }{halflife}
#'             \item{initial = }{3}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{loggamma}
#'             \item{param = }{3 1}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta3'}{
#'              \describe{
#'             \item{hyperid = }{39003}
#'             \item{name = }{logshape}
#'             \item{short.name = }{shape}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{loggamma}
#'             \item{param = }{10 10}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'log1exp'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'A nonlinear model of a covariate'}
#'               \item{constr = }{'FALSE'}
#'               \item{nrow.ncol = }{'FALSE'}
#'               \item{augmented = }{'FALSE'}
#'               \item{aug.factor = }{'1'}
#'               \item{aug.constr = }{'NULL'}
#'               \item{n.div.by = }{'NULL'}
#'               \item{n.required = }{'FALSE'}
#'               \item{set.default.values = }{'FALSE'}
#'               \item{status = }{'experimental'}
#'               \item{pdf = }{'log1exp'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 3.
#'          \describe{
#'           \item{Hyperparameter 'theta1'}{
#'              \describe{
#'             \item{hyperid = }{39011}
#'             \item{name = }{beta}
#'             \item{short.name = }{b}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 1}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta2'}{
#'              \describe{
#'             \item{hyperid = }{39012}
#'             \item{name = }{alpha}
#'             \item{short.name = }{a}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 1}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta3'}{
#'              \describe{
#'             \item{hyperid = }{39013}
#'             \item{name = }{gamma}
#'             \item{short.name = }{g}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 1}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'logdist'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'A nonlinear model of a covariate'}
#'               \item{constr = }{'FALSE'}
#'               \item{nrow.ncol = }{'FALSE'}
#'               \item{augmented = }{'FALSE'}
#'               \item{aug.factor = }{'1'}
#'               \item{aug.constr = }{'NULL'}
#'               \item{n.div.by = }{'NULL'}
#'               \item{n.required = }{'FALSE'}
#'               \item{set.default.values = }{'FALSE'}
#'               \item{status = }{'experimental'}
#'               \item{pdf = }{'logdist'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 3.
#'          \describe{
#'           \item{Hyperparameter 'theta1'}{
#'              \describe{
#'             \item{hyperid = }{39021}
#'             \item{name = }{beta}
#'             \item{short.name = }{b}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 1}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta2'}{
#'              \describe{
#'             \item{hyperid = }{39022}
#'             \item{name = }{alpha1}
#'             \item{short.name = }{a1}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{loggamma}
#'             \item{param = }{0.1 1}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta3'}{
#'              \describe{
#'             \item{hyperid = }{39023}
#'             \item{name = }{alpha2}
#'             \item{short.name = }{a2}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{loggamma}
#'             \item{param = }{0.1 1}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'          }
#'        }
#'      }
#' @section 'group':
#'
#'    Valid models in this section are:
#'      \describe{
#'       \item{Model 'exchangeable'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'Exchangeable correlations'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 1.
#'          \describe{
#'           \item{Hyperparameter 'theta'}{
#'              \describe{
#'             \item{hyperid = }{40001}
#'             \item{name = }{logit correlation}
#'             \item{short.name = }{rho}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 0.2}
#'             \item{to.theta = }{`function(x, REPLACE.ME.ngroup) log((1 + x * (ngroup - 1)) / (1 - x))`}
#'             \item{from.theta = }{`function(x, REPLACE.ME.ngroup) (exp(x) - 1) / (exp(x) + ngroup - 1)`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'exchangeablepos'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'Exchangeable positive correlations'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 1.
#'          \describe{
#'           \item{Hyperparameter 'theta'}{
#'              \describe{
#'             \item{hyperid = }{40101}
#'             \item{name = }{logit correlation}
#'             \item{short.name = }{rho}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{pc.cor0}
#'             \item{param = }{0.5 0.5}
#'             \item{to.theta = }{`function(x) log(x / (1 - x))`}
#'             \item{from.theta = }{`function(x) exp(x) / (1 + exp(x))`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'ar1'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'AR(1) correlations'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 1.
#'          \describe{
#'           \item{Hyperparameter 'theta'}{
#'              \describe{
#'             \item{hyperid = }{41001}
#'             \item{name = }{logit correlation}
#'             \item{short.name = }{rho}
#'             \item{initial = }{2}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 0.15}
#'             \item{to.theta = }{`function(x) log((1 + x) / (1 - x))`}
#'             \item{from.theta = }{`function(x) 2 * exp(x) / (1 + exp(x)) - 1`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'ar'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'AR(p) correlations'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 11.
#'          \describe{
#'           \item{Hyperparameter 'theta1'}{
#'              \describe{
#'             \item{hyperid = }{42001}
#'             \item{name = }{log precision}
#'             \item{short.name = }{prec}
#'             \item{initial = }{0}
#'             \item{fixed = }{TRUE}
#'             \item{prior = }{pc.prec}
#'             \item{param = }{3 0.01}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta2'}{
#'              \describe{
#'             \item{hyperid = }{42002}
#'             \item{name = }{pacf1}
#'             \item{short.name = }{pacf1}
#'             \item{initial = }{2}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{pc.cor0}
#'             \item{param = }{0.5 0.5}
#'             \item{to.theta = }{`function(x) log((1 + x) / (1 - x))`}
#'             \item{from.theta = }{`function(x) 2 * exp(x) / (1 + exp(x)) - 1`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta3'}{
#'              \describe{
#'             \item{hyperid = }{42003}
#'             \item{name = }{pacf2}
#'             \item{short.name = }{pacf2}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{pc.cor0}
#'             \item{param = }{0.5 0.4}
#'             \item{to.theta = }{`function(x) log((1 + x) / (1 - x))`}
#'             \item{from.theta = }{`function(x) 2 * exp(x) / (1 + exp(x)) - 1`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta4'}{
#'              \describe{
#'             \item{hyperid = }{42004}
#'             \item{name = }{pacf3}
#'             \item{short.name = }{pacf3}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{pc.cor0}
#'             \item{param = }{0.5 0.3}
#'             \item{to.theta = }{`function(x) log((1 + x) / (1 - x))`}
#'             \item{from.theta = }{`function(x) 2 * exp(x) / (1 + exp(x)) - 1`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta5'}{
#'              \describe{
#'             \item{hyperid = }{42005}
#'             \item{name = }{pacf4}
#'             \item{short.name = }{pacf4}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{pc.cor0}
#'             \item{param = }{0.5 0.2}
#'             \item{to.theta = }{`function(x) log((1 + x) / (1 - x))`}
#'             \item{from.theta = }{`function(x) 2 * exp(x) / (1 + exp(x)) - 1`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta6'}{
#'              \describe{
#'             \item{hyperid = }{42006}
#'             \item{name = }{pacf5}
#'             \item{short.name = }{pacf5}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{pc.cor0}
#'             \item{param = }{0.5 0.1}
#'             \item{to.theta = }{`function(x) log((1 + x) / (1 - x))`}
#'             \item{from.theta = }{`function(x) 2 * exp(x) / (1 + exp(x)) - 1`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta7'}{
#'              \describe{
#'             \item{hyperid = }{42007}
#'             \item{name = }{pacf6}
#'             \item{short.name = }{pacf6}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{pc.cor0}
#'             \item{param = }{0.5 0.1}
#'             \item{to.theta = }{`function(x) log((1 + x) / (1 - x))`}
#'             \item{from.theta = }{`function(x) 2 * exp(x) / (1 + exp(x)) - 1`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta8'}{
#'              \describe{
#'             \item{hyperid = }{42008}
#'             \item{name = }{pacf7}
#'             \item{short.name = }{pacf7}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{pc.cor0}
#'             \item{param = }{0.5 0.1}
#'             \item{to.theta = }{`function(x) log((1 + x) / (1 - x))`}
#'             \item{from.theta = }{`function(x) 2 * exp(x) / (1 + exp(x)) - 1`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta9'}{
#'              \describe{
#'             \item{hyperid = }{42009}
#'             \item{name = }{pacf8}
#'             \item{short.name = }{pacf8}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{pc.cor0}
#'             \item{param = }{0.5 0.1}
#'             \item{to.theta = }{`function(x) log((1 + x) / (1 - x))`}
#'             \item{from.theta = }{`function(x) 2 * exp(x) / (1 + exp(x)) - 1`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta10'}{
#'              \describe{
#'             \item{hyperid = }{42010}
#'             \item{name = }{pacf9}
#'             \item{short.name = }{pacf9}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{pc.cor0}
#'             \item{param = }{0.5 0.1}
#'             \item{to.theta = }{`function(x) log((1 + x) / (1 - x))`}
#'             \item{from.theta = }{`function(x) 2 * exp(x) / (1 + exp(x)) - 1`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta11'}{
#'              \describe{
#'             \item{hyperid = }{42011}
#'             \item{name = }{pacf10}
#'             \item{short.name = }{pacf10}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{pc.cor0}
#'             \item{param = }{0.5 0.1}
#'             \item{to.theta = }{`function(x) log((1 + x) / (1 - x))`}
#'             \item{from.theta = }{`function(x) 2 * exp(x) / (1 + exp(x)) - 1`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'rw1'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'Random walk of order 1'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 1.
#'          \describe{
#'           \item{Hyperparameter 'theta'}{
#'              \describe{
#'             \item{hyperid = }{43001}
#'             \item{name = }{log precision}
#'             \item{short.name = }{prec}
#'             \item{prior = }{loggamma}
#'             \item{param = }{1 5e-05}
#'             \item{initial = }{0}
#'             \item{fixed = }{TRUE}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'rw2'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'Random walk of order 2'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 1.
#'          \describe{
#'           \item{Hyperparameter 'theta'}{
#'              \describe{
#'             \item{hyperid = }{44001}
#'             \item{name = }{log precision}
#'             \item{short.name = }{prec}
#'             \item{prior = }{loggamma}
#'             \item{param = }{1 5e-05}
#'             \item{initial = }{0}
#'             \item{fixed = }{TRUE}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'besag'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'Besag model'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 1.
#'          \describe{
#'           \item{Hyperparameter 'theta'}{
#'              \describe{
#'             \item{hyperid = }{45001}
#'             \item{name = }{log precision}
#'             \item{short.name = }{prec}
#'             \item{prior = }{loggamma}
#'             \item{param = }{1 5e-05}
#'             \item{initial = }{0}
#'             \item{fixed = }{TRUE}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'iid'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'Independent model'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 1.
#'          \describe{
#'           \item{Hyperparameter 'theta'}{
#'              \describe{
#'             \item{hyperid = }{46001}
#'             \item{name = }{log precision}
#'             \item{short.name = }{prec}
#'             \item{prior = }{loggamma}
#'             \item{param = }{1 5e-05}
#'             \item{initial = }{0}
#'             \item{fixed = }{TRUE}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'          }
#'        }
#'      }
#' @section 'scopy':
#'
#'    Valid models in this section are:
#'      \describe{
#'       \item{Model 'rw1'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'Random walk of order 1'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 0.
#'        }
#'       \item{Model 'rw2'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'Random walk of order 2'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 0.
#'        }
#'      }
#' @section 'mix':
#'
#'    Valid models in this section are:
#'      \describe{
#'       \item{Model 'gaussian'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'Gaussian mixture'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 1.
#'          \describe{
#'           \item{Hyperparameter 'theta'}{
#'              \describe{
#'             \item{hyperid = }{47001}
#'             \item{name = }{log precision}
#'             \item{short.name = }{prec}
#'             \item{prior = }{pc.prec}
#'             \item{param = }{1 0.01}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'loggamma'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'LogGamma mixture'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 1.
#'          \describe{
#'           \item{Hyperparameter 'theta'}{
#'              \describe{
#'             \item{hyperid = }{47101}
#'             \item{name = }{log precision}
#'             \item{short.name = }{prec}
#'             \item{prior = }{pc.mgamma}
#'             \item{param = }{4.8}
#'             \item{initial = }{4}
#'             \item{fixed = }{FALSE}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'mloggamma'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'Minus-LogGamma mixture'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 1.
#'          \describe{
#'           \item{Hyperparameter 'theta'}{
#'              \describe{
#'             \item{hyperid = }{47201}
#'             \item{name = }{log precision}
#'             \item{short.name = }{prec}
#'             \item{prior = }{pc.mgamma}
#'             \item{param = }{4.8}
#'             \item{initial = }{4}
#'             \item{fixed = }{FALSE}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'          }
#'        }
#'      }
#' @section 'link':
#'
#'    Valid models in this section are:
#'      \describe{
#'       \item{Model 'default'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'The default link'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 0.
#'        }
#'       \item{Model 'cloglog'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'The complementary log-log link'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 0.
#'        }
#'       \item{Model 'ccloglog'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'The complement complementary log-log link'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 0.
#'        }
#'       \item{Model 'loglog'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'The log-log link'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 0.
#'        }
#'       \item{Model 'identity'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'The identity link'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 0.
#'        }
#'       \item{Model 'inverse'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'The inverse link'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 0.
#'        }
#'       \item{Model 'log'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'The log-link'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 0.
#'        }
#'       \item{Model 'loga'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'The loga-link'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 0.
#'        }
#'       \item{Model 'neglog'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'The negative log-link'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 0.
#'        }
#'       \item{Model 'logit'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'The logit-link'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 0.
#'        }
#'       \item{Model 'probit'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'The probit-link'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 0.
#'        }
#'       \item{Model 'cauchit'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'The cauchit-link'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 0.
#'        }
#'       \item{Model 'tan'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'The tan-link'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 0.
#'        }
#'       \item{Model 'quantile'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'The quantile-link'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 0.
#'        }
#'       \item{Model 'pquantile'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'The population quantile-link'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 0.
#'        }
#'       \item{Model 'sslogit'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'Logit link with sensitivity and specificity'}
#'               \item{status = }{'disabled'}
#'               \item{pdf = }{'NA'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 2.
#'          \describe{
#'           \item{Hyperparameter 'theta1'}{
#'              \describe{
#'             \item{hyperid = }{48001}
#'             \item{name = }{sensitivity}
#'             \item{short.name = }{sens}
#'             \item{prior = }{logitbeta}
#'             \item{param = }{10 5}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{to.theta = }{`function(x) log(x / (1 - x))`}
#'             \item{from.theta = }{`function(x) exp(x) / (1 + exp(x))`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta2'}{
#'              \describe{
#'             \item{hyperid = }{48002}
#'             \item{name = }{specificity}
#'             \item{short.name = }{spec}
#'             \item{prior = }{logitbeta}
#'             \item{param = }{10 5}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{to.theta = }{`function(x) log(x / (1 - x))`}
#'             \item{from.theta = }{`function(x) exp(x) / (1 + exp(x))`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'logoffset'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'Log-link with an offset'}
#'               \item{pdf = }{'logoffset'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 1.
#'          \describe{
#'           \item{Hyperparameter 'theta'}{
#'              \describe{
#'             \item{hyperid = }{49001}
#'             \item{name = }{beta}
#'             \item{short.name = }{b}
#'             \item{prior = }{normal}
#'             \item{param = }{0 100}
#'             \item{initial = }{0}
#'             \item{fixed = }{TRUE}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'logitoffset'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'Logit-link with an offset'}
#'               \item{status = }{'experimental'}
#'               \item{pdf = }{'logitoffset'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 1.
#'          \describe{
#'           \item{Hyperparameter 'theta'}{
#'              \describe{
#'             \item{hyperid = }{49011}
#'             \item{name = }{prob}
#'             \item{short.name = }{p}
#'             \item{prior = }{normal}
#'             \item{param = }{-1 100}
#'             \item{initial = }{-1}
#'             \item{fixed = }{FALSE}
#'             \item{to.theta = }{`function(x) log(x / (1 - x))`}
#'             \item{from.theta = }{`function(x) exp(x) / (1 + exp(x))`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'robit'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'Robit link'}
#'               \item{status = }{'experimental'}
#'               \item{pdf = }{'robit'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 1.
#'          \describe{
#'           \item{Hyperparameter 'theta'}{
#'              \describe{
#'             \item{hyperid = }{49021}
#'             \item{name = }{log degrees of freedom}
#'             \item{short.name = }{dof}
#'             \item{initial = }{1.6094379124341}
#'             \item{fixed = }{TRUE}
#'             \item{prior = }{pc.dof}
#'             \item{param = }{50 0.5}
#'             \item{to.theta = }{`function(x) log(x - 2)`}
#'             \item{from.theta = }{`function(x) 2 + exp(x)`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'sn'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'Skew-normal link'}
#'               \item{pdf = }{'linksn'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 2.
#'          \describe{
#'           \item{Hyperparameter 'theta1'}{
#'              \describe{
#'             \item{hyperid = }{49031}
#'             \item{name = }{skewness}
#'             \item{short.name = }{skew}
#'             \item{initial = }{0.00123456789}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{pc.sn}
#'             \item{param = }{10}
#'             \item{to.theta = }{`function(x, skew.max = 0.988) log((1 + x / skew.max) / (1 - x / skew.max))`}
#'             \item{from.theta = }{`function(x, skew.max = 0.988) skew.max * (2 * exp(x) / (1 + exp(x)) - 1)`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta2'}{
#'              \describe{
#'             \item{hyperid = }{49032}
#'             \item{name = }{intercept}
#'             \item{short.name = }{intercept}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{linksnintercept}
#'             \item{param = }{0 0}
#'             \item{to.theta = }{`function(x) log(x / (1 - x))`}
#'             \item{from.theta = }{`function(x) exp(x) / (1 + exp(x))`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'powerlogit'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'Power logit link'}
#'               \item{pdf = }{'powerlogit'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 2.
#'          \describe{
#'           \item{Hyperparameter 'theta1'}{
#'              \describe{
#'             \item{hyperid = }{49131}
#'             \item{name = }{power}
#'             \item{short.name = }{power}
#'             \item{initial = }{0.00123456789}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 10}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta2'}{
#'              \describe{
#'             \item{hyperid = }{49132}
#'             \item{name = }{intercept}
#'             \item{short.name = }{intercept}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{logitbeta}
#'             \item{param = }{1 1}
#'             \item{to.theta = }{`function(x) log(x / (1 - x))`}
#'             \item{from.theta = }{`function(x) exp(x) / (1 + exp(x))`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'test1'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'A test1-link function (experimental)'}
#'               \item{pdf = }{'NA'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 1.
#'          \describe{
#'           \item{Hyperparameter 'theta'}{
#'              \describe{
#'             \item{hyperid = }{50001}
#'             \item{name = }{beta}
#'             \item{short.name = }{b}
#'             \item{prior = }{normal}
#'             \item{param = }{0 100}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'special1'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'A special1-link function (experimental)'}
#'               \item{pdf = }{'NA'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 11.
#'          \describe{
#'           \item{Hyperparameter 'theta1'}{
#'              \describe{
#'             \item{hyperid = }{51001}
#'             \item{name = }{log precision}
#'             \item{short.name = }{prec}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{loggamma}
#'             \item{param = }{1 1}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta2'}{
#'              \describe{
#'             \item{hyperid = }{51002}
#'             \item{name = }{beta1}
#'             \item{short.name = }{beta1}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{mvnorm}
#'             \item{param = }{0 100}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta3'}{
#'              \describe{
#'             \item{hyperid = }{51003}
#'             \item{name = }{beta2}
#'             \item{short.name = }{beta2}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta4'}{
#'              \describe{
#'             \item{hyperid = }{51004}
#'             \item{name = }{beta3}
#'             \item{short.name = }{beta3}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta5'}{
#'              \describe{
#'             \item{hyperid = }{51005}
#'             \item{name = }{beta4}
#'             \item{short.name = }{beta4}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta6'}{
#'              \describe{
#'             \item{hyperid = }{51006}
#'             \item{name = }{beta5}
#'             \item{short.name = }{beta5}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta7'}{
#'              \describe{
#'             \item{hyperid = }{51007}
#'             \item{name = }{beta6}
#'             \item{short.name = }{beta6}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta8'}{
#'              \describe{
#'             \item{hyperid = }{51008}
#'             \item{name = }{beta7}
#'             \item{short.name = }{beta7}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta9'}{
#'              \describe{
#'             \item{hyperid = }{51009}
#'             \item{name = }{beta8}
#'             \item{short.name = }{beta8}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta10'}{
#'              \describe{
#'             \item{hyperid = }{51010}
#'             \item{name = }{beta9}
#'             \item{short.name = }{beta9}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta11'}{
#'              \describe{
#'             \item{hyperid = }{51011}
#'             \item{name = }{beta10}
#'             \item{short.name = }{beta10}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'special2'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'A special2-link function (experimental)'}
#'               \item{pdf = }{'NA'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 1.
#'          \describe{
#'           \item{Hyperparameter 'theta'}{
#'              \describe{
#'             \item{hyperid = }{52001}
#'             \item{name = }{beta}
#'             \item{short.name = }{b}
#'             \item{prior = }{normal}
#'             \item{param = }{0 10}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'          }
#'        }
#'      }
#' @section 'predictor':
#'
#'    Valid models in this section are:
#'      \describe{
#'       \item{Model 'predictor'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'(do not use)'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 1.
#'          \describe{
#'           \item{Hyperparameter 'theta'}{
#'              \describe{
#'             \item{hyperid = }{53001}
#'             \item{name = }{log precision}
#'             \item{short.name = }{prec}
#'             \item{initial = }{13.8155105579643}
#'             \item{fixed = }{TRUE}
#'             \item{prior = }{loggamma}
#'             \item{param = }{1 1e-05}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'          }
#'        }
#'      }
#' @section 'hazard':
#'
#'    Valid models in this section are:
#'      \describe{
#'       \item{Model 'rw1'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'A random walk of order 1 for the log-hazard'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 1.
#'          \describe{
#'           \item{Hyperparameter 'theta'}{
#'              \describe{
#'             \item{hyperid = }{54001}
#'             \item{name = }{log precision}
#'             \item{short.name = }{prec}
#'             \item{initial = }{4}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{loggamma}
#'             \item{param = }{1 5e-05}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'rw2'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'A random walk of order 2 for the log-hazard'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 1.
#'          \describe{
#'           \item{Hyperparameter 'theta'}{
#'              \describe{
#'             \item{hyperid = }{55001}
#'             \item{name = }{log precision}
#'             \item{short.name = }{prec}
#'             \item{initial = }{4}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{loggamma}
#'             \item{param = }{1 5e-05}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'iid'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'An iid model for the log-hazard'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 1.
#'          \describe{
#'           \item{Hyperparameter 'theta'}{
#'              \describe{
#'             \item{hyperid = }{55501}
#'             \item{name = }{log precision}
#'             \item{short.name = }{prec}
#'             \item{initial = }{4}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{loggamma}
#'             \item{param = }{1 5e-05}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'          }
#'        }
#'      }
#' @section 'likelihood':
#'
#'    Valid models in this section are:
#'      \describe{
#'       \item{Model 'poisson'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'The Poisson likelihood'}
#'               \item{survival = }{'FALSE'}
#'               \item{discrete = }{'TRUE'}
#'               \item{link = }{'default log logoffset quantile test1 special1 special2'}
#'               \item{pdf = }{'poisson'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 0.
#'        }
#'       \item{Model 'xpoisson'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'The Poisson likelihood (expert version)'}
#'               \item{survival = }{'FALSE'}
#'               \item{discrete = }{'TRUE'}
#'               \item{link = }{'default log logoffset quantile test1 special1 special2'}
#'               \item{pdf = }{'poisson'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 0.
#'        }
#'       \item{Model 'cenpoisson'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'Then censored Poisson likelihood'}
#'               \item{survival = }{'FALSE'}
#'               \item{discrete = }{'TRUE'}
#'               \item{link = }{'default log logoffset test1 special1 special2'}
#'               \item{pdf = }{'cenpoisson'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 0.
#'        }
#'       \item{Model 'cenpoisson2'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'Then censored Poisson likelihood (version 2)'}
#'               \item{survival = }{'FALSE'}
#'               \item{discrete = }{'TRUE'}
#'               \item{link = }{'default log logoffset test1 special1 special2'}
#'               \item{pdf = }{'cenpoisson2'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 0.
#'        }
#'       \item{Model 'gpoisson'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'The generalized Poisson likelihood'}
#'               \item{survival = }{'FALSE'}
#'               \item{discrete = }{'TRUE'}
#'               \item{link = }{'default log logoffset'}
#'               \item{pdf = }{'gpoisson'}
#'               \item{status = }{'experimental'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 2.
#'          \describe{
#'           \item{Hyperparameter 'theta1'}{
#'              \describe{
#'             \item{hyperid = }{56001}
#'             \item{name = }{overdispersion}
#'             \item{short.name = }{phi}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{loggamma}
#'             \item{param = }{1 1}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta2'}{
#'              \describe{
#'             \item{hyperid = }{56002}
#'             \item{name = }{p}
#'             \item{short.name = }{p}
#'             \item{initial = }{1}
#'             \item{fixed = }{TRUE}
#'             \item{prior = }{normal}
#'             \item{param = }{1 100}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'poisson.special1'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'The Poisson.special1 likelihood'}
#'               \item{survival = }{'FALSE'}
#'               \item{discrete = }{'TRUE'}
#'               \item{link = }{'default log'}
#'               \item{pdf = }{'poisson-special'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 1.
#'          \describe{
#'           \item{Hyperparameter 'theta'}{
#'              \describe{
#'             \item{hyperid = }{56100}
#'             \item{name = }{logit probability}
#'             \item{short.name = }{prob}
#'             \item{initial = }{-1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{gaussian}
#'             \item{param = }{-1 0.2}
#'             \item{to.theta = }{`function(x) log(x / (1 - x))`}
#'             \item{from.theta = }{`function(x) exp(x) / (1 + exp(x))`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model '0poisson'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'New 0-inflated Poisson'}
#'               \item{status = }{'experimental'}
#'               \item{survival = }{'FALSE'}
#'               \item{discrete = }{'TRUE'}
#'               \item{link = }{'default log quantile'}
#'               \item{link.simple = }{'default logit cauchit probit cloglog ccloglog'}
#'               \item{pdf = }{'0inflated'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 10.
#'          \describe{
#'           \item{Hyperparameter 'theta1'}{
#'              \describe{
#'             \item{hyperid = }{56201}
#'             \item{name = }{beta1}
#'             \item{short.name = }{beta1}
#'             \item{initial = }{-4}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{-4 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta2'}{
#'              \describe{
#'             \item{hyperid = }{56202}
#'             \item{name = }{beta2}
#'             \item{short.name = }{beta2}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta3'}{
#'              \describe{
#'             \item{hyperid = }{56203}
#'             \item{name = }{beta3}
#'             \item{short.name = }{beta3}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta4'}{
#'              \describe{
#'             \item{hyperid = }{56204}
#'             \item{name = }{beta4}
#'             \item{short.name = }{beta4}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta5'}{
#'              \describe{
#'             \item{hyperid = }{56205}
#'             \item{name = }{beta5}
#'             \item{short.name = }{beta5}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta6'}{
#'              \describe{
#'             \item{hyperid = }{56206}
#'             \item{name = }{beta6}
#'             \item{short.name = }{beta6}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta7'}{
#'              \describe{
#'             \item{hyperid = }{56207}
#'             \item{name = }{beta7}
#'             \item{short.name = }{beta7}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta8'}{
#'              \describe{
#'             \item{hyperid = }{56208}
#'             \item{name = }{beta8}
#'             \item{short.name = }{beta8}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta9'}{
#'              \describe{
#'             \item{hyperid = }{56209}
#'             \item{name = }{beta9}
#'             \item{short.name = }{beta9}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta10'}{
#'              \describe{
#'             \item{hyperid = }{56210}
#'             \item{name = }{beta10}
#'             \item{short.name = }{beta10}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model '0poissonS'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'New 0-inflated Poisson Swap'}
#'               \item{status = }{'experimental'}
#'               \item{survival = }{'FALSE'}
#'               \item{discrete = }{'TRUE'}
#'               \item{link = }{'default logit loga cauchit probit cloglog ccloglog loglog log sslogit logitoffset quantile pquantile robit sn powerlogit'}
#'               \item{link.simple = }{'default log'}
#'               \item{pdf = }{'0inflated'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 10.
#'          \describe{
#'           \item{Hyperparameter 'theta1'}{
#'              \describe{
#'             \item{hyperid = }{56301}
#'             \item{name = }{beta1}
#'             \item{short.name = }{beta1}
#'             \item{initial = }{-4}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{-4 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta2'}{
#'              \describe{
#'             \item{hyperid = }{56302}
#'             \item{name = }{beta2}
#'             \item{short.name = }{beta2}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta3'}{
#'              \describe{
#'             \item{hyperid = }{56303}
#'             \item{name = }{beta3}
#'             \item{short.name = }{beta3}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta4'}{
#'              \describe{
#'             \item{hyperid = }{56304}
#'             \item{name = }{beta4}
#'             \item{short.name = }{beta4}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta5'}{
#'              \describe{
#'             \item{hyperid = }{56305}
#'             \item{name = }{beta5}
#'             \item{short.name = }{beta5}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta6'}{
#'              \describe{
#'             \item{hyperid = }{56306}
#'             \item{name = }{beta6}
#'             \item{short.name = }{beta6}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta7'}{
#'              \describe{
#'             \item{hyperid = }{56307}
#'             \item{name = }{beta7}
#'             \item{short.name = }{beta7}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta8'}{
#'              \describe{
#'             \item{hyperid = }{56308}
#'             \item{name = }{beta8}
#'             \item{short.name = }{beta8}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta9'}{
#'              \describe{
#'             \item{hyperid = }{56309}
#'             \item{name = }{beta9}
#'             \item{short.name = }{beta9}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta10'}{
#'              \describe{
#'             \item{hyperid = }{56310}
#'             \item{name = }{beta10}
#'             \item{short.name = }{beta10}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'bell'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'The Bell likelihood'}
#'               \item{status = }{'experimental'}
#'               \item{survival = }{'FALSE'}
#'               \item{discrete = }{'TRUE'}
#'               \item{link = }{'default log'}
#'               \item{pdf = }{'bell'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 0.
#'        }
#'       \item{Model '0binomial'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'New 0-inflated Binomial'}
#'               \item{status = }{'experimental'}
#'               \item{survival = }{'FALSE'}
#'               \item{discrete = }{'TRUE'}
#'               \item{link = }{'default logit loga cauchit probit cloglog ccloglog loglog log'}
#'               \item{link.simple = }{'default logit cauchit probit cloglog ccloglog'}
#'               \item{pdf = }{'0inflated'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 10.
#'          \describe{
#'           \item{Hyperparameter 'theta1'}{
#'              \describe{
#'             \item{hyperid = }{56401}
#'             \item{name = }{beta1}
#'             \item{short.name = }{beta1}
#'             \item{initial = }{-4}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{-4 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta2'}{
#'              \describe{
#'             \item{hyperid = }{56402}
#'             \item{name = }{beta2}
#'             \item{short.name = }{beta2}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta3'}{
#'              \describe{
#'             \item{hyperid = }{56403}
#'             \item{name = }{beta3}
#'             \item{short.name = }{beta3}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta4'}{
#'              \describe{
#'             \item{hyperid = }{56404}
#'             \item{name = }{beta4}
#'             \item{short.name = }{beta4}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta5'}{
#'              \describe{
#'             \item{hyperid = }{56405}
#'             \item{name = }{beta5}
#'             \item{short.name = }{beta5}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta6'}{
#'              \describe{
#'             \item{hyperid = }{56406}
#'             \item{name = }{beta6}
#'             \item{short.name = }{beta6}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta7'}{
#'              \describe{
#'             \item{hyperid = }{56407}
#'             \item{name = }{beta7}
#'             \item{short.name = }{beta7}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta8'}{
#'              \describe{
#'             \item{hyperid = }{56408}
#'             \item{name = }{beta8}
#'             \item{short.name = }{beta8}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta9'}{
#'              \describe{
#'             \item{hyperid = }{56409}
#'             \item{name = }{beta9}
#'             \item{short.name = }{beta9}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta10'}{
#'              \describe{
#'             \item{hyperid = }{56410}
#'             \item{name = }{beta10}
#'             \item{short.name = }{beta10}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model '0binomialS'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'New 0-inflated Binomial Swap'}
#'               \item{status = }{'experimental'}
#'               \item{survival = }{'FALSE'}
#'               \item{discrete = }{'TRUE'}
#'               \item{link = }{'default logit loga cauchit probit cloglog ccloglog loglog log'}
#'               \item{link.simple = }{'default logit cauchit probit cloglog ccloglog'}
#'               \item{pdf = }{'0inflated'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 10.
#'          \describe{
#'           \item{Hyperparameter 'theta1'}{
#'              \describe{
#'             \item{hyperid = }{56501}
#'             \item{name = }{beta1}
#'             \item{short.name = }{beta1}
#'             \item{initial = }{-4}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{-4 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta2'}{
#'              \describe{
#'             \item{hyperid = }{56502}
#'             \item{name = }{beta2}
#'             \item{short.name = }{beta2}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta3'}{
#'              \describe{
#'             \item{hyperid = }{56503}
#'             \item{name = }{beta3}
#'             \item{short.name = }{beta3}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta4'}{
#'              \describe{
#'             \item{hyperid = }{56504}
#'             \item{name = }{beta4}
#'             \item{short.name = }{beta4}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta5'}{
#'              \describe{
#'             \item{hyperid = }{56505}
#'             \item{name = }{beta5}
#'             \item{short.name = }{beta5}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta6'}{
#'              \describe{
#'             \item{hyperid = }{56506}
#'             \item{name = }{beta6}
#'             \item{short.name = }{beta6}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta7'}{
#'              \describe{
#'             \item{hyperid = }{56507}
#'             \item{name = }{beta7}
#'             \item{short.name = }{beta7}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta8'}{
#'              \describe{
#'             \item{hyperid = }{56508}
#'             \item{name = }{beta8}
#'             \item{short.name = }{beta8}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta9'}{
#'              \describe{
#'             \item{hyperid = }{56509}
#'             \item{name = }{beta9}
#'             \item{short.name = }{beta9}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta10'}{
#'              \describe{
#'             \item{hyperid = }{56510}
#'             \item{name = }{beta10}
#'             \item{short.name = }{beta10}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'binomial'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'The Binomial likelihood'}
#'               \item{survival = }{'FALSE'}
#'               \item{discrete = }{'TRUE'}
#'               \item{link = }{'default logit loga cauchit probit cloglog ccloglog loglog log sslogit logitoffset quantile pquantile robit sn powerlogit'}
#'               \item{pdf = }{'binomial'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 0.
#'        }
#'       \item{Model 'xbinomial'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'The Binomial likelihood (expert version)'}
#'               \item{survival = }{'FALSE'}
#'               \item{discrete = }{'TRUE'}
#'               \item{link = }{'default logit loga cauchit probit cloglog ccloglog loglog log sslogit logitoffset quantile pquantile robit sn powerlogit'}
#'               \item{pdf = }{'binomial'}
#'               \item{status = }{'experimental'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 0.
#'        }
#'       \item{Model 'pom'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'Likelihood for the proportional odds model'}
#'               \item{status = }{'experimental'}
#'               \item{survival = }{'FALSE'}
#'               \item{discrete = }{'TRUE'}
#'               \item{link = }{'default identity'}
#'               \item{pdf = }{'pom'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 10.
#'          \describe{
#'           \item{Hyperparameter 'theta1'}{
#'              \describe{
#'             \item{hyperid = }{57101}
#'             \item{name = }{theta1}
#'             \item{short.name = }{theta1}
#'             \item{initial = }{NA}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{dirichlet}
#'             \item{param = }{3}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta2'}{
#'              \describe{
#'             \item{hyperid = }{57102}
#'             \item{name = }{theta2}
#'             \item{short.name = }{theta2}
#'             \item{initial = }{NA}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta3'}{
#'              \describe{
#'             \item{hyperid = }{57103}
#'             \item{name = }{theta3}
#'             \item{short.name = }{theta3}
#'             \item{initial = }{NA}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta4'}{
#'              \describe{
#'             \item{hyperid = }{57104}
#'             \item{name = }{theta4}
#'             \item{short.name = }{theta4}
#'             \item{initial = }{NA}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta5'}{
#'              \describe{
#'             \item{hyperid = }{57105}
#'             \item{name = }{theta5}
#'             \item{short.name = }{theta5}
#'             \item{initial = }{NA}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta6'}{
#'              \describe{
#'             \item{hyperid = }{57106}
#'             \item{name = }{theta6}
#'             \item{short.name = }{theta6}
#'             \item{initial = }{NA}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta7'}{
#'              \describe{
#'             \item{hyperid = }{57107}
#'             \item{name = }{theta7}
#'             \item{short.name = }{theta7}
#'             \item{initial = }{NA}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta8'}{
#'              \describe{
#'             \item{hyperid = }{57108}
#'             \item{name = }{theta8}
#'             \item{short.name = }{theta8}
#'             \item{initial = }{NA}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta9'}{
#'              \describe{
#'             \item{hyperid = }{57109}
#'             \item{name = }{theta9}
#'             \item{short.name = }{theta9}
#'             \item{initial = }{NA}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta10'}{
#'              \describe{
#'             \item{hyperid = }{57110}
#'             \item{name = }{theta10}
#'             \item{short.name = }{theta10}
#'             \item{initial = }{NA}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'bgev'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'The blended Generalized Extreme Value likelihood'}
#'               \item{status = }{'experimental'}
#'               \item{survival = }{'FALSE'}
#'               \item{discrete = }{'FALSE'}
#'               \item{link = }{'default identity log'}
#'               \item{pdf = }{'bgev'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 12.
#'          \describe{
#'           \item{Hyperparameter 'theta1'}{
#'              \describe{
#'             \item{hyperid = }{57201}
#'             \item{name = }{spread}
#'             \item{short.name = }{sd}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{loggamma}
#'             \item{param = }{1 3}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta2'}{
#'              \describe{
#'             \item{hyperid = }{57202}
#'             \item{name = }{tail}
#'             \item{short.name = }{xi}
#'             \item{initial = }{-4}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{pc.gevtail}
#'             \item{param = }{7 0 0.5}
#'             \item{to.theta = }{`function(x, interval = c(REPLACE.ME.low, REPLACE.ME.high)) log(-(interval[1] - x) / (interval[2] - x))`}
#'             \item{from.theta = }{`function(x, interval = c(REPLACE.ME.low, REPLACE.ME.high)) interval[1] + (interval[2] - interval[1]) * exp(x) / (1.0 + exp(x))`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta3'}{
#'              \describe{
#'             \item{hyperid = }{57203}
#'             \item{name = }{beta1}
#'             \item{short.name = }{beta1}
#'             \item{initial = }{NA}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 300}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta4'}{
#'              \describe{
#'             \item{hyperid = }{57204}
#'             \item{name = }{beta2}
#'             \item{short.name = }{beta2}
#'             \item{initial = }{NA}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 300}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta5'}{
#'              \describe{
#'             \item{hyperid = }{57205}
#'             \item{name = }{beta3}
#'             \item{short.name = }{beta3}
#'             \item{initial = }{NA}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 300}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta6'}{
#'              \describe{
#'             \item{hyperid = }{57206}
#'             \item{name = }{beta4}
#'             \item{short.name = }{beta4}
#'             \item{initial = }{NA}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 300}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta7'}{
#'              \describe{
#'             \item{hyperid = }{57207}
#'             \item{name = }{beta5}
#'             \item{short.name = }{beta5}
#'             \item{initial = }{NA}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 300}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta8'}{
#'              \describe{
#'             \item{hyperid = }{57208}
#'             \item{name = }{beta6}
#'             \item{short.name = }{beta6}
#'             \item{initial = }{NA}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 300}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta9'}{
#'              \describe{
#'             \item{hyperid = }{57209}
#'             \item{name = }{beta7}
#'             \item{short.name = }{beta7}
#'             \item{initial = }{NA}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 300}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta10'}{
#'              \describe{
#'             \item{hyperid = }{57210}
#'             \item{name = }{beta8}
#'             \item{short.name = }{beta8}
#'             \item{initial = }{NA}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 300}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta11'}{
#'              \describe{
#'             \item{hyperid = }{57211}
#'             \item{name = }{beta9}
#'             \item{short.name = }{beta9}
#'             \item{initial = }{NA}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 300}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta12'}{
#'              \describe{
#'             \item{hyperid = }{57212}
#'             \item{name = }{beta10}
#'             \item{short.name = }{beta}
#'             \item{initial = }{NA}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 300}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'gamma'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'The Gamma likelihood'}
#'               \item{survival = }{'FALSE'}
#'               \item{discrete = }{'FALSE'}
#'               \item{link = }{'default log quantile'}
#'               \item{pdf = }{'gamma'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 1.
#'          \describe{
#'           \item{Hyperparameter 'theta'}{
#'              \describe{
#'             \item{hyperid = }{58001}
#'             \item{name = }{precision parameter}
#'             \item{short.name = }{prec}
#'             \item{initial = }{4.60517018598809}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{loggamma}
#'             \item{param = }{1 0.01}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'gammasurv'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'The Gamma likelihood (survival)'}
#'               \item{survival = }{'TRUE'}
#'               \item{discrete = }{'FALSE'}
#'               \item{status = }{'experimental'}
#'               \item{link = }{'default log neglog quantile'}
#'               \item{pdf = }{'gammasurv'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 11.
#'          \describe{
#'           \item{Hyperparameter 'theta1'}{
#'              \describe{
#'             \item{hyperid = }{58101}
#'             \item{name = }{precision parameter}
#'             \item{short.name = }{prec}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{loggamma}
#'             \item{param = }{1 0.01}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta2'}{
#'              \describe{
#'             \item{hyperid = }{58102}
#'             \item{name = }{beta1}
#'             \item{short.name = }{beta1}
#'             \item{initial = }{-7}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{-4 100}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta3'}{
#'              \describe{
#'             \item{hyperid = }{58103}
#'             \item{name = }{beta2}
#'             \item{short.name = }{beta2}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 100}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta4'}{
#'              \describe{
#'             \item{hyperid = }{58104}
#'             \item{name = }{beta3}
#'             \item{short.name = }{beta3}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 100}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta5'}{
#'              \describe{
#'             \item{hyperid = }{58105}
#'             \item{name = }{beta4}
#'             \item{short.name = }{beta4}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 100}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta6'}{
#'              \describe{
#'             \item{hyperid = }{58106}
#'             \item{name = }{beta5}
#'             \item{short.name = }{beta5}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 100}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta7'}{
#'              \describe{
#'             \item{hyperid = }{58107}
#'             \item{name = }{beta6}
#'             \item{short.name = }{beta6}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 100}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta8'}{
#'              \describe{
#'             \item{hyperid = }{58108}
#'             \item{name = }{beta7}
#'             \item{short.name = }{beta7}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 100}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta9'}{
#'              \describe{
#'             \item{hyperid = }{58109}
#'             \item{name = }{beta8}
#'             \item{short.name = }{beta8}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 100}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta10'}{
#'              \describe{
#'             \item{hyperid = }{58110}
#'             \item{name = }{beta9}
#'             \item{short.name = }{beta9}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 100}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta11'}{
#'              \describe{
#'             \item{hyperid = }{58111}
#'             \item{name = }{beta10}
#'             \item{short.name = }{beta10}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 100}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'gammajw'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'A special case of the Gamma likelihood'}
#'               \item{survival = }{'FALSE'}
#'               \item{discrete = }{'FALSE'}
#'               \item{link = }{'default log neglog'}
#'               \item{pdf = }{'gammajw'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 0.
#'        }
#'       \item{Model 'gammajwsurv'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'A special case of the Gamma likelihood (survival)'}
#'               \item{survival = }{'TRUE'}
#'               \item{discrete = }{'FALSE'}
#'               \item{link = }{'default log'}
#'               \item{pdf = }{'gammajw'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 10.
#'          \describe{
#'           \item{Hyperparameter 'theta1'}{
#'              \describe{
#'             \item{hyperid = }{58200}
#'             \item{name = }{beta1}
#'             \item{short.name = }{beta1}
#'             \item{initial = }{-7}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{-4 100}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta2'}{
#'              \describe{
#'             \item{hyperid = }{58201}
#'             \item{name = }{beta2}
#'             \item{short.name = }{beta2}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 100}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta3'}{
#'              \describe{
#'             \item{hyperid = }{58202}
#'             \item{name = }{beta3}
#'             \item{short.name = }{beta3}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 100}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta4'}{
#'              \describe{
#'             \item{hyperid = }{58203}
#'             \item{name = }{beta4}
#'             \item{short.name = }{beta4}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 100}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta5'}{
#'              \describe{
#'             \item{hyperid = }{58204}
#'             \item{name = }{beta5}
#'             \item{short.name = }{beta5}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 100}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta6'}{
#'              \describe{
#'             \item{hyperid = }{58205}
#'             \item{name = }{beta6}
#'             \item{short.name = }{beta6}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 100}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta7'}{
#'              \describe{
#'             \item{hyperid = }{58206}
#'             \item{name = }{beta7}
#'             \item{short.name = }{beta7}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 100}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta8'}{
#'              \describe{
#'             \item{hyperid = }{58207}
#'             \item{name = }{beta8}
#'             \item{short.name = }{beta8}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 100}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta9'}{
#'              \describe{
#'             \item{hyperid = }{58208}
#'             \item{name = }{beta9}
#'             \item{short.name = }{beta9}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 100}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta10'}{
#'              \describe{
#'             \item{hyperid = }{58209}
#'             \item{name = }{beta10}
#'             \item{short.name = }{beta10}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 100}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'gammacount'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'A Gamma generalisation of the Poisson likelihood'}
#'               \item{survival = }{'FALSE'}
#'               \item{discrete = }{'FALSE'}
#'               \item{link = }{'default log'}
#'               \item{status = }{'experimental'}
#'               \item{pdf = }{'gammacount'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 1.
#'          \describe{
#'           \item{Hyperparameter 'theta'}{
#'              \describe{
#'             \item{hyperid = }{59001}
#'             \item{name = }{log alpha}
#'             \item{short.name = }{alpha}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{pc.gammacount}
#'             \item{param = }{3}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'qkumar'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'A quantile version of the Kumar likelihood'}
#'               \item{survival = }{'FALSE'}
#'               \item{discrete = }{'FALSE'}
#'               \item{link = }{'default logit loga cauchit'}
#'               \item{pdf = }{'qkumar'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 1.
#'          \describe{
#'           \item{Hyperparameter 'theta'}{
#'              \describe{
#'             \item{hyperid = }{60001}
#'             \item{name = }{precision parameter}
#'             \item{short.name = }{prec}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{loggamma}
#'             \item{param = }{1 0.1}
#'             \item{to.theta = }{`function(x, sc = 0.1) log(x) / sc`}
#'             \item{from.theta = }{`function(x, sc = 0.1) exp(sc * x)`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'qloglogistic'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'A quantile loglogistic likelihood'}
#'               \item{survival = }{'FALSE'}
#'               \item{discrete = }{'FALSE'}
#'               \item{link = }{'default log neglog'}
#'               \item{pdf = }{'qloglogistic'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 1.
#'          \describe{
#'           \item{Hyperparameter 'theta'}{
#'              \describe{
#'             \item{hyperid = }{60011}
#'             \item{name = }{log alpha}
#'             \item{short.name = }{alpha}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{loggamma}
#'             \item{param = }{25 25}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'qloglogisticsurv'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'A quantile loglogistic likelihood (survival)'}
#'               \item{survival = }{'TRUE'}
#'               \item{discrete = }{'FALSE'}
#'               \item{link = }{'default log neglog'}
#'               \item{pdf = }{'qloglogistic'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 11.
#'          \describe{
#'           \item{Hyperparameter 'theta1'}{
#'              \describe{
#'             \item{hyperid = }{60021}
#'             \item{name = }{log alpha}
#'             \item{short.name = }{alpha}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{loggamma}
#'             \item{param = }{25 25}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta2'}{
#'              \describe{
#'             \item{hyperid = }{60022}
#'             \item{name = }{beta1}
#'             \item{short.name = }{beta1}
#'             \item{initial = }{-5}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{-4 100}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta3'}{
#'              \describe{
#'             \item{hyperid = }{60023}
#'             \item{name = }{beta2}
#'             \item{short.name = }{beta2}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 100}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta4'}{
#'              \describe{
#'             \item{hyperid = }{60024}
#'             \item{name = }{beta3}
#'             \item{short.name = }{beta3}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 100}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta5'}{
#'              \describe{
#'             \item{hyperid = }{60025}
#'             \item{name = }{beta4}
#'             \item{short.name = }{beta4}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 100}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta6'}{
#'              \describe{
#'             \item{hyperid = }{60026}
#'             \item{name = }{beta5}
#'             \item{short.name = }{beta5}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 100}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta7'}{
#'              \describe{
#'             \item{hyperid = }{60027}
#'             \item{name = }{beta6}
#'             \item{short.name = }{beta6}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 100}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta8'}{
#'              \describe{
#'             \item{hyperid = }{60028}
#'             \item{name = }{beta7}
#'             \item{short.name = }{beta7}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 100}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta9'}{
#'              \describe{
#'             \item{hyperid = }{60029}
#'             \item{name = }{beta8}
#'             \item{short.name = }{beta8}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 100}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta10'}{
#'              \describe{
#'             \item{hyperid = }{60030}
#'             \item{name = }{beta9}
#'             \item{short.name = }{beta9}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 100}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta11'}{
#'              \describe{
#'             \item{hyperid = }{60031}
#'             \item{name = }{beta10}
#'             \item{short.name = }{beta10}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 100}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'beta'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'The Beta likelihood'}
#'               \item{survival = }{'FALSE'}
#'               \item{discrete = }{'FALSE'}
#'               \item{link = }{'default logit loga cauchit probit cloglog ccloglog loglog'}
#'               \item{pdf = }{'beta'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 1.
#'          \describe{
#'           \item{Hyperparameter 'theta'}{
#'              \describe{
#'             \item{hyperid = }{61001}
#'             \item{name = }{precision parameter}
#'             \item{short.name = }{phi}
#'             \item{initial = }{2.30258509299405}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{loggamma}
#'             \item{param = }{1 0.1}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'betabinomial'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'The Beta-Binomial likelihood'}
#'               \item{survival = }{'FALSE'}
#'               \item{discrete = }{'TRUE'}
#'               \item{link = }{'default logit loga cauchit probit cloglog ccloglog loglog robit sn'}
#'               \item{pdf = }{'betabinomial'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 1.
#'          \describe{
#'           \item{Hyperparameter 'theta'}{
#'              \describe{
#'             \item{hyperid = }{62001}
#'             \item{name = }{overdispersion}
#'             \item{short.name = }{rho}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{gaussian}
#'             \item{param = }{0 0.4}
#'             \item{to.theta = }{`function(x) log(x / (1 - x))`}
#'             \item{from.theta = }{`function(x) exp(x) / (1 + exp(x))`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'betabinomialna'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'The Beta-Binomial Normal approximation likelihood'}
#'               \item{survival = }{'FALSE'}
#'               \item{discrete = }{'TRUE'}
#'               \item{link = }{'default logit loga cauchit probit cloglog ccloglog loglog robit sn'}
#'               \item{pdf = }{'betabinomialna'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 1.
#'          \describe{
#'           \item{Hyperparameter 'theta'}{
#'              \describe{
#'             \item{hyperid = }{62101}
#'             \item{name = }{overdispersion}
#'             \item{short.name = }{rho}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{gaussian}
#'             \item{param = }{0 0.4}
#'             \item{to.theta = }{`function(x) log(x / (1 - x))`}
#'             \item{from.theta = }{`function(x) exp(x) / (1 + exp(x))`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'cbinomial'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'The clustered Binomial likelihood'}
#'               \item{survival = }{'FALSE'}
#'               \item{discrete = }{'TRUE'}
#'               \item{link = }{'default logit loga cauchit probit cloglog ccloglog loglog robit sn'}
#'               \item{status = }{'experimental'}
#'               \item{pdf = }{'cbinomial'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 0.
#'        }
#'       \item{Model 'nbinomial'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'The negBinomial likelihood'}
#'               \item{survival = }{'FALSE'}
#'               \item{discrete = }{'TRUE'}
#'               \item{link = }{'default log logoffset quantile'}
#'               \item{pdf = }{'nbinomial'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 1.
#'          \describe{
#'           \item{Hyperparameter 'theta'}{
#'              \describe{
#'             \item{hyperid = }{63001}
#'             \item{name = }{size}
#'             \item{short.name = }{size}
#'             \item{initial = }{2.30258509299405}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{pc.mgamma}
#'             \item{param = }{7}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'nbinomial2'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'The negBinomial2 likelihood'}
#'               \item{survival = }{'FALSE'}
#'               \item{discrete = }{'TRUE'}
#'               \item{link = }{'default logit loga cauchit probit cloglog ccloglog loglog'}
#'               \item{pdf = }{'nbinomial'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 0.
#'        }
#'       \item{Model 'cennbinomial2'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'The CenNegBinomial2 likelihood (similar to cenpoisson2)'}
#'               \item{status = }{'experimental'}
#'               \item{survival = }{'FALSE'}
#'               \item{discrete = }{'TRUE'}
#'               \item{link = }{'default log logoffset quantile'}
#'               \item{pdf = }{'cennbinomial2'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 1.
#'          \describe{
#'           \item{Hyperparameter 'theta'}{
#'              \describe{
#'             \item{hyperid = }{63101}
#'             \item{name = }{size}
#'             \item{short.name = }{size}
#'             \item{initial = }{2.30258509299405}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{pc.mgamma}
#'             \item{param = }{7}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'simplex'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'The simplex likelihood'}
#'               \item{survival = }{'FALSE'}
#'               \item{discrete = }{'FALSE'}
#'               \item{link = }{'default logit loga cauchit probit cloglog ccloglog loglog'}
#'               \item{pdf = }{'simplex'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 1.
#'          \describe{
#'           \item{Hyperparameter 'theta'}{
#'              \describe{
#'             \item{hyperid = }{64001}
#'             \item{name = }{log precision}
#'             \item{short.name = }{prec}
#'             \item{initial = }{4}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{loggamma}
#'             \item{param = }{1 5e-05}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'gaussian'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'The Gaussian likelihoood'}
#'               \item{survival = }{'FALSE'}
#'               \item{discrete = }{'FALSE'}
#'               \item{link = }{'default identity logit loga cauchit log logoffset'}
#'               \item{pdf = }{'gaussian'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 2.
#'          \describe{
#'           \item{Hyperparameter 'theta1'}{
#'              \describe{
#'             \item{hyperid = }{65001}
#'             \item{name = }{log precision}
#'             \item{short.name = }{prec}
#'             \item{initial = }{4}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{loggamma}
#'             \item{param = }{1 5e-05}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta2'}{
#'              \describe{
#'             \item{hyperid = }{65002}
#'             \item{name = }{log precision offset}
#'             \item{short.name = }{precoffset}
#'             \item{initial = }{72.0873067782343}
#'             \item{fixed = }{TRUE}
#'             \item{prior = }{none}
#'             \item{param = }{}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'gaussianjw'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'The GaussianJW likelihoood'}
#'               \item{status = }{'experimental'}
#'               \item{survival = }{'FALSE'}
#'               \item{discrete = }{'FALSE'}
#'               \item{link = }{'default logit probit'}
#'               \item{pdf = }{'gaussianjw'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 3.
#'          \describe{
#'           \item{Hyperparameter 'theta1'}{
#'              \describe{
#'             \item{hyperid = }{65101}
#'             \item{name = }{beta1}
#'             \item{short.name = }{beta1}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 100}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta2'}{
#'              \describe{
#'             \item{hyperid = }{65102}
#'             \item{name = }{beta2}
#'             \item{short.name = }{beta2}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{1 100}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta3'}{
#'              \describe{
#'             \item{hyperid = }{65103}
#'             \item{name = }{beta3}
#'             \item{short.name = }{beta3}
#'             \item{initial = }{-1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{-1 100}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'agaussian'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'The aggregated Gaussian likelihoood'}
#'               \item{status = }{'experimental'}
#'               \item{survival = }{'FALSE'}
#'               \item{discrete = }{'FALSE'}
#'               \item{link = }{'default identity logit loga cauchit log logoffset'}
#'               \item{pdf = }{'agaussian'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 1.
#'          \describe{
#'           \item{Hyperparameter 'theta'}{
#'              \describe{
#'             \item{hyperid = }{66001}
#'             \item{name = }{log precision}
#'             \item{short.name = }{prec}
#'             \item{initial = }{4}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{loggamma}
#'             \item{param = }{1 5e-05}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'circularnormal'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'The circular Gaussian likelihoood'}
#'               \item{survival = }{'FALSE'}
#'               \item{discrete = }{'FALSE'}
#'               \item{link = }{'default tan'}
#'               \item{pdf = }{'circular-normal'}
#'               \item{status = }{'experimental'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 1.
#'          \describe{
#'           \item{Hyperparameter 'theta'}{
#'              \describe{
#'             \item{hyperid = }{67001}
#'             \item{name = }{log precision parameter}
#'             \item{short.name = }{prec}
#'             \item{initial = }{2}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{loggamma}
#'             \item{param = }{1 0.01}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'wrappedcauchy'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'The wrapped Cauchy likelihoood'}
#'               \item{survival = }{'FALSE'}
#'               \item{discrete = }{'FALSE'}
#'               \item{link = }{'default tan'}
#'               \item{pdf = }{'wrapped-cauchy'}
#'               \item{status = }{'disabled'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 1.
#'          \describe{
#'           \item{Hyperparameter 'theta'}{
#'              \describe{
#'             \item{hyperid = }{68001}
#'             \item{name = }{log precision parameter}
#'             \item{short.name = }{prec}
#'             \item{initial = }{2}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{loggamma}
#'             \item{param = }{1 0.005}
#'             \item{to.theta = }{`function(x) log(x / (1 - x))`}
#'             \item{from.theta = }{`function(x) exp(x) / (1 + exp(x))`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'iidgamma'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'(experimental)'}
#'               \item{survival = }{'FALSE'}
#'               \item{discrete = }{'FALSE'}
#'               \item{link = }{'default identity'}
#'               \item{pdf = }{'iidgamma'}
#'               \item{status = }{'experimental'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 2.
#'          \describe{
#'           \item{Hyperparameter 'theta1'}{
#'              \describe{
#'             \item{hyperid = }{69001}
#'             \item{name = }{logshape}
#'             \item{short.name = }{shape}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{loggamma}
#'             \item{param = }{100 100}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta2'}{
#'              \describe{
#'             \item{hyperid = }{69002}
#'             \item{name = }{lograte}
#'             \item{short.name = }{rate}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{loggamma}
#'             \item{param = }{100 100}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'iidlogitbeta'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'(experimental)'}
#'               \item{survival = }{'FALSE'}
#'               \item{discrete = }{'FALSE'}
#'               \item{link = }{'default logit loga'}
#'               \item{pdf = }{'iidlogitbeta'}
#'               \item{status = }{'experimental'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 2.
#'          \describe{
#'           \item{Hyperparameter 'theta1'}{
#'              \describe{
#'             \item{hyperid = }{70001}
#'             \item{name = }{log.a}
#'             \item{short.name = }{a}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{loggamma}
#'             \item{param = }{1 1}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta2'}{
#'              \describe{
#'             \item{hyperid = }{70002}
#'             \item{name = }{log.b}
#'             \item{short.name = }{b}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{loggamma}
#'             \item{param = }{1 1}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'loggammafrailty'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'(experimental)'}
#'               \item{survival = }{'FALSE'}
#'               \item{discrete = }{'FALSE'}
#'               \item{link = }{'default identity'}
#'               \item{pdf = }{'loggammafrailty'}
#'               \item{status = }{'experimental'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 1.
#'          \describe{
#'           \item{Hyperparameter 'theta'}{
#'              \describe{
#'             \item{hyperid = }{71001}
#'             \item{name = }{log precision}
#'             \item{short.name = }{prec}
#'             \item{initial = }{4}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{loggamma}
#'             \item{param = }{1 5e-05}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'logistic'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'The Logistic likelihoood'}
#'               \item{survival = }{'FALSE'}
#'               \item{discrete = }{'FALSE'}
#'               \item{link = }{'default identity'}
#'               \item{pdf = }{'logistic'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 1.
#'          \describe{
#'           \item{Hyperparameter 'theta'}{
#'              \describe{
#'             \item{hyperid = }{72001}
#'             \item{name = }{log precision}
#'             \item{short.name = }{prec}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{loggamma}
#'             \item{param = }{1 5e-05}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'sn'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'The Skew-Normal likelihoood'}
#'               \item{status = }{'experimental'}
#'               \item{survival = }{'FALSE'}
#'               \item{discrete = }{'FALSE'}
#'               \item{link = }{'default identity'}
#'               \item{pdf = }{'sn'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 2.
#'          \describe{
#'           \item{Hyperparameter 'theta1'}{
#'              \describe{
#'             \item{hyperid = }{74001}
#'             \item{name = }{log precision}
#'             \item{short.name = }{prec}
#'             \item{initial = }{4}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{loggamma}
#'             \item{param = }{1 5e-05}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta2'}{
#'              \describe{
#'             \item{hyperid = }{74002}
#'             \item{name = }{logit skew}
#'             \item{short.name = }{skew}
#'             \item{initial = }{0.00123456789}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{pc.sn}
#'             \item{param = }{10}
#'             \item{to.theta = }{`function(x, skew.max = 0.988) log((1 + x / skew.max) / (1 - x / skew.max))`}
#'             \item{from.theta = }{`function(x, skew.max = 0.988) skew.max * (2 * exp(x) / (1 + exp(x)) - 1)`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'gev'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'The Generalized Extreme Value likelihood'}
#'               \item{survival = }{'FALSE'}
#'               \item{discrete = }{'FALSE'}
#'               \item{link = }{'default identity'}
#'               \item{status = }{'disabled: Use likelihood model 'bgev' instead; see inla.doc('bgev')'}
#'               \item{pdf = }{'gev'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 2.
#'          \describe{
#'           \item{Hyperparameter 'theta1'}{
#'              \describe{
#'             \item{hyperid = }{76001}
#'             \item{name = }{log precision}
#'             \item{short.name = }{prec}
#'             \item{initial = }{4}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{loggamma}
#'             \item{param = }{1 5e-05}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta2'}{
#'              \describe{
#'             \item{hyperid = }{76002}
#'             \item{name = }{tail parameter}
#'             \item{short.name = }{tail}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{gaussian}
#'             \item{param = }{0 25}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'lognormal'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'The log-Normal likelihood'}
#'               \item{survival = }{'FALSE'}
#'               \item{discrete = }{'FALSE'}
#'               \item{link = }{'default identity'}
#'               \item{pdf = }{'lognormal'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 1.
#'          \describe{
#'           \item{Hyperparameter 'theta'}{
#'              \describe{
#'             \item{hyperid = }{77101}
#'             \item{name = }{log precision}
#'             \item{short.name = }{prec}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{loggamma}
#'             \item{param = }{1 5e-05}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'lognormalsurv'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'The log-Normal likelihood (survival)'}
#'               \item{survival = }{'TRUE'}
#'               \item{discrete = }{'FALSE'}
#'               \item{link = }{'default identity'}
#'               \item{pdf = }{'lognormal'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 11.
#'          \describe{
#'           \item{Hyperparameter 'theta'}{
#'              \describe{
#'             \item{hyperid = }{78001}
#'             \item{name = }{log precision}
#'             \item{short.name = }{prec}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{loggamma}
#'             \item{param = }{1 5e-05}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta2'}{
#'              \describe{
#'             \item{hyperid = }{78002}
#'             \item{name = }{beta1}
#'             \item{short.name = }{beta1}
#'             \item{initial = }{-7}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{-4 100}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta3'}{
#'              \describe{
#'             \item{hyperid = }{78003}
#'             \item{name = }{beta2}
#'             \item{short.name = }{beta2}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 100}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta4'}{
#'              \describe{
#'             \item{hyperid = }{78004}
#'             \item{name = }{beta3}
#'             \item{short.name = }{beta3}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 100}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta5'}{
#'              \describe{
#'             \item{hyperid = }{78005}
#'             \item{name = }{beta4}
#'             \item{short.name = }{beta4}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 100}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta6'}{
#'              \describe{
#'             \item{hyperid = }{78006}
#'             \item{name = }{beta5}
#'             \item{short.name = }{beta5}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 100}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta7'}{
#'              \describe{
#'             \item{hyperid = }{78007}
#'             \item{name = }{beta6}
#'             \item{short.name = }{beta6}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 100}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta8'}{
#'              \describe{
#'             \item{hyperid = }{78008}
#'             \item{name = }{beta7}
#'             \item{short.name = }{beta7}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 100}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta9'}{
#'              \describe{
#'             \item{hyperid = }{78009}
#'             \item{name = }{beta8}
#'             \item{short.name = }{beta8}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 100}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta10'}{
#'              \describe{
#'             \item{hyperid = }{78010}
#'             \item{name = }{beta9}
#'             \item{short.name = }{beta9}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 100}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta11'}{
#'              \describe{
#'             \item{hyperid = }{78011}
#'             \item{name = }{beta10}
#'             \item{short.name = }{beta10}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 100}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'exponential'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'The Exponential likelihood'}
#'               \item{survival = }{'FALSE'}
#'               \item{discrete = }{'FALSE'}
#'               \item{link = }{'default log'}
#'               \item{pdf = }{'exponential'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 0.
#'        }
#'       \item{Model 'exponentialsurv'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'The Exponential likelihood (survival)'}
#'               \item{survival = }{'TRUE'}
#'               \item{discrete = }{'FALSE'}
#'               \item{link = }{'default log neglog'}
#'               \item{pdf = }{'exponential'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 10.
#'          \describe{
#'           \item{Hyperparameter 'theta1'}{
#'              \describe{
#'             \item{hyperid = }{78020}
#'             \item{name = }{beta1}
#'             \item{short.name = }{beta1}
#'             \item{initial = }{-4}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{-1 100}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta2'}{
#'              \describe{
#'             \item{hyperid = }{78021}
#'             \item{name = }{beta2}
#'             \item{short.name = }{beta2}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 100}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta3'}{
#'              \describe{
#'             \item{hyperid = }{78022}
#'             \item{name = }{beta3}
#'             \item{short.name = }{beta3}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 100}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta4'}{
#'              \describe{
#'             \item{hyperid = }{78023}
#'             \item{name = }{beta4}
#'             \item{short.name = }{beta4}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 100}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta5'}{
#'              \describe{
#'             \item{hyperid = }{78024}
#'             \item{name = }{beta5}
#'             \item{short.name = }{beta5}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 100}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta6'}{
#'              \describe{
#'             \item{hyperid = }{78025}
#'             \item{name = }{beta6}
#'             \item{short.name = }{beta6}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 100}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta7'}{
#'              \describe{
#'             \item{hyperid = }{78026}
#'             \item{name = }{beta7}
#'             \item{short.name = }{beta7}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 100}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta8'}{
#'              \describe{
#'             \item{hyperid = }{78027}
#'             \item{name = }{beta8}
#'             \item{short.name = }{beta8}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 100}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta9'}{
#'              \describe{
#'             \item{hyperid = }{78028}
#'             \item{name = }{beta9}
#'             \item{short.name = }{beta9}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 100}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta10'}{
#'              \describe{
#'             \item{hyperid = }{78029}
#'             \item{name = }{beta10}
#'             \item{short.name = }{beta10}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 100}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'coxph'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'Cox-proportional hazard likelihood'}
#'               \item{survival = }{'TRUE'}
#'               \item{discrete = }{'TRUE'}
#'               \item{link = }{'default log neglog'}
#'               \item{pdf = }{'coxph'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 0.
#'        }
#'       \item{Model 'weibull'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'The Weibull likelihood'}
#'               \item{survival = }{'FALSE'}
#'               \item{discrete = }{'FALSE'}
#'               \item{link = }{'default log neglog quantile'}
#'               \item{pdf = }{'weibull'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 1.
#'          \describe{
#'           \item{Hyperparameter 'theta'}{
#'              \describe{
#'             \item{hyperid = }{79001}
#'             \item{name = }{log alpha}
#'             \item{short.name = }{alpha}
#'             \item{initial = }{-2}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{pc.alphaw}
#'             \item{param = }{5}
#'             \item{to.theta = }{`function(x, sc = 0.1) log(x) / sc`}
#'             \item{from.theta = }{`function(x, sc = 0.1) exp(sc * x)`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'weibullsurv'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'The Weibull likelihood (survival)'}
#'               \item{survival = }{'TRUE'}
#'               \item{discrete = }{'FALSE'}
#'               \item{link = }{'default log neglog quantile'}
#'               \item{pdf = }{'weibull'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 11.
#'          \describe{
#'           \item{Hyperparameter 'theta'}{
#'              \describe{
#'             \item{hyperid = }{79101}
#'             \item{name = }{log alpha}
#'             \item{short.name = }{alpha}
#'             \item{initial = }{-2}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{pc.alphaw}
#'             \item{param = }{5}
#'             \item{to.theta = }{`function(x, sc = 0.1) log(x) / sc`}
#'             \item{from.theta = }{`function(x, sc = 0.1) exp(sc * x)`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta2'}{
#'              \describe{
#'             \item{hyperid = }{79102}
#'             \item{name = }{beta1}
#'             \item{short.name = }{beta1}
#'             \item{initial = }{-7}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{-4 100}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta3'}{
#'              \describe{
#'             \item{hyperid = }{79103}
#'             \item{name = }{beta2}
#'             \item{short.name = }{beta2}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 100}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta4'}{
#'              \describe{
#'             \item{hyperid = }{79104}
#'             \item{name = }{beta3}
#'             \item{short.name = }{beta3}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 100}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta5'}{
#'              \describe{
#'             \item{hyperid = }{79105}
#'             \item{name = }{beta4}
#'             \item{short.name = }{beta4}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 100}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta6'}{
#'              \describe{
#'             \item{hyperid = }{79106}
#'             \item{name = }{beta5}
#'             \item{short.name = }{beta5}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 100}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta7'}{
#'              \describe{
#'             \item{hyperid = }{79107}
#'             \item{name = }{beta6}
#'             \item{short.name = }{beta6}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 100}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta8'}{
#'              \describe{
#'             \item{hyperid = }{79108}
#'             \item{name = }{beta7}
#'             \item{short.name = }{beta7}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 100}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta9'}{
#'              \describe{
#'             \item{hyperid = }{79109}
#'             \item{name = }{beta8}
#'             \item{short.name = }{beta8}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 100}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta10'}{
#'              \describe{
#'             \item{hyperid = }{79110}
#'             \item{name = }{beta9}
#'             \item{short.name = }{beta9}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 100}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta11'}{
#'              \describe{
#'             \item{hyperid = }{79111}
#'             \item{name = }{beta10}
#'             \item{short.name = }{beta10}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 100}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'loglogistic'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'The loglogistic likelihood'}
#'               \item{survival = }{'FALSE'}
#'               \item{discrete = }{'FALSE'}
#'               \item{link = }{'default log neglog'}
#'               \item{pdf = }{'loglogistic'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 1.
#'          \describe{
#'           \item{Hyperparameter 'theta'}{
#'              \describe{
#'             \item{hyperid = }{80001}
#'             \item{name = }{log alpha}
#'             \item{short.name = }{alpha}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{loggamma}
#'             \item{param = }{25 25}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'loglogisticsurv'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'The loglogistic likelihood (survival)'}
#'               \item{survival = }{'TRUE'}
#'               \item{discrete = }{'FALSE'}
#'               \item{link = }{'default log neglog'}
#'               \item{pdf = }{'loglogistic'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 11.
#'          \describe{
#'           \item{Hyperparameter 'theta1'}{
#'              \describe{
#'             \item{hyperid = }{80011}
#'             \item{name = }{log alpha}
#'             \item{short.name = }{alpha}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{loggamma}
#'             \item{param = }{25 25}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta2'}{
#'              \describe{
#'             \item{hyperid = }{80012}
#'             \item{name = }{beta1}
#'             \item{short.name = }{beta1}
#'             \item{initial = }{-5}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{-4 100}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta3'}{
#'              \describe{
#'             \item{hyperid = }{80013}
#'             \item{name = }{beta2}
#'             \item{short.name = }{beta2}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 100}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta4'}{
#'              \describe{
#'             \item{hyperid = }{80014}
#'             \item{name = }{beta3}
#'             \item{short.name = }{beta3}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 100}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta5'}{
#'              \describe{
#'             \item{hyperid = }{80015}
#'             \item{name = }{beta4}
#'             \item{short.name = }{beta4}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 100}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta6'}{
#'              \describe{
#'             \item{hyperid = }{80016}
#'             \item{name = }{beta5}
#'             \item{short.name = }{beta5}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 100}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta7'}{
#'              \describe{
#'             \item{hyperid = }{80017}
#'             \item{name = }{beta6}
#'             \item{short.name = }{beta6}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 100}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta8'}{
#'              \describe{
#'             \item{hyperid = }{80018}
#'             \item{name = }{beta7}
#'             \item{short.name = }{beta7}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 100}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta9'}{
#'              \describe{
#'             \item{hyperid = }{80019}
#'             \item{name = }{beta8}
#'             \item{short.name = }{beta8}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 100}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta10'}{
#'              \describe{
#'             \item{hyperid = }{80020}
#'             \item{name = }{beta9}
#'             \item{short.name = }{beta9}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 100}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta11'}{
#'              \describe{
#'             \item{hyperid = }{80021}
#'             \item{name = }{beta10}
#'             \item{short.name = }{beta10}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 100}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'stochvol'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'The Gaussian stochvol likelihood'}
#'               \item{survival = }{'FALSE'}
#'               \item{discrete = }{'FALSE'}
#'               \item{link = }{'default log'}
#'               \item{pdf = }{'stochvolgaussian'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 1.
#'          \describe{
#'           \item{Hyperparameter 'theta'}{
#'              \describe{
#'             \item{hyperid = }{82001}
#'             \item{name = }{log precision}
#'             \item{short.name = }{prec}
#'             \item{initial = }{500}
#'             \item{fixed = }{TRUE}
#'             \item{prior = }{loggamma}
#'             \item{param = }{1 0.005}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'stochvolsn'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'The SkewNormal stochvol likelihood'}
#'               \item{status = }{'experimental'}
#'               \item{survival = }{'FALSE'}
#'               \item{discrete = }{'FALSE'}
#'               \item{link = }{'default log'}
#'               \item{pdf = }{'stochvolsn'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 2.
#'          \describe{
#'           \item{Hyperparameter 'theta1'}{
#'              \describe{
#'             \item{hyperid = }{82101}
#'             \item{name = }{logit skew}
#'             \item{short.name = }{skew}
#'             \item{initial = }{0.00123456789}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{pc.sn}
#'             \item{param = }{10}
#'             \item{to.theta = }{`function(x, skew.max = 0.988) log((1 + x / skew.max) / (1 - x / skew.max))`}
#'             \item{from.theta = }{`function(x, skew.max = 0.988) skew.max * (2 * exp(x) / (1 + exp(x)) - 1)`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta2'}{
#'              \describe{
#'             \item{hyperid = }{82102}
#'             \item{name = }{log precision}
#'             \item{short.name = }{prec}
#'             \item{initial = }{500}
#'             \item{fixed = }{TRUE}
#'             \item{prior = }{loggamma}
#'             \item{param = }{1 0.005}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'stochvolt'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'The Student-t stochvol likelihood'}
#'               \item{survival = }{'FALSE'}
#'               \item{discrete = }{'FALSE'}
#'               \item{link = }{'default log'}
#'               \item{pdf = }{'stochvolt'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 1.
#'          \describe{
#'           \item{Hyperparameter 'theta'}{
#'              \describe{
#'             \item{hyperid = }{83001}
#'             \item{name = }{log degrees of freedom}
#'             \item{short.name = }{dof}
#'             \item{initial = }{4}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{pc.dof}
#'             \item{param = }{15 0.5}
#'             \item{to.theta = }{`function(x) log(x - 2)`}
#'             \item{from.theta = }{`function(x) 2 + exp(x)`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'stochvolnig'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'The Normal inverse Gaussian stochvol likelihood'}
#'               \item{survival = }{'FALSE'}
#'               \item{discrete = }{'FALSE'}
#'               \item{link = }{'default log'}
#'               \item{pdf = }{'stochvolnig'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 2.
#'          \describe{
#'           \item{Hyperparameter 'theta1'}{
#'              \describe{
#'             \item{hyperid = }{84001}
#'             \item{name = }{skewness}
#'             \item{short.name = }{skew}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{gaussian}
#'             \item{param = }{0 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta2'}{
#'              \describe{
#'             \item{hyperid = }{84002}
#'             \item{name = }{shape}
#'             \item{short.name = }{shape}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{loggamma}
#'             \item{param = }{1 0.5}
#'             \item{to.theta = }{`function(x) log(x - 1)`}
#'             \item{from.theta = }{`function(x) 1 + exp(x)`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'zeroinflatedpoisson0'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'Zero-inflated Poisson, type 0'}
#'               \item{survival = }{'FALSE'}
#'               \item{discrete = }{'FALSE'}
#'               \item{link = }{'default log'}
#'               \item{pdf = }{'zeroinflated'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 1.
#'          \describe{
#'           \item{Hyperparameter 'theta'}{
#'              \describe{
#'             \item{hyperid = }{85001}
#'             \item{name = }{logit probability}
#'             \item{short.name = }{prob}
#'             \item{initial = }{-1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{gaussian}
#'             \item{param = }{-1 0.2}
#'             \item{to.theta = }{`function(x) log(x / (1 - x))`}
#'             \item{from.theta = }{`function(x) exp(x) / (1 + exp(x))`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'zeroinflatedpoisson1'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'Zero-inflated Poisson, type 1'}
#'               \item{survival = }{'FALSE'}
#'               \item{discrete = }{'FALSE'}
#'               \item{link = }{'default log'}
#'               \item{pdf = }{'zeroinflated'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 1.
#'          \describe{
#'           \item{Hyperparameter 'theta'}{
#'              \describe{
#'             \item{hyperid = }{86001}
#'             \item{name = }{logit probability}
#'             \item{short.name = }{prob}
#'             \item{initial = }{-1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{gaussian}
#'             \item{param = }{-1 0.2}
#'             \item{to.theta = }{`function(x) log(x / (1 - x))`}
#'             \item{from.theta = }{`function(x) exp(x) / (1 + exp(x))`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'zeroinflatedpoisson2'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'Zero-inflated Poisson, type 2'}
#'               \item{survival = }{'FALSE'}
#'               \item{discrete = }{'FALSE'}
#'               \item{link = }{'default log'}
#'               \item{pdf = }{'zeroinflated'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 1.
#'          \describe{
#'           \item{Hyperparameter 'theta'}{
#'              \describe{
#'             \item{hyperid = }{87001}
#'             \item{name = }{log alpha}
#'             \item{short.name = }{a}
#'             \item{initial = }{0.693147180559945}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{gaussian}
#'             \item{param = }{0.693147180559945 1}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'zeroinflatedcenpoisson0'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'Zero-inflated censored Poisson, type 0'}
#'               \item{status = }{'experimental'}
#'               \item{survival = }{'FALSE'}
#'               \item{discrete = }{'FALSE'}
#'               \item{link = }{'default log'}
#'               \item{pdf = }{'zeroinflated'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 1.
#'          \describe{
#'           \item{Hyperparameter 'theta'}{
#'              \describe{
#'             \item{hyperid = }{87101}
#'             \item{name = }{logit probability}
#'             \item{short.name = }{prob}
#'             \item{initial = }{-1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{gaussian}
#'             \item{param = }{-1 0.2}
#'             \item{to.theta = }{`function(x) log(x / (1 - x))`}
#'             \item{from.theta = }{`function(x) exp(x) / (1 + exp(x))`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'zeroinflatedcenpoisson1'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'Zero-inflated censored Poisson, type 1'}
#'               \item{status = }{'experimental'}
#'               \item{survival = }{'FALSE'}
#'               \item{discrete = }{'FALSE'}
#'               \item{link = }{'default log'}
#'               \item{pdf = }{'zeroinflated'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 1.
#'          \describe{
#'           \item{Hyperparameter 'theta'}{
#'              \describe{
#'             \item{hyperid = }{87201}
#'             \item{name = }{logit probability}
#'             \item{short.name = }{prob}
#'             \item{initial = }{-1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{gaussian}
#'             \item{param = }{-1 0.2}
#'             \item{to.theta = }{`function(x) log(x / (1 - x))`}
#'             \item{from.theta = }{`function(x) exp(x) / (1 + exp(x))`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'zeroinflatedbetabinomial0'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'Zero-inflated Beta-Binomial, type 0'}
#'               \item{survival = }{'FALSE'}
#'               \item{discrete = }{'TRUE'}
#'               \item{link = }{'default logit loga cauchit probit cloglog ccloglog loglog robit sn'}
#'               \item{pdf = }{'zeroinflated'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 2.
#'          \describe{
#'           \item{Hyperparameter 'theta1'}{
#'              \describe{
#'             \item{hyperid = }{88001}
#'             \item{name = }{overdispersion}
#'             \item{short.name = }{rho}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{gaussian}
#'             \item{param = }{0 0.4}
#'             \item{to.theta = }{`function(x) log(x / (1 - x))`}
#'             \item{from.theta = }{`function(x) exp(x) / (1 + exp(x))`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta2'}{
#'              \describe{
#'             \item{hyperid = }{88002}
#'             \item{name = }{logit probability}
#'             \item{short.name = }{prob}
#'             \item{initial = }{-1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{gaussian}
#'             \item{param = }{-1 0.2}
#'             \item{to.theta = }{`function(x) log(x / (1 - x))`}
#'             \item{from.theta = }{`function(x) exp(x) / (1 + exp(x))`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'zeroinflatedbetabinomial1'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'Zero-inflated Beta-Binomial, type 1'}
#'               \item{survival = }{'FALSE'}
#'               \item{discrete = }{'TRUE'}
#'               \item{link = }{'default logit loga cauchit probit cloglog ccloglog loglog robit sn'}
#'               \item{pdf = }{'zeroinflated'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 2.
#'          \describe{
#'           \item{Hyperparameter 'theta1'}{
#'              \describe{
#'             \item{hyperid = }{89001}
#'             \item{name = }{overdispersion}
#'             \item{short.name = }{rho}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{gaussian}
#'             \item{param = }{0 0.4}
#'             \item{to.theta = }{`function(x) log(x / (1 - x))`}
#'             \item{from.theta = }{`function(x) exp(x) / (1 + exp(x))`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta2'}{
#'              \describe{
#'             \item{hyperid = }{89002}
#'             \item{name = }{logit probability}
#'             \item{short.name = }{prob}
#'             \item{initial = }{-1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{gaussian}
#'             \item{param = }{-1 0.2}
#'             \item{to.theta = }{`function(x) log(x / (1 - x))`}
#'             \item{from.theta = }{`function(x) exp(x) / (1 + exp(x))`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'zeroinflatedbinomial0'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'Zero-inflated Binomial, type 0'}
#'               \item{survival = }{'FALSE'}
#'               \item{discrete = }{'FALSE'}
#'               \item{link = }{'default logit loga cauchit probit cloglog ccloglog loglog robit sn'}
#'               \item{pdf = }{'zeroinflated'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 1.
#'          \describe{
#'           \item{Hyperparameter 'theta'}{
#'              \describe{
#'             \item{hyperid = }{90001}
#'             \item{name = }{logit probability}
#'             \item{short.name = }{prob}
#'             \item{initial = }{-1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{gaussian}
#'             \item{param = }{-1 0.2}
#'             \item{to.theta = }{`function(x) log(x / (1 - x))`}
#'             \item{from.theta = }{`function(x) exp(x) / (1 + exp(x))`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'zeroinflatedbinomial1'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'Zero-inflated Binomial, type 1'}
#'               \item{survival = }{'FALSE'}
#'               \item{discrete = }{'FALSE'}
#'               \item{link = }{'default logit loga cauchit probit cloglog ccloglog loglog robit sn'}
#'               \item{pdf = }{'zeroinflated'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 1.
#'          \describe{
#'           \item{Hyperparameter 'theta'}{
#'              \describe{
#'             \item{hyperid = }{91001}
#'             \item{name = }{logit probability}
#'             \item{short.name = }{prob}
#'             \item{initial = }{-1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{gaussian}
#'             \item{param = }{-1 0.2}
#'             \item{to.theta = }{`function(x) log(x / (1 - x))`}
#'             \item{from.theta = }{`function(x) exp(x) / (1 + exp(x))`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'zeroinflatedbinomial2'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'Zero-inflated Binomial, type 2'}
#'               \item{survival = }{'FALSE'}
#'               \item{discrete = }{'FALSE'}
#'               \item{link = }{'default logit loga cauchit probit cloglog ccloglog loglog robit sn'}
#'               \item{pdf = }{'zeroinflated'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 1.
#'          \describe{
#'           \item{Hyperparameter 'theta'}{
#'              \describe{
#'             \item{hyperid = }{92001}
#'             \item{name = }{alpha}
#'             \item{short.name = }{alpha}
#'             \item{initial = }{-1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{gaussian}
#'             \item{param = }{-1 0.2}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'zeroninflatedbinomial2'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'Zero and N inflated binomial, type 2'}
#'               \item{survival = }{'FALSE'}
#'               \item{discrete = }{'FALSE'}
#'               \item{link = }{'default logit loga cauchit probit cloglog ccloglog loglog robit sn'}
#'               \item{pdf = }{'NA'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 2.
#'          \describe{
#'           \item{Hyperparameter 'theta1'}{
#'              \describe{
#'             \item{hyperid = }{93001}
#'             \item{name = }{alpha1}
#'             \item{short.name = }{alpha1}
#'             \item{initial = }{-1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{gaussian}
#'             \item{param = }{-1 0.2}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta2'}{
#'              \describe{
#'             \item{hyperid = }{93002}
#'             \item{name = }{alpha2}
#'             \item{short.name = }{alpha2}
#'             \item{initial = }{-1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{gaussian}
#'             \item{param = }{-1 0.2}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'zeroninflatedbinomial3'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'Zero and N inflated binomial, type 3'}
#'               \item{status = }{'experimental'}
#'               \item{survival = }{'FALSE'}
#'               \item{discrete = }{'FALSE'}
#'               \item{link = }{'default logit loga cauchit probit cloglog ccloglog loglog robit sn'}
#'               \item{pdf = }{'zeroinflated'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 2.
#'          \describe{
#'           \item{Hyperparameter 'theta1'}{
#'              \describe{
#'             \item{hyperid = }{93101}
#'             \item{name = }{alpha0}
#'             \item{short.name = }{alpha0}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{loggamma}
#'             \item{param = }{1 1}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta2'}{
#'              \describe{
#'             \item{hyperid = }{93102}
#'             \item{name = }{alphaN}
#'             \item{short.name = }{alphaN}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{loggamma}
#'             \item{param = }{1 1}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'zeroinflatedbetabinomial2'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'Zero inflated Beta-Binomial, type 2'}
#'               \item{survival = }{'FALSE'}
#'               \item{discrete = }{'FALSE'}
#'               \item{link = }{'default logit loga cauchit probit cloglog ccloglog loglog robit sn'}
#'               \item{pdf = }{'zeroinflated'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 2.
#'          \describe{
#'           \item{Hyperparameter 'theta1'}{
#'              \describe{
#'             \item{hyperid = }{94001}
#'             \item{name = }{log alpha}
#'             \item{short.name = }{a}
#'             \item{initial = }{0.693147180559945}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{gaussian}
#'             \item{param = }{0.693147180559945 1}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta2'}{
#'              \describe{
#'             \item{hyperid = }{94002}
#'             \item{name = }{beta}
#'             \item{short.name = }{b}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{gaussian}
#'             \item{param = }{0 1}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'zeroinflatednbinomial0'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'Zero inflated negBinomial, type 0'}
#'               \item{survival = }{'FALSE'}
#'               \item{discrete = }{'FALSE'}
#'               \item{link = }{'default log'}
#'               \item{pdf = }{'zeroinflated'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 2.
#'          \describe{
#'           \item{Hyperparameter 'theta1'}{
#'              \describe{
#'             \item{hyperid = }{95001}
#'             \item{name = }{log size}
#'             \item{short.name = }{size}
#'             \item{initial = }{2.30258509299405}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{pc.mgamma}
#'             \item{param = }{7}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta2'}{
#'              \describe{
#'             \item{hyperid = }{95002}
#'             \item{name = }{logit probability}
#'             \item{short.name = }{prob}
#'             \item{initial = }{-1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{gaussian}
#'             \item{param = }{-1 0.2}
#'             \item{to.theta = }{`function(x) log(x / (1 - x))`}
#'             \item{from.theta = }{`function(x) exp(x) / (1 + exp(x))`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'zeroinflatednbinomial1'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'Zero inflated negBinomial, type 1'}
#'               \item{survival = }{'FALSE'}
#'               \item{discrete = }{'FALSE'}
#'               \item{link = }{'default log'}
#'               \item{pdf = }{'zeroinflated'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 2.
#'          \describe{
#'           \item{Hyperparameter 'theta1'}{
#'              \describe{
#'             \item{hyperid = }{96001}
#'             \item{name = }{log size}
#'             \item{short.name = }{size}
#'             \item{initial = }{2.30258509299405}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{pc.mgamma}
#'             \item{param = }{7}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta2'}{
#'              \describe{
#'             \item{hyperid = }{96002}
#'             \item{name = }{logit probability}
#'             \item{short.name = }{prob}
#'             \item{initial = }{-1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{gaussian}
#'             \item{param = }{-1 0.2}
#'             \item{to.theta = }{`function(x) log(x / (1 - x))`}
#'             \item{from.theta = }{`function(x) exp(x) / (1 + exp(x))`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'zeroinflatednbinomial1strata2'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'Zero inflated negBinomial, type 1, strata 2'}
#'               \item{status = }{'experimental'}
#'               \item{survival = }{'FALSE'}
#'               \item{discrete = }{'FALSE'}
#'               \item{link = }{'default log'}
#'               \item{pdf = }{'zeroinflated'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 11.
#'          \describe{
#'           \item{Hyperparameter 'theta1'}{
#'              \describe{
#'             \item{hyperid = }{97001}
#'             \item{name = }{log size}
#'             \item{short.name = }{size}
#'             \item{initial = }{2.30258509299405}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{pc.mgamma}
#'             \item{param = }{7}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta2'}{
#'              \describe{
#'             \item{hyperid = }{97002}
#'             \item{name = }{logit probability 1}
#'             \item{short.name = }{prob1}
#'             \item{initial = }{-1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{gaussian}
#'             \item{param = }{-1 0.2}
#'             \item{to.theta = }{`function(x) log(x / (1 - x))`}
#'             \item{from.theta = }{`function(x) exp(x) / (1 + exp(x))`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta3'}{
#'              \describe{
#'             \item{hyperid = }{97003}
#'             \item{name = }{logit probability 2}
#'             \item{short.name = }{prob2}
#'             \item{initial = }{-1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{gaussian}
#'             \item{param = }{-1 0.2}
#'             \item{to.theta = }{`function(x) log(x / (1 - x))`}
#'             \item{from.theta = }{`function(x) exp(x) / (1 + exp(x))`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta4'}{
#'              \describe{
#'             \item{hyperid = }{97004}
#'             \item{name = }{logit probability 3}
#'             \item{short.name = }{prob3}
#'             \item{initial = }{-1}
#'             \item{fixed = }{TRUE}
#'             \item{prior = }{gaussian}
#'             \item{param = }{-1 0.2}
#'             \item{to.theta = }{`function(x) log(x / (1 - x))`}
#'             \item{from.theta = }{`function(x) exp(x) / (1 + exp(x))`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta5'}{
#'              \describe{
#'             \item{hyperid = }{97005}
#'             \item{name = }{logit probability 4}
#'             \item{short.name = }{prob4}
#'             \item{initial = }{-1}
#'             \item{fixed = }{TRUE}
#'             \item{prior = }{gaussian}
#'             \item{param = }{-1 0.2}
#'             \item{to.theta = }{`function(x) log(x / (1 - x))`}
#'             \item{from.theta = }{`function(x) exp(x) / (1 + exp(x))`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta6'}{
#'              \describe{
#'             \item{hyperid = }{97006}
#'             \item{name = }{logit probability 5}
#'             \item{short.name = }{prob5}
#'             \item{initial = }{-1}
#'             \item{fixed = }{TRUE}
#'             \item{prior = }{gaussian}
#'             \item{param = }{-1 0.2}
#'             \item{to.theta = }{`function(x) log(x / (1 - x))`}
#'             \item{from.theta = }{`function(x) exp(x) / (1 + exp(x))`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta7'}{
#'              \describe{
#'             \item{hyperid = }{97007}
#'             \item{name = }{logit probability 6}
#'             \item{short.name = }{prob6}
#'             \item{initial = }{-1}
#'             \item{fixed = }{TRUE}
#'             \item{prior = }{gaussian}
#'             \item{param = }{-1 0.2}
#'             \item{to.theta = }{`function(x) log(x / (1 - x))`}
#'             \item{from.theta = }{`function(x) exp(x) / (1 + exp(x))`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta8'}{
#'              \describe{
#'             \item{hyperid = }{97008}
#'             \item{name = }{logit probability 7}
#'             \item{short.name = }{prob7}
#'             \item{initial = }{-1}
#'             \item{fixed = }{TRUE}
#'             \item{prior = }{gaussian}
#'             \item{param = }{-1 0.2}
#'             \item{to.theta = }{`function(x) log(x / (1 - x))`}
#'             \item{from.theta = }{`function(x) exp(x) / (1 + exp(x))`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta9'}{
#'              \describe{
#'             \item{hyperid = }{97009}
#'             \item{name = }{logit probability 8}
#'             \item{short.name = }{prob8}
#'             \item{initial = }{-1}
#'             \item{fixed = }{TRUE}
#'             \item{prior = }{gaussian}
#'             \item{param = }{-1 0.2}
#'             \item{to.theta = }{`function(x) log(x / (1 - x))`}
#'             \item{from.theta = }{`function(x) exp(x) / (1 + exp(x))`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta10'}{
#'              \describe{
#'             \item{hyperid = }{97010}
#'             \item{name = }{logit probability 9}
#'             \item{short.name = }{prob9}
#'             \item{initial = }{-1}
#'             \item{fixed = }{TRUE}
#'             \item{prior = }{gaussian}
#'             \item{param = }{-1 0.2}
#'             \item{to.theta = }{`function(x) log(x / (1 - x))`}
#'             \item{from.theta = }{`function(x) exp(x) / (1 + exp(x))`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta11'}{
#'              \describe{
#'             \item{hyperid = }{97011}
#'             \item{name = }{logit probability 10}
#'             \item{short.name = }{prob10}
#'             \item{initial = }{-1}
#'             \item{fixed = }{TRUE}
#'             \item{prior = }{gaussian}
#'             \item{param = }{-1 0.2}
#'             \item{to.theta = }{`function(x) log(x / (1 - x))`}
#'             \item{from.theta = }{`function(x) exp(x) / (1 + exp(x))`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'zeroinflatednbinomial1strata3'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'Zero inflated negBinomial, type 1, strata 3'}
#'               \item{status = }{'experimental'}
#'               \item{survival = }{'FALSE'}
#'               \item{discrete = }{'FALSE'}
#'               \item{link = }{'default log'}
#'               \item{pdf = }{'zeroinflated'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 11.
#'          \describe{
#'           \item{Hyperparameter 'theta1'}{
#'              \describe{
#'             \item{hyperid = }{98001}
#'             \item{name = }{logit probability}
#'             \item{short.name = }{prob}
#'             \item{initial = }{-1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{gaussian}
#'             \item{param = }{-1 0.2}
#'             \item{to.theta = }{`function(x) log(x / (1 - x))`}
#'             \item{from.theta = }{`function(x) exp(x) / (1 + exp(x))`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta2'}{
#'              \describe{
#'             \item{hyperid = }{98002}
#'             \item{name = }{log size 1}
#'             \item{short.name = }{size1}
#'             \item{initial = }{2.30258509299405}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{pc.mgamma}
#'             \item{param = }{7}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta3'}{
#'              \describe{
#'             \item{hyperid = }{98003}
#'             \item{name = }{log size 2}
#'             \item{short.name = }{size2}
#'             \item{initial = }{2.30258509299405}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{pc.mgamma}
#'             \item{param = }{7}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta4'}{
#'              \describe{
#'             \item{hyperid = }{98004}
#'             \item{name = }{log size 3}
#'             \item{short.name = }{size3}
#'             \item{initial = }{2.30258509299405}
#'             \item{fixed = }{TRUE}
#'             \item{prior = }{pc.mgamma}
#'             \item{param = }{7}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta5'}{
#'              \describe{
#'             \item{hyperid = }{98005}
#'             \item{name = }{log size 4}
#'             \item{short.name = }{size4}
#'             \item{initial = }{2.30258509299405}
#'             \item{fixed = }{TRUE}
#'             \item{prior = }{pc.mgamma}
#'             \item{param = }{7}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta6'}{
#'              \describe{
#'             \item{hyperid = }{98006}
#'             \item{name = }{log size 5}
#'             \item{short.name = }{size5}
#'             \item{initial = }{2.30258509299405}
#'             \item{fixed = }{TRUE}
#'             \item{prior = }{pc.mgamma}
#'             \item{param = }{7}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta7'}{
#'              \describe{
#'             \item{hyperid = }{98007}
#'             \item{name = }{log size 6}
#'             \item{short.name = }{size6}
#'             \item{initial = }{2.30258509299405}
#'             \item{fixed = }{TRUE}
#'             \item{prior = }{pc.mgamma}
#'             \item{param = }{7}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta8'}{
#'              \describe{
#'             \item{hyperid = }{98008}
#'             \item{name = }{log size 7}
#'             \item{short.name = }{size7}
#'             \item{initial = }{2.30258509299405}
#'             \item{fixed = }{TRUE}
#'             \item{prior = }{pc.mgamma}
#'             \item{param = }{7}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta9'}{
#'              \describe{
#'             \item{hyperid = }{98009}
#'             \item{name = }{log size 8}
#'             \item{short.name = }{size8}
#'             \item{initial = }{2.30258509299405}
#'             \item{fixed = }{TRUE}
#'             \item{prior = }{pc.mgamma}
#'             \item{param = }{7}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta10'}{
#'              \describe{
#'             \item{hyperid = }{98010}
#'             \item{name = }{log size 9}
#'             \item{short.name = }{size9}
#'             \item{initial = }{2.30258509299405}
#'             \item{fixed = }{TRUE}
#'             \item{prior = }{pc.mgamma}
#'             \item{param = }{7}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta11'}{
#'              \describe{
#'             \item{hyperid = }{98011}
#'             \item{name = }{log size 10}
#'             \item{short.name = }{size10}
#'             \item{initial = }{2.30258509299405}
#'             \item{fixed = }{TRUE}
#'             \item{prior = }{pc.mgamma}
#'             \item{param = }{7}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'zeroinflatednbinomial2'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'Zero inflated negBinomial, type 2'}
#'               \item{survival = }{'FALSE'}
#'               \item{discrete = }{'FALSE'}
#'               \item{link = }{'default log'}
#'               \item{pdf = }{'zeroinflated'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 2.
#'          \describe{
#'           \item{Hyperparameter 'theta1'}{
#'              \describe{
#'             \item{hyperid = }{99001}
#'             \item{name = }{log size}
#'             \item{short.name = }{size}
#'             \item{initial = }{2.30258509299405}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{pc.mgamma}
#'             \item{param = }{7}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta2'}{
#'              \describe{
#'             \item{hyperid = }{99002}
#'             \item{name = }{log alpha}
#'             \item{short.name = }{a}
#'             \item{initial = }{0.693147180559945}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{gaussian}
#'             \item{param = }{2 1}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 't'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'Student-t likelihood'}
#'               \item{survival = }{'FALSE'}
#'               \item{discrete = }{'FALSE'}
#'               \item{link = }{'default identity'}
#'               \item{pdf = }{'student-t'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 2.
#'          \describe{
#'           \item{Hyperparameter 'theta1'}{
#'              \describe{
#'             \item{hyperid = }{100001}
#'             \item{name = }{log precision}
#'             \item{short.name = }{prec}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{loggamma}
#'             \item{param = }{1 5e-05}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta2'}{
#'              \describe{
#'             \item{hyperid = }{100002}
#'             \item{name = }{log degrees of freedom}
#'             \item{short.name = }{dof}
#'             \item{initial = }{5}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{pc.dof}
#'             \item{param = }{15 0.5}
#'             \item{to.theta = }{`function(x) log(x - 2)`}
#'             \item{from.theta = }{`function(x) 2 + exp(x)`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'tstrata'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'A stratified version of the Student-t likelihood'}
#'               \item{survival = }{'FALSE'}
#'               \item{discrete = }{'FALSE'}
#'               \item{link = }{'default identity'}
#'               \item{pdf = }{'tstrata'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 11.
#'          \describe{
#'           \item{Hyperparameter 'theta1'}{
#'              \describe{
#'             \item{hyperid = }{101001}
#'             \item{name = }{log degrees of freedom}
#'             \item{short.name = }{dof}
#'             \item{initial = }{4}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{pc.dof}
#'             \item{param = }{15 0.5}
#'             \item{to.theta = }{`function(x) log(x - 5)`}
#'             \item{from.theta = }{`function(x) 5 + exp(x)`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta2'}{
#'              \describe{
#'             \item{hyperid = }{101002}
#'             \item{name = }{log precision1}
#'             \item{short.name = }{prec1}
#'             \item{initial = }{2}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{loggamma}
#'             \item{param = }{1 5e-05}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta3'}{
#'              \describe{
#'             \item{hyperid = }{101003}
#'             \item{name = }{log precision2}
#'             \item{short.name = }{prec2}
#'             \item{initial = }{2}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{loggamma}
#'             \item{param = }{1 5e-05}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta4'}{
#'              \describe{
#'             \item{hyperid = }{101004}
#'             \item{name = }{log precision3}
#'             \item{short.name = }{prec3}
#'             \item{initial = }{2}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{loggamma}
#'             \item{param = }{1 5e-05}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta5'}{
#'              \describe{
#'             \item{hyperid = }{101005}
#'             \item{name = }{log precision4}
#'             \item{short.name = }{prec4}
#'             \item{initial = }{2}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{loggamma}
#'             \item{param = }{1 5e-05}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta6'}{
#'              \describe{
#'             \item{hyperid = }{101006}
#'             \item{name = }{log precision5}
#'             \item{short.name = }{prec5}
#'             \item{initial = }{2}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{loggamma}
#'             \item{param = }{1 5e-05}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta7'}{
#'              \describe{
#'             \item{hyperid = }{101007}
#'             \item{name = }{log precision6}
#'             \item{short.name = }{prec6}
#'             \item{initial = }{2}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{loggamma}
#'             \item{param = }{1 5e-05}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta8'}{
#'              \describe{
#'             \item{hyperid = }{101008}
#'             \item{name = }{log precision7}
#'             \item{short.name = }{prec7}
#'             \item{initial = }{2}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{loggamma}
#'             \item{param = }{1 5e-05}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta9'}{
#'              \describe{
#'             \item{hyperid = }{101009}
#'             \item{name = }{log precision8}
#'             \item{short.name = }{prec8}
#'             \item{initial = }{2}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{loggamma}
#'             \item{param = }{1 5e-05}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta10'}{
#'              \describe{
#'             \item{hyperid = }{101010}
#'             \item{name = }{log precision9}
#'             \item{short.name = }{prec9}
#'             \item{initial = }{2}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{loggamma}
#'             \item{param = }{1 5e-05}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta11'}{
#'              \describe{
#'             \item{hyperid = }{101011}
#'             \item{name = }{log precision10}
#'             \item{short.name = }{prec10}
#'             \item{initial = }{2}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{loggamma}
#'             \item{param = }{1 5e-05}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'nmix'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'Binomial-Poisson mixture'}
#'               \item{status = }{'experimental'}
#'               \item{survival = }{'FALSE'}
#'               \item{discrete = }{'TRUE'}
#'               \item{link = }{'default logit loga probit'}
#'               \item{pdf = }{'nmix'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 15.
#'          \describe{
#'           \item{Hyperparameter 'theta1'}{
#'              \describe{
#'             \item{hyperid = }{101101}
#'             \item{name = }{beta1}
#'             \item{short.name = }{beta1}
#'             \item{initial = }{2.30258509299405}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 0.5}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta2'}{
#'              \describe{
#'             \item{hyperid = }{101102}
#'             \item{name = }{beta2}
#'             \item{short.name = }{beta2}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 1}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta3'}{
#'              \describe{
#'             \item{hyperid = }{101103}
#'             \item{name = }{beta3}
#'             \item{short.name = }{beta3}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 1}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta4'}{
#'              \describe{
#'             \item{hyperid = }{101104}
#'             \item{name = }{beta4}
#'             \item{short.name = }{beta4}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 1}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta5'}{
#'              \describe{
#'             \item{hyperid = }{101105}
#'             \item{name = }{beta5}
#'             \item{short.name = }{beta5}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 1}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta6'}{
#'              \describe{
#'             \item{hyperid = }{101106}
#'             \item{name = }{beta6}
#'             \item{short.name = }{beta6}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 1}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta7'}{
#'              \describe{
#'             \item{hyperid = }{101107}
#'             \item{name = }{beta7}
#'             \item{short.name = }{beta7}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 1}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta8'}{
#'              \describe{
#'             \item{hyperid = }{101108}
#'             \item{name = }{beta8}
#'             \item{short.name = }{beta8}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 1}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta9'}{
#'              \describe{
#'             \item{hyperid = }{101109}
#'             \item{name = }{beta9}
#'             \item{short.name = }{beta9}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 1}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta10'}{
#'              \describe{
#'             \item{hyperid = }{101110}
#'             \item{name = }{beta10}
#'             \item{short.name = }{beta10}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 1}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta11'}{
#'              \describe{
#'             \item{hyperid = }{101111}
#'             \item{name = }{beta11}
#'             \item{short.name = }{beta11}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 1}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta12'}{
#'              \describe{
#'             \item{hyperid = }{101112}
#'             \item{name = }{beta12}
#'             \item{short.name = }{beta12}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 1}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta13'}{
#'              \describe{
#'             \item{hyperid = }{101113}
#'             \item{name = }{beta13}
#'             \item{short.name = }{beta13}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 1}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta14'}{
#'              \describe{
#'             \item{hyperid = }{101114}
#'             \item{name = }{beta14}
#'             \item{short.name = }{beta14}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 1}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta15'}{
#'              \describe{
#'             \item{hyperid = }{101115}
#'             \item{name = }{beta15}
#'             \item{short.name = }{beta15}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 1}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'nmixnb'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'NegBinomial-Poisson mixture'}
#'               \item{status = }{'experimental'}
#'               \item{survival = }{'FALSE'}
#'               \item{discrete = }{'TRUE'}
#'               \item{link = }{'default logit loga probit'}
#'               \item{pdf = }{'nmixnb'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 16.
#'          \describe{
#'           \item{Hyperparameter 'theta1'}{
#'              \describe{
#'             \item{hyperid = }{101121}
#'             \item{name = }{beta1}
#'             \item{short.name = }{beta1}
#'             \item{initial = }{2.30258509299405}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 0.5}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta2'}{
#'              \describe{
#'             \item{hyperid = }{101122}
#'             \item{name = }{beta2}
#'             \item{short.name = }{beta2}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 1}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta3'}{
#'              \describe{
#'             \item{hyperid = }{101123}
#'             \item{name = }{beta3}
#'             \item{short.name = }{beta3}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 1}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta4'}{
#'              \describe{
#'             \item{hyperid = }{101124}
#'             \item{name = }{beta4}
#'             \item{short.name = }{beta4}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 1}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta5'}{
#'              \describe{
#'             \item{hyperid = }{101125}
#'             \item{name = }{beta5}
#'             \item{short.name = }{beta5}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 1}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta6'}{
#'              \describe{
#'             \item{hyperid = }{101126}
#'             \item{name = }{beta6}
#'             \item{short.name = }{beta6}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 1}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta7'}{
#'              \describe{
#'             \item{hyperid = }{101127}
#'             \item{name = }{beta7}
#'             \item{short.name = }{beta7}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 1}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta8'}{
#'              \describe{
#'             \item{hyperid = }{101128}
#'             \item{name = }{beta8}
#'             \item{short.name = }{beta8}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 1}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta9'}{
#'              \describe{
#'             \item{hyperid = }{101129}
#'             \item{name = }{beta9}
#'             \item{short.name = }{beta9}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 1}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta10'}{
#'              \describe{
#'             \item{hyperid = }{101130}
#'             \item{name = }{beta10}
#'             \item{short.name = }{beta10}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 1}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta11'}{
#'              \describe{
#'             \item{hyperid = }{101131}
#'             \item{name = }{beta11}
#'             \item{short.name = }{beta11}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 1}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta12'}{
#'              \describe{
#'             \item{hyperid = }{101132}
#'             \item{name = }{beta12}
#'             \item{short.name = }{beta12}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 1}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta13'}{
#'              \describe{
#'             \item{hyperid = }{101133}
#'             \item{name = }{beta13}
#'             \item{short.name = }{beta13}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 1}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta14'}{
#'              \describe{
#'             \item{hyperid = }{101134}
#'             \item{name = }{beta14}
#'             \item{short.name = }{beta14}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 1}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta15'}{
#'              \describe{
#'             \item{hyperid = }{101135}
#'             \item{name = }{beta15}
#'             \item{short.name = }{beta15}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 1}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta16'}{
#'              \describe{
#'             \item{hyperid = }{101136}
#'             \item{name = }{overdispersion}
#'             \item{short.name = }{overdispersion}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{pc.gamma}
#'             \item{param = }{7}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'gp'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'Generalized Pareto likelihood'}
#'               \item{status = }{'experimental'}
#'               \item{survival = }{'FALSE'}
#'               \item{discrete = }{'TRUE'}
#'               \item{link = }{'default quantile'}
#'               \item{pdf = }{'genPareto'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 1.
#'          \describe{
#'           \item{Hyperparameter 'theta'}{
#'              \describe{
#'             \item{hyperid = }{101201}
#'             \item{name = }{tail}
#'             \item{short.name = }{xi}
#'             \item{initial = }{-4}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{pc.gevtail}
#'             \item{param = }{7 0 0.5}
#'             \item{to.theta = }{`function(x, interval = c(REPLACE.ME.low, REPLACE.ME.high)) log(-(interval[1] - x) / (interval[2] - x))`}
#'             \item{from.theta = }{`function(x, interval = c(REPLACE.ME.low, REPLACE.ME.high)) interval[1] + (interval[2] - interval[1]) * exp(x) / (1.0 + exp(x))`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'dgp'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'Discrete generalized Pareto likelihood'}
#'               \item{status = }{'experimental'}
#'               \item{survival = }{'FALSE'}
#'               \item{discrete = }{'TRUE'}
#'               \item{link = }{'default quantile'}
#'               \item{pdf = }{'dgp'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 1.
#'          \describe{
#'           \item{Hyperparameter 'theta'}{
#'              \describe{
#'             \item{hyperid = }{101301}
#'             \item{name = }{tail}
#'             \item{short.name = }{xi}
#'             \item{initial = }{2}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{pc.gevtail}
#'             \item{param = }{7 0 0.5}
#'             \item{to.theta = }{`function(x, interval = c(REPLACE.ME.low, REPLACE.ME.high)) log(-(interval[1] - x) / (interval[2] - x))`}
#'             \item{from.theta = }{`function(x, interval = c(REPLACE.ME.low, REPLACE.ME.high)) interval[1] + (interval[2] - interval[1]) * exp(x) / (1.0 + exp(x))`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'logperiodogram'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'Likelihood for the log-periodogram'}
#'               \item{survival = }{'FALSE'}
#'               \item{discrete = }{'FALSE'}
#'               \item{link = }{'default identity'}
#'               \item{pdf = }{'NA'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 0.
#'        }
#'       \item{Model 'tweedie'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'Tweedie distribution'}
#'               \item{status = }{'experimental'}
#'               \item{survival = }{'FALSE'}
#'               \item{discrete = }{'FALSE'}
#'               \item{link = }{'default log'}
#'               \item{pdf = }{'tweedie'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 2.
#'          \describe{
#'           \item{Hyperparameter 'theta1'}{
#'              \describe{
#'             \item{hyperid = }{102101}
#'             \item{name = }{p}
#'             \item{short.name = }{p}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 100}
#'             \item{to.theta = }{`function(x, interval = c(1.0, 2.0)) log(-(interval[1] - x) / (interval[2] - x))`}
#'             \item{from.theta = }{`function(x, interval = c(1.0, 2.0)) interval[1] + (interval[2] - interval[1]) * exp(x) / (1.0 + exp(x))`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta2'}{
#'              \describe{
#'             \item{hyperid = }{102201}
#'             \item{name = }{dispersion}
#'             \item{short.name = }{phi}
#'             \item{initial = }{-4}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{loggamma}
#'             \item{param = }{100 100}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'fmri'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'fmri distribution (special nc-chi)'}
#'               \item{status = }{'experimental'}
#'               \item{survival = }{'FALSE'}
#'               \item{discrete = }{'FALSE'}
#'               \item{link = }{'default log'}
#'               \item{pdf = }{'fmri'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 2.
#'          \describe{
#'           \item{Hyperparameter 'theta1'}{
#'              \describe{
#'             \item{hyperid = }{103101}
#'             \item{name = }{precision}
#'             \item{short.name = }{prec}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{loggamma}
#'             \item{param = }{10 10}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta2'}{
#'              \describe{
#'             \item{hyperid = }{103202}
#'             \item{name = }{dof}
#'             \item{short.name = }{df}
#'             \item{initial = }{4}
#'             \item{fixed = }{TRUE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 1}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'fmrisurv'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'fmri distribution (special nc-chi)'}
#'               \item{status = }{'experimental'}
#'               \item{survival = }{'TRUE'}
#'               \item{discrete = }{'FALSE'}
#'               \item{link = }{'default log'}
#'               \item{pdf = }{'fmri'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 2.
#'          \describe{
#'           \item{Hyperparameter 'theta1'}{
#'              \describe{
#'             \item{hyperid = }{104101}
#'             \item{name = }{precision}
#'             \item{short.name = }{prec}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{loggamma}
#'             \item{param = }{10 10}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta2'}{
#'              \describe{
#'             \item{hyperid = }{104201}
#'             \item{name = }{dof}
#'             \item{short.name = }{df}
#'             \item{initial = }{4}
#'             \item{fixed = }{TRUE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 1}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'gompertz'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'gompertz distribution'}
#'               \item{status = }{'experimental'}
#'               \item{survival = }{'FALSE'}
#'               \item{discrete = }{'FALSE'}
#'               \item{link = }{'default log neglog'}
#'               \item{pdf = }{'gompertz'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 1.
#'          \describe{
#'           \item{Hyperparameter 'theta'}{
#'              \describe{
#'             \item{hyperid = }{105101}
#'             \item{name = }{shape}
#'             \item{short.name = }{alpha}
#'             \item{initial = }{-1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 1}
#'             \item{to.theta = }{`function(x, sc = 0.1) log(x) / sc`}
#'             \item{from.theta = }{`function(x, sc = 0.1) exp(sc * x)`}
#'             }
#'           }
#'          }
#'        }
#'       \item{Model 'gompertzsurv'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'gompertz distribution'}
#'               \item{status = }{'experimental'}
#'               \item{survival = }{'TRUE'}
#'               \item{discrete = }{'FALSE'}
#'               \item{link = }{'default log neglog'}
#'               \item{pdf = }{'gompertz'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 11.
#'          \describe{
#'           \item{Hyperparameter 'theta1'}{
#'              \describe{
#'             \item{hyperid = }{106101}
#'             \item{name = }{shape}
#'             \item{short.name = }{alpha}
#'             \item{initial = }{-10}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 1}
#'             \item{to.theta = }{`function(x, sc = 0.1) log(x) / sc`}
#'             \item{from.theta = }{`function(x, sc = 0.1) exp(sc * x)`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta2'}{
#'              \describe{
#'             \item{hyperid = }{106102}
#'             \item{name = }{beta1}
#'             \item{short.name = }{beta1}
#'             \item{initial = }{-5}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{-4 100}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta3'}{
#'              \describe{
#'             \item{hyperid = }{106103}
#'             \item{name = }{beta2}
#'             \item{short.name = }{beta2}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 100}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta4'}{
#'              \describe{
#'             \item{hyperid = }{106104}
#'             \item{name = }{beta3}
#'             \item{short.name = }{beta3}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 100}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta5'}{
#'              \describe{
#'             \item{hyperid = }{106105}
#'             \item{name = }{beta4}
#'             \item{short.name = }{beta4}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 100}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta6'}{
#'              \describe{
#'             \item{hyperid = }{106106}
#'             \item{name = }{beta5}
#'             \item{short.name = }{beta5}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 100}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta7'}{
#'              \describe{
#'             \item{hyperid = }{106107}
#'             \item{name = }{beta6}
#'             \item{short.name = }{beta6}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 100}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta8'}{
#'              \describe{
#'             \item{hyperid = }{106108}
#'             \item{name = }{beta7}
#'             \item{short.name = }{beta7}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 100}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta9'}{
#'              \describe{
#'             \item{hyperid = }{106109}
#'             \item{name = }{beta8}
#'             \item{short.name = }{beta8}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 100}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta10'}{
#'              \describe{
#'             \item{hyperid = }{106110}
#'             \item{name = }{beta9}
#'             \item{short.name = }{beta9}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 100}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta11'}{
#'              \describe{
#'             \item{hyperid = }{106111}
#'             \item{name = }{beta10}
#'             \item{short.name = }{beta10}
#'             \item{initial = }{0}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{0 100}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'          }
#'        }
#'      }
#' @section 'prior':
#'
#'    Valid models in this section are:
#'      \describe{
#'       \item{Model 'normal'.}{
#'          Number of parameters in the prior = 2
#'        }
#'       \item{Model 'gaussian'.}{
#'          Number of parameters in the prior = 2
#'        }
#'       \item{Model 'linksnintercept'.}{
#'          Number of parameters in the prior = 2
#'        }
#'       \item{Model 'wishart1d'.}{
#'          Number of parameters in the prior = 2
#'        }
#'       \item{Model 'wishart2d'.}{
#'          Number of parameters in the prior = 4
#'        }
#'       \item{Model 'wishart3d'.}{
#'          Number of parameters in the prior = 7
#'        }
#'       \item{Model 'wishart4d'.}{
#'          Number of parameters in the prior = 11
#'        }
#'       \item{Model 'wishart5d'.}{
#'          Number of parameters in the prior = 16
#'        }
#'       \item{Model 'loggamma'.}{
#'          Number of parameters in the prior = 2
#'        }
#'       \item{Model 'gamma'.}{
#'          Number of parameters in the prior = 2
#'        }
#'       \item{Model 'minuslogsqrtruncnormal'.}{
#'          Number of parameters in the prior = 2
#'        }
#'       \item{Model 'logtnormal'.}{
#'          Number of parameters in the prior = 2
#'        }
#'       \item{Model 'logtgaussian'.}{
#'          Number of parameters in the prior = 2
#'        }
#'       \item{Model 'flat'.}{
#'          Number of parameters in the prior = 0
#'        }
#'       \item{Model 'logflat'.}{
#'          Number of parameters in the prior = 0
#'        }
#'       \item{Model 'logiflat'.}{
#'          Number of parameters in the prior = 0
#'        }
#'       \item{Model 'mvnorm'.}{
#'          Number of parameters in the prior = -1
#'        }
#'       \item{Model 'pc.alphaw'.}{
#'          Number of parameters in the prior = 1
#'        }
#'       \item{Model 'pc.ar'.}{
#'          Number of parameters in the prior = 1
#'        }
#'       \item{Model 'dirichlet'.}{
#'          Number of parameters in the prior = 1
#'        }
#'       \item{Model 'none'.}{
#'          Number of parameters in the prior = 0
#'        }
#'       \item{Model 'invalid'.}{
#'          Number of parameters in the prior = 0
#'        }
#'       \item{Model 'betacorrelation'.}{
#'          Number of parameters in the prior = 2
#'        }
#'       \item{Model 'logitbeta'.}{
#'          Number of parameters in the prior = 2
#'        }
#'       \item{Model 'pc.prec'.}{
#'          Number of parameters in the prior = 2
#'        }
#'       \item{Model 'pc.dof'.}{
#'          Number of parameters in the prior = 2
#'        }
#'       \item{Model 'pc.cor0'.}{
#'          Number of parameters in the prior = 2
#'        }
#'       \item{Model 'pc.cor1'.}{
#'          Number of parameters in the prior = 2
#'        }
#'       \item{Model 'pc.fgnh'.}{
#'          Number of parameters in the prior = 2
#'        }
#'       \item{Model 'pc.spde.GA'.}{
#'          Number of parameters in the prior = 4
#'        }
#'       \item{Model 'pc.matern'.}{
#'          Number of parameters in the prior = 3
#'        }
#'       \item{Model 'pc.range'.}{
#'          Number of parameters in the prior = 2
#'        }
#'       \item{Model 'pc.sn'.}{
#'          Number of parameters in the prior = 1
#'        }
#'       \item{Model 'pc.gamma'.}{
#'          Number of parameters in the prior = 1
#'        }
#'       \item{Model 'pc.mgamma'.}{
#'          Number of parameters in the prior = 1
#'        }
#'       \item{Model 'pc.gammacount'.}{
#'          Number of parameters in the prior = 1
#'        }
#'       \item{Model 'pc.gevtail'.}{
#'          Number of parameters in the prior = 3
#'        }
#'       \item{Model 'pc'.}{
#'          Number of parameters in the prior = 2
#'        }
#'       \item{Model 'ref.ar'.}{
#'          Number of parameters in the prior = 0
#'        }
#'       \item{Model 'pom'.}{
#'          Number of parameters in the prior = 0
#'        }
#'       \item{Model 'jeffreystdf'.}{
#'          Number of parameters in the prior = 0
#'        }
#'       \item{Model 'wishartkd'.}{
#'          Number of parameters in the prior = 211
#'        }
#'       \item{Model 'expression:'.}{
#'          Number of parameters in the prior = -1
#'        }
#'       \item{Model 'table:'.}{
#'          Number of parameters in the prior = -1
#'        }
#'      }
#' @section 'wrapper':
#'
#'    Valid models in this section are:
#'      \describe{
#'       \item{Model 'joint'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{doc = }{'(experimental)'}
#'               \item{constr = }{'FALSE'}
#'               \item{nrow.ncol = }{'FALSE'}
#'               \item{augmented = }{'FALSE'}
#'               \item{aug.factor = }{'1'}
#'               \item{aug.constr = }{'NULL'}
#'               \item{n.div.by = }{'NULL'}
#'               \item{n.required = }{'FALSE'}
#'               \item{set.default.values = }{'FALSE'}
#'               \item{pdf = }{'NA'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 1.
#'          \describe{
#'           \item{Hyperparameter 'theta'}{
#'              \describe{
#'             \item{hyperid = }{102001}
#'             \item{name = }{log precision}
#'             \item{short.name = }{prec}
#'             \item{initial = }{0}
#'             \item{fixed = }{TRUE}
#'             \item{prior = }{loggamma}
#'             \item{param = }{1 5e-05}
#'             \item{to.theta = }{`function(x) log(x)`}
#'             \item{from.theta = }{`function(x) exp(x)`}
#'             }
#'           }
#'          }
#'        }
#'      }
#' @section 'lp.scale':
#'
#'    Valid models in this section are:
#'      \describe{
#'       \item{Model 'lp.scale'.}{
#'          \describe{
#'           \item{Properties:}{
#'             \describe{
#'               \item{pdf = }{'lp.scale'}
#'              }
#'            }
#'          }
#'         Number of hyperparmeters is 100.
#'          \describe{
#'           \item{Hyperparameter 'theta1'}{
#'              \describe{
#'             \item{hyperid = }{103001}
#'             \item{name = }{beta1}
#'             \item{short.name = }{b1}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{1 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta2'}{
#'              \describe{
#'             \item{hyperid = }{103002}
#'             \item{name = }{beta2}
#'             \item{short.name = }{b2}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{1 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta3'}{
#'              \describe{
#'             \item{hyperid = }{103003}
#'             \item{name = }{beta3}
#'             \item{short.name = }{b3}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{1 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta4'}{
#'              \describe{
#'             \item{hyperid = }{103004}
#'             \item{name = }{beta4}
#'             \item{short.name = }{b4}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{1 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta5'}{
#'              \describe{
#'             \item{hyperid = }{103005}
#'             \item{name = }{beta5}
#'             \item{short.name = }{b5}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{1 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta6'}{
#'              \describe{
#'             \item{hyperid = }{103006}
#'             \item{name = }{beta6}
#'             \item{short.name = }{b6}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{1 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta7'}{
#'              \describe{
#'             \item{hyperid = }{103007}
#'             \item{name = }{beta7}
#'             \item{short.name = }{b7}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{1 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta8'}{
#'              \describe{
#'             \item{hyperid = }{103008}
#'             \item{name = }{beta8}
#'             \item{short.name = }{b8}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{1 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta9'}{
#'              \describe{
#'             \item{hyperid = }{103009}
#'             \item{name = }{beta9}
#'             \item{short.name = }{b9}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{1 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta10'}{
#'              \describe{
#'             \item{hyperid = }{103010}
#'             \item{name = }{beta10}
#'             \item{short.name = }{b10}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{1 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta11'}{
#'              \describe{
#'             \item{hyperid = }{103011}
#'             \item{name = }{beta11}
#'             \item{short.name = }{b11}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{1 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta12'}{
#'              \describe{
#'             \item{hyperid = }{103012}
#'             \item{name = }{beta12}
#'             \item{short.name = }{b12}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{1 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta13'}{
#'              \describe{
#'             \item{hyperid = }{103013}
#'             \item{name = }{beta13}
#'             \item{short.name = }{b13}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{1 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta14'}{
#'              \describe{
#'             \item{hyperid = }{103014}
#'             \item{name = }{beta14}
#'             \item{short.name = }{b14}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{1 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta15'}{
#'              \describe{
#'             \item{hyperid = }{103015}
#'             \item{name = }{beta15}
#'             \item{short.name = }{b15}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{1 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta16'}{
#'              \describe{
#'             \item{hyperid = }{103016}
#'             \item{name = }{beta16}
#'             \item{short.name = }{b16}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{1 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta17'}{
#'              \describe{
#'             \item{hyperid = }{103017}
#'             \item{name = }{beta17}
#'             \item{short.name = }{b17}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{1 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta18'}{
#'              \describe{
#'             \item{hyperid = }{103018}
#'             \item{name = }{beta18}
#'             \item{short.name = }{b18}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{1 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta19'}{
#'              \describe{
#'             \item{hyperid = }{103019}
#'             \item{name = }{beta19}
#'             \item{short.name = }{b19}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{1 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta20'}{
#'              \describe{
#'             \item{hyperid = }{103020}
#'             \item{name = }{beta20}
#'             \item{short.name = }{b20}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{1 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta21'}{
#'              \describe{
#'             \item{hyperid = }{103021}
#'             \item{name = }{beta21}
#'             \item{short.name = }{b21}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{1 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta22'}{
#'              \describe{
#'             \item{hyperid = }{103022}
#'             \item{name = }{beta22}
#'             \item{short.name = }{b22}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{1 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta23'}{
#'              \describe{
#'             \item{hyperid = }{103023}
#'             \item{name = }{beta23}
#'             \item{short.name = }{b23}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{1 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta24'}{
#'              \describe{
#'             \item{hyperid = }{103024}
#'             \item{name = }{beta24}
#'             \item{short.name = }{b24}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{1 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta25'}{
#'              \describe{
#'             \item{hyperid = }{103025}
#'             \item{name = }{beta25}
#'             \item{short.name = }{b25}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{1 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta26'}{
#'              \describe{
#'             \item{hyperid = }{103026}
#'             \item{name = }{beta26}
#'             \item{short.name = }{b26}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{1 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta27'}{
#'              \describe{
#'             \item{hyperid = }{103027}
#'             \item{name = }{beta27}
#'             \item{short.name = }{b27}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{1 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta28'}{
#'              \describe{
#'             \item{hyperid = }{103028}
#'             \item{name = }{beta28}
#'             \item{short.name = }{b28}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{1 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta29'}{
#'              \describe{
#'             \item{hyperid = }{103029}
#'             \item{name = }{beta29}
#'             \item{short.name = }{b29}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{1 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta30'}{
#'              \describe{
#'             \item{hyperid = }{103030}
#'             \item{name = }{beta30}
#'             \item{short.name = }{b30}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{1 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta31'}{
#'              \describe{
#'             \item{hyperid = }{103031}
#'             \item{name = }{beta31}
#'             \item{short.name = }{b31}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{1 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta32'}{
#'              \describe{
#'             \item{hyperid = }{103032}
#'             \item{name = }{beta32}
#'             \item{short.name = }{b32}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{1 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta33'}{
#'              \describe{
#'             \item{hyperid = }{103033}
#'             \item{name = }{beta33}
#'             \item{short.name = }{b33}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{1 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta34'}{
#'              \describe{
#'             \item{hyperid = }{103034}
#'             \item{name = }{beta34}
#'             \item{short.name = }{b34}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{1 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta35'}{
#'              \describe{
#'             \item{hyperid = }{103035}
#'             \item{name = }{beta35}
#'             \item{short.name = }{b35}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{1 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta36'}{
#'              \describe{
#'             \item{hyperid = }{103036}
#'             \item{name = }{beta36}
#'             \item{short.name = }{b36}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{1 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta37'}{
#'              \describe{
#'             \item{hyperid = }{103037}
#'             \item{name = }{beta37}
#'             \item{short.name = }{b37}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{1 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta38'}{
#'              \describe{
#'             \item{hyperid = }{103038}
#'             \item{name = }{beta38}
#'             \item{short.name = }{b38}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{1 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta39'}{
#'              \describe{
#'             \item{hyperid = }{103039}
#'             \item{name = }{beta39}
#'             \item{short.name = }{b39}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{1 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta40'}{
#'              \describe{
#'             \item{hyperid = }{103040}
#'             \item{name = }{beta40}
#'             \item{short.name = }{b40}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{1 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta41'}{
#'              \describe{
#'             \item{hyperid = }{103041}
#'             \item{name = }{beta41}
#'             \item{short.name = }{b41}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{1 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta42'}{
#'              \describe{
#'             \item{hyperid = }{103042}
#'             \item{name = }{beta42}
#'             \item{short.name = }{b42}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{1 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta43'}{
#'              \describe{
#'             \item{hyperid = }{103043}
#'             \item{name = }{beta43}
#'             \item{short.name = }{b43}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{1 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta44'}{
#'              \describe{
#'             \item{hyperid = }{103044}
#'             \item{name = }{beta44}
#'             \item{short.name = }{b44}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{1 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta45'}{
#'              \describe{
#'             \item{hyperid = }{103045}
#'             \item{name = }{beta45}
#'             \item{short.name = }{b45}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{1 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta46'}{
#'              \describe{
#'             \item{hyperid = }{103046}
#'             \item{name = }{beta46}
#'             \item{short.name = }{b46}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{1 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta47'}{
#'              \describe{
#'             \item{hyperid = }{103047}
#'             \item{name = }{beta47}
#'             \item{short.name = }{b47}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{1 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta48'}{
#'              \describe{
#'             \item{hyperid = }{103048}
#'             \item{name = }{beta48}
#'             \item{short.name = }{b48}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{1 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta49'}{
#'              \describe{
#'             \item{hyperid = }{103049}
#'             \item{name = }{beta49}
#'             \item{short.name = }{b49}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{1 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta50'}{
#'              \describe{
#'             \item{hyperid = }{103050}
#'             \item{name = }{beta50}
#'             \item{short.name = }{b50}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{1 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta51'}{
#'              \describe{
#'             \item{hyperid = }{103051}
#'             \item{name = }{beta51}
#'             \item{short.name = }{b51}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{1 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta52'}{
#'              \describe{
#'             \item{hyperid = }{103052}
#'             \item{name = }{beta52}
#'             \item{short.name = }{b52}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{1 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta53'}{
#'              \describe{
#'             \item{hyperid = }{103053}
#'             \item{name = }{beta53}
#'             \item{short.name = }{b53}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{1 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta54'}{
#'              \describe{
#'             \item{hyperid = }{103054}
#'             \item{name = }{beta54}
#'             \item{short.name = }{b54}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{1 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta55'}{
#'              \describe{
#'             \item{hyperid = }{103055}
#'             \item{name = }{beta55}
#'             \item{short.name = }{b55}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{1 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta56'}{
#'              \describe{
#'             \item{hyperid = }{103056}
#'             \item{name = }{beta56}
#'             \item{short.name = }{b56}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{1 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta57'}{
#'              \describe{
#'             \item{hyperid = }{103057}
#'             \item{name = }{beta57}
#'             \item{short.name = }{b57}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{1 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta58'}{
#'              \describe{
#'             \item{hyperid = }{103058}
#'             \item{name = }{beta58}
#'             \item{short.name = }{b58}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{1 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta59'}{
#'              \describe{
#'             \item{hyperid = }{103059}
#'             \item{name = }{beta59}
#'             \item{short.name = }{b59}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{1 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta60'}{
#'              \describe{
#'             \item{hyperid = }{103060}
#'             \item{name = }{beta60}
#'             \item{short.name = }{b60}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{1 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta61'}{
#'              \describe{
#'             \item{hyperid = }{103061}
#'             \item{name = }{beta61}
#'             \item{short.name = }{b61}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{1 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta62'}{
#'              \describe{
#'             \item{hyperid = }{103062}
#'             \item{name = }{beta62}
#'             \item{short.name = }{b62}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{1 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta63'}{
#'              \describe{
#'             \item{hyperid = }{103063}
#'             \item{name = }{beta63}
#'             \item{short.name = }{b63}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{1 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta64'}{
#'              \describe{
#'             \item{hyperid = }{103064}
#'             \item{name = }{beta64}
#'             \item{short.name = }{b64}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{1 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta65'}{
#'              \describe{
#'             \item{hyperid = }{103065}
#'             \item{name = }{beta65}
#'             \item{short.name = }{b65}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{1 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta66'}{
#'              \describe{
#'             \item{hyperid = }{103066}
#'             \item{name = }{beta66}
#'             \item{short.name = }{b66}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{1 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta67'}{
#'              \describe{
#'             \item{hyperid = }{103067}
#'             \item{name = }{beta67}
#'             \item{short.name = }{b67}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{1 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta68'}{
#'              \describe{
#'             \item{hyperid = }{103068}
#'             \item{name = }{beta68}
#'             \item{short.name = }{b68}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{1 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta69'}{
#'              \describe{
#'             \item{hyperid = }{103069}
#'             \item{name = }{beta69}
#'             \item{short.name = }{b69}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{1 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta70'}{
#'              \describe{
#'             \item{hyperid = }{103070}
#'             \item{name = }{beta70}
#'             \item{short.name = }{b70}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{1 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta71'}{
#'              \describe{
#'             \item{hyperid = }{103071}
#'             \item{name = }{beta71}
#'             \item{short.name = }{b71}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{1 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta72'}{
#'              \describe{
#'             \item{hyperid = }{103072}
#'             \item{name = }{beta72}
#'             \item{short.name = }{b72}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{1 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta73'}{
#'              \describe{
#'             \item{hyperid = }{103073}
#'             \item{name = }{beta73}
#'             \item{short.name = }{b73}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{1 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta74'}{
#'              \describe{
#'             \item{hyperid = }{103074}
#'             \item{name = }{beta74}
#'             \item{short.name = }{b74}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{1 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta75'}{
#'              \describe{
#'             \item{hyperid = }{103075}
#'             \item{name = }{beta75}
#'             \item{short.name = }{b75}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{1 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta76'}{
#'              \describe{
#'             \item{hyperid = }{103076}
#'             \item{name = }{beta76}
#'             \item{short.name = }{b76}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{1 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta77'}{
#'              \describe{
#'             \item{hyperid = }{103077}
#'             \item{name = }{beta77}
#'             \item{short.name = }{b77}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{1 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta78'}{
#'              \describe{
#'             \item{hyperid = }{103078}
#'             \item{name = }{beta78}
#'             \item{short.name = }{b78}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{1 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta79'}{
#'              \describe{
#'             \item{hyperid = }{103079}
#'             \item{name = }{beta79}
#'             \item{short.name = }{b79}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{1 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta80'}{
#'              \describe{
#'             \item{hyperid = }{103080}
#'             \item{name = }{beta80}
#'             \item{short.name = }{b80}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{1 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta81'}{
#'              \describe{
#'             \item{hyperid = }{103081}
#'             \item{name = }{beta81}
#'             \item{short.name = }{b81}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{1 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta82'}{
#'              \describe{
#'             \item{hyperid = }{103082}
#'             \item{name = }{beta82}
#'             \item{short.name = }{b82}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{1 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta83'}{
#'              \describe{
#'             \item{hyperid = }{103083}
#'             \item{name = }{beta83}
#'             \item{short.name = }{b83}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{1 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta84'}{
#'              \describe{
#'             \item{hyperid = }{103084}
#'             \item{name = }{beta84}
#'             \item{short.name = }{b84}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{1 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta85'}{
#'              \describe{
#'             \item{hyperid = }{103085}
#'             \item{name = }{beta85}
#'             \item{short.name = }{b85}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{1 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta86'}{
#'              \describe{
#'             \item{hyperid = }{103086}
#'             \item{name = }{beta86}
#'             \item{short.name = }{b86}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{1 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta87'}{
#'              \describe{
#'             \item{hyperid = }{103087}
#'             \item{name = }{beta87}
#'             \item{short.name = }{b87}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{1 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta88'}{
#'              \describe{
#'             \item{hyperid = }{103088}
#'             \item{name = }{beta88}
#'             \item{short.name = }{b88}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{1 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta89'}{
#'              \describe{
#'             \item{hyperid = }{103089}
#'             \item{name = }{beta89}
#'             \item{short.name = }{b89}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{1 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta90'}{
#'              \describe{
#'             \item{hyperid = }{103090}
#'             \item{name = }{beta90}
#'             \item{short.name = }{b90}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{1 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta91'}{
#'              \describe{
#'             \item{hyperid = }{103091}
#'             \item{name = }{beta91}
#'             \item{short.name = }{b91}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{1 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta92'}{
#'              \describe{
#'             \item{hyperid = }{103092}
#'             \item{name = }{beta92}
#'             \item{short.name = }{b92}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{1 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta93'}{
#'              \describe{
#'             \item{hyperid = }{103093}
#'             \item{name = }{beta93}
#'             \item{short.name = }{b93}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{1 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta94'}{
#'              \describe{
#'             \item{hyperid = }{103094}
#'             \item{name = }{beta94}
#'             \item{short.name = }{b94}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{1 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta95'}{
#'              \describe{
#'             \item{hyperid = }{103095}
#'             \item{name = }{beta95}
#'             \item{short.name = }{b95}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{1 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta96'}{
#'              \describe{
#'             \item{hyperid = }{103096}
#'             \item{name = }{beta96}
#'             \item{short.name = }{b96}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{1 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta97'}{
#'              \describe{
#'             \item{hyperid = }{103097}
#'             \item{name = }{beta97}
#'             \item{short.name = }{b97}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{1 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta98'}{
#'              \describe{
#'             \item{hyperid = }{103098}
#'             \item{name = }{beta98}
#'             \item{short.name = }{b98}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{1 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta99'}{
#'              \describe{
#'             \item{hyperid = }{103099}
#'             \item{name = }{beta99}
#'             \item{short.name = }{b99}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{1 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'           \item{Hyperparameter 'theta100'}{
#'              \describe{
#'             \item{hyperid = }{103100}
#'             \item{name = }{beta100}
#'             \item{short.name = }{b100}
#'             \item{initial = }{1}
#'             \item{fixed = }{FALSE}
#'             \item{prior = }{normal}
#'             \item{param = }{1 10}
#'             \item{to.theta = }{`function(x) x`}
#'             \item{from.theta = }{`function(x) x`}
#'             }
#'           }
#'          }
#'        }
#'      }
#' @examples
#' ## How to set hyperparameters to pass as the argument 'hyper'. This
#' ## format is compatible with the old style (using 'initial', 'fixed',
#' ## 'prior', 'param'), but the new style using 'hyper' takes precedence
#' ## over the old style. The two styles can also be mixed. The old style
#' ## might be removed from the code in the future...
#'
#' ## Only a subset need to be given
#' hyper <- list(theta = list(initial = 2))
#' ## The `name' can be used instead of 'theta', or 'theta1', 'theta2',...
#' hyper <- list(precision = list(initial = 2))
#' hyper <- list(precision = list(prior = "flat", param = numeric(0)))
#' hyper <- list(theta2 = list(initial = 3), theta1 = list(prior = "gaussian"))
#' ## The 'short.name' can be used instead of 'name'
#' hyper <- list(rho = list(param = c(0, 1)))
NULL
