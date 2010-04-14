\name{inla}
\alias{inla}
\alias{xinla}

\title{Bayesian analysis of structured additive models}

\description{\code{inla} performs a full Bayesian analysis of  additive models using Integrated Nested Laplace approximation
 }

\usage{
inla =
    function (formula,
              family = "gaussian", 
              data = data.frame(),
              quantiles=c(0.025,0.975),
              E = NULL,
              offset=NULL,
              scale = NULL,
              Ntrials = NULL,
              verbose = FALSE,
              control.compute = list(),
              control.predictor = list(),
              control.data = list(),
              control.inla = list(),
              control.results = list(),
              control.fixed = list(),
              control.mode = list(),
              control.expert = list(),
              only.hyperparam = FALSE,
              inla.call = inla.getOption("inla.call"),
              inla.arg = inla.getOption("inla.arg"),
              num.threads = inla.getOption("num.threads"),
              keep = inla.getOption("keep"),
              working.directory = inla.getOption("working.directory"),
              silent = inla.getOption("silent"),
              disable.check = inla.getOption("disable.check"),
              debug = inla.getOption("debug"),
              user.hook = NULL,
              user.hook.arg = NULL
              )
}

\arguments{

  \item{formula}{A \code{inla} formula (see also
    \code{\link{inla.models}}), like

    \code{y ~ 1 + z + f(ind, model="iid")} + f(ind2, weights, model="ar1")

    This is much like the formula for a \code{glm} except that smooth or
    spatial terms can be added to the right hand side of the formula.
    See \code{\link{f}} for full details. Each smooth or spatial term
    specified through \code{f} should correspond to separate column of
    the data frame \code{data}.

    The response variable, \code{y} in this example, is usually a
    uni-variate response variable. The exception is for survival-models
    (see variable \link{family}), where the response is a data.frame or
    list, with elements

    \code{y = list(truncation=, event=, lower=, upper=, time=)}.

    Here, \code{truncation} is the truncation (like conditioning on
    \code{time} > \code{truncation}), \code{event} is 1 if the failure
    time is observed, 0 for right censoring, 2 for left censoring and 3
    for interval censoring. \code{lower} and \code{upper} are the lower
    and upper limits of the failure \code{time}. Note that none of the
    events use all variables, and those who are not used should be set
    to zero.
}
  
  \item{family}{A string indicating the likelihood family. The default
    is \code{gaussian} with identity link. Other possible choices are:
    \code{poisson}, \code{binomial}, \code{zeroinflated_binomial_0},
    \code{zeroinflated_binomial_1}, \code{T}, \code{stochvol},
    \code{stochvol_t}, \code{stochvol_nig}, \code{logperiodogram},
    \code{exponential(*)}, \code{weibull(*)} and \code{ps(*)}. *-marked
    models are for survival, and there is a special requirement for the
    response variable in the formula.
    }
  
  \item{data}{An optional data frame or a list containing the variables
    in the model.  If not found in \code{data}, the variables are taken
    from \code{environment(formula)}, typically the environment from
    which \code{inla} is called.}
  
  \item{quantiles}{A vector of quantiles,
    \eqn{p(0),p(1),\dots}{p(0),p(1),\ldots} to compute for each
    posterior marginal. The function returns, for each posterior
    marginal, the values \eqn{x(0),x(1),\dots}{x(0),x(1),\ldots} such
    that \deqn{\mbox{Prob}(X<x(p))=p}{Prob(X<x)=p} }
  
  \item{E}{Known component in the mean for the Poisson likelihoods
    defined as \deqn{E_i\exp(\eta_i)}{E exp(eta)} where
    \deqn{\eta_i}{eta} is the linear predictor. Should be \code{NULL} or
    a numeric vector of the same length of the data vector. See
    \cite{inla.model}. (Default 1)}

  \item{offset}{This can be used to specify an a-priori known component
          to be included in the linear predictor during fitting.  This
          should be \code{NULL} or a numeric vector of length either one
          or equal to the number of cases. One or more \code{offset()}
          terms can be included in the formula instead or as well, and
          if both are specified their sum is used.}

   \item{scale}{Fixed (optional) scale parameters of the precision for
          Gaussian and Student-T response models. Default values are all
          1.}
  
  \item{Ntrials}{A vector containing the number of trials for the
    \code{binomial} likelihood. (Default 1)}
  
  \item{verbose}{Boolean indicating if the \code{inla}-program should
    run in a verbose mode (default \code{FALSE}). }

  \item{control.compute}{see \code{?control.compute}}

  \item{control.inla}{see \code{?control.inla}}

  \item{control.results}{see \code{?control.results}}

  \item{control.fixed}{see \code{?control.fixed}}
  \item{control.mode}{see \code{?control.mode}}
  \item{control.expert}{see \code{?control.expert}}

  \item{only.hyperparam}{A boolean variable saying if only the
    hyperparameters are to be computed. Mainly internal use.}

  \item{inla.call}{The path to, or the name of, the
    \code{inla}-program. This is program is installed together with the
    \code{R}-package, but, for example, a native compiled version can be
    used instead to improve the performance.}

  \item{inla.arg}{A string indicating ALL arguments to the 'inla'
    program and do not include default arguments. (OOPS: This is an
    expert option!)}

  \item{num.threads}{Maximum number of threads the \code{inla}-program
    will use. For Windows this defaults to 1, otherwise its defaults to
    \code{NULL} (for which the system takes over control).}

  \item{keep}{A boolean variable indicating the working files (ini file,
    data files and results files) should be kept. If TRUE and no
    \code{working.directory} is specified the working files are stored
    in a directory called "inla"}

  \item{working.directory}{A string giving the name of an alternative
    directory where to store the working files}

  \item{silent}{A boolean variable defining how the
    \code{inla}-program should be ``silent''.}

  \item{disable.check}{Disable check for presence of the \code{inla}-program?}

  \item{debug}{If \code{TRUE}, then enable some debug output.}

  \item{user.hook}{This defines an optional user-defined function, which
    can be called just after the .ini-file is created, usually, to add
    extra information to the .ini-file. See the function
    \code{inla.user.hook} for an example of a such function, and the
    arguments to it.}

  \item{user.hook.args}{This defines an optional argument to \code{user.hook}.}
}


\value{%%

  \code{inla} returns an object of class \code{"inla"}. This is a list
  containing the following arguments:
  
  \item{summary.fixed}{Mean and standard deviation (plus, possibly
    quantiles and cdf) of the the fixed effects of the model.}

  \item{marginals.fixed}{Posterior marginal densities of the fixed
    effects of the model.}
  
  \item{summary.random}{List of mean and standard deviation (plus,
    possibly quantiles and cdf) of the the smooth or spatial effects
    defined through \code{f}.}
  
  \item{marginals.random}{If \code{return.marginals.random}=\code{TRUE}
    in \code{control.results} (default), posterior marginal densities of
    the random effects defined through \code{f}.}

  \item{model.random}{List of models for the smooth or spatial effects
    defined through \code{f}.}         

  \item{random.levels}{List of levels for the smooth or spatial effects
    defined through \code{f}.}

  \item{summary.linear.predictor}{Mean and standard deviation (plus,
    possibly quantiles and cdf) of the the linear fit on link scale.}
  
  \item{marginals.linear.predictor}{If
    \code{return.marginals.predictor}=\code{TRUE} in
    \code{control.predictor} (default), the posterior marginal densities of the
    linear fit on link scale.}

  \item{summary.fitted.values}{Mean and standard deviation (plus,
    possibly quantiles and cdf) of the fitted values obtained by
    transforming the linear predictors by the inverse of the link
    function.}

  \item{marginals.fitted.values}{If
    \code{return.marginals.random}=\code{TRUE} in
    \code{control.predictor} (default), the posterior marginal densities
    of the fitted values obtained by transforming the linear predictors
    by the inverse of the link function. }
  
  \item{nhyper}{Total number of hyperparameters in the model}
  
  \item{summary.hyperpar}{If \code{hyperpar}=\code{TRUE} in
    \code{control.compute}, mean and standard deviation of the posterior
    marginals for the hyperparameters.}
  
  \item{marginals.hyperpar}{If \code{hyperpar}=\code{TRUE} in
    \code{control.compute}, a list including the posterior marginals for
    the hyperparameters of the model, otherwise \code{NULL}. }

  \item{cpo}{If \code{cpo}=\code{TRUE} in \code{control.compute}, the
    values of conditional predictive ordinate (CPO) and the probability
    integral transform (PIT) for the model, otherwise \code{NULL}.}

   
  \item{mlik}{If \code{mlik}=\code{TRUE} in \code{control.compute}, the
    marginal likelihood of the model, otherwise \code{NULL}}
  
  \item{dic}{If \code{dic}=\code{TRUE} in \code{control.compute}, the
    deviance information criteria and effective number of parameters,
    otherwise \code{NULL}}
  
  \item{neffp}{Expected effective number of parameters in the model. The
    standard deviation of the expected number of parameters and the
    number of replicas for parameter are also returned}

  \item{control.predictor}{The value of the '\code{control.predictor}'
  argument used.}
  
  \item{control.inla}{The value of the '\code{control.inla}' argument
  used.}
  
  \item{control.data}{The value of the '\code{control.data}' argument
  used.}
  
  \item{call}{The matched call.}                   
  
  \item{formula}{The formula supplied}
  
  \item{family}{The '\code{family}' object used.}
%%
}
%%

\references{ Rue, H. and Martino, S. and Chopin, N. (2009)
\emph{Approximate Bayesian Inference for latent Gaussian models using
Integrated Nested Laplace Approximations, JRSS-series B (with
discussion)}, vol 71, no 2, pp 319-392.

  Rue, H and Held, L. (2005) \emph{Gaussian Markov Random Fields -
  Theory and Applications} Chapman and Hall
  
  Martino, S. and Rue, H. (2008) \emph{Implementing Approximate Bayesian
  Inference using Integrated Nested Laplace Approximation: a manual for
  the inla program} Preprint N.2 from Dep. of Mathematical Sciences
  (NTNU Norway)

}

\author{Sara Martino and Havard Rue \email{hrue@math.ntnu.no} }

\seealso{\code{\link{f}}, \code{\link{inla.models}},
  \code{\link{inla.hyperpar}} }

\examples{

### *******
### Example Surg. From WinBugs/OpenBugs Examples vol I: Binomial
### Regression with Random Effects.

data(Surg)
formula = r ~ f(hospital,model="iid",param=c(0.001,0.001))
mod.surg = inla(formula,data=Surg,family="binomial",Ntrials=n)

# improved estimation of the hyperparameters
h.surg = inla.hyperpar(mod.surg)

### *******
### Example Seeds. From WinBugs/OpenBugs Examples vol I: Random effect
### logistic regression

data(Seeds)
formula = r ~ x1*x2+f(plate,model="iid")
mod.seeds = inla(formula,data=Seeds,family="binomial",Ntrials=n)

## improved estimation of the hyperparameters
h.seeds = inla.hyperpar(mod.seeds)

### *******
### Example Salm. From WinBugs/OpenBugs Examples vol I: extra - Poisson
### variation in dose - response study

data(Salm)
mod.salm = inla(y~log(dose+10)+dose+f(rand,model="iid",param=c(0.001,0.001)),
                family="poisson",data=Salm)

## improved estimation of the hyperparameters
h.salm = inla.hyperpar(mod.salm)

### *******
### Example Epil. From WinBugs/OpenBugs Examples vol I: Mixed effect
### model with Poisson counts

data(Epil)

## build centered covariates
M.Trt  =  Epil$Trt-mean(Epil$Trt)
M.lBase4 = log(Epil$Base/4)-mean(log(Epil$Base/4))
M.V4 = Epil$V4-mean(Epil$V4)
M.lAge = log(Epil$Age)-mean(log(Epil$Age))
M.BT = Epil$Trt*M.lBase4-mean(Epil$Trt*M.lBase4)

Epil.data = list(y=Epil$y, Ind = Epil$Ind, rand = Epil$rand,
                 M.Trt=M.Trt, M.lBase4 = M.lBase4, M.V4 = M.V4,
                 M.lAge = M.lAge, M.BT = M.BT)
formula = y ~ M.lBase4 + M.Trt + M.BT+ M.lAge + M.V4 + 
              f(Ind,model="iid") + f(rand,model="iid")

model = inla(formula,family="poisson", data = Epil.data)


### *******
### Tokyo example. from Rue&Held "Gaussian Markov Random Fields"(Ch
### 4.3.4)

data(Tokyo)

## NB in this model there is no common mean, this has to be specified
## using -1 in the formula. This example will give a warning due to
## the flag cyclic=TRUE and the (general) use of localtions. In this
## example this is ok.
formula = y ~ f(time,model="rw2",cyclic=TRUE,param=c(1,0.0001)) - 1
mod.tokyo = inla(formula,family="binomial",Ntrials=n,data=Tokyo)


### *******
### Drivers example. From Rue&Held "Gaussian Markov Random Fields"(Ch
### 4.2.1). Time series model with seasonal and trend, and with Gaussian
### resposes.

data(Drivers)
formula  = sqrt(y) ~ belt + f(trend,model="rw2",param=c(1,0.0005)) +
                     f(seasonal,model="seasonal",season.length=12,param=c(1,0.1))
mod = inla(formula, family="gaussian", data=Drivers,
           control.data=list(param=c(4,4)))

### *******
### Disease Mapping. Oral cavity data example. From Rue & Martino 2007
### ``Approximate Inference for Hierarchical Gaussian Markov Random
### Field Models'', JSPI, 2007.

data(Oral)
g = system.file("demodata/germany.graph", package="INLA")

## the data set has to contain one separate column for each term
## specified through f() so in this case it is necessary to add one
## column to the data set

Oral = cbind(Oral,region.struct=Oral$region)
formula = Y ~ f(region.struct,model="besag",graph.file=g,param=c(1,0.01)) +
              f(region,model="iid",param=c(1,0.01))
mod = inla(formula,family="poisson", data=Oral, E=E)

### *******
### Disease mapping with covariate. Larynx data example, from Rue &
### Martino 2007 ``Approximate Inference for Hierarchical Gaussian
### Markov Random Field Models'' JSPI.

data(Germany)
g = system.file("demodata/germany.graph", package="INLA")

## the data set has to contain one separate column for each term
## specified through f() so in this case it is necessary to add one
## column to the data set

Germany = cbind(Germany,region.struct=Germany$region)
formula = Y ~ f(region.struct,model="besag", graph.file=g,
                param=c(1,0.00005)) + f(x,model="rw2",param=c(1,0.05)) +
                f(region,model="iid")
mod  =  inla(formula,family="poisson",data=Germany,E=E)


### *******
### Example 4.3.5 from the GMRF-book by Rue & Held. 
data(Cancer)
g = system.file("demodata/ddr.graph", package="INLA")

## the data set has to contain one separate column for each term
## specified through f() so in this case it is necessary to add one
## column to the data set

Cancer = cbind(Cancer,region.struct=Cancer$region)
formula  =  Y ~ f(region.struct,model="besag",graph.file=g,param=c(1,0.01)) + 
                f(Age,model="rw2",param=c(1,0.001)) +
                f(region,model="iid",constr=TRUE)

mod  =  inla(formula,family="binomial",data=Cancer,Ntrial=N, control.inla=list(h=0.05))



### *******
### Example 4.2.2 from the GMRF-book by Rue & Held: Munich rent data

data(Munich)
g = system.file("demodata/munich.graph", package="INLA")

## Note that here we what to have an estimator of the effect of `year`
## also the for years where we have no observation, therefore we give a
## vector with all possible values assumed by the covariate `year`, that
## is seq(1918,2001)

formula = rent ~ f(location,model="besag",graph.file=g,initial=1,param=c(1,0.01)) +
                 f(year,model="rw2",values=seq(1918,2001),param=c(1,0.01)) +
                 f(floor.size,model="rw2",param=c(1,0.01)) +
                 Gute.Wohnlage + Beste.Wohnlage + Keine.Wwv + Keine.Zh +
                 Kein.Badkach  + Besond.Bad + Gehobene.Kueche +
                 zim1 + zim2 + zim3 + zim4 + zim5 + zim6 -1
mod  =  inla(formula,data=Munich, verbose = TRUE,
             control.data=list(initial=-1),
             control.predictor = list(initial = 12),
             control.inla = list(h=1e-4))
 

### *******
### Zambia example. From the BayesX manual; see
### http://www.stat.uni-muenchen.de/~bayesx/bayesx.html

data(Zambia)
g = system.file("demodata/zambia.graph", package="INLA")
Zambia$distr.unstruct  =  Zambia$district

formula = hazstd ~ f(bmi.distr,model="rw2") + f(agc,model="rw2") +
                   f(district,model="besag",graph.file=g, param=c(1,0.01)) +
                   f(distr.unstruct,model="iid",param=c(1,0.01)) +
                   rcw+edu1+ edu2+ tpr+ sex
mod.zambia  =  inla(formula,data=Zambia, control.data=list(initial = 1),
                    control.inla=list(h=1e-4),
                    control.predictor=list(initial = 10),
                    verbose=TRUE)

### *******
#### Lip-cancer Example. From the OpenBUGS\GeoBUGS Examples Manual.

data(Scotland)
g = system.file("demodata/scotland.graph", package="INLA")
formula = Counts~f(Region,model="besag", graph.file = g,
                     param=c(0.5,0.0005)) + I(X/10)
mod.scotland = inla(formula,family="poisson",E=E,data=Scotland)
}



