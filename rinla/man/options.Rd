\name{inla.option}
\alias{inla.option}
\alias{inla.options}
\alias{inla.setOption}
\alias{inla.getOption}

\title{Set and get global options for INLA}

\description{Set and get global options for INLA}
\usage{
inla.setOption(option, value)
inla.getOption(option)
}

\arguments{

  \item{option}{The option to set. If \code{option = NULL} then
    \code{inla.getOption} then \code{inla.getOption} will display the
    current defaults, otherwise, \code{option} must be one of

    inla.call: The path to the inla-program.

    inla.arg: Additional arguments to \code{inla.call}

    num.threads: Number of threads to use.
    
    silent: Run the inla-program in a silent mode?

    disable.check: Check the inla-program path upfront?

    keep: Keep temporary files?

    working.directory: The name of the working directory.

    debug: Print some debug information in the R-interface

    value: The value of \code{option}

    cygwin: The home of the CYGWIN installation (default
    C:/cygwin). [Remote computing for Windows only]

    cygwin.home: The $HOME directory in the cygwin-shell. [Remote
    computing for Windows only]

    ssh.auth.sock: The ssh bind-adress (value of $SSH_AUTH_SOCK int he
    cygwin-shell). [Remote computing for Windows only]

The options are stored in the variable \code{inla.options} in the
    \code{.GlobalEnv}-environment.
  }
}

\author{Havard Rue \email{hrue@math.ntnu.no}}

\examples{

## Display all default options
inla.getOption()
## set number of threads
inla.setOption("num.threads", 2)
## check it
inla.getOption("num.threads")
## this is where its stored
inla.options
## remove it
rm("inla.options")
## and we are back to the default value
inla.getOption("num.threads")
}


