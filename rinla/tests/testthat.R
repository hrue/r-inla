# For development, run with devtools::test() from the rinla/ directory
# Requires the binaries to be locally available via the relevant
#   rinla/inst/bin/... directory

library(testthat)
library(INLA)

test_check("INLA")
