This directory contains code for the Leukaemia example used in
F Lindgren, H Rue, J Lindström
"An explicit link between Gaussian fields and Gaussian Markov
 random fields: The stochastic partial differential equation approach"
J.R.S.S. series B, to be read March 16th, 2011
http://www.r-inla.org/events/readpaperatrss16thmarch2011

Unpack, and run
 source("leuk-demo.R")
The comments in the code should hopefully explain what the different
parts do, with further explanations below.

You will need the "testing" version of INLA;
 inla.upgrade(testing=TRUE)

The spde/gmrf/inla interface is being updated from the current
(2011-02-16) slightly ad hoc version into something more streamlined,
in particular for plotting.  We will try to keep backwards
compatibility and/or update the contents of this demonstration code.
Update (2011-03-08): The code has been updated to use the latest
version new inla.mesh/spde interface.

In the code, element  Q[i,j]  of the precision matrix corresponds to
the link between the weights for triangulation basis nodes  i  and  j
with coordinates stored in   mesh$loc[i,]   and   mesh$loc[j,]

The  proj  object contains information for mapping between the
triangulation and a dense lattice of points; used in the code for
plotting purposes. An example of how to use this is in the plotting
part of the code.

For each data point,  Leuk$district, is the region index, so that
mesh.district = c()
mesh.district[mesh$idx$loc] = Leuk$district
should give the mesh index mapped to the triangulation vertices.
No such mapping is available for the dense lattice in  proj.

Contact Finn Lindgren, finn.lindgren@math.ntnu.no, if you have
questions about the code.

Finn Lindgren, 2011-03-08
