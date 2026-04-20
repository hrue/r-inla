SHELL=/bin/bash

all : inla.remote inla.q inla.submit 

inla.remote : inla.remote.core
	rm -f $@
	{ echo '#!/usr/bin/env bash'; echo; cat $^; } > $@
	chmod 755 $@

inla.q : inla.q.core
	rm -f $@
	{ echo '#!/usr/bin/env bash'; echo; cat $^; } > $@
	chmod 755 $@

inla.submit : inla.submit.core
	rm -f $@
	{ echo '#!/usr/bin/env bash'; echo; cat $^; } > $@
	chmod 755 $@

.PHONY : all
