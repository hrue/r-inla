SHELL=/bin/bash

all : inla.remote inla.remote.cygwin

inla.remote : inla.remote.core
	rm -f $@
	{ echo '#!/usr/bin/env bash'; echo; echo "CYGWIN=0"; echo; cat $^; } > $@
	chmod 755 $@

inla.remote.cygwin : inla.remote.core
	rm -f $@
	{ echo '#!/usr/bin/env.exe bash.exe'; echo; echo "CYGWIN=1"; echo; cat $^; } > $@
	chmod 755 $@

.PHONY : all
