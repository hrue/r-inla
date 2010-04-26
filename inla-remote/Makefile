SHELL=/bin/bash

all : inla.remote inla.remote.cygwin

inla.remote : inla.remote.core
	rm -f $@
	{ echo '#!/usr/bin/env bash'; cat $^; } > $@
	chmod 755 $@

inla.remote.cygwin : inla.remote.core
	rm -f $@
	{ echo '#!/usr/bin/env.exe bash.exe'; cat $^; } > $@
	chmod 755 $@

.PHONY : all
