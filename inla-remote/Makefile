SHELL=/bin/bash
# Only try to update the hgstamp if hg is available.  Should also
# check if we are in a repository.  To use existing HG file, use "make
# HG= target" Would be prettier and more robust with a mercurial
# update hook.
HG:= $(shell which hg)
HGVERSION:= $(shell if [ -x "$(HG)" ]; then $(HG) parents --template 'hgid: {node|short}  date: {date|date}'; elif [ -e hgstamp ]; then cat hgstamp; else echo 'hgid: unknown  date: unknown' ; fi)


all : inla.remote inla.remote.cygwin

inla.remote : inla.remote.core
	rm -f $@
	{ echo '#!/usr/bin/env bash'; echo; echo -e HGVERSION=\""$(HGVERSION)"\"; echo; cat $^; } > $@
	chmod 755 $@

inla.remote.cygwin : inla.remote.core
	rm -f $@
	{ echo '#!/usr/bin/env.exe bash.exe'; echo; echo -e HGVERSION=\""$(HGVERSION)"\"; echo; cat $^; } > $@
	chmod 755 $@

.PHONY : all
