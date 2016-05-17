
all: TAGS make-man-pages.timestamp

make-man-pages.timestamp: $(wildcard *.R)
	./make-man-pages
	touch $@

clean :; rm make-man-pages.timestamp

TAGS : $(wildcard *.R)
	R --vanilla -e 'library(utils); rtags(path=".", ofile="TAGS")'

.PHONY : all clean

