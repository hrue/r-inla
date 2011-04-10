
all: TAGS make-man-pages.timestamp

make-man-pages.timestamp: $(wildcard *.R)
	./make-man-pages
	touch $@

TAGS : $(wildcard *.R)
	etags $^

