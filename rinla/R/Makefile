
all: TAGS make-man-pages.timestamp

make-man-pages.timestamp: $(wildcard *.R)
	./make-man-pages
	touch $@

clean :; rm make-man-pages.timestamp

TAGS : $(wildcard *.R)
	etags $^

.PHONY : all clean TAGS all

