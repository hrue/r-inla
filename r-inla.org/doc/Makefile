all:; @(for d in *; do test -d $$d &&  make -C $$d; done)
clean:; @(for d in *; do test -d $$d && make -C $$d $@; done)
.PHONY: all clean


