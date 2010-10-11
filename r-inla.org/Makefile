all:; 
	@ find . -type d -mindepth 1| while read d; do if [ -f "$$d/Makefile" ]; then $(MAKE) -C "$$d" ;  fi; done
