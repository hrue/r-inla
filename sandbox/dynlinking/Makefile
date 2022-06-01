CFLAGS= 
test : test.o
	$(CC) -o $@ $^ -ldl

test2 : test2.o
	$(CC) -o $@ $^ -lltdl

test3 : test3.o libelias.so
	$(CC) -o $@ test3.o -lltdl -lm

libelias.so : elias.c
	$(CC) -c -fPIC elias.c
	$(CC) -shared -o $@ elias.o
