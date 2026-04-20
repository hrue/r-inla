CC = gcc
CFLAGS = -O3 -march=native -DSHOW_TEST_OUTPUT

test : test.o sparse-vec.o
	$(CC) -o $@ $^ -lblas -lgomp -lgsl -lgslcblas -lm
