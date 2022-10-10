CC = gcc
CFLAGS = -O3 -march=native -DSHOW_TEST_OUTPUT

test : test.o dot-product.o
	$(CC) -o $@ $^ -lblas -lgomp -lgsl -lgslcblas -lm
