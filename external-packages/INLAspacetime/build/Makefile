CFLAGS = $(FLAGS) -O3 -mtune=generic -fopenmp -pipe -fpic -I.
LDFLAGS = 
TARGET = libINLAspacetime.so
OBJ = $(patsubst %.c,%.o,$(wildcard *.c))

$(TARGET) : $(OBJ)
	$(CC) $(LDFLAGS) -o $@ -shared $^ $(LIBBLAS) -lgfortran -lgomp -lm

clean :; rm -f $(TARGET) $(OBJ)
