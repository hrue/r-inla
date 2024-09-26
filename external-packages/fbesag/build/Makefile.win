CFLAGS = $(FLAGS) -O3 -mtune=generic -fopenmp  -fpic -I.
LDFLAGS = 
TARGET = libfbesag.so
OBJ = $(patsubst %.c,%.o,$(wildcard *.c))

$(TARGET) : $(OBJ)
	$(CC) $(LDFLAGS) -o $@ -shared $^ $(LIBBLAS) -lgfortran -lgomp -lm

clean :; rm -f $(TARGET) $(OBJ)
