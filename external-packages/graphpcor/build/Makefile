CFLAGS = $(FLAGS) $(RINC) -O3 -mtune=generic -fopenmp  -fpic -I.
LDFLAGS = $(LLIB)
TARGET = libgraphpcor.so
OBJ = $(patsubst %.c,%.o,$(wildcard *.c))

$(TARGET) : $(OBJ)
	$(CC) $(LDFLAGS) -o $@ -shared $^ $(LIBBLAS) -lgfortran -lgomp -lltdl -lm

clean :; rm -f $(TARGET) $(OBJ)
