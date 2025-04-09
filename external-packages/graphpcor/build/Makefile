CFLAGS = $(FLAGS) $(INC) -O3 -mtune=generic -fopenmp  -fpic -I.
LDFLAGS = $(LIB)
TARGET = libgraphpcor.so
OBJ = $(patsubst %.c,%.o,$(wildcard *.c))

$(TARGET) : $(OBJ)
	$(CC) $(LDFLAGS) -o $@ -shared $^ $(LIBBLAS) -lgfortran -lgomp $(LIBR) $(LIBLTDL) -lm

clean :; rm -f $(TARGET) $(OBJ)
