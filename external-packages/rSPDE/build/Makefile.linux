CFLAGS = $(FLAGS) -O3 -mtune=generic -fopenmp  -fpic -I.
CXXFLAGS = $(FLAGS) -O3 -mtune=generic -fopenmp  -fpic -I. -I$(EIGEN)
LDFLAGS = 
TARGET = librSPDE.so
OBJ = $(patsubst %.c,%.o,$(wildcard *.c)) $(patsubst %.cpp,%.o,$(wildcard *.cpp))

$(TARGET) : $(OBJ)
	$(CXX) $(LDFLAGS) -o $@ -shared $^ $(LIBBLAS) -lgfortran -lgomp -lm

clean :; rm -f $(TARGET) $(OBJ)
