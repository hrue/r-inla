CXXFLAGS = $(FLAGS) -O3 -mtune=generic -fopenmp -pipe -fpic -I.
LDFLAGS = 
TARGET = libfbesag.so
OBJ = $(patsubst %.cpp,%.o,$(wildcard *.cpp))

$(TARGET) : $(OBJ)
	$(CXX) $(LDFLAGS) -o $@ -shared $^ $(LIBBLAS) -lgfortran -lgomp -lm

clean :; rm -f $(TARGET) $(OBJ)
