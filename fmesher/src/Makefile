OBJ = xtmpl.o predicates.o mesh.o meshc.o vector.o ioutils.o
FLAGS = -g -Wall
CXXFLAGS = -g -Wall
#CXXFLAGS = -g -Wall -DDEBUG
LDFLAGS = -lX11

all: fmesher main

$(OBJ): Makefile
ioutils.o: ioutils.tpp ioutils.h
vector.o: vector.h
mesh.o: mesh.h vector.h
meshc.o: meshc.h mesh.h vector.h
main.o: main.cpp fmesher.h mesh.h meshc.h predicates.h vector.h ioutils.h ioutils.tpp
fmesher.o: fmesher.cpp fmesher.h mesh.h meshc.h predicates.h vector.h ioutils.h ioutils.tpp
predicates.o: predicates.h

fmesher: $(OBJ) Makefile fmesher.o
	g++ $(FLAGS) $(CXXFLAGS) $(LDFLAGS) $(OBJ) fmesher.o -o fmesher

main: $(OBJ) Makefile main.o
	g++ $(FLAGS) $(CXXFLAGS) $(LDFLAGS) $(OBJ) main.o -o main
