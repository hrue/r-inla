OBJ = xtmpl.o predicates.o mesh.o meshc.o main.o vector.o ioutils.o
FLAGS = -g -Wall
CXXFLAGS = -g -Wall
#CXXFLAGS = -g -Wall -DDEBUG
LDFLAGS = -lX11

all: main

$(OBJ): Makefile
ioutils.o: ioutils.tpp ioutils.h
vector.o: vector.h
mesh.o: mesh.h vector.h
meshc.o: meshc.h mesh.h vector.h
main.o: fmesher.h mesh.h meshc.h predicates.h vector.h ioutils.h ioutils.tpp
predicates.o: predicates.h

main: $(OBJ) Makefile
	g++ $(FLAGS) $(CXXFLAGS) $(LDFLAGS) $(OBJ) -o main
