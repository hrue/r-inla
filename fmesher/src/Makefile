OBJ = xtmpl.o predicates.o mesh.o meshc.o main.o
FLAGS = -g -Wall
CXXFLAGS = -g -pg -Wall
LDFLAGS = -lX11

all: main

$(OBJ): Makefile $(wildcard *.h)
mesh.o: mesh.h
meshc.o: meshc.h mesh.h
main.o: fmesher.h mesh.h meshc.h

main: $(OBJ) Makefile
	g++ $(FLAGS) $(CXXFLAGS) $(LDFLAGS) $(OBJ) -o main
