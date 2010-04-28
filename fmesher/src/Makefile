OBJ = mesh.o predicates.o main.o xtmpl.o
FLAGS = -g -pg -Wall
CXXFLAGS = -g -pg -Wall
LDFLAGS = -lX11

all: main

mesh.o: mesh.h Makefile
main.o: mesh.h Makefile

main: $(OBJ) Makefile
	g++ $(FLAGS) $(CXXFLAGS) $(LDFLAGS) $(OBJ) -o main
