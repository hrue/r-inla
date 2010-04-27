OBJ = mesh.o predicates.o main.o xtmpl.o
FLAGS = -g -pg -Wall
CXXFLAGS = -g -pg -Wall
LDFLAGS = -lX11

all: main

main: $(OBJ) mesh.h Makefile
	g++ $(FLAGS) $(CXXFLAGS) $(LDFLAGS) $(OBJ) -o main
