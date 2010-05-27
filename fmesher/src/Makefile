OBJ = xtmpl.o predicates.o mesh.o meshc.o main.o vector.o
FLAGS = -g -Wall
CXXFLAGS = -g -Wall
#CXXFLAGS = -g -Wall -DDEBUG
LDFLAGS = -lX11

all: main

$(OBJ): Makefile
mesh.o: mesh.h
meshc.o: meshc.h mesh.h
main.o: fmesher.h mesh.h meshc.h predicates.h vector.h
predicates.o: predicates.h

main: $(OBJ) Makefile
	g++ $(FLAGS) $(CXXFLAGS) $(LDFLAGS) $(OBJ) -o main
