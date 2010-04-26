OBJ = mesh.o predicates.o main.o xtmpl.o
FLAGS = -g -Wall
LDFLAGS = -lX11

all: main

main: $(OBJ)
	g++ $(LDFLAGS) $(OBJ) -o main
