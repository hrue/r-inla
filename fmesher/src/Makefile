OBJ = mesh.o predicates.o main.o xtmpl.o
FLAGS = -g -Wall
CXXFLAGS = -g -Weffc++ -Wall
LDFLAGS = -lX11

all: main

main: $(OBJ)
	g++ $(FLAGS) $(CXXFLAGS) $(LDFLAGS) $(OBJ) -o main
