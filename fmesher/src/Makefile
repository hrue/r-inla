OBJ = xtmpl.o predicates.o mesh.o meshc.o vector.o ioutils.o
FLAGS = -g -Wall
CFLAGS = $(FLAGS)
CXXFLAGS = -g -Wall
#CXXFLAGS = -g -Wall -DDEBUG
LDFLAGS = -lX11
CC = gcc
CXX = g++

ALL = fmesher main

all: $(ALL)

$(OBJ): Makefile
ioutils.o: ioutils.tcc ioutils.hh
vector.o: vector.hh
mesh.o: mesh.hh vector.hh
meshc.o: meshc.hh mesh.hh vector.hh
main.o: main.cc fmesher.hh mesh.hh meshc.hh predicates.hh vector.hh ioutils.hh ioutils.tcc
fmesher.o: fmesher.cc fmesher.hh mesh.hh meshc.hh predicates.hh vector.hh ioutils.hh ioutils.tcc
predicates.o: predicates.hh

fmesher: $(OBJ) Makefile fmesher.o
	g++ $(FLAGS) $(CXXFLAGS) $(LDFLAGS) $(OBJ) fmesher.o -o fmesher

main: $(OBJ) Makefile main.o
	g++ $(FLAGS) $(CXXFLAGS) $(LDFLAGS) $(OBJ) main.o -o main

clean :;
	rm $(ALL) *.o

.PHONY : clean
