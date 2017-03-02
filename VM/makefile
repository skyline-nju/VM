all: main.out

main.out: main.o output.o grid.o node.o disorder.o
	icc -o main.out main.o output.o grid.o node.o disorder.o

main.o: main.cpp disorder.h output.h
	icc -c main.cpp

output.o: output.cpp output.h
	icc -c output.cpp

grid.o: grid.cpp grid.h
	icc -c grid.cpp

node.o: node.cpp node.h rand.h
	icc -c node.cpp

disorder.o: disorder.cpp disorder.h
	icc -c disorder.cpp

clean:
	rm *o
