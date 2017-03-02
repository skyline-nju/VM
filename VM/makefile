all: main.out

main.out: main.o output.o grid.o node.o run.o
	icc -o main.out main.o output.o grid.o node.o run.o

main.o: main.cpp run.h
	icc -c main.cpp

run.o: run.cpp grid.h output.h
	icc -c run.cpp

output.o: output.cpp output.h
	icc -c output.cpp

grid.o: grid.cpp grid.h
	icc -c grid.cpp

node.o: node.cpp node.h rand.h
	icc -c node.cpp

clean:
	rm *o
