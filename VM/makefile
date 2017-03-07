CPPFLAGS += -O2
# source
SRC = $(wildcard *.cpp)
# source object
SRC_OBJ = $(SRC:.cpp=.o)
# source binary
SRC_BIN =main.out

$(SRC_BIN): $(SRC_OBJ)
	g++ -o $@ $^ $(CPPFLAGS)

clean:
	rm *o
