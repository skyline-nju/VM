CXX = g++
CPPFLAGS += -O2 -std=c++11
# source
SRC = $(wildcard *.cpp)
# source object
SRC_OBJ = $(SRC:.cpp=.o)
# source binary
SRC_BIN =main.out

$(SRC_BIN): $(SRC_OBJ)
	$(CXX) $(CPPFLAGS) -o $@ $^ 

clean:
	rm *o
