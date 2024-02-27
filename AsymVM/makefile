# choose compiler: icc or g++
CXX = g++
CPPFLAGS += -O2 -std=c++11

# include path
CFLAGS += -I/home-yw/users/nsyw449_YK/dy/Program/fftw3/include
# lib path
LDFLAGS += -L/home-yw/users/nsyw449/dy/Program/fftw3/lib
LIBS += -lfftw3 -lm
# source
SRC = $(wildcard *.cpp)
# source object
SRC_OBJ = $(SRC:.cpp=.o)
# source binary
SRC_BIN =main.out

$(SRC_BIN): $(SRC_OBJ)
	$(CXX) $(CPPFLAGS) $(CFLAGS) $(LDFLAGS)-o $@ $^ $(LIBS)

clean:
	rm *o
