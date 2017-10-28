CXX = icc
CPPFLAGS += -O2 -std=c++11 -fPIC

#FFTW
FFTW_INC := -I/home-yw/users/nsyw449_YK/dy/Program/fftw3/include
FFTW_LIB := -L/home-yw/users/nsyw449/dy/Program/fftw3/lib
#python
PYTHON_INC := -I/home-yw/users/nsyw449_YK/dy/Program/anaconda3/include/python3.6m
NUMPY_INC := -I/home-yw/users/nsyw449_YK/dy/Program/anaconda3/pkgs/numpy-1.11.3-py36_0/lib/python3.6/site-packages/numpy/core/include
PYTHON_LIB := -L/home-yw/users/nsyw449_YK/dy/Program/anaconda3/lib
LIBS += -lfftw3 -lm -lpython3

_VMpy.so: comn.o grid.o io_data.o node.o run.o spatial_corr.o VMpy.o VMpy_wrap.o
	$(CXX) -shared comn.o grid.o io_data.o node.o run.o spatial_corr.o VMpy.o VMpy_wrap.o -o _VMpy.so $(PYTHON_LIB) $(FFTW_LIB) $(LIBS)

VMpy_wrap.o: VMpy_wrap.cxx
	$(CXX) $(CPPFLAGS) -c VMpy_wrap.cxx $(PYTHON_INC) $(NUMPY_INC)
VMpy_wrap.cxx: VMpy.h VMpy.cpp VMpy.i
	swig -python -c++ VMpy.i
VMpy.o: VMpy.cpp VMpy.h
	$(CXX) $(CPPFLAGS) -c VMpy.cpp
run.o: run.h run.cpp grid.h io_data.h
	$(CXX) $(CPPFLAGS) -c run.cpp
io_data.o: io_data.cpp io_data.h node.h spatial_corr.h
	$(CXX) $(CPPFLAGS) -c io_data.cpp
spatial_corr.o: spatial_corr.h spatial_corr.cpp
	$(CXX) $(CPPFLAGS) -c spatial_corr.cpp -lfftw3 -lm $(FFTW_INC) $(FFTW_LIB)
grid.o: node.h grid.h grid.cpp
	$(CXX) $(CPPFLAGS) -c grid.cpp
node.o: node.cpp node.h rand.h comn.h
	$(CXX) $(CPPFLAGS) -c node.cpp
comn.o: comn.cpp comn.h
	$(CXX) $(CPPFLAGS) -c comn.cpp


clean:
	rm *o *cxx
