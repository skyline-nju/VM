CXX = icc
CPPFLAGS += -O2 -std=c++11 -fPIC

#python
PYTHON_INC := -I/home-yw/users/nsyw449_YK/dy/Program/anaconda3/include/python3.6m
NUMPY_INC := -I/home-yw/users/nsyw449_YK/dy/Program/anaconda3/pkgs/numpy-1.11.3-py36_0/lib/python3.6/site-packages/numpy/core/include
PYTHON_LIB := -L/home-yw/users/nsyw449_YK/dy/Program/anaconda3/lib
LIBS += -lpython3

_XYmodel.so: XYmodel.o XYmodel_wrap.o
	$(CXX) -shared XYmodel.o XYmodel_wrap.o -o _XYmodel.so $(PYTHON_LIB) $(LIBS)

XYmodel.o: XYmodel.h XYmodel.cpp rand.h
	$(CXX) $(CPPFLAGS) -c XYmodel.cpp

XYmodel_wrap.o: XYmodel_wrap.cxx
	$(CXX) $(CPPFLAGS) -c XYmodel_wrap.cxx $(PYTHON_INC) $(NUMPY_INC)

XYmodel_wrap.cxx: XYmodel.i XYmodel.cpp XYmodel.h
	swig -python -c++ XYmodel.i


clean:
	rm *o *cxx