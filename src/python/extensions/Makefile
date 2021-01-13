PYINC=-I/usr/include/python2.4/
NINC=-I$(HOME)/usr/lib64/python2.4/site-packages/numpy/core/include/numpy


all: _slice_update.so _getalpha.so
clean:
	rm -f _getalpha.so _slice_update.so
	rm -f alpha.o get_alpha.o slice_update.o

# ---- Link --------------------------- 
_slice_update.so: slice_update.o
	gcc -m64 -shared -o _slice_update.so slice_update.o -fPIC

_getalpha.so: alpha.o get_alpha.o
	gcc -m64 -shared -o _getalpha.so get_alpha.o alpha.o -fPIC -lgfortran

alpha.o: alpha.f90
	gfortran -m64 -O3 -c alpha.f90 $(QINC) -fPIC 

get_alpha.o: get_alpha.c alpha.o
	gcc -O3 -m64 -c get_alpha.c $(NINC) $(PYINC) -fPIC


# ---- gcc C compile ------------------
slice_update.o: slice_update.c 
	gcc -m64 -O3 -c slice_update.c $(NINC) $(PYINC) -fPIC

