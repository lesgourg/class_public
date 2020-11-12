CXX = g++
CFLAGS = -O2 -fopenmp -I../include -I../external/HyRec2020 -I../external/RecfastCLASS -I../external/heating
CLASSMODULES = ../build/arrays.o ../build/background.o ../build/common.o \
	../build/dei_rkck.o ../build/distortions.o ../build/energy_injection.o \
	../build/evolver_ndf15.o ../build/evolver_rkck.o ../build/growTable.o \
	../build/helium.o ../build/history.o ../build/hydrogen.o \
	../build/hyperspherical.o ../build/hyrectools.o \
	../build/injection.o ../build/input.o ../build/lensing.o \
	../build/noninjection.o ../build/nonlinear.o ../build/output.o \
	../build/parser.o ../build/perturbations.o ../build/primordial.o \
	../build/quadrature.o ../build/sparse.o ../build/spectra.o \
	../build/thermodynamics.o ../build/transfer.o \
	../build/trigonometric_integrals.o ../build/wrap_hyrec.o ../build/wrap_recfast.o

all: testKlass Makefile

testKlass: testKlass.o Engine.o ClassEngine.o
	$(CXX) $(CFLAGS) $(CLASSMODULES) ClassEngine.o Engine.o testKlass.o -o testKlass

testKlass.o: testKlass.cc
	$(CXX) $(CFLAGS) -c testKlass.cc -o testKlass.o

ClassEngine.o: ClassEngine.cc ClassEngine.hh
	$(CXX) $(CFLAGS) -c ClassEngine.cc -o ClassEngine.o

Engine.o: Engine.cc Engine.hh
	$(CXX) $(CFLAGS) -c Engine.cc -o Engine.o

clean:
	rm -rf *.o testKlass
