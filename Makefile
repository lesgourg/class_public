#Some Makefile for CLASS.
#Julien Lesgourgues, 28.11.2011

MDIR := $(shell pwd)
WRKDIR = $(MDIR)/build

.base:
	if ! [ -e $(WRKDIR) ]; then mkdir $(WRKDIR) ; mkdir $(WRKDIR)/lib; fi;
	touch build/.base

vpath %.h source:tools:main:include
vpath %.c source:tools:main:test
vpath %.cpp source:tools:main:test
vpath %.o build
vpath %.opp build
vpath .base build

########################################################
###### LINES TO ADAPT TO YOUR PLATEFORM ################
########################################################

# your C compiler:
CC        = gcc
CXX       = g++

# your tool for creating static libraries:
AR        = ar rv

# Your python interpreter.
# In order to use Python 3, you can manually
# substitute python3 to python in the line below, or you can simply
# add a compilation option on the terminal command line:
# "PYTHON=python3 make all" (THanks to Marius Millea for pyhton3
# compatibility)
PYTHON ?= python

# your optimization flag
OPTFLAG = -O3 -ffast-math #-march=native

# all other compilation flags
CCFLAG = -g -fPIC
CXXFLAG = $(CCFLAG) -std=c++11 -Wno-write-strings
LDFLAG = -g -fPIC
LIBRARIES = -lm -lpthread

# leave blank to compile without HyRec, or put path to HyRec directory
# (with no slash at the end: e.g. hyrec or ../hyrec)
HYREC = hyrec

########################################################
###### IN PRINCIPLE THE REST SHOULD BE LEFT UNCHANGED ##
########################################################

# pass current working directory to the code
CCFLAG += -D__CLASSDIR__='"$(MDIR)"'

# where to find include files *.h
INCLUDES = -I../include -I../tools -I../source

# automatically add external programs if needed. First, initialize to blank.
EXTERNAL =

# eventually update flags for including HyRec
ifneq ($(HYREC),)
vpath %.c $(HYREC)
CCFLAG += -DHYREC
#LDFLAGS += -DHYREC
INCLUDES += -I../hyrec
EXTERNAL += hyrectools.o helium.o hydrogen.o history.o
endif
.SUFFIXES: .c .cpp .o .opp .h

# We could let gcc generate dependency information automatically, see this link:
# https://make.mad-scientist.net/papers/advanced-auto-dependency-generation/
# However, a clean build of CLASS is so fast that we just rebuild everything if *any*
# .h-file changed.
H_ALL = $(notdir $(wildcard include/*.h) $(wildcard tools/*.h) $(wildcard source/*.h))

%.o: %.c .base $(H_ALL)
	cd $(WRKDIR);$(CC) $(OPTFLAG) $(CCFLAG) $(INCLUDES) -c ../$< -o $*.o

%.opp: %.cpp .base $(H_ALL)
	cd $(WRKDIR); $(CXX) $(OPTFLAG) $(CXXFLAG) $(INCLUDES) -c ../$< -o $*.opp

TOOLS_O = growTable.o dei_rkck.o sparse.o evolver_rkck.o evolver_ndf15.o arrays.o parser.opp quadrature.o hyperspherical.o common.o trigonometric_integrals.o

TOOLS_OPP = non_cold_dark_matter.opp exceptions.opp

TOOLS = $(TOOLS_O) $(TOOLS_OPP)

SOURCE = input_module.opp background_module.opp thermodynamics_module.opp perturbations_module.opp primordial_module.opp nonlinear_module.opp transfer_module.opp spectra_module.opp lensing_module.opp cosmology.opp

OUTPUT = output_module.opp

CLASS = class.opp

TEST_LOOPS = test_loops.o

TEST_LOOPS_OMP = test_loops_omp.o

TEST_DEGENERACY = test_degeneracy.o

TEST_SPECTRA = test_spectra.o

TEST_TRANSFER = test_transfer.o

TEST_NONLINEAR = test_nonlinear.o

TEST_PERTURBATIONS = test_perturbations.o

TEST_THERMODYNAMICS = test_thermodynamics.o

TEST_BACKGROUND = test_background.o

TEST_SIGMA = test_sigma.o

TEST_HYPERSPHERICAL = test_hyperspherical.o

TEST_STEPHANE = test_stephane.o

all: class libclass.a classy

libclass.a: $(TOOLS) $(SOURCE) $(EXTERNAL)
	$(AR)  $@ $(addprefix build/, $(TOOLS) $(SOURCE) $(EXTERNAL))

class: $(TOOLS) $(SOURCE) $(EXTERNAL) $(OUTPUT) $(CLASS)
	$(CXX) $(OPTFLAG) $(LDFLAG) -o class $(addprefix build/,$(notdir $^)) $(LIBRARIES)

test_sigma: $(TOOLS) $(SOURCE) $(EXTERNAL) $(OUTPUT) $(TEST_SIGMA)
	$(CC) $(OPTFLAG) $(LDFLAG) -o test_sigma $(addprefix build/,$(notdir $^)) $(LIBRARIES)

test_loops: $(TOOLS) $(SOURCE) $(EXTERNAL) $(OUTPUT) $(TEST_LOOPS)
	$(CXX) $(OPTFLAG) $(LDFLAG) -o $@ $(addprefix build/,$(notdir $^)) $(LIBRARIES)

test_loops_omp: $(TOOLS) $(SOURCE) $(EXTERNAL) $(OUTPUT) $(TEST_LOOPS_OMP)
	$(CC) $(OPTFLAG) $(LDFLAG) -o $@ $(addprefix build/,$(notdir $^)) $(LIBRARIES)

test_stephane: $(TOOLS) $(SOURCE) $(EXTERNAL) $(OUTPUT) $(TEST_STEPHANE)
	$(CC) $(OPTFLAG) $(LDFLAG) -o $@ $(addprefix build/,$(notdir $^)) $(LIBRARIES)

test_degeneracy: $(TOOLS) $(SOURCE) $(EXTERNAL) $(OUTPUT) $(TEST_DEGENERACY)
	$(CC) $(OPTFLAG) $(LDFLAG) -o $@ $(addprefix build/,$(notdir $^)) $(LIBRARIES)

test_spectra: $(TOOLS) $(SOURCE) $(EXTERNAL) $(TEST_SPECTRA)
	$(CC) $(OPTFLAG) $(LDFLAG) -o  $@ $(addprefix build/,$(notdir $^)) $(LIBRARIES)

test_transfer: $(TOOLS) $(SOURCE) $(EXTERNAL) $(TEST_TRANSFER)
	$(CC) $(OPTFLAG) $(LDFLAG) -o  $@ $(addprefix build/,$(notdir $^)) $(LIBRARIES)

test_nonlinear: $(TOOLS) $(SOURCE) $(EXTERNAL) $(TEST_NONLINEAR)
	$(CC) $(OPTFLAG) $(LDFLAG) -o  $@ $(addprefix build/,$(notdir $^)) $(LIBRARIES)

test_perturbations: $(TOOLS) $(SOURCE) $(EXTERNAL) $(TEST_PERTURBATIONS)
	$(CC) $(OPTFLAG) $(LDFLAG) -o  $@ $(addprefix build/,$(notdir $^)) $(LIBRARIES)

test_thermodynamics: $(TOOLS) $(SOURCE) $(EXTERNAL) $(TEST_THERMODYNAMICS)
	$(CC) $(OPTFLAG) $(LDFLAG) -o  $@ $(addprefix build/,$(notdir $^)) $(LIBRARIES)

test_background: $(TOOLS) $(SOURCE) $(EXTERNAL) $(TEST_BACKGROUND)
	$(CC) $(OPTFLAG) $(LDFLAG) -o  $@ $(addprefix build/,$(notdir $^)) $(LIBRARIES)

test_hyperspherical: $(TOOLS) $(TEST_HYPERSPHERICAL)
	$(CC) $(OPTFLAG) $(LDFLAG) -o test_hyperspherical $(addprefix build/,$(notdir $^)) $(LIBRARIES)

python/cclassy.pxd: $(H_ALL)
	$(PYTHON) python/generate_wrapper.py

classy: libclass.a python/classy.pyx python/cclassy.pxd
	rm -f python/classy.cpp
	cd python; export CC=$(CC); export CXX=$(CXX); $(PYTHON) setup.py install || $(PYTHON) setup.py install --user

clean: .base
	rm -rf $(WRKDIR);
	rm -f libclass.a
	rm -f $(MDIR)/python/classy.cpp
	rm -f $(MDIR)/python/cclassy.pxd
	rm -rf $(MDIR)/python/build
