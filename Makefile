#Some Makefile for CLASS.
#Julien Lesgourgues, 28.11.2011
#Nils Schöneberg, Matteo Lucca, 27.02.2019

MDIR := $(shell pwd)
WRKDIR = $(MDIR)/build

.base:
	if ! [ -e $(WRKDIR) ]; then mkdir $(WRKDIR) ; mkdir $(WRKDIR)/lib; fi;
	touch build/.base

vpath %.c source:tools:main:test
vpath %.o build
vpath %.opp build
vpath .base build

########################################################
###### LINES TO ADAPT TO YOUR PLATEFORM ################
########################################################

# Your C and C++ compilers.
#
# CC and CXX are chosen, in decreasing order of priority, from:
#   1. the make command line:   make CC=gcc-14 CXX=g++-14
#   2. the environment:         export CC=... CXX=...
#      (conda's clang_osx-*/gcc_linux-* packages set these when the
#      environment is activated, and pip passes the environment on to make)
#   3. the defaults chosen below ('cc' and 'c++').
# Always set CC and CXX together, from the same compiler family: the C++
# standard library used by libclass.a and by the python extension must match
# (libc++ with clang, libstdc++ with gcc), otherwise linking fails.
#
# CHANGED: CC used to be hardcoded as 'gcc' ("CC = gcc"), and the C++ compiler
# as "CPP = g++ ...". A plain assignment in a Makefile overrides the
# environment, so a user's CC/CXX (and conda's) were silently ignored and
# /usr/bin/gcc and /usr/bin/g++ were always used, even inside a conda
# environment with its own compilers. Only the make command line could change
# them, and that is not reachable through 'pip install'.
#
# CC: make has a built-in default (CC = cc), which the environment and the
# command line already override, so CC is deliberately NOT assigned here.
# We default to the generic name 'cc' rather than 'gcc' because:
#   - on macOS 'gcc' is always Apple clang (cc is too, so nothing is lost);
#   - in a conda environment created with the 'compilers' metapackage
#     'gcc' is a real GNU gcc while 'cc' and 'c++' are conda's clang, so
#     'cc' and 'c++' are the pair that belongs together;
#   - on Linux 'cc' is the system compiler, normally gcc.
# Other compilers (icc, pgcc, ...) can be selected with CC=icc, CXX=icpc, etc.
#
# CXX: make's built-in default is 'g++', which would pair a conda 'cc' (clang)
# with a system 'g++' (see above). Override that default only if the user has
# not set CXX ('origin' is 'default' only when neither the environment nor
# the command line defined it).
ifeq ($(origin CXX),default)
CXX = c++
endif

# C++ language flags for the C++ compiler.
# CHANGED: these were previously part of the compiler variable ("CPP = g++
# --std=c++11 ..."). They are split out so that a user-supplied CXX (which
# is a bare compiler name, e.g. conda's arm64-apple-darwin20.0.0-clang++)
# still gets them. They are used as "$(CXX) $(CXXSTD)" below.
#   --std=c++11          include/parallel.h needs C++11 (std::future etc.)
#   -fpermissive         CLASS is C-style code compiled as C++ (see .opp)
#   -Wno-write-strings   string literals are assigned to char *
CXXSTD   = --std=c++11 -fpermissive -Wno-write-strings

# CHANGED: the variable called CPP was removed. In make's and the
# environment's convention CPP is the C *preprocessor* (conda sets
# CPP=...-clang-cpp, make's built-in is "$(CC) -E"); the C++ compiler is CXX.
# Using CPP for the C++ compiler made an environment-provided CPP unusable
# and could not be combined with 'make -e'.

# your tool for creating static libraries:
AR        = ar rv

# Your python interpreter.
# In order to use Python 3, you can manually
# substitute python3 to python in the line below, or you can simply
# add a compilation option on the terminal command line:
# "PYTHON=python3 make all" (Thanks to Marius Millea for python3 compatibility)
PYTHON ?= python

# your optimization flag
OPTFLAG = -O3
#OPTFLAG = -Ofast -ffast-math #-march=native
#OPTFLAG = -fast

# your openmp flag (comment for compiling without openmp)
OMPFLAG   = -pthread #-fopenmp
#OMPFLAG   = -mp -mp=nonuma -mp=allcores -g
#OMPFLAG   = -openmp

# all other compilation flags
CCFLAG = -g -fPIC
LDFLAG = -g -fPIC

# leave blank to compile without HyRec, or put path to HyRec directory
# (with no slash at the end: e.g. "external/RecfastCLASS")
HYREC = external/HyRec2020
RECFAST = external/RecfastCLASS
HEATING = external/heating
HALOFIT = external/Halofit
HMCODE = external/HMcode

########################################################
###### IN PRINCIPLE THE REST SHOULD BE LEFT UNCHANGED ##
########################################################

# pass current working directory to the code
CLASSDIR ?= $(MDIR)
CCFLAG += -D__CLASSDIR__='"$(CLASSDIR)"'

# where to find include files *.h
INCLUDES = -I../include
HEADERFILES = $(wildcard ./include/*.h)

# automatically add external programs if needed. First, initialize to blank.
EXTERNAL =

vpath %.c $(RECFAST)
#CCFLAG += -DRECFAST
INCLUDES += -I../$(RECFAST)
EXTERNAL += wrap_recfast.o
HEADERFILES += $(wildcard ./$(RECFAST)/*.h)

vpath %.c $(HEATING)
#CCFLAG += -DHEATING
INCLUDES += -I../$(HEATING)
EXTERNAL += injection.o noninjection.o
HEADERFILES += $(wildcard ./$(HEATING)/*.h)

# update flags for including HyRec
ifneq ($(HYREC),)
vpath %.c $(HYREC)
CCFLAG += -DHYREC
#LDFLAGS += -DHYREC
INCLUDES += -I../$(HYREC)
EXTERNAL += hyrectools.o helium.o hydrogen.o history.o wrap_hyrec.o energy_injection.o
HEADERFILES += $(wildcard ./$(HYREC)/*.h)
endif

vpath %.c $(HALOFIT)
INCLUDES += -I../$(HALOFIT)
EXTERNAL += halofit.o
HEADERFILES += $(wildcard ./$(HALOFIT)/*.h)

vpath %.c $(HMCODE)
INCLUDES += -I../$(HMCODE)
EXTERNAL += hmcode.opp
HEADERFILES += $(wildcard ./$(HMCODE)/*.h)

%.o:  %.c .base $(HEADERFILES)
	cd $(WRKDIR);$(CC) $(OPTFLAG) $(OMPFLAG) $(CCFLAG) $(INCLUDES) -c ../$< -o $*.o

# The .opp rule compiles the C sources that include parallel.h as C++ (the
# object files are named .opp so that they are built by this rule rather than
# by the C rule above; the list of which files are .opp is in TOOLS, SOURCE
# and EXTERNAL below).
# CHANGED: "-x c++" was added before the source file, and $(CPP) became
# "$(CXX) $(CXXSTD)". The sources have a .c extension; g++ silently treats
# them as C++ in C++ mode, but clang++ does so with a deprecation warning
# ("treating 'c' input as 'c++' when in C++ mode, this behavior is
# deprecated"), which is expected to become an error. "-x c++" states the
# language explicitly for every compiler.
%.opp:  %.c .base $(HEADERFILES)
	cd $(WRKDIR);$(CXX) $(CXXSTD) $(OPTFLAG) $(OMPFLAG) $(CCFLAG) $(INCLUDES) -c -x c++ ../$< -o $*.opp

TOOLS = growTable.o dei_rkck.o sparse.o evolver_rkck.o evolver_ndf15.o evolver_rkdp45.o arrays.opp parser.o quadrature.o hyperspherical.opp common.o trigonometric_integrals.o

SOURCE = input.o background.o thermodynamics.o perturbations.opp primordial.opp fourier.o transfer.opp harmonic.opp lensing.opp distortions.o

INPUT = input.o

PRECISION = precision.o

BACKGROUND = background.o

THERMO = thermodynamics.o

PERTURBATIONS = perturbations.o

TRANSFER = transfer.o

PRIMORDIAL = primordial.o

HARMONIC = harmonic.o

FOURIER = fourier.o

LENSING = lensing.o

DISTORTIONS = distortions.o

OUTPUT = output.o

CLASS = class.o

TEST_LOOPS = test_loops.o

TEST_LOOPS_OMP = test_loops_omp.o

TEST_HARMONIC = test_harmonic.o

TEST_TRANSFER = test_transfer.o

TEST_FOURIER = test_fourier.o

TEST_PERTURBATIONS = test_perturbations.o

TEST_THERMODYNAMICS = test_thermodynamics.o

TEST_BACKGROUND = test_background.o

TEST_HYPERSPHERICAL = test_hyperspherical.o

C_TOOLS =  $(addprefix tools/, $(addsuffix .c,$(basename $(TOOLS))))
C_SOURCE = $(addprefix source/, $(addsuffix .c,$(basename $(SOURCE) $(OUTPUT))))
C_TEST = $(addprefix test/, $(addsuffix .c,$(basename $(TEST_DEGENERACY) $(TEST_LOOPS) $(TEST_TRANSFER) $(TEST_FOURIER) $(TEST_PERTURBATIONS) $(TEST_THERMODYNAMICS))))
C_MAIN = $(addprefix main/, $(addsuffix .c,$(basename $(CLASS))))
C_ALL = $(C_MAIN) $(C_TOOLS) $(C_SOURCE)
H_ALL = $(addprefix include/, common.h svnversion.h $(addsuffix .h, $(basename $(notdir $(C_ALL)))))
PRE_ALL = cl_ref.pre clt_permille.pre
INI_ALL = explanatory.ini lcdm.ini
MISC_FILES = Makefile CPU psd_FD_single.dat myselection.dat myselection_bins.dat myevolution.dat README bbn/sBBN.dat external_Pk/* cpp
PYTHON_FILES = python/classy.pyx python/setup.py python/cclassy.pxd python/test_class.py

all: class libclass.a classy

libclass.a: $(TOOLS) $(SOURCE) $(EXTERNAL)
	$(AR)  $@ $(addprefix build/, $(TOOLS) $(SOURCE) $(EXTERNAL))

# CHANGED (this rule and the test_* rules below): the executables are linked
# with "$(CXX) $(CXXSTD)" instead of "$(CPP)", because the .opp objects need
# the C++ runtime library. test_hyperspherical has no C++ objects and is
# linked with $(CC).
class: $(TOOLS) $(SOURCE) $(EXTERNAL) $(OUTPUT) $(CLASS)
	$(CXX) $(CXXSTD) $(OPTFLAG) $(OMPFLAG) $(LDFLAG) -o class $(addprefix build/,$(notdir $^)) -lm

test_loops: $(TOOLS) $(SOURCE) $(EXTERNAL) $(OUTPUT) $(TEST_LOOPS)
	$(CXX) $(CXXSTD) $(OPTFLAG) $(OMPFLAG) $(LDFLAG) -o $@ $(addprefix build/,$(notdir $^)) -lm

test_loops_omp: $(TOOLS) $(SOURCE) $(EXTERNAL) $(OUTPUT) $(TEST_LOOPS_OMP)
	$(CXX) $(CXXSTD) $(OPTFLAG) $(OMPFLAG) $(LDFLAG) -o $@ $(addprefix build/,$(notdir $^)) -lm

test_harmonic: $(TOOLS) $(SOURCE) $(EXTERNAL) $(TEST_HARMONIC)
	$(CXX) $(CXXSTD) $(OPTFLAG) $(OMPFLAG) $(LDFLAG) -o  $@ $(addprefix build/,$(notdir $^)) -lm

test_transfer: $(TOOLS) $(SOURCE) $(EXTERNAL) $(TEST_TRANSFER)
	$(CXX) $(CXXSTD) $(OPTFLAG) $(OMPFLAG) $(LDFLAG) -o  $@ $(addprefix build/,$(notdir $^)) -lm

test_fourier: $(TOOLS) $(SOURCE) $(EXTERNAL) $(TEST_FOURIER)
	$(CXX) $(CXXSTD) $(OPTFLAG) $(OMPFLAG) $(LDFLAG) -o  $@ $(addprefix build/,$(notdir $^)) -lm

test_perturbations: $(TOOLS) $(SOURCE) $(EXTERNAL) $(TEST_PERTURBATIONS)
	$(CXX) $(CXXSTD) $(OPTFLAG) $(OMPFLAG) $(LDFLAG) -o  $@ $(addprefix build/,$(notdir $^)) -lm

test_thermodynamics: $(TOOLS) $(SOURCE) $(EXTERNAL) $(TEST_THERMODYNAMICS)
	$(CXX) $(CXXSTD) $(OPTFLAG) $(OMPFLAG) $(LDFLAG) -o  $@ $(addprefix build/,$(notdir $^)) -lm

test_background: $(TOOLS) $(SOURCE) $(EXTERNAL) $(TEST_BACKGROUND)
	$(CXX) $(CXXSTD) $(OPTFLAG) $(OMPFLAG) $(LDFLAG) -o  $@ $(addprefix build/,$(notdir $^)) -lm

test_hyperspherical: $(TOOLS) $(TEST_HYPERSPHERICAL)
	$(CC) $(OPTFLAG) $(OMPFLAG) $(LDFLAG) -o test_hyperspherical $(addprefix build/,$(notdir $^)) -lm

tar: $(C_ALL) $(C_TEST) $(H_ALL) $(PRE_ALL) $(INI_ALL) $(MISC_FILES) $(HYREC) $(PYTHON_FILES)
	tar czvf class.tar.gz $(C_ALL) $(H_ALL) $(PRE_ALL) $(INI_ALL) $(MISC_FILES) $(HYREC) $(PYTHON_FILES)

# CHANGED: "export CC=$(CC)" became 'export CC="$(CC)" CXX="$(CXX)"'.
# The python build (setup.py) runs 'make libclass.a' again and compiles the
# python extension, and both must see the same compilers as this Makefile.
# Exporting the values explicitly passes on the choice made here even when it
# came from a default (make does not export its built-in defaults) or from
# the make command line. CXX is new; it was not passed on before. The values
# are quoted because a compiler variable may contain spaces (e.g. "ccache cc").
classy: libclass.a python/classy.pyx python/cclassy.pxd
	export CC="$(CC)" CXX="$(CXX)"; output=$$($(PYTHON) -m pip install . 2>&1); \
    echo "$$output"; \
    if echo "$$output" | grep -q "ERROR: Cannot uninstall"; then \
        site_packages=$$($(PYTHON) -c "import distutils.sysconfig; print(distutils.sysconfig.get_python_lib())" || $(PYTHON) -c "import site; print(site.getsitepackages()[0])") && \
        echo "Cleaning up previous installation in: $$site_packages" && \
        rm -rf $$site_packages/classy $$site_packages/classy.* $$site_packages/classy-* && \
        $(PYTHON) -m pip install .; \
    fi

clean: .base
	rm -rf $(WRKDIR);
	rm -f libclass.a
	rm -f $(MDIR)/python/classy.c
	rm -rf $(MDIR)/python/build
