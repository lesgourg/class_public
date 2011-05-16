#Some Makefile for CLASS.
#Julien Lesgourgues, 18.04.2010

MDIR := $(shell pwd)
WRKDIR = $(MDIR)/build
#PMCLIB = $(MDIR)/../../../pmclib

.base:
	if ! [ -a $(WRKDIR) ]; then mkdir $(WRKDIR) ; mkdir $(WRKDIR)/lib; fi;
	touch build/.base

vpath %.c source:tools:main:test
vpath %.o build
vpath .base build

CC       = gcc
#CC       = pgcc
AR        = ar rv

#CCFLAG   = -fast -mp -mp=nonuma -mp=allcores -g
#LDFLAG   = -fast -mp -mp=nonuma -mp=allcores -g
#CCFLAG   = -O4 -Wall -pg
#LDFLAG   = -O4 -Wall -pg
#CCFLAG   = -O0 -Wall -ggdb
#LDFLAG   = -O0 -Wall -ggdb
CCFLAG   = -O4 -Wall -fopenmp
LDFLAG   = -O4 -Wall -fopenmp
#LCCFLAG = -O0
#LDFLAG = -O0
#CCFLAG   = -fast -fopenmp -Wall
#LDFLAG   = -fast -fopenmp -Wall
#CCFLAG = -O0 -ggdb -g -Wall
#LDFLAG = -O0 -ggdb -g -Wall
#CCFLAG = -O4 -arch i386 -pg
#LDFLAG = -O4 -arch i386 -pg
#CCFLAG   = -O2 -Wall -g
#LDFLAG   = -O2 -Wall -g

#-L$(PMCLIB)/lib -lerrorio -lreadConf -lgsl -lgslcblas -llua

INCLUDES = ../include
#-I../../../../pmclib/include/pmclib/tools -I../../../../pmclib/include

%.o:  %.c .base
	cd $(WRKDIR);$(CC) $(CCFLAG) -I$(INCLUDES) -c ../$< -o $*.o

#TOOLS = growTable.o dei_rkck.o evolver_rkck.o arrays.o parser.o quadrature.o
TOOLS = growTable.o dei_rkck.o sparse.o evolver_rkck.o  evolver_ndf15.o arrays.o parser.o quadrature.o

INPUT = input.o

PRECISION = precision.o

BACKGROUND = background.o

THERMO = thermodynamics.o

PERTURBATIONS = perturbations.o 

BESSEL = bessel.o

TRANSFER = transfer.o

PRIMORDIAL = primordial.o

SPECTRA = spectra.o

NONLINEAR = trg.o nonlinear.o

LENSING = lensing.o

OUTPUT = output.o

CLASS = class.o

TEST_OPTIMIZE_1D = test_optimize_1D.o

TEST_OPTIMIZE = test_optimize.o

CHI2 = chi2_camb.o

TEST_TIMING = test_timing.o

TEST_LOOPS = test_loops.o

TEST_SPECTRA = test_spectra.o

TEST_LENSING = test_lensing.o

TEST_TRANSFER = test_transfer.o

TEST_BESSEL = test_bessel.o

TEST_PERTURBATIONS = test_perturbations.o

TEST_THERMODYNAMICS = test_thermodynamics.o

TEST_BACKGROUND = test_background.o

TEST_TRG = test_trg.o

TEST_KARIM = test_karim.o

TEST_PBC = test_pbc.o

all: class libclass.a

libclass.a: $(TOOLS) $(INPUT) $(BACKGROUND) $(THERMO) $(PERTURBATIONS) $(BESSEL) $(TRANSFER) $(PRIMORDIAL) $(SPECTRA) $(NONLINEAR) $(LENSING) $(OUTPUT)
	$(AR)  $@ $(addprefix build/,$(notdir $^))

class: $(TOOLS) $(INPUT) $(BACKGROUND) $(THERMO) $(PERTURBATIONS) $(BESSEL) $(TRANSFER) $(PRIMORDIAL) $(SPECTRA) $(NONLINEAR) $(LENSING) $(OUTPUT) $(CLASS)
	$(CC) $(LDFLAG) -o  $@ $(addprefix build/,$(notdir $^)) -lm

test_timing: $(TOOLS) $(INPUT) $(BACKGROUND) $(THERMO) $(PERTURBATIONS) $(BESSEL) $(TRANSFER) $(PRIMORDIAL) $(SPECTRA) $(NONLINEAR) $(LENSING) $(OUTPUT) $(TEST_TIMING)
	$(CC) $(LDFLAG) -o  $@ $(addprefix build/,$(notdir $^)) -lm

test_lensing: $(TOOLS) $(INPUT) $(BACKGROUND) $(THERMO) $(PERTURBATIONS) $(BESSEL) $(TRANSFER) $(PRIMORDIAL) $(SPECTRA) $(NONLINEAR) $(LENSING) $(OUTPUT) $(TEST_LENSING)
	$(CC) $(LDFLAG) -o  $@ $(addprefix build/,$(notdir $^)) -lm

test_pbc: $(TOOLS) $(INPUT) $(BACKGROUND) $(THERMO) $(PERTURBATIONS) $(BESSEL) $(TRANSFER) $(PRIMORDIAL) $(SPECTRA) $(LENSING) $(OUTPUT) $(TEST_PBC)
	$(CC) $(LDFLAG) -o  $@ $(addprefix build/,$(notdir $^)) -lm

test_optimize_1D: $(TOOLS) $(INPUT) $(BACKGROUND) $(THERMO) $(PERTURBATIONS) $(BESSEL) $(TRANSFER) $(PRIMORDIAL) $(SPECTRA) $(LENSING) $(OUTPUT) $(TEST_OPTIMIZE_1D)
	$(CC) $(LDFLAG) -o  $@ $(addprefix build/,$(notdir $^)) -lm

test_optimize: $(TOOLS) $(INPUT) $(BACKGROUND) $(THERMO) $(PERTURBATIONS) $(BESSEL) $(TRANSFER) $(PRIMORDIAL) $(SPECTRA) $(NONLINEAR) $(LENSING) $(OUTPUT) $(TEST_OPTIMIZE)
	$(CC) $(LDFLAG) -o  $@ $(addprefix build/,$(notdir $^)) -lm

chi2: $(TOOLS) $(INPUT) $(BACKGROUND) $(THERMO) $(PERTURBATIONS) $(BESSEL) $(TRANSFER) $(PRIMORDIAL) $(SPECTRA) $(NONLINEAR) $(LENSING) $(OUTPUT) $(CHI2)
	$(CC) $(LDFLAG) -o  $@ $(addprefix build/,$(notdir $^)) -lm

test_loops: $(TOOLS) $(INPUT) $(BACKGROUND) $(THERMO) $(PERTURBATIONS) $(BESSEL) $(TRANSFER) $(PRIMORDIAL) $(SPECTRA) $(NONLINEAR) $(LENSING) $(OUTPUT) $(TEST_LOOPS)
	$(CC) $(LDFLAG) -o  $@ $(addprefix build/,$(notdir $^)) -lm

test_trg: $(TOOLS) $(INPUT) $(BACKGROUND) $(THERMO) $(PERTURBATIONS) $(BESSEL) $(TRANSFER) $(PRIMORDIAL) $(SPECTRA) $(LENSING) $(TRG) $(TEST_TRG)
	$(CC) $(LDFLAG) -o  $@ $(addprefix build/,$(notdir $^)) -lm

test_spectra: $(TOOLS) $(INPUT) $(BACKGROUND) $(THERMO) $(PERTURBATIONS) $(BESSEL) $(TRANSFER) $(PRIMORDIAL) $(SPECTRA) $(TEST_SPECTRA)
	$(CC) $(LDFLAG) -o  $@ $(addprefix build/,$(notdir $^)) -lm

test_transfer: $(TOOLS) $(INPUT) $(BACKGROUND) $(THERMO) $(PERTURBATIONS) $(BESSEL) $(TRANSFER) $(TEST_TRANSFER)
	$(CC) $(LDFLAG) -o  $@ $(addprefix build/,$(notdir $^)) -lm

test_bessel: $(TOOLS) $(INPUT) $(BACKGROUND) $(THERMO) $(PERTURBATIONS) $(BESSEL) $(TEST_BESSEL)
	$(CC) $(LDFLAG) -o  $@ $(addprefix build/,$(notdir $^)) -lm

test_perturbations: $(TOOLS) $(INPUT) $(BACKGROUND) $(THERMO) $(PERTURBATIONS) $(TEST_PERTURBATIONS)
	$(CC) $(LDFLAG) -o  $@ $(addprefix build/,$(notdir $^)) -lm

test_thermodynamics: $(TOOLS) $(INPUT) $(BACKGROUND) $(THERMO) $(TEST_THERMODYNAMICS)
	$(CC) $(LDFLAG) -o  $@ $(addprefix build/,$(notdir $^)) -lm

test_background: $(TOOLS) $(INPUT) $(BACKGROUND) $(TEST_BACKGROUND)
	$(CC) $(LDFLAG) -o  $@ $(addprefix build/,$(notdir $^)) -lm

test_karim: $(TOOLS) $(INPUT) $(BACKGROUND) $(THERMO) $(PERTURBATIONS) $(BESSEL) $(TRANSFER) $(PRIMORDIAL) $(SPECTRA) $(LENSING) $(OUTPUT) $(TEST_KARIM)
	$(CC) $(LDFLAG) -o  $@ $(addprefix build/,$(notdir $^)) -lm

clean: .base
	rm -rf $(WRKDIR);
