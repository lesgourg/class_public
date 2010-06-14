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

CCFLAG   = -O2 -fopenmp
LDFLAG   = -O2 -fopenmp
#CCFLAG = -O2
#LDFLAG = -O2 

#-L$(PMCLIB)/lib -lerrorio -lreadConf -lgsl -lgslcblas -llua

INCLUDES = ../include
#-I../../../../pmclib/include/pmclib/tools -I../../../../pmclib/include

%.o:  %.c .base
	cd $(WRKDIR);$(CC) $(CCFLAG) -I$(INCLUDES) -c ../$< -o $*.o

TOOLS = growTable.o dei_rkck.o tools_arrays.o parser.o

INPUT = input.o

PRECISION = precision.o

BACKGROUND = background.o

THERMO = thermodynamics.o

PERTURBATIONS = perturbations.o 

BESSEL = bessel.o

TRANSFER = transfer.o

PRIMORDIAL = primordial.o

SPECTRA = spectra.o

OUTPUT = output.o

CLASS = class.o

TRG = trg.o

TEST_TRANSFER = test_transfer.o

TEST_BESSEL = test_bessel.o

TEST_PERTURBATIONS = test_perturbations.o

TEST_THERMODYNAMICS = test_thermodynamics.o

TEST_BACKGROUND = test_background.o

TEST_TRG = test_trg.o

class: $(TOOLS) $(INPUT) $(PRECISION) $(BACKGROUND) $(THERMO) $(PERTURBATIONS) $(BESSEL) $(TRANSFER) $(PRIMORDIAL) $(SPECTRA) $(OUTPUT) $(CLASS)
	$(CC) $(LDFLAG) -o  $@ $(addprefix build/,$(notdir $^)) -lm

test_trg: $(TOOLS) $(PRECISION) $(INPUT) $(BACKGROUND) $(THERMO) $(PERTURBATIONS) $(BESSEL) $(TRANSFER) $(PRIMORDIAL) $(SPECTRA) $(TRG) $(TEST_TRG)
	$(CC) $(LDFLAG) -o  $@ $(addprefix build/,$(notdir $^)) -lm

test_transfer: $(TOOLS) $(PRECISION) $(INPUT) $(BACKGROUND) $(THERMO) $(PERTURBATIONS) $(BESSEL) $(TRANSFER) $(TEST_TRANSFER)
	$(CC) $(LDFLAG) -o  $@ $(addprefix build/,$(notdir $^)) -lm

test_bessel: $(TOOLS) $(PRECISION) $(INPUT) $(BACKGROUND) $(THERMO) $(PERTURBATIONS) $(BESSEL) $(TEST_BESSEL)
	$(CC) $(LDFLAG) -o  $@ $(addprefix build/,$(notdir $^)) -lm

test_perturbations: $(TOOLS) $(PRECISION) $(INPUT) $(BACKGROUND) $(THERMO) $(PERTURBATIONS) $(TEST_PERTURBATIONS)
	$(CC) $(LDFLAG) -o  $@ $(addprefix build/,$(notdir $^)) -lm

test_thermodynamics: $(TOOLS) $(PRECISION) $(INPUT) $(BACKGROUND) $(THERMO) $(TEST_THERMODYNAMICS)
	$(CC) $(LDFLAG) -o  $@ $(addprefix build/,$(notdir $^)) -lm

test_background: $(TOOLS) $(PRECISION) $(INPUT) $(BACKGROUND) $(TEST_BACKGROUND)
	$(CC) $(LDFLAG) -o  $@ $(addprefix build/,$(notdir $^)) -lm

clean: .base
	rm -rf $(WRKDIR);
