PYTHON ?= python3
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

# your C compiler:
CC       = gcc
#CC       = icc
#CC       = pgcc
CPP      = g++ --std=c++11 -fpermissive -Wno-write-strings

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

%.opp:  %.c .base $(HEADERFILES)
	cd $(WRKDIR);$(CPP) $(OPTFLAG) $(OMPFLAG) $(CCFLAG) $(INCLUDES) -c ../$< -o $*.opp

TOOLS = growTable.o dei_rkck.o sparse.o evolver_rkck.o  evolver_ndf15.o arrays.opp parser.o quadrature.o hyperspherical.opp common.o trigonometric_integrals.o

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
MISC_FILES = Makefile CPU psd_FD_single.dat myselection.dat myevolution.dat README bbn/sBBN.dat external_Pk/* cpp
PYTHON_FILES = python/classy.pyx python/setup.py python/cclassy.pxd python/test_class.py

all: class libclass.a classy

libclass.a: $(TOOLS) $(SOURCE) $(EXTERNAL)
	$(AR)  $@ $(addprefix build/, $(TOOLS) $(SOURCE) $(EXTERNAL))

class: $(TOOLS) $(SOURCE) $(EXTERNAL) $(OUTPUT) $(CLASS)
	$(CPP) $(OPTFLAG) $(OMPFLAG) $(LDFLAG) -o class $(addprefix build/,$(notdir $^)) -lm

test_loops: $(TOOLS) $(SOURCE) $(EXTERNAL) $(OUTPUT) $(TEST_LOOPS)
	$(CPP) $(OPTFLAG) $(OMPFLAG) $(LDFLAG) -o $@ $(addprefix build/,$(notdir $^)) -lm

test_loops_omp: $(TOOLS) $(SOURCE) $(EXTERNAL) $(OUTPUT) $(TEST_LOOPS_OMP)
	$(CPP) $(OPTFLAG) $(OMPFLAG) $(LDFLAG) -o $@ $(addprefix build/,$(notdir $^)) -lm

test_harmonic: $(TOOLS) $(SOURCE) $(EXTERNAL) $(TEST_HARMONIC)
	$(CPP) $(OPTFLAG) $(OMPFLAG) $(LDFLAG) -o  $@ $(addprefix build/,$(notdir $^)) -lm

test_transfer: $(TOOLS) $(SOURCE) $(EXTERNAL) $(TEST_TRANSFER)
	$(CPP) $(OPTFLAG) $(OMPFLAG) $(LDFLAG) -o  $@ $(addprefix build/,$(notdir $^)) -lm

test_fourier: $(TOOLS) $(SOURCE) $(EXTERNAL) $(TEST_FOURIER)
	$(CPP) $(OPTFLAG) $(OMPFLAG) $(LDFLAG) -o  $@ $(addprefix build/,$(notdir $^)) -lm

test_perturbations: $(TOOLS) $(SOURCE) $(EXTERNAL) $(TEST_PERTURBATIONS)
	$(CPP) $(OPTFLAG) $(OMPFLAG) $(LDFLAG) -o  $@ $(addprefix build/,$(notdir $^)) -lm

test_thermodynamics: $(TOOLS) $(SOURCE) $(EXTERNAL) $(TEST_THERMODYNAMICS)
	$(CPP) $(OPTFLAG) $(OMPFLAG) $(LDFLAG) -o  $@ $(addprefix build/,$(notdir $^)) -lm

test_background: $(TOOLS) $(SOURCE) $(EXTERNAL) $(TEST_BACKGROUND)
	$(CPP) $(OPTFLAG) $(OMPFLAG) $(LDFLAG) -o  $@ $(addprefix build/,$(notdir $^)) -lm

test_hyperspherical: $(TOOLS) $(TEST_HYPERSPHERICAL)
	$(CC) $(OPTFLAG) $(OMPFLAG) $(LDFLAG) -o test_hyperspherical $(addprefix build/,$(notdir $^)) -lm

tar: $(C_ALL) $(C_TEST) $(H_ALL) $(PRE_ALL) $(INI_ALL) $(MISC_FILES) $(HYREC) $(PYTHON_FILES)
	tar czvf class.tar.gz $(C_ALL) $(H_ALL) $(PRE_ALL) $(INI_ALL) $(MISC_FILES) $(HYREC) $(PYTHON_FILES)

classy: libclass.a python/classy.pyx python/cclassy.pxd
	export CC=$(CC); output=$$($(PYTHON) -m pip install . 2>&1); \
    echo "$$output"; \
    if echo "$$output" | grep -q "ERROR: Cannot uninstall"; then \
        site_packages=$$($(PYTHON) -c "import distutils.sysconfig; print(distutils.sysconfig.get_python_lib())" || $(PYTHON) -c "import site; print(site.getsitepackages()[0])") && \
        echo "Cleaning up previous installation in: $$site_packages" && \
        rm -rf $$site_packages/classy* && \
        $(PYTHON) -m pip install .; \
    fi

clean: .base
	rm -rf $(WRKDIR);
	rm -f libclass.a
	rm -f $(MDIR)/python/classy.c
	rm -rf $(MDIR)/python/build
.PHONY: hs-run hs-clean
hs-run: class
	./class hs_parse_test.ini | tee hs_all.log
	@echo "Artifacts:"
	@ls -lh hs_cmb_summary.tsv hs_bao_summary.tsv hs_distances.tsv hs_distances_plus.tsv hs_summary.csv hs_meta.txt || true
hs-clean:
	rm -f hs_*.log hs_*.tsv hs_summary.csv
.PHONY: hs-summary hs-run-all
hs-summary:
	@set -e; \
	if [ -f hs_cmb_summary.tsv ] && [ -f hs_bao_summary.tsv ]; then \
	  { \
	    head -n 1 hs_cmb_summary.tsv; \
	    echo ',,,,'; \
	    head -n 1 hs_bao_summary.tsv; \
	  } | paste -d, - - - > hs_summary.csv; \
	  { \
	    sed -n '2p' hs_cmb_summary.tsv; \
	    echo ',,,,'; \
	    sed -n '2p' hs_bao_summary.tsv | head -n 1; \
	  } | paste -d, - - - >> hs_summary.csv; \
	else \
	  echo "hs-summary: missing hs_cmb_summary.tsv or hs_bao_summary.tsv; run \`make hs-run\` first" >&2; \
	  exit 1; \
	fi
	@echo "Built hs_summary.csv"
hs-run-all: class
	./class hs_parse_test.ini | tee hs_all.log
	$(MAKE) hs-summary
	$(MAKE) hs-meta
	$(MAKE) hs-dist-plus
	$(MAKE) hs-bundle
	@echo "Artifacts:"
	@ls -lh hs_meta.txt || true
	@ls -lh hs_meta.txt || true
	@ls -lh hs_cmb_summary.tsv hs_bao_summary.tsv hs_distances.tsv hs_distances_plus.tsv hs_summary.csv hs_artifacts.zip hs_meta.txt || true
.PHONY: hs-cls hs-cls-tsv hs-cls-all
hs-cls: class
	./class hs_cls.ini | tee hs_cls.log
	@echo "CLASS wrote:" && ls -lh output/*_cl*.dat || true
hs-cls-all: hs-cls hs-cls-tsv hs-cls-dl-tsv
	@zip -9 hs_cls_artifacts.zip hs_cls.ini hs_cls.log hs_cls_*.tsv >/dev/null || true
	@echo "Artifacts:" && ls -lh hs_cls_*.tsv hs_cls_artifacts.zip || true
	@ls -lh hs_meta.txt || true
	@ls -lh hs_meta.txt || true
hs-cls-tsv:
	./scripts/hs_cls_tsv.sh
.PHONY: hs-cls-dl-tsv hs-cls-bundle
hs-cls-dl-tsv:
	@rm -f hs_cls__Dl.tsv
	@set -e; \
	for b in TT EE TE PP; do \
	  in="hs_cls_$${b}.tsv"; out="hs_cls_$${b}_Dl.tsv"; \
	  if [ -f "$$in" ]; then \
	    awk 'BEGIN{pi=4*atan2(1,1); print "ell\tD_ell"} NR>1{ell=$$1; Cl=$$2; print ell"\t"(ell*(ell+1)*Cl/(2*pi))}' "$$in" > "$$out"; \
	  fi; \
	done; \
	echo "Wrote:" && ls -lh hs_cls_*_Dl.tsv 2>/dev/null || true
hs-cls-bundle: hs-cls-tsv hs-cls-dl-tsv
	@zip -9 hs_cls_artifacts.zip hs_cls.ini hs_cls.log hs_cls_*.tsv >/dev/null || true
	@echo "Artifacts:" && ls -lh hs_cls_*.tsv hs_cls_artifacts.zip || true
	@ls -lh hs_meta.txt || true
	@ls -lh hs_meta.txt || true

.PHONY: hs-sha
hs-sha: hs-release
	@f=$$(ls -1t hs_release_*.zip | head -n1); shasum -a 256 "$$f" hs_artifacts.zip hs_cls_artifacts.zip > SHA256SUMS.txt; echo "Wrote SHA256SUMS.txt"
.PHONY: hs-dist-plus
hs-dist-plus:
	./scripts/hs_dist_plus.sh
.PHONY: hs-release
hs-release: hs-run-all hs-cls-all
	@rev=$$(git rev-parse --short HEAD); \
	zip -9 "hs_release_$${rev}.zip" hs_artifacts.zip hs_cls_artifacts.zip hs_meta.txt README-HS.md >/dev/null || true; \
	echo "Built hs_release_$${rev}.zip"
.PHONY: hs-verify
hs-verify:
	@set -e; for f in hs_cmb_summary.tsv hs_bao_summary.tsv hs_distances.tsv hs_distances_plus.tsv hs_summary.csv hs_all.log hs_meta.txt; do \
		[ -s $$f ] || { echo "Missing $$f"; exit 1; }; \
	done
	@grep -Eri "nan|inf" -- hs_*.tsv || echo "No NaNs/Infs ✅"
.PHONY: hs-publish
hs-publish: hs-run-all hs-cls-all hs-verify
	@tag=$$(git rev-parse --short HEAD); ts=$$(date -u +%Y%m%dT%H%M%SZ); \
	name="hs_artifacts_$${tag}_$${ts}.zip"; \
	rm -f "$$name"; \
	zip -9 "$$name" hs_artifacts.zip hs_cls_artifacts.zip hs_meta.txt README-HS.md >/dev/null || true; \
	echo "Wrote $$name"
.PHONY: hs-meta
hs-meta:
	@./scripts/hs_meta.sh
.PHONY: hs-bundle
hs-bundle:
	@rm -f hs_artifacts.zip
	@zip -9 hs_artifacts.zip hs_*.tsv hs_summary.csv hs_all.log hs_meta.txt >/dev/null || true
	@echo "Built hs_artifacts.zip"
