CC = gcc
AR = ar rv
CCFLAG = -O3
LDFLAG = -O3

%.o: %.c 
	$(CC) $(CCFLAG) -c $*.c -o $*.o

HYREC_SRC = hyrectools.o helium.o hydrogen.o history.o 
HYREC_EXE = hyrec.o 

clean: 
	rm *.o

hyrec: $(HYREC_SRC) $(HYREC_EXE)
	$(CC) $(LDFLAG) -o hyrec $(HYREC_SRC) $(HYREC_EXE) -lm

libhyrec.a: $(HYREC_SRC)
	$(AR) $@ $(HYREC_SRC)
