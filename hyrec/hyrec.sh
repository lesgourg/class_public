rm *.o
gcc -03 -c arrays.c -o arrays.o
gcc -03 -c helium.c -o helium.o
gcc -03 -c hydrogen.c -o hydrogen.o
gcc -03 -c history.c -o history.o
gcc -03 -c hyrec.c -o hyrec.o
gcc -03 arrays.o helium.o hydrogen.o history.o hyrec.o -o hyrec
ar rv libhyrec.a arrays.o helium.o hydrogen.o history.o
