#!/bin/bash

# go to the directory doc/manual/latex and compile the latex files to
# produce the pdf doc
cd ../manual/latex

make

make

cp refman.pdf ../CLASS_MANUAL.pdf
cd ../../input
