#!/bin/bash

cd ../manual/latex

make

make

cp refman.pdf ../CLASS_MANUAL.pdf
cd ../../input
