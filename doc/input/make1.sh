#!/bin/bash

# Let doxygen parse the code, with settings described in doxyconf.
# This will produce the html doc, as well as the .tex files necessary for the
# pdf doc.
doxygen doxyconf

# The repository contains a doxygen.sty file generated automatically
# with the command 'doxygen -w latex header.tex footer.tex doxygen.sty
# doxyconf'. As long as the latex compilation works, this file can
# just be copied to the directory where the latex compilation is
# done. However, a new doxygen.sty should be generated after each
# change of doxyconf or with each new version of doxygen.
cp doxygen.sty ../manual/latex/doxygen.sty
