Updating the manual  
===================
Author: D. C. Hooper (Hooper@physik.rwth-aachen.de)

This pdf manual and accompanying web version have been generated using the doxygen software (http://www.doxygen.org). This software directly reads the code and extracts the necessary comments to form the manual, meaning it is very easy to generate newer versions of the manual as desired.

To keep the manual up-to-date, a new version should be generated after any major upgrade to `CLASS`. To keep track of how updated the manual is, the title page also displays the last modification date.

To generate a new version of this manual, one should install the doxygen software. Once installed, doxygen uses a specific configuration file to know how to read the code. The configuration file for this project can be found in /class/doc/input/doxyconf. To run doxygen, navigate in terminal to the above-mentioned folder containing the configuration file and type

`doxygen doxyconf`

This will generate a new version of the html manual and the necessary files to make the pdf version. Note that any changes in the `doxyconf` file can dramatically impact the outcome, so the configuration file should only be modified with great care.

Currently doxygen does not offer the option to order the output chapters in the pdf version of the manual. Hence, before compiling the pdf one must check that the manual is ordered correctly. To do this, navigate to /class/doc/output/latex. From here, the `refman.tex` file can be easily modified to obtain the desired order. Once the `refman.tex file is correct, the pdf can be created in the same directory by typing

`make`

in the terminal. This will result in the generation of the pdf manual. It is often useful to run `make` twice consecutively, to insure all the references and links have been generated correctly. The updated version of the manual should now be ready. For convenience, one can copy the final pdf to /class/doc/output.

