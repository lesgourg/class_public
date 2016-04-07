Updating the manual  
===================
Author: D. C. Hooper (hooper@physik.rwth-aachen.de)

This pdf manual and accompanying web version have been generated using the `doxygen` software (http://www.doxygen.org). This software directly reads the code and extracts the necessary comments to form the manual, meaning it is very easy to generate newer versions of the manual as desired.

To maintain the usefulness of the manual, a new version should be generated after any major upgrade to `CLASS`. To keep track of how up-to-date the manual is the title page also displays the last modification date.

Generating a new version of this manual is straightforward. First, you need to install the `doxygen` software, which can be done by following the instructions on the software's webpage. The location where you install this software is irrelevant; it doesn't need to be in the same folder as `CLASS`.

Once installed, navigate to the class/doc/input directory and run the first script 

` . make1.sh`

This will generate a new version of the html manual and the necessary files to make the pdf version. Unfortunately, `doxygen` does not yet offer the option to automatically order the output chapters in the pdf version of the manual. Hence, before compiling the pdf, this must be done manually. To do this you need to find the `refman.tex` file in class/doc/manual/latex. With this file you can modify the title page, headers, footers, and chapter ordering for the final pdf. Once you have this file with your desired configuration, navigate back to the class/doc/input directory, and run the second script 

` . make2.sh`

You should now be able to find the finished pdf in the class/doc/manual/CLASS_MANUAL.pdf.

As a final comment, doxygen uses two main configuration files: `doxyconf` and `doxygen.sty`, both located in class/doc/input. Changes to these files can dramatically impact the outcome, so any modifications to these files should be done with great care.
