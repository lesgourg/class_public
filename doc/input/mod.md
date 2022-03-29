Updating the manual
===================
Author: D. C. Hooper (hooper@physik.rwth-aachen.de)

This pdf manual and accompanying web version have been generated using the `doxygen` software (http://www.doxygen.org). This software directly reads the code and extracts the necessary comments to form the manual, meaning it is very easy to generate newer versions of the manual as desired.

### For CLASS developpers: ###

To maintain the usefulness of the manual, a new version should be generated after any major upgrade to `CLASS`. To keep track of how up-to-date the manual is the title page also displays the last modification date. The manual is generated	automatically from the code, excepted a few chapters written manually in the files

    README.md
    doc/input/chap2.md
    doc/input/chap3.md
    doc/input/mod.md
    external_Pk/README.md

You can	update these files, or add new ones that should	be declared in the `INPUT=` field of `doc/input/doxyconf`.

Generating a new version of this manual is straightforward. First, you need to install the `doxygen` software, which can be done by following the instructions on the software's webpage. The location where you install this software is irrelevant; it doesn't need to be in the same folder as `CLASS`. For Mac OSX, homebrew users can install the software with `brew install doxygen --with-graphviz`.

Once installed, navigate to the class/doc/input directory and run the first script

` . make1.sh`

This will generate a new version of the html manual and the necessary files to make the pdf version. Unfortunately, `doxygen` does not yet offer the option to automatically order the output chapters in the pdf version of the manual. Hence, before compiling the pdf, this must be done manually. To do this you need to find the `refman.tex` file in class/doc/manual/latex. With this file you can modify the title page, headers, footers, and chapter ordering for the final pdf. Usually we just make two things: add manually the line

    \vspace*{1cm}
    {\large Last updated \today}\\

after

    {\Large CLASS MANUAL }\\

and move manually the chapters `"The external Pk mode"` and `"Updating the manual"` to the end, after the automatically generated part. Once you have this file with your desired configuration, navigate back to the class/doc/input directory, and run the second script

` . make2.sh`

You should now be able to find the finished pdf in `class/doc/manual/CLASS_MANUAL.pdf`. Finally you can commit the changes to git, but not all the content of `doc/` is necessary: only `doc/README`, `doc/input/` and `doc/manual/CLASS_MANUAL.pdf`. Since version 2.8, we are not committing anymore `doc/manual/html/` because it was too big (and complicating the version history): users only get the PDF manual from git.

As a final comment, doxygen uses two main configuration files: `doxyconf` and `doxygen.sty`, both located in class/doc/input. Changes to these files can dramatically impact the outcome, so any modifications to these files should be done with great care.
