# Script to change the names of CLASS modules (by Nils Schöneberg & Julien Lesgourgues)
#
# Can be used to:
#  - rename module files, module prefixes, module structures, module structure acronyms
#  - undo renaming
#  - clean the generated log and backup files
#
# usage: CLASS_rename.py [-h] --method {rename,undo,clean} [-v | -q]
#
# optional arguments:
#  -h, --help            show this help message and exit
#  --method {rename,undo,clean}   rename / undo renaming / clean
#  -v, --verbose         Increase the verbosity of the program for more detailed output
#  -q, --quiet           Make the program entirely quiet, setting the verbosity to 0.
#                        Also disables the user confirmation, so use it carefully
#
# The actual renaming to be performed has to be set beforehand in the section below.
# Currently this is set for the transformation
# of CLASS v2.10.8 into CLASS v3.0.0 and backwards.

### EDIT ONLY BELOW ###
### EDIT ONLY BELOW ###

module_filename = ["thermodynamics","perturbations","nonlinear","transfer","spectra"]
module_prefix = ["thermodynamics","perturb","nonlinear","transfer","spectra"]
structure_longname = ["thermo","perturbs","nonlinear","transfers","spectra"]
structure_shortname = ["th","pt","nl","tr","sp"]

newmodule_filename = ["thermodynamics","perturbations","fourier","transfer","harmonic"]
newmodule_prefix = ["thermodynamics","perturbations","fourier","transfer","harmonic"]
newstructure_longname = ["thermodynamics","perturbations","fourier","transfer","harmonic"]
newstructure_shortname= ["th","pt","fo","tr","hr"]

# Potential problem: structure short names are just two
# letters. Combinations of the same two letters may
# appear casually. Thus in some sub-cases we first
# check for exceptions.

# to identify these exception, for each short name (e.g. 'nl'), run:
#
# > grep "nl\." */*.c */*.h */*.py */*.pyx */*.pxd */*.ipynb */*.ini *.ini
#
# > grep "\&nl" */*.c */*.h */*.py */*.pyx */*.pxd */*.ipynb */*.ini *.ini
#
# and check whether some of the lines feature an nl that has nothing
# to do with the stucture short name. If yes, write the exception in
# the dictionary below.

exceptions = {"th":[],
              "pt":[],
              "nl":["nl_corr","R_nl"],
              "tr":[],
              "sp":["osp.","resp"]}

prefix_exceptions = {"thermodynamics":[],
                     "perturbations":[],
                     "nonlinear":["nonlinear_method","nonlinear_scale","nonlinear_min_k_max"],
                     "transfer":[],
                     "spectra":[]}

src_folder = "source"
incl_folder = "include"
test_folder = "test"

### EDIT ONLY ABOVE ###
### EDIT ONLY ABOVE ###




import os
import argparse

# parse the arguments of the command line
parser = argparse.ArgumentParser(description='Change the names of CLASS modules')
parser.add_argument("--method",choices=["rename","undo","clean"], required=True,help="Do you want to rename / undo renaming / clean the generated log and backup files? Type 'rename','undo', or 'clean'")
# default verbosity is 1, can be increased with -v or decreased with -q
group = parser.add_mutually_exclusive_group()
group.add_argument("-v", "--verbose", action="count",default=1,help="Increase the verbosity of the program for more detailed output")
group.add_argument("-q", "--quiet", action="store_true",help="Make the program entirely quiet, setting the verbosity to 0. This also disables the confirmation required by the user, so use it carefully")
parse_dict = parser.parse_args()

if parse_dict.quiet:
  parse_dict.verbose = 0

# Inform the user about the starting of the actual routine
if parse_dict.verbose>0:
  print("START RENAMING ROUTINE v.0.3 (credits Nils Schöneberg & Julien Lesgourgues)")
  print("CHECKING ALL FILES IN DIRECTORY : "+os.path.abspath("."))

# Find the list of all the directories that we will parse and in which we will do changes in some files
#
# Get the list of all folders and subfolders in the local folder
# After this step, each element x is such that x[0] contains a folder name 'folder/subfolder/.../'
folder_list = [x for x in os.walk(".")]
# remove .git, doc, build folders
folder_list = [x for x in folder_list if not (".git" in x[0])]
#folder_list = [x for x in folder_list if not ("doc" in x[0])]
folder_list = [x for x in folder_list if not ("doc/manual" in x[0])]
folder_list = [x for x in folder_list if not ("doc/input/latex" in x[0])]
folder_list = [x for x in folder_list if not ("build" in x[0])]
# remove the folder of the RealSpaceInterface containing cached data
folder_list = [x for x in folder_list if not ("RealSpaceInterface/static" in x[0])]
# keep only the list of all folders, not the files they contain
folder_list = [x[0] for x in folder_list]

if parse_dict.verbose > 0:
  # show the list of 'folder/subfolder/.../'
  print("FOLDER LIST : "+" ".join(folder_list))

  # let the user confirm or abort
  read = input("Continue? (y/n)")
  if not read.startswith("y"):
    quit()

###############
# 'undo' mode #
###############

if parse_dict.method == "undo":

  # For each changed module/file name, go back to old file names (e.g. 'fourier.c' -> 'nonlinear.c') such that they can be overwritten by the corresponding .old files
  for i in range(len(module_filename)):
    xf = module_filename[i]
    yf = newmodule_filename[i]
    os.rename(os.path.join(src_folder,yf+".c"),os.path.join(src_folder,xf+".c"))
    os.rename(os.path.join(incl_folder,yf+".h"),os.path.join(incl_folder,xf+".h"))
    os.rename(os.path.join(test_folder,"test_"+yf+".c"),os.path.join(test_folder,"test_"+xf+".c"))
    if parse_dict.verbose > 0:
      print("REVERTED TO MODULE NAME "+xf)

  # find all folders containing .old and/or .unchanged files
  for fldername in folder_list:
    # First, get name of all files and subfolders in this folder
    filelist_all = os.listdir(fldername)
    filelist = []
    for fname in filelist_all:
      tmp_name = os.path.join(fldername, fname)
      # remove subfolder names, keep only file names
      if os.path.isdir(tmp_name):
        continue
      if tmp_name.endswith(".old"):
        filelist.append(fname)
      elif tmp_name.endswith(".unchanged"):
        filelist.append(fname)

    if parse_dict.verbose > 2:
      print (fldername, filelist)

    for filename in filelist:
      # remove the log files *.unchanged
      # (this can be done safely, as all relevant information is in the .old files)
      if(".unchanged" in filename):
        os.remove(os.path.join(fldername,filename))

      # remove the .old extensions, thus overwriting the changed files with the old files
      if(".old" in filename):
        os.rename(os.path.join(fldername,filename),os.path.join(fldername,filename.replace(".old","")))
        if parse_dict.verbose > 1:
          print ("mv "+os.path.join(fldername,filename)+" "+os.path.join(fldername,filename.replace(".old",""))+"!")

    if parse_dict.verbose > 0:
      print ("IN "+fname+", DELETED .unchanged AND RESTORED ORIGINAL FROM .old FILES")


################
# 'clean' mode #
################

elif parse_dict.method == "clean":

  # find all folders containing .old and/or .unchanged files
  for fldername in folder_list:
    # First, get name of all files and subfolders in this folder
    filelist_all = os.listdir(fldername)
    filelist = []
    for fname in filelist_all:
      tmp_name = os.path.join(fldername, fname)
      # remove subfolder names, keep only file names
      if os.path.isdir(tmp_name):
        continue
      if tmp_name.endswith(".old"):
        filelist.append(fname)
      elif tmp_name.endswith(".unchanged"):
        filelist.append(fname)

    if parse_dict.verbose > 2:
      print (fldername,filelist)

    for filename in filelist:
      # just remove any .unchanged or .old files
      if(".unchanged" in filename or ".old" in filename):
        os.remove(os.path.join(fldername,filename))

    if parse_dict.verbose > 0:
      print ("IN "+fname+", DELETED .unchanged AND .old FILES")

    try:
      # remove log files Makefile.old and possibly autostep.py
      os.remove("Makefile.old")
      if parse_dict.verbose > 0:
        print ("REMOVED Makefile.old")
      os.remove(os.path.join("python","autosetup.py"))
      if parse_dict.verbose > 0:
        print ("REMOVED python/autosetup.py")
    except:
      pass

#################
# 'rename' mode #
#################

elif parse_dict.method == "rename":

  # Some operations only have to be done for the first iteration over all files
  # One example of this is the generation of the backup .old files
  # Thus, we keep track if this is our first iteration
  first_loop = True

  # loop over each module to be renamed/modified
  for i in range(len(module_filename)):

    xf = module_filename[i]
    xp = module_prefix[i]
    xsl = structure_longname[i]
    xss = structure_shortname[i]

    yf = newmodule_filename[i]
    yp = newmodule_prefix[i]
    ysl = newstructure_longname[i]
    yss = newstructure_shortname[i]

    if parse_dict.verbose > 0:
      print("BEGIN RENAMING {} -> {}".format(xsl,ysl))

    # Parse and possibly do changes in each file of the folder fldername
    for fldername in folder_list:

      # Establish the list of file to be parsed and possibly modified in this folder

      # First, get name of all files and subfolders in this folder
      filelist_all = os.listdir(fldername)
      filelist = []
      for fname in filelist_all:
        tmp_name = os.path.join(fldername, fname)
        # remove subfolder names, keep only file names
        if os.path.isdir(tmp_name):
          continue
        # ignore the automatically generated python setup file
        if "autosetup.py" in tmp_name:
          continue
        # take into account other files with extension .c, .py, .pyx, .pxd, .ipynb, .h, .ini
        # but not the .py of the local (root) folder (and thus e.g. not this script!)
        if tmp_name.endswith(".c"):
          filelist.append(fname)
        elif tmp_name.endswith(".py"):
          if fldername != '.':
            filelist.append(fname)
        elif tmp_name.endswith(".pyx"):
          filelist.append(fname)
        elif tmp_name.endswith(".pxd"):
          filelist.append(fname)
        elif tmp_name.endswith(".ipynb"):
          filelist.append(fname)
        elif tmp_name.endswith(".h"):
          filelist.append(fname)
        elif tmp_name.endswith(".ini"):
          filelist.append(fname)
        elif tmp_name.endswith(".md"):
          filelist.append(fname)

      # show the list of file to be parsed and possibly modified in this folder
      if parse_dict.verbose > 1:
        print("WILL MODIFY ALL FILES IN FOLDER '{}': [".format(fldername)+",".join(filelist)+"]")

      # iterate over all files in the current folder
      for filename in filelist:
        # open input file (with old names)
        with open(os.path.join(fldername,filename),"r") as inf:
          # open temporary output file (where we will subsititue the new names)
          with open(os.path.join(fldername,filename+".tmp"),"w") as outf:
            # open a log file with extension .unchanged where we will store lines that were not changed but should have, potentially (for visual inspection)
            with open(os.path.join(fldername,filename+".unchanged"),"a") as unchf:
              # iterate over each line in the input file
              line = inf.readline()
              while line:

                # I. Treat lines where the full structure name appears, e.g. 'nonlinear'
                if "struct "+xsl in line:
                  if "struct "+xsl+" "+xss in line:
                    # replace each structure declaration (e.g. 'struct nonlinear nl' -> 'struct fourier fo')
                    # we isolate this case because it is very useful to catch many occurences of the structure short name (e.g. 'nl') already here
                    line = line.replace("struct "+xsl+" "+xss,"struct "+ysl+" "+yss)
                  else:
                    # replace other occurences (e.g. 'struct nonlinear' -> 'struct fourier')
                    # Special care is needed here! Check that the next character is not a letter
                    # Thus we only allow for a small selection of relevant possibilities
                    for char in ['\t','\n',' ','*','`',';',':']:
                      line = line.replace("struct "+xsl+char,"struct "+ysl+char)
                if "cdef "+xsl in line:
                  if "cdef "+xsl+" "+xss in line:
                    # replace each structure declaration (e.g. 'cdef nonlinear nl' -> 'cdef fourier fo')
                    # we isolate this case because it is very useful to catch many occurences of the structure short name (e.g. 'nl') already here
                    line = line.replace("cdef "+xsl+" "+xss,"cdef "+ysl+" "+yss)
                  else:
                    # replace other occurences (e.g. 'cdef nonlinear' -> 'cdef fourier')
                    line = line.replace("cdef "+xsl,"cdef "+ysl)
                if xsl+" structure" in line:
                  line = line.replace(xsl+" structure",ysl+" structure")

                # II. Treat lines where the module (= file) name appears
                if xf.upper() in line:
                  # replace capitalized module name (e.g. '__NONLINEAR__' -> '__FOURIER__')
                  line = line.replace(xf.upper(),yf.upper())
                if xf+".c" in line:
                  # replace full filename in the comments (e.g. 'nonlinear.c' --> 'fourier.c')
                  line = line.replace(xf+".c",yf+".c")
                if xf+".h" in line:
                  # replace full filename in the comments (e.g. 'nonlinear.h' --> 'fourier.h')
                  line = line.replace(xf+".h",yf+".h")
                if xf+" module" in line:
                  # replace full filename in the comments (e.g. 'nonlinear module' --> 'fourier module')
                  line = line.replace(xf+" module",yf+" module")
                if "\""+xf+"\"" in line:
                  # replace full filename in quotation marks (e.g. '"nonlinear"' --> '"fourier"')
                  line = line.replace("\""+xf+"\"","\""+yf+"\"")


                # III. Treat lines where the prefix appears
                if xp+"_" in line:
                   # replace all prefix names (e.g. 'nonlinear' -> 'fourier')
                   # For all prefix exceptions, substitute the problematic string with 'xx'
                   for i,x in enumerate(prefix_exceptions[xp]):
                     if x in line:
                       line = line.replace(x,prefix_exceptions[xp][i].replace(xp,'xx'))
                   # Now replace all the corresponding names where the prefix appears
                   line = line.replace(xp+"_",yp+"_")
                   # Finally, re-substitute the original exception string instead of the 'xx'
                   for i,x in enumerate(prefix_exceptions[xp]):
                     if x.replace(xp,'xx') in line:
                       line = line.replace(x.replace(xp,'xx'),prefix_exceptions[xp][i])

                # IV. Treat line where short structure name appears, e.g. 'nl'
                if xss in line:

                  # replace pointers towards structure (e.g. 'pnl' -> 'pfo')
                  if "p"+xss in line:
                    line = line.replace("p"+xss,"p"+yss)

                  # replace structure addresses (e.g. '&nl' -> '&fo') and structure short names before dots (e.g. 'nl.error_message' -> 'fo.error_message')
                  if "&"+xss in line or  xss+"." in line:
                    # For all exceptions, substitute the problematic string with 'xx'
                    for i,x in enumerate(exceptions[xss]):
                      if x in line:
                        line = line.replace(x,exceptions[xss][i].replace(xss,'xx'))
                    # Now replace all structure short names before dots and addresses
                    line = line.replace("&"+xss,"&"+yss)
                    line = line.replace(xss+".",yss+".")
                    # Finally, re-substitute the original exception string instead of the 'xx'
                    for i,x in enumerate(exceptions[xss]):
                      if x.replace(xss,'xx') in line:
                        line = line.replace(x.replace(xss,'xx'),exceptions[xss][i])

                  # replace structures as fields of bigger structures in python (e.g. 'self.nl' -> 'self.fo')
                  if "self."+xss in line:
                    line = line.replace("self."+xss,"self."+yss)

                  # if the line did contain the short name in another circumstances, print it in the log file .unchanged
                  if xss in line:
                    # Mark the occurence of the short name by arrows (e.g. 'only' -> 'o-->nl<--y')
                    unchf.write(line.replace(xss,"-->"+xss+"<--"))

                # write the line (changed or not) in the temporary output file
                outf.write(line)
                line = inf.readline()

        # keep the input file but add to it an extension .old (so we keep it as a backup, if something goes wrong) e.g. nonlinear.c -> nonlinear.c.old
        # This is done only in the first loop over modules.
        if first_loop == True:
          os.rename(os.path.join(fldername,filename),os.path.join(fldername,filename+".old"))
        # give to the temporary output file name its final extension (e.g. 'nonlinear.c.tmp' -> 'nonlinear.c')
        os.rename(os.path.join(fldername,filename+".tmp"),os.path.join(fldername,filename))

      if parse_dict.verbose > 1:
        print("SUCCESS IN FOLDER {}".format(fldername))

    # work on the Makefile
    if parse_dict.verbose>1:
      print("MODIFY MAKEFILE")
    with open("Makefile","r") as inf:
      # implement the changes in Makefile.tmp
      with open("Makefile.tmp","w") as outf:
        line = inf.readline()
        while line:
          # replace long names (e.g. 'nonlinear' -> 'fourier')
          if xf in line:
            line = line.replace(xf,yf)
          # replace long names when capitalized
          if xf.upper() in line:
            line = line.replace(xf.upper(),yf.upper())
          outf.write(line)
          line = inf.readline()
    # keep old version with additional .old extension
    if first_loop == True:
      os.rename("Makefile","Makefile.old")
    # rename Makefile.tmp -> Makefile
    os.rename("Makefile.tmp","Makefile")
    if parse_dict.verbose>1:
      print("SUCCESS IN MODIFYING MAKEFILE")

    # change actual file names (e.g. 'nonlinear.c' -> 'fourier.c')
    if parse_dict.verbose>1:
      print("RENAME MODULE "+yf)
    os.rename(os.path.join(src_folder,xf+".c"),os.path.join(src_folder,yf+".c"))
    os.rename(os.path.join(incl_folder,xf+".h"),os.path.join(incl_folder,yf+".h"))
    os.rename(os.path.join(test_folder,"test_"+xf+".c"),os.path.join(test_folder,"test_"+yf+".c"))
    if parse_dict.verbose>1:
      print("SUCCESS IN RENAMING MODULE "+yf)


    if parse_dict.verbose > 0:
      print("SUCCESS FOR RENAMING {} -> {}".format(xf,yf))

    # done for this particular module
    first_loop = False

# end of loop over modulea
if parse_dict.verbose>0:
  print("SUCCESS!")
