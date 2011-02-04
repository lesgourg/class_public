#!/bin/sh
# CPU: a CLASS Plotting Utility

### PLOTTING FUNCTIONS ###

function plot_cl_lin {
local output=$1.lin.plt
local field=2
if [ "$2" != "" ]
then
  local field=$2
fi
( echo "set term wxt enhanced"
echo "set xlabel 'l'"
echo "set ylabel 'l(l+1)C_l / 2{/Symbol p}'"
echo "set key right"
echo "set title 'CMB computed with presicion, etc, date of the run ?"
echo "plot '$PWD/$1' u 1:$field title 'dimensionless C_l' w l"
) > $PWD/$output
echo "File '$output' generated, launching gnuplot now"
gnuplot -persist $PWD/$output
}

function plot_cl_log {
local output=$1.log.plt
local field=2
if [ "$2" != "" ]
then
  local field=$2
fi
local lmax=$(more $PWD/$1 | grep multipoles | sed 's/# i.e. number of multipoles equal to //')
if [ "$3" != "" ]
then
  let lmax=$3
fi
( echo "set term wxt enhanced"
echo "set logscale x"
echo "set xr [:$lmax]"
echo "set xlabel 'l'"
echo "set ylabel 'l(l+1)C_l / 2{/Symbol p}'"
echo "set key left"
echo "set title 'CMB computed with presicion, etc, date of the run ?"
echo "plot '$PWD/$1' u 1:$field title 'dimensionless C_l' w l"
) > $PWD/$output
echo "File '$output' generated, launching gnuplot now"
gnuplot -persist $PWD/$output
}

function plot_cl_loglinear {
local output=$1.log_linear.plt
local field=2
if [ "$2" != "" ]
then
  local field=$2
fi
local lr=1.
local lmax=$(more $PWD/$1 | grep "multipoles" | sed 's/# i.e. number of multipoles equal to //')
if [ "$3" != "" ]
then
  let lmax=$3
fi
( echo "set term wxt enhanced"
echo "set xlabel 'l'"
echo "set ylabel 'l(l+1)C_l / 2{/Symbol p}'"
echo "set key right"
echo "set title 'CMB multipoles computed..."
echo "set xtics ('2' (sqrt(2/$lr)),'100' (sqrt(100/$lr)), '500' (sqrt(500/$lr)), '1000' (sqrt(1000/$lr)), '1500' (sqrt(1500/$lr)), '2000' (sqrt(2000/$lr)),'2500' (sqrt(2500/$lr)), '$lmax' (sqrt($lmax/$lr)))"
echo "set xr [sqrt(2):sqrt($lmax)]"
echo "plot '$PWD/$1' u (sqrt(\$1/$lr)):$field title 'dimensionless C_l' w l"
) > $PWD/$output
echo "File '$output' generated, launching gnuplot now"
gnuplot -persist $PWD/$output
}


function plot_pk {
local output=$1.plt
local field=2
if [ "$2" != "" ]
then
  local field=$2
fi
local z
if [ "$3" != "" ]
then
  let z=$3
else
  let z=$(more $PWD/$1 | grep redshift | sed 's/# Matter power spectrum P(k) at redshift z=//')
fi
( echo "set term wxt enhanced"
echo "set logscale"
echo "set xlabel 'k (h/Mpc)'"
echo "set ylabel 'P_k (Mpc/h)^3'"
echo "set key right"
echo "set title 'Power spectrum at z=$z, run $1"
echo "plot '$PWD/$1' u 1:$field title 'P_k' w l"
) > $PWD/$output
echo "File '$output' generated, launching gnuplot now"
gnuplot -persist $PWD/$output
}

### MERGING FUNCTIONS ###

function merge_cl {
if [ "$1" = "$2" ]
then
  echo "  You are comparing two identical files, merging aborted"
  exit
fi
local temp=${#1}-14
local base_name_1=${1:7:$temp}
local short_base_name_1=$PWD/"output/"$base_name_1"_parameters.ini"
local lmax=$(more $short_base_name_1 | grep l_max_scalars | sed 's/l_max_scalars = //') # Extract the l_max from the parameter file
let temp=${#2}-14
local base_name_2=${2:7:$temp}
local output="output/"$base_name_1"_"$base_name_2"_cl.dat"
paste $1 $2 > $PWD/$output
echo "  File '$output' created: merging successful"
local field="((\$2/\$7-1.)*100.)"
if [ "$3" = "-lin" ] || [ "$3" = "" ]
then
  plot_cl_lin $output $field 
fi
if [ "$3" = "-log" ]
then
  plot_cl_log $output $field $lmax
fi
if [ "$3" = "-ll" ]
then
  plot_cl_loglinear $output $field $lmax
fi
}

function merge_pk {
if [ "$1" = "$2" ]
then
  echo "  You are comparing two identical files, merging aborted"
  exit
fi
local z1=$(more $1 | grep redshift | sed 's/# Matter power spectrum P(k) at redshift z=//')
local z2=$(more $2 | grep redshift | set 's/# Matter power spectrum P(k) at redshift z=//')
if [ "$z1" != "$z2" ]
then
  echo "  You are comparing the spectrum at two different redshifts, proceed ? (y or n)"
  read
  if [ "$REPLY" = "n" ]
  then
    exit
  fi
fi
local temp
if ls $PWD/$1 | grep -q "z"
then
  let temp=${#1}-17
else
  let temp=${#1}-14
fi
local base_name_1=${1:7:$temp}
if ls $PWD/$2 | grep -q "z"
then
  let temp=${#1}-17
else
  let temp=${#1}-14
fi
local base_name_2=${2:7:$temp}
local output="output/"$base_name_1"_"$base_name_2"_pk.dat"
paste $1 $2 > $PWD/$output
echo "  File '$output' created: merging successful"
local field="(\$2/\$3)"
plot_pk $output $field $z1
}

### MISCELLANEOUS ###

function help {
if [ "$1" = "--help" ]
then
  clear
fi
echo "CPU, a CLASS Plotting Utility"
echo ""
echo "Usage: CPU.sh [FILE...] [OPTIONS...]"
echo ""
echo " CPU.sh plots what you want in a pretty gnuplot graph,"
echo " by extracting the information directly from the .dat files."
echo ""
echo " It also generates a gnuplot script to let you modify"
echo " the details to create a more specific graph."
echo ""
echo " This program also supports the comparison between two files,"
echo " just input their two names as an argument, and voila !"
echo ""
echo " Note that you have to be in the main directory of 'class' for the script to run. "
echo " Indicate file name with the full relative path, i.e. 'output/base_zi_pk.dat'"
echo ""
echo ""
echo "    File:"
echo ""
echo " For now, you can specify at most two files."
echo " *If you only specify one, from the main directory of CLASS, i.e. output/test_cl.dat"
echo "  the program simply writes a gnuplot script, and then read it with gnuplot"
echo " *If you specify two files, the program merge the two files if relevant, and then plot the result"
echo " *For a _cl.dat input, by default the program output only the linear display"
echo ""
echo "    Option (without any file):"
echo ""
echo " -h,   --help    : display this help"
echo " -c,   --clean   : clean all .plt previously generated by the program"
echo ""
echo "    Options after a file:"
echo ""
echo " -lin            : select only the linear     plot (for Cl only) "
echo " -log            : select only the log        plot (for Cl only)"
echo " -ll             : select only the log-linear plot (for Cl only)"
echo ""
echo "    Things to add: an option for printing it to a png file with a quick -p"
echo ""
echo "  To remove the annoying use of ./ in front of the script, put a copy in your ~/bin/ !"
echo ""
echo "  If the program is run without argument,"
echo "  it might one day launch a fancy graphical interface. Maybe"
}

### MAIN ###

# Runs in command line mode if the script is called with an argument.
# It is the fastest way to run this,the drawback is that you can specify less options: plots both linear and log for Cl's
# Typically useful when one wants to quickly check the output of a run, can be run in a script on a cluster for instance

# If no argument is provided, runs into the interactive plotting mode, where you can be more specific on the type of display.

if [[ "$1" = "--help" || "$1" = "-h" ]]
then
  help
  exit
elif [[ "$1" = "--clean" || "$1" = "-c" ]]
then
  echo Cleaning...
  if ls $PWD/output | grep -q ".plt" 
  then
    rm $PWD/output/*.plt
    echo "  output directory is now free from any .plt files"
  else 
    echo "  No file to remove: the output folder does not contain any .plt"
  fi
  exit
elif [ "$1" = "" ]
then
  clear
  echo Class Plotting Utility, interactive mode, yet to find something satisfying
  echo "For a detailled help on the usage of this script, try to run CPU.sh --help"
  exit
elif [ "$#" -lt 4 ]
then
  if echo "$1" | grep -q "cl.dat"
  then
    if echo "$2" | grep -q "cl.dat"
    then
      if echo "$3" | grep -q "cl.dat"
      then
	echo "--> Please try to compare only two files at a time"
	exit
      else
	merge_cl $1 $2 $3
	exit
      fi
    else
      if [ "$2" = "-lin" ] || [ "$2" = "" ]
      then
	plot_cl_lin $1
      fi
      if [ "$2" = "-log" ]
      then
	plot_cl_log $1
      fi
      if [ "$2" = "-ll" ]
      then
	plot_cl_loglinear $1
	exit
      fi
    fi
  elif echo "$1" | grep -q "pk.dat"
  then
    if echo "$2" | grep -q "pk.dat"
    then
      merge_pk $1 $2
      exit
    else
      plot_pk $1
      exit
    fi
  else
    echo "--> Please specify a valid file name (try to run CPU.sh --help)"
    exit
  fi
else
  echo Unrecognizied option, try to run CPU.sh --help for detailed instructions
fi

