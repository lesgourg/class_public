"""
.. module:: CPU
    :synopsis: CPU, a CLASS Plotting Utility
.. moduleauthor:: Benjamin Audren <benj_audren@yahoo.fr>
.. version:: 2.0

This is a small python program aimed to gain time when comparing two spectra,
i.e. from CAMB and CLASS, or a non-linear spectrum to a linear one.  It is
designed to be used in a command line fashion, not being restricted to your
CLASS directory, though it recognized mainly CLASS output format.  Far from
perfect, or complete, it could use any suggestion for enhancing it, just to
avoid losing time on useless matters for others.  Be warned that, when
comparing with other format, the following is assumed: there are no empty line
(especially at the end of file). Gnuplot comment lines (starting with a # are
allowed). This issue will cause a non-very descriptive error in CPU, any
suggestion for testing it is welcome.  Example of use: To superimpose two
different spectra and see their global shape :
python CPU.py output/lcdm_z2_pk.dat output/lncdm_z2_pk.dat
To see in details their ratio:
python CPU.py output/lcdm_z2_pk.dat output/lncdm_z2_pk.dat -r

"""
import numpy as np
import os
import matplotlib.pyplot as plt
import sys
import argparse
import textwrap

START_LINE = {}
START_LINE['error'] = [r' /|\   ',
                       r'/_o_\  ',
                       r'       ']
START_LINE['warning'] = [r' /!\ ',
                         r'     ']
START_LINE['info'] = [r' /!\ ',
                      r'     ']

STANDARD_LENGTH = 80  # standard, increase if you have a big screen


def create_parser():
    parser = argparse.ArgumentParser(
        description=(
            'CPU, a CLASS Plotting Utility, specify wether you want\n'
            'to superimpose, or plot the ratio of different files.'),
        epilog=(
            'A standard usage would be, for instance:\n'
            'python CPU.py output/test_pk.dat output/test_pk_nl_density.dat'
            ' -r\npython CPU.py output/wmap_cl.dat output/planck_cl.dat'),
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument(
        'files', type=str, nargs='*', help='Files to plot')
    parser.add_argument('-r', '--ratio', dest='ratio', action='store_true',
                        help='Plot the ratio of the spectra')
    parser.add_argument('-s', '--selection', dest='selection',
                        nargs='+',
                        help='specify the fields you want to plot.')
    parser.add_argument('--scale', choices=['lin', 'loglog', 'loglin'],
                        type=str,
                        help='Specify the scale to use for the plot')
    parser.add_argument(
        '-p, --print',
        dest='printfile', action='store_true', default=False,
        help='print the graph directly in a .png file')
    parser.add_argument(
        '-r, --repeat',
        dest='repeat', action='store_true', default=False,
        help='repeat the step for all redshifts with same base name')
    return parser


def plot_CLASS_output(files, selection, ratio=False, output_name='',
                      extension='', x_variable='', scale=''):
    """
    Load the data to numpy arrays, write a Python script and plot them.

    Inspired heavily by the matlab version by Thomas Tram

    Parameters
    ----------
    files : list
        List of files to plot
    selection : list, or string
        List of items to plot, which should match the way they appear in the
        file, for instance: ['TT', 'BB]

    Keyword Arguments
    -----------------
    ratio : bool
        If set to yes, plots the ratio of the files, taking as a reference the
        first one
    output_name : str
        Specify a different name for the produced figure (by default, it takes
        the name of the first file, and replace the .dat by .pdf)
    extension : str


    """
    # Load all the graphs
    data = []
    for data_file in files:
        data.append(np.loadtxt(data_file))

    # Create the python script, and initialise it
    python_script_path = files[0]+'.py'
    text = '''
import matplotlib.pyplot as plt
import numpy as np\n'''

    # Create the full_path_files list, that contains the absolute path, so that
    # the future python script can import them directly.
    full_path_files = [os.path.abspath(elem) for elem in files]
    text += '''files = %s\n''' % full_path_files
    text += '''
data = []
for data_file in files:
    data.append(np.loadtxt(data_file))
selection = %s"\n''' % selection

    # Recover the number of columns in the first file, as well as their title.
    with open(files[0], 'r') as header_file:
        header = [line for line in header_file if line[0] == '#']
        header = header[-1]

    # Count the number of columns in the file, and recover their name. Thanks
    # Thomas Tram for the trick
    indices = [i+1 for i in range(len(header)) if
               header.startswith(':', i)]
    num_columns = len(indices)
    long_names = [header[indices[i]:indices[(i+1)]-2].strip() if i < num_columns-1
             else header[indices[i]:].strip()
             for i in range(num_columns)]


    # Process long_names further to handle special cases, and extract names,
    # which will correspond to the tags specified in "selection".
    names = []
    for index, name in enumerate(long_names):
        # This can happen in the background file
        if name.startswith('(.)', 0):
            temp_name = name[3:]
            long_names[index] = replace_scale(name)
        # Otherwise, we simply
        else:
            names.append(name)
            pass

    print long_names
    if num_columns == 2:
        selection = None
    elif num_columns > 2:
        # in case selection was only a string, cast it to a list
        if isinstance(selection, str):
            selection = [selection]
        for elem in selection:
            print elem
            if elem not in names:
                raise InputError(
                    "Expected in selection a name of a field in the"
                    " specified files.")
    # Create the figure and ax objects
    fig, ax = plt.subplots()
    text += '\nfig, ax = plt.subplots()\n'

    # if ratio is not set, then simply plot them all
    if not ratio:
        text += 'for curve in data:'
        for curve in data:
            text += 'ax.'
            #if spectrum_type == 'cl_lin':
                #text += 'plot'
                #ax.plot(curve[:, 0], curve[:, colnum])
            #elif spectrum_type == 'cl_log':
                #text += 'loglog'
                #ax.loglog(curve[:, 0], curve[:, colnum])
            #elif spectrum_type == 'pk':
                #text += 'loglog'
                #ax.loglog(curve[:, 0], curve[:, colnum])
            #text += '(curve[:, 0], curve[:, colnum])\n'
    else:
        ref = data[0]
        #for index in range(1, len(data)):
            #current = data[index]
            #if np.allclose(current[0], ref[0]):
                #ax.plot(current[0], current[colnum]/ref[colnum])
    text += 'plt.show()\n'
    plt.show()

    # Write to the python file all the issued commands. You can then reproduce
    # the plot by running "python output/something_cl.dat.py"
    with open(python_script_path, 'w') as python_script:
        python_script.write(text)


# Errors ###########################
class MyError(Exception):
    """
    Base class defining the general presentation of error messages
    """
    def __init__(self, message):
        """Reformat the name of the class for easier reading"""
        Exception.__init__(self)
        self.message = message
        name = self.__class__.__name__
        self.name = ''
        # Extract the name, and add spaces between the capital letters
        for index, letter in enumerate(name):
            if letter.isupper():
                if index > 0:
                    self.name += ' ' + letter
                else:
                    self.name += letter
            else:
                self.name += letter

    def __str__(self):
        """Define the behaviour under the print statement"""
        return '\n\n' + self.name + ':' + pretty_print(
            self.message, "error", True)


def pretty_print(string, status, return_string=False):
    """
    Return the string formatted according to its status

    The input is a potentially long message, describing the problem.
    According to the severity of its status (so far, 'error' will exit the
    program, whereas 'warning' and 'info' will go through anyway).

    Standard length has been defined globally, as well as the ascii-art
    dictionary of arrays START_LINE.

    """

    if return_string:
        output = ''
    length = STANDARD_LENGTH-len(START_LINE[status][0])
    # Remove unwanted spaces (coming from carriage returns in the input string)
    # and handle voluntary carriage returns specified with \n
    first_cleanup = [' '.join(elem.lstrip(' ').split())
                     for elem in string.split('\n')]
    splitted = []
    # Recover the lines splitted at correct length
    for elem in first_cleanup:
        splitted.extend(textwrap.wrap(elem, length))

    if status == 'error':
        # Add a blank line so that the error displays better
        if return_string:
            output += '\n'
        else:
            print

    # Add in front the appropriate fancy display
    index = 0
    for line in splitted:
        # If the number of needed lines is bigger than the ascii-art, the last
        # line of it (empty) will be used.
        if index < len(START_LINE[status]):
            start_index = index
        else:
            start_index = len(START_LINE[status])-1
        if return_string:
            output += START_LINE[status][start_index]+line+'\n'
        else:
            print START_LINE[status][start_index]+line
        index += 1
    if return_string:
        return output
    else:
        return


class FormatError(MyError):
    """Format not recognised"""
    pass


class TypeError(MyError):
    """Spectrum type not recognised"""
    pass


class NumberOfFilesError(MyError):
    """Invalid number of files"""
    pass


class InputError(MyError):
    """Incompatible input requirements"""
    pass


def replace_scale(string):
    """
    This assumes that the string starts with "(.)", which will be replaced by
    (8piG/3)

    >>> print replace_scale('(.)toto')
    >>> '(8\\pi G/3)toto'
    """
    string_list = list(string)
    string_list.pop(1)
    string_list[1:1] = list('8\\pi G/3')
    return ''.join(string_list)


def main():
    print '~~~ Running CPU, a CLASS Plotting Utility ~~~'
    parser = create_parser()
    # Parse the command line arguments
    args = parser.parse_args()

    # if there are no argument in the input, print usage
    if len(args.files) == 0:
        parser.print_usage()
        return

    # if the first file name contains cl or pk, infer the type of desired
    # spectrum
    if not args.selection:
        if args.files[0].rfind('cl') != -1:
            selection = 'TT'
            scale = 'loglog'
        elif args.files[0].rfind('pk') != -1:
            selection = 'P'
            scale = 'loglog'
        else:
            raise TypeError(
                "Please specify a field to plot")
        if not args.scale:
            args.scale = scale
        args.selection = selection

    # If ratio is asked, but only one file was passed in argument, politely
    # complain
    if args.ratio:
        if len(args.files) < 2:
            raise NumberOfFilesError(
                "If you want me to compute a ratio between two files, "
                "I strongly encourage you to give me at least two of them.")
    # actual plotting. By default, a simple superposition of the graph is
    # performed. If asked to be divided, the ratio is shown - whether a need
    # for interpolation arises or not.
    plot_CLASS_output(args.files, args.selection,
                      ratio=args.ratio, scale=args.scale)

if __name__ == '__main__':
    sys.exit(main())
