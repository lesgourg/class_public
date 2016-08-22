#!/usr/bin/env python
"""
.. module:: CPU
    :synopsis: CPU, a CLASS Plotting Utility
.. moduleauthor:: Benjamin Audren <benjamin.audren@gmail.com>
.. credits:: Benjamin Audren, Jesus Torrado
.. version:: 2.0

This is a small python program aimed to gain time when comparing two spectra,
e.g. from CAMB and CLASS, or a non-linear spectrum to a linear one.

It is designed to be used in a command line fashion, not being restricted to
your CLASS directory, though it recognizes mainly CLASS output format. Far from
perfect, or complete, it could use any suggestion for enhancing it,
just to avoid losing time on useless matters for others.

Be warned that, when comparing with other format, the following is assumed:
there are no empty line (especially at the end of file). Gnuplot comment lines
(starting with a # are allowed). This issue will cause a non-very descriptive
error in CPU, any suggestion for testing it is welcome.

Example of use:
- To superimpose two different spectra and see their global shape :
python CPU.py output/lcdm_z2_pk.dat output/lncdm_z2_pk.dat
- To see in details their ratio:
python CPU.py output/lcdm_z2_pk.dat output/lncdm_z2_pk.dat -r

The "PlanckScale" is taken with permission from Jesus Torrado's:
cosmo_mini_toolbox, available under GPLv3 at
https://github.com/JesusTorrado/cosmo_mini_toolbox

"""

from __future__ import unicode_literals, print_function

# System imports
import os
import sys
import argparse

# Numerics
import numpy as np
from numpy import ma
from scipy.interpolate import InterpolatedUnivariateSpline
from math import floor

# Plotting
import matplotlib.pyplot as plt
from matplotlib import scale as mscale
from matplotlib.transforms import Transform
from matplotlib.ticker import FixedLocator


def CPU_parser():
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
    parser.add_argument('-y', '--y-axis', dest='y_axis', nargs='+',
                        help='specify the fields you want to plot.')
    parser.add_argument('-x', '--x-axis', dest='x_axis', type=str,
                        help='specify the field to be used on the x-axis')
    parser.add_argument('--scale', type=str,
                        choices=['lin', 'loglog', 'loglin', 'george'],
                        help='Specify the scale to use for the plot')
    parser.add_argument('--xlim', dest='xlim', nargs='+', type=float,
                        default=[], help='Specify the x range')
    parser.add_argument('--ylim', dest='ylim', nargs='+', type=float,
                        default=[], help='Specify the y range')
    parser.add_argument(
        '-p, --print',
        dest='printfile', default='',
        help=('print the graph directly in a file. If no name is specified, it'
              'uses the name of the first input file'))
    parser.add_argument(
        '--repeat',
        dest='repeat', action='store_true', default=False,
        help='repeat the step for all redshifts with same base name')
    return parser


def plot_CLASS_output(files, x_axis, y_axis, ratio=False, printing='',
                      output_name='', extension='', x_variable='',
                      scale='lin', xlim=[], ylim=[]):
    """
    Load the data to numpy arrays, write all the commands for plotting to a
    Python script for further refinment, and display them.

    Inspired heavily by the matlab version by Thomas Tram

    Parameters
    ----------
    files : list
        List of files to plot
    x-axis : string
        name of the column to use as the x coordinate
    y-axis : list, str
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
    # Define the python script name, and the pdf path
    python_script_path = os.path.splitext(files[0])[0]+'.py'

    # The variable text will contain all the lines to be printed in the end to
    # the python script path, joined with newline characters. Beware of the
    # indentation.
    text = ['import matplotlib.pyplot as plt',
            'import numpy as np',
            'import itertools', '']

    # Load all the graphs
    data = []
    for data_file in files:
        data.append(np.loadtxt(data_file))

    # Create the full_path_files list, that contains the absolute path, so that
    # the future python script can import them directly.
    full_path_files = [os.path.abspath(elem) for elem in files]

    text += ['files = %s' % full_path_files]
    text += ['data = []',
             'for data_file in files:',
             '    data.append(np.loadtxt(data_file))']

    # Recover the base name of the files, everything before the dot
    roots = [elem.split(os.path.sep)[-1].split('.')[0] for elem in files]
    text += ['roots = [%s]' % ', '.join(["'%s'" % root for root in roots])]

    # Create the figure and ax objects
    fig, ax = plt.subplots()
    text += ['', 'fig, ax = plt.subplots()']

    # if ratio is not set, then simply plot them all
    original_y_axis = y_axis
    legend = []
    if not ratio:
        for index, curve in enumerate(data):
            # Recover the number of columns in the first file, as well as their
            # title.
            num_columns, names, tex_names = extract_headers(files[index])

            text += ['', 'index, curve = %i, data[%i]' % (index, index)]
            # Check if everything is in order
            if num_columns == 2:
                y_axis = [names[1]]
            elif num_columns > 2:
                # in case y_axis was only a string, cast it to a list
                if isinstance(original_y_axis, str):
                    y_axis = [original_y_axis]
                else:
                    y_axis = original_y_axis

            # Store the selected text and tex_names to the script
            selected = []
            for elem in y_axis:
                selected.extend(
                    [name for name in names if name.find(elem) != -1 and
                     name not in selected])
            if not y_axis:
                selected = names[1:]
            y_axis = selected

            # Decide for the x_axis, by default the index will be set to zero
            x_index = 0
            if x_axis:
                for index_name, name in enumerate(names):
                    if name.find(x_axis) != -1:
                        x_index = index_name
                        break
            # Store to text
            text += ['y_axis = %s' % selected]
            text += ['tex_names = %s' % [elem for (elem, name) in
                     zip(tex_names, names) if name in selected]]
            text += ["x_axis = '%s'" % tex_names[x_index]]
            text += ["ylim = %s" % ylim]
            text += ["xlim = %s" % xlim]

            for selec in y_axis:
                index_selec = names.index(selec)
                plot_line = 'ax.'
                if scale == 'lin':
                    plot_line += 'plot(curve[:, %i], curve[:, %i])' % (
                        x_index, index_selec)
                    ax.plot(curve[:, x_index], curve[:, index_selec])
                elif scale == 'loglog':
                    plot_line += 'loglog(curve[:, %i], abs(curve[:, %i]))' % (
                        x_index, index_selec)
                    ax.loglog(curve[:, x_index], abs(curve[:, index_selec]))
                elif scale == 'loglin':
                    plot_line += 'semilogx(curve[:, %i], curve[:, %i])' % (
                        x_index, index_selec)
                    ax.semilogx(curve[:, x_index], curve[:, index_selec])
                elif scale == 'george':
                    plot_line += 'plot(curve[:, %i], curve[:, %i])' % (
                        x_index, index_selec)
                    ax.plot(curve[:, x_index], curve[:, index_selec])
                    ax.set_xscale('planck')
                text += [plot_line]

            legend.extend([roots[index]+': '+elem for elem in y_axis])

        ax.legend(legend, loc='best')
        text += ["",
                 "ax.legend([root+': '+elem for (root, elem) in",
                 "    itertools.product(roots, y_axis)], loc='best')",
                 ""]
    else:
        ref = data[0]
        num_columns, ref_curve_names, ref_tex_names = extract_headers(files[0])
        # Check if everything is in order
        if num_columns == 2:
            y_axis_ref = [ref_curve_names[1]]
        elif num_columns > 2:
            # in case y_axis was only a string, cast it to a list
            if isinstance(original_y_axis, str):
                y_axis_ref = [original_y_axis]
            else:
                y_axis_ref = original_y_axis

        # Store the selected text and tex_names to the script
        selected = []
        for elem in y_axis_ref:
            selected.extend([name for name in ref_curve_names if name.find(elem) != -1 and
                             name not in selected])
        y_axis_ref = selected

        # Decide for the x_axis, by default the index will be set to zero
        x_index_ref = 0
        if x_axis:
            for index_name, name in enumerate(ref_curve_names):
                if name.find(x_axis) != -1:
                    x_index_ref = index_name
                    break

        for idx in range(1, len(data)):
            current = data[idx]
            num_columns, names, tex_names = extract_headers(files[idx])

            # Check if everything is in order
            if num_columns == 2:
                y_axis = [names[1]]
            elif num_columns > 2:
                # in case y_axis was only a string, cast it to a list
                if isinstance(original_y_axis, str):
                    y_axis = [original_y_axis]
                else:
                    y_axis = original_y_axis

            # Store the selected text and tex_names to the script
            selected = []
            for elem in y_axis:
                selected.extend([name for name in names if name.find(elem) != -1 and
                                 name not in selected])
            y_axis = selected

            text += ['y_axis = %s' % selected]
            text += ['tex_names = %s' % [elem for (elem, name) in
                                         zip(tex_names, names) if name in selected]]

            # Decide for the x_axis, by default the index will be set to zero
            x_index = 0
            if x_axis:
                for index_name, name in enumerate(names):
                    if name.find(x_axis) != -1:
                        x_index = index_name
                        break

            text += ["x_axis = '%s'" % tex_names[x_index]]
            for selec in y_axis:
                # Do the interpolation
                axis = ref[:, x_index_ref]
                reference = ref[:, ref_curve_names.index(selec)]
                #plt.loglog(current[:, x_index], current[:, names.index(selec)])
                #plt.show()
                #interpolated = splrep(current[:, x_index],
                                      #current[:, names.index(selec)])
                interpolated = InterpolatedUnivariateSpline(current[:, x_index],
                                      current[:, names.index(selec)])
                if scale == 'lin':
                    #ax.plot(axis, splev(ref[:, x_index_ref],
                                        #interpolated)/reference-1)
                    ax.plot(axis, interpolated(ref[:, x_index_ref])/reference-1)
                elif scale == 'loglin':
                    #ax.semilogx(axis, splev(ref[:, x_index_ref],
                                            #interpolated)/reference-1)
                    ax.semilogx(axis, interpolated(ref[:, x_index_ref])/reference-1)
                elif scale == 'loglog':
                    raise InputError(
                        "loglog plot is not available for ratios")

    if 'TT' in names:
        ax.set_xlabel('$\ell$', fontsize=16)
        text += ["ax.set_xlabel('$\ell$', fontsize=16)"]
    elif 'P' in names:
        ax.set_xlabel('$k$ [$h$/Mpc]', fontsize=16)
        text += ["ax.set_xlabel('$k$ [$h$/Mpc]', fontsize=16)"]
    else:
        ax.set_xlabel(tex_names[x_index], fontsize=16)
        text += ["ax.set_xlabel('%s', fontsize=16)" % tex_names[x_index]]
    if xlim:
        if len(xlim) > 1:
            ax.set_xlim(xlim)
            text += ["ax.set_xlim(xlim)"]
        else:
            ax.set_xlim(xlim[0])
            text += ["ax.set_xlim(xlim[0])"]
        ax.set_ylim()
        text += ["ax.set_ylim()"]
    if ylim:
        if len(ylim) > 1:
            ax.set_ylim(ylim)
            text += ["ax.set_ylim(ylim)"]
        else:
            ax.set_ylim(ylim[0])
            text += ["ax.set_ylim(ylim[0])"]
    text += ['plt.show()']
    plt.show()

    # If the use wants to print the figure to a file
    if printing:
        fig.savefig(printing)
        text += ["fig.savefig('%s')" % printing]

    # Write to the python file all the issued commands. You can then reproduce
    # the plot by running "python output/something_cl.dat.py"
    with open(python_script_path, 'w') as python_script:
        print('Creating a python script to reproduce the figure')
        print('--> stored in %s' % python_script_path)
        python_script.write('\n'.join(text))

    # If the use wants to print the figure to a file
    if printing:
        fig.savefig(printing)


class FormatError(Exception):
    """Format not recognised"""
    pass


class TypeError(Exception):
    """Spectrum type not recognised"""
    pass


class NumberOfFilesError(Exception):
    """Invalid number of files"""
    pass


class InputError(Exception):
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


def process_long_names(long_names):
    """
    Given the names extracted from the header, return two arrays, one with the
    short version, and one tex version

    >>> names, tex_names = process_long_names(['(.)toto', 'proper time [Gyr]'])
    >>> print names
    >>> ['toto', 'proper time']
    >>> print tex_names
    >>> ['(8\\pi G/3)toto, 'proper time [Gyr]']

    """
    names = []
    tex_names = []
    # First pass, to remove the leading scales
    for name in long_names:
        # This can happen in the background file
        if name.startswith('(.)', 0):
            temp_name = name[3:]
            names.append(temp_name)
            tex_names.append(replace_scale(name))
        # Otherwise, we simply
        else:
            names.append(name)
            tex_names.append(name)

    # Finally, remove any extra spacing
    names = [''.join(elem.split()) for elem in names]
    return names, tex_names


def extract_headers(header_path):
    with open(header_path, 'r') as header_file:
        header = [line for line in header_file if line[0] == '#']
        header = header[-1]

    # Count the number of columns in the file, and recover their name. Thanks
    # Thomas Tram for the trick
    indices = [i+1 for i in range(len(header)) if
               header.startswith(':', i)]
    num_columns = len(indices)
    long_names = [header[indices[i]:indices[(i+1)]-3].strip() if i < num_columns-1
                  else header[indices[i]:].strip()
                  for i in range(num_columns)]

    # Process long_names further to handle special cases, and extract names,
    # which will correspond to the tags specified in "y_axis".
    names, tex_names = process_long_names(long_names)

    return num_columns, names, tex_names


def main():
    print('~~~ Running CPU, a CLASS Plotting Utility ~~~')
    parser = CPU_parser()
    # Parse the command line arguments
    args = parser.parse_args()

    # if there are no argument in the input, print usage
    if len(args.files) == 0:
        parser.print_usage()
        return

    # if the first file name contains cl or pk, infer the type of desired
    # spectrum
    if not args.y_axis:
        if args.files[0].rfind('cl') != -1:
            scale = 'loglog'
        elif args.files[0].rfind('pk') != -1:
            scale = 'loglog'
        else:
            scale = 'lin'
        args.y_axis = []
    else:
        scale = ''
    if not args.scale:
        if scale:
            args.scale = scale
        else:
            args.scale = 'lin'

    # Remove extra spacing in the y_axis list
    args.y_axis = [''.join(elem.split()) for elem in args.y_axis]
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
    if args.ratio and args.scale == 'loglog':
        print("Defaulting to loglin scale")
        args.scale = 'loglin'

    plot_CLASS_output(args.files, args.x_axis, args.y_axis,
                      ratio=args.ratio, printing=args.printfile,
                      scale=args.scale, xlim=args.xlim, ylim=args.ylim)


# Helper code from cosmo_mini_toolbox, by Jesus Torrado, available fully at
# https://github.com/JesusTorrado/cosmo_mini_toolbox, to use the log then
# linear scale for the multipole axis when plotting Cl.
nonpos = "mask"
change = 50.0
factor = 500.


def _mask_nonpos(a):
    """
    Return a Numpy masked array where all non-positive 1 are
    masked. If there are no non-positive, the original array
    is returned.
    """
    mask = a <= 0.0
    if mask.any():
        return ma.MaskedArray(a, mask=mask)
    return a


def _clip_smaller_than_one(a):
    a[a <= 0.0] = 1e-300
    return a


class PlanckScale(mscale.ScaleBase):
    """
    Scale used by the Planck collaboration to plot Temperature power spectra:
    base-10 logarithmic up to l=50, and linear from there on.

    Care is taken so non-positive values are not plotted.
    """
    name = 'planck'

    def __init__(self, axis, **kwargs):
        pass

    def set_default_locators_and_formatters(self, axis):
        axis.set_major_locator(
            FixedLocator(
                np.concatenate((np.array([2, 10, change]),
                                np.arange(500, 2500, 500)))))
        axis.set_minor_locator(
            FixedLocator(
                np.concatenate((np.arange(2, 10),
                                np.arange(10, 50, 10),
                                np.arange(floor(change/100), 2500, 100)))))

    def get_transform(self):
        """
        Return a :class:`~matplotlib.transforms.Transform` instance
        appropriate for the given logarithm base.
        """
        return self.PlanckTransform(nonpos)

    def limit_range_for_scale(self, vmin, vmax, minpos):
        """
        Limit the domain to positive values.
        """
        return (vmin <= 0.0 and minpos or vmin,
                vmax <= 0.0 and minpos or vmax)

    class PlanckTransform(Transform):
        input_dims = 1
        output_dims = 1
        is_separable = True
        has_inverse = True

        def __init__(self, nonpos):
            Transform.__init__(self)
            if nonpos == 'mask':
                self._handle_nonpos = _mask_nonpos
            else:
                self._handle_nonpos = _clip_nonpos

        def transform_non_affine(self, a):
            lower = a[np.where(a<=change)]
            greater = a[np.where(a> change)]
            if lower.size:
                lower = self._handle_nonpos(lower * 10.0)/10.0
                if isinstance(lower, ma.MaskedArray):
                    lower = ma.log10(lower)
                else:
                    lower = np.log10(lower)
                lower = factor*lower
            if greater.size:
                greater = (factor*np.log10(change) + (greater-change))
            # Only low
            if not(greater.size):
                return lower
            # Only high
            if not(lower.size):
                return greater
            return np.concatenate((lower, greater))

        def inverted(self):
            return PlanckScale.InvertedPlanckTransform()

    class InvertedPlanckTransform(Transform):
        input_dims = 1
        output_dims = 1
        is_separable = True
        has_inverse = True

        def transform_non_affine(self, a):
            lower = a[np.where(a<=factor*np.log10(change))]
            greater = a[np.where(a> factor*np.log10(change))]
            if lower.size:
                if isinstance(lower, ma.MaskedArray):
                    lower = ma.power(10.0, lower/float(factor))
                else:
                    lower = np.power(10.0, lower/float(factor))
            if greater.size:
                greater = (greater + change - factor*np.log10(change))
            # Only low
            if not(greater.size):
                return lower
            # Only high
            if not(lower.size):
                return greater
            return np.concatenate((lower, greater))

        def inverted(self):
            return PlanckTransform()

# Finished. Register the scale!
mscale.register_scale(PlanckScale)

if __name__ == '__main__':
    sys.exit(main())
