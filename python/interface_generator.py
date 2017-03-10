"""
Automatically reads header files to generate an interface
"""
from __future__ import division, print_function
import sys
import logging
try:
    from collections import OrderedDict as od
except ImportError:
    try:
        from ordereddict import OrderedDict as od
    except ImportError:
        raise ImportError(
            "If you are running with Python v2.5 or 2.6"
            " you need to manually install the ordereddict"
            " package.")
try:
    import colorlog
except ImportError:
    raise ImportError(
        "You have to install the colorlog module"
        " with pip, or easy-install.")

SPACING = '    '
NAMING_CONVENTION = {
    'precision': {'python': 'precision',
                  'function': 'precision'},
    'background': {'python': 'background',
                   'function': 'background'},
    'thermo': {'python': 'thermodynamics',
               'function': 'thermodynamics'},
    'perturbs': {'python': 'perturbations',
                 'function': 'perturb'},
    'transfers': {'python': 'transfer',
                  'function': 'transfer'},
    'primordial': {'python': 'primordial',
                   'function': 'primordial'},
    'spectra': {'python': 'spectra',
                'function': 'spectra'},
    'lensing': {'python': 'lensing',
                'function': 'lensing'},
    'nonlinear': {'python': 'nonlinear',
                  'function': 'nonlinear'},
    'output': {'python': 'output',
               'function': 'output'},
    }


def main():
    # create logger
    logger = create_logger()

    # Recover all sub-header files
    main_header = '../include/class.h'
    headers = []

    with open(main_header, 'r') as header_file:
        in_modules = False
        for line in header_file:
            if in_modules:
                if line.strip() == '':
                    in_modules = False
                    continue
                if line.find('common') == -1 and line.find('input') == -1:
                    headers.append(
                        '../include/%s' % line.split()[-1].strip('"'))
            if line.find('class modules') != -1:
                in_modules = True

    logger.info('Extracted the following headers: %s', ', '.join(headers))
    output = 'classy.pyx'
    logger.info('Creating %s', output)
    structs = od()
    output_file = open(output, 'w')
    write_imports(output_file)
    output_file.write('cdef extern from "class.h":\n')

    # First write the first non automatic bits
    output_file.write(
        SPACING+'ctypedef char FileArg[40]\n' +
        SPACING+'ctypedef char* ErrorMsg\n' +
        SPACING+'cdef struct precision:\n' +
        2*SPACING+'ErrorMsg error_message\n\n' +
        SPACING+'cdef int _FAILURE_\n' +
        SPACING+'cdef int _FALSE_\n' +
        SPACING+'cdef int _TRUE_\n')

    for header in headers:
        extract_headers(header, structs, output_file, logger)
    logger.info("Finished extracting headers")
    for struct_name, struct in structs.items():
        create_wrapper_class(struct_name, struct, output_file, logger)

    return


def extract_headers(header, structs, output_file, logger):
    """toto"""
    # Initialise the two flags controlling the exploration of the main
    # structure
    in_struct, main_struct_finished = False, False
    # Flags for exploring enums (only the ones before the struct)
    in_enum = False
    # flag dealing with extracting docstrings
    comment_partially_recovered = False
    # Flag keeping track of multiple variables
    multiple_var = False
    # Flag recovering the functions
    in_function_definitions, in_function, in_init = False, False, False
    with open(header, 'r') as header_file:
        logger.info("reading %s" % header)
        for line in header_file:
            # First case, recover the enums
            if not main_struct_finished and not in_struct:
                if line.find("enum ") != -1 and line.find("{") != -1:
                    enum_members = []
                    if line.find(";") == -1:
                        in_enum = True
                        enum_name = line.strip("enum").strip().strip('{')
                    else:
                        in_enum = False
                        line = line.strip("enum").strip().strip(';')
                        enum_name, enum_sign = line.split(' ', 1)
                        enum_sign = enum_sign.strip('}').strip('{')
                        for elem in enum_sign.split(','):
                            enum_members.append(elem.strip())
                        output_file.write(
                            SPACING + 'cdef enum %s:\n' % enum_name)
                        for elem in enum_members:
                            output_file.write(2*SPACING + elem + '\n')
                        output_file.write('\n')
                elif in_enum:
                    if line.find('};') != -1:
                        in_enum = False
                        output_file.write(
                            SPACING + 'cdef enum %s:\n' % enum_name)
                        for elem in enum_members:
                            output_file.write(2*SPACING + elem + '\n')
                        output_file.write('\n')
                    else:
                        if line.strip() != '':
                            enum_members.append(line.split()[0].strip().strip(','))
            if line.find("struct ") != -1 and not main_struct_finished:
                in_struct = True
                # Recover the name
                logger.debug("in struct: %s" % line)
                struct_name = line.strip().split()[1]
                logger.debug("struct name: %s" % struct_name)
                structs[struct_name] = {}
                structs[struct_name].update(
                    NAMING_CONVENTION[struct_name])
                output_file.write("%scdef struct %s:\n" % (
                    SPACING, struct_name))
                continue
            elif in_struct:
                if line.find("};\n") != -1:
                    output_file.write('\n')
                    in_struct, main_struct_finished = False, True
                else:
                    # if the line is not empty or does not contain only a
                    # comment:
                    if line.strip() == '' or line.strip()[:2] == '/*':
                        continue
                    logger.debug(
                        "potentially non empty line: %s" % line.strip())
                    #elif line.find('/**') != -1 or line.find('*/') != -1:
                        #continue
                    if line.find(';') == -1 and not comment_partially_recovered:
                        logger.debug("--> Discarded")
                        continue
                    elif line.find(';') != -1 and not comment_partially_recovered:
                        var_doc = ''
                        var_part, begin_comment = line.strip().split(';', 1)
                        var_doc += begin_comment.strip()[4:].strip()
                        # 2 things can happen: there can be arrays, and there
                        # can be several variables defined in one line...
                        # If array, slightly more complex
                        if var_part.find('*') != -1:
                            # if no comma is found, it means it is a single
                            # variable: good !
                            if var_part.find(',') == -1:
                                # remove if commented (starts with /*)
                                if var_part[:2] in ['/*', '//']:
                                    continue
                                multiple_var = False
                                var_type, var_stars, var_name = var_part.strip().split()
                                structs[struct_name][var_name] = [
                                    var_type, var_stars]
                            else:
                                # Count how many variables are defined
                                multiple_var = True
                                all_vars = [elem.strip() for elem in
                                            var_part.split('*')[-1].split(',')]
                                var_type, var_stars = (var_part.strip().
                                                       split()[:2])
                                for var in all_vars:
                                    structs[struct_name][var] = [
                                        var_type, var_stars]
                        else:
                            # Again check for more than one variable
                            var_stars = ''
                            if var_part.find(',') == -1:
                                multiple_var = False
                                var_type, var_name = var_part.strip().split(' ', 1)
                                # Check if enum
                                if var_type == 'enum':
                                    enum_name, var_name = var_name.split()
                                    var_type += ' '+enum_name
                                structs[struct_name][var_name] = [
                                    var_type, var_stars]
                            else:
                                multiple_var = True
                                all_vars = [elem.strip() for elem in
                                            var_part.split()[2:].split(',')]
                                var_type = (var_part.strip().split()[0])
                                for var in all_vars:
                                    structs[struct_name][var] = [
                                        var_type, var_stars]
                        # If the comment is finished, pass
                        if var_doc[-2:] != '*/':
                            comment_partially_recovered = True
                        else:
                            var_doc = var_doc[:-2].replace('\\f$', '$').strip()
                            structs[struct_name][var_name].append(var_doc)
                            logger.debug(
                                "extracted the variable %s, " % var_name +
                                "of type %s, with docstring: %s" % (
                                    ''.join([var_stars, var_type]), var_doc))

                        if not multiple_var:
                            output_file.write(2*SPACING+' '.join(
                                [elem for elem in [var_type, var_stars, var_name]
                                 if elem])+'\n')
                        else:
                            for var in all_vars:
                                output_file.write(2*SPACING+' '.join(
                                    [elem for elem in [var_type, var_stars, var]
                                     if elem])+'\n')

                    if comment_partially_recovered:
                        logger.debug("--> Accepted")
                        var_doc += ' '+line.strip()
                        if var_doc[-2:] == '*/':
                            comment_partially_recovered = False
                            var_doc = var_doc[:-2].replace('\\f$', '$').strip()
                            structs[struct_name][var_name].append(var_doc)
                            logger.debug(
                                "extracted the variable %s, " % var_name +
                                "of type %s, with docstring: %s" % (
                                    ''.join([var_stars, var_type]), var_doc))

            elif main_struct_finished:
                if line.find('extern "C"') != -1:
                    in_function_definitions = True
                if not in_function_definitions:
                    continue
                else:
                    if line.find('(') != -1:
                        in_function = True
                        logger.debug("Found a function")
                        func_type, func_name = line.split('(')[0].strip().split()
                        logger.debug('%s %s' % (func_name, func_type))
                        func_param = []

                        if func_name == structs[struct_name]['function']+'_init':
                            logger.info("found the init function")
                            in_init = True
                            structs[struct_name]['init'] = [func_name]
                        output_file.write(SPACING+'%s %s(' % (
                            func_type, func_name))
                    elif in_function:
                        # recover the signature of the function
                        line = line.strip().strip(',')
                        if line.find('struct') != -1:
                            if in_init:
                                name = line.split('*')[0].strip()[7:]
                                structs[struct_name]['init'].append(name)
                            func_param.append('void *')
                        elif line.find('*') != -1:
                            # Taking into account with or without spaces
                            temp = ''.join(line.strip(',').split())
                            last_star = len(temp)-temp[::-1].find('*')
                            func_param.append(temp[:last_star])
                        elif line.find(')') == -1:
                            if line != '':
                                func_param.append(line.split()[0])
                        else:
                            logger.debug('signature extracted')
                            in_function = False
                            if in_init:
                                in_init = False
                            output_file.write(', '.join(func_param) + ')\n')

                    elif line.find('}') != -1:
                        output_file.write('\n')
                        in_function_definitions = False
            #print line.strip()


def create_wrapper_class(struct_name, struct, of, logger):
    """TODO"""
    of.write('# Defining wrapper around struct %s\n' % struct_name)
    of.write('cdef class %s:\n' % (
        NAMING_CONVENTION[struct_name]['python'].capitalize()))

    ## recover the number of additional arguments:
    init_name, argument_names = struct['init'][0], struct['init'][1:]

    for companion in argument_names:
        of.write(SPACING+'cdef %s _%s\n' % (companion, companion))
        #logger.info("structure: %s, python name: %s" % (
            #companion, NAMING_CONVENTION[companion]['python']))
    of.write('\n')

    # Define the array variables for all needed
    array_variables = []
    variables = []
    for key, value in struct.items():
        if key != 'init':
            if value[1]:
                array_variables.append(key)
                variables.append(key)
                of.write(SPACING+'cdef np.ndarray %s_arr\n' % key)
            else:
                variables.append(key)
    of.write('\n')

    # write the init
    of.write(SPACING+'def __init__(self')
    for companion in argument_names:
        of.write(", %s py_%s" % (
            NAMING_CONVENTION[companion]['python'].capitalize(), companion))
    of.write('):\n\n')

    # pointing the pointers where they belong
    for companion in argument_names:
        of.write(2*SPACING+"self._%s = py_%s._%s\n" % (
            companion, companion, companion))

    # Writing the call to structname_init()
    of.write(2*SPACING+'%s_init(\n' % struct_name)
    for companion in argument_names:
        of.write(3*SPACING+'&(self._%s),\n' % companion)
    of.write(3*SPACING+'&(self._%s))\n\n' % struct_name)

    #of.write(2*SPACING+'%s_init(&(self._%s))\n\n' % (
        #struct_name, struct_name))
    for array in array_variables:
        of.write(2*SPACING+'# Wrapping %s\n' % array)
        of.write(2*SPACING+'%s_wrapper = ArrayWrapper()\n' % array)
        of.write(
            2*SPACING+"%s_wrapper.set_data(%d, '%s', "
            "<void*> self._%s.%s)\n" % (
                array, 2, struct[array].strip('*'), struct_name, array))
        of.write(
            2*SPACING+'self.%s_arr = np.array(%s_wrapper, '
            'copy=False)\n' % (
                array, array))
        of.write(2*SPACING+'self.%s_arr.base = '
                 '<PyObject*> %s_wrapper\n' % (
                     array, array))
        of.write(2*SPACING+'Py_INCREF(%s_wrapper)\n\n' % array)
    #raise NotImplementedError('multiple init are not supported')

    # Write the properties
    for key in variables:
        of.write(SPACING+'property %s:\n' % key)
        if key not in array_variables:
            of.write(2*SPACING+'def __get__(self):\n')
            of.write(3*SPACING+'return self._%s.%s\n' % (struct_name, key))
            of.write(2*SPACING+'def __set__(self, rhs):\n')
            of.write(3*SPACING+'self._%s.%s = rhs\n' % (struct_name, key))
        else:
            of.write(2*SPACING+'def __get__(self):\n')
            of.write(3*SPACING+'return self.%s_arr\n' % key)
            of.write(2*SPACING+'def __set__(self, rhs):\n')
            of.write(3*SPACING+'self.%s_arr[:] = rhs\n' % key)
        of.write('\n')

    # Add blank lines
    of.write('\n\n')


def write_imports(output_file):
    """TODO"""
    a = '''# Author: Gael Varoquaux
# License: BSD
from libc.stdlib cimport free
from cpython cimport PyObject, Py_INCREF

# Import the Python-level symbols of numpy
import numpy as np

# Import the C-level symbols of numpy
cimport numpy as np

# Numpy must be initialized. When using numpy from C or Cython you must
# _always_ do that, or you will have segfaults
np.import_array()

cdef class ArrayWrapper:
    cdef void* data_ptr
    cdef int size
    cdef int type

    cdef set_data(self, int size, char* type, void* data_ptr):
        """ Set the data of the array

        This cannot be done in the constructor as it must recieve C-level
        arguments.

        Parameters:
        -----------
        size: int
            Length of the array.
        data_ptr: void*
            Pointer to the data

        """
        self.data_ptr = data_ptr
        self.size = size
        if type.find('int') != -1:
            self.type = np.NPY_INT
        elif type.find('float') != -1:
            self.type = np.NPY_FLOAT
        elif type.find('double') != -1:
            self.type = np.NPY_DOUBLE
        elif type.find('long') != -1:
            self.type = np.NPY_LONG

    def __array__(self):
        """ Here we use the __array__ method, that is called when numpy
        tries to get an array from the object."""
        cdef np.npy_intp shape[1]
        shape[0] = <np.npy_intp> self.size
        # Create a 1D array, of length 'size'
        ndarray = np.PyArray_SimpleNewFromData(1, shape,
            self.type, self.data_ptr)
        return ndarray

    def __dealloc__(self):
        """ Frees the array. This is called by Python when all the
        references to the object are gone. """
        free(<void*>self.data_ptr)\n\n'''
    output_file.write(a)


def create_logger():
    """Nothing"""

    logger = logging.getLogger('simple_example')
    #logger.setLevel(logging.DEBUG)
    logger.setLevel(logging.INFO)

    # create console handler and set level to debug
    console_handler = logging.StreamHandler()
    #console_handler.setLevel(logging.DEBUG)
    console_handler.setLevel(logging.INFO)

    # create formatter
    #formatter = logging.Formatter(
        #"%(asctime)s %(module)s: L%(lineno) 4s %(funcName) 15s"
        #" | %(levelname) -10s  --> %(message)s")
    formatter = colorlog.ColoredFormatter(
        "%(asctime)s %(module)s: L%(lineno) 4s %(blue)s%(funcName) 15s%(reset)s"
        " | %(log_color)s%(levelname) -10s  --> %(message)s%(reset)s",
        datefmt=None,
        reset=True,
        log_colors={
            'DEBUG':    'cyan',
            'INFO':     'green',
            'WARNING':  'yellow',
            'ERROR':    'red',
            'CRITICAL': 'red',
        })

    # add formatter to console_handler
    console_handler.setFormatter(formatter)

    # add console_handler to logger
    logger.addHandler(console_handler)

    return logger

if __name__ == "__main__":
    sys.exit(main())
