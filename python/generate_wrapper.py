#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
import subprocess
rootdir = os.path.dirname(os.path.abspath(__file__)) + '/..'

h_files = []
for subdir, dirs, files in os.walk(rootdir):
    for file in files:
        ext = os.path.splitext(file)[-1].lower()
        if ext != '.h':
            continue
        h_files.append(os.path.join(subdir, file))


# In[2]:


preample = """
from libcpp cimport bool
from libcpp.map cimport map
from libcpp.memory cimport shared_ptr
from libcpp.pair cimport pair
from libcpp.string cimport string
from libcpp.vector cimport vector

"""


# In[3]:


definition_names = [
    '_MAX_NUMBER_OF_K_FILES_',
    '_MAXTITLESTRINGLENGTH_',
    '_FILENAMESIZE_',
    '_LINE_LENGTH_MAX_',
    '_Z_PK_NUM_MAX_',
    '_SELECTION_NUM_MAX_',
    '_ARGUMENT_LENGTH_MAX_',
    '_ERRORMSGSIZE_',
    '_TRUE_',
    '_FALSE_',
    '_SUCCESS_',
    '_FAILURE_',
]
definitions = []
for file in h_files:
    with open(file) as fid:
        for line in fid:
            if "#define" not in line:
                continue
            words = line.split()
            if len(words) < 3:
                continue
            if words[0].strip() != "#define":
                continue
            for d in definition_names:
                if d == words[1].strip():
                    definitions.append("DEF " + d + " = " + words[2])
                    continue


# In[4]:


enums = []
enums.append('cdef extern from "class.h":')
enums.append('    pair[string, string] get_my_py_error_message()')
enums.append('')
enums.append('    ctypedef char FileArg[_ARGUMENT_LENGTH_MAX_]')
enums.append('    ctypedef char ErrorMsg[_ERRORMSGSIZE_]')
enums.append('    ctypedef char FileName[_FILENAMESIZE_]')
enums.append('')
enum_names = [
    'equation_of_state',
    'file_format',
    'halofit_integral_type',
    'hmcode_baryonic_feedback_model',
    'linear_or_logarithmic',
    'non_linear_method',
    'out_sigmas',
    'pk_outputs',
    'primordial_spectrum_type',
    'selection_type',
    'source_extrapolation',
    'spatial_curvature',
]

for file in h_files:
    with open(file) as fid:
        enum_found = False
        for line in fid:
            if not enum_found:
                #Check for struct
                if 'enum' not in line:
                    continue
                elif '(' in line:
                    # Probably a function definition
                    continue
                elif '{' not in line:
                    # Probably a declaration of an enum variable in a struct.
                    continue
                for s in enum_names:
                    if 'enum ' + s + ' ' in line or 'enum ' + s + '\n' in line:
                        if '}' not in line:
                            #Assume multi-line enum:
                            enum_found = True
                        line = line.strip()
                        line = line.replace('{',':').replace('}','').replace(';','')
                        enums.append('    ctypedef ' + line)
                        break
            else:
                # Check for multiline enum ended
                if '};' in line:
                    enum_found = False
                elif line:
                    words = line.strip().split()
                    enums.append('       ' + words[0])


# In[5]:


structs = []
struct_names = [
    'background',
    'lensing',
    'nonlinear',
    'output',
    'perturbs',
    'precision',
    'primordial',
    'spectra',
    'thermo',
    'transfers',
]
allowed_types = ['double', 'int', 'short', 'FileArg']

for file in h_files:
    # Special treatment of common.h
    if file == 'common.h':
        subprocess.run(['gcc', '-E', rootdir+'/include/common.h','-o', rootdir+'/include/tmp'])
        file = rootdir+'/include/tmp'
    with open(file) as fid:
        struct_found = False
        for line in fid:
            if not struct_found:
                #Check for struct
                if 'struct' not in line:
                    continue
                for s in struct_names:
                    if 'struct ' + s + ' ' in line or 'struct ' + s + '\n' in line:
                        struct_found = True
                        structs.append('    cdef struct ' + s + ':')
            else:
                # Check for struct ended
                if '};' in line:
                    structs.append('')
                    struct_found = False
                else:
                    words = line.strip().split()
                    if len(words)>1 and words[0] in allowed_types:
                        if '*' in words[1]:
                            # It is a pointer, like this: double * a_pointer.
                            continue
                        structs.append('        ' + words[0] + ' ' + words[1].strip(';'))
                    elif len(words)>2 and words[0] == 'enum' and words[1] in enum_names:
                        structs.append('        ' + words[1] + ' ' + words[2].strip(';'))




# In[6]:

classes = []

class_names = ['FileContent','InputModule', 'BackgroundModule', 'ThermodynamicsModule', 'PerturbationsModule',
                'PrimordialModule', 'NonlinearModule', 'TransferModule', 'SpectraModule', 'LensingModule']
allowed_types = ['double', 'int', 'short', 'char', 'bool', 'void', 'ErrorMsg', 'FileArg',
                 'std::map<std::string, std::vector<double>>',
                 'std::map<std::string, int>'] + struct_names

for file in h_files:
    with open(file) as fid:
        class_name = ''
        for line in fid:
            if not class_name:
                # Check for struct
                if 'class' not in line and 'struct' not in line:
                    continue
                for s in class_names:
                    if (('struct ' + s + ' ' in line or 'struct ' + s + '\n' in line) or
                        ('class ' + s + ' ' in line or 'class ' + s + '\n' in line)):
                        class_name = s
                        classes.append('cdef extern from "' + file + '":')
                        classes.append('    cdef cppclass ' + s + ':')
            elif 'private:' in line or '};' in line:
                # Class just ended ended or we have reached the private section of the class
                if 'Module' in class_name:
                    # If module, add the inherited variable error_message_
                    classes.append('        ErrorMsg error_message_')
                class_name = ''
                classes.append('')
            else:
                line = line.strip()
                typename = ''
                for sometype in allowed_types:
                    if line.startswith(sometype):
                        typename = sometype
                        break
                if not typename:
                    continue

                variable_name_begin = -1
                for index in range(len(typename), len(line)):
                    if line[index] == '*':
                        typename += '*'
                    elif line[index] != ' ':
                        variable_name_begin = index
                        break
                if variable_name_begin == -1:
                    continue

                variable_name_end = -1
                parantheses_count = 0
                is_function = False
                for index in range(variable_name_begin, len(line)):
                    if ((line[index] == ' ' and parantheses_count == 0) or
                        (line[index] == ';' and parantheses_count == 0)):
                        variable_name_end = index
                        break
                    if line[index] == '(':
                        parantheses_count += 1
                        is_function = True
                    elif line[index] == ')':
                        parantheses_count -= 1

                if variable_name_end == -1:
                    continue

                variable_name = line[variable_name_begin:variable_name_end]
                variable_name = variable_name.replace('enum ','')
                if is_function:
                    variable_name += ' except +'

                out_line = '        ' + typename + ' ' + variable_name
                out_line = out_line.replace('std::','').replace('<', '[').replace('>',']')
                classes.append(out_line)


# In[7]:


modules = [m for m in class_names if 'Module' in m]
cosmology_class = []
cosmology_class.append('cdef extern from "cosmology.h":')
for m in modules:
    cosmology_class.append('    ctypedef shared_ptr[const ' + m + '] ' + m + 'Ptr')
cosmology_class.append('')
cosmology_class.append('    cdef cppclass Cosmology:')
cosmology_class.append('        Cosmology(FileContent& fc) except +')
for m in modules:
    cosmology_class.append('        ' + m + 'Ptr& Get' + m + '()')


# In[8]:


with open(rootdir + '/python/cclassy.pxd', 'w') as fid:
    fid.write(preample)
    for lines in [definitions, enums, structs, classes]:
        fid.write("\n".join(lines))
        fid.write("\n\n")
