"""
.. module:: test_class
    :synopsis: python script for testing CLASS using nose
.. moduleauthor:: Benjamin Audren <benjamin.audren@gmail.com>
.. credits:: Benjamin Audren, Thomas Tram
.. version:: 1.0

This is a python script for testing CLASS and its wrapper Classy using nose.
To run the test suite, type
nosetests test_class.py
If you want to extract the problematic input parameters at a later stage,
you should type
nosetests test_class.py 1>stdoutfile 2>stderrfile
and then use the python script extract_errors.py on the stderrfile.

When adding a new input parameter to CLASS (by modifying input.c), you
should also include tests of this new input. You will be in one of the
two cases:
1:  The new input is supposed to be compatible with any existing input.
    This is the standard case when adding a new species for instance.
2:  The new input is incompatible with one of the existing inputs. This
    would be the case if you have added (or just want to test) some other
    value of an already defined parameter. (Maybe you have allowed for
    negative mass neutrinos and you want to test CLASS using a negative mass.)

In case 1, you must add an entry in the CLASS_INPUT dictionary:
CLASS_INPUT['Mnu'] = (
    [{'N_eff': 0.0, 'N_ncdm': 1, 'm_ncdm': 0.06, 'deg_ncdm': 3.0},
     {'N_eff': 1.5, 'N_ncdm': 1, 'm_ncdm': 0.03, 'deg_ncdm': 1.5}],
    'normal')
The key 'Mnu' is not being used in the code, so its purpose is just to
describe the entry to the reader.
the value is a 2-tuple where the first entry [{},{},...,{}] is an array of
dictionaries containg the actual input to CLASS. The second entry is a keyword
which can be either 'normal' or 'power'. It tells the script how this input
will be combined with other inputs.

What does 'normal' and 'power' mean?
If an entry has the 'power' keyword, it will be combined with any other entry.
If an entry has the 'normal' keyword, it will not be combined with any other
entry having the 'normal' keyword, but it will be combined with all entries
carrying the 'power keyword.
Beware that the number of tests grow a lot when using the 'power' keyword.

In case 2, you should find the relevant entry and just add a new dictionary
to the array. E.g. if you want to test some negative mass model you should add
{'N_ncdm': 1, 'm_ncdm': -0.1, 'deg_ncdm': 1.0}

How are default parameters handled?
Any input array implicitly contains the empty dictionary. That means that if
Omega_k:0.0 is the default value, writing
CLASS_INPUT['Curvature'] = (
    [{'Omega_k': 0.01},
     {'Omega_k': -0.01}],
    'normal')
will test the default value Omega_k=0.0 along with the two specified models.

How to deal with inconsistent input?
Sometimes a specific feature requires the presence of another input parameter.
For instance, if we ask for tensor modes we must have temperature and/or
polarisation in the output. If not, CLASS is supposed to fail during the
evaluation of the input module and return an error message. This fail is the
correct behaviour of CLASS. To implement such a case, modify the function
test_incompatible_input(self)

Comparing output: When the flag 'COMPARE_OUTPUT' is set to true, the code will
rerun CLASS for each case under Newtonian gauge and then compare Cl's and
matter power spectrum. If the two are not close enough, it will generate a
PDF plot of this and save it in the 'fail' folder.
"""
from __future__ import print_function
from classy import Class
from classy import CosmoSevereError
import itertools
import sys
import shutil
import os
import numpy as np
from math import log10
import matplotlib.pyplot as plt
import unittest
from nose_parameterized import parameterized

# To avoid testing for differences between synchronous and Newtonian gauge, set
# this flag to False
COMPARE_OUTPUT = True

# Dictionary of models to test the wrapper against. Each of these scenario will
# be run against all the possible output choices (nothing, tCl, mPk, etc...),
# with or without non-linearities.
# Never input the default value, as this will be **automatically** tested
# against. Indeed, when not specifying a field, CLASS takes the default input.
CLASS_INPUT = {}

#CLASS_INPUT['Mnu'] = (
#    [{'N_eff': 0.0, 'N_ncdm': 1, 'm_ncdm': 0.06, 'deg_ncdm': 3.0},
#     {'N_eff': 1.5, 'N_ncdm': 1, 'm_ncdm': 0.03, 'deg_ncdm': 1.5}],
#    'normal')

#CLASS_INPUT['Curvature'] = (
#    [{'Omega_k': 0.01},
#     {'Omega_k': -0.01}],
#    'normal')

CLASS_INPUT['Isocurvature_modes'] = (
    [{'ic': 'ad,nid,cdi', 'c_ad_cdi': -0.5}],
    'normal')

#CLASS_INPUT['Scalar_field'] = (
    #[{'Omega_scf': 0.1, 'attractor_ic_scf': 'yes',
      #'scf_parameters': '10, 0, 0, 0'}],
    #'normal')

CLASS_INPUT['Inflation'] = (
    [{'P_k_ini type': 'inflation_V'},
     {'P_k_ini type': 'inflation_H'},
     {'P_k_ini type': 'inflation_V_end'}],
    'normal')

CLASS_INPUT['modes'] = (
    [{'modes': 't'},
     {'modes': 's, t'}],
    'power')

CLASS_INPUT['Tensor_method'] = (
    [{'tensor method': 'exact'},
     {'tensor method': 'photons'}],
    'power')

CLASS_INPUT['Output_spectra'] = (
    [{'output': 'mPk', 'P_k_max_1/Mpc': 10},
     {'output': 'tCl'},
     {'output': 'tCl pCl lCl'},
     {'output': 'mPk tCl lCl', 'P_k_max_1/Mpc': 10},
     {'output': 'nCl sCl'},
     {'output': 'tCl pCl lCl nCl sCl'}],
    'power')

CLASS_INPUT['Nonlinear'] = (
    [{'non linear': 'halofit'}],
    'power')

CLASS_INPUT['Lensing'] = (
    [{'lensing': 'yes'}],
    'power')

# Let's kill the machine (replace all 'normal' flags with power', uncomment at
# you own risk)
# for k, v in CLASS_INPUT.iteritems():
#     models, state = v
#     CLASS_INPUT[k] = (models, 'power')

INPUTPOWER = []
INPUTNORMAL = [{}]
for key, value in list(CLASS_INPUT.items()):
    models, state = value
    if state == 'power':
        INPUTPOWER.append([{}]+models)
    else:
        INPUTNORMAL.extend(models)

    PRODPOWER = list(itertools.product(*INPUTPOWER))

    DICTARRAY = []
    for normelem in INPUTNORMAL:
        for powelem in PRODPOWER:  # itertools.product(*modpower):
            temp_dict = normelem.copy()
            for elem in powelem:
                temp_dict.update(elem)
            DICTARRAY.append(temp_dict)


TUPLE_ARRAY = []
for e in DICTARRAY:
    TUPLE_ARRAY.append((e, ))


def powerset(iterable):
    xs = list(iterable)
    # note we return an iterator rather than a list
    return itertools.chain.from_iterable(
        itertools.combinations(xs, n) for n in range(1, len(xs)+1))


class TestClass(unittest.TestCase):
    """
    Testing Class and its wrapper classy on different cosmologies

    To run it, do
    ~] nosetest test_class.py

    It will run many times Class, on different cosmological scenarios, and
    everytime testing for different output possibilities (none asked, only mPk,
    etc..)

    """
    @classmethod
    def setUpClass(self):
        self.faulty_figs_path = os.path.join(
            os.path.sep.join(os.path.realpath(__file__).split(
                os.path.sep)[:-1]),
            'faulty_figs')

        if os.path.isdir(self.faulty_figs_path):
            shutil.rmtree(self.faulty_figs_path)

        os.mkdir(self.faulty_figs_path)

    @classmethod
    def tearDownClass(self):
        pass

    def setUp(self):
        """
        set up data used in the tests.
        setUp is called before each test function execution.
        """
        self.cosmo = Class()
        self.cosmo_newt = Class()

        self.verbose = {
            'input_verbose': 1,
            'background_verbose': 1,
            'thermodynamics_verbose': 1,
            'perturbations_verbose': 1,
            'transfer_verbose': 1,
            'primordial_verbose': 1,
            'harmonic_verbose': 1,
            'fourier_verbose': 1,
            'lensing_verbose': 1,
            'output_verbose': 1}
        self.scenario = {}

    def tearDown(self):
        self.cosmo.struct_cleanup()
        self.cosmo.empty()
        self.cosmo_newt.struct_cleanup()
        self.cosmo_newt.empty()
        del self.scenario

    def poormansname(self, somedict):
        string = "_".join(
            [k+'='+str(v)
             for k, v in list(somedict.items())])
        string = string.replace('/', '%')
        string = string.replace(',', '')
        string = string.replace(' ', '')
        return string

    @parameterized.expand(TUPLE_ARRAY)
    def test_0wrapper_implementation(self, inputdict):
        """Create a few instances based on different cosmologies"""
        self.scenario.update(inputdict)

        self.name = self.poormansname(inputdict)

        sys.stderr.write('\n\n---------------------------------\n')
        sys.stderr.write('| Test case %s |\n' % self.name)
        sys.stderr.write('---------------------------------\n')
        for key, value in list(self.scenario.items()):
            sys.stderr.write("%s = %s\n" % (key, value))
            sys.stdout.write("%s = %s\n" % (key, value))
        sys.stderr.write("\n")

        setting = self.cosmo.set(
            dict(list(self.verbose.items())+list(self.scenario.items())))
        self.assertTrue(setting, "Class failed to initialize with input dict")

        cl_dict = {
            'tCl': ['tt'],
            'lCl': ['pp'],
            'pCl': ['ee', 'bb']}
        density_cl_list = ['nCl', 'sCl']

        # 'lensing' is always set to yes. Therefore, trying to compute 'tCl' or
        # 'pCl' will fail except if we also ask for 'lCl'. The flag
        # 'should_fail' stores this status.
        sys.stderr.write('Should')
        should_fail = self.test_incompatible_input()
        if should_fail:
            sys.stderr.write(' fail...\n')
        else:
            sys.stderr.write(' not fail...\n')

        if not should_fail:
            self.cosmo.compute()
        else:
            self.assertRaises(CosmoSevereError, self.cosmo.compute)
            return

        self.assertTrue(
            self.cosmo.state,
            "Class failed to go through all __init__ methods")
        if self.cosmo.state:
            print('--> Class is ready')
        # Depending
        if 'output' in list(self.scenario.keys()):
            # Positive tests of raw cls
            output = self.scenario['output']
            for elem in output.split():
                if elem in list(cl_dict.keys()):
                    for cl_type in cl_dict[elem]:
                        sys.stderr.write(
                            '--> testing raw_cl for %s\n' % cl_type)
                        cl = self.cosmo.raw_cl(100)
                        self.assertIsNotNone(cl, "raw_cl returned nothing")
                        self.assertEqual(
                            np.shape(cl[cl_type])[0], 101,
                            "raw_cl returned wrong size")
                    # TODO do the same for lensed if 'lCl' is there, and for
                    # density cl
                if elem == 'mPk':
                    sys.stderr.write('--> testing pk function\n')
                    pk = self.cosmo.pk(0.1, 0)
                    self.assertIsNotNone(pk, "pk returned nothing")
            # Negative tests of output functions
            if not any([elem in list(cl_dict.keys()) for elem in output.split()]):
                sys.stderr.write('--> testing absence of any Cl\n')
                self.assertRaises(CosmoSevereError, self.cosmo.raw_cl, 100)
            if 'mPk' not in output.split():
                sys.stderr.write('--> testing absence of mPk\n')
                self.assertRaises(CosmoSevereError, self.cosmo.pk, 0.1, 0)

        if COMPARE_OUTPUT:
            # Now, compute with Newtonian gauge, and compare the results
            self.cosmo_newt.set(
                dict(list(self.verbose.items())+list(self.scenario.items())))
            self.cosmo_newt.set({'gauge': 'newtonian'})
            self.cosmo_newt.compute()
            # Check that the computation worked
            self.assertTrue(
                self.cosmo_newt.state,
                "Class failed to go through all __init__ methods in Newtonian gauge")

            self.compare_output(self.cosmo, self.cosmo_newt)

    def test_incompatible_input(self):

        should_fail = False

        # If we have tensor modes, we must have one tensor observable,
        # either tCl or pCl.
        if has_tensor(self.scenario):
            if 'output' not in list(self.scenario.keys()):
                should_fail = True
            else:
                output = self.scenario['output'].split()
                if 'tCl' not in output and 'pCl' not in output:
                    should_fail = True

        # If we have specified lensing, we must have lCl in output,
        # otherwise lensing will not be read (which is an error).
        if 'lensing' in list(self.scenario.keys()):
            if 'output' not in list(self.scenario.keys()):
                should_fail = True
            else:
                output = self.scenario['output'].split()
                if 'lCl' not in output:
                    should_fail = True
                elif 'tCl' not in output and 'pCl' not in output:
                    should_fail = True

        # If we have specified a tensor method, we must have tensors.
        if 'tensor method' in list(self.scenario.keys()):
            if not has_tensor(self.scenario):
                should_fail = True

        # If we have specified non linear, we must have some form of
        # perturbations output.
        if 'non linear' in list(self.scenario.keys()):
            if 'output' not in list(self.scenario.keys()):
                should_fail = True

        # If we ask for Cl's of lensing potential, we must have scalar modes.
        if 'output' in list(self.scenario.keys()) and 'lCl' in self.scenario['output'].split():
            if 'modes' in list(self.scenario.keys()) and self.scenario['modes'].find('s') == -1:
                should_fail = True

        # If we specify initial conditions (for scalar modes), we must have
        # perturbations and scalar modes.
        if 'ic' in list(self.scenario.keys()):
            if 'modes' in list(self.scenario.keys()) and self.scenario['modes'].find('s') == -1:
                should_fail = True
            if 'output' not in list(self.scenario.keys()):
                should_fail = True

        # If we use inflation module, we must have scalar modes,
        # tensor modes, no vector modes and we should only have adiabatic IC:
        if 'P_k_ini type' in list(self.scenario.keys()) and self.scenario['P_k_ini type'].find('inflation') != -1:
            if 'modes' not in list(self.scenario.keys()):
                should_fail = True
            else:
                if self.scenario['modes'].find('s') == -1:
                    should_fail = True
                if self.scenario['modes'].find('v') != -1:
                    should_fail = True
                if self.scenario['modes'].find('t') == -1:
                    should_fail = True
            if 'ic' in list(self.scenario.keys()) and self.scenario['ic'].find('i') != -1:
                should_fail = True


        return should_fail

    def compare_output(self, reference, candidate):
        sys.stderr.write('\n\n---------------------------------\n')
        sys.stderr.write('| Comparing synch and Newt: |\n')
        sys.stderr.write('---------------------------------\n')

        for elem in ['raw_cl', 'lensed_cl', 'density_cl']:
            # Try to get the elem, but if they were not computed, a
            # CosmoComputeError should be raised. In this case, ignore the
            # whole block.
            try:
                to_test = getattr(candidate, elem)()
            except CosmoSevereError:
                continue
            ref = getattr(reference, elem)()
            for key, value in list(ref.items()):
                if key != 'ell':
                    sys.stderr.write('--> testing equality of %s %s\n' % (
                        elem, key))
                    # For all self spectra, try to compare allclose
                    if key[0] == key[1]:
                        # If it is a 'dd' or 'll', it is a dictionary.
                        if isinstance(value, dict):
                            for subkey in list(value.keys()):
                                try:
                                    np.testing.assert_allclose(
                                        value[subkey], to_test[key][subkey],
                                        rtol=1e-03, atol=1e-20)
                                except AssertionError:
                                    self.cl_faulty_plot(elem+"_"+key,
                                                        value[subkey][2:],
                                                        to_test[key][subkey][2:])
                                except TypeError:
                                    self.cl_faulty_plot(elem+"_"+key,
                                                        value[subkey][2:],
                                                        to_test[key][subkey][2:])
                        else:
                            try:
                                np.testing.assert_allclose(
                                    value, to_test[key], rtol=1e-03, atol=1e-20)
                            except AssertionError:
                                self.cl_faulty_plot(elem+"_"+key,
                                                    value[2:], to_test[key][2:])
                            except TypeError:
                                self.cl_faulty_plot(elem+"_"+key,
                                                    value[2:], to_test[key][2:])
                    # For cross-spectra, as there can be zero-crossing, we
                    # instead compare the difference.
                    else:
                        # First, we multiply each array by the biggest value
                        norm = max(
                            np.abs(value).max(), np.abs(to_test[key]).max())
                        value *= norm
                        to_test[key] *= norm
                        try:
                            np.testing.assert_array_almost_equal(
                                value, to_test[key], decimal=3)
                        except AssertionError:
                            self.cl_faulty_plot(elem+"_"+key,
                                                value[2:], to_test[key][2:])

        if 'output' in list(self.scenario.keys()):
            if self.scenario['output'].find('mPk') != -1:
                sys.stderr.write('--> testing equality of Pk')
                k = np.logspace(
                    -2, log10(self.scenario['P_k_max_1/Mpc']))
                reference_pk = np.array(
                    [reference.pk(elem, 0) for elem in k])
                candidate_pk = np.array(
                    [candidate.pk(elem, 0) for elem in k])
                try:
                    np.testing.assert_allclose(
                        reference_pk, candidate_pk, rtol=5e-03, atol=1e-20)
                except AssertionError:
                    self.pk_faulty_plot(k, reference_pk, candidate_pk)

    def cl_faulty_plot(self, cl_type, reference, candidate):
        path = os.path.join(self.faulty_figs_path, self.name)

        fig = plt.figure()
        ax_lin = plt.subplot(211)
        ax_log = plt.subplot(212)
        ell = np.arange(max(np.shape(candidate)))+2
        ax_lin.plot(ell, 1-candidate/reference)
        ax_log.loglog(ell, abs(1-candidate/reference))

        ax_lin.set_xlabel('l')
        ax_log.set_xlabel('l')
        ax_lin.set_ylabel('1-candidate/reference')
        ax_log.set_ylabel('abs(1-candidate/reference)')

        ax_lin.set_title(self.name)
        ax_log.set_title(self.name)

        ax_lin.legend([cl_type])
        ax_log.legend([cl_type])

        fig.savefig(path+'_'+cl_type+'.pdf')

        # Store parameters (contained in self.scenario) to text file
        parameters = dict(list(self.verbose.items())+list(self.scenario.items()))
        with open(path+'.ini', 'w') as param_file:
            for key, value in list(parameters.items()):
                param_file.write(key+" = "+str(value)+'\n')

    def pk_faulty_plot(self, k, reference, candidate):
        path = os.path.join(self.faulty_figs_path, self.name)

        fig = plt.figure()
        ax_lin = plt.subplot(211)
        ax_log = plt.subplot(212)
        ax_lin.plot(k, 1-candidate/reference)
        ax_log.loglog(k, abs(1-candidate/reference))

        ax_lin.set_xlabel('k')
        ax_log.set_xlabel('k')
        ax_lin.set_ylabel('1-candidate/reference')
        ax_log.set_ylabel('abs(1-candidate/reference)')

        ax_lin.set_title(self.name)
        ax_log.set_title(self.name)

        ax_lin.legend('$P_k$')
        ax_log.legend('$P_k$')

        fig.savefig(path+'_'+'pk'+'.pdf')

        # Store parameters (contained in self.scenario) to text file
        parameters = dict(list(self.verbose.items())+list(self.scenario.items()))
        with open(path+'.ini', 'w') as param_file:
            for key, value in list(parameters.items()):
                param_file.write(key+" = "+str(value)+'\n')


def has_tensor(input_dict):
    if 'modes' in list(input_dict.keys()):
        if input_dict['modes'].find('t') != -1:
            return True
    else:
        return False
    return False

if __name__ == '__main__':
    toto = TestClass()
    unittest.main()
