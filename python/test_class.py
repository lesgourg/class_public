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
has_incompatible_input(self)

Comparing output: When the flag 'COMPARE_OUTPUT_GAUGE' is set to true, the code will
rerun CLASS for each case under Newtonian gauge and then compare Cl's and
matter power spectrum. If the two are not close enough, it will generate a
PDF plot of this and save it in the 'fail' folder.
"""
from __future__ import print_function
from six import iteritems
#from future.utils import iteritems
import matplotlib as mpl
mpl.use('Agg')

import itertools
import matplotlib.pyplot as plt
import numpy as np
import os
import shutil
import unittest

from classy import Class
from classy import CosmoSevereError
from math import log10
from matplotlib.offsetbox import AnchoredText
from nose.plugins.attrib import attr
from parameterized import parameterized

# Customise test by reading environment variables
CLASS_VERBOSE = bool(int(os.getenv('CLASS_VERBOSE', '0'))) # Print output from CLASS?
COMPARE_OUTPUT_GAUGE = bool(int(os.getenv('COMPARE_OUTPUT_GAUGE', '0'))) # Compare synchronous and Newtonian gauge outputs?
COMPARE_OUTPUT_REF = bool(int(os.getenv('COMPARE_OUTPUT_REF', '0'))) # Compare classy with classyref?
POWER_ALL = bool(int(os.getenv('POWER_ALL', '0'))) # Combine every extension with each other? (Very slow!)
TEST_LEVEL = int(os.getenv('TEST_LEVEL', '0')) # 0 <= TEST_LEVEL <= 3

if COMPARE_OUTPUT_REF:
    try:
        import classyref
    except:
        COMPARE_OUTPUT_REF = False
# Define bounds on the relative and absolute errors of C(l) and P(k)
# between reference, Newtonian and Synchronous gauge
COMPARE_CL_RELATIVE_ERROR = 3e-3
COMPARE_CL_RELATIVE_ERROR_GAUGE = 5*3e-3
COMPARE_CL_ABSOLUTE_ERROR = 1e-20
COMPARE_PK_RELATIVE_ERROR = 1e-2
COMPARE_PK_RELATIVE_ERROR_GAUGE = 5*1e-2
COMPARE_PK_ABSOLUTE_ERROR = 1e-20

# Dictionary of models to test the wrapper against. Each of these scenario will
# be run against all the possible output choices (nothing, tCl, mPk, etc...),
# with or without non-linearities.
# Never input the default value, as this will be **automatically** tested
# against. Indeed, when not specifying a field, CLASS takes the default input.
CLASS_INPUT = {}

CLASS_INPUT['Output_spectra'] = (
    [{'output': 'mPk', 'P_k_max_1/Mpc': 2},
     {'output': 'tCl'},
     {'output': 'tCl pCl lCl'},
     {'output': 'mPk tCl lCl', 'P_k_max_1/Mpc': 2},
     {'output': 'nCl sCl'},
     {'output': 'tCl pCl nCl sCl'}],
    'power')

CLASS_INPUT['Nonlinear'] = (
    [{'non linear': 'halofit'}],
    'power')

CLASS_INPUT['Lensing'] = (
    [{'lensing': 'yes'}],
    'power')

if TEST_LEVEL > 0:
    CLASS_INPUT['Mnu'] = (
        [{'N_ur': 0.0, 'N_ncdm': 1, 'm_ncdm': 0.06, 'deg_ncdm': 3.0},
        {'N_ur': 1.5, 'N_ncdm': 1, 'm_ncdm': 0.03, 'deg_ncdm': 1.5}],
        'normal')

if TEST_LEVEL > 1:
    CLASS_INPUT['Curvature'] = (
        [{'Omega_k': 0.01},
        {'Omega_k': -0.01}],
        'normal')

    CLASS_INPUT['modes'] = (
        [{'modes': 't'},
        {'modes': 's, t'}],
        'power')

    CLASS_INPUT['Tensor_method'] = (
        [{'tensor method': 'exact'},
        {'tensor method': 'photons'}],
        'power')

if TEST_LEVEL > 2:
    CLASS_INPUT['Isocurvature_modes'] = (
        [{'ic': 'ad,nid,cdi', 'c_ad_cdi': -0.5}],
        'normal')

    CLASS_INPUT['Scalar_field'] = (
        [{'Omega_scf': 0.1, 'attractor_ic_scf': 'yes',
        'scf_parameters': '10, 0, 0, 0'}],
        'normal')

    CLASS_INPUT['Inflation'] = (
        [{'P_k_ini type': 'inflation_V'},
        {'P_k_ini type': 'inflation_H'}],
        'normal')
    # CLASS_INPUT['Inflation'] = (
    #     [{'P_k_ini type': 'inflation_V'},
    #     {'P_k_ini type': 'inflation_H'},
    #     {'P_k_ini type': 'inflation_V_end'}],
    #     'normal')

if POWER_ALL:
    for k, v in iteritems(CLASS_INPUT):
        models, state = v
        CLASS_INPUT[k] = (models, 'power')

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

def custom_name_func(testcase_func, param_num, param):
    special_keys = ['N_ncdm']
    somekeys = []
    for key in param.args[0].keys():
        if key in special_keys:
            somekeys.append(key)
        elif 'mega' in key:
            somekeys.append(key)

    res = '{}_{:04d}_{}'.format(
        testcase_func.__name__,
        param_num,
        parameterized.to_safe_name('_'.join(somekeys))
    )
    return res.strip('_')

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
    def setUpClass(cls):
        cls.faulty_figs_path = os.path.join(
            os.path.sep.join(os.path.realpath(__file__).split(
                os.path.sep)[:-1]),
            'faulty_figs')

        if os.path.isdir(cls.faulty_figs_path):
            shutil.rmtree(cls.faulty_figs_path)

        os.mkdir(cls.faulty_figs_path)

    @classmethod
    def tearDownClass(cls):
        pass

    def setUp(self):
        """
        set up data used in the tests.
        setUp is called before each test function execution.
        """
        self.cosmo = Class()
        self.cosmo_newt = Class()

        if CLASS_VERBOSE:
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
                'distortions_verbose': 1,
                'output_verbose': 1,
            }
        else:
            self.verbose = {}
        self.scenario = {}

    def tearDown(self):
        self.cosmo.struct_cleanup()
        self.cosmo.empty()
        self.cosmo = 0
        self.cosmo_newt.struct_cleanup()
        self.cosmo_newt.empty()
        self.cosmo_newt = 0
        del self.scenario

    def poormansname(self, somedict):
        string = "_".join(
            [k+'='+str(v)
             for k, v in list(somedict.items())])
        string = string.replace('/', '%')
        string = string.replace(',', '')
        string = string.replace(' ', '')
        return string

    @parameterized.expand(TUPLE_ARRAY, doc_func=custom_name_func, custom_name_func=custom_name_func)
    @attr('dump_ini_files')
    def test_Valgrind(self, inputdict):
        """Dump files"""
        self.scenario.update(inputdict)
        self.name = self._testMethodName
        if self.has_incompatible_input():
            return
        path = os.path.join(self.faulty_figs_path, self.name)
        self.store_ini_file(path)
        self.scenario.update({'gauge':'Newtonian'})
        self.store_ini_file(path + 'N')

    @parameterized.expand(TUPLE_ARRAY, doc_func=custom_name_func, custom_name_func=custom_name_func)
    @attr('test_scenario')
    def test_scenario(self, inputdict):
        """Test scenario"""
        self.scenario.update(inputdict)
        self.name = self._testMethodName
        self.cosmo.set(dict(itertools.chain(self.verbose.items(), self.scenario.items())))

        cl_dict = {
            'tCl': ['tt'],
            'lCl': ['pp'],
            'pCl': ['ee', 'bb'],
            'nCl': ['dd'],
            'sCl': ['ll'],
        }

        # 'lensing' is always set to yes. Therefore, trying to compute 'tCl' or
        # 'pCl' will fail except if we also ask for 'lCl'.
        if self.has_incompatible_input():
            self.assertRaises(CosmoSevereError, self.cosmo.compute)
            return
        else:
            self.cosmo.compute()

        self.assertTrue(
            self.cosmo.state,
            "Class failed to go through all __init__ methods")
        # Depending
        if 'output' in self.scenario.keys():
            # Positive tests of raw cls
            output = self.scenario['output']
            for elem in output.split():
                if elem in cl_dict.keys():
                    for cl_type in cl_dict[elem]:
                        is_density_cl = (elem == 'nCl' or elem == 'sCl')
                        if is_density_cl:
                            cl = self.cosmo.density_cl(100)
                        else:
                            cl = self.cosmo.raw_cl(100)
                        self.assertIsNotNone(cl, "raw_cl returned nothing")
                        cl_length = np.shape(cl[cl_type][0])[0] if is_density_cl else np.shape(cl[cl_type])[0]
                        self.assertEqual(cl_length, 101, "raw_cl returned wrong size")
                if elem == 'mPk':
                    pk = self.cosmo.pk(0.1, 0)
                    self.assertIsNotNone(pk, "pk returned nothing")
            # Negative tests of output functions
            if not any([elem in list(cl_dict.keys()) for elem in output.split()]):
                # testing absence of any Cl
                self.assertRaises(CosmoSevereError, self.cosmo.raw_cl, 100)
            if 'mPk' not in output.split():
                # testing absence of mPk
                self.assertRaises(CosmoSevereError, self.cosmo.pk, 0.1, 0)

        if COMPARE_OUTPUT_REF or COMPARE_OUTPUT_GAUGE:
            # Now compute same scenario in Newtonian gauge
            self.cosmo_newt.set(
                dict(list(self.verbose.items())+list(self.scenario.items())))
            self.cosmo_newt.set({'gauge': 'newtonian'})
            self.cosmo_newt.compute()

        if COMPARE_OUTPUT_GAUGE:
            # Compare synchronous and Newtonian gauge
            self.assertTrue(
                self.cosmo_newt.state,
                "Class failed to go through all __init__ methods in Newtonian gauge")

            self.compare_output(self.cosmo, "Synchronous", self.cosmo_newt, 'Newtonian', COMPARE_CL_RELATIVE_ERROR_GAUGE, COMPARE_PK_RELATIVE_ERROR_GAUGE)

        if COMPARE_OUTPUT_REF:
            # Compute reference models in both gauges and compare
            cosmo_ref = classyref.Class()
            cosmo_ref.set(self.cosmo.pars)
            cosmo_ref.compute()
            status = self.compare_output(cosmo_ref, "Reference", self.cosmo, 'Synchronous', COMPARE_CL_RELATIVE_ERROR, COMPARE_PK_RELATIVE_ERROR)
            assert status, 'Reference comparison failed in Synchronous gauge!'

            cosmo_ref = classyref.Class()
            cosmo_ref.set(self.cosmo_newt.pars)
            cosmo_ref.compute()
            self.compare_output(cosmo_ref, "Reference", self.cosmo_newt, 'Newtonian', COMPARE_CL_RELATIVE_ERROR, COMPARE_PK_RELATIVE_ERROR)
            assert status, 'Reference comparison failed in Newtonian gauge!'

    def has_incompatible_input(self):

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

        # If we ask for Cl's of lensing potential, number counts or cosmic shear, we must have scalar modes.
        # The same applies to density and velocity transfer functions and the matter power spectrum:
        if 'output' in self.scenario and 'modes' in self.scenario and self.scenario['modes'].find('s') == -1:
            requested_output_types = set(self.scenario['output'].split())
            for scalar_output_type in ['lCl', 'nCl', 'dCl', 'sCl', 'mPk', 'dTk', 'mTk', 'vTk']:
                if scalar_output_type in requested_output_types:
                    should_fail = True
                    break

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

    def compare_output(self, reference, reference_name, candidate, candidate_name, rtol_cl, rtol_pk):
        status_pass = True
        for elem in ['raw_cl', 'lensed_cl']:
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
                    # For all self spectra, try to compare allclose
                    if key[0] == key[1]:
                        # If it is a 'dd' or 'll', it is a dictionary.
                        if isinstance(value, dict):
                            for subkey in list(value.keys()):
                                try:
                                    np.testing.assert_allclose(
                                        value[subkey],
                                        to_test[key][subkey],
                                        rtol=rtol_cl,
                                        atol=COMPARE_CL_ABSOLUTE_ERROR)
                                except AssertionError:
                                    self.cl_faulty_plot(elem + "_" + key, value[subkey][2:], reference_name, to_test[key][subkey][2:], candidate_name, rtol_cl)
                                except TypeError:
                                    self.cl_faulty_plot(elem + "_" + key, value[subkey][2:], reference_name, to_test[key][subkey][2:], candidate_name, rtol_cl)
                        else:
                            try:
                                np.testing.assert_allclose(
                                    value,
                                    to_test[key],
                                    rtol=rtol_cl,
                                    atol=COMPARE_CL_ABSOLUTE_ERROR)
                            except (AssertionError, TypeError) as e:
                                self.cl_faulty_plot(elem + "_" + key, value[2:], reference_name, to_test[key][2:], candidate_name, rtol_cl)
                                status_pass = False
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
                            self.cl_faulty_plot(elem + "_" + key, value[2:], reference_name, to_test[key][2:], candidate_name, rtol_cl)
                            status_pass = False

        if 'output' in list(self.scenario.keys()):
            if self.scenario['output'].find('mPk') != -1:
                # testing equality of Pk
                k = np.logspace(-2, log10(self.scenario['P_k_max_1/Mpc']), 50)
                reference_pk = np.array([reference.pk(elem, 0) for elem in k])
                candidate_pk = np.array([candidate.pk(elem, 0) for elem in k])
                try:
                    np.testing.assert_allclose(
                        reference_pk,
                        candidate_pk,
                        rtol=rtol_pk,
                        atol=COMPARE_PK_ABSOLUTE_ERROR)
                except AssertionError:
                    self.pk_faulty_plot(k, reference_pk, reference_name, candidate_pk, candidate_name, rtol_pk)
                    status_pass = False

        return status_pass

    def store_ini_file(self, path):
        parameters = dict(list(self.verbose.items()) + list(self.scenario.items()))
        with open(path + '.ini', 'w') as param_file:
            param_file.write('# ' + str(parameters) + '\n')
            if len(parameters) == 0:
                # CLASS complains if the .ini file does not do anything.
                param_file.write('write warnings = yes\n')
            for key, value in list(parameters.items()):
                param_file.write(key + " = " + str(value)+ '\n')

    def cl_faulty_plot(self, cl_type, reference, reference_name, candidate, candidate_name, rtol):
        path = os.path.join(self.faulty_figs_path, self.name)
        fig, axes = plt.subplots(2, 1, sharex=True)
        ell = np.arange(max(np.shape(candidate))) + 2
        factor = ell*(ell + 1)/(2*np.pi) if cl_type[-2:] != 'pp' else ell**5
        axes[0].plot(ell, factor*reference, label=reference_name)
        axes[0].plot(ell, factor*candidate, label=candidate_name)
        axes[1].semilogy(ell, 100*abs(candidate/reference - 1), label=cl_type)
        axes[1].axhline(y=100*rtol, color='k', ls='--')

        axes[-1].set_xlabel(r'$\ell$')
        if cl_type[-2:] == 'pp':
            axes[0].set_ylabel(r'$\ell^5 C_\ell^\mathrm{{{_cl_type}}}$'.format(_cl_type=cl_type[-2:].upper()))
        else:
            axes[0].set_ylabel(r'$\ell(\ell + 1)/(2\pi)C_\ell^\mathrm{{{_cl_type}}}$'.format(_cl_type=cl_type[-2:].upper()))
        axes[1].set_ylabel('Relative error [%]')

        for ax in axes:
            ax.legend(loc='upper right')

        fig.tight_layout()
        fname = '{}_{}_{}_vs_{}.pdf'.format(path, cl_type, reference_name, candidate_name)
        fig.savefig(fname, bbox_inches='tight')
        plt.close(fig)

        # Store parameters (contained in self.scenario) to text file
        self.store_ini_file(path)

    def pk_faulty_plot(self, k, reference, reference_name, candidate, candidate_name, rtol):
        path = os.path.join(self.faulty_figs_path, self.name)

        fig, axes = plt.subplots(2, 1, sharex=True)
        axes[0].loglog(k, k**1.5*reference, label=reference_name)
        axes[0].loglog(k, k**1.5*candidate, label=candidate_name)
        axes[0].legend(loc='upper right')

        axes[1].loglog(k, 100*np.abs(candidate/reference - 1))
        axes[1].axhline(y=100*rtol, color='k', ls='--')

        axes[-1].set_xlabel(r'$k\quad [\mathrm{Mpc}^{-1}]$')
        axes[0].set_ylabel(r'$k^\frac{3}{2}P(k)$')
        axes[1].set_ylabel(r'Relative error [%]')

        fig.tight_layout()
        fname = path + '_pk_{}_vs_{}.pdf'.format(reference_name, candidate_name)
        fig.savefig(fname, bbox_inches='tight')
        plt.close(fig)

        # Store parameters (contained in self.scenario) to text file
        self.store_ini_file(path)

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
