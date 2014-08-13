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

# To avoid testing for differences between synchronous and Newtonian gauge, set this flag to False
COMPARE_OUTPUT = False

# Dictionary of models to test the wrapper against. Each of these scenario will
# be run against all the possible output choices (nothing, tCl, mPk, etc...),
# with or without non-linearities.
# Never input the default value, as this will be **automatically** tested
# against. Indeed, when not specifying a field, CLASS takes the default input.
CLASS_INPUT = {}

CLASS_INPUT['Mnu'] = (
    [{'N_eff': 0.0, 'N_ncdm': 1, 'm_ncdm': 0.06, 'deg_ncdm': 3.0},
     {'N_eff': 1.5, 'N_ncdm': 1, 'm_ncdm': 0.03, 'deg_ncdm': 1.5}],
    'normal')

CLASS_INPUT['Curvature'] = (
    [{'Omega_k': 0.01},
     {'Omega_k': -0.01}],
    'normal')

CLASS_INPUT['Isocurvature_modes'] = (
    [{'ic': 'ad,nid,cdi', 'c_ad_cdi': -0.5}],
    'normal')

CLASS_INPUT['Scalar_field'] = (
    [{'Omega_scf': 0.1, 'attractor_ic_scf': 'yes',
      'scf_parameters': '10, 0, 0, 0'}],
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
    [{'output': 'mPk'},
     {'output': 'tCl'},
     {'output': 'tCl pCl lCl'},
     {'output': 'mPk tCl lCl', 'P_k_max_h/Mpc':10},
     {'output': 'nCl sCl'},
     {'output': 'tCl pCl lCl nCl sCl'}],
    'power')

CLASS_INPUT['Nonlinear'] = (
    [{'non linear': 'halofit'}],
    'power')

CLASS_INPUT['Lensing'] = (
    [{'lensing': 'yes'}],
    'power')

# Let's kill the machine (replace all 'normal' flags with power', uncomment at you own risk)
# for k, v in CLASS_INPUT.iteritems():
#     models, state = v
#     CLASS_INPUT[k] = (models, 'power')

INPUTPOWER = []
INPUTNORMAL = [{}]
for key, value in CLASS_INPUT.iteritems():
    models, state = value
    if state == 'power':
        INPUTPOWER.append([{}]+models)
    else:
        INPUTNORMAL.extend(models)

    PRODPOWER = list(itertools.product(*INPUTPOWER))

    DICTARRAY=[]
    for normelem in INPUTNORMAL:
        for powelem in PRODPOWER: #itertools.product(*modpower):
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
            'spectra_verbose': 1,
            'nonlinear_verbose': 1,
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
             for k, v in somedict.iteritems()])
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
        for key, value in self.scenario.iteritems():
            sys.stderr.write("%s = %s\n" % (key, value))
        sys.stderr.write("\n")

        setting = self.cosmo.set(
            dict(self.verbose.items()+self.scenario.items()))
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
            print '--> Class is ready'
        # Depending
        if 'output' in self.scenario.keys():
            # Positive tests of raw cls
            output = self.scenario['output']
            for elem in output.split():
                if elem in cl_dict.keys():
                    for cl_type in cl_dict[elem]:
                        sys.stderr.write('--> testing raw_cl for %s' % cl_type)
                        cl = self.cosmo.raw_cl(100)
                        self.assertIsNotNone(cl, "raw_cl returned nothing")
                        self.assertEqual(
                            np.shape(cl[cl_type])[0], 101,
                            "raw_cl returned wrong size")
                    # TODO do the same for lensed if 'lCl' is there, and for density cl
                if elem == 'mPk':
                    sys.stderr.write('--> testing pk function')
                    pk = self.cosmo.pk(0.1, 0)
                    self.assertIsNotNone(pk, "pk returned nothing")
            # Negative tests of output functions
            if not any([elem in cl_dict.keys() for elem in output.split()]):
                sys.stderr.write('--> testing absence of any Cl')
                self.assertRaises(CosmoSevereError, self.cosmo.raw_cl, 100)
            if 'mPk' not in output.split():
                sys.stderr.write('--> testing absence of mPk')
                #args = (0.1, 0)
                self.assertRaises(CosmoSevereError, self.cosmo.pk, 0.1, 0)

        if COMPARE_OUTPUT:
            # Now, compute with Newtonian gauge, and compare the results
            self.cosmo_newt.set(
                dict(self.verbose.items()+self.scenario.items()))
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
            if 'output' not in self.scenario.keys():
                should_fail = True
            else:
                output = self.scenario['output'].split()
                if 'tCl' not in output and 'pCl' not in output:
                    should_fail = True

        # If we have specified lensing, we must have lCl in output,
        # otherwise lensing will not be read (which is an error).
        if 'lensing' in self.scenario.keys():
            if 'output' not in self.scenario.keys():
                should_fail = True
            else:
                output = self.scenario['output'].split()
                if 'lCl' not in output:
                    should_fail = True
                elif 'tCl' not in output and 'pCl' not in output:
                    should_fail = True

        # If we have specified a tensor method, we must have tensors.
        if 'tensor method' in self.scenario.keys():
            if not has_tensor(self.scenario):
                should_fail = True


        # If we have specified non linear, we must have some form of
        # perturbations output.
        if 'non linear' in self.scenario.keys():
            if 'output' not in self.scenario.keys():
                should_fail = True


        # If we ask for Cl's of lensing potential, we must have scalar modes.
        if 'output' in self.scenario.keys() and 'lCl' in self.scenario['output'].split():
            if 'modes' in self.scenario.keys() and self.scenario['modes'].find('s') == -1:
                should_fail = True


        return should_fail

    # @parameterized.expand(
    #     itertools.product(
    #         ('massless', 'massive', 'both'),
    #         ('photons', 'massless', 'exact'),
    #         ('t', 's, t')))
    # def test_tensors(self, scenario, method, modes):
    #     """Test the new tensor mode implementation"""
    #     self.scenario = {}
    #     if scenario == 'massless':
    #         self.scenario.update({'N_eff': 3.046, 'N_ncdm':0})
    #     elif scenario == 'massive':
    #         self.scenario.update(
    #             {'N_eff': 0, 'N_ncdm': 2, 'm_ncdm': '0.03, 0.04',
    #              'deg_ncdm': '2, 1'})
    #     elif scenario == 'both':
    #         self.scenario.update(
    #             {'N_eff': 1.5, 'N_ncdm': 2, 'm_ncdm': '0.03, 0.04',
    #              'deg_ncdm': '1, 0.5'})

    #     sys.stderr.write('\n\n---------------------------------\n')
    #     sys.stderr.write('| Test case: %s %s %s |\n' % (
    #         scenario, method, modes))
    #     sys.stderr.write('---------------------------------\n')
    #     self.scenario.update({
    #         'tensor method': method, 'modes': modes,
    #         'output': 'tCl, pCl'})
    #     for key, value in self.scenario.iteritems():
    #         sys.stderr.write("%s = %s\n" % (key, value))
    #     sys.stderr.write("\n")
    #     self.cosmo.set(
    #         dict(self.verbose.items()+self.scenario.items()))
    #     self.cosmo.compute()

    # @parameterized.expand(
    #     itertools.izip(
    #         powerset(['100*theta_s', 'Omega_dcdmdr', 'Omega_scf']),
    #         powerset([1.04, 0.20, -1]),))
    # def test_shooting_method(self, variables, values):
    #     Omega_cdm = 0.25

    #     scenario = {'Omega_b': 0.05, }

    #     for variable, value in zip(variables, values):
    #         scenario.update({variable: value})

    #     if 'Omega_dcdmdr' in variables:
    #         scenario.update({
    #             'Gamma_dcdm': 100,
    #             'Omega_cdm': Omega_cdm-scenario['Omega_dcdmdr']})
    #     else:
    #         scenario.update({
    #             'Omega_cdm': Omega_cdm})

    #     if 'Omega_scf' in variables:
    #         scenario.update({
    #             'scf_tuning_index': 0,
    #             'Omega_Lambda': 0.6,
    #             'Omega_fld': 0,
    #             # TODO, the following values fail...
    #             #'scf_parameters': '11, -1.e-40, 0.1, 0',
    #             'scf_parameters': '10, 0, 0, 0',
    #             'YHe': 0.25,  # Needed for CLASS to compute
    #             'attractor_ic_scf': 'yes'})

    #     sys.stderr.write('\n\n---------------------------------\n')
    #     sys.stderr.write('| Test shooting: %s |\n' % (
    #         ', '.join(variables)))
    #     sys.stderr.write('---------------------------------\n')
    #     for key, value in scenario.iteritems():
    #         sys.stderr.write("%s = %s\n" % (key, value))
    #     sys.stderr.write("\n")

    #     scenario.update(self.verbose)
    #     self.assertTrue(
    #         self.cosmo.set(scenario),
    #         "Class failed to initialise with this input")
    #     self.assertRaises
    #     self.cosmo.compute()

    #     # Now, check that the values are properly extracted
    #     for variable, value in zip(variables, values):
    #         if variable == '100*theta_s':
    #             computed_value = self.cosmo.get_current_derived_parameters(
    #                 [variable])[variable]
    #             self.assertAlmostEqual(
    #                 value, computed_value, places=5)
    #         # In a perfect world, we should here compare the output also for
    #         # scf and dcdm_dr, by looking into the background_table array of
    #         # the background module

    def compare_output(self, reference, candidate):
        sys.stderr.write('\n\n---------------------------------\n')
        sys.stderr.write('| Comparing synch and Newt: |\n')
        sys.stderr.write('---------------------------------\n')

        for elem in ['raw_cl', 'lensed_cl', 'density_cl']:
            to_test = getattr(candidate, elem)()
            ref = getattr(reference, elem)()
            for key, value in ref.iteritems():
                print '--> testing equality of %s %s' % (elem, key)
                try:
                    np.testing.assert_allclose(
                        value, to_test[key], rtol=1e-03, atol=1e-20)
                except AssertionError:
                    self.cl_faulty_plot(elem+"_"+key,
                                        value[2:], to_test[key][2:])

        if 'output' in self.scenario.keys():
            if self.scenario['output'].find('mPk') != -1:
                print '--> testing equality of Pk'
                k = np.logspace(
                    -2, log10(self.scenario['P_k_max_h/Mpc']))
                reference_pk = np.array(
                    [reference.pk(elem, 0) for elem in k])
                candidate_pk = np.array(
                    [candidate.pk(elem, 0) for elem in k])
                np.testing.assert_allclose(
                    reference_pk, candidate_pk, rtol=1e-03, atol=1e-20)

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

        fig.savefig(path+'.pdf')

        # Store parameters (contained in self.scenario) to text file
        parameters = dict(self.verbose.items()+self.scenario.items())
        with open(path+'.ini', 'w') as param_file:
            for key, value in parameters.iteritems():
                param_file.write(key+" = "+str(value)+'\n')


def has_tensor(input_dict):
    if 'modes' in input_dict.keys():
        if input_dict['modes'].find('t') != -1:
            return True
    else:
        return False
    return False

if __name__ == '__main__':
    toto = TestClass()
    #toto.poormansname({'output': 'tCl, mPk', 'Omega_k': 0.02})
    unittest.main()
