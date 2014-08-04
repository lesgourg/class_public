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

# Dictionary of models to test the wrapper against. Each of these scenario will
# be run against all the possible output choices (nothing, tCl, mPk, etc...),
# with or without non-linearities.
MODELS = {
    'LCDM': {},
    'Mnu': {'N_ncdm': 1, 'm_ncdm': 0.06},
    'Positive_Omega_k': {'Omega_k': 0.01},
    'Negative_Omega_k': {'Omega_k': -0.01},
    'Isocurvature_modes': {'ic': 'ad,nid,cdi', 'c_ad_cdi': -0.5},
    'Scalar_field': {'Omega_scf': 0.1, 'attractor_ic_scf': 'yes',
                     'scf_parameters': '10, 0, 0, 0'},
    }

# For extensive testing, change False to True. It will then append to MODELS
# all the possible combinations of models that are not incompatible. Obviously,
# this will be slower to test...
if False:
    INCOMPATIBLE_MODELS = ['LCDM', 'Positive_Omega_k', 'Negative_Omega_k']
    BASE_MODELS = [key for key in MODELS.keys()
                   if key not in INCOMPATIBLE_MODELS]

    ALL_MODELS = {}
    for n in range(0, len(BASE_MODELS)+1):
        for elem in itertools.combinations(BASE_MODELS, n):
            print elem
            temp = {}
            for key in elem:
                temp.update(MODELS[key])

            for key in INCOMPATIBLE_MODELS:
                ALL_MODELS['_'.join(elem+tuple([key], ))] = dict(
                    temp.items()+MODELS[key].items())
    MODELS = ALL_MODELS


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
        self.scenario = {'lensing': 'yes'}

    def tearDown(self):
        self.cosmo.struct_cleanup()
        self.cosmo.empty()
        self.cosmo_newt.struct_cleanup()
        self.cosmo_newt.empty()
        del self.scenario

    @parameterized.expand(
        itertools.product(
            MODELS.keys(),
            ({'output': ''}, {'output': 'mPk'}, {'output': 'tCl'},
             {'output': 'tCl pCl lCl'}, {'output': 'mPk tCl lCl', 'P_k_max_h/Mpc':10},
             {'output': 'nCl sCl'}, {'output': 'tCl pCl lCl nCl sCl'}),
            ({}, {'non linear': 'halofit'})))
    def test_wrapper_implementation(self, name, scenario, nonlinear):
        """Create a few instances based on different cosmologies"""
        self.scenario.update(MODELS[name])

        self.name = name
        self.nonlinear = nonlinear

        self.scenario.update(scenario)
        self.scenario.update(nonlinear)

        sys.stderr.write('\n\n---------------------------------\n')
        sys.stderr.write('| Test case %s |\n' % name)
        sys.stderr.write('---------------------------------\n')
        for key, value in self.scenario.iteritems():
            sys.stderr.write("%s = %s\n" % (key, value))
        sys.stderr.write("\n")

        setting = self.cosmo.set(
            dict(self.verbose.items()+self.scenario.items()))
        self.assertTrue(setting, "Class failed to initialize with input dict")

        cl_list = ['tCl', 'lCl', 'pCl', 'nCl', 'sCl']

        # Depending on the cases, the compute should fail or not
        should_fail = True
        output = self.scenario['output'].split()
        for elem in output:
            if elem in ['tCl', 'pCl']:
                for elem2 in output:
                    if elem2 == 'lCl':
                        should_fail = False
                        break

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
            # Positive tests
            output = self.scenario['output']
            for elem in output.split():
                if elem in cl_list:
                    print '--> testing raw_cl function'
                    cl = self.cosmo.raw_cl(100)
                    self.assertIsNotNone(cl, "raw_cl returned nothing")
                    self.assertEqual(
                        np.shape(cl['tt'])[0], 101,
                        "raw_cl returned wrong size")
                if elem == 'mPk':
                    print '--> testing pk function'
                    pk = self.cosmo.pk(0.1, 0)
                    self.assertIsNotNone(pk, "pk returned nothing")
            # Negative tests of output functions
            if not any([elem in cl_list for elem in output.split()]):
                print '--> testing absence of any Cl'
                self.assertRaises(CosmoSevereError, self.cosmo.raw_cl, 100)
            if 'mPk' not in self.scenario['output'].split():
                print '--> testing absence of mPk'
                #args = (0.1, 0)
                self.assertRaises(CosmoSevereError, self.cosmo.pk, 0.1, 0)

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

    @parameterized.expand(
        itertools.product(
            ('massless', 'massive', 'both'),
            ('photons', 'massless', 'exact'),
            ('t', 's, t')))
    def test_tensors(self, scenario, method, modes):
        """Test the new tensor mode implementation"""
        self.scenario = {}
        if scenario == 'massless':
            self.scenario.update({'N_eff': 3.046, 'N_ncdm':0})
        elif scenario == 'massive':
            self.scenario.update(
                {'N_eff': 0, 'N_ncdm': 2, 'm_ncdm': '0.03, 0.04',
                 'deg_ncdm': '2, 1'})
        elif scenario == 'both':
            self.scenario.update(
                {'N_eff': 1.5, 'N_ncdm': 2, 'm_ncdm': '0.03, 0.04',
                 'deg_ncdm': '1, 0.5'})

        sys.stderr.write('\n\n---------------------------------\n')
        sys.stderr.write('| Test case: %s %s %s |\n' % (
            scenario, method, modes))
        sys.stderr.write('---------------------------------\n')
        self.scenario.update({
            'tensor method': method, 'modes': modes,
            'output': 'tCl, pCl'})
        for key, value in self.scenario.iteritems():
            sys.stderr.write("%s = %s\n" % (key, value))
        sys.stderr.write("\n")
        self.cosmo.set(
            dict(self.verbose.items()+self.scenario.items()))
        self.cosmo.compute()

    @parameterized.expand(
        itertools.izip(
            powerset(['100*theta_s', 'Omega_dcdmdr', 'Omega_scf']),
            powerset([1.04, 0.20, -1]),))
    def test_shooting_method(self, variables, values):
        Omega_cdm = 0.25

        scenario = {'Omega_b': 0.05, }

        for variable, value in zip(variables, values):
            scenario.update({variable: value})

        if 'Omega_dcdmdr' in variables:
            scenario.update({
                'Gamma_dcdm': 100,
                'Omega_cdm': Omega_cdm-scenario['Omega_dcdmdr']})
        else:
            scenario.update({
                'Omega_cdm': Omega_cdm})

        if 'Omega_scf' in variables:
            scenario.update({
                'scf_tuning_index': 0,
                'Omega_Lambda': 0.6,
                'Omega_fld': 0,
                # TODO, the following values fail...
                #'scf_parameters': '11, -1.e-40, 0.1, 0',
                'scf_parameters': '10, 0, 0, 0',
                'YHe': 0.25,  # Needed for CLASS to compute
                'attractor_ic_scf': 'yes'})

        sys.stderr.write('\n\n---------------------------------\n')
        sys.stderr.write('| Test shooting: %s |\n' % (
            ', '.join(variables)))
        sys.stderr.write('---------------------------------\n')
        for key, value in scenario.iteritems():
            sys.stderr.write("%s = %s\n" % (key, value))
        sys.stderr.write("\n")

        scenario.update(self.verbose)
        self.assertTrue(
            self.cosmo.set(scenario),
            "Class failed to initialise with this input")
        self.assertRaises
        self.cosmo.compute()

        # Now, check that the values are properly extracted
        for variable, value in zip(variables, values):
            if variable == '100*theta_s':
                computed_value = self.cosmo.get_current_derived_parameters(
                    [variable])[variable]
                self.assertAlmostEqual(
                    value, computed_value, places=5)
            # In a perfect world, we should here compare the output also for
            # scf and dcdm_dr, by looking into the background_table array of
            # the background module

    def compare_output(self, reference, candidate):
        sys.stderr.write('\n\n---------------------------------\n')
        sys.stderr.write('| Comparing synch and Newt: |\n')
        sys.stderr.write('---------------------------------\n')

        for elem in ['raw_cl', 'lensed_cl', 'lensed_density_cl']:
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
        name = self.name+"_"
        if self.nonlinear:
            name += "NL_"
        name += cl_type

        path = os.path.join(self.faulty_figs_path, name)

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

        ax_lin.set_title(name)
        ax_log.set_title(name)

        ax_lin.legend([cl_type])
        ax_log.legend([cl_type])

        fig.savefig(path+'.pdf')

        # Store parameters (contained in self.scenario) to text file
        parameters = dict(self.verbose.items()+self.scenario.items())
        with open(path+'.ini', 'w') as param_file:
            for key, value in parameters.iteritems():
                param_file.write(key+" = "+str(value)+'\n')


if __name__ == '__main__':
    unittest.main()
