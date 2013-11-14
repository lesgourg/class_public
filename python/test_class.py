from classy import Class
import itertools
import numpy as np
import unittest
from nose_parameterized import parameterized


class TestClass(unittest.TestCase):
    """
    Testing Class and its wrapper classy on different cosmologies

    To run it, do
    ~] nosetest test_class.py

    It will run many times Class, on different cosmological scenarios, and
    everytime testing for different output possibilities (none asked, only mPk,
    etc..)

    """

    def setUp(self):
        """
        set up data used in the tests.
        setUp is called before each test function execution.
        """
        self.cosmo = Class()

        self.verbose = {
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
        self.cosmo.cleanup()
        del self.scenario

             #'Mnu',
             #'Positive_Omega_k',
             #'Negative_Omega_k'),
    #@parameterized.expand([
    @parameterized.expand(
        itertools.product(
            ('LCDM',
             'Mnu',
             'Positive_Omega_k',
             'Negative_Omega_k'),
            ({}, {'output': 'mPk'}, {'output': 'tCl'},
             {'output': 'tCl lCl'}, {'output': 'mPk tCl lCl'})))
    def test_parameters(self, name, scenario):
        """Create a few instances based on different cosmologies"""
        if name == 'Mnu':
            self.scenario.update({'N_ncdm': 1, 'm_ncdm': 0.06})
        elif name == 'Positive_Omega_k':
            self.scenario.update({'Omega_k': 0.01})
        elif name == 'Negative_Omega_k':
            self.scenario.update({'Omega_k': -0.01})

        self.scenario.update(scenario)

        print '\n\n--------------------------'
        print '| Test case %s |' % name
        print '--------------------------'
        for key, value in self.scenario.iteritems():
            print key, '=', value
        print

        setting = self.cosmo.set(
            dict(self.verbose.items()+self.scenario.items()))
        self.assertTrue(setting, "Class failed to initialize with input dict")

        self.cosmo.compute()
        self.assertTrue(
            self.cosmo.state,
            "Class failed to go through all __init__ methods")
        if self.cosmo.state:
            print '--> Class is ready'
        # Depending
        if 'output' in self.scenario.keys():
            # Positive tests
            for elem in self.scenario['output'].split():
                if elem == 'tCl':
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
            if 'tCl' not in self.scenario['output'].split():
                print '--> testing absence of tCl'
                self.assertRaises(NameError, self.cosmo.raw_cl, 100)
            #if 'mPk' not in self.scenario['output'].split():
                #print '--> testing absence of mPk'
                ##args = (1, 0)
                #pk = self.cosmo.pk(0.1, 0.)
                #print pk
                #self.assertRaises(NameError, self.cosmo.pk, *args)

        print '~~~~~~~~ passed ? '

if __name__ == '__main__':
    unittest.main()
