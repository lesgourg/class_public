from setuptools import setup

DATA_FILES = [('DarkAges/transfer_functions/original',['transfer_functions/original/Transfer_Ch{:d}.dat'.format(i+1) for i in xrange(5)])]

setup(name='DarkAges',
      version='0.99',
      description='The darkAges-package',
      author='Patrick Stoecker',
      author_email='stoecker@physik.rwth-aachen.de',
      scripts = ['bin/DarkAges'],
      license='MIT',
      data_files = DATA_FILES,
      #include_package_data=True,
      packages=['DarkAges'])
