
# Run File as: $python setup.py build_ext --inline

from setuptools import setup, Extension
from Cython.Build import cythonize
import numpy as np


setup(
    setup_requires = ['setuptools>=18.0', 'cython'],
    ext_modules = cythonize(Extension('neighbours', sources=['neighbours.pyx'])),
    include_dirs = [np.get_include()]        
)
