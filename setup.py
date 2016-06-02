#setup.py, like a python Makefile

#from distutils.core import setup
#from Cython.Build import cythonize
#import numpy as np
#
#setup(
#    ext_modules = cythonize("advection.pyx", include_path = [np.get_include()])
#)

from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import numpy as np

extensions = [
    Extension('advection', ['advection.pyx'], include_dirs=[np.get_include()]),
]

setup(
    ext_modules=cythonize(extensions)
)
