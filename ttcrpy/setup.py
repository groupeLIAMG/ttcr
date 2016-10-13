#
# python setup.py build_ext --inplace

import numpy as np
from distutils.core import setup, Extension
from Cython.Build import cythonize

setup(ext_modules = cythonize(Extension(
                                        'cgrid2d',
                                        sources=['cgrid2d.pyx', 'Grid2Dttcr.cpp'],  # additional source file(s)
                                        include_dirs=['../ttcr/',np.get_include()],
                                        language='c++',             # generate C++ code
                                        extra_compile_args=['-std=c++11'],)))
