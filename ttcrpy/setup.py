#
# python setup.py build_ext --inplace

import platform
import numpy as np
from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize

if platform.system() == 'Windows':
    extra_compile_args = ['/O2']
    include_dirs = ['C:\\Users\\berna\OneDrive\Documents\\ttcr\\ttcr',
    'C:\\Users\\berna\OneDrive\Documents\\boost_1_66_0',
    'C:\\Users\\berna\OneDrive\Documents\eigen3',np.get_include()]
elif platform.system() == 'Darwin':
    extra_compile_args = ['-std=c++11','-stdlib=libc++','-O3']
    include_dirs = ['../ttcr/','/opt/local/include','/opt/local/include/eigen3/',np.get_include()]
elif platform.system() == 'Linux':
    extra_compile_args = ['-std=c++11','-O3']
    include_dirs = ['../ttcr/','/usr/include/eigen3/',np.get_include()]



extensions = [
    Extension('cgrid2d',
              sources=['cgrid2d.pyx', 'Grid2Dttcr.cpp'],  # additional source file(s)
              include_dirs=include_dirs,
              language='c++',             # generate C++ code
              extra_compile_args=extra_compile_args,
              ),
    Extension('cgrid3d',
              sources=['cgrid3d.pyx', 'Grid3Dttcr.cpp'],  # additional source file(s)
              include_dirs=include_dirs,
              language='c++',             # generate C++ code
              extra_compile_args=extra_compile_args,
              ),
    Extension('cmesh3d',
              sources=['cmesh3d.pyx', 'Mesh3Dttcr.cpp'],
              include_dirs=include_dirs,
              language='c++',             # generate C++ code
              extra_compile_args=extra_compile_args,
              ),
]

setup(
    name='ttcrpy',
    ext_modules=cythonize(extensions, include_path=['.',]),
)
