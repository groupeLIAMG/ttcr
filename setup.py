#
# python setup.py build_ext --inplace

import glob
import platform
import os
import numpy as np
from setuptools import setup
from setuptools.extension import Extension
from Cython.Build import cythonize

if platform.system() == 'Darwin':
    os.environ['CC'] = 'clang'
    os.environ['CXX'] = 'clang++'
    os.environ['MACOSX_DEPLOYMENT_TARGET'] = '10.15'
    extra_compile_args = ['-std=c++17', '-stdlib=libc++', '-O3']
elif platform.system() == 'Windows':
    extra_compile_args = ['/O2']
elif platform.system() == 'Linux':
    extra_compile_args = ['-std=c++17', '-O3']

include_dirs = ['ttcr', 'boost_1_72_0', 'eigen-3.3.7', np.get_include()]


extensions = [
    Extension('ttcrpy.rgrid',
              sources=['ttcrpy/rgrid.pyx', 'ttcrpy/verbose.cpp'],  # additional source file(s)
              include_dirs=include_dirs,
              language='c++',             # generate C++ code
              extra_compile_args=extra_compile_args,
              ),
    Extension('ttcrpy.cgrid2d',
              sources=['ttcrpy/cgrid2d.pyx', 'ttcrpy/Grid2Dttcr.cpp'],  # additional source file(s)
              include_dirs=include_dirs,
              language='c++',             # generate C++ code
              extra_compile_args=extra_compile_args,
              ),
    Extension('ttcrpy.cgrid3d',
              sources=['ttcrpy/cgrid3d.pyx'],  # additional source file(s)
              include_dirs=include_dirs,
              language='c++',             # generate C++ code
              extra_compile_args=extra_compile_args,
              ),
    Extension('ttcrpy.cmesh2d',
              sources=['ttcrpy/cmesh2d.pyx'],
              include_dirs=include_dirs,
              language='c++',             # generate C++ code
              extra_compile_args=extra_compile_args,
              ),
    Extension('ttcrpy.cmesh3d',
              sources=['ttcrpy/cmesh3d.pyx', './ttcrpy/Mesh3Dttcr.cpp'],
              include_dirs=include_dirs,
              language='c++',             # generate C++ code
              extra_compile_args=extra_compile_args,
              ),
]

setup(
    name='ttcrpy',
    description='Code to perform raytracing for geophysical applications',
    author='Bernard Giroux',
    author_email='bernard.giroux@ete.inrs.ca',
    url='https://github.com/groupeLIAMG/ttcr',
    ext_modules=cythonize(extensions, language_level=3),
)
