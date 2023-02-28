#
# python setup.py build_ext --inplace

import platform
import os
import subprocess
import numpy as np
from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext
from Cython.Build import cythonize

if platform.system() == 'Darwin':
    extra_compile_args = ['-std=c++11', '-stdlib=libc++', '-O3', '-arch', 'x86_64',
                          '-arch', 'arm64',
                          '-Wno-unused-result', '-Wunreachable-code',
                          '-fno-common', '-dynamic', '-DNDEBUG', '-fwrapv',
                          '-pipe', '-isysroot/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk']
    extra_link_args = ['-bundle', '-undefined', 'dynamic_lookup','-arch', 'x86_64',
                       '-arch', 'arm64',
                       '-L/opt/local/lib', '-Wl,-headerpad_max_install_names',
                       '-Wl,-syslibroot,/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk']
elif platform.system() == 'Windows':
    extra_compile_args = ['/O2']
    extra_link_args = []
elif platform.system() == 'Linux':
    if 'readthedocs.org/user_builds/ttcrpy' in os.path.dirname(__file__):
        # do not optimize when building on readthedocs server
        extra_compile_args = ['-std=c++11', '-O0']
    else:
        extra_compile_args = ['-std=c++11', '-O3']
    extra_link_args = []

include_dirs = ['ttcr', 'boost_1_81_0', 'eigen-3.4.0', np.get_include()]

setup(
    ext_modules=cythonize([
    Extension('ttcrpy.rgrid',
              sources=['src/ttcrpy/rgrid.pyx', 'src/ttcrpy/verbose.cpp'],  # additional source file(s)
              include_dirs=include_dirs,
              language='c++',             # generate C++ code
              extra_compile_args=extra_compile_args,
              extra_link_args=extra_link_args,
              ),
    Extension('ttcrpy.tmesh',
              sources=['src/ttcrpy/tmesh.pyx', 'src/ttcrpy/verbose.cpp'],  # additional source file(s)
              include_dirs=include_dirs,
              language='c++',             # generate C++ code
              extra_compile_args=extra_compile_args,
              extra_link_args=extra_link_args,
              ),
    ], language_level=3)
)
