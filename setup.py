#
# python setup.py build_ext --inplace

import platform
import os
import subprocess
import numpy as np
from setuptools import setup
from setuptools.command.build_ext import build_ext
from setuptools.extension import Extension
from Cython.Build import cythonize

class custom_build_ext(build_ext):
    def build_extensions(self):
        # Override the compiler executables. Importantly, this
        # removes the "default" compiler flags that would
        # otherwise get passed on to to the compiler, i.e.,
        # distutils.sysconfig.get_var("CFLAGS").
        # Done here to remove "-Wsign-compare"
        if platform.system() == 'Darwin':
            self.compiler.set_executable("compiler_so", "/usr/bin/clang")
            self.compiler.set_executable("compiler_cxx", "/usr/bin/clang++")
            self.compiler.set_executable("linker_so", "/usr/bin/clang++")
        build_ext.build_extensions(self)


with open("README.md", "r") as fh:
    long_description = fh.read()

if platform.system() == 'Darwin':
    extra_compile_args = ['-std=c++11', '-stdlib=libc++', '-O3',
                          '-Wno-unused-result', '-Wunreachable-code',
                          '-fno-common', '-dynamic', '-DNDEBUG', '-fwrapv',
                          '-pipe', '-isysroot/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk']
    extra_link_args = ['-bundle', '-undefined', 'dynamic_lookup',
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

include_dirs = ['ttcr', 'boost_1_72_0', 'eigen-3.3.7', np.get_include()]

extensions = [
    Extension('ttcrpy.rgrid',
              sources=['ttcrpy/rgrid.pyx', 'ttcrpy/verbose.cpp'],  # additional source file(s)
              include_dirs=include_dirs,
              language='c++',             # generate C++ code
              extra_compile_args=extra_compile_args,
              extra_link_args=extra_link_args,
              ),
    Extension('ttcrpy.tmesh',
              sources=['ttcrpy/tmesh.pyx', 'ttcrpy/verbose.cpp'],  # additional source file(s)
              include_dirs=include_dirs,
              language='c++',             # generate C++ code
              extra_compile_args=extra_compile_args,
              extra_link_args=extra_link_args,
              ),
]

setup(
    name='ttcrpy',
    use_scm_version=True,
    setup_requires=['setuptools_scm'],
    description='Code to perform raytracing for geophysical applications',
    long_description=long_description,
    long_description_content_type='text/markdown',
    author='Bernard Giroux',
    author_email='bernard.giroux@ete.inrs.ca',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Intended Audience :: Education',
        'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
        'Operating System :: OS Independent',
        'Programming Language :: C++',
        'Programming Language :: Cython',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Topic :: Scientific/Engineering'
        ],
    url='https://github.com/groupeLIAMG/ttcr',
    ext_modules=cythonize(extensions, language_level=3),
    package_data={
        # If any package contains *.txt or *.rst files, include them:
        "": ["docs/*.md"],
        },
    cmdclass={"build_ext": custom_build_ext}
)
