#
# python setup.py build_ext --inplace

import platform
import os
import numpy as np
from setuptools import Extension, setup
from Cython.Build import cythonize


def _embed_kernel(cl_name, hdr_name, sym_name, delim):
    """Embed a single OpenCL .cl file as a C++ raw-string header."""
    here = os.path.dirname(os.path.abspath(__file__))
    cl_path  = os.path.join(here, 'ttcr', cl_name)
    hdr_path = os.path.join(here, 'ttcr', hdr_name)

    with open(cl_path, 'r') as f:
        src = f.read()

    assert (')' + delim + '"') not in src, \
        f'{cl_name}: kernel source contains the raw-string delimiter; choose another'

    guard = hdr_name.replace('.', '_').replace('/', '_')
    content = (
        f'// Auto-generated from {cl_name} by setup.py.  Do not edit.\n'
        f'#ifndef ttcr_{guard}\n'
        f'#define ttcr_{guard}\n'
        'namespace ttcr {\n'
        f'static const char* const {sym_name} = R"{delim}(\n{src}){delim}";\n'
        '}\n'
        '#endif\n'
    )

    existing = None
    if os.path.exists(hdr_path):
        with open(hdr_path, 'r') as f:
            existing = f.read()

    if existing != content:
        with open(hdr_path, 'w') as f:
            f.write(content)


def generate_embedded_kernels():
    """Embed all OpenCL kernel sources into C++ headers.

    An installed wheel does not need the .cl files on disk at runtime;
    the source is compiled into the extension as a string constant.
    """
    _embed_kernel('Grid3Drn_kernels.cl', 'Grid3Drn_kernels_src.h',
                  'Grid3Drn_kernels_src', 'TTCR_CL_KERNEL')
    _embed_kernel('Grid2Drn_kernels.cl', 'Grid2Drn_kernels_src.h',
                  'Grid2Drn_kernels_src', 'TTCR_CL_2D')


generate_embedded_kernels()

# Per-platform build configuration.  All platforms compile the OpenCL grids, so
# every build needs the OpenCL headers and links against the ICD loader.
libraries = []
library_dirs = []
extra_include_dirs = []

if platform.system() == 'Darwin':
    # OpenCL ships as a system framework on macOS.
    extra_compile_args = ['-std=c++17', '-stdlib=libc++', '-O3', '-arch', 'x86_64',
                          '-arch', 'arm64',
                          '-Wno-unused-result', '-Wunreachable-code',
                          '-fno-common', '-dynamic', '-DNDEBUG', '-fwrapv',
                          '-pipe', '-isysroot/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk']
    extra_link_args = ['-bundle', '-undefined', 'dynamic_lookup', '-arch', 'x86_64',
                       '-arch', 'arm64', '-framework', 'OpenCL',
                       '-L/opt/local/lib', '-Wl,-headerpad_max_install_names',
                       '-Wl,-syslibroot,/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk']
elif platform.system() == 'Windows':
    extra_compile_args = ['/O2']
    extra_link_args = []
    libraries = ['OpenCL']
    # GitHub runners expose vcpkg through VCPKG_INSTALLATION_ROOT; honour either
    # that or a user-set VCPKG_ROOT so `vcpkg install opencl` is picked up.
    vcpkg = os.environ.get('VCPKG_INSTALLATION_ROOT') or os.environ.get('VCPKG_ROOT')
    if vcpkg:
        triplet = os.environ.get('VCPKG_DEFAULT_TRIPLET', 'x64-windows')
        extra_include_dirs.append(os.path.join(vcpkg, 'installed', triplet, 'include'))
        library_dirs.append(os.path.join(vcpkg, 'installed', triplet, 'lib'))
elif platform.system() == 'Linux':
    if 'readthedocs.org/user_builds/ttcrpy' in os.path.dirname(__file__):
        # do not optimize when building on readthedocs server
        extra_compile_args = ['-std=c++17', '-O0']
    else:
        extra_compile_args = ['-std=c++17', '-O3']
    extra_link_args = []
    libraries = ['OpenCL']

include_dirs = ['ttcr', 'boost_1_91_0', 'eigen-5.0.0', np.get_include()] + extra_include_dirs

setup(
    ext_modules=cythonize([
    Extension('ttcrpy.rgrid',
              sources=['src/ttcrpy/rgrid.pyx', 'src/ttcrpy/verbose.cpp'],  # additional source file(s)
              include_dirs=include_dirs,
              libraries=libraries,
              library_dirs=library_dirs,
              language='c++',             # generate C++ code
              extra_compile_args=extra_compile_args,
              extra_link_args=extra_link_args,
              ),
    Extension('ttcrpy.tmesh',
              sources=['src/ttcrpy/tmesh.pyx', 'src/ttcrpy/verbose.cpp'],  # additional source file(s)
              include_dirs=include_dirs,
              libraries=libraries,
              library_dirs=library_dirs,
              language='c++',             # generate C++ code
              extra_compile_args=extra_compile_args,
              extra_link_args=extra_link_args,
              ),
    ], language_level=3)
)
