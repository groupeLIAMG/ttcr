#
# python setup.py build_ext --inplace

import platform
import os
import numpy as np
from setuptools import Extension, setup
from Cython.Build import cythonize


def generate_embedded_kernels():
    """Embed the OpenCL kernel source into a C++ header.

    The OpenCL grids build their kernels from ttcr::Grid3Drn_kernels_src, which
    is generated here from ttcr/Grid3Drn_kernels.cl.  Embedding the source means
    an installed wheel does not need to locate the .cl file on disk at runtime.
    """
    here = os.path.dirname(os.path.abspath(__file__))
    cl_path = os.path.join(here, 'ttcr', 'Grid3Drn_kernels.cl')
    hdr_path = os.path.join(here, 'ttcr', 'Grid3Drn_kernels_src.h')

    with open(cl_path, 'r') as f:
        src = f.read()

    delim = 'TTCR_CL_KERNEL'
    assert (')' + delim + '"') not in src, \
        'kernel source contains the raw-string delimiter; choose another'

    content = (
        '// Auto-generated from Grid3Drn_kernels.cl by setup.py.  Do not edit.\n'
        '#ifndef ttcr_Grid3Drn_kernels_src_h\n'
        '#define ttcr_Grid3Drn_kernels_src_h\n'
        'namespace ttcr {\n'
        'static const char* const Grid3Drn_kernels_src = R"%s(\n%s)%s";\n'
        '}\n'
        '#endif\n'
    ) % (delim, src, delim)

    # Only rewrite when the content changes, to avoid forcing a rebuild.
    if not os.path.exists(hdr_path) or open(hdr_path).read() != content:
        with open(hdr_path, 'w') as f:
            f.write(content)


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
