#
# python setup.py build_ext --inplace

import numpy as np
from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize


extensions = [
    Extension('cgrid2d',
              sources=['cgrid2d.pyx', 'Grid2Dttcr.cpp'],  # additional source file(s)
              include_dirs=['../ttcr/',np.get_include()],
              language='c++',             # generate C++ code
              extra_compile_args=['-std=c++11', "-mmacosx-version-min=10.9"],
              extra_link_args = ["-stdlib=libc++", "-mmacosx-version-min=10.9"],),

    Extension('cmesh3d',
              sources=['cmesh3d.pyx', 'Mesh3Dttcr.cpp'],
              include_dirs=['../ttcr/','/opt/local/include/eigen3/',np.get_include()],
              language='c++',             # generate C++ code
              extra_compile_args=['-std=c++11', "-mmacosx-version-min=10.9"],
              extra_link_args = ["-stdlib=libc++", "-mmacosx-version-min=10.9"],),
]

setup(
    name='ttcrpy',
    ext_modules=cythonize(extensions, include_path=['.',]),
)
