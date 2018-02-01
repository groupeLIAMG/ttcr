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
              extra_compile_args=['-std=c++11','-stdlib=libc++','-O3'],
              ),
    Extension('cgrid3d',
              sources=['cgrid3d.pyx', 'Grid3Dttcr.cpp'],  # additional source file(s)
              include_dirs=['../ttcr/','/opt/local/include', np.get_include()],
              language='c++',             # generate C++ code
              extra_compile_args=['-std=c++11','-stdlib=libc++','-O3'],
              ),
    Extension('cmesh3d',
              sources=['cmesh3d.pyx', 'Mesh3Dttcr.cpp'],
              include_dirs=['../ttcr/','/opt/local/include/eigen3/',np.get_include()],
              language='c++',             # generate C++ code
              extra_compile_args=['-std=c++11','-stdlib=libc++','-O3'],
              ),
]

setup(
    name='ttcrpy',
    ext_modules=cythonize(extensions, include_path=['.',]),
)
