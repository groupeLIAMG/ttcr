[![Logo](https://github.com/groupeLIAMG/ttcr/blob/master/images/ttcrpy_logo.svg)](https://github.com/groupeLIAMG/ttcr)
============================

[![pypi](https://img.shields.io/pypi/v/ttcrpy.svg)](https://pypi.org/project/ttcrpy/)
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](./01_LICENSE.txt)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1162725.svg)](https://doi.org/10.5281/zenodo.1162725)
[![Build Status](https://travis-ci.com/groupeLIAMG/ttcr.svg?branch=master)](https://travis-ci.com/groupeLIAMG/ttcr)
[![Documentation Status](https://readthedocs.org/projects/ttcrpy/badge/?version=latest)](https://ttcrpy.readthedocs.io/en/latest/?badge=latest)


This repo contains C++ and python codes for raytracing on regular and unstructured meshes.
Matlab wrappers are provided as well.

- [Python package](#heading)
- [Stand alone command-line programs](#heading)
- [Matlab mex files](#heading)
- [References](#heading)

<!-- toc -->

## Python package

ttcrpy is a package for computing traveltimes and raytracing that was
developed with geophysical applications in mind, e.g. ray-based seismic/GPR
tomography and microseismic event location (joint hypocenter-velocity
inversion).  The package contains code to perform computation on 2D and 3D
rectilinear grids, as well as 2D triangular and 3D tetrahedral meshes. Three
different algorithms have been implemented: the Fast-Sweeping Method, the
Shortest-Path Method, and the Dynamic Shortest-Path Method.  Calculations can
be run in parallel on a multi-core machine.

The core computing code is written in C++, and has been wrapped with cython.

Documentation can be found on [Read The Docs](https://ttcrpy.readthedocs.io/)

## Stand-alone command-line programs

There are three programs that can be called from the command line:

- ttcr2d : raytracing on planar 2D meshes
- ttcr2ds : raytracing on undulated surfaces
- ttcr3d : raytracing in 3D

See [documentation](https://github.com/groupeLIAMG/ttcr/blob/master/docs/command_line.md) for
command-line programs options and file formats.

### Examples

Look at the files in the examples directory for some samples.

### Compiling

The programs are coded in C++ and follow the C++11 standard.  You must have VTK
(http://vtk.org) installed on your system, to benefit from full functionalities.
Files from the eigen3 (http://eigen.tuxfamily.org) and boost
(http://www.boost.org) libraries are distributed with the source to facilitate
compilation.  These codes were compiled and tested on macs with the default
compiler (clang).  They were also tested to some extent under linux with g++
version 4.8.

## Matlab wrappers

To compile the mexfiles, you will need:
- a C++ compiler that conforms to the C++11 standard
- the source codes of the ttcr package

On my OS X machine, I use this command to compile from a terminal:
```
MATLAB=/Applications/MATLAB_R2014a.app

$MATLAB/bin/mex -O CXXFLAGS='$CXXFLAGS -std=c++11 -stdlib=libc++' \
LDFLAGS='$LDFLAGS -std=c++11 -stdlib=libc++' -largeArrayDims -v \
-I$HOME/src/ttcr/ttcr -I$HOME/src/ttcr/boost_1_72_0 \
-I$HOME/src/ttcr/eigen-3.3.7 grid2dunsp_mex.cpp
```

3D classes must be compiled with `verbose.cpp` in the list of source files, i.e.
```
$MATLAB/bin/mex -O CXXFLAGS='$CXXFLAGS -std=c++11 -stdlib=libc++' \
LDFLAGS='$LDFLAGS -std=c++11 -stdlib=libc++' -largeArrayDims -v \
-I$HOME/src/ttcr/ttcr -I$HOME/src/ttcr/boost_1_72_0 \
-I$HOME/src/ttcr/eigen-3.3.7 grid3dunfs_mex.cpp verbose.cpp
```

On a windows machine with intel compiler installed, I could compile it from the matlab prompt with:
```
mex -v -O COMPFLAGS='$COMPFLAGS /Qstd=c++11' -largeArrayDims -I../ttcr -I../boost_1_72_0 -I../eigen-3.3.7 grid2dunsp_mex.cpp
```

Unfortunately, I cannot offer extensive support for compiling on other platforms, especially windows variants.


Please report bugs to bernard dot giroux at ete.inrs.ca




## References
```

@article{doi:10.1111/1365-2478.12930,
  author = {Nasr, Maher and Giroux, Bernard and Dupuis, J. Christian},
  title = {A hybrid approach to compute seismic travel times in three-dimensional tetrahedral meshes},
  journal = {Geophysical Prospecting},
  volume = {n/a},
  number = {n/a},
  pages = {},
  keywords = {Travel time, Seismic modelling, Ray tracing, Seismics, Computing aspects},
  doi = {10.1111/1365-2478.12930},
  url = {https://onlinelibrary.wiley.com/doi/abs/10.1111/1365-2478.12930},
  eprint = {https://onlinelibrary.wiley.com/doi/pdf/10.1111/1365-2478.12930},
}

@inbook{nasr18,
  author = { Maher Nasr  and  Bernard Giroux  and  J. Christian Dupuis },
  title = {An optimized approach to compute traveltimes in 3D unstructured meshes},
  booktitle = {SEG Technical Program Expanded Abstracts 2018},
  chapter = {},
  pages = {5073-5077},
  year = {2018},
  doi = {10.1190/segam2018-2997918.1},
  URL = {https://library.seg.org/doi/abs/10.1190/segam2018-2997918.1},
  eprint = {https://library.seg.org/doi/pdf/10.1190/segam2018-2997918.1}
}

@InProceedings{giroux14,
  Title = {Comparison of grid-based methods for raytracing on unstructured meshes},
  Author = {Bernard Giroux},
  Booktitle = {SEG Technical Program Expanded Abstracts},
  Year = {2014},
  Pages = {3388-3392},
  Chapter = {649},
  DOI = {10.1190/segam2014-1197.1},
  Eprint = {http://library.seg.org/doi/pdf/10.1190/segam2014-1197.1},
  URL = {http://dx.doi.org/10.1190/segam2014-1197.1}
}


@ARTICLE{giroux13,
  author = {Bernard Giroux and Beno\^{\i}t Larouche},
  title = {Task-parallel implementation of {3D} shortest path raytracing for
	geophysical applications},
  journal = {Computers & Geosciences},
  year = {2013},
  volume = {54},
  pages = {130--141},
  number = {0},
  doi = {10.1016/j.cageo.2012.12.005}
  url = {http://dx.doi.org/10.1016/j.cageo.2012.12.005}
}

@INPROCEEDINGS{giroux13b,
  author = {Bernard Giroux},
  title = {Shortest path raytracing on tetrahedral meshes},
  booktitle = {75$^{th}$ EAGE Conference \& Exhibition},
  year = {2013},
  address = {London},
  organization = {EAGE},
  doi = {10.3997/2214-4609.20130236}
  url = {http://dx.doi.org/10.3997/2214-4609.20130236}
}

@ARTICLE{lelievre11,
  author = {Leli\`evre, Peter G. and Farquharson, Colin G. and Hurich, Charles A.},
  title = {Computing first-arrival seismic traveltimes on unstructured 3-{D}
	tetrahedral grids using the Fast Marching Method},
  journal = {Geophysical Journal International},
  year = {2011},
  volume = {184},
  pages = {885-896},
  number = {2},
  doi = {10.1111/j.1365-246X.2010.04880.x}
  url = {http://dx.doi.org/10.1111/j.1365-246X.2010.04880.x}
}

@ARTICLE{qian07,
  author = {Qian, Jianliang and Zhang, Yong-Tao and Zhao, Hong-Kai},
  title = {Fast Sweeping Methods for Eikonal Equations on Triangular Meshes},
  journal = {SIAM Journal on Numerical Analysis},
  year = {2007},
  volume = {45},
  pages = {83--107},
  number = {1},
  doi = {10.1137/050627083},
  publisher = {Society for Industrial and Applied Mathematics},
  url = {http://www.jstor.org/stable/40232919}
}

@Article{zhang06,
  Title                    = {High Order Fast Sweeping Methods for Static {H}amilton–{J}acobi Equations},
  Author                   = {Yong-Tao Zhang and Hong-Kai Zhao and Jianliang Qian},
  Journal                  = {Journal of Scientific Computing},
  Year                     = {2006},
  Number                   = {1},
  Pages                    = {25--56},
  Volume                   = {29},
  DOI                      = {10.1007/s10915-005-9014-3},
  URL                      = {http://dx.doi.org/10.1007/s10915-005-9014-3}
}

@Article{zhao05,
  Title                    = {A Fast Sweeping Method for Eikonal Equations},
  Author                   = {Zhao, Hongkai},
  Journal                  = {Mathematics of Computation},
  Year                     = {2005},
  Month                    = apr,
  Number                   = {250},
  Pages                    = {603--627},
  Volume                   = {74},
  Publisher                = {American Mathematical Society},
  URL                      = {http://www.jstor.org/stable/4100081}
}
```
