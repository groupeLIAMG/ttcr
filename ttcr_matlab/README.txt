
----------------------------------------------------
Information for compiling mexfiles used by grd2duisp
----------------------------------------------------

You need:
- a C++ compiler that conforms to the C++11 standard
- to have the eigen3 library installed on your system
		  (http://eigen.tuxfamily.org)
- as well as the boost library (http://www.boost.org)
- the source codes of the ttcr package


On my OS X machine, I use this command to compile from a terminal:

MATLAB=/Applications/MATLAB_R2014a.app

$MATLAB/bin/mex -O CXXFLAGS='$CXXFLAGS -std=c++11 -stdlib=libc++' LDFLAGS='$LDFLAGS -std=c++11 -stdlib=libc++' -largeArrayDims -v -I$HOME/src/ttcr/ttcr -I/opt/local/include -I/opt/local/include/eigen3 grid2duisp_mex.cpp


On a windows 8.1 machine with intel compiler installed, I could compile it from the matlab prompt with:

mex -O COMPFLAGS='$COMPFLAGS /Qstd=c++11' -largeArrayDims -v -I../ttcr -I'C:\libraries\boost_1_55_0' -I'C:\Program Files (x86)\Eigen\include' grid2duisp_mex.cpp

I got windows binaries of the eigen library from http://pointclouds.org/downloads/windows.html


Unfortunately, I cannot offer extensive support for compiling on other platforms, especially windows variants.

Please report bugs to bernard.giroux@ete.inrs.ca

Enjoy!

-----

Bernard Giroux
INRS-ETE
2014-04-28
