
-----------------------------------------------------
Information for compiling mexfiles 
-----------------------------------------------------

You need:
- a C++ compiler that conforms to the C++11 standard
- the source codes of the ttcr package


On my OS X machine, I use this command to compile from a terminal:

MATLAB=/Applications/MATLAB_R2014a.app

$MATLAB/bin/mex -O CXXFLAGS='$CXXFLAGS -std=c++11 -stdlib=libc++' \
LDFLAGS='$LDFLAGS -std=c++11 -stdlib=libc++' \
-largeArrayDims -v -I$HOME/src/ttcr/ttcr \
-I$HOME/src/ttcr/boost_1_72_0 -I$HOME/src/ttcr/eigen-3.3.7 grid2dunsp_mex.cpp

3D classes must be compiled with verbose.cpp in the list of source files, i.e.

$MATLAB/bin/mex -O CXXFLAGS='$CXXFLAGS -std=c++11 -stdlib=libc++' \
LDFLAGS='$LDFLAGS -std=c++11 -stdlib=libc++' -largeArrayDims -v \
-I$HOME/src/ttcr/ttcr -I$HOME/src/ttcr/boost_1_72_0 \
-I$HOME/src/ttcr/eigen-3.3.7 grid3dunfs_mex.cpp verbose.cpp


On a windows 8.1 machine with intel compiler installed, I could compile it from the matlab prompt with:

mex -v -O COMPFLAGS='$COMPFLAGS /Qstd=c++11' -largeArrayDims -I../ttcr -I../boost_1_72_0 -I../eigen-3.3.7 grid2dunsp_mex.cpp


Unfortunately, I cannot offer extensive support for compiling on other platforms, especially windows variants.

Please report bugs to bernard.giroux@ete.inrs.ca

Enjoy!

-----

Bernard Giroux
INRS-ETE
