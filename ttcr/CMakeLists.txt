cmake_minimum_required(VERSION 2.8.7)
project( ttcr )

set( VTK_DIR /usr/local/VTK/lib/cmake/vtk-7.0 )
set( EIGEN3_INCLUDE_DIR /usr/local/eigen/include/eigen3 )
set( BOOST_ROOT /usr/local/boost )

add_definitions(-DVTK)

FIND_LIBRARY(C++_LIBRARY stdc++ /usr/local/gcc/lib)

find_package(VTK 7.0 REQUIRED NO_MODULE)

set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)
find_package(Boost 1.60.0 REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

IF( NOT EIGEN3_INCLUDE_DIR )
MESSAGE( FATAL_ERROR "Please point the environment variable EIGEN3_INCLUDE_DIR to the include directory of your Eigen3 installation.")
ENDIF()
INCLUDE_DIRECTORIES ( "${EIGEN3_INCLUDE_DIR}" )

set(CMAKE_CXX_FLAGS "-std=c++11 -march=native")

set( CMAKE_INSTALL_PREFIX $ENV{HOME} )
set( ttcr3d_SRCS ttcr3d.cpp ttcr_io.cpp )
set( ttcr2d_SRCS ttcr2d.cpp ttcr_io.cpp )
set( ttcr2ds_SRCS ttcr2ds.cpp ttcr_io.cpp )
set( CMAKE_VERBOSE_MAKEFILE on )

include(${VTK_USE_FILE})

add_executable( ttcr3d ${ttcr3d_SRCS} )
add_executable( ttcr2d ${ttcr2d_SRCS} )
add_executable( ttcr2ds ${ttcr2ds_SRCS} )

target_link_libraries(ttcr3d ${VTK_LIBRARIES} ${C++_LIBRARY})
target_link_libraries(ttcr2d ${VTK_LIBRARIES} ${C++_LIBRARY})
target_link_libraries(ttcr2ds ${VTK_LIBRARIES} ${C++_LIBRARY})

set_property(TARGET ttcr3d ttcr2d ttcr2ds PROPERTY INSTALL_RPATH_USE_LINK_PATH TRUE)

install(TARGETS ttcr3d ttcr2d ttcr2ds RUNTIME DESTINATION bin)
