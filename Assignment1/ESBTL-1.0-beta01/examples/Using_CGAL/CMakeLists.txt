# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( examples_example )

CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)
 
find_package(CGAL QUIET COMPONENTS Core )

if ( CGAL_FOUND )

  include( ${CGAL_USE_FILE} )

  include( CGAL_CreateSingleSourceCGALProgram )

  include_directories (BEFORE ../../include)
#  include_directories (BEFORE /home/sloriot/CGAL/install/trunk/Skin_surface_3/include)

  create_single_source_cgal_program( "triangulation.cpp" )
  create_single_source_cgal_program( "coarse_grain_triang.cpp" )
  
  link_libraries(boost_regex reader)

else()
  
    message(STATUS "This program requires the CGAL library, and will not be compiled.")
  
endif()

