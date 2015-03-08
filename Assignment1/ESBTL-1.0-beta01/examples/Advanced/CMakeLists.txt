# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( test_example )

CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)
 
set (CMAKE_BUILD_TYPE Debug)
 
#~ find_package(BOOST)

#~ if ( BOOST_FOUND )
  include_directories(../../include)
  
  add_executable(reader reader.cpp)
  add_executable(coarse_grain coarse_grain.cpp)
  add_executable(solvation solvation.cpp)

  add_executable(to_xyzr_from_compressed to_xyzr_from_compressed.cpp)
  target_link_libraries(to_xyzr_from_compressed boost_iostreams)
  
  add_executable(create_pdbs_from_xtc create_pdbs_from_xtc.cpp)
  target_link_libraries(create_pdbs_from_xtc xdrfile)  


#~ else()
#~   
#~     message(STATUS "This program requires the BOOST library, and will not be compiled.")
#~ endif()

