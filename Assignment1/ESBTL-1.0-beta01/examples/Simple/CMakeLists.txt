# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( example )

CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)

set (CMAKE_BUILD_TYPE Debug)
#SET(CMAKE_CXX_FLAGS "-std=c++0x")

#~ find_package(BOOST)

#~ if ( BOOST_FOUND )
  include_directories(../../include)
  
  add_executable(to_xyzr to_xyzr.cpp)
  add_executable(read_models read_models.cpp)
  add_executable(write_selection write_selection.cpp)
  add_executable(select_in_system select_in_system.cpp)
  add_executable(select_in_reader select_in_reader.cpp)
  add_executable(generic_select generic_select.cpp)
  
  

#~ else()
#~   
#~     message(STATUS "This program requires the BOOST library, and will not be compiled.")
#~ endif()

