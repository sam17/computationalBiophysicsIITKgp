# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list

# Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/dementor/FOSS/ESBTL-1.0-beta01/examples/Simple

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/dementor/FOSS/ESBTL-1.0-beta01/examples/Simple

# Include any dependencies generated for this target.
include CMakeFiles/generic_select.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/generic_select.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/generic_select.dir/flags.make

CMakeFiles/generic_select.dir/generic_select.o: CMakeFiles/generic_select.dir/flags.make
CMakeFiles/generic_select.dir/generic_select.o: generic_select.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/dementor/FOSS/ESBTL-1.0-beta01/examples/Simple/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/generic_select.dir/generic_select.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/generic_select.dir/generic_select.o -c /home/dementor/FOSS/ESBTL-1.0-beta01/examples/Simple/generic_select.cpp

CMakeFiles/generic_select.dir/generic_select.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/generic_select.dir/generic_select.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/dementor/FOSS/ESBTL-1.0-beta01/examples/Simple/generic_select.cpp > CMakeFiles/generic_select.dir/generic_select.i

CMakeFiles/generic_select.dir/generic_select.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/generic_select.dir/generic_select.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/dementor/FOSS/ESBTL-1.0-beta01/examples/Simple/generic_select.cpp -o CMakeFiles/generic_select.dir/generic_select.s

CMakeFiles/generic_select.dir/generic_select.o.requires:
.PHONY : CMakeFiles/generic_select.dir/generic_select.o.requires

CMakeFiles/generic_select.dir/generic_select.o.provides: CMakeFiles/generic_select.dir/generic_select.o.requires
	$(MAKE) -f CMakeFiles/generic_select.dir/build.make CMakeFiles/generic_select.dir/generic_select.o.provides.build
.PHONY : CMakeFiles/generic_select.dir/generic_select.o.provides

CMakeFiles/generic_select.dir/generic_select.o.provides.build: CMakeFiles/generic_select.dir/generic_select.o

# Object files for target generic_select
generic_select_OBJECTS = \
"CMakeFiles/generic_select.dir/generic_select.o"

# External object files for target generic_select
generic_select_EXTERNAL_OBJECTS =

generic_select: CMakeFiles/generic_select.dir/generic_select.o
generic_select: CMakeFiles/generic_select.dir/build.make
generic_select: CMakeFiles/generic_select.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable generic_select"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/generic_select.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/generic_select.dir/build: generic_select
.PHONY : CMakeFiles/generic_select.dir/build

CMakeFiles/generic_select.dir/requires: CMakeFiles/generic_select.dir/generic_select.o.requires
.PHONY : CMakeFiles/generic_select.dir/requires

CMakeFiles/generic_select.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/generic_select.dir/cmake_clean.cmake
.PHONY : CMakeFiles/generic_select.dir/clean

CMakeFiles/generic_select.dir/depend:
	cd /home/dementor/FOSS/ESBTL-1.0-beta01/examples/Simple && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/dementor/FOSS/ESBTL-1.0-beta01/examples/Simple /home/dementor/FOSS/ESBTL-1.0-beta01/examples/Simple /home/dementor/FOSS/ESBTL-1.0-beta01/examples/Simple /home/dementor/FOSS/ESBTL-1.0-beta01/examples/Simple /home/dementor/FOSS/ESBTL-1.0-beta01/examples/Simple/CMakeFiles/generic_select.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/generic_select.dir/depend
