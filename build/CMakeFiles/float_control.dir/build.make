# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.16

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


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
CMAKE_SOURCE_DIR = /home/xiao/proj/openloong-dyn-control

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/xiao/proj/openloong-dyn-control/build

# Include any dependencies generated for this target.
include CMakeFiles/float_control.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/float_control.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/float_control.dir/flags.make

CMakeFiles/float_control.dir/demo/float_control.cpp.o: CMakeFiles/float_control.dir/flags.make
CMakeFiles/float_control.dir/demo/float_control.cpp.o: ../demo/float_control.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/xiao/proj/openloong-dyn-control/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/float_control.dir/demo/float_control.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/float_control.dir/demo/float_control.cpp.o -c /home/xiao/proj/openloong-dyn-control/demo/float_control.cpp

CMakeFiles/float_control.dir/demo/float_control.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/float_control.dir/demo/float_control.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/xiao/proj/openloong-dyn-control/demo/float_control.cpp > CMakeFiles/float_control.dir/demo/float_control.cpp.i

CMakeFiles/float_control.dir/demo/float_control.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/float_control.dir/demo/float_control.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/xiao/proj/openloong-dyn-control/demo/float_control.cpp -o CMakeFiles/float_control.dir/demo/float_control.cpp.s

# Object files for target float_control
float_control_OBJECTS = \
"CMakeFiles/float_control.dir/demo/float_control.cpp.o"

# External object files for target float_control
float_control_EXTERNAL_OBJECTS =

float_control: CMakeFiles/float_control.dir/demo/float_control.cpp.o
float_control: CMakeFiles/float_control.dir/build.make
float_control: libcore.a
float_control: CMakeFiles/float_control.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/xiao/proj/openloong-dyn-control/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable float_control"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/float_control.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/float_control.dir/build: float_control

.PHONY : CMakeFiles/float_control.dir/build

CMakeFiles/float_control.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/float_control.dir/cmake_clean.cmake
.PHONY : CMakeFiles/float_control.dir/clean

CMakeFiles/float_control.dir/depend:
	cd /home/xiao/proj/openloong-dyn-control/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/xiao/proj/openloong-dyn-control /home/xiao/proj/openloong-dyn-control /home/xiao/proj/openloong-dyn-control/build /home/xiao/proj/openloong-dyn-control/build /home/xiao/proj/openloong-dyn-control/build/CMakeFiles/float_control.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/float_control.dir/depend

