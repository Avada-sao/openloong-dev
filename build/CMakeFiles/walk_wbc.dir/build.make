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
include CMakeFiles/walk_wbc.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/walk_wbc.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/walk_wbc.dir/flags.make

CMakeFiles/walk_wbc.dir/demo/walk_wbc.cpp.o: CMakeFiles/walk_wbc.dir/flags.make
CMakeFiles/walk_wbc.dir/demo/walk_wbc.cpp.o: ../demo/walk_wbc.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/xiao/proj/openloong-dyn-control/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/walk_wbc.dir/demo/walk_wbc.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/walk_wbc.dir/demo/walk_wbc.cpp.o -c /home/xiao/proj/openloong-dyn-control/demo/walk_wbc.cpp

CMakeFiles/walk_wbc.dir/demo/walk_wbc.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/walk_wbc.dir/demo/walk_wbc.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/xiao/proj/openloong-dyn-control/demo/walk_wbc.cpp > CMakeFiles/walk_wbc.dir/demo/walk_wbc.cpp.i

CMakeFiles/walk_wbc.dir/demo/walk_wbc.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/walk_wbc.dir/demo/walk_wbc.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/xiao/proj/openloong-dyn-control/demo/walk_wbc.cpp -o CMakeFiles/walk_wbc.dir/demo/walk_wbc.cpp.s

# Object files for target walk_wbc
walk_wbc_OBJECTS = \
"CMakeFiles/walk_wbc.dir/demo/walk_wbc.cpp.o"

# External object files for target walk_wbc
walk_wbc_EXTERNAL_OBJECTS =

walk_wbc: CMakeFiles/walk_wbc.dir/demo/walk_wbc.cpp.o
walk_wbc: CMakeFiles/walk_wbc.dir/build.make
walk_wbc: libcore.a
walk_wbc: CMakeFiles/walk_wbc.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/xiao/proj/openloong-dyn-control/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable walk_wbc"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/walk_wbc.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/walk_wbc.dir/build: walk_wbc

.PHONY : CMakeFiles/walk_wbc.dir/build

CMakeFiles/walk_wbc.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/walk_wbc.dir/cmake_clean.cmake
.PHONY : CMakeFiles/walk_wbc.dir/clean

CMakeFiles/walk_wbc.dir/depend:
	cd /home/xiao/proj/openloong-dyn-control/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/xiao/proj/openloong-dyn-control /home/xiao/proj/openloong-dyn-control /home/xiao/proj/openloong-dyn-control/build /home/xiao/proj/openloong-dyn-control/build /home/xiao/proj/openloong-dyn-control/build/CMakeFiles/walk_wbc.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/walk_wbc.dir/depend

