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
include CMakeFiles/walk_mpc_wbc_joystick.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/walk_mpc_wbc_joystick.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/walk_mpc_wbc_joystick.dir/flags.make

CMakeFiles/walk_mpc_wbc_joystick.dir/demo/walk_mpc_wbc_joystick.cpp.o: CMakeFiles/walk_mpc_wbc_joystick.dir/flags.make
CMakeFiles/walk_mpc_wbc_joystick.dir/demo/walk_mpc_wbc_joystick.cpp.o: ../demo/walk_mpc_wbc_joystick.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/xiao/proj/openloong-dyn-control/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/walk_mpc_wbc_joystick.dir/demo/walk_mpc_wbc_joystick.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/walk_mpc_wbc_joystick.dir/demo/walk_mpc_wbc_joystick.cpp.o -c /home/xiao/proj/openloong-dyn-control/demo/walk_mpc_wbc_joystick.cpp

CMakeFiles/walk_mpc_wbc_joystick.dir/demo/walk_mpc_wbc_joystick.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/walk_mpc_wbc_joystick.dir/demo/walk_mpc_wbc_joystick.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/xiao/proj/openloong-dyn-control/demo/walk_mpc_wbc_joystick.cpp > CMakeFiles/walk_mpc_wbc_joystick.dir/demo/walk_mpc_wbc_joystick.cpp.i

CMakeFiles/walk_mpc_wbc_joystick.dir/demo/walk_mpc_wbc_joystick.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/walk_mpc_wbc_joystick.dir/demo/walk_mpc_wbc_joystick.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/xiao/proj/openloong-dyn-control/demo/walk_mpc_wbc_joystick.cpp -o CMakeFiles/walk_mpc_wbc_joystick.dir/demo/walk_mpc_wbc_joystick.cpp.s

# Object files for target walk_mpc_wbc_joystick
walk_mpc_wbc_joystick_OBJECTS = \
"CMakeFiles/walk_mpc_wbc_joystick.dir/demo/walk_mpc_wbc_joystick.cpp.o"

# External object files for target walk_mpc_wbc_joystick
walk_mpc_wbc_joystick_EXTERNAL_OBJECTS =

walk_mpc_wbc_joystick: CMakeFiles/walk_mpc_wbc_joystick.dir/demo/walk_mpc_wbc_joystick.cpp.o
walk_mpc_wbc_joystick: CMakeFiles/walk_mpc_wbc_joystick.dir/build.make
walk_mpc_wbc_joystick: libcore.a
walk_mpc_wbc_joystick: CMakeFiles/walk_mpc_wbc_joystick.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/xiao/proj/openloong-dyn-control/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable walk_mpc_wbc_joystick"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/walk_mpc_wbc_joystick.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/walk_mpc_wbc_joystick.dir/build: walk_mpc_wbc_joystick

.PHONY : CMakeFiles/walk_mpc_wbc_joystick.dir/build

CMakeFiles/walk_mpc_wbc_joystick.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/walk_mpc_wbc_joystick.dir/cmake_clean.cmake
.PHONY : CMakeFiles/walk_mpc_wbc_joystick.dir/clean

CMakeFiles/walk_mpc_wbc_joystick.dir/depend:
	cd /home/xiao/proj/openloong-dyn-control/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/xiao/proj/openloong-dyn-control /home/xiao/proj/openloong-dyn-control /home/xiao/proj/openloong-dyn-control/build /home/xiao/proj/openloong-dyn-control/build /home/xiao/proj/openloong-dyn-control/build/CMakeFiles/walk_mpc_wbc_joystick.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/walk_mpc_wbc_joystick.dir/depend

