# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.18

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Disable VCS-based implicit rules.
% : %,v


# Disable VCS-based implicit rules.
% : RCS/%


# Disable VCS-based implicit rules.
% : RCS/%,v


# Disable VCS-based implicit rules.
% : SCCS/s.%


# Disable VCS-based implicit rules.
% : s.%


.SUFFIXES: .hpux_make_needs_suffix_list


# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/local/bin/cmake

# The command to remove a file.
RM = /usr/local/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/haitan/Codes/RAZOR

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/haitan/Codes/RAZOR/build

# Include any dependencies generated for this target.
include CMakeFiles/RAZOR-main.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/RAZOR-main.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/RAZOR-main.dir/flags.make

CMakeFiles/RAZOR-main.dir/src/main.cpp.o: CMakeFiles/RAZOR-main.dir/flags.make
CMakeFiles/RAZOR-main.dir/src/main.cpp.o: ../src/main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/haitan/Codes/RAZOR/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/RAZOR-main.dir/src/main.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/RAZOR-main.dir/src/main.cpp.o -c /home/haitan/Codes/RAZOR/src/main.cpp

CMakeFiles/RAZOR-main.dir/src/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/RAZOR-main.dir/src/main.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/haitan/Codes/RAZOR/src/main.cpp > CMakeFiles/RAZOR-main.dir/src/main.cpp.i

CMakeFiles/RAZOR-main.dir/src/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/RAZOR-main.dir/src/main.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/haitan/Codes/RAZOR/src/main.cpp -o CMakeFiles/RAZOR-main.dir/src/main.cpp.s

# Object files for target RAZOR-main
RAZOR__main_OBJECTS = \
"CMakeFiles/RAZOR-main.dir/src/main.cpp.o"

# External object files for target RAZOR-main
RAZOR__main_EXTERNAL_OBJECTS =

RAZOR-main: CMakeFiles/RAZOR-main.dir/src/main.cpp.o
RAZOR-main: CMakeFiles/RAZOR-main.dir/build.make
RAZOR-main: libRAZOR.so.0.0.1
RAZOR-main: /usr/local/lib/smart-math/libsmart-math.so
RAZOR-main: /usr/local/lib/smart-uq/libsmart-uq.so
RAZOR-main: /usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5_cpp.so
RAZOR-main: /usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5.so
RAZOR-main: /usr/lib/x86_64-linux-gnu/libpthread.so
RAZOR-main: /usr/lib/x86_64-linux-gnu/libsz.so
RAZOR-main: /usr/lib/x86_64-linux-gnu/libz.so
RAZOR-main: /usr/lib/x86_64-linux-gnu/libdl.so
RAZOR-main: /usr/lib/x86_64-linux-gnu/libm.so
RAZOR-main: /usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5_cpp.so
RAZOR-main: /usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5.so
RAZOR-main: /usr/lib/x86_64-linux-gnu/libpthread.so
RAZOR-main: /usr/lib/x86_64-linux-gnu/libsz.so
RAZOR-main: /usr/lib/x86_64-linux-gnu/libz.so
RAZOR-main: /usr/lib/x86_64-linux-gnu/libdl.so
RAZOR-main: /usr/lib/x86_64-linux-gnu/libm.so
RAZOR-main: CMakeFiles/RAZOR-main.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/haitan/Codes/RAZOR/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable RAZOR-main"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/RAZOR-main.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/RAZOR-main.dir/build: RAZOR-main

.PHONY : CMakeFiles/RAZOR-main.dir/build

CMakeFiles/RAZOR-main.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/RAZOR-main.dir/cmake_clean.cmake
.PHONY : CMakeFiles/RAZOR-main.dir/clean

CMakeFiles/RAZOR-main.dir/depend:
	cd /home/haitan/Codes/RAZOR/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/haitan/Codes/RAZOR /home/haitan/Codes/RAZOR /home/haitan/Codes/RAZOR/build /home/haitan/Codes/RAZOR/build /home/haitan/Codes/RAZOR/build/CMakeFiles/RAZOR-main.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/RAZOR-main.dir/depend
