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
include CMakeFiles/RAZOR.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/RAZOR.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/RAZOR.dir/flags.make

CMakeFiles/RAZOR.dir/libraries/acquire_training_set.cpp.o: CMakeFiles/RAZOR.dir/flags.make
CMakeFiles/RAZOR.dir/libraries/acquire_training_set.cpp.o: ../libraries/acquire_training_set.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/haitan/Codes/RAZOR/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/RAZOR.dir/libraries/acquire_training_set.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/RAZOR.dir/libraries/acquire_training_set.cpp.o -c /home/haitan/Codes/RAZOR/libraries/acquire_training_set.cpp

CMakeFiles/RAZOR.dir/libraries/acquire_training_set.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/RAZOR.dir/libraries/acquire_training_set.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/haitan/Codes/RAZOR/libraries/acquire_training_set.cpp > CMakeFiles/RAZOR.dir/libraries/acquire_training_set.cpp.i

CMakeFiles/RAZOR.dir/libraries/acquire_training_set.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/RAZOR.dir/libraries/acquire_training_set.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/haitan/Codes/RAZOR/libraries/acquire_training_set.cpp -o CMakeFiles/RAZOR.dir/libraries/acquire_training_set.cpp.s

CMakeFiles/RAZOR.dir/libraries/low_dimensional_model.cpp.o: CMakeFiles/RAZOR.dir/flags.make
CMakeFiles/RAZOR.dir/libraries/low_dimensional_model.cpp.o: ../libraries/low_dimensional_model.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/haitan/Codes/RAZOR/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/RAZOR.dir/libraries/low_dimensional_model.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/RAZOR.dir/libraries/low_dimensional_model.cpp.o -c /home/haitan/Codes/RAZOR/libraries/low_dimensional_model.cpp

CMakeFiles/RAZOR.dir/libraries/low_dimensional_model.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/RAZOR.dir/libraries/low_dimensional_model.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/haitan/Codes/RAZOR/libraries/low_dimensional_model.cpp > CMakeFiles/RAZOR.dir/libraries/low_dimensional_model.cpp.i

CMakeFiles/RAZOR.dir/libraries/low_dimensional_model.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/RAZOR.dir/libraries/low_dimensional_model.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/haitan/Codes/RAZOR/libraries/low_dimensional_model.cpp -o CMakeFiles/RAZOR.dir/libraries/low_dimensional_model.cpp.s

CMakeFiles/RAZOR.dir/libraries/low_dimensional_solution.cpp.o: CMakeFiles/RAZOR.dir/flags.make
CMakeFiles/RAZOR.dir/libraries/low_dimensional_solution.cpp.o: ../libraries/low_dimensional_solution.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/haitan/Codes/RAZOR/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/RAZOR.dir/libraries/low_dimensional_solution.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/RAZOR.dir/libraries/low_dimensional_solution.cpp.o -c /home/haitan/Codes/RAZOR/libraries/low_dimensional_solution.cpp

CMakeFiles/RAZOR.dir/libraries/low_dimensional_solution.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/RAZOR.dir/libraries/low_dimensional_solution.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/haitan/Codes/RAZOR/libraries/low_dimensional_solution.cpp > CMakeFiles/RAZOR.dir/libraries/low_dimensional_solution.cpp.i

CMakeFiles/RAZOR.dir/libraries/low_dimensional_solution.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/RAZOR.dir/libraries/low_dimensional_solution.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/haitan/Codes/RAZOR/libraries/low_dimensional_solution.cpp -o CMakeFiles/RAZOR.dir/libraries/low_dimensional_solution.cpp.s

CMakeFiles/RAZOR.dir/libraries/modal_identification.cpp.o: CMakeFiles/RAZOR.dir/flags.make
CMakeFiles/RAZOR.dir/libraries/modal_identification.cpp.o: ../libraries/modal_identification.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/haitan/Codes/RAZOR/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object CMakeFiles/RAZOR.dir/libraries/modal_identification.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/RAZOR.dir/libraries/modal_identification.cpp.o -c /home/haitan/Codes/RAZOR/libraries/modal_identification.cpp

CMakeFiles/RAZOR.dir/libraries/modal_identification.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/RAZOR.dir/libraries/modal_identification.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/haitan/Codes/RAZOR/libraries/modal_identification.cpp > CMakeFiles/RAZOR.dir/libraries/modal_identification.cpp.i

CMakeFiles/RAZOR.dir/libraries/modal_identification.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/RAZOR.dir/libraries/modal_identification.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/haitan/Codes/RAZOR/libraries/modal_identification.cpp -o CMakeFiles/RAZOR.dir/libraries/modal_identification.cpp.s

CMakeFiles/RAZOR.dir/libraries/read_data.cpp.o: CMakeFiles/RAZOR.dir/flags.make
CMakeFiles/RAZOR.dir/libraries/read_data.cpp.o: ../libraries/read_data.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/haitan/Codes/RAZOR/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object CMakeFiles/RAZOR.dir/libraries/read_data.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/RAZOR.dir/libraries/read_data.cpp.o -c /home/haitan/Codes/RAZOR/libraries/read_data.cpp

CMakeFiles/RAZOR.dir/libraries/read_data.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/RAZOR.dir/libraries/read_data.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/haitan/Codes/RAZOR/libraries/read_data.cpp > CMakeFiles/RAZOR.dir/libraries/read_data.cpp.i

CMakeFiles/RAZOR.dir/libraries/read_data.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/RAZOR.dir/libraries/read_data.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/haitan/Codes/RAZOR/libraries/read_data.cpp -o CMakeFiles/RAZOR.dir/libraries/read_data.cpp.s

CMakeFiles/RAZOR.dir/libraries/write_data.cpp.o: CMakeFiles/RAZOR.dir/flags.make
CMakeFiles/RAZOR.dir/libraries/write_data.cpp.o: ../libraries/write_data.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/haitan/Codes/RAZOR/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object CMakeFiles/RAZOR.dir/libraries/write_data.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/RAZOR.dir/libraries/write_data.cpp.o -c /home/haitan/Codes/RAZOR/libraries/write_data.cpp

CMakeFiles/RAZOR.dir/libraries/write_data.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/RAZOR.dir/libraries/write_data.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/haitan/Codes/RAZOR/libraries/write_data.cpp > CMakeFiles/RAZOR.dir/libraries/write_data.cpp.i

CMakeFiles/RAZOR.dir/libraries/write_data.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/RAZOR.dir/libraries/write_data.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/haitan/Codes/RAZOR/libraries/write_data.cpp -o CMakeFiles/RAZOR.dir/libraries/write_data.cpp.s

# Object files for target RAZOR
RAZOR_OBJECTS = \
"CMakeFiles/RAZOR.dir/libraries/acquire_training_set.cpp.o" \
"CMakeFiles/RAZOR.dir/libraries/low_dimensional_model.cpp.o" \
"CMakeFiles/RAZOR.dir/libraries/low_dimensional_solution.cpp.o" \
"CMakeFiles/RAZOR.dir/libraries/modal_identification.cpp.o" \
"CMakeFiles/RAZOR.dir/libraries/read_data.cpp.o" \
"CMakeFiles/RAZOR.dir/libraries/write_data.cpp.o"

# External object files for target RAZOR
RAZOR_EXTERNAL_OBJECTS =

libRAZOR.so.0.0.1: CMakeFiles/RAZOR.dir/libraries/acquire_training_set.cpp.o
libRAZOR.so.0.0.1: CMakeFiles/RAZOR.dir/libraries/low_dimensional_model.cpp.o
libRAZOR.so.0.0.1: CMakeFiles/RAZOR.dir/libraries/low_dimensional_solution.cpp.o
libRAZOR.so.0.0.1: CMakeFiles/RAZOR.dir/libraries/modal_identification.cpp.o
libRAZOR.so.0.0.1: CMakeFiles/RAZOR.dir/libraries/read_data.cpp.o
libRAZOR.so.0.0.1: CMakeFiles/RAZOR.dir/libraries/write_data.cpp.o
libRAZOR.so.0.0.1: CMakeFiles/RAZOR.dir/build.make
libRAZOR.so.0.0.1: /usr/local/lib/smart-math/libsmart-math.so
libRAZOR.so.0.0.1: /usr/local/lib/smart-uq/libsmart-uq.so
libRAZOR.so.0.0.1: /usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5_cpp.so
libRAZOR.so.0.0.1: /usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5.so
libRAZOR.so.0.0.1: /usr/lib/x86_64-linux-gnu/libpthread.so
libRAZOR.so.0.0.1: /usr/lib/x86_64-linux-gnu/libsz.so
libRAZOR.so.0.0.1: /usr/lib/x86_64-linux-gnu/libz.so
libRAZOR.so.0.0.1: /usr/lib/x86_64-linux-gnu/libdl.so
libRAZOR.so.0.0.1: /usr/lib/x86_64-linux-gnu/libm.so
libRAZOR.so.0.0.1: /usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5_cpp.so
libRAZOR.so.0.0.1: /usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5.so
libRAZOR.so.0.0.1: /usr/lib/x86_64-linux-gnu/libpthread.so
libRAZOR.so.0.0.1: /usr/lib/x86_64-linux-gnu/libsz.so
libRAZOR.so.0.0.1: /usr/lib/x86_64-linux-gnu/libz.so
libRAZOR.so.0.0.1: /usr/lib/x86_64-linux-gnu/libdl.so
libRAZOR.so.0.0.1: /usr/lib/x86_64-linux-gnu/libm.so
libRAZOR.so.0.0.1: CMakeFiles/RAZOR.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/haitan/Codes/RAZOR/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Linking CXX shared library libRAZOR.so"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/RAZOR.dir/link.txt --verbose=$(VERBOSE)
	$(CMAKE_COMMAND) -E cmake_symlink_library libRAZOR.so.0.0.1 libRAZOR.so.0.0.1 libRAZOR.so

libRAZOR.so: libRAZOR.so.0.0.1
	@$(CMAKE_COMMAND) -E touch_nocreate libRAZOR.so

# Rule to build all files generated by this target.
CMakeFiles/RAZOR.dir/build: libRAZOR.so

.PHONY : CMakeFiles/RAZOR.dir/build

CMakeFiles/RAZOR.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/RAZOR.dir/cmake_clean.cmake
.PHONY : CMakeFiles/RAZOR.dir/clean

CMakeFiles/RAZOR.dir/depend:
	cd /home/haitan/Codes/RAZOR/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/haitan/Codes/RAZOR /home/haitan/Codes/RAZOR /home/haitan/Codes/RAZOR/build /home/haitan/Codes/RAZOR/build /home/haitan/Codes/RAZOR/build/CMakeFiles/RAZOR.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/RAZOR.dir/depend

