# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.17

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
CMAKE_COMMAND = /home/garen_lee/jialiang.li/software/clion-2020.3.4/bin/cmake/linux/bin/cmake

# The command to remove a file.
RM = /home/garen_lee/jialiang.li/software/clion-2020.3.4/bin/cmake/linux/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/garen_lee/jialiang.li/code/test/ctest_pro

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/garen_lee/jialiang.li/code/test/ctest_pro/cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/ctest_pro.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/ctest_pro.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/ctest_pro.dir/flags.make

CMakeFiles/ctest_pro.dir/main.cpp.o: CMakeFiles/ctest_pro.dir/flags.make
CMakeFiles/ctest_pro.dir/main.cpp.o: ../main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/garen_lee/jialiang.li/code/test/ctest_pro/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/ctest_pro.dir/main.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/ctest_pro.dir/main.cpp.o -c /home/garen_lee/jialiang.li/code/test/ctest_pro/main.cpp

CMakeFiles/ctest_pro.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/ctest_pro.dir/main.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/garen_lee/jialiang.li/code/test/ctest_pro/main.cpp > CMakeFiles/ctest_pro.dir/main.cpp.i

CMakeFiles/ctest_pro.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/ctest_pro.dir/main.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/garen_lee/jialiang.li/code/test/ctest_pro/main.cpp -o CMakeFiles/ctest_pro.dir/main.cpp.s

# Object files for target ctest_pro
ctest_pro_OBJECTS = \
"CMakeFiles/ctest_pro.dir/main.cpp.o"

# External object files for target ctest_pro
ctest_pro_EXTERNAL_OBJECTS =

ctest_pro: CMakeFiles/ctest_pro.dir/main.cpp.o
ctest_pro: CMakeFiles/ctest_pro.dir/build.make
ctest_pro: CMakeFiles/ctest_pro.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/garen_lee/jialiang.li/code/test/ctest_pro/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ctest_pro"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/ctest_pro.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/ctest_pro.dir/build: ctest_pro

.PHONY : CMakeFiles/ctest_pro.dir/build

CMakeFiles/ctest_pro.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/ctest_pro.dir/cmake_clean.cmake
.PHONY : CMakeFiles/ctest_pro.dir/clean

CMakeFiles/ctest_pro.dir/depend:
	cd /home/garen_lee/jialiang.li/code/test/ctest_pro/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/garen_lee/jialiang.li/code/test/ctest_pro /home/garen_lee/jialiang.li/code/test/ctest_pro /home/garen_lee/jialiang.li/code/test/ctest_pro/cmake-build-debug /home/garen_lee/jialiang.li/code/test/ctest_pro/cmake-build-debug /home/garen_lee/jialiang.li/code/test/ctest_pro/cmake-build-debug/CMakeFiles/ctest_pro.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/ctest_pro.dir/depend

