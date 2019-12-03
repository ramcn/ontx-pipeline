# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.5

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
CMAKE_SOURCE_DIR = /home/chakenal/ontx-pipeline/SPAdes-3.13.1/src

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/chakenal/ontx-pipeline/SPAdes-3.13.1/build_spades

# Include any dependencies generated for this target.
include common/sequence/CMakeFiles/sequence.dir/depend.make

# Include the progress variables for this target.
include common/sequence/CMakeFiles/sequence.dir/progress.make

# Include the compile flags for this target's objects.
include common/sequence/CMakeFiles/sequence.dir/flags.make

common/sequence/CMakeFiles/sequence.dir/sequence_tools.cpp.o: common/sequence/CMakeFiles/sequence.dir/flags.make
common/sequence/CMakeFiles/sequence.dir/sequence_tools.cpp.o: /home/chakenal/ontx-pipeline/SPAdes-3.13.1/src/common/sequence/sequence_tools.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/chakenal/ontx-pipeline/SPAdes-3.13.1/build_spades/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object common/sequence/CMakeFiles/sequence.dir/sequence_tools.cpp.o"
	cd /home/chakenal/ontx-pipeline/SPAdes-3.13.1/build_spades/common/sequence && /usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/sequence.dir/sequence_tools.cpp.o -c /home/chakenal/ontx-pipeline/SPAdes-3.13.1/src/common/sequence/sequence_tools.cpp

common/sequence/CMakeFiles/sequence.dir/sequence_tools.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/sequence.dir/sequence_tools.cpp.i"
	cd /home/chakenal/ontx-pipeline/SPAdes-3.13.1/build_spades/common/sequence && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/chakenal/ontx-pipeline/SPAdes-3.13.1/src/common/sequence/sequence_tools.cpp > CMakeFiles/sequence.dir/sequence_tools.cpp.i

common/sequence/CMakeFiles/sequence.dir/sequence_tools.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/sequence.dir/sequence_tools.cpp.s"
	cd /home/chakenal/ontx-pipeline/SPAdes-3.13.1/build_spades/common/sequence && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/chakenal/ontx-pipeline/SPAdes-3.13.1/src/common/sequence/sequence_tools.cpp -o CMakeFiles/sequence.dir/sequence_tools.cpp.s

common/sequence/CMakeFiles/sequence.dir/sequence_tools.cpp.o.requires:

.PHONY : common/sequence/CMakeFiles/sequence.dir/sequence_tools.cpp.o.requires

common/sequence/CMakeFiles/sequence.dir/sequence_tools.cpp.o.provides: common/sequence/CMakeFiles/sequence.dir/sequence_tools.cpp.o.requires
	$(MAKE) -f common/sequence/CMakeFiles/sequence.dir/build.make common/sequence/CMakeFiles/sequence.dir/sequence_tools.cpp.o.provides.build
.PHONY : common/sequence/CMakeFiles/sequence.dir/sequence_tools.cpp.o.provides

common/sequence/CMakeFiles/sequence.dir/sequence_tools.cpp.o.provides.build: common/sequence/CMakeFiles/sequence.dir/sequence_tools.cpp.o


# Object files for target sequence
sequence_OBJECTS = \
"CMakeFiles/sequence.dir/sequence_tools.cpp.o"

# External object files for target sequence
sequence_EXTERNAL_OBJECTS =

common/sequence/libsequence.a: common/sequence/CMakeFiles/sequence.dir/sequence_tools.cpp.o
common/sequence/libsequence.a: common/sequence/CMakeFiles/sequence.dir/build.make
common/sequence/libsequence.a: common/sequence/CMakeFiles/sequence.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/chakenal/ontx-pipeline/SPAdes-3.13.1/build_spades/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX static library libsequence.a"
	cd /home/chakenal/ontx-pipeline/SPAdes-3.13.1/build_spades/common/sequence && $(CMAKE_COMMAND) -P CMakeFiles/sequence.dir/cmake_clean_target.cmake
	cd /home/chakenal/ontx-pipeline/SPAdes-3.13.1/build_spades/common/sequence && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/sequence.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
common/sequence/CMakeFiles/sequence.dir/build: common/sequence/libsequence.a

.PHONY : common/sequence/CMakeFiles/sequence.dir/build

common/sequence/CMakeFiles/sequence.dir/requires: common/sequence/CMakeFiles/sequence.dir/sequence_tools.cpp.o.requires

.PHONY : common/sequence/CMakeFiles/sequence.dir/requires

common/sequence/CMakeFiles/sequence.dir/clean:
	cd /home/chakenal/ontx-pipeline/SPAdes-3.13.1/build_spades/common/sequence && $(CMAKE_COMMAND) -P CMakeFiles/sequence.dir/cmake_clean.cmake
.PHONY : common/sequence/CMakeFiles/sequence.dir/clean

common/sequence/CMakeFiles/sequence.dir/depend:
	cd /home/chakenal/ontx-pipeline/SPAdes-3.13.1/build_spades && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/chakenal/ontx-pipeline/SPAdes-3.13.1/src /home/chakenal/ontx-pipeline/SPAdes-3.13.1/src/common/sequence /home/chakenal/ontx-pipeline/SPAdes-3.13.1/build_spades /home/chakenal/ontx-pipeline/SPAdes-3.13.1/build_spades/common/sequence /home/chakenal/ontx-pipeline/SPAdes-3.13.1/build_spades/common/sequence/CMakeFiles/sequence.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : common/sequence/CMakeFiles/sequence.dir/depend

