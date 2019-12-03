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
include ext/gfa1/CMakeFiles/gfa1.dir/depend.make

# Include the progress variables for this target.
include ext/gfa1/CMakeFiles/gfa1.dir/progress.make

# Include the compile flags for this target's objects.
include ext/gfa1/CMakeFiles/gfa1.dir/flags.make

ext/gfa1/CMakeFiles/gfa1.dir/gfa.c.o: ext/gfa1/CMakeFiles/gfa1.dir/flags.make
ext/gfa1/CMakeFiles/gfa1.dir/gfa.c.o: /home/chakenal/ontx-pipeline/SPAdes-3.13.1/ext/src/gfa1/gfa.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/chakenal/ontx-pipeline/SPAdes-3.13.1/build_spades/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object ext/gfa1/CMakeFiles/gfa1.dir/gfa.c.o"
	cd /home/chakenal/ontx-pipeline/SPAdes-3.13.1/build_spades/ext/gfa1 && /usr/bin/cc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/gfa1.dir/gfa.c.o   -c /home/chakenal/ontx-pipeline/SPAdes-3.13.1/ext/src/gfa1/gfa.c

ext/gfa1/CMakeFiles/gfa1.dir/gfa.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/gfa1.dir/gfa.c.i"
	cd /home/chakenal/ontx-pipeline/SPAdes-3.13.1/build_spades/ext/gfa1 && /usr/bin/cc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/chakenal/ontx-pipeline/SPAdes-3.13.1/ext/src/gfa1/gfa.c > CMakeFiles/gfa1.dir/gfa.c.i

ext/gfa1/CMakeFiles/gfa1.dir/gfa.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/gfa1.dir/gfa.c.s"
	cd /home/chakenal/ontx-pipeline/SPAdes-3.13.1/build_spades/ext/gfa1 && /usr/bin/cc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/chakenal/ontx-pipeline/SPAdes-3.13.1/ext/src/gfa1/gfa.c -o CMakeFiles/gfa1.dir/gfa.c.s

ext/gfa1/CMakeFiles/gfa1.dir/gfa.c.o.requires:

.PHONY : ext/gfa1/CMakeFiles/gfa1.dir/gfa.c.o.requires

ext/gfa1/CMakeFiles/gfa1.dir/gfa.c.o.provides: ext/gfa1/CMakeFiles/gfa1.dir/gfa.c.o.requires
	$(MAKE) -f ext/gfa1/CMakeFiles/gfa1.dir/build.make ext/gfa1/CMakeFiles/gfa1.dir/gfa.c.o.provides.build
.PHONY : ext/gfa1/CMakeFiles/gfa1.dir/gfa.c.o.provides

ext/gfa1/CMakeFiles/gfa1.dir/gfa.c.o.provides.build: ext/gfa1/CMakeFiles/gfa1.dir/gfa.c.o


# Object files for target gfa1
gfa1_OBJECTS = \
"CMakeFiles/gfa1.dir/gfa.c.o"

# External object files for target gfa1
gfa1_EXTERNAL_OBJECTS =

ext/gfa1/libgfa1.a: ext/gfa1/CMakeFiles/gfa1.dir/gfa.c.o
ext/gfa1/libgfa1.a: ext/gfa1/CMakeFiles/gfa1.dir/build.make
ext/gfa1/libgfa1.a: ext/gfa1/CMakeFiles/gfa1.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/chakenal/ontx-pipeline/SPAdes-3.13.1/build_spades/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking C static library libgfa1.a"
	cd /home/chakenal/ontx-pipeline/SPAdes-3.13.1/build_spades/ext/gfa1 && $(CMAKE_COMMAND) -P CMakeFiles/gfa1.dir/cmake_clean_target.cmake
	cd /home/chakenal/ontx-pipeline/SPAdes-3.13.1/build_spades/ext/gfa1 && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/gfa1.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
ext/gfa1/CMakeFiles/gfa1.dir/build: ext/gfa1/libgfa1.a

.PHONY : ext/gfa1/CMakeFiles/gfa1.dir/build

ext/gfa1/CMakeFiles/gfa1.dir/requires: ext/gfa1/CMakeFiles/gfa1.dir/gfa.c.o.requires

.PHONY : ext/gfa1/CMakeFiles/gfa1.dir/requires

ext/gfa1/CMakeFiles/gfa1.dir/clean:
	cd /home/chakenal/ontx-pipeline/SPAdes-3.13.1/build_spades/ext/gfa1 && $(CMAKE_COMMAND) -P CMakeFiles/gfa1.dir/cmake_clean.cmake
.PHONY : ext/gfa1/CMakeFiles/gfa1.dir/clean

ext/gfa1/CMakeFiles/gfa1.dir/depend:
	cd /home/chakenal/ontx-pipeline/SPAdes-3.13.1/build_spades && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/chakenal/ontx-pipeline/SPAdes-3.13.1/src /home/chakenal/ontx-pipeline/SPAdes-3.13.1/ext/src/gfa1 /home/chakenal/ontx-pipeline/SPAdes-3.13.1/build_spades /home/chakenal/ontx-pipeline/SPAdes-3.13.1/build_spades/ext/gfa1 /home/chakenal/ontx-pipeline/SPAdes-3.13.1/build_spades/ext/gfa1/CMakeFiles/gfa1.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : ext/gfa1/CMakeFiles/gfa1.dir/depend

