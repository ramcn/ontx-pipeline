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
include common/utils/CMakeFiles/utils.dir/depend.make

# Include the progress variables for this target.
include common/utils/CMakeFiles/utils.dir/progress.make

# Include the compile flags for this target's objects.
include common/utils/CMakeFiles/utils.dir/flags.make

common/utils/CMakeFiles/utils.dir/memory_limit.cpp.o: common/utils/CMakeFiles/utils.dir/flags.make
common/utils/CMakeFiles/utils.dir/memory_limit.cpp.o: /home/chakenal/ontx-pipeline/SPAdes-3.13.1/src/common/utils/memory_limit.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/chakenal/ontx-pipeline/SPAdes-3.13.1/build_spades/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object common/utils/CMakeFiles/utils.dir/memory_limit.cpp.o"
	cd /home/chakenal/ontx-pipeline/SPAdes-3.13.1/build_spades/common/utils && /usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/utils.dir/memory_limit.cpp.o -c /home/chakenal/ontx-pipeline/SPAdes-3.13.1/src/common/utils/memory_limit.cpp

common/utils/CMakeFiles/utils.dir/memory_limit.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/utils.dir/memory_limit.cpp.i"
	cd /home/chakenal/ontx-pipeline/SPAdes-3.13.1/build_spades/common/utils && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/chakenal/ontx-pipeline/SPAdes-3.13.1/src/common/utils/memory_limit.cpp > CMakeFiles/utils.dir/memory_limit.cpp.i

common/utils/CMakeFiles/utils.dir/memory_limit.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/utils.dir/memory_limit.cpp.s"
	cd /home/chakenal/ontx-pipeline/SPAdes-3.13.1/build_spades/common/utils && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/chakenal/ontx-pipeline/SPAdes-3.13.1/src/common/utils/memory_limit.cpp -o CMakeFiles/utils.dir/memory_limit.cpp.s

common/utils/CMakeFiles/utils.dir/memory_limit.cpp.o.requires:

.PHONY : common/utils/CMakeFiles/utils.dir/memory_limit.cpp.o.requires

common/utils/CMakeFiles/utils.dir/memory_limit.cpp.o.provides: common/utils/CMakeFiles/utils.dir/memory_limit.cpp.o.requires
	$(MAKE) -f common/utils/CMakeFiles/utils.dir/build.make common/utils/CMakeFiles/utils.dir/memory_limit.cpp.o.provides.build
.PHONY : common/utils/CMakeFiles/utils.dir/memory_limit.cpp.o.provides

common/utils/CMakeFiles/utils.dir/memory_limit.cpp.o.provides.build: common/utils/CMakeFiles/utils.dir/memory_limit.cpp.o


common/utils/CMakeFiles/utils.dir/filesystem/copy_file.cpp.o: common/utils/CMakeFiles/utils.dir/flags.make
common/utils/CMakeFiles/utils.dir/filesystem/copy_file.cpp.o: /home/chakenal/ontx-pipeline/SPAdes-3.13.1/src/common/utils/filesystem/copy_file.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/chakenal/ontx-pipeline/SPAdes-3.13.1/build_spades/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object common/utils/CMakeFiles/utils.dir/filesystem/copy_file.cpp.o"
	cd /home/chakenal/ontx-pipeline/SPAdes-3.13.1/build_spades/common/utils && /usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/utils.dir/filesystem/copy_file.cpp.o -c /home/chakenal/ontx-pipeline/SPAdes-3.13.1/src/common/utils/filesystem/copy_file.cpp

common/utils/CMakeFiles/utils.dir/filesystem/copy_file.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/utils.dir/filesystem/copy_file.cpp.i"
	cd /home/chakenal/ontx-pipeline/SPAdes-3.13.1/build_spades/common/utils && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/chakenal/ontx-pipeline/SPAdes-3.13.1/src/common/utils/filesystem/copy_file.cpp > CMakeFiles/utils.dir/filesystem/copy_file.cpp.i

common/utils/CMakeFiles/utils.dir/filesystem/copy_file.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/utils.dir/filesystem/copy_file.cpp.s"
	cd /home/chakenal/ontx-pipeline/SPAdes-3.13.1/build_spades/common/utils && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/chakenal/ontx-pipeline/SPAdes-3.13.1/src/common/utils/filesystem/copy_file.cpp -o CMakeFiles/utils.dir/filesystem/copy_file.cpp.s

common/utils/CMakeFiles/utils.dir/filesystem/copy_file.cpp.o.requires:

.PHONY : common/utils/CMakeFiles/utils.dir/filesystem/copy_file.cpp.o.requires

common/utils/CMakeFiles/utils.dir/filesystem/copy_file.cpp.o.provides: common/utils/CMakeFiles/utils.dir/filesystem/copy_file.cpp.o.requires
	$(MAKE) -f common/utils/CMakeFiles/utils.dir/build.make common/utils/CMakeFiles/utils.dir/filesystem/copy_file.cpp.o.provides.build
.PHONY : common/utils/CMakeFiles/utils.dir/filesystem/copy_file.cpp.o.provides

common/utils/CMakeFiles/utils.dir/filesystem/copy_file.cpp.o.provides.build: common/utils/CMakeFiles/utils.dir/filesystem/copy_file.cpp.o


common/utils/CMakeFiles/utils.dir/filesystem/path_helper.cpp.o: common/utils/CMakeFiles/utils.dir/flags.make
common/utils/CMakeFiles/utils.dir/filesystem/path_helper.cpp.o: /home/chakenal/ontx-pipeline/SPAdes-3.13.1/src/common/utils/filesystem/path_helper.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/chakenal/ontx-pipeline/SPAdes-3.13.1/build_spades/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object common/utils/CMakeFiles/utils.dir/filesystem/path_helper.cpp.o"
	cd /home/chakenal/ontx-pipeline/SPAdes-3.13.1/build_spades/common/utils && /usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/utils.dir/filesystem/path_helper.cpp.o -c /home/chakenal/ontx-pipeline/SPAdes-3.13.1/src/common/utils/filesystem/path_helper.cpp

common/utils/CMakeFiles/utils.dir/filesystem/path_helper.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/utils.dir/filesystem/path_helper.cpp.i"
	cd /home/chakenal/ontx-pipeline/SPAdes-3.13.1/build_spades/common/utils && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/chakenal/ontx-pipeline/SPAdes-3.13.1/src/common/utils/filesystem/path_helper.cpp > CMakeFiles/utils.dir/filesystem/path_helper.cpp.i

common/utils/CMakeFiles/utils.dir/filesystem/path_helper.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/utils.dir/filesystem/path_helper.cpp.s"
	cd /home/chakenal/ontx-pipeline/SPAdes-3.13.1/build_spades/common/utils && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/chakenal/ontx-pipeline/SPAdes-3.13.1/src/common/utils/filesystem/path_helper.cpp -o CMakeFiles/utils.dir/filesystem/path_helper.cpp.s

common/utils/CMakeFiles/utils.dir/filesystem/path_helper.cpp.o.requires:

.PHONY : common/utils/CMakeFiles/utils.dir/filesystem/path_helper.cpp.o.requires

common/utils/CMakeFiles/utils.dir/filesystem/path_helper.cpp.o.provides: common/utils/CMakeFiles/utils.dir/filesystem/path_helper.cpp.o.requires
	$(MAKE) -f common/utils/CMakeFiles/utils.dir/build.make common/utils/CMakeFiles/utils.dir/filesystem/path_helper.cpp.o.provides.build
.PHONY : common/utils/CMakeFiles/utils.dir/filesystem/path_helper.cpp.o.provides

common/utils/CMakeFiles/utils.dir/filesystem/path_helper.cpp.o.provides.build: common/utils/CMakeFiles/utils.dir/filesystem/path_helper.cpp.o


common/utils/CMakeFiles/utils.dir/filesystem/temporary.cpp.o: common/utils/CMakeFiles/utils.dir/flags.make
common/utils/CMakeFiles/utils.dir/filesystem/temporary.cpp.o: /home/chakenal/ontx-pipeline/SPAdes-3.13.1/src/common/utils/filesystem/temporary.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/chakenal/ontx-pipeline/SPAdes-3.13.1/build_spades/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object common/utils/CMakeFiles/utils.dir/filesystem/temporary.cpp.o"
	cd /home/chakenal/ontx-pipeline/SPAdes-3.13.1/build_spades/common/utils && /usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/utils.dir/filesystem/temporary.cpp.o -c /home/chakenal/ontx-pipeline/SPAdes-3.13.1/src/common/utils/filesystem/temporary.cpp

common/utils/CMakeFiles/utils.dir/filesystem/temporary.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/utils.dir/filesystem/temporary.cpp.i"
	cd /home/chakenal/ontx-pipeline/SPAdes-3.13.1/build_spades/common/utils && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/chakenal/ontx-pipeline/SPAdes-3.13.1/src/common/utils/filesystem/temporary.cpp > CMakeFiles/utils.dir/filesystem/temporary.cpp.i

common/utils/CMakeFiles/utils.dir/filesystem/temporary.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/utils.dir/filesystem/temporary.cpp.s"
	cd /home/chakenal/ontx-pipeline/SPAdes-3.13.1/build_spades/common/utils && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/chakenal/ontx-pipeline/SPAdes-3.13.1/src/common/utils/filesystem/temporary.cpp -o CMakeFiles/utils.dir/filesystem/temporary.cpp.s

common/utils/CMakeFiles/utils.dir/filesystem/temporary.cpp.o.requires:

.PHONY : common/utils/CMakeFiles/utils.dir/filesystem/temporary.cpp.o.requires

common/utils/CMakeFiles/utils.dir/filesystem/temporary.cpp.o.provides: common/utils/CMakeFiles/utils.dir/filesystem/temporary.cpp.o.requires
	$(MAKE) -f common/utils/CMakeFiles/utils.dir/build.make common/utils/CMakeFiles/utils.dir/filesystem/temporary.cpp.o.provides.build
.PHONY : common/utils/CMakeFiles/utils.dir/filesystem/temporary.cpp.o.provides

common/utils/CMakeFiles/utils.dir/filesystem/temporary.cpp.o.provides.build: common/utils/CMakeFiles/utils.dir/filesystem/temporary.cpp.o


common/utils/CMakeFiles/utils.dir/logger/logger_impl.cpp.o: common/utils/CMakeFiles/utils.dir/flags.make
common/utils/CMakeFiles/utils.dir/logger/logger_impl.cpp.o: /home/chakenal/ontx-pipeline/SPAdes-3.13.1/src/common/utils/logger/logger_impl.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/chakenal/ontx-pipeline/SPAdes-3.13.1/build_spades/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object common/utils/CMakeFiles/utils.dir/logger/logger_impl.cpp.o"
	cd /home/chakenal/ontx-pipeline/SPAdes-3.13.1/build_spades/common/utils && /usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/utils.dir/logger/logger_impl.cpp.o -c /home/chakenal/ontx-pipeline/SPAdes-3.13.1/src/common/utils/logger/logger_impl.cpp

common/utils/CMakeFiles/utils.dir/logger/logger_impl.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/utils.dir/logger/logger_impl.cpp.i"
	cd /home/chakenal/ontx-pipeline/SPAdes-3.13.1/build_spades/common/utils && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/chakenal/ontx-pipeline/SPAdes-3.13.1/src/common/utils/logger/logger_impl.cpp > CMakeFiles/utils.dir/logger/logger_impl.cpp.i

common/utils/CMakeFiles/utils.dir/logger/logger_impl.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/utils.dir/logger/logger_impl.cpp.s"
	cd /home/chakenal/ontx-pipeline/SPAdes-3.13.1/build_spades/common/utils && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/chakenal/ontx-pipeline/SPAdes-3.13.1/src/common/utils/logger/logger_impl.cpp -o CMakeFiles/utils.dir/logger/logger_impl.cpp.s

common/utils/CMakeFiles/utils.dir/logger/logger_impl.cpp.o.requires:

.PHONY : common/utils/CMakeFiles/utils.dir/logger/logger_impl.cpp.o.requires

common/utils/CMakeFiles/utils.dir/logger/logger_impl.cpp.o.provides: common/utils/CMakeFiles/utils.dir/logger/logger_impl.cpp.o.requires
	$(MAKE) -f common/utils/CMakeFiles/utils.dir/build.make common/utils/CMakeFiles/utils.dir/logger/logger_impl.cpp.o.provides.build
.PHONY : common/utils/CMakeFiles/utils.dir/logger/logger_impl.cpp.o.provides

common/utils/CMakeFiles/utils.dir/logger/logger_impl.cpp.o.provides.build: common/utils/CMakeFiles/utils.dir/logger/logger_impl.cpp.o


common/utils/CMakeFiles/utils.dir/autocompletion.cpp.o: common/utils/CMakeFiles/utils.dir/flags.make
common/utils/CMakeFiles/utils.dir/autocompletion.cpp.o: /home/chakenal/ontx-pipeline/SPAdes-3.13.1/src/common/utils/autocompletion.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/chakenal/ontx-pipeline/SPAdes-3.13.1/build_spades/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object common/utils/CMakeFiles/utils.dir/autocompletion.cpp.o"
	cd /home/chakenal/ontx-pipeline/SPAdes-3.13.1/build_spades/common/utils && /usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/utils.dir/autocompletion.cpp.o -c /home/chakenal/ontx-pipeline/SPAdes-3.13.1/src/common/utils/autocompletion.cpp

common/utils/CMakeFiles/utils.dir/autocompletion.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/utils.dir/autocompletion.cpp.i"
	cd /home/chakenal/ontx-pipeline/SPAdes-3.13.1/build_spades/common/utils && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/chakenal/ontx-pipeline/SPAdes-3.13.1/src/common/utils/autocompletion.cpp > CMakeFiles/utils.dir/autocompletion.cpp.i

common/utils/CMakeFiles/utils.dir/autocompletion.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/utils.dir/autocompletion.cpp.s"
	cd /home/chakenal/ontx-pipeline/SPAdes-3.13.1/build_spades/common/utils && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/chakenal/ontx-pipeline/SPAdes-3.13.1/src/common/utils/autocompletion.cpp -o CMakeFiles/utils.dir/autocompletion.cpp.s

common/utils/CMakeFiles/utils.dir/autocompletion.cpp.o.requires:

.PHONY : common/utils/CMakeFiles/utils.dir/autocompletion.cpp.o.requires

common/utils/CMakeFiles/utils.dir/autocompletion.cpp.o.provides: common/utils/CMakeFiles/utils.dir/autocompletion.cpp.o.requires
	$(MAKE) -f common/utils/CMakeFiles/utils.dir/build.make common/utils/CMakeFiles/utils.dir/autocompletion.cpp.o.provides.build
.PHONY : common/utils/CMakeFiles/utils.dir/autocompletion.cpp.o.provides

common/utils/CMakeFiles/utils.dir/autocompletion.cpp.o.provides.build: common/utils/CMakeFiles/utils.dir/autocompletion.cpp.o


# Object files for target utils
utils_OBJECTS = \
"CMakeFiles/utils.dir/memory_limit.cpp.o" \
"CMakeFiles/utils.dir/filesystem/copy_file.cpp.o" \
"CMakeFiles/utils.dir/filesystem/path_helper.cpp.o" \
"CMakeFiles/utils.dir/filesystem/temporary.cpp.o" \
"CMakeFiles/utils.dir/logger/logger_impl.cpp.o" \
"CMakeFiles/utils.dir/autocompletion.cpp.o"

# External object files for target utils
utils_EXTERNAL_OBJECTS =

common/utils/libutils.a: common/utils/CMakeFiles/utils.dir/memory_limit.cpp.o
common/utils/libutils.a: common/utils/CMakeFiles/utils.dir/filesystem/copy_file.cpp.o
common/utils/libutils.a: common/utils/CMakeFiles/utils.dir/filesystem/path_helper.cpp.o
common/utils/libutils.a: common/utils/CMakeFiles/utils.dir/filesystem/temporary.cpp.o
common/utils/libutils.a: common/utils/CMakeFiles/utils.dir/logger/logger_impl.cpp.o
common/utils/libutils.a: common/utils/CMakeFiles/utils.dir/autocompletion.cpp.o
common/utils/libutils.a: common/utils/CMakeFiles/utils.dir/build.make
common/utils/libutils.a: common/utils/CMakeFiles/utils.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/chakenal/ontx-pipeline/SPAdes-3.13.1/build_spades/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Linking CXX static library libutils.a"
	cd /home/chakenal/ontx-pipeline/SPAdes-3.13.1/build_spades/common/utils && $(CMAKE_COMMAND) -P CMakeFiles/utils.dir/cmake_clean_target.cmake
	cd /home/chakenal/ontx-pipeline/SPAdes-3.13.1/build_spades/common/utils && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/utils.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
common/utils/CMakeFiles/utils.dir/build: common/utils/libutils.a

.PHONY : common/utils/CMakeFiles/utils.dir/build

common/utils/CMakeFiles/utils.dir/requires: common/utils/CMakeFiles/utils.dir/memory_limit.cpp.o.requires
common/utils/CMakeFiles/utils.dir/requires: common/utils/CMakeFiles/utils.dir/filesystem/copy_file.cpp.o.requires
common/utils/CMakeFiles/utils.dir/requires: common/utils/CMakeFiles/utils.dir/filesystem/path_helper.cpp.o.requires
common/utils/CMakeFiles/utils.dir/requires: common/utils/CMakeFiles/utils.dir/filesystem/temporary.cpp.o.requires
common/utils/CMakeFiles/utils.dir/requires: common/utils/CMakeFiles/utils.dir/logger/logger_impl.cpp.o.requires
common/utils/CMakeFiles/utils.dir/requires: common/utils/CMakeFiles/utils.dir/autocompletion.cpp.o.requires

.PHONY : common/utils/CMakeFiles/utils.dir/requires

common/utils/CMakeFiles/utils.dir/clean:
	cd /home/chakenal/ontx-pipeline/SPAdes-3.13.1/build_spades/common/utils && $(CMAKE_COMMAND) -P CMakeFiles/utils.dir/cmake_clean.cmake
.PHONY : common/utils/CMakeFiles/utils.dir/clean

common/utils/CMakeFiles/utils.dir/depend:
	cd /home/chakenal/ontx-pipeline/SPAdes-3.13.1/build_spades && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/chakenal/ontx-pipeline/SPAdes-3.13.1/src /home/chakenal/ontx-pipeline/SPAdes-3.13.1/src/common/utils /home/chakenal/ontx-pipeline/SPAdes-3.13.1/build_spades /home/chakenal/ontx-pipeline/SPAdes-3.13.1/build_spades/common/utils /home/chakenal/ontx-pipeline/SPAdes-3.13.1/build_spades/common/utils/CMakeFiles/utils.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : common/utils/CMakeFiles/utils.dir/depend

