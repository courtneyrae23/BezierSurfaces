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
CMAKE_SOURCE_DIR = /home/cc/cs184/fa16/class/cs184-ado/example_03

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/cc/cs184/fa16/class/cs184-ado/example_03/build

# Include any dependencies generated for this target.
include glew/CMakeFiles/glew.dir/depend.make

# Include the progress variables for this target.
include glew/CMakeFiles/glew.dir/progress.make

# Include the compile flags for this target's objects.
include glew/CMakeFiles/glew.dir/flags.make

glew/CMakeFiles/glew.dir/src/glew.c.o: glew/CMakeFiles/glew.dir/flags.make
glew/CMakeFiles/glew.dir/src/glew.c.o: ../glew/src/glew.c
	$(CMAKE_COMMAND) -E cmake_progress_report /home/cc/cs184/fa16/class/cs184-ado/example_03/build/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object glew/CMakeFiles/glew.dir/src/glew.c.o"
	cd /home/cc/cs184/fa16/class/cs184-ado/example_03/build/glew && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/glew.dir/src/glew.c.o   -c /home/cc/cs184/fa16/class/cs184-ado/example_03/glew/src/glew.c

glew/CMakeFiles/glew.dir/src/glew.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/glew.dir/src/glew.c.i"
	cd /home/cc/cs184/fa16/class/cs184-ado/example_03/build/glew && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -E /home/cc/cs184/fa16/class/cs184-ado/example_03/glew/src/glew.c > CMakeFiles/glew.dir/src/glew.c.i

glew/CMakeFiles/glew.dir/src/glew.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/glew.dir/src/glew.c.s"
	cd /home/cc/cs184/fa16/class/cs184-ado/example_03/build/glew && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -S /home/cc/cs184/fa16/class/cs184-ado/example_03/glew/src/glew.c -o CMakeFiles/glew.dir/src/glew.c.s

glew/CMakeFiles/glew.dir/src/glew.c.o.requires:
.PHONY : glew/CMakeFiles/glew.dir/src/glew.c.o.requires

glew/CMakeFiles/glew.dir/src/glew.c.o.provides: glew/CMakeFiles/glew.dir/src/glew.c.o.requires
	$(MAKE) -f glew/CMakeFiles/glew.dir/build.make glew/CMakeFiles/glew.dir/src/glew.c.o.provides.build
.PHONY : glew/CMakeFiles/glew.dir/src/glew.c.o.provides

glew/CMakeFiles/glew.dir/src/glew.c.o.provides.build: glew/CMakeFiles/glew.dir/src/glew.c.o

glew/CMakeFiles/glew.dir/src/glewinfo.c.o: glew/CMakeFiles/glew.dir/flags.make
glew/CMakeFiles/glew.dir/src/glewinfo.c.o: ../glew/src/glewinfo.c
	$(CMAKE_COMMAND) -E cmake_progress_report /home/cc/cs184/fa16/class/cs184-ado/example_03/build/CMakeFiles $(CMAKE_PROGRESS_2)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object glew/CMakeFiles/glew.dir/src/glewinfo.c.o"
	cd /home/cc/cs184/fa16/class/cs184-ado/example_03/build/glew && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/glew.dir/src/glewinfo.c.o   -c /home/cc/cs184/fa16/class/cs184-ado/example_03/glew/src/glewinfo.c

glew/CMakeFiles/glew.dir/src/glewinfo.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/glew.dir/src/glewinfo.c.i"
	cd /home/cc/cs184/fa16/class/cs184-ado/example_03/build/glew && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -E /home/cc/cs184/fa16/class/cs184-ado/example_03/glew/src/glewinfo.c > CMakeFiles/glew.dir/src/glewinfo.c.i

glew/CMakeFiles/glew.dir/src/glewinfo.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/glew.dir/src/glewinfo.c.s"
	cd /home/cc/cs184/fa16/class/cs184-ado/example_03/build/glew && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -S /home/cc/cs184/fa16/class/cs184-ado/example_03/glew/src/glewinfo.c -o CMakeFiles/glew.dir/src/glewinfo.c.s

glew/CMakeFiles/glew.dir/src/glewinfo.c.o.requires:
.PHONY : glew/CMakeFiles/glew.dir/src/glewinfo.c.o.requires

glew/CMakeFiles/glew.dir/src/glewinfo.c.o.provides: glew/CMakeFiles/glew.dir/src/glewinfo.c.o.requires
	$(MAKE) -f glew/CMakeFiles/glew.dir/build.make glew/CMakeFiles/glew.dir/src/glewinfo.c.o.provides.build
.PHONY : glew/CMakeFiles/glew.dir/src/glewinfo.c.o.provides

glew/CMakeFiles/glew.dir/src/glewinfo.c.o.provides.build: glew/CMakeFiles/glew.dir/src/glewinfo.c.o

glew/CMakeFiles/glew.dir/src/visualinfo.c.o: glew/CMakeFiles/glew.dir/flags.make
glew/CMakeFiles/glew.dir/src/visualinfo.c.o: ../glew/src/visualinfo.c
	$(CMAKE_COMMAND) -E cmake_progress_report /home/cc/cs184/fa16/class/cs184-ado/example_03/build/CMakeFiles $(CMAKE_PROGRESS_3)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object glew/CMakeFiles/glew.dir/src/visualinfo.c.o"
	cd /home/cc/cs184/fa16/class/cs184-ado/example_03/build/glew && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/glew.dir/src/visualinfo.c.o   -c /home/cc/cs184/fa16/class/cs184-ado/example_03/glew/src/visualinfo.c

glew/CMakeFiles/glew.dir/src/visualinfo.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/glew.dir/src/visualinfo.c.i"
	cd /home/cc/cs184/fa16/class/cs184-ado/example_03/build/glew && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -E /home/cc/cs184/fa16/class/cs184-ado/example_03/glew/src/visualinfo.c > CMakeFiles/glew.dir/src/visualinfo.c.i

glew/CMakeFiles/glew.dir/src/visualinfo.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/glew.dir/src/visualinfo.c.s"
	cd /home/cc/cs184/fa16/class/cs184-ado/example_03/build/glew && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -S /home/cc/cs184/fa16/class/cs184-ado/example_03/glew/src/visualinfo.c -o CMakeFiles/glew.dir/src/visualinfo.c.s

glew/CMakeFiles/glew.dir/src/visualinfo.c.o.requires:
.PHONY : glew/CMakeFiles/glew.dir/src/visualinfo.c.o.requires

glew/CMakeFiles/glew.dir/src/visualinfo.c.o.provides: glew/CMakeFiles/glew.dir/src/visualinfo.c.o.requires
	$(MAKE) -f glew/CMakeFiles/glew.dir/build.make glew/CMakeFiles/glew.dir/src/visualinfo.c.o.provides.build
.PHONY : glew/CMakeFiles/glew.dir/src/visualinfo.c.o.provides

glew/CMakeFiles/glew.dir/src/visualinfo.c.o.provides.build: glew/CMakeFiles/glew.dir/src/visualinfo.c.o

# Object files for target glew
glew_OBJECTS = \
"CMakeFiles/glew.dir/src/glew.c.o" \
"CMakeFiles/glew.dir/src/glewinfo.c.o" \
"CMakeFiles/glew.dir/src/visualinfo.c.o"

# External object files for target glew
glew_EXTERNAL_OBJECTS =

glew/libglew.a: glew/CMakeFiles/glew.dir/src/glew.c.o
glew/libglew.a: glew/CMakeFiles/glew.dir/src/glewinfo.c.o
glew/libglew.a: glew/CMakeFiles/glew.dir/src/visualinfo.c.o
glew/libglew.a: glew/CMakeFiles/glew.dir/build.make
glew/libglew.a: glew/CMakeFiles/glew.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking C static library libglew.a"
	cd /home/cc/cs184/fa16/class/cs184-ado/example_03/build/glew && $(CMAKE_COMMAND) -P CMakeFiles/glew.dir/cmake_clean_target.cmake
	cd /home/cc/cs184/fa16/class/cs184-ado/example_03/build/glew && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/glew.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
glew/CMakeFiles/glew.dir/build: glew/libglew.a
.PHONY : glew/CMakeFiles/glew.dir/build

glew/CMakeFiles/glew.dir/requires: glew/CMakeFiles/glew.dir/src/glew.c.o.requires
glew/CMakeFiles/glew.dir/requires: glew/CMakeFiles/glew.dir/src/glewinfo.c.o.requires
glew/CMakeFiles/glew.dir/requires: glew/CMakeFiles/glew.dir/src/visualinfo.c.o.requires
.PHONY : glew/CMakeFiles/glew.dir/requires

glew/CMakeFiles/glew.dir/clean:
	cd /home/cc/cs184/fa16/class/cs184-ado/example_03/build/glew && $(CMAKE_COMMAND) -P CMakeFiles/glew.dir/cmake_clean.cmake
.PHONY : glew/CMakeFiles/glew.dir/clean

glew/CMakeFiles/glew.dir/depend:
	cd /home/cc/cs184/fa16/class/cs184-ado/example_03/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/cc/cs184/fa16/class/cs184-ado/example_03 /home/cc/cs184/fa16/class/cs184-ado/example_03/glew /home/cc/cs184/fa16/class/cs184-ado/example_03/build /home/cc/cs184/fa16/class/cs184-ado/example_03/build/glew /home/cc/cs184/fa16/class/cs184-ado/example_03/build/glew/CMakeFiles/glew.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : glew/CMakeFiles/glew.dir/depend

