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
include glfw-3.2.1/tests/CMakeFiles/tearing.dir/depend.make

# Include the progress variables for this target.
include glfw-3.2.1/tests/CMakeFiles/tearing.dir/progress.make

# Include the compile flags for this target's objects.
include glfw-3.2.1/tests/CMakeFiles/tearing.dir/flags.make

glfw-3.2.1/tests/CMakeFiles/tearing.dir/tearing.c.o: glfw-3.2.1/tests/CMakeFiles/tearing.dir/flags.make
glfw-3.2.1/tests/CMakeFiles/tearing.dir/tearing.c.o: ../glfw-3.2.1/tests/tearing.c
	$(CMAKE_COMMAND) -E cmake_progress_report /home/cc/cs184/fa16/class/cs184-ado/example_03/build/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object glfw-3.2.1/tests/CMakeFiles/tearing.dir/tearing.c.o"
	cd /home/cc/cs184/fa16/class/cs184-ado/example_03/build/glfw-3.2.1/tests && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/tearing.dir/tearing.c.o   -c /home/cc/cs184/fa16/class/cs184-ado/example_03/glfw-3.2.1/tests/tearing.c

glfw-3.2.1/tests/CMakeFiles/tearing.dir/tearing.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/tearing.dir/tearing.c.i"
	cd /home/cc/cs184/fa16/class/cs184-ado/example_03/build/glfw-3.2.1/tests && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -E /home/cc/cs184/fa16/class/cs184-ado/example_03/glfw-3.2.1/tests/tearing.c > CMakeFiles/tearing.dir/tearing.c.i

glfw-3.2.1/tests/CMakeFiles/tearing.dir/tearing.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/tearing.dir/tearing.c.s"
	cd /home/cc/cs184/fa16/class/cs184-ado/example_03/build/glfw-3.2.1/tests && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -S /home/cc/cs184/fa16/class/cs184-ado/example_03/glfw-3.2.1/tests/tearing.c -o CMakeFiles/tearing.dir/tearing.c.s

glfw-3.2.1/tests/CMakeFiles/tearing.dir/tearing.c.o.requires:
.PHONY : glfw-3.2.1/tests/CMakeFiles/tearing.dir/tearing.c.o.requires

glfw-3.2.1/tests/CMakeFiles/tearing.dir/tearing.c.o.provides: glfw-3.2.1/tests/CMakeFiles/tearing.dir/tearing.c.o.requires
	$(MAKE) -f glfw-3.2.1/tests/CMakeFiles/tearing.dir/build.make glfw-3.2.1/tests/CMakeFiles/tearing.dir/tearing.c.o.provides.build
.PHONY : glfw-3.2.1/tests/CMakeFiles/tearing.dir/tearing.c.o.provides

glfw-3.2.1/tests/CMakeFiles/tearing.dir/tearing.c.o.provides.build: glfw-3.2.1/tests/CMakeFiles/tearing.dir/tearing.c.o

glfw-3.2.1/tests/CMakeFiles/tearing.dir/__/deps/getopt.c.o: glfw-3.2.1/tests/CMakeFiles/tearing.dir/flags.make
glfw-3.2.1/tests/CMakeFiles/tearing.dir/__/deps/getopt.c.o: ../glfw-3.2.1/deps/getopt.c
	$(CMAKE_COMMAND) -E cmake_progress_report /home/cc/cs184/fa16/class/cs184-ado/example_03/build/CMakeFiles $(CMAKE_PROGRESS_2)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object glfw-3.2.1/tests/CMakeFiles/tearing.dir/__/deps/getopt.c.o"
	cd /home/cc/cs184/fa16/class/cs184-ado/example_03/build/glfw-3.2.1/tests && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/tearing.dir/__/deps/getopt.c.o   -c /home/cc/cs184/fa16/class/cs184-ado/example_03/glfw-3.2.1/deps/getopt.c

glfw-3.2.1/tests/CMakeFiles/tearing.dir/__/deps/getopt.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/tearing.dir/__/deps/getopt.c.i"
	cd /home/cc/cs184/fa16/class/cs184-ado/example_03/build/glfw-3.2.1/tests && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -E /home/cc/cs184/fa16/class/cs184-ado/example_03/glfw-3.2.1/deps/getopt.c > CMakeFiles/tearing.dir/__/deps/getopt.c.i

glfw-3.2.1/tests/CMakeFiles/tearing.dir/__/deps/getopt.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/tearing.dir/__/deps/getopt.c.s"
	cd /home/cc/cs184/fa16/class/cs184-ado/example_03/build/glfw-3.2.1/tests && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -S /home/cc/cs184/fa16/class/cs184-ado/example_03/glfw-3.2.1/deps/getopt.c -o CMakeFiles/tearing.dir/__/deps/getopt.c.s

glfw-3.2.1/tests/CMakeFiles/tearing.dir/__/deps/getopt.c.o.requires:
.PHONY : glfw-3.2.1/tests/CMakeFiles/tearing.dir/__/deps/getopt.c.o.requires

glfw-3.2.1/tests/CMakeFiles/tearing.dir/__/deps/getopt.c.o.provides: glfw-3.2.1/tests/CMakeFiles/tearing.dir/__/deps/getopt.c.o.requires
	$(MAKE) -f glfw-3.2.1/tests/CMakeFiles/tearing.dir/build.make glfw-3.2.1/tests/CMakeFiles/tearing.dir/__/deps/getopt.c.o.provides.build
.PHONY : glfw-3.2.1/tests/CMakeFiles/tearing.dir/__/deps/getopt.c.o.provides

glfw-3.2.1/tests/CMakeFiles/tearing.dir/__/deps/getopt.c.o.provides.build: glfw-3.2.1/tests/CMakeFiles/tearing.dir/__/deps/getopt.c.o

glfw-3.2.1/tests/CMakeFiles/tearing.dir/__/deps/glad.c.o: glfw-3.2.1/tests/CMakeFiles/tearing.dir/flags.make
glfw-3.2.1/tests/CMakeFiles/tearing.dir/__/deps/glad.c.o: ../glfw-3.2.1/deps/glad.c
	$(CMAKE_COMMAND) -E cmake_progress_report /home/cc/cs184/fa16/class/cs184-ado/example_03/build/CMakeFiles $(CMAKE_PROGRESS_3)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object glfw-3.2.1/tests/CMakeFiles/tearing.dir/__/deps/glad.c.o"
	cd /home/cc/cs184/fa16/class/cs184-ado/example_03/build/glfw-3.2.1/tests && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/tearing.dir/__/deps/glad.c.o   -c /home/cc/cs184/fa16/class/cs184-ado/example_03/glfw-3.2.1/deps/glad.c

glfw-3.2.1/tests/CMakeFiles/tearing.dir/__/deps/glad.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/tearing.dir/__/deps/glad.c.i"
	cd /home/cc/cs184/fa16/class/cs184-ado/example_03/build/glfw-3.2.1/tests && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -E /home/cc/cs184/fa16/class/cs184-ado/example_03/glfw-3.2.1/deps/glad.c > CMakeFiles/tearing.dir/__/deps/glad.c.i

glfw-3.2.1/tests/CMakeFiles/tearing.dir/__/deps/glad.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/tearing.dir/__/deps/glad.c.s"
	cd /home/cc/cs184/fa16/class/cs184-ado/example_03/build/glfw-3.2.1/tests && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -S /home/cc/cs184/fa16/class/cs184-ado/example_03/glfw-3.2.1/deps/glad.c -o CMakeFiles/tearing.dir/__/deps/glad.c.s

glfw-3.2.1/tests/CMakeFiles/tearing.dir/__/deps/glad.c.o.requires:
.PHONY : glfw-3.2.1/tests/CMakeFiles/tearing.dir/__/deps/glad.c.o.requires

glfw-3.2.1/tests/CMakeFiles/tearing.dir/__/deps/glad.c.o.provides: glfw-3.2.1/tests/CMakeFiles/tearing.dir/__/deps/glad.c.o.requires
	$(MAKE) -f glfw-3.2.1/tests/CMakeFiles/tearing.dir/build.make glfw-3.2.1/tests/CMakeFiles/tearing.dir/__/deps/glad.c.o.provides.build
.PHONY : glfw-3.2.1/tests/CMakeFiles/tearing.dir/__/deps/glad.c.o.provides

glfw-3.2.1/tests/CMakeFiles/tearing.dir/__/deps/glad.c.o.provides.build: glfw-3.2.1/tests/CMakeFiles/tearing.dir/__/deps/glad.c.o

# Object files for target tearing
tearing_OBJECTS = \
"CMakeFiles/tearing.dir/tearing.c.o" \
"CMakeFiles/tearing.dir/__/deps/getopt.c.o" \
"CMakeFiles/tearing.dir/__/deps/glad.c.o"

# External object files for target tearing
tearing_EXTERNAL_OBJECTS =

glfw-3.2.1/tests/tearing: glfw-3.2.1/tests/CMakeFiles/tearing.dir/tearing.c.o
glfw-3.2.1/tests/tearing: glfw-3.2.1/tests/CMakeFiles/tearing.dir/__/deps/getopt.c.o
glfw-3.2.1/tests/tearing: glfw-3.2.1/tests/CMakeFiles/tearing.dir/__/deps/glad.c.o
glfw-3.2.1/tests/tearing: glfw-3.2.1/tests/CMakeFiles/tearing.dir/build.make
glfw-3.2.1/tests/tearing: glfw-3.2.1/src/libglfw3.a
glfw-3.2.1/tests/tearing: /usr/lib/x86_64-linux-gnu/librt.so
glfw-3.2.1/tests/tearing: /usr/lib/x86_64-linux-gnu/libm.so
glfw-3.2.1/tests/tearing: /usr/lib/x86_64-linux-gnu/libX11.so
glfw-3.2.1/tests/tearing: /usr/lib/x86_64-linux-gnu/libXrandr.so
glfw-3.2.1/tests/tearing: /usr/lib/x86_64-linux-gnu/libXinerama.so
glfw-3.2.1/tests/tearing: /usr/lib/x86_64-linux-gnu/libXxf86vm.so
glfw-3.2.1/tests/tearing: /usr/lib/x86_64-linux-gnu/libXcursor.so
glfw-3.2.1/tests/tearing: glfw-3.2.1/tests/CMakeFiles/tearing.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking C executable tearing"
	cd /home/cc/cs184/fa16/class/cs184-ado/example_03/build/glfw-3.2.1/tests && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/tearing.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
glfw-3.2.1/tests/CMakeFiles/tearing.dir/build: glfw-3.2.1/tests/tearing
.PHONY : glfw-3.2.1/tests/CMakeFiles/tearing.dir/build

glfw-3.2.1/tests/CMakeFiles/tearing.dir/requires: glfw-3.2.1/tests/CMakeFiles/tearing.dir/tearing.c.o.requires
glfw-3.2.1/tests/CMakeFiles/tearing.dir/requires: glfw-3.2.1/tests/CMakeFiles/tearing.dir/__/deps/getopt.c.o.requires
glfw-3.2.1/tests/CMakeFiles/tearing.dir/requires: glfw-3.2.1/tests/CMakeFiles/tearing.dir/__/deps/glad.c.o.requires
.PHONY : glfw-3.2.1/tests/CMakeFiles/tearing.dir/requires

glfw-3.2.1/tests/CMakeFiles/tearing.dir/clean:
	cd /home/cc/cs184/fa16/class/cs184-ado/example_03/build/glfw-3.2.1/tests && $(CMAKE_COMMAND) -P CMakeFiles/tearing.dir/cmake_clean.cmake
.PHONY : glfw-3.2.1/tests/CMakeFiles/tearing.dir/clean

glfw-3.2.1/tests/CMakeFiles/tearing.dir/depend:
	cd /home/cc/cs184/fa16/class/cs184-ado/example_03/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/cc/cs184/fa16/class/cs184-ado/example_03 /home/cc/cs184/fa16/class/cs184-ado/example_03/glfw-3.2.1/tests /home/cc/cs184/fa16/class/cs184-ado/example_03/build /home/cc/cs184/fa16/class/cs184-ado/example_03/build/glfw-3.2.1/tests /home/cc/cs184/fa16/class/cs184-ado/example_03/build/glfw-3.2.1/tests/CMakeFiles/tearing.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : glfw-3.2.1/tests/CMakeFiles/tearing.dir/depend

