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
CMAKE_SOURCE_DIR = /home/qjz/work/study/diffusion/Diffusion2D

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/qjz/work/study/diffusion/Diffusion2D/build

# Include any dependencies generated for this target.
include CMakeFiles/diffusion_2d.exe.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/diffusion_2d.exe.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/diffusion_2d.exe.dir/flags.make

CMakeFiles/diffusion_2d.exe.dir/src/diffusion_main.f90.o: CMakeFiles/diffusion_2d.exe.dir/flags.make
CMakeFiles/diffusion_2d.exe.dir/src/diffusion_main.f90.o: ../src/diffusion_main.f90
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/qjz/work/study/diffusion/Diffusion2D/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building Fortran object CMakeFiles/diffusion_2d.exe.dir/src/diffusion_main.f90.o"
	/home/qjz/.local/bin/gfortran  $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /home/qjz/work/study/diffusion/Diffusion2D/src/diffusion_main.f90 -o CMakeFiles/diffusion_2d.exe.dir/src/diffusion_main.f90.o

CMakeFiles/diffusion_2d.exe.dir/src/diffusion_main.f90.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/diffusion_2d.exe.dir/src/diffusion_main.f90.i"
	/home/qjz/.local/bin/gfortran  $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /home/qjz/work/study/diffusion/Diffusion2D/src/diffusion_main.f90 > CMakeFiles/diffusion_2d.exe.dir/src/diffusion_main.f90.i

CMakeFiles/diffusion_2d.exe.dir/src/diffusion_main.f90.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/diffusion_2d.exe.dir/src/diffusion_main.f90.s"
	/home/qjz/.local/bin/gfortran  $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /home/qjz/work/study/diffusion/Diffusion2D/src/diffusion_main.f90 -o CMakeFiles/diffusion_2d.exe.dir/src/diffusion_main.f90.s

CMakeFiles/diffusion_2d.exe.dir/src/diffusion_main.f90.o.requires:

.PHONY : CMakeFiles/diffusion_2d.exe.dir/src/diffusion_main.f90.o.requires

CMakeFiles/diffusion_2d.exe.dir/src/diffusion_main.f90.o.provides: CMakeFiles/diffusion_2d.exe.dir/src/diffusion_main.f90.o.requires
	$(MAKE) -f CMakeFiles/diffusion_2d.exe.dir/build.make CMakeFiles/diffusion_2d.exe.dir/src/diffusion_main.f90.o.provides.build
.PHONY : CMakeFiles/diffusion_2d.exe.dir/src/diffusion_main.f90.o.provides

CMakeFiles/diffusion_2d.exe.dir/src/diffusion_main.f90.o.provides.build: CMakeFiles/diffusion_2d.exe.dir/src/diffusion_main.f90.o


CMakeFiles/diffusion_2d.exe.dir/src/scheme_mod.f90.o: CMakeFiles/diffusion_2d.exe.dir/flags.make
CMakeFiles/diffusion_2d.exe.dir/src/scheme_mod.f90.o: ../src/scheme_mod.f90
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/qjz/work/study/diffusion/Diffusion2D/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building Fortran object CMakeFiles/diffusion_2d.exe.dir/src/scheme_mod.f90.o"
	/home/qjz/.local/bin/gfortran  $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /home/qjz/work/study/diffusion/Diffusion2D/src/scheme_mod.f90 -o CMakeFiles/diffusion_2d.exe.dir/src/scheme_mod.f90.o

CMakeFiles/diffusion_2d.exe.dir/src/scheme_mod.f90.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/diffusion_2d.exe.dir/src/scheme_mod.f90.i"
	/home/qjz/.local/bin/gfortran  $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /home/qjz/work/study/diffusion/Diffusion2D/src/scheme_mod.f90 > CMakeFiles/diffusion_2d.exe.dir/src/scheme_mod.f90.i

CMakeFiles/diffusion_2d.exe.dir/src/scheme_mod.f90.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/diffusion_2d.exe.dir/src/scheme_mod.f90.s"
	/home/qjz/.local/bin/gfortran  $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /home/qjz/work/study/diffusion/Diffusion2D/src/scheme_mod.f90 -o CMakeFiles/diffusion_2d.exe.dir/src/scheme_mod.f90.s

CMakeFiles/diffusion_2d.exe.dir/src/scheme_mod.f90.o.requires:

.PHONY : CMakeFiles/diffusion_2d.exe.dir/src/scheme_mod.f90.o.requires

CMakeFiles/diffusion_2d.exe.dir/src/scheme_mod.f90.o.provides: CMakeFiles/diffusion_2d.exe.dir/src/scheme_mod.f90.o.requires
	$(MAKE) -f CMakeFiles/diffusion_2d.exe.dir/build.make CMakeFiles/diffusion_2d.exe.dir/src/scheme_mod.f90.o.provides.build
.PHONY : CMakeFiles/diffusion_2d.exe.dir/src/scheme_mod.f90.o.provides

CMakeFiles/diffusion_2d.exe.dir/src/scheme_mod.f90.o.provides.build: CMakeFiles/diffusion_2d.exe.dir/src/scheme_mod.f90.o


# Object files for target diffusion_2d.exe
diffusion_2d_exe_OBJECTS = \
"CMakeFiles/diffusion_2d.exe.dir/src/diffusion_main.f90.o" \
"CMakeFiles/diffusion_2d.exe.dir/src/scheme_mod.f90.o"

# External object files for target diffusion_2d.exe
diffusion_2d_exe_EXTERNAL_OBJECTS =

diffusion_2d.exe: CMakeFiles/diffusion_2d.exe.dir/src/diffusion_main.f90.o
diffusion_2d.exe: CMakeFiles/diffusion_2d.exe.dir/src/scheme_mod.f90.o
diffusion_2d.exe: CMakeFiles/diffusion_2d.exe.dir/build.make
diffusion_2d.exe: /usr/local/mpich2/lib/libmpifort.so
diffusion_2d.exe: /usr/local/mpich2/lib/libmpi.so
diffusion_2d.exe: CMakeFiles/diffusion_2d.exe.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/qjz/work/study/diffusion/Diffusion2D/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking Fortran executable diffusion_2d.exe"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/diffusion_2d.exe.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/diffusion_2d.exe.dir/build: diffusion_2d.exe

.PHONY : CMakeFiles/diffusion_2d.exe.dir/build

CMakeFiles/diffusion_2d.exe.dir/requires: CMakeFiles/diffusion_2d.exe.dir/src/diffusion_main.f90.o.requires
CMakeFiles/diffusion_2d.exe.dir/requires: CMakeFiles/diffusion_2d.exe.dir/src/scheme_mod.f90.o.requires

.PHONY : CMakeFiles/diffusion_2d.exe.dir/requires

CMakeFiles/diffusion_2d.exe.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/diffusion_2d.exe.dir/cmake_clean.cmake
.PHONY : CMakeFiles/diffusion_2d.exe.dir/clean

CMakeFiles/diffusion_2d.exe.dir/depend:
	cd /home/qjz/work/study/diffusion/Diffusion2D/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/qjz/work/study/diffusion/Diffusion2D /home/qjz/work/study/diffusion/Diffusion2D /home/qjz/work/study/diffusion/Diffusion2D/build /home/qjz/work/study/diffusion/Diffusion2D/build /home/qjz/work/study/diffusion/Diffusion2D/build/CMakeFiles/diffusion_2d.exe.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/diffusion_2d.exe.dir/depend
