# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.16

# Note that incremental build could trigger a call to cmake_copy_f90_mod on each re-build

CMakeFiles/diffusion_1d.exe.dir/src/history_mod.F90.o: CMakeFiles/diffusion_1d.exe.dir/mesh_mod.mod.stamp
CMakeFiles/diffusion_1d.exe.dir/src/history_mod.F90.o: CMakeFiles/diffusion_1d.exe.dir/namelist_mod.mod.stamp
CMakeFiles/diffusion_1d.exe.dir/src/history_mod.F90.o: /usr/local/starman/clang_10.0.0/Packages/netcdf/4.6.0/include/netcdf.mod
CMakeFiles/diffusion_1d.exe.dir/src/history_mod.F90.o: CMakeFiles/diffusion_1d.exe.dir/state_mod.mod.stamp
CMakeFiles/diffusion_1d.exe.dir/src/history_mod.F90.o.provides.build: CMakeFiles/diffusion_1d.exe.dir/history_mod.mod.stamp
CMakeFiles/diffusion_1d.exe.dir/history_mod.mod.stamp: CMakeFiles/diffusion_1d.exe.dir/src/history_mod.F90.o
	$(CMAKE_COMMAND) -E cmake_copy_f90_mod .mods/history_mod.mod CMakeFiles/diffusion_1d.exe.dir/history_mod.mod.stamp GNU
CMakeFiles/diffusion_1d.exe.dir/src/history_mod.F90.o.provides.build:
	$(CMAKE_COMMAND) -E touch CMakeFiles/diffusion_1d.exe.dir/src/history_mod.F90.o.provides.build
CMakeFiles/diffusion_1d.exe.dir/build: CMakeFiles/diffusion_1d.exe.dir/src/history_mod.F90.o.provides.build

CMakeFiles/diffusion_1d.exe.dir/src/main.F90.o: CMakeFiles/diffusion_1d.exe.dir/history_mod.mod.stamp
CMakeFiles/diffusion_1d.exe.dir/src/main.F90.o: CMakeFiles/diffusion_1d.exe.dir/mesh_mod.mod.stamp
CMakeFiles/diffusion_1d.exe.dir/src/main.F90.o: CMakeFiles/diffusion_1d.exe.dir/namelist_mod.mod.stamp
CMakeFiles/diffusion_1d.exe.dir/src/main.F90.o: CMakeFiles/diffusion_1d.exe.dir/scheme_mod.mod.stamp
CMakeFiles/diffusion_1d.exe.dir/src/main.F90.o: CMakeFiles/diffusion_1d.exe.dir/state_mod.mod.stamp

CMakeFiles/diffusion_1d.exe.dir/src/mesh_mod.F90.o: CMakeFiles/diffusion_1d.exe.dir/namelist_mod.mod.stamp
CMakeFiles/diffusion_1d.exe.dir/src/mesh_mod.F90.o.provides.build: CMakeFiles/diffusion_1d.exe.dir/mesh_mod.mod.stamp
CMakeFiles/diffusion_1d.exe.dir/mesh_mod.mod.stamp: CMakeFiles/diffusion_1d.exe.dir/src/mesh_mod.F90.o
	$(CMAKE_COMMAND) -E cmake_copy_f90_mod .mods/mesh_mod.mod CMakeFiles/diffusion_1d.exe.dir/mesh_mod.mod.stamp GNU
CMakeFiles/diffusion_1d.exe.dir/src/mesh_mod.F90.o.provides.build:
	$(CMAKE_COMMAND) -E touch CMakeFiles/diffusion_1d.exe.dir/src/mesh_mod.F90.o.provides.build
CMakeFiles/diffusion_1d.exe.dir/build: CMakeFiles/diffusion_1d.exe.dir/src/mesh_mod.F90.o.provides.build

CMakeFiles/diffusion_1d.exe.dir/src/namelist_mod.F90.o.provides.build: CMakeFiles/diffusion_1d.exe.dir/namelist_mod.mod.stamp
CMakeFiles/diffusion_1d.exe.dir/namelist_mod.mod.stamp: CMakeFiles/diffusion_1d.exe.dir/src/namelist_mod.F90.o
	$(CMAKE_COMMAND) -E cmake_copy_f90_mod .mods/namelist_mod.mod CMakeFiles/diffusion_1d.exe.dir/namelist_mod.mod.stamp GNU
CMakeFiles/diffusion_1d.exe.dir/src/namelist_mod.F90.o.provides.build:
	$(CMAKE_COMMAND) -E touch CMakeFiles/diffusion_1d.exe.dir/src/namelist_mod.F90.o.provides.build
CMakeFiles/diffusion_1d.exe.dir/build: CMakeFiles/diffusion_1d.exe.dir/src/namelist_mod.F90.o.provides.build

CMakeFiles/diffusion_1d.exe.dir/src/scheme_mod.F90.o: CMakeFiles/diffusion_1d.exe.dir/mesh_mod.mod.stamp
CMakeFiles/diffusion_1d.exe.dir/src/scheme_mod.F90.o: CMakeFiles/diffusion_1d.exe.dir/namelist_mod.mod.stamp
CMakeFiles/diffusion_1d.exe.dir/src/scheme_mod.F90.o.provides.build: CMakeFiles/diffusion_1d.exe.dir/scheme_mod.mod.stamp
CMakeFiles/diffusion_1d.exe.dir/scheme_mod.mod.stamp: CMakeFiles/diffusion_1d.exe.dir/src/scheme_mod.F90.o
	$(CMAKE_COMMAND) -E cmake_copy_f90_mod .mods/scheme_mod.mod CMakeFiles/diffusion_1d.exe.dir/scheme_mod.mod.stamp GNU
CMakeFiles/diffusion_1d.exe.dir/src/scheme_mod.F90.o.provides.build:
	$(CMAKE_COMMAND) -E touch CMakeFiles/diffusion_1d.exe.dir/src/scheme_mod.F90.o.provides.build
CMakeFiles/diffusion_1d.exe.dir/build: CMakeFiles/diffusion_1d.exe.dir/src/scheme_mod.F90.o.provides.build

CMakeFiles/diffusion_1d.exe.dir/src/state_mod.F90.o: CMakeFiles/diffusion_1d.exe.dir/mesh_mod.mod.stamp
CMakeFiles/diffusion_1d.exe.dir/src/state_mod.F90.o: CMakeFiles/diffusion_1d.exe.dir/namelist_mod.mod.stamp
CMakeFiles/diffusion_1d.exe.dir/src/state_mod.F90.o: CMakeFiles/diffusion_1d.exe.dir/scheme_mod.mod.stamp
CMakeFiles/diffusion_1d.exe.dir/src/state_mod.F90.o.provides.build: CMakeFiles/diffusion_1d.exe.dir/state_mod.mod.stamp
CMakeFiles/diffusion_1d.exe.dir/state_mod.mod.stamp: CMakeFiles/diffusion_1d.exe.dir/src/state_mod.F90.o
	$(CMAKE_COMMAND) -E cmake_copy_f90_mod .mods/state_mod.mod CMakeFiles/diffusion_1d.exe.dir/state_mod.mod.stamp GNU
CMakeFiles/diffusion_1d.exe.dir/src/state_mod.F90.o.provides.build:
	$(CMAKE_COMMAND) -E touch CMakeFiles/diffusion_1d.exe.dir/src/state_mod.F90.o.provides.build
CMakeFiles/diffusion_1d.exe.dir/build: CMakeFiles/diffusion_1d.exe.dir/src/state_mod.F90.o.provides.build
