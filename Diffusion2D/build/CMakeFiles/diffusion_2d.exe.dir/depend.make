# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.5


CMakeFiles/diffusion_2d.exe.dir/src/diffusion_main.f90.o.requires: CMakeFiles/diffusion_2d.exe.dir/diffusion_scheme.mod.proxy
CMakeFiles/diffusion_2d.exe.dir/src/diffusion_main.f90.o: CMakeFiles/diffusion_2d.exe.dir/diffusion_scheme.mod.stamp

CMakeFiles/diffusion_2d.exe.dir/diffusion_scheme.mod.proxy: CMakeFiles/diffusion_2d.exe.dir/src/scheme_mod.f90.o.provides
CMakeFiles/diffusion_2d.exe.dir/src/scheme_mod.f90.o.provides.build:
	$(CMAKE_COMMAND) -E cmake_copy_f90_mod .mods/diffusion_scheme CMakeFiles/diffusion_2d.exe.dir/diffusion_scheme.mod.stamp GNU
	$(CMAKE_COMMAND) -E touch CMakeFiles/diffusion_2d.exe.dir/src/scheme_mod.f90.o.provides.build
CMakeFiles/diffusion_2d.exe.dir/build: CMakeFiles/diffusion_2d.exe.dir/src/scheme_mod.f90.o.provides.build