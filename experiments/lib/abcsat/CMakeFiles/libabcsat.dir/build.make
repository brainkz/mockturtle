# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.25

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
CMAKE_COMMAND = /opt/homebrew/Cellar/cmake/3.25.2/bin/cmake

# The command to remove a file.
RM = /opt/homebrew/Cellar/cmake/3.25.2/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/brainkz/Documents/GitHub/mockturtle

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/brainkz/Documents/GitHub/mockturtle/experiments

# Include any dependencies generated for this target.
include lib/abcsat/CMakeFiles/libabcsat.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include lib/abcsat/CMakeFiles/libabcsat.dir/compiler_depend.make

# Include the progress variables for this target.
include lib/abcsat/CMakeFiles/libabcsat.dir/progress.make

# Include the compile flags for this target's objects.
include lib/abcsat/CMakeFiles/libabcsat.dir/flags.make

lib/abcsat/CMakeFiles/libabcsat.dir/AbcGlucose.cpp.o: lib/abcsat/CMakeFiles/libabcsat.dir/flags.make
lib/abcsat/CMakeFiles/libabcsat.dir/AbcGlucose.cpp.o: /Users/brainkz/Documents/GitHub/mockturtle/lib/abcsat/AbcGlucose.cpp
lib/abcsat/CMakeFiles/libabcsat.dir/AbcGlucose.cpp.o: lib/abcsat/CMakeFiles/libabcsat.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/brainkz/Documents/GitHub/mockturtle/experiments/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object lib/abcsat/CMakeFiles/libabcsat.dir/AbcGlucose.cpp.o"
	cd /Users/brainkz/Documents/GitHub/mockturtle/experiments/lib/abcsat && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT lib/abcsat/CMakeFiles/libabcsat.dir/AbcGlucose.cpp.o -MF CMakeFiles/libabcsat.dir/AbcGlucose.cpp.o.d -o CMakeFiles/libabcsat.dir/AbcGlucose.cpp.o -c /Users/brainkz/Documents/GitHub/mockturtle/lib/abcsat/AbcGlucose.cpp

lib/abcsat/CMakeFiles/libabcsat.dir/AbcGlucose.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/libabcsat.dir/AbcGlucose.cpp.i"
	cd /Users/brainkz/Documents/GitHub/mockturtle/experiments/lib/abcsat && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/brainkz/Documents/GitHub/mockturtle/lib/abcsat/AbcGlucose.cpp > CMakeFiles/libabcsat.dir/AbcGlucose.cpp.i

lib/abcsat/CMakeFiles/libabcsat.dir/AbcGlucose.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/libabcsat.dir/AbcGlucose.cpp.s"
	cd /Users/brainkz/Documents/GitHub/mockturtle/experiments/lib/abcsat && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/brainkz/Documents/GitHub/mockturtle/lib/abcsat/AbcGlucose.cpp -o CMakeFiles/libabcsat.dir/AbcGlucose.cpp.s

lib/abcsat/CMakeFiles/libabcsat.dir/Glucose.cpp.o: lib/abcsat/CMakeFiles/libabcsat.dir/flags.make
lib/abcsat/CMakeFiles/libabcsat.dir/Glucose.cpp.o: /Users/brainkz/Documents/GitHub/mockturtle/lib/abcsat/Glucose.cpp
lib/abcsat/CMakeFiles/libabcsat.dir/Glucose.cpp.o: lib/abcsat/CMakeFiles/libabcsat.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/brainkz/Documents/GitHub/mockturtle/experiments/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object lib/abcsat/CMakeFiles/libabcsat.dir/Glucose.cpp.o"
	cd /Users/brainkz/Documents/GitHub/mockturtle/experiments/lib/abcsat && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT lib/abcsat/CMakeFiles/libabcsat.dir/Glucose.cpp.o -MF CMakeFiles/libabcsat.dir/Glucose.cpp.o.d -o CMakeFiles/libabcsat.dir/Glucose.cpp.o -c /Users/brainkz/Documents/GitHub/mockturtle/lib/abcsat/Glucose.cpp

lib/abcsat/CMakeFiles/libabcsat.dir/Glucose.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/libabcsat.dir/Glucose.cpp.i"
	cd /Users/brainkz/Documents/GitHub/mockturtle/experiments/lib/abcsat && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/brainkz/Documents/GitHub/mockturtle/lib/abcsat/Glucose.cpp > CMakeFiles/libabcsat.dir/Glucose.cpp.i

lib/abcsat/CMakeFiles/libabcsat.dir/Glucose.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/libabcsat.dir/Glucose.cpp.s"
	cd /Users/brainkz/Documents/GitHub/mockturtle/experiments/lib/abcsat && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/brainkz/Documents/GitHub/mockturtle/lib/abcsat/Glucose.cpp -o CMakeFiles/libabcsat.dir/Glucose.cpp.s

lib/abcsat/CMakeFiles/libabcsat.dir/SimpSolver.cpp.o: lib/abcsat/CMakeFiles/libabcsat.dir/flags.make
lib/abcsat/CMakeFiles/libabcsat.dir/SimpSolver.cpp.o: /Users/brainkz/Documents/GitHub/mockturtle/lib/abcsat/SimpSolver.cpp
lib/abcsat/CMakeFiles/libabcsat.dir/SimpSolver.cpp.o: lib/abcsat/CMakeFiles/libabcsat.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/brainkz/Documents/GitHub/mockturtle/experiments/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object lib/abcsat/CMakeFiles/libabcsat.dir/SimpSolver.cpp.o"
	cd /Users/brainkz/Documents/GitHub/mockturtle/experiments/lib/abcsat && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT lib/abcsat/CMakeFiles/libabcsat.dir/SimpSolver.cpp.o -MF CMakeFiles/libabcsat.dir/SimpSolver.cpp.o.d -o CMakeFiles/libabcsat.dir/SimpSolver.cpp.o -c /Users/brainkz/Documents/GitHub/mockturtle/lib/abcsat/SimpSolver.cpp

lib/abcsat/CMakeFiles/libabcsat.dir/SimpSolver.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/libabcsat.dir/SimpSolver.cpp.i"
	cd /Users/brainkz/Documents/GitHub/mockturtle/experiments/lib/abcsat && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/brainkz/Documents/GitHub/mockturtle/lib/abcsat/SimpSolver.cpp > CMakeFiles/libabcsat.dir/SimpSolver.cpp.i

lib/abcsat/CMakeFiles/libabcsat.dir/SimpSolver.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/libabcsat.dir/SimpSolver.cpp.s"
	cd /Users/brainkz/Documents/GitHub/mockturtle/experiments/lib/abcsat && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/brainkz/Documents/GitHub/mockturtle/lib/abcsat/SimpSolver.cpp -o CMakeFiles/libabcsat.dir/SimpSolver.cpp.s

lib/abcsat/CMakeFiles/libabcsat.dir/satSolver.cpp.o: lib/abcsat/CMakeFiles/libabcsat.dir/flags.make
lib/abcsat/CMakeFiles/libabcsat.dir/satSolver.cpp.o: /Users/brainkz/Documents/GitHub/mockturtle/lib/abcsat/satSolver.cpp
lib/abcsat/CMakeFiles/libabcsat.dir/satSolver.cpp.o: lib/abcsat/CMakeFiles/libabcsat.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/brainkz/Documents/GitHub/mockturtle/experiments/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object lib/abcsat/CMakeFiles/libabcsat.dir/satSolver.cpp.o"
	cd /Users/brainkz/Documents/GitHub/mockturtle/experiments/lib/abcsat && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT lib/abcsat/CMakeFiles/libabcsat.dir/satSolver.cpp.o -MF CMakeFiles/libabcsat.dir/satSolver.cpp.o.d -o CMakeFiles/libabcsat.dir/satSolver.cpp.o -c /Users/brainkz/Documents/GitHub/mockturtle/lib/abcsat/satSolver.cpp

lib/abcsat/CMakeFiles/libabcsat.dir/satSolver.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/libabcsat.dir/satSolver.cpp.i"
	cd /Users/brainkz/Documents/GitHub/mockturtle/experiments/lib/abcsat && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/brainkz/Documents/GitHub/mockturtle/lib/abcsat/satSolver.cpp > CMakeFiles/libabcsat.dir/satSolver.cpp.i

lib/abcsat/CMakeFiles/libabcsat.dir/satSolver.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/libabcsat.dir/satSolver.cpp.s"
	cd /Users/brainkz/Documents/GitHub/mockturtle/experiments/lib/abcsat && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/brainkz/Documents/GitHub/mockturtle/lib/abcsat/satSolver.cpp -o CMakeFiles/libabcsat.dir/satSolver.cpp.s

lib/abcsat/CMakeFiles/libabcsat.dir/satStore.cpp.o: lib/abcsat/CMakeFiles/libabcsat.dir/flags.make
lib/abcsat/CMakeFiles/libabcsat.dir/satStore.cpp.o: /Users/brainkz/Documents/GitHub/mockturtle/lib/abcsat/satStore.cpp
lib/abcsat/CMakeFiles/libabcsat.dir/satStore.cpp.o: lib/abcsat/CMakeFiles/libabcsat.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/brainkz/Documents/GitHub/mockturtle/experiments/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object lib/abcsat/CMakeFiles/libabcsat.dir/satStore.cpp.o"
	cd /Users/brainkz/Documents/GitHub/mockturtle/experiments/lib/abcsat && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT lib/abcsat/CMakeFiles/libabcsat.dir/satStore.cpp.o -MF CMakeFiles/libabcsat.dir/satStore.cpp.o.d -o CMakeFiles/libabcsat.dir/satStore.cpp.o -c /Users/brainkz/Documents/GitHub/mockturtle/lib/abcsat/satStore.cpp

lib/abcsat/CMakeFiles/libabcsat.dir/satStore.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/libabcsat.dir/satStore.cpp.i"
	cd /Users/brainkz/Documents/GitHub/mockturtle/experiments/lib/abcsat && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/brainkz/Documents/GitHub/mockturtle/lib/abcsat/satStore.cpp > CMakeFiles/libabcsat.dir/satStore.cpp.i

lib/abcsat/CMakeFiles/libabcsat.dir/satStore.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/libabcsat.dir/satStore.cpp.s"
	cd /Users/brainkz/Documents/GitHub/mockturtle/experiments/lib/abcsat && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/brainkz/Documents/GitHub/mockturtle/lib/abcsat/satStore.cpp -o CMakeFiles/libabcsat.dir/satStore.cpp.s

# Object files for target libabcsat
libabcsat_OBJECTS = \
"CMakeFiles/libabcsat.dir/AbcGlucose.cpp.o" \
"CMakeFiles/libabcsat.dir/Glucose.cpp.o" \
"CMakeFiles/libabcsat.dir/SimpSolver.cpp.o" \
"CMakeFiles/libabcsat.dir/satSolver.cpp.o" \
"CMakeFiles/libabcsat.dir/satStore.cpp.o"

# External object files for target libabcsat
libabcsat_EXTERNAL_OBJECTS =

lib/abcsat/liblibabcsat.a: lib/abcsat/CMakeFiles/libabcsat.dir/AbcGlucose.cpp.o
lib/abcsat/liblibabcsat.a: lib/abcsat/CMakeFiles/libabcsat.dir/Glucose.cpp.o
lib/abcsat/liblibabcsat.a: lib/abcsat/CMakeFiles/libabcsat.dir/SimpSolver.cpp.o
lib/abcsat/liblibabcsat.a: lib/abcsat/CMakeFiles/libabcsat.dir/satSolver.cpp.o
lib/abcsat/liblibabcsat.a: lib/abcsat/CMakeFiles/libabcsat.dir/satStore.cpp.o
lib/abcsat/liblibabcsat.a: lib/abcsat/CMakeFiles/libabcsat.dir/build.make
lib/abcsat/liblibabcsat.a: lib/abcsat/CMakeFiles/libabcsat.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/brainkz/Documents/GitHub/mockturtle/experiments/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Linking CXX static library liblibabcsat.a"
	cd /Users/brainkz/Documents/GitHub/mockturtle/experiments/lib/abcsat && $(CMAKE_COMMAND) -P CMakeFiles/libabcsat.dir/cmake_clean_target.cmake
	cd /Users/brainkz/Documents/GitHub/mockturtle/experiments/lib/abcsat && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/libabcsat.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
lib/abcsat/CMakeFiles/libabcsat.dir/build: lib/abcsat/liblibabcsat.a
.PHONY : lib/abcsat/CMakeFiles/libabcsat.dir/build

lib/abcsat/CMakeFiles/libabcsat.dir/clean:
	cd /Users/brainkz/Documents/GitHub/mockturtle/experiments/lib/abcsat && $(CMAKE_COMMAND) -P CMakeFiles/libabcsat.dir/cmake_clean.cmake
.PHONY : lib/abcsat/CMakeFiles/libabcsat.dir/clean

lib/abcsat/CMakeFiles/libabcsat.dir/depend:
	cd /Users/brainkz/Documents/GitHub/mockturtle/experiments && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/brainkz/Documents/GitHub/mockturtle /Users/brainkz/Documents/GitHub/mockturtle/lib/abcsat /Users/brainkz/Documents/GitHub/mockturtle/experiments /Users/brainkz/Documents/GitHub/mockturtle/experiments/lib/abcsat /Users/brainkz/Documents/GitHub/mockturtle/experiments/lib/abcsat/CMakeFiles/libabcsat.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : lib/abcsat/CMakeFiles/libabcsat.dir/depend

