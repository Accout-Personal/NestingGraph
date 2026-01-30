## Building
This project uses CMake for building the binaries.

This project requires a compiler compatible with the C++23 feature set.

To build the project first clone the repository and pull all the submodules:
```
git clone https://github.com/Accout-Personal/NestingGraph
cd NestingGraph
git submodule update --init --recursive

```

For linux
(yet to be updated)

Then create a build folder and initialize CMake files:
```
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
```
use `cmake -DCMAKE_BUILD_TYPE=Debug ..` for Debug builds (on linux).

Then to build the executable run:
 - on Linux -> `cmake --build .`
   The binaries is located in the following path: `build/nesting_graph`, `build/nfpGen`.
 - on Windows -> `cmake --build . --config Release`. The binaries are located in the following path: `build/Release/nesting_graph.exe`, `build/Release/nfpGen.exe`
   
   For debug mode use `cmake --build . --config Debug` (Or without --config option ). The binaries are located in the following path: `build/Debug/nesting_graph.exe`, `build/Debug/nfpGen.exe`
