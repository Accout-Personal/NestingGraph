set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_EXECUTABLE_ENABLE_EXPORTS OFF CACHE BOOL "" FORCE)
set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS OFF CACHE BOOL "" FORCE)
cmake_minimum_required(VERSION 3.20)
if(WIN32)
    # Windows-specific static settings
    # This bakes the C++ runtime (msvcrt) into the exe
    set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
    add_definitions(-DCGAL_USE_STATIC_LIBS=1)
    set(VCPKG_TARGET_TRIPLET "x64-windows-static" CACHE STRING "")
    set(VCPKG_HOST_TRIPLET "x64-windows-static" CACHE STRING "")
    
else()
    # Linux/macOS static settings (vcpkg usually uses 'x64-linux')
    # Note: Static linking on Linux is handled slightly differently via compiler flags
    set(VCPKG_TARGET_TRIPLET "x64-linux" CACHE STRING "")
endif()

project(nestab LANGUAGES CXX)
set(Boost_USE_STATIC_LIBS ON)
set(Boost_USE_STATIC_RUNTIME ON)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


add_library(utils STATIC
  src/common/geometry_util.cpp
  src/common/nfp_lib.cpp
)

find_package(CGAL CONFIG REQUIRED)

add_library(nfp_cgal SHARED
  src/nfp/NFPGenerator.cpp
)

# export symbols automatically for now (step-1 quick win)
set_target_properties(nfp_cgal PROPERTIES
  WINDOWS_EXPORT_ALL_SYMBOLS ON
)

target_include_directories(nfp_cgal PUBLIC
  src
)

# IMPORTANT: export when building DLL
target_compile_definitions(nfp_cgal PRIVATE NFP_BRIDGE_BUILD)

target_link_libraries(nfp_cgal PRIVATE CGAL::CGAL)

if(TARGET CGAL::CGAL_Core)
  target_link_libraries(nfp_cgal PRIVATE CGAL::CGAL_Core)
endif()

target_include_directories(utils PRIVATE
  src
  ${CMAKE_SOURCE_DIR}/third_party/json/include
)

set(JSON_BuildTests OFF CACHE INTERNAL "")
set(JSON_Install OFF CACHE INTERNAL "")
add_subdirectory(third_party/json)

add_executable(nfpGen
  src/NFPGen.cpp
)

add_executable(nesting_graph
  src/NestingGraph.cpp
)

add_executable(nfpTest
  src/NFPTest.cpp
)

add_executable(graphCount
  src/graphCount.cpp
)

set_target_properties(nfpGen PROPERTIES
  ENABLE_EXPORTS OFF
)

target_include_directories(nfpTest PRIVATE
  src
  ${CMAKE_SOURCE_DIR}/third_party/json/include
)

target_include_directories(nfpGen PRIVATE
  src
  ${CMAKE_SOURCE_DIR}/third_party/json/include
)

target_include_directories(nesting_graph PRIVATE
  src
  ${CMAKE_SOURCE_DIR}/third_party/json/include
)

target_include_directories(graphCount PRIVATE
  src
  ${CMAKE_SOURCE_DIR}/third_party/json/include
)


target_link_libraries(nfpTest PRIVATE
  utils
  nfp_cgal
  nlohmann_json::nlohmann_json
)

target_link_libraries(nfpGen PRIVATE
  utils
  nfp_cgal
  nlohmann_json::nlohmann_json
)

target_link_libraries(nesting_graph PRIVATE
  utils
  nfp_cgal
  nlohmann_json::nlohmann_json
)

target_link_libraries(graphCount PRIVATE
  utils
  nlohmann_json::nlohmann_json
)
