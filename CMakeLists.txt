set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_EXECUTABLE_ENABLE_EXPORTS OFF CACHE BOOL "" FORCE)
set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS OFF CACHE BOOL "" FORCE)
cmake_minimum_required(VERSION 3.20)
project(nfp_refactor LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

file(GLOB CLIPPER2_SOURCES
  ${CMAKE_SOURCE_DIR}/third_party/Clipper2/CPP/Clipper2Lib/src/*.cpp
)

# Some versions also have utils in a sibling folder.
# If this folder exists in your tree, enable it:
file(GLOB CLIPPER2_UTILS
  ${CMAKE_SOURCE_DIR}/third_party/Clipper2/CPP/Clipper2Lib/utils/*.cpp
)

add_library(clipper2 STATIC
  ${CLIPPER2_SOURCES}
  ${CLIPPER2_UTILS}
)

add_library(utils STATIC
  src/common/geometry_util.cpp
  src/common/nfp_lib.cpp
)

target_include_directories(utils PRIVATE
  src
  ${CMAKE_SOURCE_DIR}/third_party/Clipper2/CPP/Clipper2Lib/include
  ${CMAKE_SOURCE_DIR}/third_party/json/include
)

target_include_directories(clipper2 PUBLIC
  ${CMAKE_SOURCE_DIR}/third_party/Clipper2/CPP/Clipper2Lib/include
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

set_target_properties(nfpGen PROPERTIES
  ENABLE_EXPORTS OFF
)

target_include_directories(nfpGen PRIVATE
  src
  ${CMAKE_SOURCE_DIR}/third_party/Clipper2/CPP/Clipper2Lib/include
  ${CMAKE_SOURCE_DIR}/third_party/json/include
)

target_include_directories(nesting_graph PRIVATE
  src
  ${CMAKE_SOURCE_DIR}/third_party/Clipper2/CPP/Clipper2Lib/include
  ${CMAKE_SOURCE_DIR}/third_party/json/include
)

target_link_libraries(nfpGen PRIVATE
  utils
  clipper2
  nlohmann_json::nlohmann_json
)

target_link_libraries(nesting_graph PRIVATE
  utils
  clipper2
  nlohmann_json::nlohmann_json
)

