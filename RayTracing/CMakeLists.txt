cmake_minimum_required(VERSION 3.12)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

project(RayTracer)

# Add executable
add_executable(main main.cpp parser/pugixml/src/pugixml.cpp)

# Add include directories
target_include_directories(main PRIVATE source)
target_include_directories(main PRIVATE parser/pugixml/src)

# Set output directory
set_target_properties(main PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/build
)

