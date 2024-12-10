set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(PYBIND11_NEWPYTHON ON)

project(hermespy-rt)

find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
find_package(pybind11 CONFIG REQUIRED)

# Configure the C++ library
add_library(hermespy-rt-static STATIC compute_paths.c)
target_include_directories(hermespy-rt-static PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

# Configure the python module
pybind11_add_module(hermespy-rt-binding MODULE compute_paths_pybind11.cpp)
set_target_properties(hermespy-rt-binding PROPERTIES OUTPUT_NAME "rt")
target_link_libraries(hermespy-rt-binding PRIVATE hermespy-rt-static pybind11::headers)
target_include_directories(hermespy-rt-binding PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
