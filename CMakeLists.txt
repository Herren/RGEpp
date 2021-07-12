cmake_minimum_required(VERSION 3.10)
project(RGE++)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Somehow compilers struggle with the templates...
add_compile_options(-ftemplate-depth=4096)

# Add .h files
include_directories(include)
include_directories(models)

# Add external libraries
find_package(PkgConfig)
pkg_search_module(Eigen3 REQUIRED eigen3)
include_directories(${Eigen3_INCLUDE_DIRS})

find_package(Boost 1.70.0 REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

# Define all examples
add_executable(nsm_example examples/nsm_example.cpp models/nsm.cpp src/ckm.cpp)
add_executable(sm_example examples/sm_example.cpp models/sm.cpp src/ckm.cpp src/pmns.cpp)
add_executable(smtest examples/smtest.cpp models/sm.cpp src/ckm.cpp src/pmns.cpp)
add_executable(thdm_example examples/thdm_example.cpp models/thdmi.cpp models/thdmii.cpp models/thdmx.cpp models/thdmy.cpp src/ckm.cpp src/pmns.cpp)
add_executable(nusmtest examples/nusmtest.cpp models/nusm.cpp examples/rundown.cpp src/ckm.cpp src/pmns.cpp)
add_executable(running_plot examples/running_plot.cpp models/sm.cpp src/ckm.cpp src/pmns.cpp)
add_executable(landau_example examples/landau_example.cpp models/sm.cpp src/ckm.cpp src/pmns.cpp)
add_executable(mssmtest examples/mssmtest.cpp models/mssm.cpp src/ckm.cpp src/pmns.cpp)
add_executable(numssmtest examples/numssmtest.cpp examples/rundown.cpp models/numssm.cpp src/ckm.cpp src/pmns.cpp)
add_executable(numssm_example examples/numssm_example.cpp examples/rundown.cpp models/numssm.cpp src/ckm.cpp src/pmns.cpp)
