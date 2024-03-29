cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

project(move)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# libtorch
set(CMAKE_PREFIX_PATH ~/Software/Programming/libtorch-cuda10.1-1.7.1)
find_package(Torch REQUIRED)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${TORCH_CXX_FLAGS}")

# Cpp-Library
include_directories(~/Library/Cpp-Library/include)
link_directories(~/Library/Cpp-Library/lib)

# Torch-Chemistry
include_directories(~/Library/Torch-Chemistry/include)
link_directories(~/Library/Torch-Chemistry/lib)

include_directories(include)
add_executable(move.exe
    source/PES.cpp
    source/basic.cpp
    source/PESintegral.cpp
    source/expectation.cpp
    source/window.cpp
    source/EOM.cpp
    source/main.cpp
)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -qopenmp -parallel -ipo -static-intel")
target_link_libraries(move.exe tchem CL "${TORCH_LIBRARIES}")
