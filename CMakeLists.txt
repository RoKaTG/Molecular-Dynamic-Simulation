cmake_minimum_required(VERSION 3.10)
project(md_simulator)

# Force C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# Aggressive flag
set(CMAKE_CXX_FLAGS_RELEASE
    "-g3 -Ofast -march=native -funroll-loops -flto -ffast-math -fstrict-aliasing -DNDEBUG ${CMAKE_CXX_FLAGS_RELEASE}"
)

# Dependencies
find_package(MPI REQUIRED)
find_package(OpenMP REQUIRED)
find_package(Qt5 COMPONENTS Widgets REQUIRED)

# Includes
include_directories(${MPI_INCLUDE_PATH})
include_directories(${CMAKE_SOURCE_DIR}/src)
include_directories(${CMAKE_SOURCE_DIR}/include)
include_directories(${Qt5Widgets_INCLUDE_DIRS})

# Source files
set(SOURCES
    src/main.cpp
    src/particle.cpp
    src/system.cpp
    src/lennard_jones.cpp
    src/verlet_integrator.cpp
    src/benchmark.cpp
    src/utils.cpp
    src/md_utils.cpp
    # Qt file (soon)
)

# Executable
add_executable(md_simulator ${SOURCES})

# Linker
target_link_libraries(md_simulator 
    PRIVATE 
        MPI::MPI_CXX 
        OpenMP::OpenMP_CXX 
        Qt5::Widgets
)
