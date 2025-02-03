cmake_minimum_required(VERSION 3.29)
project(LunarSim_v0_2)

set(CMAKE_CXX_STANDARD 20)

# Include directories
include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${PROJECT_SOURCE_DIR}/src)
# Include Eigen directory
include_directories("external/eigen-3.4.0")
# Include FFTW3 directory
include_directories("external/fftw-3.3.5-dll64")
link_directories("external/fftw-3.3.5-dll64")

# Source files (excluding main.cpp from src/)
file(GLOB_RECURSE SOURCES src/*.cpp)
file(GLOB_RECURSE HEADERS include/*.h)

# Add main.cpp explicitly
set(SOURCES ${SOURCES} main.cpp)

# Executable target
add_executable(LunarSim_v0_2 ${SOURCES} ${HEADERS})

# Link FFTW3
link_directories(external/fftw-3.3.5-dll64)
target_link_libraries(LunarSim_v0_2 PUBLIC libfftw3-3)

# Link yaml-cpp
include(FetchContent)
set(GIT_EXECUTABLE "C:/Program Files/Git/cmd/git.exe")

FetchContent_Declare(
    yaml-cpp
    GIT_REPOSITORY https://github.com/jbeder/yaml-cpp.git
    GIT_TAG master
)
FetchContent_MakeAvailable(yaml-cpp)

# Set the library directory
target_link_libraries(LunarSim_v0_2 PUBLIC yaml-cpp::yaml-cpp)

target_link_libraries(LunarSim_v0_2 PUBLIC libfftw3-3)

