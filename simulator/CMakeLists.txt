# Minimum required version of CMake
cmake_minimum_required(VERSION 3.10)

# Project name and version
project(fds VERSION 1.0)

# Set the C++ standard
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_BUILD_TYPE Debug)

# Specify the output directory for compiled binaries
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Add the include directory to the compiler search path
include_directories(${PROJECT_SOURCE_DIR}/include)

# Collect all .cpp files in the src directory
file(GLOB SOURCES "${PROJECT_SOURCE_DIR}/src/*.cpp")

find_package(Eigen3 3.3 REQUIRED NO_MODULE)

# Define the executable target
add_executable(${PROJECT_NAME} ${SOURCES})
target_link_libraries(${PROJECT_NAME} Eigen3::Eigen)

# Copy resource files to the output directory after build
set(RESOURCE_FILES
    "${PROJECT_SOURCE_DIR}/resources/tle.txt"
    "${PROJECT_SOURCE_DIR}/resources/config.json")

foreach(RESOURCE_FILE ${RESOURCE_FILES})
    add_custom_command(
        TARGET ${PROJECT_NAME} POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy ${RESOURCE_FILE} $<TARGET_FILE_DIR:${PROJECT_NAME}>
    )
endforeach()

target_compile_definitions(${PROJECT_NAME} PRIVATE BUILD_OUTPUT_PATH="${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")

# Set up testing
add_executable(test_OrbitalMechanics tests/test_OrbitalMechanics.cpp src/OrbitalMechanics.cpp)
enable_testing()
add_test(NAME OrbitalMechanicsTest COMMAND test_OrbitalMechanics)
