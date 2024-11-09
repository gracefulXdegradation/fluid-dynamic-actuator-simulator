# Set the minimum version of CMake that can be used
# To find the cmake version run
# $ cmake --version
cmake_minimum_required(VERSION 3.10)

# Set the project name
project (flight_dynamics_simulator VERSION 1.0)

# Create a sources variable with a link to all cpp files to compile
set(SOURCES
    src/Hello.cpp
    src/main.cpp
)

# Add an executable with the above sources
add_executable(${PROJECT_NAME} ${SOURCES})

# Specify the resource directory and static file path
set(RESOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/resources")
set(RESOURCE_FILES "${RESOURCE_DIR}/tle.txt")

# Copy static files to the build directory
foreach(RESOURCE_FILE ${RESOURCE_FILES})
    add_custom_command(
        TARGET ${PROJECT_NAME} POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy ${RESOURCE_FILE} $<TARGET_FILE_DIR:${PROJECT_NAME}>
    )
endforeach()

target_compile_definitions(${PROJECT_NAME} PRIVATE TLE_PATH="${CMAKE_CURRENT_BINARY_DIR}/tle.txt")

# Set the directories that should be included in the build command for this target
# when running g++ these will be included as -I/directory/path/
target_include_directories(${PROJECT_NAME}
    PRIVATE 
        ${PROJECT_SOURCE_DIR}/include
)