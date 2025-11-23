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

# Collect all .cpp files in the src directory, excluding job_processor_main.cpp
file(GLOB ALL_SOURCES "${PROJECT_SOURCE_DIR}/src/*.cpp")
list(FILTER ALL_SOURCES EXCLUDE REGEX ".*job_processor_main\\.cpp$")
set(SOURCES ${ALL_SOURCES})

find_package(Eigen3 3.3 REQUIRED NO_MODULE)
include_directories(${EIGEN3_INCLUDE_DIR})

find_package(Boost 1.74.0 COMPONENTS) 
include_directories(${Boost_INCLUDE_DIRS})

# Find PostgreSQL
find_package(PostgreSQL REQUIRED)

# Find pqxx using pkg-config
find_package(PkgConfig REQUIRED)
pkg_check_modules(PQXX REQUIRED libpqxx)

# Find hiredis
pkg_check_modules(HIREDIS REQUIRED hiredis)

# Define the main simulator executable
add_executable(${PROJECT_NAME} ${SOURCES})
target_link_libraries(${PROJECT_NAME} 
    Eigen3::Eigen 
    ${Boost_LIBRARIES}
    PostgreSQL::PostgreSQL
    ${PQXX_LIBRARIES}
    ${HIREDIS_LIBRARIES}
)
target_include_directories(${PROJECT_NAME} PRIVATE 
    ${PQXX_INCLUDE_DIRS}
    ${HIREDIS_INCLUDE_DIRS}
)
target_compile_options(${PROJECT_NAME} PRIVATE ${PQXX_CFLAGS_OTHER})

# Define the job processor executable
add_executable(job_processor 
    src/job_processor_main.cpp
    src/JobProcessor.cpp
    src/RedisClient.cpp
    src/Database.cpp
)
target_link_libraries(job_processor
    Eigen3::Eigen
    ${Boost_LIBRARIES}
    PostgreSQL::PostgreSQL
    ${PQXX_LIBRARIES}
    ${HIREDIS_LIBRARIES}
)
target_include_directories(job_processor PRIVATE 
    ${PROJECT_SOURCE_DIR}/include
    ${PQXX_INCLUDE_DIRS}
    ${HIREDIS_INCLUDE_DIRS}
)
target_compile_options(job_processor PRIVATE ${PQXX_CFLAGS_OTHER})

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
