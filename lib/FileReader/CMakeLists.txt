cmake_minimum_required(VERSION 3.5)
set(project_name FileReader)
project(${project_name} LANGUAGES CXX)

set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

SET(MAJOR_VERSION 0)
SET(MINOR_VERSION 0)
SET(MICRO_VERSION 1)
SET(VERSION "${MAJOR_VERSION}.${MINOR_VERSION}.${MICRO_VERSION}")

add_definitions(-DVERSION="${VERSION}")
add_definitions(-DSOFTWARE_NAME="${project_name}")
add_definitions(-DCORPORATION_NAME="Damavan-imaging")
add_definitions(-DNO_YAML)

option(FILEREADER_BUILD_TESTS "Enable testing" OFF)

find_package(Qt5 COMPONENTS Core Gui Test REQUIRED)
include_directories(./ ./lib/CommonQuadratorHeader/510)

if(FILEREADER_BUILD_TESTS)
	add_compile_options(-fprofile-arcs -ftest-coverage -fPIC -O0)
	enable_testing()
	add_subdirectory(tests)

endif()

add_library(${project_name} STATIC filereader.cpp)
target_link_libraries(${project_name} Qt5::Core Qt5::Gui)

