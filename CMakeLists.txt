######################################################################
# @author      : lopezmar (lopezmar@$HOSTNAME)
# @modified by : o.allegrini (o.allegrini@ipnl.in2p3.fr)
# @file        : CMakeLists
# @created     : jeudi oct. 15, 2020 16:45:17 CEST
######################################################################

cmake_minimum_required(VERSION 3.5)
project(Test)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


#####################################################################
#
#                           ROOT
#
#####################################################################

#---Locate the ROOT package and defines a number of variables (e.g. ROOT_INCLUDE_DIRS)

find_package(ROOT REQUIRED COMPONENTS Spectrum)
include(${ROOT_USE_FILE})


#####################################################################
#
#                           DAMAVAN
#
#####################################################################

add_subdirectory(lib/FileReader)
include_directories(/Users/allegrini/Documents/11.CameraTemporal/TestFileReader/lib/FileReader 
/Users/allegrini/Documents/11.CameraTemporal/TestFileReader/lib/FileReader/lib/CommonQuadratorHeader/510
/Users/allegrini/Documents/11.CameraTemporal/TestFileReader/functions
${ROOT_INCLUDE_DIRS})


add_definitions(-DNO_YAML ${ROOT_CXX_FLAGS})

add_executable(Calibration main_calibration.cpp)
add_executable(Analysis main.cpp)

target_link_libraries(Calibration FileReader ${ROOT_LIBRARIES})
target_link_libraries(Analysis FileReader ${ROOT_LIBRARIES})
