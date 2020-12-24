include_directories(..)
include_directories(.)

add_executable(FileReader-test
    tst_filereader.cpp
#    ../src/Configuration.cpp
)

add_definitions(-DSOURCES_DIR=\"${CMAKE_CURRENT_SOURCE_DIR}\")

set_target_properties(FileReader-test PROPERTIES
    CXX_STANDARD 11
    CXX_STANDARD_REQUIRED ON
)

target_link_libraries(FileReader-test FileReader -lgcov Qt5::Test)

add_test(FileReader-test ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/FileReader-test)

