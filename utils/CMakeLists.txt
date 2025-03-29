cmake_minimum_required (VERSION 3.12)
project(project-utils VERSION 1.0.0 LANGUAGES C CXX)

find_package(MPI REQUIRED)

add_library(utils "src/databasics.cpp" "include/databasics.hpp")
target_include_directories(utils PUBLIC "include")
target_link_libraries(utils PUBLIC MPI::MPI_CXX)

add_library(orderedio "src/orderedio.c" "include/orderedio.h")
target_include_directories(orderedio PUBLIC "include")
target_link_libraries(orderedio PUBLIC MPI::MPI_CXX)


add_library(datageneration OBJECT "src/datageneration.cpp"
    "include/datageneration.hpp"
    "../src/basic_defs.h"
)
target_include_directories(datageneration PUBLIC "../src")
target_link_libraries(datageneration utils)


####### TESTS ######
find_package(GTest)
if(GTEST_FOUND)
include(GoogleTest)

add_executable(utiltests "tests/utiltests.cpp" "tests/gtest-mpi-listener.hpp")
target_include_directories(utiltests PUBLIC "src" "include")
target_link_libraries(utiltests utils)
target_link_libraries(utiltests GTest::GTest GTest::Main MPI::MPI_CXX)

endif(GTEST_FOUND)
