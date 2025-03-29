cmake_minimum_required (VERSION 3.12)
project(project VERSION 1.0.0 LANGUAGES C CXX)

# Location of additional cmake modules
#
#add_definitions(-DDATA_SIZE_MULTIPLE_OF_RANKS)#TODO: Remove this and the dependent code (sp19 leftover)

set(CMAKE_MODULE_PATH
    ${CMAKE_MODULE_PATH}
    ${PROJECT_SOURCE_DIR}/cmake
    )

include("cs484common")

##### Project debugging #####
#set(CMAKE_CXX_CLANG_TIDY "clang-tidy;-checks=*") #use a linter for the entire project.


find_package(MPI REQUIRED)
find_package(Threads REQUIRED)

add_subdirectory("utils")

####### SOLUTION CODE ########
#set(CMAKE_CXX_CLANG_TIDY "clang-tidy;-checks=mpi*") #Use a linter for solution code.

add_library(student_solution OBJECT "src/solution.cpp"
	"src/solution.hpp"
	"src/basic_defs.h"
)
target_link_libraries(student_solution utils) #The utility code is available to students
target_link_libraries(student_solution MPI::MPI_CXX) #Students will use MPI
target_link_libraries(student_solution Threads::Threads) #Students might use pthreads

unset(CMAKE_CXX_CLANG_TIDY)#Do not use the linter for provided code.
##### END SOLUTION ######


##### DRIVER / BENCHMARK CODE ########
add_executable(sorter
		"src/main_sort.cpp"
		"src/basic_defs.h"
)
target_include_directories(sorter PUBLIC "${PROJECT_SOURCE_DIR}/src")
target_link_libraries(sorter student_solution)
target_link_libraries(sorter datageneration)
target_link_libraries(sorter utils orderedio)
target_link_libraries(sorter m MPI::MPI_CXX)

###### TESTS #######
add_subdirectory("tests")
