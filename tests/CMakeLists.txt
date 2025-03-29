cmake_minimum_required (VERSION 3.12)
project(project-tests VERSION 1.0.0 LANGUAGES C CXX)


#Tests and grading
find_package(GTest REQUIRED)

include(CTest)
include(GoogleTest)

######## Testing ##########
add_library(test_helpers OBJECT "test_helpers.cpp" "./include/test_helpers.hpp")
target_include_directories(test_helpers PUBLIC "./include")
target_link_libraries(test_helpers datageneration)

add_executable(main_tests  "main_tests.cpp"
						"./include/gtest-mpi-listener.hpp"
)
target_include_directories(main_tests PUBLIC "../src" "../utils/include" "./include")
target_link_libraries(main_tests datageneration test_helpers utils student_solution)
target_link_libraries(main_tests m GTest::GTest GTest::Main MPI::MPI_CXX Threads::Threads)

gtest_discover_tests(main_tests)




if(DEFINED GRADING)
	#Don't compile student tests at grading time.
	message(STATUS "BUILDING FOR GRADING: ${GRADING}.")
	add_definitions(-DGRADING="${GRADING}")

else(DEFINED GRADING)
	#Student tests only get compiled when _not_ grading.
	#add_definitions(-DGRADING="${GRADING}")

	######## Student written tests ##########
	add_executable(student_tests "student_tests.cpp"
							"./include/gtest-mpi-listener.hpp")
	include_directories(student_tests "../src" "../utils/include" "./include")
	target_link_libraries(student_tests datageneration test_helpers utils student_solution)
	target_link_libraries(student_tests m GTest::GTest GTest::Main MPI::MPI_CXX Threads::Threads)
endif(DEFINED GRADING)
