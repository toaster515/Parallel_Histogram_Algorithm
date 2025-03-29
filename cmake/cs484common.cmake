
include("PreventInTreeBuilds") #Very important to keep students from building in wrong place


set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CXX_FLAGS "-O3 -Wfatal-errors")
set(CMAKE_CXX_FLAGS_DEBUG "-g -Wall -Wextra -Wfatal-errors")

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_DATA_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/data)

find_package(Threads REQUIRED)

#we need an openmp compiler
find_package(OpenMP)
if(OPENMP_FOUND)
	message("Found OpenMP, building parallel part.")

	#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
	#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
	#set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")

else(OPENMP_FOUND)
		message("ERROR: OpenMP not found.")
endif(OPENMP_FOUND)



#we're usually going to need mpi
find_package(MPI)
if(MPI_CXX_FOUND)
	message("Found MPI")
else(MPI_CXX_FOUND)
	message(FATAL_ERROR "Could not find MPI")
endif(MPI_CXX_FOUND)

#some assignments require Charm++ or aMPI.

#find_package(Charm)
#if(CHARM_FOUND)
#	message("Charm found: ${CHARM_COMPILER}")
#else(CHARM_FOUND)
#	message(FATAL_ERROR "Could not find Charm++")
#endif(CHARM_FOUND)
#Tests and grading

enable_testing()
find_package(GTest REQUIRED)#comes with CMake
if(GTEST_FOUND)
	#include(GoogleTest)

else(GTEST_FOUND)
	message(FATAL_ERROR "DANGER: Could not find GoogleTest")
endif(GTEST_FOUND)


#Benchmarking
find_package(GBench)
if(GBENCH_FOUND)

else(GBENCH_FOUND)
	message(FATAL_ERROR "DANGER: Could not find GoogleBenchmark")
endif(GBENCH_FOUND)
