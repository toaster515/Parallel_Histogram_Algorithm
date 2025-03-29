
#include <cstdlib>
#include <ctime>

#include <random>
#include <string>
#include <sstream>
#include <vector>

#include <gtest/gtest.h>
#include "gtest-mpi-listener.hpp"
#include "gtest-timeout.hpp"
#include "mpi.h"


#include "databasics.hpp"
#include "datageneration.hpp"
#include "basic_defs.h"
#include "solution.hpp"

#include "test_helpers.hpp"

#define REBALANCE_TIMEOUT_S 30
#define SPLITTERS_TIMEOUT_S 30
#define MOVE_TIMEOUT_S 30
#define SORT_TIMEOUT_S 5



TEST(rebalance, correct){
	const int num_iterations = 1;
	// Get number of processes
	int rank, nProcs;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nProcs);


	for(int i = 0;i<num_iterations;i++){
		MPI_Barrier(MPI_COMM_WORLD);
		dist_sort_t *my_data;
		dist_sort_size_t my_data_size;
		dist_sort_size_t global_data_size;

		dist_sort_t *rebalanced_data = nullptr;
		dist_sort_size_t rebalanced_data_size = 0;

		test_rebalance_datagen(&my_data, &my_data_size, &global_data_size, 0);

		uint64_t confirm_mask = randuint64();
		uint64_t confirm_modulus = randuint64() % 0xffffffff;



		uint64_t gen_checksum = distributed_checksum(my_data, my_data_size, confirm_mask, confirm_modulus, MPI_COMM_WORLD);

		TEST_TIMEOUT_BEGIN
		rebalance(my_data, my_data_size, &rebalanced_data, &rebalanced_data_size);
		TEST_TIMEOUT_FAIL_END(REBALANCE_TIMEOUT_S*1000)

		ASSERT_NE(rebalanced_data, (dist_sort_t*)nullptr) << "No data was malloc()ed/returned.";
		ASSERT_GT(rebalanced_data_size, 0) << "No data was returned.";

		if(nullptr != rebalanced_data){
			uint64_t rebalanced_checksum = distributed_checksum(rebalanced_data, rebalanced_data_size, confirm_mask, confirm_modulus, MPI_COMM_WORLD);
			EXPECT_EQ(gen_checksum,rebalanced_checksum) << "Mismatched checksums mean data was lost, gained, or corrupted. (Order invariant.)";
		}

		free(my_data);
		if(nullptr != rebalanced_data){
			free(rebalanced_data);
		}
	}
}

TEST(rebalance, balanced){
	const int num_iterations = 1;
	const double allowable_size_deviation = 0.01;
	// Get number of processes
	int rank, nProcs;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nProcs);


	for(int i = 0;i<num_iterations;i++){
		MPI_Barrier(MPI_COMM_WORLD);
		dist_sort_t *my_data;
		dist_sort_size_t my_data_size;
		dist_sort_size_t global_data_size;

		dist_sort_t *rebalanced_data = nullptr;
		dist_sort_size_t rebalanced_data_size = 0;

		test_rebalance_datagen(&my_data, &my_data_size, &global_data_size, 1);

		uint64_t confirm_mask = randuint64();
		uint64_t confirm_modulus = randuint64() % 0xffffffff;

		uint64_t gen_checksum = distributed_checksum(my_data, my_data_size, confirm_mask, confirm_modulus, MPI_COMM_WORLD);

		TEST_TIMEOUT_BEGIN
		rebalance(my_data, my_data_size, &rebalanced_data, &rebalanced_data_size);
		TEST_TIMEOUT_FAIL_END(REBALANCE_TIMEOUT_S*1000)

		ASSERT_NE(rebalanced_data, (dist_sort_t*)nullptr) << "No data was malloc()ed/returned.";
		ASSERT_GT(rebalanced_data_size, 0) << "No data was returned.";

		if(nullptr != rebalanced_data){
			uint64_t rebalanced_checksum = distributed_checksum(rebalanced_data, rebalanced_data_size, confirm_mask, confirm_modulus, MPI_COMM_WORLD);
			//We don't check this because we want to have the option at grading time to separate correctness and balance.
			//EXPECT_EQ(gen_checksum,rebalanced_checksum) << "Mismatched checksums mean data was lost, gained, or corrupted. (Order invariant.)";
		}

		//Test that data was rebalanced within the required thresholds.
		long double average_size = (long double)(global_data_size / nProcs);

		long double my_rebalanced_size_deviation = fabs((long double)(rebalanced_data_size - average_size) / average_size);

		EXPECT_LE(my_rebalanced_size_deviation, allowable_size_deviation) << "This process (" << rank << ") has a deviation of " << my_rebalanced_size_deviation << " .";

		free(my_data);
		if(nullptr != rebalanced_data){
			free(rebalanced_data);
		}
	}// End of for:num_iterations
}


TEST(findSplitters, balanced){

	const int num_iterations = 1;
	const double allowable_bin_deviation = 0.01;
	// Get number of processes
	int rank, nProcs;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nProcs);


	dist_sort_t *test_data;
	dist_sort_size_t test_data_size;
	dist_sort_size_t global_size;

	dist_sort_t *test_splitters;
	dist_sort_size_t *test_counts;
	int num_test_splitters;


	{//may become a loop later, ensure scoping is good now.
		num_test_splitters = nProcs;
		test_splitters = (dist_sort_t*)malloc(num_test_splitters*sizeof(dist_sort_t));
		test_counts = (dist_sort_size_t*)malloc(num_test_splitters*sizeof(dist_sort_size_t));

		test_splitters_datagen(&test_data, &test_data_size, &global_size, 0 );//Uniform

		//findSpliters expects sorted data, at least according to our API.
		std::sort(test_data,test_data+test_data_size);

		TEST_TIMEOUT_BEGIN
		findSplitters(test_data, test_data_size, test_splitters, test_counts, num_test_splitters);
		TEST_TIMEOUT_FAIL_END(SPLITTERS_TIMEOUT_S*1000)

		long double avg_size = (long double)global_size / (long double)nProcs;
		if(0 == rank){
			for(int i = 0;i<num_test_splitters;++i){
				long double bin_deviation = fabs(test_counts[i] - avg_size)/avg_size;
				EXPECT_LE(bin_deviation,allowable_bin_deviation) << "Problem with quantiles in bin "<<i<<".";
				if(bin_deviation > allowable_bin_deviation) break;
			}
		}

		free(test_data);
		free(test_splitters);
		free(test_counts);
	}
}

TEST(moveData, correct){
	const int num_iterations = 1;
	const double allowable_bin_deviation = 0.01;
	// Get number of processes
	int rank, nProcs;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nProcs);


	dist_sort_t *test_data;
	dist_sort_size_t test_data_size;
	dist_sort_size_t global_size;

	dist_sort_t *moved_data = nullptr;
	dist_sort_size_t moved_data_size = 0;
	dist_sort_size_t post_move_global_size;

	dist_sort_t *test_splitters;
	dist_sort_size_t *test_counts;
	int num_test_splitters;


	{//may become a loop later, ensure scoping is good now.
		num_test_splitters = nProcs;

		test_splitters = (dist_sort_t*)malloc(num_test_splitters*sizeof(dist_sort_t));
		test_counts = (dist_sort_size_t*)malloc(num_test_splitters*sizeof(dist_sort_size_t));

		dist_sort_t *test_regions = (dist_sort_t*)malloc(num_test_splitters*sizeof(dist_sort_t));
		dist_sort_size_t *test_weights = (dist_sort_size_t*)malloc(num_test_splitters*sizeof(dist_sort_size_t));

		if(0 == rank){
			for(int i = 0;i<num_test_splitters;++i){
				test_regions[i] = randuint64();
				test_weights[i] = randuint64() % 1234 + 10;
			}
			std::sort(test_regions, test_regions+num_test_splitters);
			std::sort(test_weights, test_weights+num_test_splitters);
		}
		MPI_Bcast(test_regions,num_test_splitters,MPI_UINT64_T, 0, MPI_COMM_WORLD);
		MPI_Bcast(test_weights,num_test_splitters,MPI_UINT64_T, 0, MPI_COMM_WORLD);

		//test_splitters_datagen(&test_data, &test_data_size, &global_size, 1 );//Known histogram

		splitter_correctness_datagen(&test_data, &test_data_size, &global_size,
			test_regions, test_weights, num_test_splitters,
			0 );

		//findSpliters expects sorted data, at least according to our API.
		std::sort(test_data,test_data+test_data_size);

		//student may or may not have actually used the counts in moveData, need to be ready if they did.
		for(int i = 0;i<num_test_splitters;++i){ test_weights[i] *= nProcs;}

		uint64_t confirm_mask = randuint64();
		uint64_t confirm_modulus = randuint64() % 0xffffffff;
		uint64_t gen_checksum = distributed_checksum(test_data, test_data_size, confirm_mask, confirm_modulus, MPI_COMM_WORLD);

		TEST_TIMEOUT_BEGIN
		moveData(test_data, test_data_size, &moved_data, &moved_data_size,test_regions,test_weights, nProcs);
		TEST_TIMEOUT_FAIL_END(MOVE_TIMEOUT_S*1000)

		ASSERT_NE(moved_data,(dist_sort_t*)nullptr) << "No output memory was malloc()ed";
		ASSERT_GT(moved_data_size,0) << "There should probably be more than 0 data on every rank";

		if(nullptr != moved_data){
			uint64_t moved_checksum = distributed_checksum(moved_data, moved_data_size, confirm_mask, confirm_modulus, MPI_COMM_WORLD);
			if(0 == rank){
				EXPECT_EQ(gen_checksum, moved_checksum) << "Data checksum wrong after moveData.";
			}

			//sort after move
			std::sort(moved_data,moved_data+moved_data_size);
			if(rank > 0){
				EXPECT_LT(test_regions[rank-1], moved_data[0]) << "This rank (" << rank << ") holds data too small for its assigned region.";
			}
			if(rank < (nProcs -1)){
				//Last splitter should be equal to the maximum in the data, but if it's not, that's a problem in findSplitters, not in moveData.
				//Someone has to take all the data, so it might go to the last processor even if the last splitter is wrong.
				//Rememeber that the last splitter should just contain the max anyway.
				EXPECT_LE(moved_data[moved_data_size-1], test_regions[rank]) << "This rank (" << rank << ") holds data too large for its assigned region.";
			}
		}

		free(test_data);
		free(test_splitters);
		free(test_counts);
		free(test_regions);
		free(test_weights);
		if(nullptr != moved_data){
			free(moved_data);
		}
	}
}

TEST(sort, correct){
	const dist_sort_size_t sort_test_size = 1234;
	// Get number of processes
	int rank, nProcs;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nProcs);

	dist_sort_t *thedata;// = (dist_sort_t*)malloc(sort_test_size*sizeof(dist_sort_t));

	if(0 == rank){//One test is sufficient.
		generateData(0, &thedata, sort_test_size, 0, INT64_MAX);

		TEST_TIMEOUT_BEGIN
		sort(thedata, sort_test_size);
		TEST_TIMEOUT_FAIL_END(SORT_TIMEOUT_S*1000)

		for(int i = 1;i<sort_test_size;++i){
			ASSERT_LE(thedata[i-1] , thedata[i]) << "Data not monotonically increasing.";
			//if(thedata[i-1] > thedata[i]) break;
		}

		//This happens even if the test is failed.
		free(thedata);
	}
}


int main(int argc, char** argv) {
	// Filter out Google Test arguments
	::testing::InitGoogleTest(&argc, argv);

	// Initialize MPI
	int mpi_thread_level;
	MPI_Init_thread(&argc, &argv,MPI_THREAD_MULTIPLE,&mpi_thread_level);
	if(mpi_thread_level < MPI_THREAD_MULTIPLE){
		std::cerr << "We require MPI_THREAD_MULTIPLE in case students used threads." << std::endl;
		exit(1);
	}

	int world_rank;
	int world_size;

	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);

	//srand(time(nullptr)+world_rank);
	srand(1337+world_rank);

	// Add object that will finalize MPI on exit; Google Test owns this pointer
	::testing::AddGlobalTestEnvironment(new MPIEnvironment);

	// Get the event listener list.
	::testing::TestEventListeners& listeners =
		::testing::UnitTest::GetInstance()->listeners();

	//only get json/xml output from rank 0
	if(0 != world_rank){
		// Remove default listener
		delete listeners.Release(listeners.default_result_printer());
		// Remove default file output
		delete listeners.Release(listeners.default_xml_generator());
	}
	listeners.Append(new MPICollectTestResults);

	/*
	if(4 > world_size){
		std::cerr << " Can't run these tests on fewer than 4 ranks " << std::endl;
		exit(1);
	}
	*/

	// Run tests, then clean up and exit
	int foo = RUN_ALL_TESTS();

	return foo;
}
