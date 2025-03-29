#include <cstdlib>
#include <ctime>

#include <random>
#include <string>
#include <sstream>
#include <vector>

#include "mpi.h"

#include "databasics.hpp"
#include "datageneration.hpp"
#include "basic_defs.h"


#include "test_helpers.hpp"


dist_sort_size_t test_rebalance_datagen(dist_sort_t **thedata, dist_sort_size_t *mycount, dist_sort_size_t *globalcount, int mode ){
	// Get number of processes
	int rank, num_processors;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &num_processors);

	dist_sort_size_t local_N;
	dist_sort_size_t global_N;
	//Round up to the nearest multiple of the number of processors.
	//We promised that. Should remove later.
	//TODO(bjlunt2): Create an enum to specify the mode more clearly.
	switch(mode){
		case 0:
			local_N = data_randsize(0, 30000, 20000);
			break;
		case 1:
			local_N = 50000 / (1 + rank);
			break;
		default:
			local_N = data_randsize(0, 50000, 10000);
	}
	
	#ifdef DATA_SIZE_MULTIPLE_OF_RANKS
	/*
	Ensure that whatever the datasize, it is a multiple of the number of ranks.
	The easy way is to make it a multiple of the number of ranks _on each rank_.
	*/
	local_N += num_processors - (local_N % num_processors);
	#endif

	MPI_Allreduce(&local_N, &global_N, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);

	generateData(0, thedata, local_N, 0, INT64_MAX);

	//TODO: Actually allocate data.
	*mycount = local_N;
	*globalcount = global_N;
}

dist_sort_size_t test_splitters_datagen(dist_sort_t **thedata, dist_sort_size_t *mycount, dist_sort_size_t *globalcount, int mode ){
	// Get number of processes
	int rank, num_processors;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &num_processors);

	dist_sort_size_t local_N;
	dist_sort_size_t global_N;
	//Round up to the nearest multiple of the number of processors.
	//We promised that. Should remove later.

	local_N = 1000 + data_randsize(0, 2000, 5000);
	//local_N += num_processors - (local_N % num_processors);

	if(0 == mode){//Uniform Distribution over uint64
		dist_sort_t *tmp = (dist_sort_t*)malloc(local_N*sizeof(dist_sort_t));
		int i;
		for(i = 0;i< local_N;i++){
			tmp[i] = randuint64();
		}
		*thedata = tmp;
	}else if(1 == mode){
		//Sawtooth.

	}else{
		generateData(0, thedata, local_N, 0, INT64_MAX);
	}


	MPI_Allreduce(&local_N, &global_N, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);

	//TODO: Actually allocate data.
	*mycount = local_N;
	*globalcount = global_N;
}

dist_sort_size_t splitter_correctness_datagen(dist_sort_t **thedata, dist_sort_size_t *mycount, dist_sort_size_t *globalcount,
	const dist_sort_t *regions, const dist_sort_size_t *weights, const int num_regions,
	const int mode ){
	// Get number of processes
	int rank, num_processors;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &num_processors);

	dist_sort_size_t local_N;
	dist_sort_size_t global_N;
	//Round up to the nearest multiple of the number of processors.
	//We promised that. Should remove later.

	local_N = 1000 + data_randsize(0, 2000, 5000);
	//local_N += num_processors - (local_N % num_processors);
	if(0 == mode){
		local_N = 0; for(int i = 0;i<num_regions;++i){ local_N += weights[i];}
		dist_sort_t *tmp = (dist_sort_t*)malloc(local_N*sizeof(dist_sort_t));

		dist_sort_t region_minimum = 0;
		dist_sort_t *tmp_region = tmp;

		for(int i = 0;i<num_regions;i++){
			dist_sort_t region_size = regions[i] - region_minimum;
			for(int j = 0;j<(weights[i]-1);j++){
				//random values guaranteed to be within region.
				tmp_region[j] = randuint64() % region_size + region_minimum;
			}
			//To ensure that the splitters come out just right
			tmp_region[weights[i]-1] = regions[i];

			region_minimum = regions[i];
			tmp_region = tmp_region + weights[i];
		}


		*thedata = tmp;
	}else if(1 == mode){
		unsigned seed1 = rand(); //+std::chrono::system_clock::now().time_since_epoch().count();
		std::mt19937 g1 (seed1);

		//std::default_random_engine generator;
		std::vector<double> intervals(regions,regions+num_regions);
		std::vector<double> use_weights(weights,weights+num_regions);
		std::piecewise_constant_distribution<double>distribution(intervals.begin(),intervals.end(),use_weights.begin());

		dist_sort_t *tmp = (dist_sort_t*)malloc(local_N*sizeof(dist_sort_t));
		int i;
		for(i = 0;i< local_N;i++){
			tmp[i] = (dist_sort_t)distribution(g1);
		}
		*thedata = tmp;
	}else{
		generateData(0, thedata, local_N, 0, INT64_MAX);
	}


	MPI_Allreduce(&local_N, &global_N, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);

	//TODO: Actually allocate data.
	*mycount = local_N;
	*globalcount = global_N;
}
