#ifndef _TEST_HELPERS_HPP_
#define _TEST_HELPERS_HPP_

#include "basic_defs.h"

dist_sort_size_t test_rebalance_datagen(dist_sort_t **thedata, dist_sort_size_t *mycount, dist_sort_size_t *globalcount, int mode );

dist_sort_size_t test_splitters_datagen(dist_sort_t **thedata, dist_sort_size_t *mycount, dist_sort_size_t *globalcount, int mode );

dist_sort_size_t splitter_correctness_datagen(dist_sort_t **thedata, dist_sort_size_t *mycount, dist_sort_size_t *globalcount,
	const dist_sort_t *regions, const dist_sort_size_t *weights, const int num_regions,
	const int mode );

#endif // _TEST_HELPERS_HPP_
