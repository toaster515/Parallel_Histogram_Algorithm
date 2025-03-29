#include <cmath>
#include <algorithm>
#include <iostream>
#include <utility>
#include <queue>
#include <vector>

#include <cassert>
#include <cstring>

#include "basic_defs.h"
#include "databasics.hpp"
#include "solution.hpp"

//***************************************************************************************************************************************************************************************************************
//__________________________________________________________Rebalance_____________________________________________________________________
//***************************************************************************************************************************************************************************************************************
void rebalance(const dist_sort_t *data, const dist_sort_size_t myDataCount, dist_sort_t **rebalancedData, dist_sort_size_t *rCount) {

//****************************************************************************************************************************************
	MPI_Datatype swap_instr_t;
    MPI_Type_contiguous(3, MPI_INT, &swap_instr_t);
    MPI_Type_commit(&swap_instr_t);
	dist_sort_size_t global_count;
	MPI_Allreduce(&myDataCount, &global_count, 1, MPI_TYPE_DIST_SORT_SIZE_T, MPI_SUM, MPI_COMM_WORLD);
	int g_count = (int)global_count;
	int nProcs;
	MPI_Comm_size(MPI_COMM_WORLD, &nProcs);
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Barrier(MPI_COMM_WORLD);

//****************************************************************************************************************************************
	//Open a window, all ranks put their data counts in an 0's array indexed to their ranks, essentially a gather
	//Rank 0 sends instructions for distribution
	int *counts=(int*)malloc(sizeof(int)*nProcs);
	int myCount = (int)myDataCount;
	MPI_Win window;	
	MPI_Win_create(counts, nProcs*sizeof(int), sizeof(int), MPI_INFO_NULL, MPI_COMM_WORLD, &window);
	MPI_Win_fence(0, window);
		MPI_Put(&myCount, 1, MPI_INT, 0, rank, 1, MPI_INT, window);
	MPI_Win_fence(0, window);
	MPI_Win_free(&window);
	
	MPI_Barrier(MPI_COMM_WORLD);
//****************************************************************************************************************************************
	//Rank 0 uses counts to determine optimal counts
	//creates window for other ranks to GET their optimal count
	int *bal_counts=(int*)malloc(sizeof(int)*nProcs);
	if (rank==0){		
		int *bal_diff=(int*)malloc(sizeof(int)*nProcs);
		int bal = g_count / nProcs;
		int rem = g_count % nProcs;
		//Determine Optimal Counts
		for (int i=0; i<nProcs; i++){
			bal_counts[i]=bal;
			if (rem>0){
				bal_counts[i]+=1;
				rem-=1;
			}
		}
		//Determine Node counts and difference from optimal
		for (int i =0; i<nProcs; i++){
			bal_diff[i]=counts[i]-bal_counts[i];
		}
	}
	MPI_Win_create(bal_counts, nProcs*sizeof(int), sizeof(int), MPI_INFO_NULL, MPI_COMM_WORLD, &window);
	MPI_Win_fence(0, window);
		if (rank!=0){MPI_Get(&bal_counts[0], nProcs, MPI_INT, 0, 0, nProcs, MPI_INT, window);}
	MPI_Win_fence(0, window);
	MPI_Win_free(&window);
	MPI_Barrier(MPI_COMM_WORLD);
	int myDif = myCount - bal_counts[rank];

	//Ranks allocate appropriate space for their balanced data arrays
	*rebalancedData = (dist_sort_t*)malloc(sizeof(dist_sort_t)*bal_counts[rank]);
	*rCount = (dist_sort_size_t)bal_counts[rank];
	//*rCount = (dist_sort_size_t*)bal_counts[rank];
	MPI_Barrier(MPI_COMM_WORLD);
//****************************************************************************************************************************************
	//Open Window
	//ranks fill Rank0 window with their differences from optimal count
	int *bal_diff=(int*)malloc(sizeof(int)*nProcs);
	MPI_Win_create(bal_diff, nProcs*sizeof(int), sizeof(int), MPI_INFO_NULL, MPI_COMM_WORLD, &window);
	MPI_Win_fence(0, window);
		MPI_Put(&myDif, 1, MPI_INT, 0, rank, 1, MPI_INT, window);
	MPI_Win_fence(0, window);
	MPI_Win_free(&window);
	MPI_Barrier(MPI_COMM_WORLD);
//****************************************************************************************************************************************
	//swap data
	std::queue<swap_instr> instructions;
	if (rank==0){//I am the captain
		//Create and Fill Queues
		std::queue<std::pair<int, int>> surplus_queue;
		std::queue<std::pair<int, int>> deficit_queue;
		std::pair<int, int> tmp, sq, dq;

		for (int i=0; i<nProcs; i++){
			if (bal_diff[i]>0){
				tmp.first = i;
				tmp.second = bal_diff[i];
				surplus_queue.push(tmp);
			}else if(bal_diff[i] < 0){
				tmp.first = i;
				tmp.second = bal_diff[i];
				deficit_queue.push(tmp);
			}
		}
		
		//Master uses queues to send messages: 
		//command instructions for nodes to send their data to each other
		MPI_Request req;
		sq = surplus_queue.front();
		dq = deficit_queue.front();
		while(true){
			if (surplus_queue.size()==0){
				break;
			}
			if (deficit_queue.size()==0){
				break;
			}
			if (sq.second > abs(dq.second)){//surplus greater than deficit
				swap_instr snd_msg={0, dq.first, abs(dq.second)};
				if (sq.first==0){//I'm sending
					instructions.push(snd_msg);
				}else{
					MPI_Isend(&snd_msg, 1, swap_instr_t, sq.first, 0, MPI_COMM_WORLD, &req);
				}
				swap_instr rec_msg={1, sq.first, abs(dq.second)};
				if (dq.first==0){//I'm receiving
					instructions.push(rec_msg);
				}else{
					MPI_Isend(&rec_msg, 1, swap_instr_t, dq.first, 0, MPI_COMM_WORLD, &req);
				}
				bal_diff[sq.first]-=abs(dq.second);
				bal_diff[dq.first]+=abs(dq.second);
				sq.second += dq.second;
				deficit_queue.pop();
				dq=deficit_queue.front();

			}else if (sq.second < abs(dq.second)){//deficit greater than surplus
				swap_instr snd_msg={0, dq.first, sq.second};
				if (sq.first==0){//I'm sending
					instructions.push(snd_msg);
				}else{				
					MPI_Isend(&snd_msg, 1, swap_instr_t, sq.first, 0, MPI_COMM_WORLD, &req);
				}
				swap_instr rec_msg={1, sq.first, sq.second};
				if(dq.first==0){//Im receiving
					instructions.push(rec_msg);
				}else{
					MPI_Isend(&rec_msg, 1, swap_instr_t, dq.first, 0, MPI_COMM_WORLD, &req);
				}
				bal_diff[sq.first]-=sq.second;
				bal_diff[dq.first]+=sq.second;
				dq.second += sq.second;
				surplus_queue.pop();
				sq=surplus_queue.front();

			}else{//same size
				swap_instr snd_msg={0, dq.first, sq.second};
				if (sq.first==0){//I'm sending
					instructions.push(snd_msg);
				}else{
					MPI_Isend(&snd_msg, 1, swap_instr_t, sq.first, 0, MPI_COMM_WORLD, &req);
				}
				swap_instr rec_msg={1, sq.first, sq.second};
				if (dq.first==0){//I'm receiving
					instructions.push(rec_msg);
				}else{
					MPI_Isend(&rec_msg, 1, swap_instr_t, dq.first, 0, MPI_COMM_WORLD, &req);
				}
				bal_diff[sq.first]-=sq.second;
				bal_diff[dq.first]+=sq.second;
				surplus_queue.pop();
				deficit_queue.pop();
				sq=surplus_queue.front();
				dq=deficit_queue.front();
			}
		}
		//After loop break, all instructions have been sent, send break instr to all nodes
		swap_instr end_msg={-1,-1,-1};
		for (int i=1; i<nProcs; i++){
			MPI_Isend(&end_msg, 1, swap_instr_t, i, 0, MPI_COMM_WORLD, &req);
		}
	}else{//Wait for the master precious

		MPI_Status recv_status, msg_status;
		MPI_Request req;
		swap_instr msg;

		while (true){//Receive instructions from the master until told to stop 
			MPI_Irecv(&msg, 1, swap_instr_t, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &req);
			MPI_Wait(&req, &recv_status);
			if (msg.action==-1){break;}
			else{
				instructions.push(msg);
			}
		}
	}

//****************************************************************************************************************************************
	//fill rebalanced arrays up to optimal or mycount
	int pos;
	if (myDif >= 0){//I've got enough or more
		for (int i=0; i<bal_counts[rank]; i++){
			(*rebalancedData)[i]=data[i];
		}  
		pos=myCount;                             
		
	}else if (myDif < 0){//I've got a deficit
		for (int i=0; i<myCount; i++){
			(*rebalancedData)[i]=data[i];
		}
		pos=0;
	}
	MPI_Barrier(MPI_COMM_WORLD);

//****************************************************************************************************************************************
	//Loop through queues of message instructions 
	//indexer pos approaches optimal, recv adds, sends subtract

	MPI_Status recv_status;
	MPI_Request recv_req, send_req;
	dist_sort_t *tmp = (dist_sort_t*)malloc(sizeof(dist_sort_t) * abs(myDif));
	swap_instr msg;
	while(!instructions.empty()){
		msg = instructions.front();
		instructions.pop();

		if (msg.action==1){//receiving
			//recv into tmp, update position
			MPI_Irecv(&tmp[pos], msg.num, MPI_TYPE_DIST_SORT_T, msg.idx, 2, MPI_COMM_WORLD, &recv_req);
			MPI_Wait(&recv_req, &recv_status);
			pos+=msg.num;
		}else if (msg.action==0){//sending
			//send, update pos
			pos-=msg.num;
			MPI_Isend(&data[pos], msg.num, MPI_TYPE_DIST_SORT_T, msg.idx, 2, MPI_COMM_WORLD, &send_req);
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);

	//Update those with deficits from tmp arrays
	int c=0;
	if (myDif < 0){
		for (int i=0; i<abs(myDif); i++){
			(*rebalancedData)[i+myCount] = tmp[i];
		}
	}
	
	free(tmp);
	free(bal_counts);
	free(bal_diff);
	free(counts);

	MPI_Barrier(MPI_COMM_WORLD);
//****************************************************************************************************************************************
}

//***************************************************************************************************************************************************************************************************************
//_______________________________________________________Find Splitters___________________________________________________________________
//***************************************************************************************************************************************************************************************************************

void findSplitters(const dist_sort_t *data, const dist_sort_size_t data_size, dist_sort_t *splitters, dist_sort_size_t *counts, int numSplitters) {

//****************************************************************************************************************************************
	//Setup

	int nProcs, rank;
	MPI_Comm_size(MPI_COMM_WORLD, &nProcs);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	MPI_Datatype probe_bounds_t;
    MPI_Type_contiguous(2, MPI_TYPE_DIST_SORT_T, &probe_bounds_t);
    MPI_Type_commit(&probe_bounds_t);

	MPI_Barrier(MPI_COMM_WORLD);

	//Global Item Count
	dist_sort_size_t g_count;
	MPI_Allreduce(&data_size, &g_count, 1, MPI_TYPE_DIST_SORT_SIZE_T, MPI_SUM, MPI_COMM_WORLD);
	int global_count = (int)g_count;
	
	dist_sort_t *sorted_data = (dist_sort_t*)malloc(sizeof(dist_sort_t) * data_size);
	memcpy(sorted_data, data, sizeof(dist_sort_t)*data_size);
	sort(sorted_data, data_size);

	dist_sort_t my_min = sorted_data[0];
	dist_sort_t my_max = sorted_data[data_size-1];
	MPI_Barrier(MPI_COMM_WORLD);
	//Approximate real number of items that should be in each bin
	int b_size = global_count / numSplitters;
	float cutoff = (float)b_size*0.001;
//****************************************************************************************************************************************
	//Ranks Share Lower and Upper bounds, Root finds global bounds and broadcasts

	probe_bounds *pbs = (probe_bounds*)malloc(sizeof(probe_bounds) * nProcs);
	probe_bounds my_bounds={my_min, my_max};
	MPI_Win window;	
	MPI_Win_create(pbs, nProcs*sizeof(probe_bounds), sizeof(probe_bounds), MPI_INFO_NULL, MPI_COMM_WORLD, &window);
	MPI_Win_fence(0, window);
		MPI_Put(&my_bounds, 1, probe_bounds_t, 0, rank, 1, probe_bounds_t, window);
	MPI_Win_fence(0, window);
	MPI_Win_free(&window);

	probe_bounds global_bounds;
	if (rank==0){
		global_bounds={my_min, my_max};
		for (int i=0; i<nProcs; i++){
			if (global_bounds.lower > pbs[i].lower){
				global_bounds.lower = pbs[i].lower;
			}
			if (global_bounds.upper < pbs[i].upper){
				global_bounds.upper = pbs[i].upper;
			}
		}
	}
	MPI_Bcast(&global_bounds, 1, probe_bounds_t, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
//****************************************************************************************************************************************
	//Big Histograming Loop
	//It's like a blender. 
	
	dist_sort_t *probe_keys=(dist_sort_t*)malloc(sizeof(dist_sort_t) * numSplitters);
	dist_sort_t *probe_counts=(dist_sort_t*)malloc(sizeof(dist_sort_t) * numSplitters * nProcs);
	dist_sort_size_t *totalcounts=(dist_sort_size_t*)malloc(sizeof(dist_sort_size_t) * numSplitters);
	float *splits=(float*)malloc(sizeof(float) * numSplitters);
	float *opt_split=(float*)malloc(sizeof(float) * numSplitters);
	dist_sort_t bound = (global_bounds.upper - global_bounds.lower);
	for (int i=0; i<numSplitters; i++){
		opt_split[i]=(1/(float)numSplitters);
	}
	
	int hist_iter=0;
	int iter_max=5;
	int big_gate=0;
	MPI_Barrier(MPI_COMM_WORLD);

	while(true){
		//Loop Process:
		//Root Calculates key intervals, Fills Keys, Shares window, Broadcasts to everyone
		//Essentially, Root sends out beacons of probe_keys
		//Nodes check the new stuff and report back, Root returns or rerolls

		if (rank==0){//Rank 0 calc initial Keys, interval of even spaced keys from min to max
			if (hist_iter==0){
				//First Iteration: Fill splits with Even Distribution over bound
				for (int i=0; i<numSplitters; i++){
					splits[i]=opt_split[i];
				}
			}//All other iterations: Rank 0 adjusts splits based on results
			dist_sort_t splitpoint;
			for (int i=0; i<numSplitters-1; i++){
				splitpoint=(float)bound*splits[i];
				if (i==0){probe_keys[i] = global_bounds.lower+splitpoint;}
				else{probe_keys[i] = probe_keys[i-1]+splitpoint;}
			}
			probe_keys[numSplitters-1]=global_bounds.upper;
		}

		//Rank 0 Broadcast probe_keys
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Bcast(probe_keys, numSplitters, MPI_TYPE_DIST_SORT_T, 0, MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);

		//All Ranks Count Bins from probe keys
		dist_sort_size_t *my_counts=(dist_sort_size_t*)malloc(sizeof(dist_sort_size_t) * numSplitters);
		for (int i=0; i<numSplitters; i++){my_counts[i]=0;}

		int pos=0;
		for (int i=0; i<data_size; i++){
			if (sorted_data[i] > probe_keys[pos]){pos+=1;}
			my_counts[pos]+=1;
		}

		//All Ranks count and update probe_counts at 0 for their rank index
		MPI_Win_create(probe_counts, sizeof(dist_sort_t) * numSplitters * nProcs, sizeof(dist_sort_t), MPI_INFO_NULL, MPI_COMM_WORLD, &window);
		MPI_Win_fence(0, window);
		MPI_Put(&my_counts[0], numSplitters, MPI_TYPE_DIST_SORT_SIZE_T, 0, (numSplitters * rank), numSplitters, MPI_TYPE_DIST_SORT_SIZE_T, window);
		MPI_Win_fence(0, window);
		MPI_Win_free(&window);
		MPI_Barrier(MPI_COMM_WORLD);
		
		if (rank==0){
			//LOOP CUT OFF VALVE
			//Root determines if cutoff, else adjusts keys, rinse and repeat
			
			for (int i=0; i<nProcs; i++){totalcounts[i]=0;}
			for (int i=0; i<nProcs; i++){
				for (int j=0; j<numSplitters; j++){
					totalcounts[j]+=probe_counts[(i*nProcs) + j];
				}
			}
			//if all within the cutoff margin, break
			int countdiff;
			int lilgate=1;
			for (int i=0; i<nProcs; i++){
				countdiff=abs(totalcounts[i]-b_size);
				if (countdiff>cutoff){
					lilgate=0;
				}
			}
			if (lilgate==1){//Kill Loop
				big_gate=1;
			}
			
			//Next iteration splits
			float new_split;	
			for (int i=0; i<numSplitters; i++){
				new_split = splits[i]+(opt_split[i]-((float)totalcounts[i]/(float)global_count));
				splits[i]=new_split;
			}
		}

		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Bcast(&big_gate, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);
		hist_iter++;
		if (big_gate==1){
			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Bcast(totalcounts, numSplitters, MPI_TYPE_DIST_SORT_SIZE_T, 0, MPI_COMM_WORLD);
			MPI_Barrier(MPI_COMM_WORLD);		
			memcpy(counts, totalcounts, sizeof(dist_sort_size_t) * numSplitters);
			memcpy(splitters, probe_keys, sizeof(dist_sort_t) * numSplitters);	
			break;}
		if (hist_iter>iter_max){break;}
	}//Hist Loop
	
	free(sorted_data);
	free(probe_keys);
	free(probe_counts);
	free(totalcounts);
	free(splits);
	free(opt_split);

	MPI_Barrier(MPI_COMM_WORLD);
}


//***************************************************************************************************************************************************************************************************************
//________________________________________________________Move Data_______________________________________________________________________
//***************************************************************************************************************************************************************************************************************
void moveData(const dist_sort_t *const sendData, const dist_sort_size_t sDataCount,
		dist_sort_t **recvData, dist_sort_size_t *rDataCount,
		const dist_sort_t *const splitters, const dist_sort_t *const counts, int numSplitters) {
//****************************************************************************************************************************************
	//Setup
	int nProcs, rank;
	MPI_Comm_size(MPI_COMM_WORLD, &nProcs);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Win window;

	dist_sort_t *split_keys=(dist_sort_t*)malloc(sizeof(dist_sort_t) * numSplitters);
	dist_sort_size_t *bin_counts=(dist_sort_size_t*)malloc(sizeof(dist_sort_size_t) * numSplitters);

	if (rank==0){
		memcpy(split_keys, splitters, sizeof(dist_sort_t) * numSplitters);
		memcpy(bin_counts, counts, sizeof(dist_sort_size_t) * numSplitters);	
	}

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Bcast(bin_counts, numSplitters, MPI_TYPE_DIST_SORT_SIZE_T, 0, MPI_COMM_WORLD);
	MPI_Bcast(split_keys, numSplitters, MPI_TYPE_DIST_SORT_T, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
//****************************************************************************************************************************************
	//Create indexes to move data correctly, everyone get on the same page

	int *split_indexes=(int*)malloc(sizeof(int) * numSplitters);
	int *mybins=(int*)malloc(sizeof(int)*numSplitters);
	memset(split_indexes, 0, numSplitters);
	
	//Index my bins relative to my sendData
	int key=0;
	for (int i=0; i<sDataCount; i++){
		if(sendData[i] > split_keys[key]){
			split_indexes[key]=i;
			if (key==0){
				mybins[key]=i;
			}else{
				mybins[key]=split_indexes[key]-split_indexes[key-1];
			}
			key++;
			if(key>=numSplitters-1){
				split_indexes[key]=sDataCount;//last cutoff is last index
				mybins[key]=split_indexes[key]-split_indexes[key-1];
				break;
			}
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	int *allbins=(int*)malloc(sizeof(int) * numSplitters * nProcs);
	MPI_Win_create(allbins, sizeof(int)*numSplitters*nProcs, sizeof(int), MPI_INFO_NULL, MPI_COMM_WORLD, &window);
	MPI_Win_fence(0, window);
	MPI_Put(&mybins[0], numSplitters, MPI_INT, 0, (rank*numSplitters), numSplitters, MPI_INT, window);
	MPI_Win_fence(0, window);
	MPI_Win_free(&window);

	MPI_Barrier(MPI_COMM_WORLD);
	int *allindexes=(int*)malloc(sizeof(int) * numSplitters * nProcs);
	MPI_Win_create(allindexes, sizeof(int)*numSplitters*nProcs, sizeof(int), MPI_INFO_NULL, MPI_COMM_WORLD, &window);
	MPI_Win_fence(0, window);
	MPI_Put(&split_indexes[0], numSplitters, MPI_INT, 0, (rank*numSplitters), numSplitters, MPI_INT, window);
	MPI_Win_fence(0, window);
	MPI_Win_free(&window);

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Bcast(allbins, numSplitters*nProcs, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(allindexes, numSplitters*nProcs, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	
	int toget=0;
	for (int i=0; i<nProcs; i++){
		toget+=allbins[(i*nProcs)+rank];
	}

	MPI_Barrier(MPI_COMM_WORLD);
//****************************************************************************************************************************************
	//Data Movement
	//Two communication options implemented, One-Sided and Point-to-Point

	int bin_pos=0;
	int amt, offset;
	*rDataCount=(dist_sort_size_t)toget;
	*recvData=(dist_sort_t*)malloc(sizeof(dist_sort_t)*toget);

	int move_option=2;
	if (rank==0){
		if (move_option==1){
			printf("\nRESULTS: One-Sided\n");
		}else if(move_option==2){
			printf("\nRESULTS: Point-Point\n");
		}	
	}

	if (move_option==1){//One-Sided Communication
		MPI_Win_create((void*)sendData, sizeof(dist_sort_t)*sDataCount, sizeof(dist_sort_t), MPI_INFO_NULL, MPI_COMM_WORLD, &window);
		MPI_Win_fence(0, window);
		for (int i=0; i<nProcs; i++){
				amt = allbins[(i*nProcs)+rank];
				if (rank==0){
					offset = 0;
				}else{
					offset = allindexes[(i*nProcs)+(rank-1)];
				}
				MPI_Get(&(*recvData)[bin_pos], amt, MPI_TYPE_DIST_SORT_T, i, offset, amt, MPI_TYPE_DIST_SORT_T, window);
				bin_pos+=allbins[(i*nProcs)+rank];
		}
		MPI_Win_fence(0, window);
		MPI_Win_free(&window);
		MPI_Barrier(MPI_COMM_WORLD);
	
	}else if(move_option==2){//Point to Point Communication
		MPI_Status status;
		MPI_Request req;
		MPI_Barrier(MPI_COMM_WORLD);
		for (int r=0; r<nProcs; r++){
			if (rank==r){//recving
				for (int i=0; i<nProcs; i++){
					amt = allbins[(i*nProcs)+rank];
					if (i!=rank){
						MPI_Recv(&(*recvData)[bin_pos], amt, MPI_TYPE_DIST_SORT_T, i, 1, MPI_COMM_WORLD, &status);
						bin_pos+=amt;
					}else{//My send data to my recv data
						if (r==0){
							offset = 0;
						}else{
							offset = split_indexes[rank-1];
						}
						for (int j=0; j<mybins[rank]; j++){
							(*recvData)[bin_pos+j]=sendData[offset+j];
						}
						bin_pos+=amt;
					}
				}
			}else{//sending
				if (r==0){
					offset = 0;
				}else{
					offset = split_indexes[r-1];
				}
				amt = mybins[r];
				MPI_Isend(&sendData[offset], amt, MPI_TYPE_DIST_SORT_T, r, 1, MPI_COMM_WORLD, &req);
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}

	free(split_indexes);
	free(split_keys);
	free(mybins);
	free(allbins);
	free(allindexes);

	MPI_Barrier(MPI_COMM_WORLD);	
}

//***************************************************************************************************************************************************************************************************************
//_____________________________________________________________Sort_______________________________________________________________________
//****************************************************************************************************************************************
void sort(dist_sort_t *data, dist_sort_size_t size) {
	// You are welcome to use this sort function.
	// If you wish to implement your own, you can do that too.
	// Don't use bubblesort.
	std::sort(data,data+size);
}

