
#include <stdlib.h>
#include <math.h>

#include "utils.h"
#include "MCMCcontrol.h"
#include "AlignmentProcessor.h"
#include "GenericTree.h"
#include "PopulationTree.h"
#include "LocusDataLikelihood.h"
#include "patch.h"

#include "CombStats.h"



// --- FUNCTION IMPLEMENTATIONS -----------------------------------------------


void calculateCombStats() {
	for (int comb = 0 ; comb < dataSetup.popTree->numPops ; comb++){
		if (isFeasableComb(comb)){
			initSpecificCombStats(comb);
			combNumCoals(comb);
//			computeCombCoalStats(comb);
		}
	}
}




void combNumCoals(int comb){
	combNumCoals_rec(comb, comb);
}

void combNumCoals_rec(int comb, int currentPop){
	if (isLeaf(currentPop)){
		countLeafCoals(comb, currentPop);
	} else{
		countNonLeafCoals(comb, currentPop);
	}
}

void countLeafCoals(int comb, int leaf) {
	for(int gene=0; gene<dataSetup.numLoci; gene++) {
		countLeafGeneCoals(comb, leaf, gene);
	}
	if (comb_stats[comb].leaves[leaf].debug_original_num_coals !=
			genetree_stats_total.num_coals[leaf]) { // TODO - debugggggggggg
				printf("ERROR: comb leaf num coal not correct!");
				exit(-1);
	}
}

void countLeafGeneCoals(int comb, int leaf, int gene){

	float combAge = comb_stats[comb].age;
	int belowCombAge = 0;
	int aboveCombAge = 0;
	int eventId = event_chains[gene].first_event[leaf];
	float eventAge = 0.0;

	while (eventId >= 0) {
		eventAge += event_chains[gene].events[eventId].elapsed_time;

		if (event_chains[gene].events[eventId].type == COAL){
			if (eventAge <= combAge){
				belowCombAge++;
			} else {
				aboveCombAge++; //TODO - verify with Ilan - if a coal-event occurs on comb-age, to which Pop does it belong? - - - currently to the LEAF
			}
		}

		eventId = event_chains[gene].events[eventId].next;
	}

	comb_stats[comb].leaves[leaf].num_coals_total += belowCombAge;
	comb_stats[comb].num_coals_total += aboveCombAge;
	comb_stats[comb].leaves[leaf].debug_original_num_coals += (belowCombAge + aboveCombAge); // TODO - remove debug statement
}


void countNonLeafCoals(int comb, int currentPop) {

	int leftSon = dataSetup.popTree->pops[currentPop]->sons[LEFT]->id;
	int rightSon = dataSetup.popTree->pops[currentPop]->sons[RIGHT]->id;

	combNumCoals_rec(comb, leftSon);
	combNumCoals_rec(comb, rightSon);

	comb_stats[comb].num_coals_total += genetree_stats_total.num_coals[currentPop];
}


//void computeCombCoalStats(){
//	for(int gen=0; gen<dataSetup.numLoci; gen++) {
//		computeCombCoalStats_rec(dataSetup.popTree->rootPop, gen);
//	}
//}

//void computeCombCoalStats_rec(int comb, int gen) {
//	int leftSon, rightSon;
//
//	if (isLeaf(comb)){
//		fillupLeafCombStats(comb, gen);
//	} else{
//		leftSon = dataSetup.popTree->pops[comb]->sons[LEFT]->id;
//		rightSon = dataSetup.popTree->pops[comb]->sons[RIGHT]->id;
//
//		computeCombCoalStats_rec(leftSon, gen);
//		computeCombCoalStats_rec(rightSon, gen);
//
//		fillupCombStats(comb, gen);
//
//	}
////	debug_print_combstats(comb, gen); //CAUTION! uncommenting this increases running time tenfold! (because of IO)
//}

//void fillupLeafCombStats(int comb, int gen){
//	appendPopToComb(comb, gen, 0); // since this is a leaf, starting point of the comb_stats arrays should be 0
//}

//void appendPopToComb(int comb, int gen, int startingPoint){
//	int i = startingPoint;
//	int event = event_chains[gen].first_event[comb];
//	double combStartTime = dataSetup.popTree->pops[comb]->age;
//	double eventAge = combStartTime;
//
//	for ( ; event >= 0 ; i++, event = event_chains[gen].events[event].next){
//		eventAge += event_chains[gen].events[event].elapsed_time;
//		comb_stats[comb].sorted_ages[i] = eventAge;
//		comb_stats[comb].elapsed_times[i] = event_chains[gen].events[event].elapsed_time;
//		comb_stats[comb].num_lineages[i] = event_chains[gen].events[event].num_lineages;
//		comb_stats[comb].event_types[i] = event_chains[gen].events[event].type;
//	}
//	comb_stats[comb].num_events = i;
////	comb_stats[comb].coal_stats_total += genetree_stats[gen].coal_stats_total[comb]; // TODO - decide which of these lines is better
//	comb_stats[comb].coal_stats_total += getCoalStats(comb_stats[comb].elapsed_times + startingPoint, comb_stats[comb].num_lineages + startingPoint, comb_stats[comb].num_events - startingPoint);
//
//#ifdef CHECKCOMB
//	test_compareGphocsVsCombPopCoalStats(comb, gen, startingPoint);
//#endif
//}


//void fillupCombStats(int comb, int gen){
//	addChildenIntoCombStats(comb, gen);
//	addCurrentPopIntoCombStats(comb, gen);
//}

//void addChildenIntoCombStats(int comb, int gen){
//	mergeChildern(comb, gen);
//	addChildrenCombStats(comb, gen);
//}



//void mergeChildern(int comb, int gen){
//	int i, j = 0, k = 0;
//	double leftAge, rightAge;
//	int leftSon, rightSon;
//	leftSon = dataSetup.popTree->pops[comb]->sons[LEFT]->id;
//	rightSon = dataSetup.popTree->pops[comb]->sons[RIGHT]->id;
//	int m = comb_stats[leftSon].num_events;
//	int n = comb_stats[rightSon].num_events;
//
//	for (i = 0 ; i < m + n; ) {
//		if ( j < m && k < n) {
//			leftAge = comb_stats[leftSon].sorted_ages[j];
//			rightAge = comb_stats[rightSon].sorted_ages[k];
//			if (leftAge < rightAge){
//				comb_stats[comb].event_types[i] = comb_stats[leftSon].event_types[j];
//				comb_stats[comb].sorted_ages[i] = comb_stats[leftSon].sorted_ages[j];
//				comb_stats[comb].num_lineages[i] =
//						comb_stats[leftSon].num_lineages[j] + comb_stats[rightSon].num_lineages[k];
//				j++;
//			}
//			else{
//				comb_stats[comb].event_types[i] = comb_stats[rightSon].event_types[k];
//				comb_stats[comb].sorted_ages[i] = comb_stats[rightSon].sorted_ages[k];
//				comb_stats[comb].num_lineages[i] =
//						comb_stats[leftSon].num_lineages[j] + comb_stats[rightSon].num_lineages[k];
//				k++;
//			}
//			i++;
//		}
//		else if (j == m) {
//			for (; i < m + n ;) {
//				comb_stats[comb].event_types[i] = comb_stats[rightSon].event_types[k];
//				comb_stats[comb].sorted_ages[i] = comb_stats[rightSon].sorted_ages[k];
//				comb_stats[comb].num_lineages[i] =
//						comb_stats[leftSon].num_lineages[j-1] + comb_stats[rightSon].num_lineages[k];
//				k++;
//				i++;
//			}
//		}
//		else {
//			for (; i < m + n;) {
//				comb_stats[comb].event_types[i] = comb_stats[leftSon].event_types[j];
//				comb_stats[comb].sorted_ages[i] = comb_stats[leftSon].sorted_ages[j];
//				comb_stats[comb].num_lineages[i] =
//						comb_stats[leftSon].num_lineages[j] + comb_stats[rightSon].num_lineages[k-1];
//				j++;
//				i++;
//			}
//		}
//	}
//	for (i = 0 ; i < m + n; i++ ) {
//		comb_stats[comb].elapsed_times[i] = (i==0) ? 0 :
//				comb_stats[comb].sorted_ages[i] - comb_stats[comb].sorted_ages[i-1];
//	}
//
//	comb_stats[comb].num_events = m + n;
//}

//void addChildrenCombStats(int comb, int gen){
//	comb_stats[comb].coal_stats_total +=
//			getCoalStats(comb_stats[comb].elapsed_times, comb_stats[comb].num_lineages, comb_stats[comb].num_events);
//}

//void addCurrentPopIntoCombStats(int comb, int gen){
//	appendPopToComb(comb, gen, comb_stats[comb].num_events); // start filling the comb_stats arrays from the last known event
//}


//double getCoalStats(double* elapsed_times, int* num_lineages, int size){
//	int n;
//	double t;
//	double result = 0.0;
//	for( int i = 0 ; i < size ; i++) {
//		n = num_lineages[i];
//		t = elapsed_times[i];
//		result += n*(n-1)*t;
//	}
//	return result;
//}













double getCombAge(int comb){

	if (isLeaf(comb)){
		return 99999.9; // TODO - replace with MAX_DOUBLE const
	}

	if (areChildrenLeaves(comb)){
		return dataSetup.popTree->pops[comb]->age;
	}

	if (isFeasableComb(comb)){
		double left_min = getCombAge(dataSetup.popTree->pops[comb]->sons[LEFT]->id);
		double right_min = getCombAge(dataSetup.popTree->pops[comb]->sons[RIGHT]->id);
		return fmin(left_min, right_min);
	}

	printf("ERROR: bug in combAge algorithm. Should not reach here!");
	exit(-1);
	return 99999.9;
}
int	isLeaf(int pop){
	Population *population, *left_son, *right_son;

	population = dataSetup.popTree->pops[pop];

	left_son = population->sons[LEFT];
	right_son = population->sons[RIGHT];

	if (left_son || right_son){ // TODO - validate this is a valid boolean test
		return FALSE;
	} else {
		return TRUE;
	}
}
int areChildrenLeaves(int pop){
	if (isLeaf(pop)){
		return FALSE;
	}
	int left_son = dataSetup.popTree->pops[pop]->sons[LEFT]->id;
	int right_son = dataSetup.popTree->pops[pop]->sons[RIGHT]->id;
	return (isLeaf(left_son) && isLeaf(right_son));
}
int isFeasableComb(int pop){
	if (isLeaf(pop)){
		return FALSE;
	}
	if (areChildrenLeaves(pop)) { //TODO - ask Ilan: should AB:A,B be a feasable comb?
		return FALSE;
	}
	return TRUE;
}
void initSpecificCombStats(int comb){
	comb_stats[comb].coal_stats_total   = 0.0;
	comb_stats[comb].mig_stats_total 	= 0.0;
	comb_stats[comb].num_coals_total 	= 0;
	comb_stats[comb].num_migs_total 	= 0;
	comb_stats[comb].num_events 		= 0;
	comb_stats[comb].debug_total_error  = 0;
	comb_stats[comb].age 				= getCombAge(comb);

	for (int leaf = 0 ; leaf < dataSetup.popTree->numCurPops ; leaf++){
		comb_stats[comb].leaves[leaf].num_migs_total = 0;
		comb_stats[comb].leaves[leaf].num_coals_total = 0;
		comb_stats[comb].leaves[leaf].mig_stats = 0.0;
		comb_stats[comb].leaves[leaf].coal_stats = 0.0;
		comb_stats[comb].leaves[leaf].debug_original_num_coals = 0; // TODO - remove if not in-use
	}

}
void allocateCombMem(){
	comb_stats=(struct COMB_STATS*)malloc(dataSetup.popTree->numPops*sizeof(struct COMB_STATS));
	int max_events = 2*dataSetup.numSamples+ 4*MAX_MIGS + 3*dataSetup.popTree->numMigBands + dataSetup.popTree->numPops + 10;

	for (int comb = 0; comb < dataSetup.popTree->numPops; comb++) { // TODO - consider not iterating over non-comb-feasable pops
		comb_stats[comb].sorted_ages = (double*)malloc(max_events*sizeof(double));
		comb_stats[comb].elapsed_times = (double*)malloc(max_events*sizeof(double));
		comb_stats[comb].num_lineages = (int*)malloc(max_events*sizeof(int));
		comb_stats[comb].event_types = (int*)malloc(max_events*sizeof(int));
		comb_stats[comb].coal_stats_total = 0.0;

		comb_stats[comb].mig_stats_total = 0.0;

		comb_stats[comb].leaves=(struct COMB_LEAF*)malloc(dataSetup.popTree->numCurPops*sizeof(struct COMB_LEAF));
		for (int leaf = 0 ; leaf < dataSetup.popTree->numCurPops ; leaf++){
			comb_stats[comb].leaves[leaf].coal_stats = 0.0;
			comb_stats[comb].leaves[leaf].mig_stats = 0.0;
			comb_stats[comb].leaves[leaf].num_migs_total = 0;
			comb_stats[comb].leaves[leaf].num_coals_total = 0;
		}


		if((comb_stats == NULL) || // TODO - add error checks for comb_stats_leaves
				(!comb_stats[comb].sorted_ages) || (!comb_stats[comb].elapsed_times) ||
				(!comb_stats[comb].num_lineages)|| (!comb_stats[comb].event_types)) {
			fprintf(stderr, "\nError: Out Of Memory comb_stats\n");
			exit(-1);
		}
	}

}
void freeCombMem(){ // TODO - implement
}
void debug_print_errors(){
	for (int comb = 0; comb < dataSetup.popTree->numPops; comb++) {
		fprintf(stderr, "comb:%s, total_error=%0.35f, relative_error=%0.35f\n",
				dataSetup.popTree->pops[comb]->name , comb_stats[comb].debug_total_error,
				comb_stats[comb].debug_total_error / comb_stats[comb].coal_stats_total);
	}
	printf("\n");
}
void debug_print_combstats(int comb, int gen){
	printf("=== comb:%s, gen:%d, num_events:%d, num_coals_total:%d === >>\n",
			dataSetup.popTree->pops[comb]->name, gen, comb_stats[comb].num_events, comb_stats[comb].num_coals_total);
	for (int i = 0 ; i < comb_stats[comb].num_events ; i++){
		printf("age:%0.25f, elapsed_time:%0.25f, num_lin:%d, type:%s \n",
				comb_stats[comb].sorted_ages[i], comb_stats[comb].elapsed_times[i],
				comb_stats[comb].num_lineages[i],	/*getEventTypeName(comb_stats[comb].event_types[i])*/ "UNIMPLEMENTED");
	}
	printf("====================================================================================\n");
	fflush(stdout);
}
