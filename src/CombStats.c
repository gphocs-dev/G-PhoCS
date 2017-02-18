
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
	initCombStats();
	for (int comb = 0 ; comb < dataSetup.popTree->numPops ; comb++){
		if (isFeasibleComb(comb)){
			for(int gene=0; gene<dataSetup.numLoci; gene++) {
				calculateSufficientStats(comb, gene);
			}

		}
	}
	assertRootNumCoals();
}


void calculateSufficientStats(int comb, int gene){
	coalescence(comb, gene);
//	migrations(comb, gene);
}

void coalescence(int comb, int gene){
	coalescence_rec(comb, comb, gene);

}

void coalescence_rec(int comb, int currentPop, int gene){
	if (isLeaf(currentPop)){
		handleLeafCoals(comb, currentPop, gene);
	} else{
		coalescence_rec(comb, getSon(currentPop, LEFT), gene);
		coalescence_rec(comb, getSon(currentPop, RIGHT), gene);

		handleNonLeafCoals(comb, currentPop, gene);
	}
}


void handleLeafCoals(int comb, int leaf, int gene) {

	double elapsedTime, eventAge = 0.0;
	int numLineages, eventType, numEventsAboveComb = 0;
	double combAge = comb_stats[comb].age;

	Stats* belowCombLeafStats = &comb_stats[comb].leaves[leaf].below_comb;
	Stats* aboveCombLeafStats = &comb_stats[comb].leaves[leaf].above_comb;
	Stats* combTotalStats = &comb_stats[comb].total;

	int eventId = event_chains[gene].first_event[leaf];
	while (eventId >= 0) {
		elapsedTime = event_chains[gene].events[eventId].elapsed_time;
		eventType = event_chains[gene].events[eventId].type;
		numLineages = event_chains[gene].events[eventId].num_lineages;
		eventAge += elapsedTime;

		if (eventAge <= combAge){
			if (eventType == COAL){
				belowCombLeafStats->num_coals++;
			}
			// you've got all the data you need so just update coal_stats -
			belowCombLeafStats->coal_stats += numLineages*(numLineages-1)*elapsedTime;  // TODO - discuss this with Ilan (how I update coal stats on EVERY event and not just numLin-changing events)
		}
		else {
			if (eventType == COAL){
				aboveCombLeafStats->coal_stats++;
				combTotalStats->num_coals++;
			}
			// you're a tiny part of the comb so supply future methods the data they need -
			aboveCombLeafStats->sorted_ages[numEventsAboveComb] 	= eventAge;
			aboveCombLeafStats->elapsed_times[numEventsAboveComb] 	= elapsedTime;
			aboveCombLeafStats->num_lineages[numEventsAboveComb] 	= numLineages;
			aboveCombLeafStats->event_types[numEventsAboveComb] 	= eventType;
			numEventsAboveComb++;
		}
		eventId = event_chains[gene].events[eventId].next;
	}
	aboveCombLeafStats->num_events = numEventsAboveComb;
	combTotalStats->num_events += numEventsAboveComb;
}

void handleNonLeafCoals(int comb, int currentPop, int gene) {
	handleNonLeafNumCoals(comb, currentPop, gene);
//	handleNonLeafCoalStats(comb, currentPop, gene); // TODO - reactivate coal stats
}

void handleNonLeafNumCoals(int comb, int currentPop, int gene) {
	comb_stats[comb].total.num_coals += genetree_stats[gene].num_coals[currentPop];
}

void handleNonLeafCoalStats(int comb, int currentPop, int gene){
	mergeChildernIntoCurrent(comb, currentPop, gene);
	appendCurrent(comb, currentPop, gene); // start filling the comb_stats arrays from the last known event
}

void mergeChildernIntoCurrent(int comb, int currentPop, int gen){
	int i, j = 0, k = 0;
	double leftAge, rightAge;
	int leftSon, rightSon;
	Stats *leftStats, *rightStats, *currentStats;

	leftSon = getSon(currentPop, LEFT);
	rightSon = getSon(currentPop, RIGHT);

	currentStats = getStats(comb, currentPop);
	leftStats = getStats(comb, leftSon);
	rightStats = getStats(comb, rightSon);


	int m = leftStats->num_events;
	int n = rightStats->num_events;

	for (i = 0 ; i < m + n; ) {
		if ( j < m && k < n) {

			leftAge = leftStats->sorted_ages[j];
			rightAge = rightStats->sorted_ages[k];
			if (leftAge < rightAge){
				currentStats->event_types[i]  = leftStats->event_types[j];
				currentStats->sorted_ages[i]  = leftStats->sorted_ages[j];
				currentStats->num_lineages[i] = leftStats->num_lineages[j] + rightStats->num_lineages[k];
				j++;
			}
			else{
				currentStats->event_types[i]  = rightStats->event_types[j];
				currentStats->sorted_ages[i]  = rightStats->sorted_ages[j];
				currentStats->num_lineages[i] = leftStats->num_lineages[j] + rightStats->num_lineages[k];
				k++;
			}
			i++;
		}
		else if (j == m) {
			for (; i < m + n ;) {
				currentStats->event_types[i] = rightStats->event_types[k];
				currentStats->sorted_ages[i] = rightStats->sorted_ages[k];
				currentStats->num_lineages[i] =
						leftStats->num_lineages[j-1] + rightStats->num_lineages[k];
				k++;
				i++;
			}
		}
		else {
			for (; i < m + n;) {
				currentStats->event_types[i] = leftStats->event_types[j];
				currentStats->sorted_ages[i] = leftStats->sorted_ages[j];
				currentStats->num_lineages[i] =
						leftStats->num_lineages[j] + rightStats->num_lineages[k-1];
				j++;
				i++;
			}
		}
	}
	for (i = 0 ; i < m + n; i++ ) {
		currentStats->elapsed_times[i] = (i==0) ? 0 :
				currentStats->sorted_ages[i] - currentStats->sorted_ages[i-1];
	}

	currentStats->num_events = m + n;
}

void appendCurrent(int comb, int currentPop, int gene){
	Stats *currentStats = getStats(comb, currentPop);

	int startingPoint = currentStats->num_events;
	int i = startingPoint;
	int event = event_chains[gene].first_event[currentPop];
	double combStartTime = dataSetup.popTree->pops[currentPop]->age;
	double eventAge = combStartTime;

	for ( ; event >= 0 ; i++, event = event_chains[gene].events[event].next){
		eventAge += event_chains[gene].events[event].elapsed_time;
		currentStats->sorted_ages[i] = eventAge;
		currentStats->elapsed_times[i] = event_chains[gene].events[event].elapsed_time;
		currentStats->num_lineages[i] = event_chains[gene].events[event].num_lineages;
		currentStats->event_types[i] = event_chains[gene].events[event].type;
	}
	currentStats->num_events = i;
}

double calculateCoalStats(double* elapsed_times, int* num_lineages, int size){
	int n;
	double t;
	double result = 0.0;
	for( int i = 0 ; i < size ; i++) {
		n = num_lineages[i];
		t = elapsed_times[i];
		result += n*(n-1)*t;
	}
	return result;
}












//TODO - find places for these functions
double getCombAge(int comb){
	if (isLeaf(comb)){
		return 99999.9; // TODO - replace with MAX_DOUBLE const
	} else if (areChildrenLeaves(comb)){
		return dataSetup.popTree->pops[comb]->age;
	} else if (isFeasibleComb(comb)){
		double left_min = getCombAge(dataSetup.popTree->pops[comb]->sons[LEFT]->id);
		double right_min = getCombAge(dataSetup.popTree->pops[comb]->sons[RIGHT]->id);
		return fmin(left_min, right_min);
	} else {
		printf("ERROR: bug in combAge algorithm. Should not reach here!");
		exit(-1);
		return 99999.9;
	}
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
int isFeasibleComb(int pop){
	if (isLeaf(pop)){
		return FALSE;
	} else if (areChildrenLeaves(pop)){
		return FALSE;
	} else {
		return TRUE;
	}
}

int getSon(int pop, int SON){
	return dataSetup.popTree->pops[pop]->sons[SON]->id;
}
Stats* getStats(int comb, int pop){ //TODO - rename to something more descriptive
	if (isLeaf(pop)){
		return &comb_stats[comb].leaves[pop].above_comb;
	} else {
		return &comb_stats[comb].clades[pop];
	}
}


void initCombStats(){
	for (int comb = 0 ; comb < dataSetup.popTree->numPops ; comb++){
		if (isFeasibleComb(comb)){

			initStats(&comb_stats[comb].total);
			comb_stats[comb].age = getCombAge(comb);

			for (int pop = 0 ; pop < dataSetup.popTree->numCurPops ; pop++){
				if (isLeaf(pop)){
					initStats(&comb_stats[comb].leaves[pop].above_comb);
					initStats(&comb_stats[comb].leaves[pop].below_comb);
				} else {
					initStats(&comb_stats[comb].clades[pop]);
				}
			}
		}
	}
}
void initStats(Stats* stats){
	stats->coal_stats = 0.0;
	stats->mig_stats = 0.0;
	stats->num_coals  = 0;
	stats->num_events = 0;
	stats->num_migs = 0;
}
void allocateCombMem(){
	comb_stats=malloc(dataSetup.popTree->numPops*sizeof(struct COMB_STATS));

	for (int comb = 0; comb < dataSetup.popTree->numPops; comb++) {
		allocateStats(&comb_stats[comb].total);
		comb_stats[comb].leaves=malloc(dataSetup.popTree->numCurPops*sizeof(LeafStats));
		comb_stats[comb].clades=malloc(dataSetup.popTree->numCurPops*sizeof(Stats));
		for (int pop = 0 ; pop < dataSetup.popTree->numCurPops ; pop++){
			if (isLeaf(pop)){
				allocateStats(&comb_stats[comb].leaves[pop].below_comb);
				allocateStats(&comb_stats[comb].leaves[pop].above_comb);
			} else {
				allocateStats(&comb_stats[comb].clades[pop]);
			}
		}
		if(comb_stats == NULL){ // TODO - add memory allocation test for all of comb_stats
			fprintf(stderr, "\nError: Out Of Memory comb_stats\n");
			exit(-1);
		}
	}
}
void allocateStats(Stats* stats){
	int max_events = 2*dataSetup.numSamples+ 4*MAX_MIGS + 3*dataSetup.popTree->numMigBands + dataSetup.popTree->numPops + 10;
	stats->sorted_ages = (double*)malloc(max_events*sizeof(double));
	stats->elapsed_times = (double*)malloc(max_events*sizeof(double));
	stats->num_lineages = (int*)malloc(max_events*sizeof(int));
	stats->event_types = (int*)malloc(max_events*sizeof(int));
}

void freeCombMem(){ // TODO - implement
}

//void debug_print_errors(){
//	for (int comb = 0; comb < dataSetup.popTree->numPops; comb++) {
//		fprintf(stderr, "comb:%s, total_error=%0.35f, relative_error=%0.35f\n",
//				dataSetup.popTree->pops[comb]->name , comb_stats[comb].debug_total_error,
//				comb_stats[comb].debug_total_error / comb_stats[comb].coal_stats_total);
//	}
//	printf("\n");
//}
//void debug_print_combstats(int comb, int gen){
//	printf("=== comb:%s, gen:%d, num_events:%d, num_coals_below_comb:%d === >>\n",
//			dataSetup.popTree->pops[comb]->name, gen, comb_stats[comb].num_events, comb_stats[comb].num_coals_total);
//	for (int i = 0 ; i < comb_stats[comb].num_events ; i++){
//		printf("age:%0.25f, elapsed_time:%0.25f, num_lin:%d, type:%s \n",
//				comb_stats[comb].sorted_ages[i], comb_stats[comb].elapsed_times[i],
//				comb_stats[comb].num_lineages[i],	/*getEventTypeName(comb_stats[comb].event_types[i])*/ "UNIMPLEMENTED");
//	}
//	printf("====================================================================================\n");
//	fflush(stdout);
//}

// TODO - extract tests to different source file

void assertRootNumCoals(){
	int root = getPopIdByName(dataSetup.popTree, "root");

	int maxCoals = dataSetup.numLoci*(dataSetup.numSamples-1);

	int actualCoals = comb_stats[root].total.num_coals;
	for (int pop = 0 ; pop < dataSetup.popTree->numPops ; pop++){
		if (isLeaf(pop)){
			actualCoals += comb_stats[root].leaves[pop].below_comb.num_coals;
		}
	}
	if (maxCoals != actualCoals) { // TODO - rewrite this test
		printf("comb %s: Expected coalescence - %d. Actual coalescence - %d",
				"root", maxCoals, actualCoals);
		exit(-1);
	}
}

//void test_all(int comb){ // TODO - rewrite tests
//	test_validateCountCoals(comb);
//	test_validateCoalStats(comb);
//}
//void test_validateCountCoals(int comb){
//	int sumCoalsInComb = test_countCoalsEntireComb(comb);
//	int sumCoalsInClade = test_countCoalsEntireClade(comb);
//	if (sumCoalsInClade != sumCoalsInComb){
//		printf("Error in test_validateCountCoals(int comb): Comb had %d coals, however Clade had %d!", sumCoalsInComb, sumCoalsInClade);
//		exit(-1);
//	}
//}
//int test_countCoalsEntireComb(int comb){
//	int count = 0;
//	count += comb_stats[comb].num_coals_total;
//	count += test_countCoalsCombLeaves_rec(comb, comb);
//	return count;
//}
//int test_countCoalsCombLeaves_rec(int comb, int pop){
//	if (isLeaf(pop)){
//		return comb_stats[comb].leaves[pop].num_coals_below_comb;
//	} else {
//		int left_son = dataSetup.popTree->pops[pop]->sons[LEFT]->id;
//		int right_son = dataSetup.popTree->pops[pop]->sons[RIGHT]->id;
//		return (test_countCoalsCombLeaves_rec(comb, left_son) + test_countCoalsCombLeaves_rec(comb, right_son));
//	}
//}
//int test_countCoalsEntireClade(int clade){
//	if (isLeaf(clade)){
//		return genetree_stats_total.num_coals[clade];
//	} else {
//		int left_son  = dataSetup.popTree->pops[clade]->sons[ LEFT ]->id;
//		int right_son = dataSetup.popTree->pops[clade]->sons[ RIGHT]->id;
//		return genetree_stats_total.num_coals[clade] +
//				test_countCoalsEntireClade(left_son) +
//				test_countCoalsEntireClade(right_son);
//	}
//}
//void test_validateCoalStats(int comb){
//}
//
