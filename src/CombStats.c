
#include <stdlib.h>

#include "utils.h"
#include "MCMCcontrol.h"
#include "AlignmentProcessor.h"
#include "GenericTree.h"
#include "PopulationTree.h"
#include "LocusDataLikelihood.h"
#include "patch.h"

#include "CombStats.h"





// --- CONSTANTS --------------------------------------------------------------


#define COMB_PERCISION					0.0000000001



// --- GLOBAL DATA STRUCTURES -------------------------------------------------



// --- FUNCTION IMPLEMENTATIONS -----------------------------------------------


void computeCombStats() {

	return;

	initCombStats();
	computeCombNumCoals();
	computeCombCoalStats();
}

void initCombStats(){
	for (int comb = 0; comb < dataSetup.popTree->numPops; comb++) {
		initSpecificCombStats(comb);
	}
}


void initSpecificCombStats(int comb){
	comb_stats[comb].coal_stats_total   = 0.0;
	comb_stats[comb].mig_stats_total 	= 0.0;
	comb_stats[comb].num_coals_total 	= 0;
	comb_stats[comb].num_migs_total 	= 0;
	comb_stats[comb].num_events 		= 0;
	comb_stats[comb].debug_total_error  = 0;
}



void computeCombNumCoals(){
	computeCombNumCoals_rec(dataSetup.popTree->rootPop);
}

void computeCombNumCoals_rec(int pop){
	int leftSon, rightSon;

	if (isLeafPopulation(pop)){
		comb_stats[pop].num_coals_total = genetree_stats_total.num_coals[pop];
	} else{
		leftSon = dataSetup.popTree->pops[pop]->sons[LEFT]->id;
		rightSon = dataSetup.popTree->pops[pop]->sons[RIGHT]->id;

		computeCombNumCoals_rec(leftSon);
		computeCombNumCoals_rec(rightSon);

		comb_stats[pop].num_coals_total = genetree_stats_total.num_coals[pop] +
				comb_stats[leftSon].num_coals_total + comb_stats[rightSon].num_coals_total;
	}
}


void computeCombCoalStats(){
	for(int gen=0; gen<dataSetup.numLoci; gen++) {
		computeCombCoalStats_rec(dataSetup.popTree->rootPop, gen);
	}
}

void computeCombCoalStats_rec(int comb, int gen) {
	int leftSon, rightSon;

	if (isLeafPopulation(comb)){
		fillupLeafCombStats(comb, gen);
	} else{
		leftSon = dataSetup.popTree->pops[comb]->sons[LEFT]->id;
		rightSon = dataSetup.popTree->pops[comb]->sons[RIGHT]->id;

		computeCombCoalStats_rec(leftSon, gen);
		computeCombCoalStats_rec(rightSon, gen);

		fillupCombStats(comb, gen);

	}
//	debug_print_combstats(comb, gen); //CAUTION! uncommenting this increases running time tenfold! (because of IO)
}

void fillupLeafCombStats(int comb, int gen){
	appendPopToComb(comb, gen, 0); // since this is a leaf, starting point of the comb_stats arrays should be 0
}

void appendPopToComb(int comb, int gen, int startingPoint){
	int i = startingPoint;
	int event = event_chains[gen].first_event[comb];
	double combStartTime = dataSetup.popTree->pops[comb]->age;
	double eventAge = combStartTime;

	for ( ; event >= 0 ; i++, event = event_chains[gen].events[event].next){
		eventAge += event_chains[gen].events[event].elapsed_time;
		comb_stats[comb].sorted_ages[i] = eventAge;
		comb_stats[comb].elapsed_times[i] = event_chains[gen].events[event].elapsed_time;
		comb_stats[comb].num_lineages[i] = event_chains[gen].events[event].num_lineages;
		comb_stats[comb].event_types[i] = event_chains[gen].events[event].type;
	}
	comb_stats[comb].num_events = i;
//	comb_stats[comb].coal_stats_total += genetree_stats[gen].coal_stats[comb]; // TODO - decide which of these lines is better
	comb_stats[comb].coal_stats_total += getCoalStats(comb_stats[comb].elapsed_times + startingPoint, comb_stats[comb].num_lineages + startingPoint, comb_stats[comb].num_events - startingPoint);

#ifdef CHECKCOMB
	test_compareGphocsVsCombPopCoalStats(comb, gen, startingPoint);
#endif
}


void fillupCombStats(int comb, int gen){
	addChildenIntoCombStats(comb, gen);
	addCurrentPopIntoCombStats(comb, gen);
}

void addChildenIntoCombStats(int comb, int gen){
	mergeChildern(comb, gen);
	addChildrenCombStats(comb, gen);
}



void mergeChildern(int comb, int gen){
	int i, j = 0, k = 0;
	double leftAge, rightAge;
	int leftSon, rightSon;
	leftSon = dataSetup.popTree->pops[comb]->sons[LEFT]->id;
	rightSon = dataSetup.popTree->pops[comb]->sons[RIGHT]->id;
	int m = comb_stats[leftSon].num_events;
	int n = comb_stats[rightSon].num_events;

	for (i = 0 ; i < m + n; ) {
		if ( j < m && k < n) {
			leftAge = comb_stats[leftSon].sorted_ages[j];
			rightAge = comb_stats[rightSon].sorted_ages[k];
			if (leftAge < rightAge){
				comb_stats[comb].event_types[i] = comb_stats[leftSon].event_types[j];
				comb_stats[comb].sorted_ages[i] = comb_stats[leftSon].sorted_ages[j];
				comb_stats[comb].num_lineages[i] =
						comb_stats[leftSon].num_lineages[j] + comb_stats[rightSon].num_lineages[k];
				j++;
			}
			else{
				comb_stats[comb].event_types[i] = comb_stats[rightSon].event_types[k];
				comb_stats[comb].sorted_ages[i] = comb_stats[rightSon].sorted_ages[k];
				comb_stats[comb].num_lineages[i] =
						comb_stats[leftSon].num_lineages[j] + comb_stats[rightSon].num_lineages[k];
				k++;
			}
			i++;
		}
		else if (j == m) {
			for (; i < m + n ;) {
				comb_stats[comb].event_types[i] = comb_stats[rightSon].event_types[k];
				comb_stats[comb].sorted_ages[i] = comb_stats[rightSon].sorted_ages[k];
				comb_stats[comb].num_lineages[i] =
						comb_stats[leftSon].num_lineages[j-1] + comb_stats[rightSon].num_lineages[k];
				k++;
				i++;
			}
		}
		else {
			for (; i < m + n;) {
				comb_stats[comb].event_types[i] = comb_stats[leftSon].event_types[j];
				comb_stats[comb].sorted_ages[i] = comb_stats[leftSon].sorted_ages[j];
				comb_stats[comb].num_lineages[i] =
						comb_stats[leftSon].num_lineages[j] + comb_stats[rightSon].num_lineages[k-1];
				j++;
				i++;
			}
		}
	}
	for (i = 0 ; i < m + n; i++ ) {
		comb_stats[comb].elapsed_times[i] = (i==0) ? 0 :
				comb_stats[comb].sorted_ages[i] - comb_stats[comb].sorted_ages[i-1];
	}

	comb_stats[comb].num_events = m + n;
}

void addChildrenCombStats(int comb, int gen){
	comb_stats[comb].coal_stats_total +=
			getCoalStats(comb_stats[comb].elapsed_times, comb_stats[comb].num_lineages, comb_stats[comb].num_events);
}

void addCurrentPopIntoCombStats(int comb, int gen){
	appendPopToComb(comb, gen, comb_stats[comb].num_events); // start filling the comb_stats arrays from the last known event
}


double getCoalStats(double* elapsed_times, int* num_lineages, int size){
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





/***
 * This method is conjoined to "appendPopToComb".
 * It must run inside it, so don't move it.
 */
void test_compareGphocsVsCombPopCoalStats(int comb, int gen, int startingPoint){
	double gphocsValue = genetree_stats[gen].coal_stats[comb];
	double combValue = getCoalStats(comb_stats[comb].elapsed_times + startingPoint,
			comb_stats[comb].num_lineages + startingPoint, comb_stats[comb].num_events - startingPoint);

	double absoluteError = gphocsValue - combValue;
	double relativeError = absoluteError/combValue;
	if (COMB_PERCISION < fabs(relativeError)) {
		debug_print_combstats(comb, gen);
		comb_stats[comb].debug_total_error += absoluteError;
		fprintf(stderr, "Failed test_compareGphocsVsCombPopCoalStats:\n");
		fprintf(stderr, "comb:%s, gen:%d, startingPoint:%d\n", dataSetup.popTree->pops[comb]->name, gen, startingPoint);
		fprintf(stderr, "gphocs value  :%0.45f\n", gphocsValue);
		fprintf(stderr, "comb  value  :%0.45f\n", combValue);
		fprintf(stderr, "absolute error:%0.45f\n", absoluteError);
		fprintf(stderr, "relative error:%0.45f\n", relativeError);
		fprintf(stderr, "\nrandon seed   :%d\n", mcmcSetup.randomSeed);
		fflush(stderr);
		exit(-1);
	}
}
void test_validateRootCombVsFlatStats(){
	int rootPop = dataSetup.popTree->rootPop;
	double rootCombCoalStats = comb_stats[rootPop].coal_stats_total;
	double flatCoalStats = genetree_stats_flat.coal_stats_flat;
	double absoluteError = rootCombCoalStats - flatCoalStats;
	double relativeError = absoluteError/rootCombCoalStats;
	if (COMB_PERCISION*0.0001 < fabs(relativeError)) {
		fprintf(stderr, "Failed test_validateRootCombVsFlatStats:\n");
		fprintf(stderr, "flat		:%0.45f\n", flatCoalStats);
		fprintf(stderr, "rootComb	:%0.45f\n", rootCombCoalStats);
		fprintf(stderr, "absoluteError	:%0.45f\n", absoluteError);
		fprintf(stderr, "relativeError	:%0.45f\n", relativeError);
		fflush(stderr);
		exit(-1);
	}
}


#define TRUE 1 //TODO - where should these consts be?!
#define FALSE 0
int	isLeafPopulation(int pop){
	Population *population, *left_son, *right_son;

	population = dataSetup.popTree->pops[pop];

	left_son = population->sons[LEFT];
	right_son = population->sons[RIGHT];

	if (left_son || right_son){
		return FALSE;
	} else {
		return TRUE;
	}
}

