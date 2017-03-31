#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <string.h>

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
				finalizeCombCoalStats(comb);
//				debug_printCombGene(comb);
			}
		}
	}
	assertRootNumCoals();
//	assertRootCoalStats();
//	assertBottomCombs();
	assertCombLeaves();
}


void calculateSufficientStats(int comb, int gene){
	coalescence(comb, gene);
	migrations(comb, gene);
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
		updateEventVars(gene, &eventId, &elapsedTime, &eventType, &numLineages, &eventAge);

		if (eventAge <= combAge){
			if (eventType == COAL){
				belowCombLeafStats->num_coals++;
			}
			// you've got all the data you need so just update coal_stats -
			belowCombLeafStats->coal_stats += numLineages*(numLineages-1)*elapsedTime;
		}
		else if (eventAge > combAge){
			if (eventType == COAL){
				aboveCombLeafStats->num_coals++;
				combTotalStats->num_coals++;
			}
			// you're a tiny part of the comb so supply future methods the data they need -
			aboveCombLeafStats->sorted_ages[numEventsAboveComb] 	= eventAge;
			aboveCombLeafStats->elapsed_times[numEventsAboveComb] 	= elapsedTime;
			aboveCombLeafStats->num_lineages[numEventsAboveComb] 	= numLineages;
			aboveCombLeafStats->event_types[numEventsAboveComb] 	= eventType;

			numEventsAboveComb++;
		}
	}

	aboveCombLeafStats->num_events = numEventsAboveComb;
	combTotalStats->num_events += numEventsAboveComb;

	// TODO - this stat is not outputed. consider removing it performancewise
	aboveCombLeafStats->coal_stats += calculateCoalStats(aboveCombLeafStats->elapsed_times, aboveCombLeafStats->num_lineages, numEventsAboveComb);
}

void handleNonLeafCoals(int comb, int currentPop, int gene) {
	handleNonLeafNumCoals(comb, currentPop, gene);
	handleNonLeafCoalStats(comb, currentPop, gene);
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

	currentStats = getCombPopStats(comb, currentPop);
	leftStats = getCombPopStats(comb, leftSon);
	rightStats = getCombPopStats(comb, rightSon);


	int m = leftStats->num_events;
	int n = rightStats->num_events;

	for (i = 0 ; i < m + n; i++) {
		if ( j < m && k < n) {
			leftAge = leftStats->sorted_ages[j];
			rightAge = rightStats->sorted_ages[k];

			currentStats->num_lineages[i] = leftStats->num_lineages[j] + rightStats->num_lineages[k];
			if (leftAge < rightAge){
				currentStats->event_types[i]  = leftStats->event_types[j];
				currentStats->sorted_ages[i]  = leftStats->sorted_ages[j];
				j++;
			}
			else{
				currentStats->event_types[i]  = rightStats->event_types[k];
				currentStats->sorted_ages[i]  = rightStats->sorted_ages[k];
				k++;
			}
		}
		else if (j == m) {
				currentStats->event_types[i]  = rightStats->event_types[k];
				currentStats->sorted_ages[i]  = rightStats->sorted_ages[k];
				currentStats->num_lineages[i] = leftStats->num_lineages[j-1] + rightStats->num_lineages[k];
				k++;
		}
		else if (k == n){
				currentStats->event_types[i]  = leftStats->event_types[j];
				currentStats->sorted_ages[i]  = leftStats->sorted_ages[j];
				currentStats->num_lineages[i] = leftStats->num_lineages[j] + rightStats->num_lineages[k-1];
				j++;
		}
	}
	for (i = 0 ; i < m + n; i++ ) {
		currentStats->elapsed_times[i] = (i==0) ? 0 :
				currentStats->sorted_ages[i] - currentStats->sorted_ages[i-1];
	}

	currentStats->num_events = m + n;
}

void appendCurrent(int comb, int currentPop, int gene){
	Stats *currentStats = getCombPopStats(comb, currentPop);

	int startingPoint = currentStats->num_events;
	int i = startingPoint;
	int event = event_chains[gene].first_event[currentPop];
	double startTime = dataSetup.popTree->pops[currentPop]->age; // do I need to start from pop age or from last event age (or are they equal)?
	double eventAge = startTime;

	for ( ; event >= 0 ; i++, event = event_chains[gene].events[event].next){
		eventAge += event_chains[gene].events[event].elapsed_time;
		currentStats->sorted_ages[i] = eventAge;
		currentStats->elapsed_times[i] = event_chains[gene].events[event].elapsed_time;
		currentStats->num_lineages[i] = event_chains[gene].events[event].num_lineages;
		currentStats->event_types[i] = event_chains[gene].events[event].type;
	}
	currentStats->num_events = i;
}

void finalizeCombCoalStats(int comb){
	double* elapsedTimes = comb_stats[comb].clades[comb].elapsed_times;
	int* numLineages = comb_stats[comb].clades[comb].num_lineages;
	int size = comb_stats[comb].clades[comb].num_events;
	comb_stats[comb].total.coal_stats += calculateCoalStats(elapsedTimes, numLineages, size);
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


void migrations(int comb, int gene){
	for (int mig = 0 ; mig < dataSetup.popTree->numMigBands ; mig++){
		if (isMigOfComb(mig, comb)){
			if (isLeafMigBand(mig, comb)){
				handleLeafMigStats(comb, mig, gene);
			}
			if (isMigBandInternal(mig, comb)) {
				// ignore internal migbands. their stats aren't used
			}
			if (isMigBandExternal(mig, comb)){
				// handle migbands from outside to the comb body
			}
		}
	}
}

void handleLeafMigStats(int comb, int mig, int gene){
	double elapsedTime, eventAge = 0.0;
	int numLineages, eventType;
	double combAge = comb_stats[comb].age;

	int targetLeaf = getTargetPop(mig);

	MigStats* migLeafStats = &comb_stats[comb].leafMigs[targetLeaf];

	int eventId = event_chains[gene].first_event[targetLeaf];
	while (eventAge <= combAge && eventId > 0) {
		updateEventVars(gene, &eventId, &elapsedTime, &eventType, &numLineages, &eventAge);

		migLeafStats->mig_stats += numLineages*elapsedTime;
		if (eventType == IN_MIG) migLeafStats->num_migs++;
	}
}

void updateEventVars(int gene, int* eventId, double*elapsedTime, int* eventType, int* numLineages, double* eventAge) {
	*elapsedTime = event_chains[gene].events[*eventId].elapsed_time;
	*eventType = event_chains[gene].events[*eventId].type;
	*numLineages = event_chains[gene].events[*eventId].num_lineages;
	*eventAge += *elapsedTime;

	*eventId = event_chains[gene].events[*eventId].next;
}










//TODO - find places for these functions
double getCombAge(int comb){
	if (isLeaf(comb)){
		return DBL_MAX;
	} else if (areChildrenLeaves(comb)){
		return dataSetup.popTree->pops[comb]->age;
	} else if (isFeasibleComb(comb)){
		double left_min = getCombAge(dataSetup.popTree->pops[comb]->sons[LEFT]->id);
		double right_min = getCombAge(dataSetup.popTree->pops[comb]->sons[RIGHT]->id);
		return fmin(left_min, right_min);
	} else {
		printf("ERROR: bug in combAge algorithm. Should not reach here!");
		exit(-1);
		return DBL_MAX;
	}
}
int	isLeaf(int pop){
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
	} else if (areChildrenLeaves(pop)){ // a population whos two children are leaves is considered a "trivial comb" and is usually not interesting (except for algorithm test purposes)
		return FALSE; // Set TRUE if you want to run "assertBottomCombs()" test // TODO - always set back to FALSE after testing
	} else {
		return TRUE;
	}
}
int isAncestralTo(int father, int son){
	return dataSetup.popTree->pops[father]->isAncestralTo[son];
}
int isMigOfComb(int mig, int comb){ // if migrations flow into the comb
	int target = getTargetPop(mig);
	return isAncestralTo(comb, target) || target == comb;
}
int isMigBandExternal(int mig, int comb){
	int source = getSourcePop(mig);
	int target = getTargetPop(mig);
	return isAncestralTo(comb, target)  && !isAncestralTo(comb, source);
}
int isMigBandInternal(int mig, int comb){
	int source = getSourcePop(mig);
	int target = getTargetPop(mig);

	return isAncestralTo(comb, source) && isAncestralTo(comb, target) && (!isLeaf(source) || !isLeaf(target));
}

int isLeafMigBand(int mig, int comb){
	int target = getTargetPop(mig);
	int source = getSourcePop(mig);
	return isLeaf(target) && (isLeaf(source) || !isAncestralTo(comb, source));
}

char* getEventTypeName(int eventType){
	switch(eventType){
		case COAL:
			return "COAL";
		case IN_MIG:
			return "IN_MIG";
		case OUT_MIG:
			return "OUT_MIG";
		case MIG_BAND_START:
			return "MIG_START";
		case MIG_BAND_END:
			return "MIG_END";
		case SAMPLES_START:
			return "SAM_START";
		case END_CHAIN:
			return "END_CHAIN";
		case DUMMY:
			return "DUMMY";
		default:
			return "UNDEFINED";
	}
}
int getSon(int pop, int SON){
	return dataSetup.popTree->pops[pop]->sons[SON]->id;
}

int getSourcePop(int mig){
	return dataSetup.popTree->migBands[mig].sourcePop;
}
int getTargetPop(int mig){
	return dataSetup.popTree->migBands[mig].targetPop;
}
char* getPopName(int pop){
	return dataSetup.popTree->pops[pop]->name;
}
char* getMigName(int mig){
	char* sourceName = getPopName(getSourcePop(mig));
	char* targetName = getPopName(getTargetPop(mig));
	return concat(sourceName, targetName);
}
char* concat(const char *s1, const char *s2){ // TODO - THIS SHOULD NOT BE USED  IN PRODUCTION SINCE THE MEMORY ISN'T RELEASED
    char *result = malloc(strlen(s1)+strlen(s2)+3);//+3 for the zero-terminator and arrow
    strcpy(result, s1);
    strcat(result, "->");
    strcat(result, s2);
    return result;
}


Stats* getCombPopStats(int comb, int pop){
	if (isLeaf(pop)){
		return &comb_stats[comb].leaves[pop].above_comb;
	} else {
		return &comb_stats[comb].clades[pop];
	}
}

void initCombStats(){
	initPopStats();
	initMigStats();
}
void initPopStats() {
	for (int comb = 0; comb < dataSetup.popTree->numPops; comb++) {
		if (isFeasibleComb(comb)) {
			initStats(&comb_stats[comb].total);
			comb_stats[comb].age = getCombAge(comb);
			for (int pop = 0; pop < dataSetup.popTree->numPops; pop++) {
				if (isLeaf(pop)) {
					initStats(&comb_stats[comb].leaves[pop].above_comb);
					initStats(&comb_stats[comb].leaves[pop].below_comb);
				} else {
					initStats(&comb_stats[comb].clades[pop]);
				}
			}
		}
	}
}
void initStats(Stats* stats){ // TODO - rename signature to include "pop"
	stats->coal_stats = 0.0;
	stats->num_coals  = 0;
	stats->num_events = 0;
}
void initMigStats() {
	for (int comb = 0; comb < dataSetup.popTree->numPops; comb++) {
		if (isFeasibleComb(comb)) {
			for (int mig = 0; mig < dataSetup.popTree->numMigBands; mig++) {
				if (isMigOfComb(mig, comb)) {
					comb_stats[comb].leafMigs[mig].mig_stats = 0.0;
					comb_stats[comb].leafMigs[mig].num_migs = 0;
				}
			}
		}
	}
}
void allocateCombMem(){
	comb_stats=malloc(dataSetup.popTree->numPops*sizeof(struct COMB_STATS));

	allocatePopsMem();
	allocateMigBandsMem();
}
void allocatePopsMem() {
	for (int comb = 0; comb < dataSetup.popTree->numPops; comb++) {
		allocateStats(&comb_stats[comb].total);
		comb_stats[comb].leaves = malloc(
				dataSetup.popTree->numCurPops * sizeof(LeafStats));
		comb_stats[comb].clades = malloc(
				dataSetup.popTree->numCurPops * sizeof(Stats));
		for (int pop = 0; pop < dataSetup.popTree->numPops; pop++) {
			if (isLeaf(pop)) {
				allocateStats(&comb_stats[comb].leaves[pop].below_comb);
				allocateStats(&comb_stats[comb].leaves[pop].above_comb);
			} else {
				allocateStats(&comb_stats[comb].clades[pop]);
			}
		}
		if (comb_stats == NULL) {
			// TODO - add memory allocation test for all of comb_stats
			fprintf(stderr, "\nError: Out Of Memory comb_stats\n");
			exit(-1);
		}
	}
}
void allocateStats(Stats* stats){ // TODO - rename signature to include "pop"
	int max_events = 2*dataSetup.numSamples+ 4*MAX_MIGS + 3*dataSetup.popTree->numMigBands + dataSetup.popTree->numPops + 10;
	stats->sorted_ages   = (double*)malloc(max_events*sizeof(double));
	stats->elapsed_times = (double*)malloc(max_events*sizeof(double));
	stats->num_lineages  = (int*)malloc(max_events*sizeof(int));
	stats->event_types   = (int*)malloc(max_events*sizeof(int));
}
void allocateMigBandsMem() {
	int maxMigBands = dataSetup.popTree->numMigBands;
	for (int comb = 0; comb < dataSetup.popTree->numPops; comb++) {
		if (isFeasibleComb(comb)) {
			comb_stats[comb].leafMigs = malloc(maxMigBands * sizeof(MigStats));
		}
	}
}
void freeCombMem(){ // TODO - implement
}


// TODO - extract tests to different source file
double COMB_RELATIVE_PERCISION = 	0.000000000001;
void debug_printCombGene(int comb){
	char* combName = dataSetup.popTree->pops[comb]->name;
	double combAge = comb_stats[comb].age;
	double* elapsedTimes = comb_stats[comb].clades[comb].elapsed_times;
	int* numLineages = comb_stats[comb].clades[comb].num_lineages;
	int* eventTypes = comb_stats[comb].clades[comb].event_types;
	int size = comb_stats[comb].clades[comb].num_events;
	double currentAge = combAge + elapsedTimes[0];

	printf("\ncomb:%s, size:%d\n",combName, size);
	for (int i = 0 ; i< size ; i++){
		printf("type:%s\tnumLins:%d\tage:%0.35f\n",getEventTypeName(eventTypes[i]), numLineages[i], currentAge);
		currentAge += elapsedTimes[i];
	}
}
void assertRootNumCoals(){
	int root = getPopIdByName(dataSetup.popTree, "root");

	int maxCoals = dataSetup.numLoci*(dataSetup.numSamples-1);

	int actualCoals = comb_stats[root].total.num_coals;
	for (int pop = 0 ; pop < dataSetup.popTree->numPops ; pop++){
		if (isLeaf(pop)){
			actualCoals += comb_stats[root].leaves[pop].below_comb.num_coals;
		}
	}

	if (actualCoals != maxCoals) {
		printf("comb %s: Expected coalescence - %d. Actual coalescence - %d",
				"root", maxCoals, actualCoals);
		exit(-1);
	}
}
/**
 * Tests whether root comb with mocked comb_age:=0.0 gives the same coal_stats as a flat model.
 *  To enable this test, you must hard-code getCombAge() to return 0.0 always
 */
void assertRootCoalStats(){
	int root = getPopIdByName(dataSetup.popTree, "root");
	double actualCoalStats = comb_stats[root].total.coal_stats;
	double expectedCoalStats = genetree_stats_flat.coal_stats_flat;

	double error = fabs(actualCoalStats - expectedCoalStats);
	double relativeError = error/expectedCoalStats;

	if (relativeError > COMB_RELATIVE_PERCISION){
		printf("Error while checking root coal_stats:\nExpected:%0.35f\tActual:%0.35f\tRelative Error:%0.35f",
				expectedCoalStats, actualCoalStats, relativeError);
		exit(-1);
	}
}
/**
 * Compares trivial combs (pops directly above two leaves) with regular pops.
 * Can only be used when isFeasibleComb() allows trivial combs (and it usually shouldn't)
 */
void assertBottomCombs(){
	for (int comb = 0 ; comb < dataSetup.popTree->numPops ; comb++){
		if (isFeasibleComb(comb) && areChildrenLeaves(comb)){
			assertBottomCombsNumCoals(comb);
			assertBottomCombsCoalStats(comb);
		}
	}
}
void assertBottomCombsNumCoals(int comb){
	int expected = genetree_stats_total.num_coals[comb];
	int actual = comb_stats[comb].total.num_coals;

	if (expected != actual){
		printf("Error while checking comb %s num_coals:\nExpected num_coals %d. actual is %d",
				dataSetup.popTree->pops[comb]->name, expected, actual);
		exit(-1);
	}
}
void assertBottomCombsCoalStats(int comb){
	double expected = genetree_stats_total.coal_stats[comb];
	double actual = comb_stats[comb].total.coal_stats;
	double error = fabs(actual - expected);
	double relativeError = error/expected;
	if (relativeError > COMB_RELATIVE_PERCISION){
		printf("Error while checking comb %s coal_stats:\nExpected:%0.35f\tActual:%0.35f\tRelative Error:%0.35f\tAbsolute Error:%0.35f",
				dataSetup.popTree->pops[comb]->name,
				expected, actual, relativeError, error);
		exit(-1);
	}
}
void assertCombLeaves(){
	for (int comb = 0; comb < dataSetup.popTree->numPops ; comb++){
		if (isFeasibleComb(comb)){
			for (int leaf = 0; leaf < dataSetup.popTree->numPops ; leaf++){
				if (isLeaf(leaf) && isAncestralTo(comb, leaf)){
					assertCombLeafNumCoals(comb, leaf);
					assertCombLeafCoalStats(comb, leaf);
				}
			}
		}
	}
}
void assertCombLeafNumCoals(int comb, int leaf){
	int expectedNumCoals = genetree_stats_total.num_coals[leaf];
	int actualNumCoals = comb_stats[comb].leaves[leaf].above_comb.num_coals
			+ comb_stats[comb].leaves[leaf].below_comb.num_coals;
	if (expectedNumCoals != actualNumCoals){
		printf("Error while checking leaf %s num_coals:\nExpected num_coals %d. actual is %d",
				dataSetup.popTree->pops[leaf]->name, expectedNumCoals, actualNumCoals);
		exit(-1);
	}
}
void assertCombLeafCoalStats(int comb, int leaf){
	double expectedCoalStats = genetree_stats_total.coal_stats[leaf];
	double actualCoalStats = comb_stats[comb].leaves[leaf].above_comb.coal_stats
			+ comb_stats[comb].leaves[leaf].below_comb.coal_stats;
	double error = fabs(actualCoalStats - expectedCoalStats);
	double relativeError = error/expectedCoalStats;
	if (relativeError > COMB_RELATIVE_PERCISION){
		printf("Error while checking leaf %s coal_stats:\nExpected:%0.35f\tActual:%0.35f\tRelative Error:%0.35f\tAbsolute Error:%0.35f",
				dataSetup.popTree->pops[leaf]->name,
				expectedCoalStats, actualCoalStats, relativeError, error);
		exit(-1);
	}
}
