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
#include "DataLayer.h"
#include "MemoryMng.h"
#include "CombStats.h"
#include "McRefCommon.h"

COMB_STATS* comb_stats;

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
	//	assertRootNumCoals();
	//	assertRootCoalStats();
	//	assertBottomCombs();
	assertCombLeaves();
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

	double elapsedTime = 0.0, eventAge = 0.0, previousAge = 0.0;
	int numLineages, eventType, eventId, nextId;
	double combAge = comb_stats[comb].age;

	Stats* belowCombLeafStats = &comb_stats[comb].leaves[leaf].below_comb;
	Stats* aboveCombLeafStats = &comb_stats[comb].leaves[leaf].above_comb;
	Stats* combTotalStats = &comb_stats[comb].total;

	setupFirstEventVars(gene, leaf, &eventId, &nextId, &elapsedTime, &eventType, &numLineages);


	handleBelowCombAge(gene, combAge, &eventId, &nextId, &elapsedTime, &eventType, &numLineages, &eventAge, &previousAge, belowCombLeafStats);
	handleAboveCombAge(gene, combAge, &eventId, &nextId, &elapsedTime, &eventType, &numLineages, &eventAge, &previousAge, belowCombLeafStats, aboveCombLeafStats, combTotalStats);

}

void handleBelowCombAge(int gene, double combAge, int* eventId, int* nextId, double* elapsedTime, int* eventType, int* numLineages, double* eventAge, double* previousAge, Stats* belowCombLeafStats) {
	while (*eventAge < combAge && *nextId >=0){
		if (*eventType == COAL) belowCombLeafStats->num_coals++;
		belowCombLeafStats->coal_stats += (*numLineages) * (*numLineages-1) * (*elapsedTime);
		incrementEventVars(gene, eventId, nextId, elapsedTime, eventType, numLineages, eventAge, previousAge);
	}
}

void handleAboveCombAge(int gene, double combAge, int* eventId, int* nextId, double* elapsedTime, int* eventType, int* numLineages, double* eventAge, double* previousAge, Stats* belowCombLeafStats, Stats* aboveCombLeafStats, Stats* combTotalStats){
	handleFirstEventAboveCombAge  (gene, combAge, eventId, nextId, elapsedTime, eventType, numLineages, eventAge, previousAge, belowCombLeafStats, aboveCombLeafStats);
	handleRestOfEventsAboveCombAge(gene,          eventId, nextId, elapsedTime, eventType, numLineages, eventAge, previousAge,                     aboveCombLeafStats, combTotalStats);
}

void handleFirstEventAboveCombAge(int gene, double combAge, int* eventId, int* nextId, double*elapsedTime, int* eventType, int* numLineages, double* eventAge, double* previousAge, Stats* belowCombLeafStats, Stats* aboveCombLeafStats){

	double pseudoEventAge, pseudoElapsedTimeBelow, pseudoElapsedTimeAbove;

	if (*eventAge == combAge) {
		if (*eventType == COAL) belowCombLeafStats->num_coals++;
		belowCombLeafStats->coal_stats += (*numLineages)*(*numLineages-1)*(*elapsedTime);

		pseudoElapsedTimeAbove	= 0.0;

	} else {
		pseudoEventAge = combAge;
		pseudoElapsedTimeBelow = pseudoEventAge - *previousAge;
		belowCombLeafStats->coal_stats += (*numLineages)*(*numLineages-1)*(pseudoElapsedTimeBelow);

		pseudoElapsedTimeAbove = *eventAge - pseudoEventAge;
	}

	aboveCombLeafStats->elapsed_times[0]	= pseudoElapsedTimeAbove;
	aboveCombLeafStats->sorted_ages[0] 		= *eventAge;
	aboveCombLeafStats->num_lineages[0] 	= *numLineages;
	aboveCombLeafStats->event_types[0] 		= *eventType;
	aboveCombLeafStats->event_ids[0] 	    = *eventId;
	aboveCombLeafStats->num_events 			= 1;

	incrementEventVars(gene, eventId, nextId, elapsedTime, eventType, numLineages, eventAge, previousAge);
}

void handleRestOfEventsAboveCombAge(int gene, int* eventId, int* nextId, double*elapsedTime, int* eventType, int* numLineages, double* eventAge, double* previousAge, Stats* aboveCombLeafStats, Stats* combTotalStats){

	int numEventsAboveComb = aboveCombLeafStats->num_events;

	while (*eventId >= 0) {

		if (*eventType == COAL){
			aboveCombLeafStats->num_coals++;
			combTotalStats->num_coals++;
		}
		// you're a tiny part of the comb so supply future methods the data they need -
		aboveCombLeafStats->sorted_ages[numEventsAboveComb] 	= *eventAge;
		aboveCombLeafStats->elapsed_times[numEventsAboveComb] 	= *elapsedTime;
		aboveCombLeafStats->num_lineages[numEventsAboveComb] 	= *numLineages;
		aboveCombLeafStats->event_types[numEventsAboveComb] 	= *eventType;
		aboveCombLeafStats->event_ids[numEventsAboveComb] 	    = *eventId;

		numEventsAboveComb++;
		incrementEventVars(gene, eventId, nextId, elapsedTime, eventType, numLineages, eventAge, previousAge);
	}
	aboveCombLeafStats->num_events = numEventsAboveComb;
	combTotalStats->num_events += numEventsAboveComb;

	// TODO - this stat is not outputed. consider removing it performancewise
	aboveCombLeafStats->coal_stats += calculateCoalStats(aboveCombLeafStats->elapsed_times, aboveCombLeafStats->num_lineages, numEventsAboveComb);
}

void setupFirstEventVars(int gene, int leaf, int* eventId, int *nextId,
		double* elapsedTime, int* eventType, int* numLineages) {

	*eventId = event_chains[gene].first_event[leaf];
	*nextId = event_chains[gene].events[*eventId].getNextIdx();
	*elapsedTime = event_chains[gene].events[*eventId].getElapsedTime(); // TODO - assert this is zero
	*eventType = event_chains[gene].events[*eventId].getType();
	*numLineages = event_chains[gene].events[*eventId].getNumLineages();
}

void incrementEventVars(int gene, int* eventId, int *nextId,
		double*elapsedTime, int* eventType, int* numLineages,
		double* eventAge, double* previousAge) {

	*eventId = event_chains[gene].events[*eventId].getNextIdx();
	*nextId = event_chains[gene].events[*eventId].getNextIdx();
	*elapsedTime = event_chains[gene].events[*eventId].getElapsedTime();
	*eventType = event_chains[gene].events[*eventId].getType();
	*numLineages = event_chains[gene].events[*eventId].getNumLineages();
	*previousAge = *eventAge;
	*eventAge += *elapsedTime;
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
	Stats *leftStats, *rightStats, *currentStats;

	currentStats = getCombPopStats(comb, currentPop);
	leftStats = getCombPopStats(comb, getSon(currentPop, LEFT));
	rightStats = getCombPopStats(comb, getSon(currentPop, RIGHT));


	int m = leftStats->num_events;
	int n = rightStats->num_events;

	for (i = 0 ; i < m + n; i++) {
		if (m == 0){ // there are no events in left son, so choose the event from the right son
			copyStaticEventStats(rightStats, k, currentStats, i);
			k++;
		}
		else if (n == 0){ // there are no events in left son, so choose the event from the left son
			copyStaticEventStats(leftStats, j, currentStats, i);
			j++;
		}
		else if ( j < m && k < n) { // both sons have more events
			currentStats->num_lineages[i] = leftStats->num_lineages[j] + rightStats->num_lineages[k];

			leftAge = leftStats->sorted_ages[j];
			rightAge = rightStats->sorted_ages[k];
			if (leftAge < rightAge){ // since left event is more recent, choose it
				copyStaticEventStats(leftStats, j, currentStats, i);
				j++;
			}
			else{  // since right event is more recent, choose it
				copyStaticEventStats(rightStats, k, currentStats, i);
				k++;
			}
		}
		else if (j == m) { // all events in left son were handled
			copyStaticEventStats(rightStats, k, currentStats, i);
			currentStats->num_lineages[i] = leftStats->num_lineages[j-1] + rightStats->num_lineages[k];
			k++;
		}
		else if (k == n){ // all events in right son were handled
			copyStaticEventStats(leftStats, j, currentStats, i);
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
void copyStaticEventStats(Stats* sourceStats, int n, Stats* targetStats, int m) {
	targetStats->event_ids[m]   = sourceStats->event_ids[n];
	targetStats->event_types[m] = sourceStats->event_types[n];
	targetStats->sorted_ages[m] = sourceStats->sorted_ages[n];
}

void appendCurrent(int comb, int currentPop, int gene){
	Stats *currentStats = getCombPopStats(comb, currentPop);

	int startingPoint = currentStats->num_events;
	int i = startingPoint;
	int eventId = event_chains[gene].first_event[currentPop];
	double startTime = dataSetup.popTree->pops[currentPop]->age; // do I need to start from pop age or from last eventId age (or are they equal)?
	double eventAge = startTime;

	for ( ; eventId >= 0 ; i++, eventId = event_chains[gene].events[eventId].getNextIdx()){
		eventAge += event_chains[gene].events[eventId].getElapsedTime();
		currentStats->sorted_ages[i]   = eventAge;
		currentStats->elapsed_times[i] = event_chains[gene].events[eventId].getElapsedTime();
		currentStats->num_lineages[i]  = event_chains[gene].events[eventId].getNumLineages();
		currentStats->event_types[i]   = event_chains[gene].events[eventId].getType();
		currentStats->event_ids[i]     = eventId;
	}
	currentStats->num_events = i;
}

void finalizeCombCoalStats(int comb){
	double* elapsedTimes = comb_stats[comb].combs[comb].elapsed_times;
	int* numLineages = comb_stats[comb].combs[comb].num_lineages;
	int size = comb_stats[comb].combs[comb].num_events;
	comb_stats[comb].total.coal_stats += calculateCoalStats(elapsedTimes, numLineages, size);
}



void migrations(int comb, int gene){
	for (int migband = 0 ; migband < dataSetup.popTree->numMigBands ; migband++){
		if (isMigOfComb(migband, comb)){
			if (isLeafMigBand(migband, comb)){
				handleLeafMigStats(comb, migband, gene);
			}
			if (isMigBandInternal(migband, comb)) {
				// ignore internal migbands. their stats aren't used
			}
			if (isMigBandExternal(migband, comb)){
				// handleExternalMigStats(comb, mig, gene);
			}
		}
	}
}

void handleLeafMigStats(int comb, int mig, int gene){
	double elapsedTime, eventAge = 0.0, previousAge;
	double combAge = comb_stats[comb].age;
	int numLineages;
	int targetLeaf = getTargetPop(mig);
	int eventId = event_chains[gene].first_event[targetLeaf];
	int nextId = event_chains[gene].events[eventId].getNextIdx();
	int eventType = event_chains[gene].events[eventId].getType();
	MigStats* migLeafStats = &comb_stats[comb].leafMigs[targetLeaf];

	fastFwdPastMigBandStart(gene, &eventId, &nextId, &elapsedTime, &eventType, &numLineages, &eventAge, &previousAge);

	while (eventAge <= combAge && eventId > 0) {
		updateLeafMigStats(numLineages, elapsedTime, eventType, migLeafStats);

		incrementEventVars(gene, &eventId, &nextId, &elapsedTime, &eventType, &numLineages, &eventAge, &previousAge);

		if (eventType == MIG_BAND_END) {
			updateLeafMigStats(numLineages, elapsedTime, eventType, migLeafStats);
			break;
		}
	}
}

void fastFwdPastMigBandStart(int gene, int* eventId, int* nextId, double*elapsedTime, int* eventType, int* numLineages, double* eventAge, double* previousAge) {
	while (*eventType != MIG_BAND_START && *eventId > 0){
		incrementEventVars(gene, eventId, nextId, elapsedTime, eventType, numLineages, eventAge, previousAge);
	}
	incrementEventVars(gene, eventId, nextId, elapsedTime, eventType, numLineages, eventAge, previousAge);
}

void updateLeafMigStats(int numLineages, double elapsedTime, int eventType, MigStats* migLeafStats) {
	migLeafStats->mig_stats += numLineages * elapsedTime;
	if (eventType == IN_MIG) migLeafStats->num_migs++;// TODO - make sure the event is on the right migband
}


void handleExternalMigStats(int comb, int mig, int gene){
	int* eventIds = comb_stats[comb].combs[comb].event_ids;
	int* eventTypes = comb_stats[comb].combs[comb].event_types;
	double* sortedAges = comb_stats[comb].combs[comb].sorted_ages;
	int numEvents = comb_stats[comb].combs[comb].num_events;

	for (int i = 0 ; i < numEvents ; i++){
		printf("%s\n", getEventTypeName(eventTypes[i]));
		if (eventTypes[i] == IN_MIG || eventTypes[i] == OUT_MIG){
			int eventId = eventIds[i];
			double migAge = genetree_migs[gene].mignodes[eventId].age;
			int migSourceId = genetree_migs[gene].mignodes[eventId].source_event;
			int migTargetId = genetree_migs[gene].mignodes[eventId].target_event;
			double eventAge = sortedAges[i];

			printf("eventAge:%f, migAge:%f", eventAge , migAge);
			printf("eventId:%d, migSourceId:%d, migTargetId:%d", eventId, migSourceId , migTargetId);
		}
	}
}










//TODO - find places for these functions
Stats* getCombPopStats(int comb, int pop){
	if (isLeaf(pop)){
		return &comb_stats[comb].leaves[pop].above_comb;
	} else {
		return &comb_stats[comb].combs[pop];
	}
}

void initCombStats(){
	initPopStats();
	initMigStats();
}
void initPopStats() {
	for (int comb = 0; comb < dataSetup.popTree->numPops; comb++) {
		if (!isLeaf(comb)) {
			initStats(&comb_stats[comb].total);
			comb_stats[comb].age = getCombAge(comb);
			for (int pop = 0; pop < dataSetup.popTree->numPops; pop++) {
				if (isLeaf(pop)) {
					initStats(&comb_stats[comb].leaves[pop].above_comb);
					initStats(&comb_stats[comb].leaves[pop].below_comb);
				} else {
					initStats(&comb_stats[comb].combs[pop]);
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
	comb_stats= (COMB_STATS*) malloc(dataSetup.popTree->numPops*sizeof(COMB_STATS));

	allocatePopsMem();
	allocateMigBandsMem();
}
void allocatePopsMem() {
	for (int comb = 0; comb < dataSetup.popTree->numPops; comb++) {
		allocateStats(&comb_stats[comb].total);
		comb_stats[comb].leaves = (LeafStats*) malloc(dataSetup.popTree->numPops*sizeof(LeafStats));
		comb_stats[comb].combs  = (Stats*)     malloc(dataSetup.popTree->numPops*sizeof(Stats));
		for (int pop = 0 ; pop < dataSetup.popTree->numPops ; pop++){
			if (isLeaf(pop)){
				allocateStats(&comb_stats[comb].leaves[pop].below_comb);
				allocateStats(&comb_stats[comb].leaves[pop].above_comb);
			} else {
				allocateStats(&comb_stats[comb].combs[pop]);
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
	stats->sorted_ages   = (double*) malloc(max_events*sizeof(double));
	stats->elapsed_times = (double*) malloc(max_events*sizeof(double));
	stats->num_lineages  = (int*) malloc(max_events*sizeof(int));
	stats->event_types   = (int*) malloc(max_events*sizeof(int));
	stats->event_ids     = (int*) malloc(max_events*sizeof(int));
}
void allocateMigBandsMem() {
	int maxMigBands = dataSetup.popTree->numMigBands;
	for (int comb = 0; comb < dataSetup.popTree->numPops; comb++) {
		if (isFeasibleComb(comb)) {
			comb_stats[comb].leafMigs = (MigStats*) malloc(maxMigBands * sizeof(MigStats));
		}
	}
}


int isFeasibleComb(int pop){
	return !isLeaf(pop);
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




void freeCombMem(){ // TODO - implement
}


// TODO - extract tests to different source file
double COMB_RELATIVE_PERCISION = 	0.0000000000001;

void debug_printCombGene(int comb){
	char* combName = dataSetup.popTree->pops[comb]->name;
	double combAge = comb_stats[comb].age;
	double* elapsedTimes = comb_stats[comb].combs[comb].elapsed_times;
	int* numLineages = comb_stats[comb].combs[comb].num_lineages;
	int* eventTypes = comb_stats[comb].combs[comb].event_types;
	int size = comb_stats[comb].combs[comb].num_events;
	double currentAge = combAge + elapsedTimes[0];

	printf("\ncomb:%s, size:%d\n",combName, size);
	for (int i = 0 ; i< size ; i++){
		printf("type:%s\tnumLins:%d\tage:%0.35f\n",getEventTypeName(eventTypes[i]), numLineages[i], currentAge);
		currentAge += elapsedTimes[i];
	}
}
void assertRootNumCoals(){ // this test only works when I hard-codedly force combAge to Zero
	int root = getPopIdByName(dataSetup.popTree, "root");

	int maxCoals = dataSetup.numLoci*(dataSetup.numSamples-1);

	int actualCoals = comb_stats[root].total.num_coals;
	for (int pop = 0 ; pop < dataSetup.popTree->numPops ; pop++){
		if (isLeaf(pop)){
			actualCoals += comb_stats[root].leaves[pop].below_comb.num_coals;
		}
	}

	if (actualCoals != maxCoals) {
		printf("comb %s: Expected coalescence events - %d. Actual coalescence events - %d",
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
		printf("\nError while checking comb %s num_coals:\nExpected num_coals %d. actual is %d",
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
		printf("\nError while checking comb %s coal_stats:\nExpected:%0.35f\tActual:%0.35f\tRelative Error:%0.35f\tAbsolute Error:%0.35f",
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
					//					assertCombLeafNumCoals(comb, leaf);
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
		printf("\nError while checking leaf %s num_coals:\nExpected num_coals %d. actual is %d",
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
		printf("\nError while checking leaf %s coal_stats:\nExpected:%0.35f\tActual:%0.35f\tRelative Error:%0.35f\tAbsolute Error:%0.35f",
				dataSetup.popTree->pops[leaf]->name,
				expectedCoalStats, actualCoalStats, relativeError, error);
		exit(-1);
	}
}

