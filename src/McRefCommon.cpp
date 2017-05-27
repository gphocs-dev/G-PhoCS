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
#include "DataLayer.h"

#include "McRefCommon.h"
#include "patch.h"


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
const char* getEventTypeName(int eventType)
{
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
//char* getMigName(int mig){
//	char* sourceName = getPopName(getSourcePop(mig));
//	char* targetName = getPopName(getTargetPop(mig));
//	return concat(sourceName, targetName);
//}
//char* concat(const char *s1, const char *s2){ // TODO - THIS SHOULD NOT BE USED  IN PRODUCTION SINCE THE MEMORY ISN'T RELEASED
//    char *result = malloc(strlen(s1)+strlen(s2)+3);//+3 for the zero-terminator and arrow
//    strcpy(result, s1);
//    strcat(result, "->");
//    strcat(result, s2);
//    return result;
//}


