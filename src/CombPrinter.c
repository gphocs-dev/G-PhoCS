#include "CombPrinter.h"
#include "CombStats.h"
#include "MCMCcontrol.h"
#include "GPhoCS.h"
#include "patch.h"
#include <stdio.h>

void printCombStatsHeader(FILE* file){
	fprintf(file, "iteration\t");
	for (int comb = 0; comb < dataSetup.popTree->numPops; comb++) {
		if (isFeasibleComb(comb)){
			printOneCombHeader(comb, file);
		}
	}

	fprintf(file, "\n");
}
void printOneCombHeader(int comb, FILE* file){
	char* combName = getPopName(comb);

	fprintf(file, "C_%s cs\tC_%s nc\t", combName, combName);

	printCombPopHeaders(comb, combName, file);
	printCombMigHeaders(comb, combName, file);

}
void printCombPopHeaders(int comb, char* combName, FILE* file) {
	for (int pop = 0; pop < dataSetup.popTree->numPops; pop++) {
		if (isLeaf(pop) && isAncestralTo(comb, pop)) {
			char* leafName = getPopName(pop);
			fprintf(file, "C_%s_%s cs\tC_%s_%s nc\t",
					combName, leafName, combName, leafName);
		}
	}
}
void printCombMigHeaders(int comb, char* combName, FILE* file) {
	for (int mig = 0; mig < dataSetup.popTree->numMigBands; mig++) {
		if (isMigOfComb(mig, comb)) {
			int source = getSourcePop(mig);
			int target = getTargetPop(mig);
			char* sourceName = getPopName(source);
			char* targetName = getPopName(target);
			fprintf(file, "C_%s_%s->%s ms\tC_%s_%s->%s nm\t",
					combName, sourceName, targetName, combName, sourceName, targetName);
		}
	}
}

void printCombStats(int iteration, FILE* file){

	fprintf(file, "%d\t", iteration);
	for (int comb = 0; comb < dataSetup.popTree->numPops; comb++) {
		if (isFeasibleComb(comb)){
			printOneCombStats(comb, file);
		}
	}

	fprintf(file, "\n");
	fflush(file);
}
void printOneCombStats(int comb, FILE* file){
	fprintf(file, "%0.35f\t%d\t",  //TODO - turn float percision to a MACRO/CONSTS
			comb_stats[comb].total.coal_stats,
			comb_stats[comb].total.num_coals);
	printCombCoalStats(comb, file);
	printCombMigStats(comb, file);
}
void printCombCoalStats(int comb, FILE* file) {
	for (int pop = 0; pop < dataSetup.popTree->numPops; pop++) {
		if (isLeaf(pop) && isAncestralTo(comb, pop)) {
			fprintf(file, "%0.35f\t%d\t",
					comb_stats[comb].leaves[pop].below_comb.coal_stats,
					comb_stats[comb].leaves[pop].below_comb.num_coals);
		}
	}
}
void printCombMigStats(int comb, FILE* file) {
	for (int mig = 0; mig < dataSetup.popTree->numMigBands; mig++) {
		if (isMigOfComb(mig, comb)) {
			fprintf(file, "%0.35f\t%d\t",
					comb_stats[comb].leafMigs[mig].mig_stats,
					comb_stats[comb].leafMigs[mig].num_migs);
		}
	}
}
