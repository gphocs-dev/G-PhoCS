#include "CladePrinter.h"
#include "MemoryMng.h"
#include "McRefCommon.h"
#include "CladeStats.h"
#include <stdio.h>

void printCladeStatsHeader(FILE* file){
	fprintf(file, "iteration\t");
	for (int clade = 0; clade < dataSetup.popTree->numPops; clade++) {
		printSpecificCladeHeader(clade, file);
	}
	for (int pop = 0; pop < dataSetup.popTree->numPops; pop++) {
		printSpecificPopHeader(pop, file);
	}
	for(int mig_band=0; mig_band<dataSetup.popTree->numMigBands; mig_band++) {
			printSpecificMigHeader(mig_band, file);
	}

	fprintf(file, "\n");
}
void printSpecificCladeHeader(int clade, FILE* file){
	char* cladeName = dataSetup.popTree->popArray[clade].name;
	fprintf(file, "%s_%s\t%s_%s\t",
			cladeName, "coal_stats_total",
			cladeName, "num_coals_total");
}
void printSpecificPopHeader(int pop, FILE* file){
	char* popName = dataSetup.popTree->popArray[pop].name;
	fprintf(file, "%s_%s\t%s_%s\t",
			popName, "_pop_coal_stats_total",
			popName, "_pop_num_coals_total");
}
void printSpecificMigHeader(int mig_band, FILE* file){
	char* sourceName = getPopName(dataSetup.popTree->migBands[mig_band].sourcePop);
	char* targetName = getPopName(dataSetup.popTree->migBands[mig_band].targetPop);
	fprintf(file, "%s->%s_%s\t%s->%s_%s\t",
			sourceName, targetName, "mig_stats",
			sourceName, targetName, "num_migs");
}

void printCladeStats(int iteration, FILE* file){
	fprintf(file, "%d\t", iteration);
	for (int clade = 0; clade < dataSetup.popTree->numPops; clade++) {
		printSpecificCladeStats(clade, file);
	}
	for (int pop = 0; pop < dataSetup.popTree->numPops; pop++) {
		printSpecificPopStats(pop, file);
	}
	for(int mig_band=0; mig_band<dataSetup.popTree->numMigBands; mig_band++) {
		printSpecificMigBandStats(mig_band, file);
	}
	fprintf(file, "\n");
	fflush(file);
}

void printSpecificCladeStats(int clade, FILE* file){
	fprintf(file, "%0.35f\t%d\t",
			clade_stats[clade].coal_stats_total,
			clade_stats[clade].num_coals_total);
}
void printSpecificPopStats(int pop, FILE* file){
	fprintf(file, "%0.35f\t%d\t",
			genetree_stats_total.coal_stats[pop],
			genetree_stats_total.num_coals[pop]);
}
void printSpecificMigBandStats(int mig_band, FILE* file){
	fprintf(file, "%0.35f\t%d\t",
			genetree_stats_total.mig_stats[mig_band],
			genetree_stats_total.num_migs[mig_band]);
}
