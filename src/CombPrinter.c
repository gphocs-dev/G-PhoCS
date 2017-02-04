/***********************************************************************************
 *	printCoalStats
 * 	- prints coalescent statistics to file
 *  - this includes stats needed to compute null likelihood as well as stats on coal times
 *	- returns 0
 ***********************************************************************************/

#include "CombPrinter.h"
#include "CombStats.h"

#include "MCMCcontrol.h"

#include "GPhoCS.h"
#include "patch.h"


#include <stdio.h>

void printCombStats(int iteration, FILE* file){

	return;

	if(iteration<0) {
		printCombStatsHeader(file);

		return;
	}

	fprintf(file, "%d\t", iteration);
	for (int comb = 0; comb < dataSetup.popTree->numPops; comb++) {
		printSpecificCombStats(comb, file);
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

void printCombStatsHeader(FILE* file){
	fprintf(file, "iteration\t");
	for (int comb = 0; comb < dataSetup.popTree->numPops; comb++) {
		printSpecificCombHeader(comb, file);
	}
	for (int pop = 0; pop < dataSetup.popTree->numPops; pop++) {
		printSpecificPopHeader(pop, file);
	}
	for(int mig_band=0; mig_band<dataSetup.popTree->numMigBands; mig_band++) {
		printSpecificMigHeader(mig_band, file);
	}

	fprintf(file, "\n");
}
void printSpecificCombHeader(int comb, FILE* file){
	char* combName = dataSetup.popTree->popArray[comb].name;
	fprintf(file, "%s_%s\t%s_%s\t",
			combName, "coal_stats_total",
			combName, "num_coals_total");
}
void printSpecificPopHeader(int pop, FILE* file){
	char* popName = dataSetup.popTree->popArray[pop].name;
	fprintf(file, "%s_%s\t%s_%s\t",
			popName, "_pop_coal_stats_total",
			popName, "_pop_num_coals_total");
}
void printSpecificMigHeader(int mig_band, FILE* file){
//	char* migBandName = dataSetup.popTree->migBands[mig_band].name;
	char* migBandName = "unimplemented";
	fprintf(file, "%s_%s\t%s_%s\t",
			migBandName, "_mig_stats",
			migBandName, "_num_migs");
}


void printSpecificCombStats(int comb, FILE* file){
	fprintf(file, "%0.35f\t%d\t",
			comb_stats[comb].coal_stats_total,
			comb_stats[comb].num_coals_total);
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
