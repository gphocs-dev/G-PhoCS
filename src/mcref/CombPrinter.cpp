#include "CombPrinter.h"
#include "CombStats.h"
#include "../MCMCcontrol.h"
#include "../patch.h"
#include "McRefCommon.h"

void printCombStatsHeader(FILE *file) {
  fprintf(file, "iteration");
  for (int comb = 0; comb < dataSetup.popTree->numPops; comb++) {
    if (isFeasibleComb(comb)) {
      printOneCombHeader(comb, file);
    }
  }

  fprintf(file, "\n");
}

void printOneCombHeader(int comb, FILE *file) {
  char *combName = getPopName(comb);

  fprintf(file, "\tC_%s cs\tC_%s nc", combName, combName);

  printCombPopHeaders(comb, combName, file);
//  printCombMigHeaders(comb, combName, file);

}

void printCombPopHeaders(int comb, char *combName, FILE *file) {
  for (int pop = 0; pop < dataSetup.popTree->numPops; pop++) {
    if (isLeafPop(pop) && isAncestralTo(comb, pop)) {
      char *leafName = getPopName(pop);
      fprintf(file, "\tC_%s_%s cs\tC_%s_%s nc",
              combName, leafName, combName, leafName);
    }
  }
}


void printCombStats(int iteration, FILE *file) {

  fprintf(file, "%d", iteration);
  for (int comb = 0; comb < dataSetup.popTree->numPops; comb++) {
    if (isFeasibleComb(comb)) {
      printOneCombStats(comb, file);
    }
  }
  fprintf(file, "\n");
  fflush(file);
}

void printOneCombStats(int comb, FILE *file) {
  fprintf(file, "\t%0.35f\t%d",  //TODO - turn float precision to a MACRO/CONST
          comb_stats[comb].total.coal_stats,
          comb_stats[comb].total.num_coals);
  printCombCoalStats(comb, file);
  printCombMigStats(comb, file);
}

void printCombCoalStats(int comb, FILE *file) {
  for (int pop = 0; pop < dataSetup.popTree->numPops; pop++) {
    if (isLeafPop(pop) && isAncestralTo(comb, pop)) {
      fprintf(file, "\t%0.35f\t%d",
              comb_stats[comb].leaves[pop].below_comb.coal_stats,
              comb_stats[comb].leaves[pop].below_comb.num_coals);
    }
  }
}

void printCombMigStats(int comb, FILE *file) {
  for (int mig = 0; mig < dataSetup.popTree->numMigBands; mig++) {
    if (isCombLeafMigBand(mig, comb)) {
      fprintf(file, "\t%0.35f\t%d",
              comb_stats[comb].leafMigs[mig].mig_stats,
              comb_stats[comb].leafMigs[mig].num_migs);
    }
  }
}

