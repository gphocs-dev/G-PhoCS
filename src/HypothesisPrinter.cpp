#include "CombPrinter.h"
#include "HypothesisPrinter.h"
#include "patch.h"
#include "McRefCommon.h"
#include "MemoryMng.h"

void printHypStatsHeader(FILE *file) {

  fprintf(file, "iteration\t");
  for (int pop = 0; pop < dataSetup.popTree_->numPops; pop++) {
    printOnePopHeader(pop, file);
  }

  for (int migband = 0; migband < dataSetup.popTree_->numMigBands; migband++) {
    printOneMigBandHeader(migband, file);
  }

  fprintf(file, "\n");

}

void printOnePopHeader(int pop, FILE *file) {
  char *popName = dataSetup.popTree_->popArray[pop].name;
  fprintf(file, "P_%s %s\tP_%s %s\t",
          popName, "cs",
          popName, "nc");
}

void printOneMigBandHeader(int migband, FILE *file) {
  char *source = getPopName(dataSetup.popTree_->migBands[migband].sourcePop);
  char *target = getPopName(dataSetup.popTree_->migBands[migband].targetPop);
  fprintf(file, "MB_%s->%s ms\tMB_%s->%s nm\t",
          source, target,
          source, target);
}

void printHypStats(int iteration, FILE *file) {

  fprintf(file, "%d\t", iteration);

  for (int pop = 0; pop < dataSetup.popTree_->numPops; pop++) {
    printOnePopStats(pop, file);
  }

  for (int migband = 0; migband < dataSetup.popTree_->numMigBands; migband++) {
    printOneMigBandStats(migband, file);
  }

  fprintf(file, "\n");
  fflush(file);
}

void printOnePopStats(int pop, FILE *file) {
  fprintf(file, "%0.35f\t%d\t",
          genetree_stats_total.coal_stats[pop],
          genetree_stats_total.num_coals[pop]);
}

void printOneMigBandStats(int migband, FILE *file) {
  fprintf(file, "%0.35f\t%d\t",
          genetree_stats_total.mig_stats[migband],
          genetree_stats_total.num_migs[migband]);
}