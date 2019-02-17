#include "CladePrinter.h"
#include "MemoryMng.h"
#include "CladeStats.h"

void printCladeStatsHeader(FILE *file) {
  fprintf(file, "iteration\t");
  for (int clade = 0; clade < dataSetup.popTree_->numPops; clade++) {
    printSpecificCladeHeader(clade, file);
  }
  fprintf(file, "\n");
}

void printSpecificCladeHeader(int clade, FILE *file) {
  char *cladeName = dataSetup.popTree_->popArray[clade].name;
  fprintf(file, "%s_%s\t%s_%s\t",
          cladeName, "coal_stats_total",
          cladeName, "num_coals_total");
}

void printCladeStats(int iteration, FILE *file) {
  fprintf(file, "%d\t", iteration);
  for (int clade = 0; clade < dataSetup.popTree_->numPops; clade++) {
    printSpecificCladeStats(clade, file);
  }
  fprintf(file, "\n");
  fflush(file);
}

void printSpecificCladeStats(int clade, FILE *file) {
  fprintf(file, "%0.35f\t%d\t",
          clade_stats[clade].coal_stats_total,
          clade_stats[clade].num_coals_total);
}
