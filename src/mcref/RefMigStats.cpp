#include <cstdio>
#include "RefMigStats.h"
#include "../patch.h"

double *migBandRefStats;


void allocateRefMigStatsMem() {
  migBandRefStats = (double *) malloc(dataSetup.popTree->numMigBands * sizeof(double));
}

void initRefMigStats() {
  for (int band = 0; band < dataSetup.popTree->numMigBands; band++) {
    migBandRefStats[band] = 0.0;
  }
}

void printRefMigStatsHeader(FILE *file) {
  printf("printing migrefstats header\n");
}

void calculateReferenceMigrationStats() {
  initRefMigStats()
  printf("calculating migrefstats\n");
}

void printReferenceMigrationStats(int iteration, FILE *file) {
  printf("printing migrefstats\n");
}
