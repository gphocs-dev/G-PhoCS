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
  fprintf(file, "iteration");
  for (int band = 0; band < dataSetup.popTree->numMigBands; band++) {
    int sourceId = dataSetup.popTree->migBands[band].sourcePop;
    int targetId = dataSetup.popTree->migBands[band].targetPop;
    Population &source = dataSetup.popTree->popArray[sourceId];
    Population &target = dataSetup.popTree->popArray[targetId];
    fprintf(file, "\t%s->%s", source.name, target.name); // TODO - refactor migbandName to some member of `MigrationBand`
  }
  fprintf(file, "\n");
}

void calculateReferenceMigrationStats() {
  initRefMigStats();
  printf("calculating migrefstats\n");
}

void printReferenceMigrationStats(int iteration, FILE *file) {
  printf("printing migrefstats\n");
}
