#include <cstdio>
#include "RefMigStats.h"
#include "../patch.h"
#include "../GPhoCS.h"
#include "../McRefCommon.h"

double *migBandRefStats;


void calculateReferenceMigrationStats() {
  initRefMigStats();
  for (int gen = 0; gen < dataSetup.numLoci; gen++) {
    int rootNodeId = dataState.lociData[gen]->root;
    calculateReferenceGenMigStats(rootNodeId, gen);
  }
}

int calculateReferenceGenMigStats(int nodeId, int gen) {
  int currPop = getNodePop(nodeId, gen);
  if (isLeafNode(nodeId, gen))
    return migLcaPop(nodeId, gen, currPop);

  LikelihoodNode *currentNode = getNode(nodeId, gen);
  int leftLca = calculateReferenceGenMigStats(currentNode->leftSon, gen);
  int rightLca = calculateReferenceGenMigStats(currentNode->rightSon, gen);

  // TODO - calculate migstats of both child nodes

  int lca_pop = lca_pops[leftLca][rightLca];
  return migLcaPop(nodeId, gen, lca_pop);
}

void initRefMigStats() {
  for (int band = 0; band < dataSetup.popTree->numMigBands; band++) {
    migBandRefStats[band] = 0.0;
  }
}







// ==================================================

void allocateRefMigStatsMem() {
  migBandRefStats = (double *) malloc(dataSetup.popTree->numMigBands * sizeof(double));
}

void printRefMigStatsHeader(FILE *file) {
  fprintf(file, "iteration");
  for (int band = 0; band < dataSetup.popTree->numMigBands; band++) {
    int sourceId = dataSetup.popTree->migBands[band].sourcePop;
    int targetId = dataSetup.popTree->migBands[band].targetPop;
    Population &source = dataSetup.popTree->popArray[sourceId];
    Population &target = dataSetup.popTree->popArray[targetId];
    fprintf(file, "\t%s->%s", source.name,
            target.name); // TODO - refactor migbandName to some member of `MigrationBand`
  }
  fprintf(file, "\n");
}

void printReferenceMigrationStats(int iteration, FILE *file) {
  fprintf(file, "%d", iteration);
  for (int band = 0; band < dataSetup.popTree->numMigBands; band++) {
    fprintf(file, "\t%.40f", migBandRefStats[band]);
  }
  fprintf(file, "\n");
  fflush(file);
}
