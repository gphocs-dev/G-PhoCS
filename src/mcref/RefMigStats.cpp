#include <cstdio>
#include "RefMigStats.h"
#include "../patch.h"
#include "../TauBounds.h"
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
  int lcaPop = lca_pops[leftLca][rightLca];

  addChildEdgesToMigStats(currentNode, gen, lcaPop);

  return migLcaPop(nodeId, gen, lcaPop);
}

void addChildEdgesToMigStats(const LikelihoodNode *currentNode, int gen, int lcaPop) {
  for (int mb = 0; mb < dataSetup.popTree->numMigBands; mb++) {
    int target = dataSetup.popTree->migBands[mb].targetPop;
    if (isAncestralOrEqual(target, lcaPop)) {
      LikelihoodNode *leftSon = getNode(currentNode->leftSon, gen);
      migBandRefStats[mb] += migBandIntersection(mb, leftSon->age, currentNode->age);
      LikelihoodNode *rightSon = getNode(currentNode->rightSon, gen);
      migBandRefStats[mb] += migBandIntersection(mb, rightSon->age, currentNode->age);
    }
  }
}

double migBandIntersection(int mb, double fromAge, double toAge) {
  int target = dataSetup.popTree->migBands[mb].targetPop;
  int source = dataSetup.popTree->migBands[mb].sourcePop;
  double mbUbound = fmin(tau_ubounds[getFather(source)], tau_ubounds[getFather(target)]);

  return fromAge < mbUbound ? fmin(mbUbound, toAge) - fromAge : 0.0;
}

void initRefMigStats() {
  for (int mb = 0; mb < dataSetup.popTree->numMigBands; mb++) {
    migBandRefStats[mb] = 0.0;
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
