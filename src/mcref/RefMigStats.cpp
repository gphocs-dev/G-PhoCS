#include <cstdio>
#include "RefMigStats.h"
#include "../patch.h"
#include "../TauBounds.h"
#include "../GPhoCS.h"
#include "../McRefCommon.h"
#include "../MemoryMng.h"

double *migBandRefStats;


void calculateReferenceMigrationStats() {
  initRefMigStats();
  for (int gen = 0; gen < dataSetup.numLoci; gen++) {
    int rootNodeId = dataState.lociData[gen]->root;
    double rootNodeAge = getNode(rootNodeId, gen)->age; // TODO - make sure this workaround makes sense
    calculateReferenceGenMigStats(rootNodeId, rootNodeAge, gen);
  }
}

int calculateReferenceGenMigStats(int nodeId, double parentNodeAge, int gen) {
  int nodeLcaPop, leftLcaPop, rightLcaPop, migNodeId;
  double currentNodeAge = 0.0, bottomAge, topAge;
  LikelihoodNode *currentNode;

  /*** determine population and age of current node ***/
  if (isLeafNode(nodeId, gen)) {
    nodeLcaPop = getNodePop(nodeId, gen);
  } else {
    currentNode = getNode(nodeId, gen);
    currentNodeAge = currentNode->age;
    leftLcaPop = calculateReferenceGenMigStats(currentNode->leftSon, currentNodeAge, gen);
    rightLcaPop = calculateReferenceGenMigStats(currentNode->rightSon, currentNodeAge, gen);
    nodeLcaPop = lca_pops[leftLcaPop][rightLcaPop];
  }

  /*** iterate through all migration events above current node ***/
  bottomAge = currentNodeAge;
  while (TRUE) {
    migNodeId = getMigNodeAbove(nodeId, gen, bottomAge);
    topAge = migNodeId >= 0 ? genetree_migs[gen].mignodes[migNodeId].age : parentNodeAge;

    /*** add contribution of edge fragment to all relevant mig bands - code copied from addChildEdgesToMigStats() ***/
    for (int mb = 0; mb < dataSetup.popTree->numMigBands; mb++) {
      int target = dataSetup.popTree->migBands[mb].targetPop;
      if (isAncestralTo(target, nodeLcaPop)) {
        migBandRefStats[mb] += migBandIntersection(mb, bottomAge, topAge);
      }
    }

    /*** if no more mig nodes, then nothing more to do ***/
    if (migNodeId < 0) break;

    /*** Otherwise, update pop and bottomAge ***/
    nodeLcaPop = genetree_migs[gen].mignodes[migNodeId].source_pop;
    bottomAge = topAge;
  }

  return nodeLcaPop;
}


double migBandIntersection(int mb, double bottomAge, double topAge) {
  int target = dataSetup.popTree->migBands[mb].targetPop;
  int source = dataSetup.popTree->migBands[mb].sourcePop;
  double mbUbound = fmin(tau_ubounds[getFather(source)], tau_ubounds[getFather(target)]);

  return fmax(0.0, fmin(mbUbound, topAge) - bottomAge);
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
    fprintf(file, "\t%s->%s", source.name, target.name); // TODO - refactor migbandName to a member of `MigrationBand`
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
