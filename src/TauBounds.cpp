#include <map>

#include "TauBounds.h"
#include "patch.h"
#include "GPhoCS.h"
#include "McRefCommon.h"

double *tau_bounds;
int *lca_pops;


void calculateTauBounds() {
  initializeTauBounds();
  for (int gen = 0; gen < dataSetup.numLoci; gen++) {
    initializeLcaPops(gen);
    int rootNodeId = dataState.lociData[gen]->root;
    calculateLociTauBounds(rootNodeId, gen);
  }
}

void calculateLociTauBounds(int nodeId, int gen) {
  LikelihoodNode *currentNode = getNode(nodeId, gen);

  if (isLeafNode(currentNode)) {
    return;
  }

  calculateLociTauBounds(currentNode->leftSon, gen);
  calculateLociTauBounds(currentNode->rightSon, gen);

  lca_pops[nodeId] = lca(lca_pops[currentNode->leftSon], lca_pops[currentNode->rightSon]);

  updateTauBoundsOfDescendants(lca_pops[nodeId], currentNode->age);
}

void updateTauBoundsOfDescendants(int pop, double bound) {
  if (isLeaf(pop)) return;

  if (bound < dataSetup.popTree->pops[pop]->age) {
    //nocommit
    printf("error in tau_bounds");
    exit(-1);
  }

  tau_bounds[pop] = fmin(tau_bounds[pop], bound);

  updateTauBoundsOfDescendants(getSon(pop, LEFT), bound);
  updateTauBoundsOfDescendants(getSon(pop, RIGHT), bound);
}

void initializeLcaPops(int gen) {
  for (int nodeId = 0; nodeId < numNodes(); nodeId++) {
    LikelihoodNode *currentNode = getNode(nodeId, gen);
    if (isLeafNode(currentNode)) {
      lca_pops[nodeId] = getLeafNodePop(nodeId, gen);
    }
  }
}


void printTauBoundsHeader(FILE *file) {
  fprintf(file, "iteration");
  for (int pop = 0; pop < dataSetup.popTree->numPops; pop++) {
    if (isLeaf(pop)) continue;
    fprintf(file, "\t%s_bound\t%s_tau", getPopName(pop), getPopName(pop));
  }
  fprintf(file, "\n");
}

void printTauBounds(int iteration, FILE *file) {
  fprintf(file, "%d", iteration);
  for (int pop = 0; pop < dataSetup.popTree->numPops; pop++) {
    if (isLeaf(pop)) continue;
    fprintf(file, "\t%.40f\t%.40f", tau_bounds[pop], dataSetup.popTree->pops[pop]->age);
  }
  fprintf(file, "\n");
}

void allocateTauBoundsMem() {
  tau_bounds = (double *) malloc(dataSetup.popTree->numPops * sizeof(double));


  lca_pops = (int *) malloc(numNodes() * sizeof(int));
}

void initializeTauBounds() {
  for (int pop = 0; pop < dataSetup.popTree->numPops; pop++) {
    if (!isLeaf(pop)) {
      tau_bounds[pop] = 100.0; // TODO - init with proper value (tau of father pop? MAX/INF value?)
    }
  }
}
