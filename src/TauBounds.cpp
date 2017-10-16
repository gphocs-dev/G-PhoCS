
#include <map>

#include "TauBounds.h"
#include "patch.h"
#include "GPhoCS.h"
#include "McRefCommon.h"

double *tau_bounds; // Temporary array. Reinitialized per iteration. Holds final tau bounds
int *lca_pops; // Temporary array. Reinitialized per genealogy. Holds lca_pop of nodes of current genealogy


void calculateTauBounds() {
  initializeTauBounds();
  for (int gen = 0; gen < dataSetup.numLoci; gen++) {
    initializeLcaPops(gen);
    int rootNodeId = dataState.lociData[gen]->root;
    calculateLociTauBounds(rootNodeId, gen);
  }
  if (DEBUG_TAU_BOUNDS) runTauBoundsAssertions();
}

void calculateLociTauBounds(int nodeId, int gen) {
  LikelihoodNode *currentNode = getNode(nodeId, gen);

  if (isLeafNode(currentNode))
    return;

  calculateLociTauBounds(currentNode->leftSon, gen);
  calculateLociTauBounds(currentNode->rightSon, gen);

  lca_pops[nodeId] = lca(lca_pops[currentNode->leftSon], lca_pops[currentNode->rightSon]);

  updateTauBoundsOfDescendants(lca_pops[nodeId], currentNode->age);
}

void updateTauBoundsOfDescendants(int pop, double bound) {
  if (isLeafPop(pop)) return;

  tau_bounds[pop] = fmin(tau_bounds[pop], bound);

  updateTauBoundsOfDescendants(getSon(pop, LEFT), bound);
  updateTauBoundsOfDescendants(getSon(pop, RIGHT), bound);
}


void printTauBoundsHeader(FILE *file) {
  fprintf(file, "iteration");
  for (int pop = 0; pop < dataSetup.popTree->numPops; pop++) {
    if (isLeafPop(pop)) continue;
    fprintf(file, "\t%s_bound\t%s_tau", getPopName(pop), getPopName(pop));
  }
  fprintf(file, "\n");
}

void printTauBounds(int iteration, FILE *file) {
  fprintf(file, "%d", iteration);
  for (int pop = 0; pop < dataSetup.popTree->numPops; pop++) {
    if (isLeafPop(pop)) continue;
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
    if (isLeafPop(pop)) {
      tau_bounds[pop] = 0.0;
    } else {
      tau_bounds[pop] = 100.0; // TODO - init with proper value (tau of father pop? MAX/INF value?)
    }
  }
}

void initializeLcaPops(int gen) {
  for (int nodeId = 0; nodeId < numNodes(); nodeId++) {
    LikelihoodNode *currentNode = getNode(nodeId, gen);
    if (isLeafNode(currentNode)) {
      lca_pops[nodeId] = getNodePop(nodeId, gen);
    }
  }
}


void runTauBoundsAssertions() {
  assertBoundsAreMonotonousAscending();
  assertTausAreSmallerThanBounds();
}

void assertTausAreSmallerThanBounds() {
  for (int pop = 0; pop < dataSetup.popTree->numPops; pop++) {
    if (getPopAge(pop) > tau_bounds[pop]) {
      printf("age of %s is over its tau bound.\n", getPopName(pop));
      printf("%.10f is the age and-\n %.10f is the bound", getPopAge(pop), tau_bounds[pop]);
      exit(-1);
    }
  }
}

void assertBoundsAreMonotonousAscending() {
  for (int pop = 0; pop < dataSetup.popTree->numPops; pop++) {
    for (int ancestor = 0; ancestor < dataSetup.popTree->numPops; ancestor++) {
      if (isAncestralTo(ancestor, pop)) {
        if (tau_bounds[pop] > tau_bounds[ancestor]) {
          printf("\nTau bound of pop %s was larger than of his ancestor %s", getPopName(pop), getPopName(ancestor));
          printf("\n%.10f is the bound of %s and - ", tau_bounds[pop], getPopName(pop));
          printf("\n%.10f is the bound of %s.", tau_bounds[ancestor], getPopName(ancestor));
          exit(-1);
        }
      }
    }
  }
}
