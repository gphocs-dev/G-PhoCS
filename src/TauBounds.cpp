
#include <map>

#include "TauBounds.h"
#include "patch.h"
#include "GPhoCS.h"
#include "McRefCommon.h"

double *tau_bounds; // Temporary array. Reinitialized per iteration. Holds final tau bounds
int **lca_pops; // "cache" of the lca of every pair of pops. i.e. - lca_pops[1][2] is the lca of pops 1 & 2


void calculateTauBounds() {
  initializeBounds();

  for (int gen = 0; gen < dataSetup.numLoci; gen++) {
    int rootNodeId = dataState.lociData[gen]->root;
    calculateLociTauBounds(rootNodeId, gen);
  }

  propagateBoundsDownPopTree();

  if (DEBUG_TAU_BOUNDS) runTauBoundsAssertions();
}

int calculateLociTauBounds(int nodeId, int gen) {
  LikelihoodNode *currentNode = getNode(nodeId, gen);
  int pop = getNodePop(nodeId, gen);

  if (isLeafPop(pop))
    return pop;

  int leftLca = calculateLociTauBounds(currentNode->leftSon, gen);
  int rightLca = calculateLociTauBounds(currentNode->rightSon, gen);

  int lca_pop = lca_pops[leftLca][rightLca];

  tau_bounds[lca_pop] = fmin(tau_bounds[lca_pop], currentNode->age);

  return lca_pop;
}

void propagateBoundsDownPopTree() {
  for (int pop = 0; pop < dataSetup.popTree->numPops; pop++) {
    updateTauBoundsOfDescendants(pop, tau_bounds[pop]);
  }
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
  lca_pops = (int **) malloc(dataSetup.popTree->numPops * sizeof(int *));
  for (int pop = 0; pop < dataSetup.popTree->numPops; pop++) {
    lca_pops[pop] = (int *) malloc(dataSetup.popTree->numPops * sizeof(int));
  }
  computeLcas();
}

void computeLcas() {
  for (int pop1 = 0; pop1 < dataSetup.popTree->numPops; pop1++) {
    for (int pop2 = 0; pop2 < dataSetup.popTree->numPops; pop2++) {
      lca_pops[pop1][pop2] = lca(pop1, pop2);
    }
  }
}

void initializeBounds() {
  for (int pop = 0; pop < dataSetup.popTree->numPops; pop++) {
    if (isLeafPop(pop)) {
      tau_bounds[pop] = 0.0;
    } else {
      tau_bounds[pop] = 100.0; // TODO - init with proper value (tau of father pop? MAX/INF value?)
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
