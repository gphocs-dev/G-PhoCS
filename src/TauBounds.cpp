
#include <map>

#include "TauBounds.h"
#include "patch.h"
#include "GPhoCS.h"
#include "McRefCommon.h"
#include "MemoryMng.h"

double *tau_ubounds; // Temporary array. Reinitialized per iteration. Holds upper tau bounds
double *tau_lbounds; // Ditto, but with lower tau bounds

void calculateTauBounds() {
  initializeBounds();
  calculateUpperAndLowerBounds();
  propagateBoundsAcrossPopTree();
  if (DEBUG_TAU_BOUNDS) runTauBoundsAssertions();
}

void initializeBounds() {
  for (int pop = 0; pop < dataSetup.popTree->numPops; pop++) {
    tau_lbounds[pop] = 0.0;
    if (isLeafPop(pop)) {
      tau_ubounds[pop] = 0.0;
    } else {
      tau_ubounds[pop] = OLDAGE;
    }
  }
}

void calculateUpperAndLowerBounds() {
  int mig, target, targetFather, source, sourceFather, rootNodeId, i;
  double age;

  for (int gen = 0; gen < dataSetup.numLoci; gen++) {

    rootNodeId = dataState.lociData[gen]->root;
    calculateLocusTauUpperBounds(rootNodeId, gen);  // ubound1

    for (i = 0; i < genetree_migs[gen].num_migs; i++) {
      mig = genetree_migs[gen].living_mignodes[i];
      target = genetree_migs[gen].mignodes[mig].target_pop;
      targetFather = getFather(target);
      source = genetree_migs[gen].mignodes[mig].source_pop;
      sourceFather = getFather(source);
      age = genetree_migs[gen].mignodes[mig].age;

      tau_ubounds[source] = fmin(tau_ubounds[source], age); // ubound2
      tau_ubounds[target] = fmin(tau_ubounds[target], age); // ubound2

      tau_lbounds[sourceFather] = fmax(tau_lbounds[sourceFather], age); // lbound
      tau_lbounds[targetFather] = fmax(tau_lbounds[targetFather], age); // lbound
    }
  }
}

int calculateLocusTauUpperBounds(int nodeId, int gen) {
  LikelihoodNode *currentNode = getNode(nodeId, gen);
  int pop = getNodePop(nodeId, gen);

  if (isLeafPop(pop))
    return migLcaPop(nodeId, gen, pop);

  int leftLca = calculateLocusTauUpperBounds(currentNode->leftSon, gen);
  int rightLca = calculateLocusTauUpperBounds(currentNode->rightSon, gen);

  int lca_pop = lca_pops[leftLca][rightLca];
  tau_ubounds[lca_pop] = fmin(tau_ubounds[lca_pop], currentNode->age);
  return migLcaPop(nodeId, gen, lca_pop);
}



void propagateBoundsAcrossPopTree() {
  int root = dataSetup.popTree->rootPop;
  updateLowerBoundsOfDescendants(root);
  updateUpperBoundsOfDescendants(root, tau_ubounds[root]);
}

double updateLowerBoundsOfDescendants(int pop) {
  if (isLeafPop(pop)) return 0.0;

  double left_lbound = updateLowerBoundsOfDescendants(getSon(pop, LEFT));
  double right_lbound = updateLowerBoundsOfDescendants(getSon(pop, RIGHT));
  tau_lbounds[pop] = fmax(tau_lbounds[pop], fmax(left_lbound, right_lbound));
  return tau_lbounds[pop];
}

void updateUpperBoundsOfDescendants(int pop, double bound) {
  if (isLeafPop(pop)) return;

  tau_ubounds[pop] = fmin(tau_ubounds[pop], bound);
  updateUpperBoundsOfDescendants(getSon(pop, LEFT), tau_ubounds[pop]);
  updateUpperBoundsOfDescendants(getSon(pop, RIGHT), tau_ubounds[pop]);
}


/*  ===========================================================================  */

void printTauBoundsHeader(FILE *file) {
  fprintf(file, "iteration");
  for (int pop = 0; pop < dataSetup.popTree->numPops; pop++) {
    if (isLeafPop(pop)) continue;
    fprintf(file, "\t%s_ubound\t%s_lbound\t%s_tau", getPopName(pop), getPopName(pop), getPopName(pop));
  }
  fprintf(file, "\n");
}

void printTauBounds(int iteration, FILE *file) {
  fprintf(file, "%d", iteration);
  for (int pop = 0; pop < dataSetup.popTree->numPops; pop++) {
    if (isLeafPop(pop)) continue;
    fprintf(file, "\t%.40f\t%.40f\t%.40f", tau_ubounds[pop], tau_lbounds[pop], dataSetup.popTree->pops[pop]->age);
  }
  fprintf(file, "\n");
  fflush(file);
}

void allocateTauBoundsMem() {
  tau_ubounds = (double *) malloc(dataSetup.popTree->numPops * sizeof(double));
  tau_lbounds = (double *) malloc(dataSetup.popTree->numPops * sizeof(double));
  lca_pops = (int **) malloc(dataSetup.popTree->numPops * sizeof(int *));
  for (int pop = 0; pop < dataSetup.popTree->numPops; pop++) {
    lca_pops[pop] = (int *) malloc(dataSetup.popTree->numPops * sizeof(int));
  }
  computeLcas();
}

void runTauBoundsAssertions() {
  assertBoundsAreMonotonousAscending();
  assertTausAreBetweenBounds();
}

void assertTausAreBetweenBounds() {
  for (int pop = 0; pop < dataSetup.popTree->numPops; pop++) {
    double tau = getPopAge(pop);
    if (tau > tau_ubounds[pop] || tau < tau_lbounds[pop]) {
      printf("age of %s isnt inside its tau bound.\n", getPopName(pop));
      printf("%.10f is the age and-\n [%.10f, %.10f] is the bound", getPopAge(pop), tau_lbounds[pop], tau_ubounds[pop]);
      exit(-1);
    }
  }
}

void assertBoundsAreMonotonousAscending() {
  for (int pop = 0; pop < dataSetup.popTree->numPops; pop++) {
    for (int ancestor = 0; ancestor < dataSetup.popTree->numPops; ancestor++) {
      if (isAncestralTo(ancestor, pop)) {
        if (tau_ubounds[pop] > tau_ubounds[ancestor] or tau_lbounds[pop] > tau_lbounds[ancestor]) {
          printf("\nTau bound of pop %s was larger than that of his ancestor %s", getPopName(pop),
                 getPopName(ancestor));
          printf("\n[%.10f, %.10f] is the bound of %s and - ", tau_lbounds[pop], tau_ubounds[pop], getPopName(pop));
          printf("\n[%.10f, %.10f] is the bound of %s.", tau_lbounds[ancestor], tau_ubounds[ancestor],
                 getPopName(ancestor));
          exit(-1);
        }
      }
    }
  }
}
