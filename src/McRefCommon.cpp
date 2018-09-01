#include "MCMCcontrol.h"
#include "McRefCommon.h"
#include "patch.h"
#include "GPhoCS.h"
#include "MemoryMng.h"


double calculateCoalStats(double *elapsed_times, int *num_lineages, int size) {
  int n;
  double t;
  double result = 0.0;
  for (int i = 0; i < size; i++) {
    n = num_lineages[i];
    t = elapsed_times[i];
    result += n * (n - 1) * t;
  }
  return result;
}

int isLeafNode(int nodeId, int gen) {
  LocusData *locus = dataState.lociData[gen];
  return nodeId < locus->numLeaves;
}

int isLeafPop(int pop) {
  Population *population, *left_son, *right_son;

  population = dataSetup.popTree->pops[pop];

  left_son = population->sons[LEFT];
  right_son = population->sons[RIGHT];

  if (left_son || right_son) {
    return FALSE;
  } else {
    return TRUE;
  }
}

int areChildrenLeaves(int pop) {
  if (isLeafPop(pop)) {
    return FALSE;
  }
  int left_son = getSon(pop, LEFT);
  int right_son = getSon(pop, RIGHT);
  return (isLeafPop(left_son) && isLeafPop(right_son));
}

int isAncestralTo(int ancestor, int descendant) {
  return dataSetup.popTree->pops[ancestor]->isAncestralTo[descendant];
}
int isAncestralOrEqual(int ancestor, int descendant) {
  return ancestor == descendant || isAncestralTo(ancestor, descendant);
}

int getSon(int pop, int SON) {
  return dataSetup.popTree->pops[pop]->sons[SON]->id;
}


double getPopAge(int pop) {
  return dataSetup.popTree->pops[pop]->age;
}

int getSourcePop(int mig) {
  return dataSetup.popTree->migBands[mig].sourcePop;
}

int getTargetPop(int mig) {
  return dataSetup.popTree->migBands[mig].targetPop;
}

char *getPopName(int popId) {
  return dataSetup.popTree->pops[popId]->name;
}

int lca(int pop1, int pop2) {
  int ancestor1 = pop1;
  while (ancestor1 != -1) {
    int ancestor2 = pop2;
    while (ancestor2 != -1) {
      if (ancestor1 == ancestor2)
        return ancestor1;
      ancestor2 = getFather(ancestor2);
    }
    ancestor1 = getFather(ancestor1);
  }
  return -1;
}

int getFather(int popId) {
  Population *pop = dataSetup.popTree->pops[popId];
  if (pop && pop->father) return pop->father->id;
  else return -1;
}


LikelihoodNode *getNode(int nodeId, int gen) {
  LocusData *loci = dataState.lociData[gen];
  return loci->nodeArray[nodeId];
}

int getNodePop(int nodeId, int gen) {
  return nodePops[gen][nodeId];
}

bool hasNextEvent(EventChain chain, int event) {
  int next = chain.events[event].getNextIdx();
  return next >= 0;
}

/**
 *  if nodeId is right below a migration event, returns the source population of the oldest such migration event.
 *  Otherwise, returns defaultLcaPop
 */
int migLcaPop(int nodeId, int gen, int defaultLcaPop) { // TODO - refactor to utilize getMigNodeAbove
  int oldestMigSource = -1;
  double oldestMigAge = -1.0;
  GENETREE_MIGS &genMigs = genetree_migs[gen];

  for (int i = 0; i < genMigs.num_migs; i++) {
    int mig = genMigs.living_mignodes[i];
    if ((genMigs.mignodes[mig].gtree_branch == nodeId) && (genMigs.mignodes[mig].age > oldestMigAge)) {
      oldestMigSource = genMigs.mignodes[mig].source_pop;
      oldestMigAge = genMigs.mignodes[mig].age;
    }
  }
  return oldestMigSource == -1 ? defaultLcaPop : oldestMigSource;
}

int getMigNodeAbove(int nodeId, int gen, double requiredAge) {
  int i, mig;
  GENETREE_MIGS &genMigs = genetree_migs[gen];

  for (i = 0, mig = genMigs.living_mignodes[i]; i < genMigs.num_migs; i++) {
    if ((genMigs.mignodes[mig].gtree_branch == nodeId) && (genMigs.mignodes[mig].age > requiredAge)) {
      return mig;
    }
  }
  return -1;
}


bool areAlmostEqual(double eventAge, double combAge) {
  return relativeDistance(eventAge, combAge) <= requiredRelativePrecision();
}

double relativeDistance(double dbl1, double dbl2) {
  double absDist = fabs(dbl1 - dbl2);
  double sum = dbl1 + dbl2;
  if (sum == 0.0) return 0.0;
  else return (absDist * 2) / sum;
}

double RELATIVE_PRECISION = 0.000000000001;

double requiredRelativePrecision() {
  return RELATIVE_PRECISION;
}

int **lca_pops;
void computeLcas() {
  for (int pop1 = 0; pop1 < dataSetup.popTree->numPops; pop1++) {
    for (int pop2 = 0; pop2 < dataSetup.popTree->numPops; pop2++) {
      lca_pops[pop1][pop2] = lca(pop1, pop2);
    }
  }
}