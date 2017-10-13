#include <unordered_set>
#include "CombStats.h"
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

int isAncestralTo(int father, int son) {
  return dataSetup.popTree->pops[father]->isAncestralTo[son];
}

const char *getEventTypeName(int eventType) {
  switch (eventType) {
    case COAL:
      return "COAL";
    case IN_MIG:
      return "IN_MIG";
    case OUT_MIG:
      return "OUT_MIG";
    case MIG_BAND_START:
      return "MIG_START";
    case MIG_BAND_END:
      return "MIG_END";
    case SAMPLES_START:
      return "SAM_START";
    case END_CHAIN:
      return "END_CHAIN";
    case DUMMY:
      return "DUMMY";
    default:
      return "UNDEFINED";
  }
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
//char* getMigName(int mig){
//	char* sourceName = getPopName(getSourcePop(mig));
//	char* targetName = getPopName(getTargetPop(mig));
//	return concat(sourceName, targetName);
//}
//char* concat(const char *s1, const char *s2){ // THIS SHOULD NOT BE USED IN PRODUCTION! This memory isn't released
//    char *result = malloc(strlen(s1)+strlen(s2)+3);//+3 for the zero-terminator and arrow
//    strcpy(result, s1);
//    strcat(result, "->");
//    strcat(result, s2);
//    return result;
//}



int lca(int pop1, int pop2) {
  std::unordered_set<int> pop1_ancestors = {}; //TODO - do I need to release this set?
  int pop1_ancestor = pop1;
  while (pop1_ancestor != -1) {
    pop1_ancestors.insert(pop1_ancestor);
    pop1_ancestor = getPopFather(pop1_ancestor);
  }


  int pop2_ancestor = pop2;
  while (pop2_ancestor != -1) {
    if (pop1_ancestors.count(pop2_ancestor)) {
      return pop2_ancestor;
    }
    pop2_ancestor = getPopFather(pop2_ancestor);
  }

  return -1;
}

int getPopFather(int popId) {
  Population *pop = dataSetup.popTree->pops[popId];
  if (pop && pop->father) return pop->father->id;
  else return -1;
}



LikelihoodNode *getNode(int nodeId, int gen) { // TODO - move to util
  LocusData *loci = dataState.lociData[gen];
  return loci->nodeArray[nodeId];
}

int getNodePop(int nodeId, int gen) {
  return nodePops[gen][nodeId];
}

bool isLeafNode(LikelihoodNode *node) {
  return (node->leftSon == -1) || (node->rightSon == -1);  //TOASK - make sure this is correct
}

int numNodes() { return (2 * dataSetup.numSamples) - 1; } // TODO - make MACRO?



bool hasNextEvent(EventChain chain, int event) {
  int next = chain.events[event].getNextIdx();
  return next >= 0;
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