#include "CombStats.h"
#include "MCMCcontrol.h"
#include "McRefCommon.h"
#include "patch.h"


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

int isLeaf(int pop) {
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
  if (isLeaf(pop)) {
    return FALSE;
  }
  int left_son = getSon(pop, LEFT);
  int right_son = getSon(pop, RIGHT);
  return (isLeaf(left_son) && isLeaf(right_son));
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