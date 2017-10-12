#include <map>

#include "TauBounds.h"
#include "patch.h"
#include "GPhoCS.h"
#include "DataLayer.h"
#include "McRefCommon.h"
#include "MemoryMng.h"
#include <unordered_set>

double *tau_bounds;
int *lca_pops;


void calculateTauBounds1() {
  for (int gen = 0; gen < dataSetup.numLoci; gen++) {
    initializeLcaPops(gen);
    int rootNodeId = dataState.lociData[gen]->root;
    calculateLociTauBounds(rootNodeId, gen);
  }
}

void initializeLcaPops(int gen) {

  for (int nodeId = 0; nodeId < numNodes(); nodeId++) {

    LikelihoodNode *currentNode = getNode(nodeId, gen);

    if (isLeafNode(currentNode)) {
      lca_pops[nodeId] = getLeafNodePop(nodeId, gen);
    } else {
      lca_pops[nodeId] = -1;
    }
  }
}

int getLeafNodePop(int nodeId, int gen) {
  return nodePops[gen][nodeId]; // TODO - implement
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

void updateTauBoundsOfDescendants(int pop, double bound) {
  if (isLeaf(pop)) return;

  tau_bounds[pop] = fmin(tau_bounds[pop], bound);

  updateTauBoundsOfDescendants(getSon(pop, LEFT), bound);
  updateTauBoundsOfDescendants(getSon(pop, RIGHT), bound);
}

LikelihoodNode *getNode(int nodeId, int gen) { // TODO - move to util
  LocusData *loci = dataState.lociData[gen];
  return loci->nodeArray[nodeId];
}

bool isLeafNode(LikelihoodNode *node) {
  return (node->leftSon == -1) || (node->rightSon == -1);  //TOASK - make sure this is correct
}

void calculateTauBounds2() {
  map<int, int> eventToPopMap;

  for (int gen = 0; gen < dataSetup.numLoci; gen++) {
    for (int currentPop = 0; currentPop < dataSetup.popTree->numPops; ++currentPop) {
      if (isLeaf(currentPop)) continue;
      for (int eventIdx = event_chains[gen].first_event[currentPop];
           eventIdx >= 0;
           eventIdx = event_chains[gen].events[eventIdx].getNextIdx()) {
        int currentEventId = event_chains[gen].events[eventIdx].getId();

        eventToPopMap.insert(pair<int, int>(currentEventId, currentPop));
      }
    }

    for (int currentPop = 0; currentPop < dataSetup.popTree->numPops; ++currentPop) {
      if (isLeaf(currentPop)) continue;
      for (int eventIdx = event_chains[gen].first_event[currentPop];
           eventIdx >= 0;
           eventIdx = event_chains[gen].events[eventIdx].getNextIdx()) {
        Event &currentEvent = event_chains[gen].events[eventIdx];
        int currentEventId = currentEvent.getId();
        EventType eventType = currentEvent.getType();
        LocusData *currentLocus = dataState.lociData[gen];

        if (eventType == COAL) {
          int leftSonId = getNodeSon(currentLocus, currentEventId, LEFT);
          int leftSonPop = eventToPopMap[leftSonId];
          int rightSonId = getNodeSon(currentLocus, currentEventId, RIGHT);
          int rightSonPop = eventToPopMap[rightSonId];
          if (leftSonPop != rightSonPop && rightSonPop != currentPop && leftSonPop != currentPop) {
            double eventAge = getNodeAge(currentLocus, currentEventId);
            tau_bounds[currentPop] = min(eventAge, tau_bounds[currentPop]);
          }
        }
      }
    }
    eventToPopMap.clear();
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
    if (isLeaf(pop)) {
      tau_bounds[pop] = 0.0;
    } else {
      tau_bounds[pop] = 100.0; // TODO - init with proper value (tau of father pop? MAX/INF value?)
    }
  }
}

int numNodes() { return (2 * dataSetup.numSamples) - 1; } // TODO - make MACRO?
