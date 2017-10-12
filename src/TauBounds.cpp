#include <map>
#include "LocusDataLikelihood.h"
#include "TauBounds.h"
#include "patch.h"
#include "GPhoCS.h"
#include "DataLayer.h"
#include "McRefCommon.h"

double *tau_bounds;


void calculateTauBounds1() {
  for (int gen = 0; gen < dataSetup.numLoci; gen++) {

    int rootNodeId = dataState.lociData[gen]->root;
    calculateLociTauBounds(rootNodeId, gen);
  }
}

void calculateLociTauBounds(int nodeId, int gen) {
  //TOASK - how do I get from LikelihoodNode to its population id?
  //TOASK - or even, from a leaf LikelihoodNode if it's easier

  printf("calculating tauBounds on node %d in gen %d\n", nodeId, gen);

  LikelihoodNode *currentNode = getNode(nodeId, gen);

  if (isLeafNode(currentNode)) {
    return;
  }

  calculateLociTauBounds(currentNode->leftSon, gen);
  calculateLociTauBounds(currentNode->rightSon, gen);

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
