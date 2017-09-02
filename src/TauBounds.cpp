#include <map>
#include "TauBounds.h"
#include "patch.h"
#include "GPhoCS.h"
#include "DataLayer.h"
#include "McRefCommon.h"

double *tau_bounds;

void calculateTauBounds() {
  map<int, int> eventToPopMap;

  for (int gen = 0; gen < dataSetup.numLoci; gen++) {
    for (int currentPop = 0; currentPop < dataSetup.popTree->numPops; ++currentPop) {
      if (isLeaf(currentPop)) continue;
      for (int eventIdx = event_chains[gen].first_event[currentPop];
           eventIdx >= 0;
           eventIdx = event_chains[gen].events[eventIdx].getNextIdx()) {
        Event &currentEvent = event_chains[gen].events[eventIdx];
        int eventId = currentEvent.getId();
        EventType eventType = currentEvent.getType();
        if (eventType == COAL || eventType == END_CHAIN)
          eventToPopMap.insert(pair<int, int>(eventId, currentPop));
      }
    }

    for (int currentPop = 0; currentPop < dataSetup.popTree->numPops; ++currentPop) {
      if (isLeaf(currentPop)) continue;
      for (int eventIdx = event_chains[gen].first_event[currentPop];
           eventIdx >= 0;
           eventIdx = event_chains[gen].events[eventIdx].getNextIdx()) {
        Event &currentEvent = event_chains[gen].events[eventIdx];

        int currentEventId = currentEvent.getId();
        if (COAL == currentEvent.getType()) { // and if event is COAL,
          int leftSonId = getNodeSon(dataState.lociData[gen], currentEventId, LEFT);
          int leftSonPop = eventToPopMap[leftSonId];
          int rightSonId = getNodeSon(dataState.lociData[gen], currentEventId, RIGHT);
          int rightSonPop = eventToPopMap[rightSonId];
          if (leftSonPop != rightSonPop && rightSonPop != currentPop && leftSonPop != currentPop) {
            double eventAge = getNodeAge(dataState.lociData[gen], currentEventId);
            tau_bounds[currentPop] = min(eventAge, tau_bounds[currentPop]);
          }
        }
      }
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
}

void initializeTauBounds() {
  tau_bounds[dataSetup.popTree->rootPop] = 100.0;  // TODO - init with proper value (tau of father pop? MAX/INF value?)
  for (int fatherId = 0; fatherId < dataSetup.popTree->numPops; fatherId++) {
    if (isLeaf(fatherId)) continue;
    Population *father = dataSetup.popTree->pops[fatherId];
    double tau = father->age;

    int leftSonId = father->sons[LEFT]->id;
    int rightSonId = father->sons[RIGHT]->id;

    // The maximum allowed tau, without any restricting coalescence event, is tau of father population
    tau_bounds[leftSonId] = tau;
    tau_bounds[rightSonId] = tau;
  }
}
