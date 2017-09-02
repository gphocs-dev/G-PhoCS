#include <map>
#include "TauBounds.h"
#include "patch.h"
#include "GPhoCS.h"
#include "DataLayer.h"
#include "McRefCommon.h"

double *tau_bounds;

void calculateTauBounds() {
  int nLociIdx = 0;
  int currentPop;
  int nEventIdx;
  map<int, int> CoalIDtoPopIdx;

  for (currentPop = 0; currentPop < dataSetup.popTree->numPops; ++currentPop) { // for every population,
    for (nEventIdx = event_chains[nLociIdx].first_event[currentPop];
         nEventIdx >= 0;
         nEventIdx = event_chains[nLociIdx].events[nEventIdx].getNextIdx()) { // run over it's event_chain
      Event &CurrEvent = event_chains[nLociIdx].events[nEventIdx];
      int nEventID = CurrEvent.getId(); // get THE event id,
      EventType eEventType = CurrEvent.getType();
      if (COAL == eEventType || END_CHAIN == eEventType)
        CoalIDtoPopIdx.insert(pair<int, int>(nEventID, currentPop)); // and remember for every COAL or END_CHAIN event id, what pop it occurred in.
    }
  }

  for (currentPop = 0; currentPop < dataSetup.popTree->numPops; ++currentPop) { // Now, for every population,
    for (nEventIdx = event_chains[nLociIdx].first_event[currentPop];
         nEventIdx >= 0;
         nEventIdx = event_chains[nLociIdx].events[nEventIdx].getNextIdx()) {  // run over it's event_chain,
      Event &currentEvent = event_chains[nLociIdx].events[nEventIdx];

      int currentEventID = currentEvent.getId();
//      printf("[Pop %2d] %2d %s ,", currentPop, currentEventID, getEventTypeName(currentEvent.getType())); // print curr event+pop
      if (COAL == currentEvent.getType()) { // and if event is COAL,
        int leftSonID = getNodeSon(dataState.lociData[nLociIdx], currentEventID, LEFT);
        int leftSonPop = CoalIDtoPopIdx[leftSonID];
        int rightSonID = getNodeSon(dataState.lociData[nLociIdx], currentEventID, RIGHT);
        int rightSonPop = CoalIDtoPopIdx[rightSonID];
        if (leftSonPop != rightSonPop && rightSonPop != currentPop && leftSonPop != currentPop ){
          double eventAge = getNodeAge(dataState.lociData[nLociIdx], currentEventID);
          tau_bounds[currentPop] = min(eventAge, tau_bounds[currentPop]);
        }
//        printf("COAL sons = %d [%d], %d [%d]", parentID, leftSonID, leftSonPop, rightSonID, rightSonPop);
      }
//      puts("");
    }
//    puts("");
  }
}

void printTauBoundsHeader(FILE *file) {
  fprintf(file, "iteration");
  for (int pop = 0; pop < dataSetup.popTree->numPops; pop++) {
    fprintf(file, "\t%s_bound\t%s_tau", getPopName(pop), getPopName(pop));
  }
  fprintf(file, "\n");
}

void printTauBounds(int iteration, FILE *file) {
  fprintf(file, "%d", iteration);
  for (int pop = 0; pop < dataSetup.popTree->numPops; pop++) {
    fprintf(file, "\t%.40f\t%.40f", tau_bounds[pop], dataSetup.popTree->pops[pop]->age);
  }
  fprintf(file, "\n");
}

void allocateTauBoundsMem() {
  tau_bounds = (double *) malloc(dataSetup.popTree->numPops * sizeof(double));
}

void initializeTauBounds() {
  for (int pop = 0; pop < dataSetup.popTree->numPops; pop++) {
    tau_bounds[pop] = 100.0;  // TODO - init with proper value (tau of father pop? MAX/INF value?)
  }
}
