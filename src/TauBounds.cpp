#include <map>
#include "TauBounds.h"
#include "patch.h"
#include "GPhoCS.h"
#include "DataLayer.h"
#include "McRefCommon.h"

float *current_tau_bounds;

void calculateTauBounds() {
  printf("Calculating Tau Bounds...\n");
  int nLociIdx = 0;
  int nPopIdx;
  int nEventIdx;
  map<int, int> CoalIDtoPopIdx;

  for (nPopIdx = 0; nPopIdx < dataSetup.popTree->numPops; ++nPopIdx) {
    for (nEventIdx = event_chains[nLociIdx].first_event[nPopIdx];
         nEventIdx >= 0;
         nEventIdx = event_chains[nLociIdx].events[nEventIdx].getNextIdx()) {
      Event &CurrEvent = event_chains[nLociIdx].events[nEventIdx];
      int nEventID = CurrEvent.getId();
      EventType eEventType = CurrEvent.getType();
      if (COAL == eEventType || END_CHAIN == eEventType)
        CoalIDtoPopIdx.insert(pair<int, int>(nEventID, nPopIdx));
    }
  }

  for (nPopIdx = 0; nPopIdx < dataSetup.popTree->numPops; ++nPopIdx) {
    for (nEventIdx = event_chains[nLociIdx].first_event[nPopIdx];
         nEventIdx >= 0;
         nEventIdx = event_chains[nLociIdx].events[nEventIdx].getNextIdx()) {
      Event &CurrEvent = event_chains[nLociIdx].events[nEventIdx];
      int nEventID = CurrEvent.getId();
      printf("[Pop %2d] %2d %s ", nPopIdx, nEventID, getEventTypeName(CurrEvent.getType()));
      if (COAL == event_chains[nLociIdx].events[nEventIdx].getType()) {
        int nParentID = getNodeFather(dataState.lociData[nLociIdx], nEventID);
        int nLeftSonID = getNodeSon(dataState.lociData[nLociIdx], nEventID, 0);
        int nLeftSonPopIdx = CoalIDtoPopIdx[nLeftSonID];
        int nRightSonID = getNodeSon(dataState.lociData[nLociIdx], nEventID, 1);
        int nRightSonPopIdx = CoalIDtoPopIdx[nRightSonID];
        printf("COAL father = %d, sons = %d [%d], %d [%d]",
               nParentID, nLeftSonID, nLeftSonPopIdx,
               nRightSonID, nRightSonPopIdx);
      }
      puts("");
    }
    puts("");
  }
}

void printTauBoundsHeader(FILE *file) {
  fprintf(file, "iteration");
  for (int pop = 0; pop < dataSetup.popTree->numPops; pop++) {
    fprintf(file, "\t%s", getPopName(pop));
  }
  fprintf(file, "\n");
}

void printTauBounds(int iteration, FILE *file) {
  fprintf(file, "%d", iteration);
  for (int pop = 0; pop < dataSetup.popTree->numPops; pop++) {
    fprintf(file, "\t%.10f", current_tau_bounds[pop]);
  }
  fprintf(file, "\n");
}

void allocateTauBoundsMem() {
  current_tau_bounds = (float *) malloc(dataSetup.popTree->numPops * sizeof(float));
}
