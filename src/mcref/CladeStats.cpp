#include <stdlib.h>
#include <string.h>

#include "../utils.h"
#include "../MCMCcontrol.h"
#include "../patch.h"
#include "../DataLayer.h"
#include "../MemoryMng.h"
#include "McRefCommon.h"
#include "CladeStats.h"

CLADE_STATS *clade_stats;

void calculateCladeStats() {
  initCladeStats();
  computeCladeNumMigs();
  computeCladeMigStats();
  computeCladeNumCoals();
  computeCladeCoalStats();
}

void computeCladeNumMigs() {
  for (int mig_band = 0; mig_band < dataSetup.popTree->numMigBands; mig_band++) {
    clade_stats[mig_band].num_migs_total = genetree_stats_total.num_migs[mig_band];
  }
}

void computeCladeMigStats() {
  for (int mig_band = 0; mig_band < dataSetup.popTree->numMigBands; mig_band++) {
    clade_stats[mig_band].mig_stats_total = genetree_stats_total.mig_stats[mig_band];
  }
}

void computeCladeNumCoals() {
  computeCladeNumCoals_rec(dataSetup.popTree->rootPop);
}

void computeCladeNumCoals_rec(int pop) {
  int leftSon, rightSon;

  if (isLeafPop(pop)) {
    clade_stats[pop].num_coals_total = genetree_stats_total.num_coals[pop];
  } else {
    leftSon = dataSetup.popTree->pops[pop]->sons[LEFT]->id;
    rightSon = dataSetup.popTree->pops[pop]->sons[RIGHT]->id;

    computeCladeNumCoals_rec(leftSon);
    computeCladeNumCoals_rec(rightSon);

    clade_stats[pop].num_coals_total = genetree_stats_total.num_coals[pop] +
                                       clade_stats[leftSon].num_coals_total + clade_stats[rightSon].num_coals_total;
  }
}


void computeCladeCoalStats() {
  for (int gen = 0; gen < dataSetup.numLoci; gen++) {
    computeCladeCoalStats_rec(dataSetup.popTree->rootPop, gen);
  }
}

void computeCladeCoalStats_rec(int clade, int gen) {
  int leftSon, rightSon;

  if (isLeafPop(clade)) {
    fillUpLeafCladeStats(clade, gen);
  } else {
    leftSon = dataSetup.popTree->pops[clade]->sons[LEFT]->id;
    rightSon = dataSetup.popTree->pops[clade]->sons[RIGHT]->id;

    computeCladeCoalStats_rec(leftSon, gen);
    computeCladeCoalStats_rec(rightSon, gen);

    fillUpCladeStats(clade, gen);

  }
}

void fillUpLeafCladeStats(int clade, int gen) {
  appendPopToClade(clade, gen, 0); // since this is a leaf, starting point of the clade_stats arrays should be 0
}

void appendPopToClade(int clade, int gen, int startingPoint) {
  int i = startingPoint;
  int event = event_chains[gen].first_event[clade];
  double cladeStartTime = dataSetup.popTree->pops[clade]->age;
  double eventAge = cladeStartTime;

  for (; event >= 0; i++, event = event_chains[gen].events[event].getNextIdx()) {
    eventAge += event_chains[gen].events[event].getElapsedTime();  // TODO - replace with getNodeAge(dataState.lociData[nLociIdx], currentEventID);
    clade_stats[clade].sorted_ages[i] = eventAge;
    clade_stats[clade].elapsed_times[i] = event_chains[gen].events[event].getElapsedTime();
    clade_stats[clade].num_lineages[i] = event_chains[gen].events[event].getNumLineages();
    clade_stats[clade].event_types[i] = event_chains[gen].events[event].getType();
  }
  clade_stats[clade].num_events = i;
//	clade_stats[clade].coal_stats_total += genetree_stats[gen].coal_stats[clade];
  clade_stats[clade].coal_stats_total += calculateCoalStats(clade_stats[clade].elapsed_times + startingPoint, clade_stats[clade].num_lineages + startingPoint,
                                                            clade_stats[clade].num_events - startingPoint);
}

void fillUpCladeStats(int clade, int gen) {
  addChildrenIntoCladeStats(clade);
  addCurrentPopIntoCladeStats(clade, gen);
}

void addChildrenIntoCladeStats(int clade) {
  mergeChildren(clade);
  addChildrenCladeStats(clade);
}


void mergeChildren(int clade) {
  int i, j = 0, k = 0;
  double leftAge, rightAge;
  int leftSon, rightSon;
  leftSon = dataSetup.popTree->pops[clade]->sons[LEFT]->id;
  rightSon = dataSetup.popTree->pops[clade]->sons[RIGHT]->id;
  int m = clade_stats[leftSon].num_events;
  int n = clade_stats[rightSon].num_events;

  for (i = 0; i < m + n;) {
    if (j < m && k < n) {
      leftAge = clade_stats[leftSon].sorted_ages[j];
      rightAge = clade_stats[rightSon].sorted_ages[k];
      if (leftAge < rightAge) {
        clade_stats[clade].event_types[i] = clade_stats[leftSon].event_types[j];
        clade_stats[clade].sorted_ages[i] = clade_stats[leftSon].sorted_ages[j];
        clade_stats[clade].num_lineages[i] =
               clade_stats[leftSon].num_lineages[j] + clade_stats[rightSon].num_lineages[k];
        j++;
      } else {
        clade_stats[clade].event_types[i] = clade_stats[rightSon].event_types[k];
        clade_stats[clade].sorted_ages[i] = clade_stats[rightSon].sorted_ages[k];
        clade_stats[clade].num_lineages[i] =
               clade_stats[leftSon].num_lineages[j] + clade_stats[rightSon].num_lineages[k];
        k++;
      }
      i++;
    } else if (j == m) {
      for (; i < m + n;) {
        clade_stats[clade].event_types[i] = clade_stats[rightSon].event_types[k];
        clade_stats[clade].sorted_ages[i] = clade_stats[rightSon].sorted_ages[k];
        clade_stats[clade].num_lineages[i] =
               clade_stats[leftSon].num_lineages[j - 1] + clade_stats[rightSon].num_lineages[k];
        k++;
        i++;
      }
    } else {
      for (; i < m + n;) {
        clade_stats[clade].event_types[i] = clade_stats[leftSon].event_types[j];
        clade_stats[clade].sorted_ages[i] = clade_stats[leftSon].sorted_ages[j];
        clade_stats[clade].num_lineages[i] =
               clade_stats[leftSon].num_lineages[j] + clade_stats[rightSon].num_lineages[k - 1];
        j++;
        i++;
      }
    }
  }
  for (i = 0; i < m + n; i++) {
    clade_stats[clade].elapsed_times[i] = (i == 0) ? 0 :
                                          clade_stats[clade].sorted_ages[i] - clade_stats[clade].sorted_ages[i - 1];
  }

  clade_stats[clade].num_events = m + n;
}

void addChildrenCladeStats(int clade) {
  clade_stats[clade].coal_stats_total +=
         calculateCoalStats(clade_stats[clade].elapsed_times, clade_stats[clade].num_lineages, clade_stats[clade].num_events);
}

void addCurrentPopIntoCladeStats(int clade, int gen) {
  appendPopToClade(clade, gen, clade_stats[clade].num_events); // start filling the clade_stats arrays from the last known event
}


void allocateCladeMem() {
  int max_events = 2 * dataSetup.numSamples + 4 * MAX_MIGS + 3 * dataSetup.popTree->numMigBands + dataSetup.popTree->numPops + 10;
  clade_stats = (CLADE_STATS *) malloc(dataSetup.popTree->numPops * sizeof(CLADE_STATS));

  for (int clade = 0; clade < dataSetup.popTree->numPops; clade++) {
    clade_stats[clade].sorted_ages = (double *) malloc(max_events * sizeof(double));
    clade_stats[clade].elapsed_times = (double *) malloc(max_events * sizeof(double));
    clade_stats[clade].num_lineages = (int *) malloc(max_events * sizeof(int));
    clade_stats[clade].event_types = (int *) malloc(max_events * sizeof(int));
    clade_stats[clade].coal_stats_total = 0.0;

    clade_stats[clade].mig_stats_total = 0.0;


    if ((clade_stats == NULL) ||
        (!clade_stats[clade].sorted_ages) || (!clade_stats[clade].elapsed_times) ||
        (!clade_stats[clade].num_lineages) || (!clade_stats[clade].event_types)) {
      fprintf(stderr, "\nError: Out Of Memory clade_stats\n");
      exit(-1);
    }
  }
}

/***
 * cleans up all variables in clade_stats
 */
void initCladeStats() {
  for (int clade = 0; clade < dataSetup.popTree->numPops; clade++) {
    initSpecificCladeStats(clade);
  }
}

void initSpecificCladeStats(int clade) {
  clade_stats[clade].coal_stats_total = 0.0;
  clade_stats[clade].mig_stats_total = 0.0;
  clade_stats[clade].num_coals_total = 0;
  clade_stats[clade].num_migs_total = 0;
  clade_stats[clade].num_events = 0;
}

void freeCladeMem() {
}

