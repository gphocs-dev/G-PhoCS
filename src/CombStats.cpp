#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "utils.h"
#include "MCMCcontrol.h"
#include "patch.h"
#include "DataLayer.h"
#include "MemoryMng.h"
#include "CombStats.h"
#include "McRefCommon.h"
#include "CombAssertions.h"

COMB_STATS *comb_stats;

void calculateCombStats() {
  initCombStats();
  for (int comb = 0; comb < dataSetup.popTree->numPops; comb++) {
    if (isFeasibleComb(comb)) {
      for (int gene = 0; gene < dataSetup.numLoci; gene++) {
        calculateSufficientStats(comb, gene);
      }
    }
  }
  if (DEBUG_COMB_STATS) runAssertions();
}

void calculateSufficientStats(int comb, int gene) {
  coalescence(comb, gene);
  migrations(comb, gene);
  finalizeCombCoalStats(comb);
}

void coalescence(int comb, int gene) {
  coalescence_rec(comb, comb, gene);
}

void coalescence_rec(int comb, int currentPop, int gene) {
  if (isLeaf(currentPop)) {
    handleLeafCoals(comb, currentPop, gene);
  } else {
    coalescence_rec(comb, getSon(currentPop, LEFT), gene);
    coalescence_rec(comb, getSon(currentPop, RIGHT), gene);

    handleNonLeafCoals(comb, currentPop, gene);
  }
}

void handleLeafCoals(int comb, int leaf, int gene) {
  double previousAge, eventAge = 0.0;
  int eventId;
  double combAge = comb_stats[comb].age;

  EventChain chain = event_chains[gene];
  Event event;
  Stats *belowCombLeafStats = &comb_stats[comb].leaves[leaf].below_comb;
  Stats *aboveCombLeafStats = &comb_stats[comb].leaves[leaf].above_comb;
  Stats *combTotalStats = &comb_stats[comb].total;

  eventId = chain.first_event[leaf];
  while (eventId >= 0) {
    event = chain.events[eventId];
    previousAge = eventAge;
    eventAge += event.getElapsedTime();

    if (isEventCompletelyBelowComb(eventAge, combAge)) {
      countCoalEventTowardsBelowComb(event, belowCombLeafStats);
    } else if (isBorderEvent(eventAge, previousAge, combAge)) {
      countCoalEventTowardsHalfAndHalf(event, eventAge, previousAge, combAge, belowCombLeafStats, aboveCombLeafStats, combTotalStats);
    } else if (isEventCompletelyInsideComb(eventAge, combAge)) {
      countCoalEventTowardsAboveComb(event, eventAge, aboveCombLeafStats, combTotalStats);
    } else {
      printErrorAndExit("leaf event wasn't counted towards any Stats");
    }
    eventId = event.getNextIdx();
  }
  aboveCombLeafStats->coal_stats += calculateCoalStats(aboveCombLeafStats->elapsed_times, aboveCombLeafStats->num_lineages, aboveCombLeafStats->num_events);
  combTotalStats->num_events += aboveCombLeafStats->num_events;
}

void countCoalEventTowardsBelowComb(Event event, Stats *belowCombLeafStats) {
  if (event.getType() == COAL) belowCombLeafStats->num_coals++;
  double elapsedTime = event.getElapsedTime();
  int numLins = event.getNumLineages();
  belowCombLeafStats->coal_stats += (numLins) * (numLins - 1) * elapsedTime;
}

void
countCoalEventTowardsHalfAndHalf(Event event, double eventAge, double previousAge, double combAge, Stats *belowCombLeafStats, Stats *aboveCombLeafStats, Stats *combTotalStats) {
  double pseudoEventAge, pseudoElapsedTimeBelow, pseudoElapsedTimeAbove;
  EventTypeIvgeny pseudoEventType;
  int eventId = event.getId();
  int numLins = event.getNumLineages();
  EventTypeIvgeny eventType = event.getType();


  if (areAlmostEqual(eventAge, combAge)) {
    if (eventType == COAL) belowCombLeafStats->num_coals++;
    belowCombLeafStats->coal_stats += (numLins) * (numLins - 1) * (event.getElapsedTime());
    pseudoElapsedTimeAbove = 0.0;
    pseudoEventType = DUMMY;
  } else { // the border event is above combAge so we need to "split it" into two events

    if (eventType == COAL) {
      aboveCombLeafStats->num_coals++;
      combTotalStats->num_coals++;
    }
    pseudoEventAge = combAge;
    pseudoElapsedTimeBelow = pseudoEventAge - previousAge;
    belowCombLeafStats->coal_stats += (numLins) * (numLins - 1) * (pseudoElapsedTimeBelow);

    pseudoElapsedTimeAbove = eventAge - pseudoEventAge;
    pseudoEventType = eventType;
  }

  aboveCombLeafStats->elapsed_times[0] = pseudoElapsedTimeAbove;
  aboveCombLeafStats->sorted_ages[0] = eventAge;
  aboveCombLeafStats->num_lineages[0] = numLins;
  aboveCombLeafStats->event_types[0] = pseudoEventType;
  aboveCombLeafStats->event_ids[0] = eventId;
  aboveCombLeafStats->num_events = 1;
}

void countCoalEventTowardsAboveComb(Event event, double eventAge, Stats *aboveCombLeafStats, Stats *combTotalStats) {
  EventTypeIvgeny eventType = event.getType();
  int numEventsAboveComb = aboveCombLeafStats->num_events;

  if (eventType == COAL) {
    aboveCombLeafStats->num_coals++;
    combTotalStats->num_coals++;
  }
  // you're a tiny part of the comb so supply future methods the data they need -
  aboveCombLeafStats->sorted_ages[numEventsAboveComb] = eventAge;
  aboveCombLeafStats->elapsed_times[numEventsAboveComb] = event.getElapsedTime();
  aboveCombLeafStats->num_lineages[numEventsAboveComb] = event.getNumLineages();
  aboveCombLeafStats->event_types[numEventsAboveComb] = event.getType();
  aboveCombLeafStats->event_ids[numEventsAboveComb] = event.getId();

  aboveCombLeafStats->num_events += 1;
}

void handleNonLeafCoals(int comb, int currentPop, int gene) {
  handleNonLeafNumCoals(comb, currentPop, gene);
  handleNonLeafCoalStats(comb, currentPop, gene);
}

void handleNonLeafNumCoals(int comb, int currentPop, int gene) {
  comb_stats[comb].total.num_coals += genetree_stats[gene].num_coals[currentPop];
}

void handleNonLeafCoalStats(int comb, int currentPop, int gene) {
  mergeChildrenIntoCurrent(comb, currentPop);
  appendCurrent(comb, currentPop, gene); // start filling the comb_stats arrays from the last known event
}

void mergeChildrenIntoCurrent(int comb, int currentPop) {
  int i, j = 0, k = 0;
  double leftAge, rightAge;
  Stats *leftStats, *rightStats, *currentStats;

  currentStats = getCombPopStats(comb, currentPop);
  leftStats = getCombPopStats(comb, getSon(currentPop, LEFT));
  rightStats = getCombPopStats(comb, getSon(currentPop, RIGHT));


  int m = leftStats->num_events;
  int n = rightStats->num_events;

  for (i = 0; i < m + n; i++) {
    if (m == 0) { // there are no events in left son, so choose the event from the right son
      copyStaticEventStats(rightStats, k, currentStats, i);
      k++;
    } else if (n == 0) { // there are no events in left son, so choose the event from the left son
      copyStaticEventStats(leftStats, j, currentStats, i);
      j++;
    } else if (j < m && k < n) { // both sons have more events
      currentStats->num_lineages[i] = leftStats->num_lineages[j] + rightStats->num_lineages[k];

      leftAge = leftStats->sorted_ages[j];
      rightAge = rightStats->sorted_ages[k];
      if (leftAge < rightAge) { // since left event is more recent, choose it
        copyStaticEventStats(leftStats, j, currentStats, i);
        j++;
      } else {  // since right event is more recent, choose it
        copyStaticEventStats(rightStats, k, currentStats, i);
        k++;
      }
    } else if (j == m) { // all events in left son were handled
      copyStaticEventStats(rightStats, k, currentStats, i);
      currentStats->num_lineages[i] = leftStats->num_lineages[j - 1] + rightStats->num_lineages[k];
      k++;
    } else if (k == n) { // all events in right son were handled
      copyStaticEventStats(leftStats, j, currentStats, i);
      currentStats->num_lineages[i] = leftStats->num_lineages[j] + rightStats->num_lineages[k - 1];
      j++;
    }
  }
  for (i = 0; i < m + n; i++) {
    currentStats->elapsed_times[i] = (i == 0) ? 0 :
                                     currentStats->sorted_ages[i] - currentStats->sorted_ages[i - 1];
  }

  currentStats->num_events = m + n;
}

void copyStaticEventStats(Stats *sourceStats, int n, Stats *targetStats, int m) {
  targetStats->event_ids[m] = sourceStats->event_ids[n];
  targetStats->event_types[m] = sourceStats->event_types[n];
  targetStats->sorted_ages[m] = sourceStats->sorted_ages[n];
}

void appendCurrent(int comb, int currentPop, int gene) {
  Stats *currentStats = getCombPopStats(comb, currentPop);
  EventChain chain = event_chains[gene];

  int startingPoint = currentStats->num_events;
  int i = startingPoint;
  int eventId = chain.first_event[currentPop];
  double startTime = dataSetup.popTree->pops[currentPop]->age; // do I need to start from pop age or from last eventId age (or are they equal)?
  double eventAge = startTime;

  for (; eventId >= 0; i++, eventId = chain.events[eventId].getNextIdx()) {
    eventAge += chain.events[eventId].getElapsedTime();
    currentStats->sorted_ages[i] = eventAge;
    currentStats->elapsed_times[i] = chain.events[eventId].getElapsedTime();
    currentStats->num_lineages[i] = chain.events[eventId].getNumLineages();
    currentStats->event_types[i] = chain.events[eventId].getType();
    currentStats->event_ids[i] = eventId;
  }
  currentStats->num_events = i;
}

void finalizeCombCoalStats(int comb) {
  double *elapsedTimes = comb_stats[comb].combs[comb].elapsed_times;
  int *numLineages = comb_stats[comb].combs[comb].num_lineages;
  int size = comb_stats[comb].combs[comb].num_events;
  comb_stats[comb].total.coal_stats += calculateCoalStats(elapsedTimes, numLineages, size);
}

void migrations(int comb, int gene) {
  for (int migband = 0; migband < dataSetup.popTree->numMigBands; migband++) {
    if (isMigOfComb(migband, comb)) {
      if (isCombLeafMigBand(migband, comb)) {
        handleCombLeavesMigBand(comb, migband, gene);
      }
      if (isMigBandInternal(migband, comb)) {
        // ignore internal migbands. their stats aren't used
      }
      if (isMigBandExternal(migband, comb)) {
      }
    }
  }
}


void handleCombLeavesMigBand(int comb, int mig, int gene) {
  double previousAge, eventAge = 0.0;
  double combAge = comb_stats[comb].age;

  EventChain chain = event_chains[gene];
  Event event;

  MigStats *leafMigStats = comb_stats[comb].leafMigs;

  int target = getTargetPop(mig);
  for (int eventId = chain.first_event[target]; eventId > 0; eventId = chain.events[eventId].getNextIdx()) {
    event = chain.events[eventId];
    previousAge = eventAge;
    eventAge += event.getElapsedTime();

    if (isEventCompletelyBelowComb(eventAge, combAge)) {
      countMigEventTowardsBelowComb(event, leafMigStats);
    } else if (isBorderEvent(eventAge, previousAge, combAge)) {
      countMigEventTowardsHalfAndHalf(event, eventAge, previousAge, combAge, leafMigStats);
    } else if (isEventCompletelyInsideComb(eventAge, combAge)) {
      countMigEventTowardsAboveComb(event, leafMigStats);
    } else {
      printErrorAndExit("leaf event wasn't counted towards any Stats");
    };
    if (event.getType() == MIG_BAND_END) break;
  }
}

void countMigEventTowardsBelowComb(Event event, MigStats *leafMigStats) {
  if (event.getType() == IN_MIG) leafMigStats->num_migs++;
  double elapsedTime = event.getElapsedTime();
  int numLins = event.getNumLineages();
  leafMigStats->mig_stats += (numLins) * elapsedTime;
}

void countMigEventTowardsHalfAndHalf(Event event, double eventAge, double previousAge, double combAge, MigStats *leafMigStats) {
  double pseudoEventAge, pseudoElapsedTimeBelow, pseudoElapsedTimeAbove;
  int numLins = event.getNumLineages();
  EventTypeIvgeny eventType = event.getType();

  if (areAlmostEqual(eventAge, combAge)) {
    if (eventType == IN_MIG) leafMigStats->num_migs++;
    leafMigStats->mig_stats += numLins * (event.getElapsedTime());
  } else { // the border event is above combAge so we need to "split it" into two events
    pseudoEventAge = combAge;
    pseudoElapsedTimeBelow = pseudoEventAge - previousAge;
    leafMigStats->mig_stats += numLins * pseudoElapsedTimeBelow;

    if (eventType == IN_MIG) leafMigStats->num_migs_above++;
    pseudoElapsedTimeAbove = eventAge - pseudoEventAge;
    leafMigStats->mig_stats_above += numLins * pseudoElapsedTimeAbove;
  }
}

void countMigEventTowardsAboveComb(Event event, MigStats *leafMigStats) {

  if (event.getType() == IN_MIG) leafMigStats->num_migs_above++;
  double elapsedTime = event.getElapsedTime();
  int numLins = event.getNumLineages();
  leafMigStats->mig_stats_above += numLins * elapsedTime;
}

Stats *getCombPopStats(int comb, int pop) {
  if (isLeaf(pop)) {
    return &comb_stats[comb].leaves[pop].above_comb;
  } else {
    return &comb_stats[comb].combs[pop];
  }
}

void initCombStats() {
  initPopStats();
  initMigStats();
}

void initPopStats() {
  for (int comb = 0; comb < dataSetup.popTree->numPops; comb++) {
    if (!isLeaf(comb)) {
      initStats(&comb_stats[comb].total);
      comb_stats[comb].age = getCombAge(comb);
      for (int pop = 0; pop < dataSetup.popTree->numPops; pop++) {
        if (isLeaf(pop)) {
          initStats(&comb_stats[comb].leaves[pop].above_comb);
          initStats(&comb_stats[comb].leaves[pop].below_comb);
        } else {
          initStats(&comb_stats[comb].combs[pop]);
        }
      }
    }
  }
}

void initStats(Stats *stats) {
  stats->coal_stats = 0.0;
  stats->num_coals = 0;
  stats->num_events = 0;
}

void initMigStats() {
  for (int comb = 0; comb < dataSetup.popTree->numPops; comb++) {
    if (isFeasibleComb(comb)) {
      for (int mig = 0; mig < dataSetup.popTree->numMigBands; mig++) {
        if (isMigOfComb(mig, comb)) {
          comb_stats[comb].leafMigs[mig].mig_stats = 0.0;
          comb_stats[comb].leafMigs[mig].num_migs = 0;
          comb_stats[comb].leafMigs[mig].mig_stats_above = 0.0;
          comb_stats[comb].leafMigs[mig].num_migs_above = 0;
        }
      }
    }
  }
}

void allocateCombMem() {
  comb_stats = (COMB_STATS *) malloc(dataSetup.popTree->numPops * sizeof(COMB_STATS));

  allocatePopsMem();
  allocateMigBandsMem();
}

void allocatePopsMem() {
  for (int comb = 0; comb < dataSetup.popTree->numPops; comb++) {
    allocateStats(&comb_stats[comb].total);
    comb_stats[comb].leaves = (LeafStats *) malloc(dataSetup.popTree->numPops * sizeof(LeafStats));
    comb_stats[comb].combs = (Stats *) malloc(dataSetup.popTree->numPops * sizeof(Stats));
    for (int pop = 0; pop < dataSetup.popTree->numPops; pop++) {
      if (isLeaf(pop)) {
        allocateStats(&comb_stats[comb].leaves[pop].below_comb);
        allocateStats(&comb_stats[comb].leaves[pop].above_comb);
      } else {
        allocateStats(&comb_stats[comb].combs[pop]);
      }
    }
    if (comb_stats == NULL) {
      // TODO - add memory allocation test for all of comb_stats
      fprintf(stderr, "\nError: Out Of Memory comb_stats\n");
      exit(-1);
    }
  }
}

void allocateStats(Stats *stats) {
  int max_events = 2 * dataSetup.numSamples + 4 * MAX_MIGS + 3 * dataSetup.popTree->numMigBands + dataSetup.popTree->numPops + 10;
  stats->sorted_ages = (double *) malloc(max_events * sizeof(double));
  stats->elapsed_times = (double *) malloc(max_events * sizeof(double));
  stats->num_lineages = (int *) malloc(max_events * sizeof(int));
  stats->event_types = (int *) malloc(max_events * sizeof(int));
  stats->event_ids = (int *) malloc(max_events * sizeof(int));
}

void allocateMigBandsMem() {
  int maxMigBands = dataSetup.popTree->numMigBands;
  for (int comb = 0; comb < dataSetup.popTree->numPops; comb++) {
    if (isFeasibleComb(comb)) {
      comb_stats[comb].leafMigs = (MigStats *) malloc(maxMigBands * sizeof(MigStats));
    }
  }
}

int isFeasibleComb(int pop) {
  return !isLeaf(pop);
}

int isMigOfComb(int mig, int comb) { // if migrations flow into the comb
  int target = getTargetPop(mig);
  return isAncestralTo(comb, target) || target == comb;
}

int isMigBandExternal(int mig, int comb) {
  int source = getSourcePop(mig);
  int target = getTargetPop(mig);
  return isAncestralTo(comb, target) && !isAncestralTo(comb, source);
}

int isMigBandInternal(int mig, int comb) {
  int source = getSourcePop(mig);
  int target = getTargetPop(mig);

  return isAncestralTo(comb, source) && isAncestralTo(comb, target) && (!isLeaf(source) || !isLeaf(target));
}

int isCombLeafMigBand(int mig, int comb) {
  int target = getTargetPop(mig);
  int source = getSourcePop(mig);
  return isLeaf(target) && isLeaf(source) && isAncestralTo(comb, target) && isAncestralTo(comb, source);
}

double getCombAge(int comb) {
  if (SET_COMB_AGE_TO_ZERO) return 0.0;
  if (isLeaf(comb)) {
    return DBL_MAX;
  } else if (areChildrenLeaves(comb)) {
    return dataSetup.popTree->pops[comb]->age;
  } else if (isFeasibleComb(comb)) {
    double left_min = getCombAge(dataSetup.popTree->pops[comb]->sons[LEFT]->id);
    double right_min = getCombAge(dataSetup.popTree->pops[comb]->sons[RIGHT]->id);
    return fmin(left_min, right_min);
  } else {
    printf("ERROR: bug in combAge algorithm. Should not reach here!");
    exit(-1);
    return DBL_MAX;
  }
}

bool isEventCompletelyBelowComb(double eventAge, double combAge) {
  return (eventAge < combAge) && !areAlmostEqual(eventAge, combAge);
}

bool isBorderEvent(double eventAge, double previousAge, double combAge) {
  return areAlmostEqual(eventAge, combAge) ||
         ((eventAge > combAge) && (previousAge < combAge));
}

bool isEventCompletelyInsideComb(double eventAge, double combAge) {
  return (eventAge > combAge) && (!areAlmostEqual(eventAge, combAge));
}


void freeCombMem() { // TODO - implement
}
