#include "CombStats.h"
#include "CombAssertions.h"
#include "utils.h"
#include "MCMCcontrol.h"
#include "patch.h"
#include "MemoryMng.h"
#include "McRefCommon.h"

void runAssertions() {
  if (SET_COMB_AGE_TO_ZERO) {
    computeFlatStats();
    assertRootNumCoals();
    assertRootCoalStats();
  } else {
    assertBottomCombs();
    assertCombLeaves();
    assertMigStats();
    assertCombLeavesEventChains();
  }
}

void assertRootNumCoals() { // this test only works when combAge is set to zero
  int root = getPopIdByName(dataSetup.popTree, "root");

  int maxCoals = dataSetup.numLoci * (dataSetup.numSamples - 1);

  int actualCoals = comb_stats[root].total.num_coals;
  for (int pop = 0; pop < dataSetup.popTree->numPops; pop++) {
    if (isLeaf(pop)) {
      actualCoals += comb_stats[root].leaves[pop].below_comb.num_coals;
    }
  }

  if (actualCoals != maxCoals) {
    printf("comb '%s' #%d: Expected coalescence events - %d. Actual coalescence events - %d",
           "root", root, maxCoals, actualCoals);
    exit(-1);
  }
}

/**
 * Tests whether root comb with mocked comb_age:=0.0 gives the same coal_stats as a flat model.
 *  To enable this test, you must hard-code getCombAge() to return 0.0 AND enable flat stats
 */
void assertRootCoalStats() {
  int root = getPopIdByName(dataSetup.popTree, "root");
  double actualCoalStats = comb_stats[root].total.coal_stats;
  double expectedCoalStats = genetree_stats_flat.coal_stats_flat;

  double error = fabs(actualCoalStats - expectedCoalStats);
  double relativeError = error / expectedCoalStats;

  if (relativeError > requiredRelativePrecision()) {
    printf("Error while checking 'root' #%d coal_stats:\nExpected:%0.35f\tActual:%0.35f\tRelative Error:%0.35f",
           root, expectedCoalStats, actualCoalStats, relativeError);
    exit(-1);
  }
}

/**
 * Compares trivial combs (pops directly above two leaves) with regular pops.
 */
void assertBottomCombs() {
  for (int comb = 0; comb < dataSetup.popTree->numPops; comb++) {
    if (isFeasibleComb(comb) && areChildrenLeaves(comb)) {
      assertBottomCombsNumCoals(comb);
      assertBottomCombsCoalStats(comb);
    }
  }
}

void assertBottomCombsNumCoals(int comb) {
  int expected = genetree_stats_total.num_coals[comb];
  int actual = comb_stats[comb].total.num_coals;

  if (expected != actual) {
    printf("\nError while checking comb '%s' #%d num_coals:\nExpected num_coals %d. actual is %d",
           getPopName(comb), comb, expected, actual);
    exit(-1);
  }
}

void assertBottomCombsCoalStats(int comb) {
  double expected = genetree_stats_total.coal_stats[comb];
  double actual = comb_stats[comb].total.coal_stats;
  double error = fabs(actual - expected);
  double relativeError = error / expected;
  if (relativeError > requiredRelativePrecision()) {
    printf("\nError while checking comb '%s' #%d coal_stats:\nExpected:%0.35f\tActual:%0.35f\tRelative Error:%0.35f\tAbsolute Error:%0.35f",
           getPopName(comb), comb,
           expected, actual, relativeError, error);
    exit(-1);
  }
}

void assertCombLeaves() {
  for (int comb = 0; comb < dataSetup.popTree->numPops; comb++) {
    if (isFeasibleComb(comb)) {
      for (int leaf = 0; leaf < dataSetup.popTree->numPops; leaf++) {
        if (isLeaf(leaf) && isAncestralTo(comb, leaf)) {
          assertCombLeafNumCoals(comb, leaf);
          assertCombLeafCoalStats(comb, leaf);
        }
      }
    }
  }
}

void assertCombLeafNumCoals(int comb, int leaf) {
  int expectedNumCoals = genetree_stats_total.num_coals[leaf];
  int actualNumCoals = comb_stats[comb].leaves[leaf].above_comb.num_coals
                       + comb_stats[comb].leaves[leaf].below_comb.num_coals;
  if (expectedNumCoals != actualNumCoals) {
    printf("\nError while checking leaf '%s' #%d num_coals:\nExpected num_coals %d. actual is %d",
           getPopName(leaf), leaf, expectedNumCoals, actualNumCoals);
    exit(-1);
  }
}

void assertCombLeafCoalStats(int comb, int leaf) {
  double expectedCoalStats = genetree_stats_total.coal_stats[leaf];
  double actualCoalStats = comb_stats[comb].leaves[leaf].above_comb.coal_stats
                           + comb_stats[comb].leaves[leaf].below_comb.coal_stats;
  double error = fabs(actualCoalStats - expectedCoalStats);
  double relativeError = error / expectedCoalStats;
  if (relativeError > requiredRelativePrecision()) {
    printf("\nError while checking leaf '%s' #%d coal_stats:\n"
                  "Expected:%0.35f\n"
                  "  Actual:%0.35f\n"
                  "   Above:%0.35f\n"
                  "   Below:%0.35f\n"
                  "RelError:%0.35f\n"
                  "AbsError:%0.35f\n",
           getPopName(leaf), leaf,
           expectedCoalStats,
           actualCoalStats,
           comb_stats[comb].leaves[leaf].above_comb.coal_stats,
           comb_stats[comb].leaves[leaf].below_comb.coal_stats,
           relativeError,
           error);
    exit(-1);
  }
}

void assertMigStats() {
  for (int migband = 0; migband < dataSetup.popTree->numMigBands; migband++) {
    for (int comb = 0; comb < dataSetup.popTree->numPops; comb++) {
      if (isFeasibleComb(comb) && isCombLeafMigBand(migband, comb)) {
        assertLeafMigStats(migband, comb);
      }
    }
  }
}

void assertLeafMigStats(int migband, int comb) {
  assertLeafMigMigStats(migband, comb);
  assertLeafMigNumMigs(migband, comb);
}

void assertLeafMigMigStats(int migband, int comb) {
  double expectedMigStats = genetree_stats_total.mig_stats[migband];
  double below_combage_migstats = comb_stats[comb].leafMigs[migband].mig_stats;
  double above_combage_migstats = comb_stats[comb].leafMigs[migband].mig_stats_above;
  double actualMigStats = below_combage_migstats + above_combage_migstats;
  double error = fabs(actualMigStats - expectedMigStats);
  double relativeError = error / expectedMigStats;
  if (relativeError > requiredRelativePrecision()) {
    int source = dataSetup.popTree->migBands[migband].sourcePop;
    int target = dataSetup.popTree->migBands[migband].targetPop;
    printf("\nError while checking migband '%s->%s'  '#%d->#%d' mig_stats for comb '%s':\n"
                  "Expected:%0.35f\n"
                  "  Actual:%0.35f\n"
                  "   Below:%0.35f\n"
                  "   Above:%0.35f\n"
                  "RelError:%0.35f\n"
                  "AbsError:%0.35f\n",
           getPopName(source), getPopName(target),
           source, target,
           getPopName(comb),
           expectedMigStats,
           actualMigStats,
           below_combage_migstats,
           above_combage_migstats,
           relativeError,
           error);
    exit(-1);
  }
}

void assertLeafMigNumMigs(int migband, int comb) {
  int expectedNumMig = genetree_stats_total.num_migs[migband];
  int actualNumMig = comb_stats[comb].leafMigs[migband].num_migs + comb_stats[comb].leafMigs[migband].num_migs_above;
  if (expectedNumMig != actualNumMig) {
    int source = dataSetup.popTree->migBands[migband].sourcePop;
    int target = dataSetup.popTree->migBands[migband].targetPop;
    printf("\nError while checking migband '%s->%s' '#%d->#%d' num migs for comb '%s'. Expected:%d. Actual:%d",
           getPopName(source), getPopName(target),
           source, target,
           getPopName(comb),
           expectedNumMig, actualNumMig);
    exit(-1);
  }
}

void assertCombLeavesEventChains() {
  for (int comb = 0; comb < dataSetup.popTree->numPops; comb++) {
    for (int leaf = 0; leaf < dataSetup.popTree->numPops; leaf++) {
      if (isFeasibleComb(comb) && isLeaf(leaf) && isAncestralTo(comb, leaf)) {
        for (int gene = 0; gene < dataSetup.numLoci; gene++) {
          assertCombLeafEventChain(comb, leaf, gene);
        }
      }
    }
  }
}

void assertCombLeafEventChain(int comb, int leaf, int gene) {

  EventChain chain = event_chains[gene];
  int firstEvent = chain.first_event[leaf];

  int lastEvent = getLastEvent(chain, firstEvent);

  assertLastEventId(lastEvent);
  assertLastEventIsSampleEnd(chain, lastEvent);
  assertLastEventIsAtLeastAsOldAsComb(comb, chain, firstEvent);

  assertChainHasAtLeastTwoEvents(chain, firstEvent);
  assertFirstEventZeroElapsedTime(chain, firstEvent);
  assertFirstEventIsSampleStart(chain, firstEvent);
}

int getLastEvent(EventChain chain, int firstEvent) {
  int lastEvent = firstEvent;
  while (hasNextEvent(chain, lastEvent)) {
    lastEvent = chain.events[lastEvent].getNextIdx();
  }
  return lastEvent;
}

void assertLastEventId(int lastEvent) {
  if (lastEvent < 0) printErrorAndExit("Last event is negative");
}

void assertLastEventIsSampleEnd(EventChain chain, int lastEvent) {
  EventType type = chain.events[lastEvent].getType();
  if (type != END_CHAIN) printErrorAndExit("Last event isn't END_CHAIN");
}

void assertLastEventIsAtLeastAsOldAsComb(int comb, EventChain chain, int firstEvent) {
  double combAge = getCombAge(comb);
  double lastEventAge = 0.0;
  int event = firstEvent;
  for (; hasNextEvent(chain, event); event = chain.events[event].getNextIdx()) {
    lastEventAge += chain.events[event].getElapsedTime();
  }
  lastEventAge += chain.events[event].getElapsedTime();

  double absError = combAge - lastEventAge;
  double relativeError = absError * 2 / (lastEventAge + combAge);

  if (lastEventAge - combAge < (-requiredRelativePrecision())) {
    fprintf(stderr, "\nlastEventAge:%0.45f\n", lastEventAge);
    fprintf(stderr, "     combAge:%0.45f\n", combAge);
    fprintf(stderr, "absolute err:%0.45f\n", absError);
    fprintf(stderr, "relative err:%0.45f\n", relativeError);
    printErrorAndExit("last event is younger than comb");
  }
}

void assertChainHasAtLeastTwoEvents(EventChain chain, int firstEvent) {
  int i = 1;
  int event = firstEvent;
  while (hasNextEvent(chain, event)) {
    event = chain.events[event].getNextIdx();
    i++;
  }
  if (i < 2) printErrorAndExit("Event chain has less than 2 events");
}

void assertFirstEventZeroElapsedTime(EventChain chain, int firstEvent) {
  double firstElapsedTime = chain.events[firstEvent].getElapsedTime();
  if (firstElapsedTime != 0.0) printErrorAndExit("first elapsed time isn't zero");
}

void assertFirstEventIsSampleStart(EventChain chain, int firstEvent) {
  EventType type = chain.events[firstEvent].getType();
  if (type != SAMPLES_START) printErrorAndExit("first event isn't SAMPLE_START");
}

void printErrorAndExit(string errorMessage) {
  fprintf(stderr, "%s", errorMessage.c_str());
  exit(-1);
}
