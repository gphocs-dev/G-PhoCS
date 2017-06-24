#ifndef G_PHOCS_COMBASSERTIONS_H
#define G_PHOCS_COMBASSERTIONS_H


void runAssertions();

void assertCombLeavesEventChains();

void assertRootNumCoals();

void assertRootCoalStats();

void assertBottomCombs();

void assertBottomCombsNumCoals(int comb);

void assertBottomCombsCoalStats(int comb);

void assertCombLeaves();

void assertCombLeafNumCoals(int comb, int leaf);

void assertCombLeafCoalStats(int comb, int leaf);

void assertMigStats();

void assertLeafMigStats(int migband, int comb);

void assertLeafMigMigStats(int migband, int comb);

void assertLeafMigNumMigs(int migband, int comb);

void assertCombLeafEventChain(int comb, int leaf, int gene);

int getLastEvent(EventChain chain, int firstEvent);

void assertLastEventId(int lastEvent);

void assertLastEventIsSampleEnd(EventChain chain, int lastEvent);

void assertLastEventIsAtLeastAsOldAsComb(int comb, EventChain chain, int lastEvent);

void assertChainHasAtLeastTwoEvents(EventChain chain, int firstEvent);

void assertFirstEventZeroElapsedTime(EventChain chain, int firstEvent);

void assertFirstEventIsSampleStart(EventChain chain, int firstEvent);

void printErrorAndExit(string errorMessage);


#endif
