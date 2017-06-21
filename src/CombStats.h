
#ifndef SRC_COMBSTATS_H_
#define SRC_COMBSTATS_H_

#include "DataLayer.h"

#define TRUE 1 // TODO - ask Ilan: where should these consts be?
#define FALSE 0


// --- GLOBAL DATA STRUCTURES -------------------------------------------------

typedef struct STATS {
    int num_coals;
    double coal_stats;
    double *sorted_ages;    // temporary array to hold sorted node ages of a specific genealogy.
    double *elapsed_times;    // temporary array to hold elapsed_time between adjacent events of a specific genealogy.
    int *num_lineages;        // temporary array to hold num of lineages of a specific genealogy.
    int *event_types;        // temporary array to hold event types
    int *event_ids;
    /** size of stats arrays (event_types, num_lineages, elapsed_times, sorted_ages)*/
    int num_events;
} Stats;
typedef struct MIGSTATS {
    int num_migs;
    double mig_stats;
    int num_migs_above;
    double mig_stats_above;
} MigStats;
typedef struct LEAF_STATS {
    Stats above_comb, below_comb;
} LeafStats;

typedef struct _COMB_STATS {
    double age;            // temporary value of comb-age
    Stats total;
    Stats *combs;
    LeafStats *leaves;
    MigStats *leafMigs;
} COMB_STATS;

extern COMB_STATS *comb_stats;


// --- FUNCTION DECLARATIONS -----------------------------------------------


// clade_stats calculation functions
void calculateCombStats();

void finalizeCombCoalStats(int comb);

void calculateSufficientStats(int comb, int gene);

void coalescence(int comb, int gene);

void coalescence_rec(int comb, int currentPop, int gene);

void handleLeafCoals(int comb, int leaf, int gene);

bool isEventCompletelyBelowComb(double eventAge, double combAge);

bool isBorderEvent(double eventAge, double previousAge, double combAge);

bool isEventCompletelyInsideComb(double eventAge, double combAge);

void countCoalEventTowardsBelowComb(Event event, Stats *belowCombLeafStats);

void countCoalEventTowardsHalfAndHalf(Event event, double eventAge, double previousAge, double combAge, Stats *belowCombLeafStats, Stats *aboveCombLeafStats);

void countCoalEventTowardsAboveComb(Event event, double eventAge, Stats *aboveCombLeafStats, Stats *combTotalStats);


void handleNonLeafCoals(int comb, int currentPop, int gene);

void handleNonLeafNumCoals(int comb, int currentPop, int gene);

void handleNonLeafCoalStats(int comb, int currentPop, int gene);

void mergeChildrenIntoCurrent(int comb, int currentPop);

void copyStaticEventStats(Stats *sourceStats, int n, Stats *targetStats, int m);

void appendCurrent(int comb, int currentPop, int gene);

void migrations(int comb, int gene);

void handleCombLeavesMigBand(int comb, int mig, int gene);

void countMigEventTowardsBelowComb(Event event, MigStats *leafMigStats);

void countMigEventTowardsHalfAndHalf(Event event, double eventAge, double previousAge, double combAge, MigStats *leafMigStats);

void assertLeafMigStats(int migband, int comb);

void assertLeafMigMigStats(int migband, int comb);

void assertLeafMigNumMigs(int migband, int comb);

void countMigEventTowardsAboveComb(Event event, MigStats *leafMigStats);

void handleExternalMigStats(int comb, int mig, int gene); // TODO - rename to better explain "O->C" type migband


Stats *getCombPopStats(int comb, int pop);

void initCombStats();

void initPopStats();

void initMigStats();

void initStats(Stats *stats);

void allocateCombMem();

void allocatePopsMem();

void allocateMigBandsMem();

void allocateStats(Stats *stats);

void freeCombMem();

double getCombAge(int comb);

int isFeasibleComb(int pop);

int isMigOfComb(int mig, int comb);

int isMigBandExternal(int mig, int comb);

int isMigBandInternal(int mig, int comb);

int isCombLeafMigBand(int mig, int comb);

bool hasNextEvent(EventChain chain, int event); // TODO - extract method to McRefCommon
bool areAlmostEqual(double eventAge, double combAge);

double relativeDistance(double dbl1, double dbl2); // TODO - extract method to McRefCommon


// TODO - extract tests to different source file
void debug_printCombGene(int comb);

void assertRootNumCoals();

void assertRootCoalStats();

void assertBottomCombs();

void assertBottomCombsCoalStats(int comb);

void assertBottomCombsNumCoals(int comb);

void assertCombLeaves();

void assertCombLeafNumCoals(int comb, int leaf);

void assertCombLeafCoalStats(int comb, int leaf);

void assertMigStats();

void assertLeafEventChain(int comb, int leaf, int gene);

int getLastEvent(EventChain chain, int firstEvent);

void assertLastEventId(int lastEvent);

void assertLastEventIsSampleEnd(EventChain chain, int lastEvent);

void assertLastEventIsAtleastAsOldAsComb(int comb, EventChain chain, int lastEvent);

void assertChainHasAtleastTwoEvents(EventChain chain, int firstEvent);

void assertFirstEventZeroElapsedTime(EventChain chain, int firstEvent);

void assertFirstEventIsSampleStart(EventChain chain, int firstEvent);

void printErrorAndExit(char *errorMessage);

#endif /* SRC_COMBSTATS_H_ */
