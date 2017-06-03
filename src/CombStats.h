
#ifndef SRC_COMBSTATS_H_
#define SRC_COMBSTATS_H_

#define TRUE 1 // TODO - ask Ilan: where should these consts be?
#define FALSE 0


// --- GLOBAL DATA STRUCTURES -------------------------------------------------

typedef struct STATS {
			int num_coals;
			double coal_stats;
			double* sorted_ages; 	// temporary array to hold sorted node ages of a specific genealogy.
			double* elapsed_times;	// temporary array to hold elapsed_time between adjacent events of a specific genealogy.
			int* num_lineages; 		// temporary array to hold num of lineages of a specific genealogy.
			int* event_types; 		// temporary array to hold event types
			int* event_ids;
			/** size of stats arrays (event_types, num_lineages, elapsed_times, sorted_ages)*/
			int num_events;
} Stats;
typedef struct MIGSTATS {
			int num_migs;
			double mig_stats;
} MigStats;
typedef struct LEAF_STATS {
	Stats above_comb, below_comb;
  } LeafStats;

typedef struct _COMB_STATS{
	double age; 			// temporary value of comb-age
	Stats total;
	Stats* combs;
	LeafStats* leaves;
	MigStats* leafMigs;
	//MigStats* combMigs;
} COMB_STATS;

extern COMB_STATS* comb_stats;


// --- FUNCTION DECLARATIONS -----------------------------------------------


// clade_stats calculation functions
void calculateCombStats();
void finalizeCombCoalStats(int comb);
void calculateSufficientStats(int comb, int gene);
void coalescence(int comb, int gene);
void coalescence_rec(int comb, int currentPop, int gene);
void handleLeafCoals(int comb, int leaf, int gene);
void handleBelowCombAge(int gene, double combAge, int* eventId, double* elapsedTime, int* eventType, int* numLineages, double* eventAge, double* previousAge, Stats* belowCombLeafStats);
void handleAboveCombAge(int gene, double combAge, int* eventId, double* elapsedTime, int* eventType, int* numLineages, double* eventAge, double* previousAge, Stats* belowCombLeafStats, Stats* aboveCombLeafStats, Stats* combTotalStats);
void handleFirstEventAboveCombAge(int gene, double combAge, int* eventId, double*elapsedTime, int* eventType, int* numLineages, double* eventAge, double* previousAge, Stats* belowCombLeafStats, Stats* aboveCombLeafStats);
void handleRestOfEventsAboveCombAge(int gene, int* eventId, double*elapsedTime, int* eventType, int* numLineages, double* eventAge, double* previousAge, Stats* aboveCombLeafStats, Stats* combTotalStats);


void handleNonLeafCoals(int comb, int currentPop, int gene);
void handleNonLeafNumCoals(int comb, int currentPop, int gene);
void handleNonLeafCoalStats(int comb, int currentPop, int gene);
void mergeChildernIntoCurrent(int comb, int currentPop, int gen);
void copyStaticEventStats(Stats* sourceStats, int n, Stats* targetStats, int m);
void appendCurrent(int comb, int currentPop, int gene);

void migrations(int comb, int gene);
void handleLeafMigStats(int comb, int mig, int gene);
void updateLeafMigStats(int numLineages, double elapsedTime, int eventType, MigStats* migLeafStats);
void fastFwdPastMigBandStart(int gene, int* eventId, double*elapsedTime, int* eventType, int* numLineages, double* eventAge, double* previousAge);
void setupFirstEventVars(int gene, int currentPop, int* eventId, double*elapsedTime, int* eventType, int* numLineages);
void incrementEventVars(int gene, int *eventId, double*elapsedTime, int* eventType, int* numLineages, double* eventAge, double* previousAge);
void handleExternalMigStats(int comb, int mig, int gene); // TODO - rename to better explain "O->C" type migband


Stats* getCombPopStats(int comb, int pop);
void initCombStats();
void initPopStats();
void initMigStats();
void initStats(Stats* stats);

void allocateCombMem();
void allocatePopsMem();
void allocateMigBandsMem();
void allocateStats(Stats* stats);

void freeCombMem();
double getCombAge(int comb);
int isFeasibleComb(int pop);
int isMigOfComb(int mig, int comb);
int isMigBandExternal(int mig, int comb);
int isMigBandInternal(int mig, int comb);
int isLeafMigBand(int mig, int comb);
bool hasNextEvent(EventChain chain, int event); // TODO - extract method to McRefCommon


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

void assertLeafEventChain(int comb, int leaf, int gene);
void assertLastEventId(int lastEvent);
void assertLastEventIsSampleEnd(EventChain chain, int lastEvent);
void assertLastEventIsAtleastAsOldAsComb(int comb, EventChain chain, int lastEvent);
void assertChainHasAtleastTwoEvents(EventChain chain, int firstEvent);
void assertFirstEventIsSampleStart(EventChain chain, int firstEvent);
void printErrorAndExit(char* errorMessage);

#endif /* SRC_COMBSTATS_H_ */
