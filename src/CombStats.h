
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
struct COMB_STATS{
	double age; 			// temporary value of comb-age
	Stats total;
	Stats* clades;
	LeafStats* leaves;
	MigStats* leafMigs;
	//MigStats* combMigs; //TODO - uncomment to see horrible bug
} *comb_stats;


// --- FUNCTION DECLARATIONS -----------------------------------------------


// clade_stats calculation functions
void calculateCombStats();
void finalizeCombCoalStats(int comb);
void calculateSufficientStats(int comb, int gene);
void coalescence(int comb, int gene);
void coalescence_rec(int comb, int currentPop, int gene);
void handleLeafCoals(int comb, int leaf, int gene);
void handleNonLeafCoals(int comb, int currentPop, int gene);
void handleNonLeafNumCoals(int comb, int currentPop, int gene);
void handleNonLeafCoalStats(int comb, int currentPop, int gene);
void mergeChildernIntoCurrent(int comb, int currentPop, int gen);
void copyStaticEventStats(Stats* sourceStats, int n, Stats* targetStats, int m);
void appendCurrent(int comb, int currentPop, int gene);
double calculateCoalStats(double* elapsed_times, int* num_lineages, int size);

void migrations(int comb, int gene);
void handleLeafMigStats(int comb, int mig, int gene);
void fastFwdPastMigBandStart(int gene, int* eventId, double*elapsedTime, int* eventType, int* numLineages, double* eventAge);
void incrementEventVars(int gene, int* eventId, double*elapsedTime, int* eventType, int* numLineages, double* eventAge);
void handleExternalMigStats(int comb, int mig, int gene); // TODO - rename to better explain "O->C" type migband


int	isLeaf(int pop);
int areChildrenLeaves(int pop);
int isFeasibleComb(int pop);
int isAncestralTo(int father, int son);
int isMigOfComb(int mig, int comb);
int isMigBandExternal(int mig, int comb);
int isMigBandInternal(int mig, int comb);
int isLeafMigBand(int mig, int comb);

char* getEventTypeName(int eventType);

double getCombAge(int comb);
int getSon(int pop, int SON);
int getSourcePop(int mig);
int getTargetPop(int mig);
char* getPopName(int pop);
char* getMigName(int mig);
char* concat(const char *s1, const char *s2); // TODO - THIS SHOULD NOT BE USED  IN PRODUCTION SINCE THE MEMORY ISN'T RELEASED

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

#endif /* SRC_COMBSTATS_H_ */
