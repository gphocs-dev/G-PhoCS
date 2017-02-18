
#ifndef SRC_COMBSTATS_H_
#define SRC_COMBSTATS_H_

#define TRUE 1 // TODO - ask Ilan: where should these consts be?
#define FALSE 0


// --- GLOBAL DATA STRUCTURES -------------------------------------------------

typedef struct STATS {
			int num_coals, num_migs;
			double coal_stats, mig_stats;
			double* sorted_ages; 	// temporary array to hold sorted node ages of a specific genealogy.
			double* elapsed_times;	// temporary array to hold elapsed_time between adjacent events of a specific genealogy.
			int* num_lineages; 		// temporary array to hold num of lineages of a specific genealogy.
			int* event_types; 		// temporary array to hold event types
			/** size of stats arrays (event_types, num_lineages, elapsed_times, sorted_ages)*/
			int num_events;
} Stats;
typedef struct LEAF_STATS {
	Stats above_comb, below_comb;
  } LeafStats;
struct COMB_STATS{
	double age; 			// temporary value of comb-age
	Stats total;
	Stats* clades;
	LeafStats* leaves;
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
void appendCurrent(int comb, int currentPop, int gene);
double calculateCoalStats(double* elapsed_times, int* num_lineages, int size);


int	isLeaf(int pop);
int areChildrenLeaves(int pop);
int isFeasibleComb(int pop);
int isAncestralTo(int father, int son);

double getCombAge(int comb);
int getSon(int pop, int SON);
Stats* getStats(int comb, int pop);
void initCombStats();
void initStats(Stats* stats);

void allocateCombMem();
void allocateStats(Stats* stats);

void freeCombMem();


// TODO - extract tests to different source file
void assertRootNumCoals();
//void test_all(int comb);
//void test_validateCoalStats(int comb);
//void test_validateCountCoals(int comb);
//int test_countCoalsCombLeaves_rec(int comb, int pop);
//int test_countCoalsEntireComb(int comb);
//int test_countCoalsEntireClade(int comb);

#endif /* SRC_COMBSTATS_H_ */
