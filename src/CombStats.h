
#ifndef SRC_COMBSTATS_H_
#define SRC_COMBSTATS_H_

#define TRUE 1 // TODO - ask Ilan: where should these consts be?
#define FALSE 0


// --- GLOBAL DATA STRUCTURES -------------------------------------------------

struct COMB_STATS{
	double coal_stats_total, mig_stats_total;
	int num_coals_total, num_migs_total;

	double* sorted_ages; 	// temporary array to hold sorted node ages of a specific genealogy.
	double* elapsed_times;	// temporary array to hold elapsed_time between adjacent events of a specific genealogy.
	int* num_lineages; 		// temporary array to hold num of lineages of a specific genealogy.
	int* event_types; 		// temporary array to hold event types

	int num_events; 		// temporary number of events in clade.
	// used as an "array length" variable for sortedAges, event_types & num_lineages arrays
	double debug_total_error; // TODO - debugggggggggg
	double age; 			// temporary value of leaves-trim-age //TODO - rewrite this comment
	struct COMB_LEAF {
	    int num_coals_total, num_migs_total;
	    double coal_stats, mig_stats;
	    int debug_original_num_coals;
	  }* leaves;

} *comb_stats;


// --- FUNCTION DECLARATIONS -----------------------------------------------


// clade_stats calculation functions
void calculateCombStats();
void initCombStats();
void initStats();
void computeCombNumMigs();
void computeCombMigStats();
void countCoals();
void countCoals_rec(int comb, int pop);
void computeCombCoalStats();
void computeCombCoalStats_rec(int clade, int gen);
void fillupLeafCombStats(int clade, int gen);
void appendPopToComb(int clade, int gen, int startingPoint);
void fillupCombStats(int clade, int gen);
void addChildenIntoCombStats(int clade, int gen);
void mergeChildern(int clade, int gen);
void addChildrenCombStats(int clade, int gen);
void addCurrentPopIntoCombStats(int clade, int gen);
double getCoalStats(double* elapsed_times, int* num_lineages, int size);



void countLeafCoals(int comb, int leaf);
void countLeafGeneCoals(int comb, int leaf, int gene);
void countNonLeafCoals(int comb, int pop); // TODO - tidy declarations
int getChildCladeNumCoal(int child);
void computeCombLeafNumCoals(int leafPop);
void countCoalsOutsideTrimmedLeafGene(int leafPop, int gene);
void freeCombMem();
void allocateCombMem();

int	isLeaf(int pop);
int areChildrenLeaves(int pop);
int isFeasableComb(int pop);

double getCombAge(int comb);

#endif /* SRC_COMBSTATS_H_ */
