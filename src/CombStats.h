
#ifndef SRC_COMBSTATS_H_
#define SRC_COMBSTATS_H_


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
	double debug_total_error;

} *comb_stats;


// --- FUNCTION DECLARATIONS -----------------------------------------------


// clade_stats calculation functions
void computeCombStats();
void initCombStats();
void initSpecificCombStats(int clade);
void computeCombNumMigs();
void computeCombMigStats();
void computeCombNumCoals();
void computeCombNumCoals_rec(int pop);
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

int	isLeafPopulation(int pop);
int printLine();


#endif /* SRC_COMBSTATS_H_ */
