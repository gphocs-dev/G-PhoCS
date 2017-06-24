#ifndef SRC_CLADESTATS_H_
#define SRC_CLADESTATS_H_

/*	clade_stats - holds statistics of all clades of tree.
 *  clade_stats[pop] holds statistics of the clade under (and including) pop.
 */
typedef struct _CLADE_STATS{
	double coal_stats_total, mig_stats_total;
	int num_coals_total, num_migs_total;

	double* sorted_ages; 	// temporary array to hold sorted node ages of a specific genealogy.
	double* elapsed_times;	// temporary array to hold elapsed_time between adjacent events of a specific genealogy.
	int* num_lineages; 		// temporary array to hold num of lineages of a specific genealogy.
	int* event_types; 		// temporary array to hold event types

	int num_events; 		// temporary number of events in clade.
	// used as an "array length" variable for sortedAges, event_types & num_lineages arrays

} CLADE_STATS;
extern CLADE_STATS* clade_stats;

void calculateCladeStats();

void allocateCladeMem();
void initCladeStats();
void initSpecificCladeStats(int clade);
void computeCladeNumMigs();
void computeCladeMigStats();
void computeCladeNumCoals();
void computeCladeNumCoals_rec(int pop);
void computeCladeCoalStats();
void computeCladeCoalStats_rec(int clade, int gen);
void fillUpLeafCladeStats(int clade, int gen);
void appendPopToClade(int clade, int gen, int startingPoint);
void fillUpCladeStats(int clade, int gen);
void addChildrenIntoCladeStats(int clade);
void mergeChildren(int clade);
void addChildrenCladeStats(int clade);
void addCurrentPopIntoCladeStats(int clade, int gen);


void freeCladeMem();

#endif
