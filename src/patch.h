/*
 * patch.h
 *
 *  Created on: Feb 4, 2017
 *      Author: ron
 */

#ifndef SRC_PATCH_H_
#define SRC_PATCH_H_


#ifndef NULL
#define NULL   ((void *) 0)
#endif


#define MAX_MIG_BANDS	100			// max migration bands in the population tree (increased in V1.3.2 to 100)
#define MAX_MIGS		10			// max migration events per genealogy
#define NSPECIES		20			// max # of species
#define NS				200			// max # of sequences
#define OLDAGE			999			// upper bound on age (can be extended...)
#define MAX_EVENTS   (NS + 2*NSPECIES + 3*MAX_MIGS)
#define NUM_DELTA_STATS_INSTANCES 2

#define DEBUG_NODE_CHANGE_NOT
#define DEBUG_RUBBERBAND_NOT
/* END OF PATCH.C DEFS */


/***************************************************************************************************************/
/******                                  GLOBAL DATA STRUCTURES                                           ******/
/***************************************************************************************************************/




typedef struct {
	int accepted;
	double datastate_dataLogLikelihood;
	double datastate_logLikelihood;
	double finetune;
	double genLogLikelihood;
	int gen;
	int id; // some debug parameters for multi thread debugging
	int result_id; // some debug parameters for multi thread debugging
} UpdateGB_InternalNode_ReturnData;

typedef struct {
  double coal_stats[2*NSPECIES-1], mig_stats[MAX_MIG_BANDS];
  int num_coals[2*NSPECIES-1], num_migs[MAX_MIG_BANDS];
} GENETREE_STATS;

/*	genetree_stats_delta
	Holds difference in stats for populations affected by changes in a genetree.
	Not generally used for changes in theta, mig_rates.
	For easy and fast update if changes are accepted.
	changed_events holds the ids of all events in the interval between original and
	new position of event. All these intervals have a change in number of lineages.
*/
typedef struct {
  int original_event;							// id of event describing original placing of node
  int updated_event;							// id of event describing updated (new) placing of node
  int num_lin_delta;							// the difference in lineage number for all events affected (typically +1 or -1)
  int num_changed_events;						// number of events affected by change
  int changed_events[MAX_EVENTS];			// an array of ids of events affected by change
  int num_pops_changed;						// number of population affected by change
  int pops_changed [2*NSPECIES-1];			// an array of populations affected by change
  int num_mig_bands_changed;					// number of migration bands affected by change
  int mig_bands_changed[MAX_MIG_BANDS];	// an array of migration bands affected by change
  double coal_stats_delta[2*NSPECIES-1];	// difference in coalescence statistics per population affected
  double mig_stats_delta[MAX_MIG_BANDS];	// difference in migration statistics per migration band affected
} GENETREE_STATS_DELTA;

/* rubberband_migs
   structure for holding data on migration events out of rubber-banded populations
   which are affected by rubber-band operation (per gen).
   num_moved_events - total number of affected events. (in/out migrations and start/end of migration bands).
   orig_events      - original copies of events.
   new_events       - new copies of events.
   pops				  - population in which each event resides.
   new_ages			  - age of each new event.

   rubberband_migs is an array of size numLoci allocated in getMem().
*/

typedef struct {

  int num_moved_events, orig_events[MAX_MIGS+MAX_MIG_BANDS], new_events[MAX_MIGS+MAX_MIG_BANDS], pops[MAX_MIGS+MAX_MIG_BANDS];
  double new_ages[MAX_MIGS+MAX_MIG_BANDS];

} RUBBERBAND_MIGS;

/*	mig_spr_stats
	holds statistics for the SPR sampling operation with migration.
	In use in UpdateGB_MigSPR and in traceLineage.
*/
typedef struct {
  int father_event_old, father_event_new;
  int father_pop_new;
  int target;
  int num_old_migs, num_new_migs;
  int old_migs[MAX_MIGS], new_migs_in[MAX_MIGS], new_migs_out[MAX_MIGS], new_migs_bands[MAX_MIGS];
  double new_migs_ages[MAX_MIGS];
  double genetree_delta_lnLd[2];
} MIG_SPR_STATS;


typedef struct {
	GENETREE_STATS genetree_stats_check;
	GENETREE_STATS_DELTA genetree_stats_delta[NUM_DELTA_STATS_INSTANCES];
	MIG_SPR_STATS mig_spr_stats;
	double genDeltaLogLikelihood;
	double genLogLikelihood;
	RUBBERBAND_MIGS rubberband_migs;
	int mig_conflict_log;
} Locus_SuperStruct;
Locus_SuperStruct *locus_data;

//double averageMigTimes[MAX_MIG_BANDS]; //Unused

/* node surrogates.
 */
int**	nodePops;			// a 2D array (numLoci X numNodes) for populations per genealogy node.
int**	nodeEvents;			// a 2D array (numLoci X numNodes) for event ids per genealogy node.

struct ADMIXTURE_STATUS {
	int 	numSampledLoci;				// number of loci for which to show admixture status
	int*	sampledLoci;				// array of sampled loci
	double** sampleLocusAdmixRate;			// array of admixture rates (up to the curent point of sampling) for each locus and each admixed sample
	int*	admixtureCounts;			// number of loci in which sample is in alternative population
	double* admixtureCoefficients;		// estimated coefficients for all samples
} admixture_status;

/* genetree_migs is a struct which contains information about
   migration events in a specific genealogy.
   Array of structs (of length numLoci is allocated in GetMem().
*/
struct GENETREE_MIGS {
  int num_migs;							// number of migration event in genetree
  int living_mignodes[MAX_MIGS];	// array of indices for mignodes for dynamic managing of migration nodes
  struct MIGNODE{
    int gtree_branch;						// id of gene tree node below the migration node
    int migration_band;					// id of migration band relevant to this node
    int target_pop, source_pop;		// source and target populations of migration event (backwards view)
    int target_event, source_event;	// ids of event corresponding to this migration node.
    double age;								// time of event
  } mignodes[MAX_MIGS];
}* genetree_migs;


/* event chain
   Each event corresponds to a time band within a population where no events
   (coalescence/migration) take place. An event is attributed with one of 5 types
   corresponding to the event taking place at the end of the interval. Events are
   sorted in a list according to chronology within a population.
   Actual array of events is allocated in getMem()
*/

typedef enum event_type {COAL, IN_MIG, OUT_MIG, MIG_BAND_START, MIG_BAND_END, SAMPLES_START, END_CHAIN, DUMMY} EventType;
typedef struct EVENT{
	EventType type;
    int node_id, next, prev;
    double elapsed_time;				// time from last event
    int num_lineages;					// number of lineages before the event
} Event;
struct EVENT_CHAIN{
  int total_events;						// total number of events pre-allocated to this chain
  int first_event[2*NSPECIES-1];		// pointers to first event for every population
  int last_event[2*NSPECIES-1];			// pointers to last event for every population
  int free_events;						// pointer to a chain of free events for use. Always have at least one free event
  Event* events;
} *event_chains;


/*	genetree stats
	holds relevant statistics for fast computation of probability of tree
	given model parameters (split times, pop sizes and migration rates).
	Holds statistics for every population is species tree.
	-- NEED TO ADD DOCUMENTATION FOR THIS --  !!!
	genetree_stats[gen] holds relevant information for genetree of genealogy 'gen'.
    array allocated in GetMem().
	genetree_stats_total holds sum of statistics for all loci. coal_stats here
    considers also all gen-specific heredity factors (but not thetas).
*/

GENETREE_STATS *genetree_stats, *genetree_stats_total_partitioned, genetree_stats_total , genetree_stats_total_check;

/*	genetree stats flat
	holds relevant statistics for fast computation of probability of tree
	given a null model (single population).
    array allocated in GetMem().
*/
struct GENETREE_STATS_FLAT{
  double coal_stats_flat, mig_stats_flat;
  int num_coals_total, num_migs_total;
  double* sortedAgesArray;
} genetree_stats_flat;

/*	genetree node stata
	holds relevant statistics for coalescent distribution across populations
    array allocated in GetMem().
*/
struct GENETREE_NODE_STATS{
  /*** memory allocation ***/
  double*   doubleArray;
  double**  doublePtrArray;
  double*** doublePtrPtrArray;
  int*      intArray;
  int**     intPtrArray;

  /*** auxilliary matrices for computation ***/
  int**     lcaMatrix;            // matrix of LCAs for all pairs of leaves
  int*      leafArray;            // array of leaves for LCA computation
  int*      firstNodesInPop;      // array with the first coal node in each pop
  double*   firstNodesAges;       // array with the age of the first coal node in each pop
  double*   nodeAges;              // array with the ages of all internal nodes

  /*** main stat matrices ***/
  double*** probCoalMatrix;       // Prob[sample pair coalesce in pop]
  double*** probFirstCoalMatrix;  // Prob[sample pair is the first coalescence in pop]
  double*** coalTimeMatrix;       // mean coal time for coalescence of pair (given they coalesce in pop)
} genetree_node_stats;









/******************************************************************************************************/
/******                                FUNCTION DECLARATIONS                                     ******/
/******************************************************************************************************/



int analyzeGenetreeFile(char* genetree_file);
int writeMScommandLine(char* outfile);
int findLastMig(int gen, int node_id, double time);
int findFirstMig(int gen, int node_id, double time);
int Coalescence1Pop (PopulationTree* popTree, GenericBinaryTree* tree, int gen, int pop, int* livingLineages);
int removeEvent(int gen, int event);

// auxiliary functions
int findInconsistency(int gen, int node);
int getSptreeNodeByName(const char* name);
int orderByAge(int subtree_root, int* ordered_nodes);
int getLineagesAtInterval(int gen, int start_event, int pop, int exc_node, int* out_array);
int getEdgesForTimePop(int gen, double time, int pop, int exc_node, int* out_array);
int populationPostOrder(int pop, int* ordered_pops);

// event chain functions
int checkAll();
int checkGtreeStructure(int gen);
int synchronizeEvents(int gen);
int printEventChains(FILE* stream, int gen);
int createEvent(int gen, int pop, double age);
int createEventBefore(int gen, int pop, int event, double elapsed_time);
int constructEventChain(int gen);
int computeFlatStats();
int computeNodeStats();
int computeTotalStats();
double recalcStats(int gen, int pop);
int recalcStats_partitioned(int gen, int pop);
int computeGenetreeStats(int gen);
int computeGenetreeStats_partitioned(void);
double gtreeLnLikelihood(int gen);
double computeDeltaLnLd(int gen, int instance);
int computeMigStatsDelta(int instance, double bottom_age, int bottom_pop, double top_age, int num_lins_delta , int gen);
int computeCoalStatsDelta(int instance, int gen, int bottom_event, int bottom_pop, int top_event, int num_lins_delta);
double considerEventMove(int gen, int instance, int event_id, int source_pop, double original_age, int target_pop, double new_age);
int acceptEventChainChanges(int gen, int instance);
int rejectEventChainChanges(int gen, int instance);
double rubberBandRipple(int gen, int do_or_redo);
double rubberBand(int gen, int pop, double static_point, double moving_point, double factor, unsigned short postORpre, int* out_num_events);
int replaceMigNodes(int gen, int node);
int traceLineage(int gen, int node, int reconnect);









#endif /* SRC_PATCH_H_ */
