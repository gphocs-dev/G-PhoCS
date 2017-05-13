#pragma once
/*============================================================================
 File: DataLayer.h

 Core data structures and functions

 Version history:
 04-Feb-2017  ronw         Created.
 12-Apr-2017  evgenyidc    Going C++.
 03-May-2017  evgenyidc    Extracting global variables definitions to CPP.
 ============================================================================*/

#include "DataLayerConstants.h"
#include "GenericTree.h"
#include "PopulationTree.h"

#include "MCMCcontrol.h"
extern DATA_SETUP dataSetup;
extern int debug;

/*=============================================================================
 *
 * GLOBAL DATA STRUCTURES
 *
 *===========================================================================*/

/*-----------------------------------------------------------------------------
 * UpdateGB_InternalNode_ReturnData
 *---------------------------------------------------------------------------*/
typedef struct _UpdateGB_InternalNode_ReturnData
{
  int    accepted;
  double datastate_dataLogLikelihood;
  double datastate_logLikelihood;
  double finetune;
  double genLogLikelihood;
  int    gen;

  // some debug parameters for multi-thread debugging
  int    id;
  int    result_id;
} UpdateGB_InternalNode_ReturnData;

/*-----------------------------------------------------------------------------
 * GENETREE_STATS
 * Holds relevant statistics for fast computation of probability of tree
 * given model parameters (split times, pop sizes and migration rates).
 * Holds statistics for every population is species tree.
 * genetree_stats[gen] holds relevant information for genetree of
 * genealogy 'gen'. Array allocated in GetMem().
 * genetree_stats_total holds sum of statistics for all loci. coal_stats here
 * considers also all gen-specific heredity factors (but not thetas).
 *
 @@TODO: NEED TO ADD DOCUMENTATION FOR THIS --  !!!
 *---------------------------------------------------------------------------*/
typedef struct _GENETREE_STATS
{
  double coal_stats[2 * NSPECIES - 1];
  double mig_stats[MAX_MIG_BANDS];
  int    num_coals[2 * NSPECIES - 1];
  int    num_migs[MAX_MIG_BANDS];
} GENETREE_STATS;

/*-----------------------------------------------------------------------------
 * GENETREE_STATS_DELTA
 * Holds difference in stats for populations affected by changes
 * in a genetree. Not generally used for changes in theta, mig_rates.
 * For easy and fast update if changes are accepted.
 * changed_events holds the ids of all events in the interval between
 * original and new position of event. All these intervals have
 * a change in number of lineages.
 *---------------------------------------------------------------------------*/
typedef struct _GENETREE_STATS_DELTA
{
  // id of event describing original placing of node
  int original_event;
  // id of event describing updated (new) placing of node
  int updated_event;
  // the difference in lineage number for all events affected (typically +/- 1)
  int num_lin_delta;
  // number of events affected by change
  int num_changed_events;
  // an array of ids of events affected by change
  int changed_events[MAX_EVENTS];
  // number of population affected by change
  int num_pops_changed;
  // an array of populations affected by change
  int pops_changed[2 * NSPECIES - 1];
  // number of migration bands affected by change
  int num_mig_bands_changed;
  // an array of migration bands affected by change
  int mig_bands_changed[MAX_MIG_BANDS];
  // difference in coalescence statistics per population affected
  double coal_stats_delta[2 * NSPECIES - 1];
  // difference in migration statistics per migration band affected
  double mig_stats_delta[MAX_MIG_BANDS];
} GENETREE_STATS_DELTA;

/*-----------------------------------------------------------------------------
 * RUBBERBAND_MIGS
 * Structure for holding data on migration events out
 * of rubber-banded populations which are affected by
 * rubber-band operation (per gen).
 * num_moved_events - total number of affected events.
 *                    (in/out migrations and start/end of migration bands).
 * orig_events      - original copies of events.
 * new_events       - new copies of events.
 * pops             - population in which each event resides.
 * new_ages         - age of each new event.
 *
 * rubberband_migs is an array of size numLoci allocated in getMem().
 *---------------------------------------------------------------------------*/
typedef struct _RUBBERBAND_MIGS
{
  int num_moved_events;
  int orig_events[MAX_MIGS + MAX_MIG_BANDS];
  int new_events[MAX_MIGS + MAX_MIG_BANDS];
  int pops[MAX_MIGS + MAX_MIG_BANDS];
  double new_ages[MAX_MIGS + MAX_MIG_BANDS];
} RUBBERBAND_MIGS;

/*-----------------------------------------------------------------------------
 * MIG_SPR_STATS
 * Holds statistics for the SPR sampling operation with migration.
 * In use in UpdateGB_MigSPR and in traceLineage.
 *---------------------------------------------------------------------------*/
typedef struct _MIG_SPR_STATS
{
  int    father_event_old;
  int    father_event_new;
  int    father_pop_new;
  int    target;
  int    num_old_migs;
  int    num_new_migs;
  int    old_migs[MAX_MIGS];
  int    new_migs_in[MAX_MIGS];
  int    new_migs_out[MAX_MIGS];
  int    new_migs_bands[MAX_MIGS];
  double new_migs_ages[MAX_MIGS];
  double genetree_delta_lnLd[2];
} MIG_SPR_STATS;

/*-----------------------------------------------------------------------------
 * Locus_SuperStruct
 *---------------------------------------------------------------------------*/
typedef struct _Locus_SuperStruct
{
  GENETREE_STATS         genetree_stats_check;
  GENETREE_STATS_DELTA   genetree_stats_delta[NUM_DELTA_STATS_INSTANCES];
  MIG_SPR_STATS          mig_spr_stats;
  double                 genDeltaLogLikelihood;
  double                 genLogLikelihood;
  RUBBERBAND_MIGS        rubberband_migs;
  int                    mig_conflict_log;
} Locus_SuperStruct;

/*-----------------------------------------------------------------------------
 * ADMIXTURE_STATUS
 *---------------------------------------------------------------------------*/
typedef struct _ADMIXTURE_STATUS
{
  // number of loci for which to show the admixture status
  int numSampledLoci;

  // array of sampled loci
  int* sampledLoci;

  // array of admixture rates (up to the curent point of sampling)
  // for each locus and each admixed sample
  double** sampleLocusAdmixRate;

  // number of loci in which sample is in alternative population
  int* admixtureCounts;

  // estimated coefficients for all samples
  double* admixtureCoefficients;

} ADMIXTURE_STATUS;

/*-----------------------------------------------------------------------------
 * GENETREE_MIGS
 * A struct which contains information about
 * migration events in a specific genealogy.
 * Array of structs (of length numLoci is allocated in GetMem().
 *---------------------------------------------------------------------------*/
typedef struct _GENETREE_MIGS
{
  // number of migration event in genetree
  int num_migs;

  // array of indices for mignodes for
  // dynamic managing of migration nodes
  int living_mignodes[MAX_MIGS];

  struct MIGNODE
  {
    // id of gene tree node below the migration node
    int gtree_branch;

    // id of migration band relevant to this node
    int migration_band;

    // source and target populations of migration event (backwards view)
    int target_pop;
    int source_pop;

    // ids of event corresponding to this migration node.
    int target_event;
    int source_event;

    // time of event
    double age;
  } mignodes[MAX_MIGS];

} GENETREE_MIGS;

/*-----------------------------------------------------------------------------
 * GENETREE_STATS_FLAT
 * Holds relevant statistics for fast computation of probability of tree
 * given a null model (single population).
 * array allocated in GetMem().
 *---------------------------------------------------------------------------*/
typedef struct _GENETREE_STATS_FLAT
{
  double  coal_stats_flat;
  double  mig_stats_flat;
  int     num_coals_total;
  int     num_migs_total;
  double* sortedAgesArray;
} GENETREE_STATS_FLAT;

/*-----------------------------------------------------------------------------
 * GENETREE_NODE_STATS
 * Holds relevant statistics for coalescent distribution across populations
 * array allocated in GetMem().
 *---------------------------------------------------------------------------*/
typedef struct _GENETREE_NODE_STATS
{
  // --- Memory allocation ---
  double*   doubleArray;
  double**  doublePtrArray;
  double*** doublePtrPtrArray;
  int*      intArray;
  int**     intPtrArray;

  // --- Auxilliary matrices for computation ---
  // matrix of LCAs for all pairs of leaves
  int**     lcaMatrix;
  // array of leaves for LCA computation
  int*      leafArray;
  // array with the first coal node in each pop
  int*      firstNodesInPop;
  // array with the age of the first coal node in each pop
  double*   firstNodesAges;
  // array with the ages of all internal nodes
  double*   nodeAges;

  // --- Main stat matrices ---
  // Prob[sample pair coalesce in pop]
  double*** probCoalMatrix;
  // Prob[sample pair is the first coalescence in pop]
  double*** probFirstCoalMatrix;
  // mean coal time for coalescence of pair (given they coalesce in pop)
  double*** coalTimeMatrix;

} GENETREE_NODE_STATS;

/*-----------------------------------------------------------------------------
 *
 * FUNCTION DECLARATIONS
 *
 *---------------------------------------------------------------------------*/
int analyzeGenetreeFile( char*               genetree_file );
int writeMScommandLine ( char*               outfile );
int findLastMig(         int                 gen,
                         int                 node_id,
                         double              time );
int findFirstMig(        int                 gen,
                         int                 node_id,
                         double              time );


int Coalescence1Pop(     PopulationTree*     popTree,
                         GenericBinaryTree*  tree,
                         int                 gen,
                         int                 pop,
                         int*                livingLineages );
int removeEvent(         int                 gen,
                         int                 event );

// auxiliary functions

int GetRandomGtree(GenericBinaryTree* tree, int gen);
int adjustRootEvents();

int findInconsistency( int gen, int node );
int getSptreeNodeByName( const char* name );
int orderByAge( int subtree_root, int* ordered_nodes );
int getLineagesAtInterval( int gen, int start_event,
                           int pop, int exc_node,
                           int* out_array );
int getEdgesForTimePop (int gen, double time,
                        int pop, int exc_node,
                        int* out_array);
int
populationPostOrder (int pop, int* ordered_pops);

// event chain functions
int
checkAll ();
int
checkGtreeStructure (int gen);
int
synchronizeEvents (int gen);
int
printEventChains (FILE* stream, int gen);
int
createEvent (int gen, int pop, double age);
int
createEventBefore (int gen, int pop, int event, double elapsed_time);
int
constructEventChain (int gen);
int
computeFlatStats ();
int
computeNodeStats ();
int
computeTotalStats ();
double
recalcStats (int gen, int pop);
int
recalcStats_partitioned (int gen, int pop);
int
computeGenetreeStats (int gen);
int
computeGenetreeStats_partitioned (void);
double
gtreeLnLikelihood (int gen);
double
computeDeltaLnLd (int gen, int instance);
int
computeMigStatsDelta (int instance, double bottom_age, int bottom_pop,
                      double top_age, int num_lins_delta, int gen);
int
computeCoalStatsDelta (int instance, int gen, int bottom_event, int bottom_pop,
                       int top_event, int num_lins_delta);
double
considerEventMove (int gen, int instance, int event_id, int source_pop,
                   double original_age, int target_pop, double new_age);
int
acceptEventChainChanges (int gen, int instance);
int
rejectEventChainChanges (int gen, int instance);
double
rubberBandRipple (int gen, int do_or_redo);
double
rubberBand (int gen, int pop, double static_point, double moving_point,
            double factor, unsigned short postORpre, int* out_num_events);
int
replaceMigNodes (int gen, int node);

//============================ END OF FILE ====================================
