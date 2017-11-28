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
#include "EventsDAG.h"

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
synchronizeEvents (int gen);
int
printEventChains (FILE* stream, int gen);
int
createEvent (int gen, int pop, double age);
int
createEventBefore (int gen, int pop, int event, double elapsed_time);
int
createEventBefore( int                   nGenIdx,
                   int                   nPopIdx,
                   EventsDAGNode<Event>* pCurrEvent,
                   double                fElapsedTime );

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
computeGenetreeStats (int gen);
//int
//recalcStats_partitioned (int gen, int pop);
//int
//computeGenetreeStats_partitioned (void);
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
