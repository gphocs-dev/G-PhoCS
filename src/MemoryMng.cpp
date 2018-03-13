/*============================================================================
 File: MemoryMng.cpp

 Basic memory management module - allocating and freeing.

 Version history:
 03-May-2017  evgenyidc      Extracting from patch.cpp
 13-Jun-2017  evgenyidc      Cleaning GetMem. Adding XDag wrapper.
 ============================================================================*/
#include "PopulationTree.h"
#include "patch.h"
#include "MemoryMng.h"
#include "DataLayer.h"
#include "MCMCcontrol.h"
#include <stdlib.h>
#include <stdio.h>
#include "EventsDAG.h"


/*-----------------------------------------------------------------------------
 *
 * Global Data Structures
 *
 *---------------------------------------------------------------------------*/

GENETREE_MIGS*       genetree_migs;
GENETREE_STATS*      genetree_stats;
GENETREE_STATS*      genetree_stats_total_partitioned;
GENETREE_STATS       genetree_stats_total;
GENETREE_STATS       genetree_stats_total_check;
GENETREE_STATS_FLAT  genetree_stats_flat;
GENETREE_NODE_STATS  genetree_node_stats;
int**                nodePops;
int**                nodeEvents;
Locus_SuperStruct*   locus_data;
extern DAGsPerLocus<Event>* pAllDAGs;


/*-----------------------------------------------------------------------------
 *
 * GetMem
 *
 *---------------------------------------------------------------------------*/

int GetMem( void )
{
  int gen, count, i, maxNodes = 2*dataSetup.numSamples-1;
  int leaf1, leaf2;

  // for debugging purposes   ELIMINATE LATER !!!!

  locus_data = (Locus_SuperStruct*)
               malloc( dataSetup.numLoci * sizeof(Locus_SuperStruct) );
  if(locus_data == NULL)
  {
    fprintf(stderr, "\nError: Out Of Memory genLogLikelihood.\n");
    exit(-1);
  }

  // get mem for auxiliary node arrays
  nodePops = (int**)
             malloc( 2 * dataSetup.numLoci * sizeof(int*) );
  if(nodePops == NULL)
  {
    fprintf(stderr, "\nError: Out Of Memory nodePop array.\n");
    exit(-1);
  }
  nodeEvents  = nodePops + dataSetup.numLoci;
  nodePops[0] = (int*)
                malloc( 2 * dataSetup.numLoci * maxNodes * sizeof(int) );
  if( nodePops[0] == NULL )
  {
    fprintf(stderr, "\nError: Out Of Memory nodePop space.\n");
    exit(-1);
  }

  for( gen = 0; gen < dataSetup.numLoci; ++gen )
  {
    nodePops[gen]   = nodePops[0]   + maxNodes * gen;
    nodeEvents[gen] = nodePops[gen] + maxNodes * dataSetup.numLoci;
  }

  // get memory for genetree_migs and initialize
  // get memory for event_chains and stats and initialize
  genetree_migs = (GENETREE_MIGS*)
                  malloc( dataSetup.numLoci * sizeof(GENETREE_MIGS) );
  if(genetree_migs == NULL)
  {
    fprintf(stderr, "\nError: Out Of Memory genetree migs.\n");
    exit(-1);
  }

  genetree_stats = (GENETREE_STATS*)
                   malloc( dataSetup.numLoci * sizeof(GENETREE_STATS) );
  if(genetree_stats == NULL)
  {
    fprintf(stderr, "\nError: Out Of Memory genetree stats.\n");
    exit(-1);
  }
  genetree_stats_total_partitioned=(GENETREE_STATS*)malloc(dataSetup.numPopPartitions*sizeof(GENETREE_STATS));
  if(genetree_stats_total_partitioned == NULL) {
    fprintf(stderr, "\nError: Out Of Memory genetree stats total partitioned.\n");
    exit(-1);
  }

  event_chains.reserve( dataSetup.numLoci );
  pAllDAGs = new DAGsPerLocus<Event>( dataSetup.numLoci, dataSetup.popTree->numPops);

  genetree_stats_flat.sortedAgesArray = (double*)
                                        malloc( maxNodes * sizeof(double) );
  if(genetree_stats_flat.sortedAgesArray == NULL)
  {
    fprintf(stderr, "\nError: Out Of Memory event chains.\n");
    exit(-1);
  }


  // max number of events should cover all possible
  // coalescences (+1 auxiliary),
  // migrations (X2), added migrations (X2) and population endings,
  // and migration bands (start + end + changed event).
  event_chains[0].total_events  =    2 * dataSetup.numSamples
                                   + 4 * MAX_MIGS
                                   + 3 * dataSetup.popTree->numMigBands
                                   + dataSetup.popTree->numPops
                                   + 10;
  for( gen = 1; gen < dataSetup.numLoci; ++gen )
    event_chains[gen].total_events = event_chains[0].total_events;

  // compute max number of events for all loci
  count = event_chains[0].total_events * dataSetup.numLoci;

  event_chains[0].events = new Event[count];
  if(event_chains[0].events == NULL)
  {
    fprintf(stderr, "\nError: Out Of Memory event space.\n");
    exit(-1);
  }

  for( gen = 0; gen < dataSetup.numLoci; ++gen )
  {
    if( gen > 0 )
    {
      event_chains[gen].events =   event_chains[gen-1].events
                                 + event_chains[gen-1].total_events;
    }
    genetree_migs[gen].num_migs = 0;
    //initialize mignodes
    for( i = 0; i < MAX_MIGS; ++i)
    {
      genetree_migs[gen].mignodes[i].age            =  0;
      genetree_migs[gen].mignodes[i].target_pop     = -1;
      genetree_migs[gen].mignodes[i].source_pop     = -1;
      genetree_migs[gen].mignodes[i].gtree_branch   = -1;
      genetree_migs[gen].mignodes[i].source_event   = -1;
      genetree_migs[gen].mignodes[i].target_event   = -1;
      genetree_migs[gen].mignodes[i].migration_band = -1;
    }


    // initialize genetree_stats_delta
    for( i = 0; i < NUM_DELTA_STATS_INSTANCES; ++i)
      locus_data[gen].genetree_stats_delta[i].init();

  }
  // initialize genetree_node_stats

  count = 3 * dataSetup.numSamples
            * dataSetup.numSamples
            * dataSetup.popTree->numPops
          + maxNodes
          + dataSetup.popTree->numPops;
  genetree_node_stats.doubleArray = (double*)
                                    malloc( count * sizeof(double) );
  if( genetree_node_stats.doubleArray == NULL )
  {
    fprintf(stderr, "\nError: Out Of Memory genetree node stats - doubles.\n");
    exit(-1);
  }

  count = 3 * dataSetup.numSamples
            * dataSetup.numSamples;
  genetree_node_stats.doublePtrArray = (double**)
                                       malloc( count * sizeof(double*) );
  if(genetree_node_stats.doublePtrArray == NULL)
  {
    fprintf( stderr,
             "\nError: Out Of Memory genetree node stats "
             "- double pointers.\n");
    exit( -1 );
  }

  count = 3 * dataSetup.numSamples;
  genetree_node_stats.doublePtrPtrArray = (double***)
                                          malloc( count * sizeof(double*) );
  if(genetree_node_stats.doublePtrPtrArray == NULL)
  {
    fprintf( stderr,
             "\nError: Out Of Memory genetree node stats "
             "- double pointers^2.\n");
    exit( -1 );
  }

  count =   dataSetup.numSamples
          * dataSetup.numSamples
          + dataSetup.numSamples
          + dataSetup.popTree->numPops;
  genetree_node_stats.intArray = (int*)
                                 malloc( count * sizeof(int) );
  if(genetree_node_stats.intArray == NULL)
  {
    fprintf(stderr, "\nError: Out Of Memory genetree node stats - ints.\n");
    exit( -1 );
  }

  count = dataSetup.numSamples;
  genetree_node_stats.intPtrArray = (int**)
                                    malloc( count * sizeof(int*) );
  if(genetree_node_stats.intPtrArray == NULL)
  {
    fprintf( stderr,
             "\nError: Out Of Memory genetree node stats - int pointers.\n" );
    exit( -1 );
  }

  genetree_node_stats.lcaMatrix           =  genetree_node_stats.intPtrArray;
  genetree_node_stats.probCoalMatrix      =  genetree_node_stats.doublePtrPtrArray;
  genetree_node_stats.probFirstCoalMatrix =  genetree_node_stats.doublePtrPtrArray +      dataSetup.numSamples;
  genetree_node_stats.coalTimeMatrix      =  genetree_node_stats.doublePtrPtrArray +  2 * dataSetup.numSamples;

  for( leaf1 = 0; leaf1 < dataSetup.numSamples; ++leaf1 )
  {
    genetree_node_stats.lcaMatrix[leaf1]           = genetree_node_stats.intArray              +     leaf1 * dataSetup.numSamples;
    genetree_node_stats.probCoalMatrix[leaf1]      = genetree_node_stats.doublePtrArray        + 3 * leaf1 * dataSetup.numSamples;
    genetree_node_stats.probFirstCoalMatrix[leaf1] = genetree_node_stats.probCoalMatrix[leaf1] +             dataSetup.numSamples;
    genetree_node_stats.coalTimeMatrix[leaf1]      = genetree_node_stats.probCoalMatrix[leaf1] + 2 *         dataSetup.numSamples;
    for( leaf2 = 0; leaf2 < dataSetup.numSamples; ++leaf2)
    {
      genetree_node_stats.probCoalMatrix     [leaf1][leaf2] = genetree_node_stats.doubleArray + 3 * dataSetup.popTree->numPops * (leaf1 * dataSetup.numSamples + leaf2);
      genetree_node_stats.probFirstCoalMatrix[leaf1][leaf2] = genetree_node_stats.probCoalMatrix[leaf1][leaf2] +     dataSetup.popTree->numPops;
      genetree_node_stats.coalTimeMatrix     [leaf1][leaf2] = genetree_node_stats.probCoalMatrix[leaf1][leaf2] + 2 * dataSetup.popTree->numPops;
    }
  }
  genetree_node_stats.leafArray       = genetree_node_stats.intArray       +     dataSetup.numSamples * dataSetup.numSamples;
  genetree_node_stats.firstNodesInPop = genetree_node_stats.leafArray      +     dataSetup.numSamples;
  genetree_node_stats.firstNodesAges  = genetree_node_stats.doubleArray    + 3 * dataSetup.popTree->numPops * dataSetup.numSamples * dataSetup.numSamples;
  genetree_node_stats.nodeAges        = genetree_node_stats.firstNodesAges +     dataSetup.popTree->numPops;

  return 0;
}
/** end of GetMem **/


/*-----------------------------------------------------------------------------
 *
 * FreeMem
 *
 *---------------------------------------------------------------------------*/
int FreeMem (void)
{
  // free(genLogLikelihood);
  free(nodePops[0]);
  free(nodePops);
  free(genetree_migs);
  free(genetree_stats);
  free(genetree_stats_total_partitioned);
  //free(rubberband_migs);
  free(event_chains[0].events);
  //free(event_chains); //done by STL
  free(genetree_stats_flat.sortedAgesArray);
  free(genetree_node_stats.doubleArray);
  free(genetree_node_stats.doublePtrArray);
  free(genetree_node_stats.doublePtrPtrArray);
  free(genetree_node_stats.intArray);
  free(genetree_node_stats.intPtrArray);

  free(locus_data);
  delete pAllDAGs;
  return 0;
}

//============================ END OF FILE ====================================

