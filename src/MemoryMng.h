#pragma once
/*============================================================================
 File: MemoryMng.h

 Basic memory management module - allocating and freeing.

 Version history:
 03-May-2017  evgenyidc      Extracting from patch.cpp
 ============================================================================*/

#include "patch.h"
#include "DataLayer.h"

/*-----------------------------------------------------------------------------
 *
 * Function declarations
 *
 *---------------------------------------------------------------------------*/

int GetMem();
int FreeMem();

/*-----------------------------------------------------------------------------
 *
 * Global data structures
 *
 *---------------------------------------------------------------------------*/
extern GENETREE_NODE_STATS genetree_node_stats;
extern GENETREE_STATS*     genetree_stats;
extern GENETREE_STATS*     genetree_stats_total_partitioned;
extern GENETREE_STATS      genetree_stats_total;
extern GENETREE_STATS      genetree_stats_total_check;
extern GENETREE_MIGS*      genetree_migs;
extern GENETREE_STATS_FLAT genetree_stats_flat;
extern Locus_SuperStruct*  locus_data;

// Node surrogates.
// a 2D array (numLoci X numNodes) for populations per genealogy node.
extern int**               nodePops;
// a 2D array (numLoci X numNodes) for event ids per genealogy node.
extern int**               nodeEvents;


//============================ END OF FILE ====================================
