#pragma once
/*============================================================================
 File: DbgErrMsgGTreeStruct.h

 Code sections dealing with error handling are moved here as macros.
 The handling is in GTree Structure validation code.

 Version history:
 10-May-2017  evgenyidc      Extracting error handling code from patch.cpp.
  ============================================================================*/

//----------------------------------------------------------------------------
#define CHECK_GTREE_STRUCT_FATAL_0034 \
if(debug)\
{\
  fprintf(stderr, "\nError: checking genetree for gen %d: ", gen);\
  fprintf(stderr, "number of lineages of event_idx %d doesn't "\
                  "match: %d, %d.\n", stack_vars.event_idx,\
                                      pCurrEvent->getNumLineages(),\
                                      num_lineages);\
}\
else\
{\
  fprintf(stderr, "Fatal Error 0034.\n");\
}\
stack_vars.res = 0;

//----------------------------------------------------------------------------
#define CHECK_GTREE_STRUCT_FATAL_0035 \
if(debug) \
{ \
  fprintf(stderr, "\nError: checking genetree for gen %d: ",gen); \
  fprintf(stderr, "prev event_idx of next event_idx" \
                  "%d doesn't match: %d, %d.\n", \
                  event_next_idx, stack_vars.event_idx, \
                  event_next_prev_idx); \
} \
else \
{ \
  fprintf(stderr, "Fatal Error 0035.\n"); \
}\
stack_vars.res = 0;


//----------------------------------------------------------------------------
#define CHECK_GTREE_STRUCT_FATAL_0036 \
if(debug)\
{\
  fprintf(stderr, "\nError: checking genetree for gen %d: ",gen);\
  fprintf(stderr, "sample age of pop %d (event_idx %d) doesn't"\
                  " match: %f, %f.\n", pop, event_idx,\
                  age,\
                  dataSetup.popTree->pops[pop]->sampleAge);\
}\
else\
{\
  fprintf(stderr, "Fatal Error 0036.\n");\
}\
res = 0;

//----------------------------------------------------------------------------
#define CHECK_GTREE_STRUCT_FATAL_0036_2 \
if(debug)\
{\
  fprintf(stderr, "\nError: checking genetree for gen %d: ",gen);\
  fprintf(stderr, "age of node %d (event_idx %d) doesn't match:"\
                  " %g, %g.\n", event_id, event_idx,\
                  age, curr_nodeAge );\
}\
else\
{\
  fprintf(stderr, "Fatal Error 0036.\n");\
}\
res = 0;

//----------------------------------------------------------------------------
#define CHECK_GTREE_STRUCT_FATAL_0037 \
if(debug)\
{\
  fprintf(stderr, "\nError: checking genetree for gen %d: ",gen);\
  fprintf(stderr, "population of node %d (event_idx %d) doesn't "\
                  "match: %d, %d.", p_stack->event_id,\
                  p_stack->event_idx, pop,\
                  nodePops[gen][p_stack->event_id]);\
}\
else\
{\
  fprintf(stderr, "Fatal Error 0037.\n");\
}\
p_stack->res = 0;

//----------------------------------------------------------------------------
#define CHECK_GTREE_STRUCT_FATAL_0038 \
if(debug)\
{\
  fprintf(stderr, "\nError: checking genetree for gen %d: ",gen);\
  fprintf(stderr, "event_idx of node %d doesn't match: %d, %d.",\
                   event_id, event_idx, nodeEvents[gen][event_id]);\
}\
else\
{\
  fprintf(stderr, "Fatal Error 0038.\n");\
}\
res = 0;

//----------------------------------------------------------------------------
#define CHECK_GTREE_STRUCT_FATAL_0039a \
if(debug)\
{\
  fprintf(stderr, "\nError: checking genetree for gen %d: ",gen);\
  fprintf(stderr, "mig band %d invalid for mig node %d .",\
                  mig_band, event_id);\
}\
else\
{\
  fprintf(stderr, "Fatal Error 0039a.\n");\
}\
p_stack->res = 0;

//----------------------------------------------------------------------------
#define CHECK_GTREE_STRUCT_FATAL_0056 \
if(debug) \
{ \
  fprintf(stderr, "\nError: checking genetree for gen %d: number of " \
                  "live mig bands %d at end of population %d.\n", \
          gen, stack_vars.num_living_mig_bands, pop); \
} \
else \
{ \
  fprintf(stderr, "Fatal Error 0056.\n"); \
}\
printGenealogyAndExit(gen,-1);

//----------------------------------------------------------------------------
#define CHECK_GTREE_STRUCT_FATAL_0057 \
if(debug)\
{\
  fprintf(stderr, "\nError: checking genetree for gen %d: ",gen);\
  fprintf(stderr, "coal stats for pop %d don't match: %g, %g "\
                  "(diff = %g)", pop,\
          genetree_stats[gen].coal_stats[pop],\
          locus_data[gen].genetree_stats_check.coal_stats[pop],\
          genetree_stats[gen].coal_stats[pop]-\
            locus_data[gen].genetree_stats_check.coal_stats[pop]);\
}\
else\
{\
  fprintf(stderr, "Fatal Error 0057.\n");\
}

//----------------------------------------------------------------------------
#define CHECK_GTREE_STRUCT_FATAL_0058 \
if(debug)\
{\
  fprintf(stderr, "\nError: checking genetree for gen %d: ",gen);\
  fprintf(stderr, "num coals for pop %d don't match: %d, %d.",\
          pop, genetree_stats[gen].num_coals[pop],\
          locus_data[gen].genetree_stats_check.num_coals[pop]);\
}\
else\
{\
  fprintf(stderr, "Fatal Error 0058.\n");\
}
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

//============================ END OF FILE ====================================
