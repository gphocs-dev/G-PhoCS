#pragma once
/*============================================================================
 File: DbgErrMsgTraceLineages.h

 Code sections dealing with error handling are moved here as macros.

 Version history:
 05-Apr-2017  evgenyidc      Extracting error handling code from traceLineage.
 19-Apr-2017  evgenyidc      Wrapping traceLineage auto vars as a struct.
 10-May-2017  evgenyidc      Name change. Only TraceLineages is covered here.
  ============================================================================*/

//----------------------------------------------------------------------------
#define TRACE_LINEAGE_INIT_FATAL_0101 \
if(debug) \
{\
  fprintf(stderr, "\nError: traceLineage: while tracing edge above"\
					        " node %d in gen %d.\n",node,gen);\
  fprintf(stderr, "  Could not find start event in pop %d.\n",p_stack->pop);\
  fprintf(stderr, ".\n");\
} \
else\
{\
  fprintf(stderr, "Fatal Error 0101.\n");\
}

//----------------------------------------------------------------------------
#define TRACE_LINEAGE_FATAL_0006 \
if(debug) \
{ \
	fprintf(stderr, "\nError: traceLineage: while tracing "\
                  "edge above node %d in gen %d.\n",\
									p_stack->node, p_stack->gen);\
	fprintf(stderr, "  Reached top event at age %f "\
					        "(reconnect == %d).\n",\
									p_stack->age, p_stack->reconnect);\
	fprintf(stderr, "  Events visited:");\
	for(i=0; \
	    i < locus_data[p_stack->gen].\
			    genetree_stats_delta[p_stack->reconnect].num_changed_events;\
	    i++)\
	{\
		fprintf(stderr, " %d",\
           locus_data[p_stack->gen].\
					 genetree_stats_delta[p_stack->reconnect].changed_events[i]);\
	}\
	fprintf(stderr, ".\n");\
}\
else\
{\
	fprintf(stderr, "Fatal Error 0006.\n");\
}\
printGenealogyAndExit(p_stack->gen,-1);

//----------------------------------------------------------------------------
#define TRACE_LINEAGE_FATAL_0007 \
if(debug) \
{\
  fprintf(stderr, "\Warning: traceLineage: while tracing"\
					        " edge above node %d in gen %d.\n", \
                  p_stack->node, p_stack->gen);\
	fprintf(stderr, "  %d living migration bands when moving to population %d"\
					        " (with total rate %f):\n",\
                  p_stack->num_live_mig_bands, \
                  p_stack->pop, p_stack->mig_rate);\
  for( i=0; i < p_stack->num_live_mig_bands; i++)\
  {\
    fprintf(stderr, " %d (rate %f)", p_stack->live_mig_bands[i],\
						dataSetup.popTree->migBands[ p_stack->live_mig_bands[i] ].migRate);\
  }\
	fprintf(stderr, ".\n");\
}\
else\
{\
  fprintf(stderr, "Warning 0007.\n");\
}
//----------------------------------------------------------------------------
#define TRACE_LINEAGE_FATAL_0008 \
if(debug)\
{\
  fprintf(stderr, "\nError traceLineage: while tracing edge above"\
					        " node %d in gen %d.\n", p_stack->node, p_stack->gen);\
  fprintf(stderr, "  age recorded at start of pop %d is %g, while age of pop"\
					        " is %g (diff = %g).\n", p_stack->pop, p_stack->age, \
									dataSetup.popTree->pops[p_stack->pop]->age,\
									p_stack->age - dataSetup.popTree->pops[p_stack->pop]->age);\
  for( i = 0; i < p_stack->num_live_mig_bands; i++ ) \
  {\
    fprintf(stderr, " %d (rate %f)",p_stack->live_mig_bands[i],\
						dataSetup.popTree->migBands[p_stack->live_mig_bands[i]].migRate);\
  }\
  fprintf(stderr, ".\n");\
}\
else\
{\
  fprintf(stderr, "Fatal Error 0008.\n");\
}\
printGenealogyAndExit(p_stack->gen,-1);

//----------------------------------------------------------------------------
#define DEBUG_TRACE_LINEAGE_ASSERT_EVENT_SAMPLE_IS_ZERO \
if(p_stack->event_sample == 0.0 && debug)\
{\
  fprintf(stderr, "\ntraceLineage() event_sample=%g when rate=%g at"\
					        " pop=%d, age=%g.\n",p_stack->event_sample,\
					        p_stack->rate,p_stack->pop,p_stack->age);\
 	printPopulationTree(dataSetup.popTree, stderr, 1);	\
  printLocusGenTree(dataState.lociData[p_stack->gen], stderr,\
										nodePops[p_stack->gen], \
                    nodeEvents[p_stack->gen]);\
  printEventChains(stderr, p_stack->gen);\
	fprintf(stderr, "\n\n*******************************"\
					        "*************************************\n\n");\
}

//----------------------------------------------------------------------------
#define TRACE_LINEAGE_FATAL_0009 \
fprintf(stderr, "\ntraceLineage() after iterating through all migration "\
				        "bands, event_sample=%g when rate=%g, mig_rate=%g at pop=%d,"\
								" age=%g.\n",	p_stack->event_sample, p_stack->rate, \
								p_stack->mig_rate, p_stack->pop, p_stack->age); \
if(p_stack->num_live_mig_bands <= 0)\
{\
  if(debug)\
  {\
	  fprintf(stderr, "\nError: traceLineage: while re-coalescing edge above "\
						        "node %d in gen %d:\n", p_stack->node, p_stack->gen);\
    fprintf(stderr, "  num_live_mig_bands=%d but migration rate is %g " \
						        "(event_sample = %g, pop=%d, age=%g).\n",\
					          p_stack->num_live_mig_bands, p_stack->mig_rate, \
					          p_stack->event_sample, p_stack->pop, p_stack->age);\
  } \
	else \
	{\
    fprintf(stderr, "Fatal Error 0009.\n");\
  }\
  printGenealogyAndExit(p_stack->gen,-1);\
}\
fprintf(stderr, "%d live mig bands:",p_stack->num_live_mig_bands);\
for( i=0, p_stack->event_sample = p_stack->mig_rate; \
     i < p_stack->num_live_mig_bands; i++) \
{\
	fprintf(stderr, " %d (rate=%g)", p_stack->live_mig_bands[i],\
          dataSetup.popTree->migBands[ p_stack->live_mig_bands[i] ].migRate);\
	p_stack->event_sample -= \
	         dataSetup.popTree->migBands[ p_stack->live_mig_bands[i] ].migRate;\
}\
fprintf(stderr, " mig_rate - all rates = %g.\n", p_stack->event_sample);\
printPopulationTree(dataSetup.popTree, stderr, 1);\
printLocusGenTree(dataState.lociData[p_stack->gen], stderr, \
									nodePops[p_stack->gen], nodeEvents[p_stack->gen]);\
printEventChains(stderr, p_stack->gen);\
fprintf(stderr, "\n\n****************************"\
				        "****************************************\n\n");
//----------------------------------------------------------------------------
#define TRACE_LINEAGE_FATAL_0009A \
if(debug) {\
	fprintf(stderr, "\nError: traceLineage: while re-coalescing edge "\
          "above node %d in gen %d:\n", p_stack->node, p_stack->gen);\
	fprintf(stderr, "  i=%d and event_sample=%g (migration rate = %g,"\
	        " pop=%d, age=%g).\n", \
		i,p_stack->event_sample,p_stack->mig_rate,p_stack->pop,p_stack->age);\
} else {\
	fprintf(stderr, "Fatal Error 0009A.\n");\
}\
printGenealogyAndExit(p_stack->gen,-1);

//----------------------------------------------------------------------------
#define TRACE_LINEAGE_FATAL_0009B \
if(debug) {\
	fprintf(stderr, "\nError: traceLineage: while re-coalescing "\
	        "edge above node %d in gen %d:\n", p_stack->node, p_stack->gen);\
	fprintf(stderr, "  selected mig_band=%d but target population %d does"\
	        " not fit current population %d (age=%g).\n",\
		      p_stack->mig_band, \
		      dataSetup.popTree->migBands[p_stack->mig_band].targetPop, \
		      p_stack->pop, p_stack->age);\
	fprintf(stderr, "%d live mig bands:", p_stack->num_live_mig_bands);\
	for( i = 0, p_stack->event_sample = p_stack->mig_rate; \
       i < p_stack->num_live_mig_bands; i++) {\
		fprintf(stderr, " %d (rate=%g)", p_stack->live_mig_bands[i], \
		       dataSetup.popTree->migBands[ p_stack->live_mig_bands[i] ].migRate);\
		p_stack->event_sample -= dataSetup.popTree->\
		                         migBands[ p_stack->live_mig_bands[i] ].migRate;\
	}\
	fprintf(stderr, " mig_rate - all rates = %g.\n", p_stack->event_sample);\
} else {\
	fprintf(stderr, "Fatal Error 0009B.\n");\
}\
printGenealogyAndExit(p_stack->gen,-1);

//----------------------------------------------------------------------------
#define TRACE_LINEAGE_FATAL_0010 \
if(debug) {\
	fprintf(stderr, "Error: problem encountered when creating new migration "\
	        "event for locus %d, mig band %d pops %d-->%d, in trace lineage.\n",\
					p_stack->gen, p_stack->mig_band, \
					dataSetup.popTree->migBands[p_stack->mig_band].sourcePop, \
					dataSetup.popTree->migBands[p_stack->mig_band].targetPop );\
} else {\
	fprintf(stderr, "Fatal Error 0010.\n");\
}\
printGenealogyAndExit(p_stack->gen,-1);
//----------------------------------------------------------------------------
#define TRACE_LINEAGE_FATAL_0011 \
if(debug) {\
fprintf(stderr, "\nError: in traceLineage: %d targets and %d lineages "\
        "(at gen %d, node %d, pop %d, , event %d, age %g)\n",\
			  p_stack->num_targets, \
			  event_chains[p_stack->gen].events[p_stack->event].getNumLineages(), \
			  p_stack->gen, p_stack->node, p_stack->pop, p_stack->event,\
			 ( p_stack->age - p_stack->t ) \
			   + event_chains[p_stack->gen].events[p_stack->event].getElapsedTime() \
       / 2.0);\
} else {\
fprintf(stderr, "Fatal Error 0011.\n");\
}\
printGenealogyAndExit(p_stack->gen,-1);

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

//============================ END OF FILE ====================================
