/*============================================================================
 File: TraceLineages.cpp

 Trace Lineages functionality

 Version history:
 07-May-2017  evgenyidc    Extracting from patch.cpp.
 ============================================================================*/

#include "DataLayerConstants.h"
#include "TraceLineages.h"
#include "MemoryMng.h"
#include "DataLayer.h"

#include <stdio.h>
#include <math.h>

#include "GPhoCS.h"

#include "DbgErrMsgTraceLineages.h"
extern int debug;


//======================== Trace Lineages Code ===============================
EventIvgeny* get_first_event_in_pop(int gen, int pop)
{
  return &(event_chains[gen].events[event_chains[gen].first_event[pop]]);
}

EventIvgeny* get_next_event(int gen, EventIvgeny* pCurrEvent)
{
  int next_idx = pCurrEvent->getNextIdx();
  if( -1 != next_idx)
    return &(event_chains[gen].events[next_idx]);
  else
    return nullptr;
}

//----------------------------------------------------------------------------
void traceLineage_init(TraceLineageAutoVars* p_stack)
{
  int i = 0;

  int& gen                 = p_stack->gen;
  int& node                = p_stack->node;
  int& reconnect           = p_stack->reconnect;
  int& pop                 = p_stack->pop;
  p_stack->age             = -1.0;
  p_stack->heredity_factor = 1.0;

  /******** initialization start  **************/
  p_stack->pop = nodePops[gen][node];

  if(node < dataSetup.numSamples)
  {
    p_stack->pEvent = get_first_event_in_pop(gen, pop);
    while( !p_stack->pEvent->isOfType(SAMPLES_START|END_CHAIN ) )
    {
      p_stack->pEvent = get_next_event(gen, p_stack->pEvent);
    }
    if( p_stack->pEvent->getType() == END_CHAIN )
    {
      TRACE_LINEAGE_INIT_FATAL_0101
    }
    p_stack->pEvent = get_next_event(gen, p_stack->pEvent);
  }
  else
  {
    p_stack->pEvent = get_next_event(gen,
                        &(event_chains[gen].events[nodeEvents[gen][node]]));
  }
  p_stack->theta =   dataSetup.popTree->pops[p_stack->pop]->theta
                   * p_stack->heredity_factor;

  p_stack->age = getNodeAge(dataState.lociData[gen], node);
  //age = nodes[node].age;
  locus_data[gen].mig_spr_stats.genetree_delta_lnLd[reconnect] = 0.0;
  if(!reconnect)
  {
    locus_data[gen].mig_spr_stats.num_old_migs = 0;
    if(node != getLocusRoot(dataState.lociData[gen]))
    {
      locus_data[gen].mig_spr_stats.father_event_old =
          nodeEvents[gen][ getNodeFather(dataState.lociData[gen], node) ];
    }
    // locus_data[gen].mig_spr_stats.father_event_old =
    //    nodes[nodes[node].father].event_id;
  }
  else
  {
    locus_data[gen].mig_spr_stats.num_new_migs = 0;
  }


  locus_data[gen].genetree_stats_delta[reconnect].clear_changed_events();
  // assume all migration bands and populations are affected
  locus_data[gen].genetree_stats_delta[reconnect].num_pops_changed =
                                                   dataSetup.popTree->numPops;
  for(i = 0; i < dataSetup.popTree->numPops; ++i)
  {
    locus_data[gen].genetree_stats_delta[reconnect].pops_changed[i] = i;
    locus_data[gen].genetree_stats_delta[reconnect].coal_stats_delta[i] = 0.0;
  }
  locus_data[gen].genetree_stats_delta[reconnect].num_mig_bands_changed =
                                                dataSetup.popTree->numMigBands;
  for(i = 0; i < dataSetup.popTree->numMigBands; ++i)
  {
    locus_data[gen].genetree_stats_delta[reconnect].mig_bands_changed[i] = i;
    locus_data[gen].genetree_stats_delta[reconnect].mig_stats_delta[i] = 0.0;
  }

  p_stack->mig_rate = 0.0;
  p_stack->num_live_mig_bands = 0;
  int& mb = p_stack->mig_band;
  for(mb = 0; mb < dataSetup.popTree->numMigBands; ++mb)
  {
    if(    dataSetup.popTree->migBands[mb].targetPop == p_stack->pop
        && dataSetup.popTree->migBands[mb].startTime  < p_stack->age
        && dataSetup.popTree->migBands[mb].endTime    > p_stack->age )
    {
      p_stack->mig_rate += dataSetup.popTree->migBands[mb].migRate;
      p_stack->live_mig_bands[p_stack->num_live_mig_bands++] = mb;
    }
  }
  p_stack->mig_source = -1;
  p_stack->proceed    =  1;
  /******** initialization end    **************/
}

//----------------------------------------------------------------------------
int
traceLineage_move_to_next_pop (TraceLineageAutoVars* p_stack)
{
  int i = 0;
  if (dataSetup.popTree->pops[p_stack->pop]->father == NULL)
  {
    // reject if trying to reconnect above OLDAGE
    if (p_stack->reconnect == 1)
    {
      return -1;
    }
    // if detaching subtree, you can reach top only if detaching the root edge
    // this has to be considered in migration scenarios
    TRACE_LINEAGE_FATAL_0006
  }
  p_stack->pop = dataSetup.popTree->pops[p_stack->pop]->father->id;
  p_stack->theta = dataSetup.popTree->pops[p_stack->pop]->theta
                   * p_stack->heredity_factor;
  p_stack->pEvent = get_first_event_in_pop(p_stack->gen, p_stack->pop);

  if (fabs (p_stack->mig_rate) > 0.00000001 || p_stack->num_live_mig_bands > 0)
  {
    TRACE_LINEAGE_FATAL_0007
  }
  p_stack->mig_rate = 0.0;

  if (fabs (p_stack->age / dataSetup.popTree->pops[p_stack->pop]->age - 1)
          > 0.01)
  {
    TRACE_LINEAGE_FATAL_0008
  }
  p_stack->age = dataSetup.popTree->pops[p_stack->pop]->age;
  return 0;
}

//----------------------------------------------------------------------------
void
traceLineage_reduce_lineages_along_pruned_edge (TraceLineageAutoVars* p_stack)
{
  // if(!reconnect) then reduce the number of lineages along pruned edge by 1
  p_stack->pEvent->decrementLineages ();
  p_stack->t = p_stack->pEvent->getElapsedTime ();
  p_stack->age += p_stack->t;
  //@@TODO: ugly
  p_stack->proceed = (p_stack->pEvent !=
                      &(event_chains[p_stack->gen]\
                         .events[ locus_data[p_stack->gen]\
                            .mig_spr_stats.father_event_old]));

  // if event is an in-migration event of edge above node,
  // record it and follow it to source population
  if( p_stack->pEvent->getType () == IN_MIG )
  {
    if (genetree_migs[p_stack->gen].mignodes[p_stack->node_id].gtree_branch
            == p_stack->node)
    {
      p_stack->mig_band =
        genetree_migs[p_stack->gen].mignodes[p_stack->node_id].migration_band;
      p_stack->mig_source =
        genetree_migs[p_stack->gen].mignodes[p_stack->node_id].source_event;
      locus_data[p_stack->gen].mig_spr_stats.\
         old_migs[locus_data[p_stack->gen].mig_spr_stats.num_old_migs++] =
           p_stack->node_id;
    }
  }
}

//----------------------------------------------------------------------------
int
traceLineage_sample_migration_event_in_interval(TraceLineageAutoVars* p_stack)
{
  int i;
  int& gen = p_stack->gen;

  // migration event - figure out where to migrate
  if(MAX_MIGS <=
        genetree_migs[gen].num_migs
      + locus_data[gen].mig_spr_stats.num_new_migs
      - locus_data[gen].mig_spr_stats.num_old_migs)
  {
    misc_stats.not_enough_migs++;
    return -1;
  }
  DEBUG_TRACE_LINEAGE_ASSERT_EVENT_SAMPLE_IS_ZERO
  for( i = 0;
       p_stack->event_sample >= 0 && i < p_stack->num_live_mig_bands;
       ++i)
  {
    p_stack->event_sample -=
      dataSetup.popTree->migBands[p_stack->live_mig_bands[i]].migRate;
  }
  if( p_stack->event_sample >= 0.0 )
  { // at this point i = num_live_mig_bands
    TRACE_LINEAGE_FATAL_0009
  }

  if(i <= 0)
  {
    TRACE_LINEAGE_FATAL_0009A
  }
  // create in and out events and follow out event
  // do not mark the events yet (unless accepted)
  locus_data[gen].mig_spr_stats.new_migs_bands[locus_data[gen].\
    mig_spr_stats.num_new_migs] = p_stack->mig_band
                                = p_stack->live_mig_bands[i-1];
  if(dataSetup.popTree->migBands[p_stack->mig_band].targetPop != p_stack->pop)
  {
      TRACE_LINEAGE_FATAL_0009B
  }

  locus_data[gen].mig_spr_stats.\
    new_migs_ages[locus_data[gen].mig_spr_stats.num_new_migs] = p_stack->age;

  int mig_event_idx = createEventBefore( gen,
                                         p_stack->pop,
                                         p_stack->pEvent,
                                         p_stack->t);
  locus_data[gen]\
    .mig_spr_stats\
      .new_migs_in[locus_data[gen].mig_spr_stats.num_new_migs] = mig_event_idx;
  p_stack->pEvent = &(event_chains[gen].events[mig_event_idx]);

  // mark source event for migration
  p_stack->mig_source =
    locus_data[gen].mig_spr_stats.\
                    new_migs_out[locus_data[gen].mig_spr_stats.num_new_migs]
      = createEvent( gen,
                     dataSetup.popTree->migBands[p_stack->mig_band].sourcePop,
                     p_stack->age);
  if( p_stack->mig_source < 0 )
  {
    TRACE_LINEAGE_FATAL_0010
  }
  locus_data[gen].mig_spr_stats.num_new_migs++;

  return 0;
}

//----------------------------------------------------------------------------
int
traceLineage_sample_coalescence_event_in_interval(TraceLineageAutoVars* p_stack)
{
  int i;
  int& gen  = p_stack->gen;
  int& node = p_stack->node;

  // coalescence event - figure out with whom to coalesce
  //getLineagesAtInterval(gen,event,pop,node,targets);  <-- UNUSED
  double tm =   (p_stack->age - p_stack->t)
           + p_stack->pEvent->getElapsedTime() / 2.0;
  p_stack->num_targets = getEdgesForTimePop( gen, tm, p_stack->pop,
                                             node, p_stack->targets);
  if(     p_stack->num_targets
      !=  p_stack->pEvent->getNumLineages())
  {
    TRACE_LINEAGE_FATAL_0011
  }
  /* printf("\nFound %d targets at gen %d for node %d, at "
            "event %d, pop %d, age %f:",
            num_targets, gen, node, event, pop, age);
     for(i=0; i< num_targets; i++)
     {
       printf(" %d", targets[i]);
     }
     printf("\n");
  */
  // this simulates choosing a lineage uniformly at random
  i = (int) ((p_stack->event_sample - p_stack->mig_rate) * p_stack->theta/2.0);
  p_stack->target = p_stack->targets[i];
  // attach lineage in genetree and create new event
  //printf("\nGrafting subtree below node %d to edge "
  //       "above node %d at time %g.",node, target, age);

  executeGenSPR( dataState.lociData[gen], node,
                 p_stack->target, p_stack->age );
  locus_data[gen].mig_spr_stats.father_pop_new = p_stack->pop;

  //GraftNode(node, target, age, pop);
  //printf("\nOriginal data log-likelihood was %g and after "
  //       "re-grafting it is %g.",data.lnpDi[gen], lnpD_gen(gen));
  locus_data[gen].mig_spr_stats.target = p_stack->target;
  int new_event_idx = createEventBefore( gen,
                                         p_stack->pop,
                                         p_stack->pEvent,
                                         p_stack->t);
  locus_data[gen].mig_spr_stats.father_event_new = new_event_idx;
  p_stack->pEvent = &(event_chains[gen].events[new_event_idx]);

  // signal to terminate
  p_stack->proceed = 0;

  return 0;
}

//----------------------------------------------------------------------------
// If reconnet is true then sample new event within interval
// (coal or mig). If coal is true then create new event, record point
// of coalescence and terminate process. If mig is true then create 2
// new events, and record source of migration for later processing.
int
traceLineage_sample_event_in_interval(TraceLineageAutoVars* p_stack)
{
  int& gen = p_stack->gen;
  int num_lngs = p_stack->pEvent->getNumLineages ();
  p_stack->rate = p_stack->mig_rate + 2.0 * num_lngs / p_stack->theta;
  // rate can be zero, if no lineages and no incoming migration bands
  double curr_elapsed_time = p_stack->pEvent->getElapsedTime ();
  if (p_stack->rate <= 0)
    p_stack->t = curr_elapsed_time;
  else
    p_stack->t = rndexp (gen, 1.0 / p_stack->rate);

  if (curr_elapsed_time <= p_stack->t)
  {
    // no event sampled in this interval
    p_stack->t = curr_elapsed_time;
    p_stack->age += p_stack->t;
  }
  else
  {
    // sample event in this interval
    p_stack->age += p_stack->t;
    p_stack->event_sample = p_stack->rate * rndu(gen);

    if(p_stack->event_sample < p_stack->mig_rate)
    {
      if (0 != traceLineage_sample_migration_event_in_interval(p_stack))
        return -1;
    }
    else
    {
      if (0 != traceLineage_sample_coalescence_event_in_interval(p_stack))
        return -1;
    }
  }      // if event in interval
  return 0;
}
//----------------------------------------------------------------------------
// if event is migration (original or new)
// follow to source and recalculate living mig_bands
// if event starts/ends migration band, change mig_rate
void
traceLineage_epilogue_recalc_mig_bands(TraceLineageAutoVars* p_stack)
{
  int& gen       = p_stack->gen;
  int& reconnect = p_stack->reconnect;

  // add log-likelihood of migration event
  locus_data[gen].mig_spr_stats.genetree_delta_lnLd[reconnect] +=
                 log(dataSetup.popTree->migBands[p_stack->mig_band].migRate);
  p_stack->pEvent = &(event_chains[gen].events[p_stack->mig_source]);
  p_stack->pop = dataSetup.popTree->migBands[p_stack->mig_band].sourcePop;
  p_stack->theta =   dataSetup.popTree->pops[p_stack->pop]->theta
                   * p_stack->heredity_factor;
  p_stack->mig_source = -1;
  p_stack->mig_rate = 0.0;
  p_stack->num_live_mig_bands = 0;
  for( p_stack->mig_band = 0;
       p_stack->mig_band < dataSetup.popTree->numMigBands;
       p_stack->mig_band++) {
    if( dataSetup.popTree->migBands[p_stack->mig_band].targetPop ==
                                                             p_stack->pop &&
        dataSetup.popTree->migBands[p_stack->mig_band].startTime <=
                                                             p_stack->age &&
        dataSetup.popTree->migBands[p_stack->mig_band].endTime > p_stack->age)
      {
        p_stack->mig_rate +=
                      dataSetup.popTree->migBands[p_stack->mig_band].migRate;
        p_stack->live_mig_bands[p_stack->num_live_mig_bands++] =
                                                           p_stack->mig_band;
      }
  }
}

//----------------------------------------------------------------------------
void
traceLineage_epilogue_record_mig_band_end(TraceLineageAutoVars* p_stack)
{
  int i = 0;
  int& gen       = p_stack->gen;

  p_stack->mig_rate -= dataSetup.popTree->migBands[p_stack->node_id].migRate;
  if(p_stack->num_live_mig_bands == 1 && p_stack->mig_rate > 0.0000001 && debug)
  {
    fprintf(stderr, "\ntraceLineage() mig_rate=%g when all mig bands "\
                    "have been emptied at pop=%d, age=%g.\n",
                    p_stack->mig_rate, p_stack->pop, p_stack->age);
    printPopulationTree(dataSetup.popTree, stderr, 1);
    printLocusGenTree(dataState.lociData[gen], stderr,
                      nodePops[gen], nodeEvents[gen]);
    printEventChains(stderr, gen);
    fprintf(stderr, "\n\n************************************"
                    "********************************\n\n");
  }
  if(p_stack->num_live_mig_bands == 1)
  {
    p_stack->mig_rate = 0.0;
  }
  for(i=0; i< p_stack->num_live_mig_bands; i++)
  {
    if( p_stack->live_mig_bands[i] == p_stack->node_id )
    {
      p_stack->live_mig_bands[i] =
                     p_stack->live_mig_bands[--p_stack->num_live_mig_bands];
      break;
    }
  }
}

//----------------------------------------------------------------------------
void
traceLineage_epilogue_record_stat_changes(TraceLineageAutoVars* p_stack)
{
  int i = 0;
  int& gen       = p_stack->gen;
  int& reconnect = p_stack->reconnect;

  locus_data[gen].genetree_stats_delta[reconnect].\
    coal_stats_delta[p_stack->pop] +=   2.0
                                      * p_stack->pEvent->getNumLineages()
                                      * p_stack->t;
  for( i = 0; i < p_stack->num_live_mig_bands; i++)
  {
    locus_data[gen].genetree_stats_delta[reconnect].\
      mig_stats_delta[p_stack->live_mig_bands[i]] += p_stack->t;
  }
  locus_data[gen]\
    .genetree_stats_delta[reconnect]\
      .push_changed_event(p_stack->pEvent);

  // add log-likelihood of no events during interval
  locus_data[gen].mig_spr_stats.genetree_delta_lnLd[reconnect] -=
    (   p_stack->mig_rate
      +    2.0
         * p_stack->pEvent->getNumLineages()
         / p_stack->theta)
    * p_stack->t;

  // if event is migration (original or new)
  // follow to source and recalculate living mig_bands
  // if event starts/ends migration band, change mig_rate
  if( p_stack->mig_source >= 0 )
  {
    traceLineage_epilogue_recalc_mig_bands(p_stack);
  }
  else if(p_stack->pEvent->getType()== MIG_BAND_START)
  {
    p_stack->mig_rate += dataSetup.popTree->migBands[p_stack->node_id].migRate;
    p_stack->live_mig_bands[p_stack->num_live_mig_bands] = p_stack->node_id;
    p_stack->num_live_mig_bands++;
  }
  else if(p_stack->pEvent->getType() == MIG_BAND_END)
  {
    traceLineage_epilogue_record_mig_band_end(p_stack);
  }
}

//----------------------------------------------------------------------------
/*  traceLineage
  Traces a lineage up a genetree.
  It reconnect == 0, then traces original edge, reduces number
  of lineages in events along the way by 1, and records the stat changes
  in genetree_stats_delta[0].
  If reconnect == 1, then re-coalesces lineage with rest of genetree according
  to coalescence with migration process determined by model parameters.
  records stat changes in genetree_stats_delta[0 ??1, maybe?? ].

  In both cases, sets the appropriate entries in mig_spr_stats.

  Returns -1 if too many migration events were sampled (otherwise 0).
*/

int traceLineage (int gen, int node, int reconnect)
{
  TraceLineageAutoVars stack_vars;
  stack_vars.gen = gen;
  stack_vars.node = node;
  stack_vars.reconnect = reconnect;

  traceLineage_init (&stack_vars);

  // loop will stop when reaching original father's event (if !reconnect)
  // or when creating new coalescent event for lineage (if reconnect)
  while (stack_vars.proceed)
  {
    // if last event in population, move to next pop
    if( nullptr == stack_vars.pEvent )//|| stack_vars.pEvent->isOfType(END_CHAIN) )
      if (0 != traceLineage_move_to_next_pop (&stack_vars))
        return -1;

    stack_vars.node_id = stack_vars.pEvent->getId ();

    if (!stack_vars.reconnect)
    {
      traceLineage_reduce_lineages_along_pruned_edge (&stack_vars);
    }
    else
    {
      if( 0 != traceLineage_sample_event_in_interval(&stack_vars))
        return -1;
    }      // if (reconnect)

    // record stat changes
    traceLineage_epilogue_record_stat_changes(&stack_vars);

    stack_vars.pEvent = get_next_event( gen, stack_vars.pEvent );
  }      // end of while

  // add log-likelihood of father coalescent
  locus_data[gen].mig_spr_stats.genetree_delta_lnLd[reconnect] +=
                                                  log (2.0 / stack_vars.theta);
  return 0;
}

//============================ END OF FILE ====================================
