#include "GTreeStructVal.h"
#include "DataLayer.h"
#include "DbgErrMsgGTreeStruct.h"
#include "LocusDataLikelihood.h"
#include "utils.h"
#include <math.h>

#include "MCMCcontrol.h"
extern DATA_SETUP  dataSetup;

#include "MemoryMng.h"
extern Locus_SuperStruct*   locus_data;

#include "GPhoCS.h"
extern DATA_STATE dataState;
//-----------------------------------------------------------------------------
int checkGtreeStructure_dispatch_check_event_chain(
                                         CheckGtreeStructureAutoVars* p_stack)
{
  int& gen                  = p_stack->gen;
  int& num_lineages         = p_stack->num_lineages;
  int& pop                  = p_stack->pop;
  int& mig_band             = p_stack->mig_band;
  int& event_idx            = p_stack->event_idx;
  int& event_id             = p_stack->event_id;
  int& num_living_mig_bands = p_stack->num_living_mig_bands;
  int& res                  = p_stack->res;
  Event* pCurrEvent         = p_stack->pCurrEvent;
  double& age               = p_stack->age;
  double PRECISION          = p_stack->PRECISION;

  int i = 0;
  double curr_nodeAge = -1.;
  switch(pCurrEvent->getType())
  {
    case(SAMPLES_START):
      num_lineages += dataSetup.numSamplesPerPop[pop];
      if( fabs( dataSetup.popTree->pops[pop]->sampleAge - age) > PRECISION )
      {
        CHECK_GTREE_STRUCT_FATAL_0036
      }
    break;
    case(COAL):
      ++(locus_data[gen].genetree_stats_check.num_coals[pop]);
      --num_lineages;
      curr_nodeAge = getNodeAge( dataState.lociData[gen], event_id );
      if( fabs( curr_nodeAge - age ) > PRECISION )
      {
        CHECK_GTREE_STRUCT_FATAL_0036_2
       }
      if(nodePops[gen][event_id] != pop)
      {
        CHECK_GTREE_STRUCT_FATAL_0037
      }
      if(nodeEvents[gen][event_id] != event_idx)
      {
        CHECK_GTREE_STRUCT_FATAL_0038
      }
    break;
    case IN_MIG:
      // figure out migration band and update its statistics
      mig_band = genetree_migs[gen].mignodes[event_id].migration_band;
      if( mig_band < 0 || mig_band > dataSetup.popTree->numMigBands )
      {
        CHECK_GTREE_STRUCT_FATAL_0039a
      }
      locus_data[gen].genetree_stats_check.num_migs[mig_band]++;
      --num_lineages;
      if( fabs(genetree_migs[gen].mignodes[event_id].age - age) > PRECISION )
      {
        if(debug)
        {
          fprintf(stderr, "\nError: checking genetree for gen %d: ",gen);
          fprintf(stderr, "age of mignode %d (target event_idx %d)"
                          " doesn't match: %g, %g.", event_id, event_idx, age,
                          genetree_migs[gen].mignodes[event_id].age);
        }
        else
        {
          fprintf(stderr, "Fatal Error 0039.\n");
        }
        res = 0;
      }
      if(dataSetup.popTree->migBands[mig_band].targetPop != pop)
      {
        if(debug)
        {
          fprintf(stderr, "\nError: checking genetree for gen %d: ",gen);
          fprintf(stderr, "target population of migration band %d "
                          "of mignode %d (target event_idx %d) "
                          "doesn't match: %d, %d.",
                          mig_band, event_id, event_idx, pop,
                          dataSetup.popTree->migBands[mig_band].targetPop);
        }
        else
        {
          fprintf(stderr, "Fatal Error 0040.\n");
        }
        res = 0;
      }
    break;
    case(OUT_MIG):
      ++num_lineages;
      mig_band = genetree_migs[gen].mignodes[event_id].migration_band;
      if(fabs(genetree_migs[gen].mignodes[event_id].age - age) > PRECISION)
      {
        if(debug)
        {
          fprintf(stderr, "\nError: checking genetree for gen %d: ",gen);
          fprintf(stderr, "age of mignode %d (source event_idx %d) doesn't"
                          " match: %g, %g.", event_id, event_idx, age,
                          genetree_migs[gen].mignodes[event_id].age);
        }
        else
        {
          fprintf(stderr, "Fatal Error 0041.\n");
        }
        res = 0;
      }
      if( dataSetup.popTree->migBands[mig_band].sourcePop != pop )
      {
        if(debug)
        {
          fprintf(stderr, "\nError: checking genetree for gen %d: ",gen);
          fprintf(stderr, "source population of migration band %d of mignode "
                          "%d (source event_idx %d) doesn't match: %d, %d.",
                          mig_band, event_id, event_idx, pop,
                          dataSetup.popTree->migBands[mig_band].sourcePop);
        }
        else
        {
          fprintf(stderr, "Fatal Error 0042.\n");
        }
        res = 0;
      }
      if(genetree_migs[gen].mignodes[event_id].source_event != event_idx)
      {
        if(debug)
        {
          fprintf(stderr, "\nError: checking genetree for gen %d: ",gen);
          fprintf(stderr, "source event_idx of node %d doesn't match: %d, %d.",
                          event_id, event_idx,
                          genetree_migs[gen].mignodes[event_id].source_event);
        }
        else
        {
          fprintf(stderr, "Fatal Error 0043.\n");
        }
        res = 0;
      }
    break;
    case(MIG_BAND_START):
      p_stack->living_mig_bands[num_living_mig_bands] = event_id;
      ++num_living_mig_bands;
      // initialize statistics for this new migration band
      locus_data[gen].genetree_stats_check.num_migs[event_id] = 0;
      locus_data[gen].genetree_stats_check.mig_stats[event_id] = 0.0;
      if(fabs(dataSetup.popTree->migBands[event_id].startTime - age) > PRECISION)
      {
        if(debug)
        {
          fprintf(stderr, "\nError: checking genetree for gen %d: ",gen);
          fprintf(stderr, "start time of migration band %d (start event_idx %d)"
                          " doesn't match: %g, %g.", event_id, event_idx, age,
                          dataSetup.popTree->migBands[event_id].startTime);
        }
        else
        {
          fprintf(stderr, "Fatal Error 0044.\n");
        }
        res = 0;
      }
      if(dataSetup.popTree->migBands[event_id].targetPop != pop)
      {
        if(debug)
        {
          fprintf(stderr, "\nError: checking genetree for gen %d: ",gen);
          fprintf(stderr, "target population of migration band %d (start "
                          "event_idx %d) doesn't match: %d, %d.",
                          event_id, event_idx, pop,
                          dataSetup.popTree->migBands[event_id].targetPop);
        }
        else
        {
          fprintf(stderr, "Fatal Error 0045.\n");
        }
        res = 0;
      }
      if( fabs(age - max2(dataSetup.popTree->pops[pop]->age,
                          dataSetup.popTree->pops[\
                            dataSetup.popTree->\
                              migBands[event_id].sourcePop]->age)) > PRECISION)
      {
        if(debug)
        {
          fprintf(stderr, "\nError: checking genetree for gen %d: ",gen);
          fprintf(stderr, "start time of migration band %d (start event %d) "
                          "doesn't match max of start times for populations "
                          "%d,%d: %g, max(%g,%g).",
                          event_id, event_idx, pop,
                          dataSetup.popTree->migBands[event_id].sourcePop,
                          age, dataSetup.popTree->pops[pop]->age,
                          dataSetup.popTree->pops[\
                            dataSetup.popTree->\
                              migBands[event_id].sourcePop]->age);
        }
        else
        {
          fprintf(stderr, "Fatal Error 0046.\n");
        }
        res = 0;
      }
    break;
    case(MIG_BAND_END):
      for(mig_band=0; mig_band<num_living_mig_bands; mig_band++)
      {
        if(p_stack->living_mig_bands[mig_band] == event_id)
          break;
      }
      if( mig_band == num_living_mig_bands )
      {
        if( debug )
        {
          fprintf(stderr, "\nError: checking genetree for gen %d: ",gen);
          fprintf(stderr, "migration band %d (start event_idx %d)"
                          " ended without starting.", event_id, event_idx);
        }
        else
        {
          fprintf(stderr, "Fatal Error 0047.\n");
        }
        res = 0;
      }
      else
      {
        --num_living_mig_bands;
        p_stack->living_mig_bands[mig_band] =
                               p_stack->living_mig_bands[num_living_mig_bands];
      }
      if(fabs(dataSetup.popTree->migBands[event_id].endTime - age) > PRECISION)
      {
        if(debug)
        {
          fprintf(stderr, "\nError: checking genetree for gen %d: ",gen);
          fprintf(stderr, "end time of migration band %d (end event_idx %d)"
                          " doesn't match: %g, %g.", event_id, event_idx, age,
                          dataSetup.popTree->migBands[event_id].endTime);
        }
        else
        {
          fprintf(stderr, "Fatal Error 0048.\n");
        }
        res = 0;
      }
      if(dataSetup.popTree->migBands[event_id].targetPop != pop)
      {
        if(debug)
        {
          fprintf(stderr, "\nError: checking genetree for gen %d: ",gen);
          fprintf(stderr, "target population of migration band %d "
                          "(end event_idx %d) doesn't match: %d, %d.",
                          event_id, event_idx, pop,
                          dataSetup.popTree->migBands[event_id].targetPop);
        }
        else
        {
          fprintf(stderr, "Fatal Error 0049.\n");
        }
        res = 0;
      }
      if(fabs( age - min2(  dataSetup.popTree->pops[pop]->father->age,
                            dataSetup.popTree->pops[\
                              dataSetup.popTree->migBands[event_id].\
                                sourcePop]->father->age) ) > PRECISION)
      {
        if(debug)
        {
          fprintf(stderr, "\nError: checking genetree for gen %d: ",gen);
          fprintf(stderr, "end time of migration band %d (end event %d) "
                          "doesn't match min of start times for parent "
                          "populations %d,%d: %g, min(%g,%g).",
                          event_id, event_idx,
                          dataSetup.popTree->pops[pop]->father->id,
                          dataSetup.popTree->pops[\
                            dataSetup.popTree->migBands[\
                              event_id].sourcePop]->father->id,
                          age, dataSetup.popTree->pops[pop]->father->age,
                          dataSetup.popTree->pops[\
                            dataSetup.popTree->migBands[\
                              event_id].sourcePop]->father->age);
        }
        else
        {
          fprintf(stderr, "Fatal Error 0050.\n");
        }
        res = 0;
      }
    break;
   case(END_CHAIN):
    if(event_id != pop)
    {
      if(debug)
      {
        fprintf(stderr, "\nError: checking genetree for gen %d: ",gen);
        fprintf(stderr, "event_id of end event_idx %d for population %d "
                        "doesn't match: %d.", event_idx, pop, event_id);
      }
      else
      {
        fprintf(stderr, "Fatal Error 0051.\n");
      }
      res = 0;
    }
    if(num_living_mig_bands != 0) {
  if(debug) {
    fprintf(stderr, "\nError: checking genetree for gen %d: ",gen);
    fprintf(stderr, "there are %d living migration bands remaining in population %d:", num_living_mig_bands,pop);
    for(i=0; i<num_living_mig_bands; i++) {
      fprintf(stderr," %d",p_stack->living_mig_bands[i]);
    }
  } else {
    fprintf(stderr, "Fatal Error 0052.\n");
  }
      res = 0;
    }
    if(event_chains[gen].events[event_idx].getNextIdx() >= 0) {
  if(debug) {
    fprintf(stderr, "\nError: checking genetree for gen %d: ",gen);
    fprintf(stderr, "remaining event_idx %d after end event_idx %d for pop %d.", event_chains[gen].events[event_idx].getNextIdx(), event_idx, pop);
  } else {
    fprintf(stderr, "Fatal Error 0053.\n");
  }
      res = 0;
    }

    //          if(pop == dataSetup.popTree->rootPop) {
    //            if(fabs(OLDAGE - age) > PERCISION) {
    //              printf("\nError checking genetree for gen %d: ",gen);
    //              printf("end time of population %d (end event %d) doesn't match OLDAGE: %f, %f.", pop, event_idx, age, OLDAGE);
    //              res = 0;
    //            }
    //          }
    if(pop != dataSetup.popTree->rootPop) {
      // add to incoming lineages of parent population
      p_stack->pop_lins_in[dataSetup.popTree->pops[pop]->father->id] +=
                                                                  num_lineages;
      if(fabs(dataSetup.popTree->pops[pop]->father->age - age) > PRECISION) {
  if(debug) {
    fprintf(stderr, "\nError: checking genetree for gen %d: ",gen);
    fprintf(stderr, "end time of population %d (end event_idx %d) doesn't match start of father pop %d: %f, %f.",
        pop, event_idx, dataSetup.popTree->pops[pop]->father->id, age, dataSetup.popTree->pops[pop]->father->age);
  } else {
    fprintf(stderr, "Fatal Error 0054.\n");
  }
        res = 0;
      }
    }
    break;
  default:
  if(debug) {
    fprintf(stderr, "\nError: checking genetree for gen %d: ",gen);
    fprintf(stderr, "Event %d of illegal type: %d.", event_idx, event_chains[gen].events[event_idx].getType());
  } else {
    fprintf(stderr, "Fatal Error 0055.\n");
  }
    res = 0;
    break;
  }// end of switch
  return res;
}

/*-----------------------------------------------------------------------------
 * checkGtreeStructure
 * Checks consistency of genetree data - nodes, mignodes, events and stats
 ----------------------------------------------------------------------------*/
int checkGtreeStructure(int gen)
{
  int i = 0;
  CheckGtreeStructureAutoVars stack_vars;
  stack_vars.gen = gen;
  int& pop                  = stack_vars.pop;
  int& num_lineages         = stack_vars.num_lineages;
  int& event_idx            = stack_vars.event_idx;
  int& event_id             = stack_vars.event_id;
  int& mig_band             = stack_vars.mig_band;
  int& num_living_mig_bands = stack_vars.num_living_mig_bands;
  int& node                 = stack_vars.node;
  int& mig                  = stack_vars.mig;
  int& num_migs             = stack_vars.num_migs;
  int& father               = stack_vars.father;
  int& res                  = stack_vars.res;
  double& age               = stack_vars.age;
  double& delta_t           = stack_vars.delta_t;
  double PRECISION          = stack_vars.PRECISION;
  Event* pCurrEvent         = stack_vars.pCurrEvent;

  // essentially follow the same path as procedure computeGenetreeStats,
  // but validates with genetree nodes

  // initialize number of in-lineages for leaf pops and order pops
  for(pop=0; pop<dataSetup.popTree->numPops; pop++)
  {
    stack_vars.pop_lins_in[pop] = 0;
  }
  // MARK CHANGHES

  // for(node=0; node<2*dataSetup.numSamples-1; node++) {
  //  if(node < dataSetup.numSamples)
  //     pop_lins_in[ nodePops[gen][node] ]++;
  // }
  populationPostOrder(dataSetup.popTree->rootPop, stack_vars.pop_queue);


  // --- Check event chain per population
  for( i = 0; i < dataSetup.popTree->numPops; ++i )
  {
    pop = stack_vars.pop_queue[i];
    locus_data[gen].genetree_stats_check.coal_stats[pop] = 0.0;
    locus_data[gen].genetree_stats_check.num_coals[pop]  = 0;
    event_idx            = event_chains[gen].first_event[pop];
    pCurrEvent           = event_chains.getEvent(gen, stack_vars.event_idx);
    num_lineages         = stack_vars.pop_lins_in[pop];
    age                  = dataSetup.popTree->pops[pop]->age;
    num_living_mig_bands = 0;


    // follow event_idx chain and check all saved data
    for( ; stack_vars.event_idx >=0;
         stack_vars.event_idx = stack_vars.pCurrEvent->getNextIdx() )
    {
      stack_vars.pCurrEvent = event_chains.getEvent(gen, stack_vars.event_idx);
      if( stack_vars.pCurrEvent->getNumLineages() != num_lineages )
      {
        CHECK_GTREE_STRUCT_FATAL_0034
      }

      int event_next_idx = stack_vars.pCurrEvent->getNextIdx();
      if( event_next_idx >=0 )
      {
        int event_next_prev_idx =
                event_chains.getEvent(gen, event_next_idx)->getPrevIdx();
        if(stack_vars.event_idx != event_next_prev_idx)
          CHECK_GTREE_STRUCT_FATAL_0035
      }
      event_id  = pCurrEvent->getId();
      delta_t   = pCurrEvent->getElapsedTime();
      age      += delta_t;

      locus_data[gen].genetree_stats_check.coal_stats[pop] +=
              num_lineages *( num_lineages - 1) * delta_t;

      for( mig_band = 0; mig_band < num_living_mig_bands; ++mig_band )
      {
        locus_data[gen].genetree_stats_check.\
          mig_stats[stack_vars.living_mig_bands[mig_band]] +=
                                                        num_lineages * delta_t;
      }
      checkGtreeStructure_dispatch_check_event_chain(&stack_vars);
    }// end of for(event_idx)
        
    if(stack_vars.num_living_mig_bands != 0)
    {
      CHECK_GTREE_STRUCT_FATAL_0056
    }
  }// end of for(pop)
  // --- End check event chain per population
  
  
  // --- Check computed stats per population ---
  for( pop = 0; pop < dataSetup.popTree->numPops; ++pop )
  {
    if( fabs(locus_data[gen].genetree_stats_check.coal_stats[pop] -
             genetree_stats[gen].coal_stats[pop]) > PRECISION )
    {
      CHECK_GTREE_STRUCT_FATAL_0057
      res = 0;
    }
    genetree_stats[gen].coal_stats[pop] =
                        locus_data[gen].genetree_stats_check.coal_stats[pop];
    if( locus_data[gen].genetree_stats_check.num_coals[pop] !=
          genetree_stats[gen].num_coals[pop] )
    {
      CHECK_GTREE_STRUCT_FATAL_0058
      res = 0;
    }
  }

  // --- Check migration bands per population ---
  for(mig_band=0; mig_band<dataSetup.popTree->numMigBands; mig_band++)
  {
    if(fabs( locus_data[gen].genetree_stats_check.mig_stats[mig_band] -
             genetree_stats[gen].mig_stats[mig_band]) > PRECISION)
    {
      if(debug)
      {
        fprintf(stderr, "\nError: checking genetree for gen %d: ",gen);
        fprintf(stderr, "mig stats for migration band %d don't match: %f, %f.",
                mig_band,
                genetree_stats[gen].mig_stats[mig_band],
                locus_data[gen].genetree_stats_check.mig_stats[mig_band]);
      }
      else
      {
        fprintf(stderr, "Fatal Error 0059.\n");
      }
      res = 0;
    }
    genetree_stats[gen].mig_stats[mig_band] =
            locus_data[gen].genetree_stats_check.mig_stats[mig_band];
    if( locus_data[gen].genetree_stats_check.num_migs[mig_band] !=
          genetree_stats[gen].num_migs[mig_band] )
    {
      if(debug)
      {
        fprintf(stderr, "\nError: checking genetree for gen %d: ",gen);
        fprintf(stderr, "num migs for migration band %d don't match: %d, %d.",
                mig_band, genetree_stats[gen].num_migs[mig_band],
                locus_data[gen].genetree_stats_check.num_migs[mig_band]);
      }
      else
      {
        fprintf(stderr, "Fatal Error 0060.\n");
      }
      res = 0;
    }
  }

  // --- check tree structure ---
  num_migs = 0;
  for(node=0; node<2*dataSetup.numSamples-1; node++)
  {
    pop = nodePops[gen][node];
    mig = findFirstMig(gen,node,-1);
    age = getNodeAge(dataState.lociData[gen],node);
    //age = gnodes[gen][node].age;
    while(mig >= 0)
    {
      ++num_migs;
      if( genetree_migs[gen].mignodes[mig].gtree_branch != node )
      {
        if(debug)
        {
          fprintf(stderr, "\nError: checking genetree for gen %d: ",gen);
          fprintf(stderr, "migration node %d for node %d points "
                          "to wrong edge: %d.",
                  mig, node, genetree_migs[gen].mignodes[mig].gtree_branch);
        }
        else
        {
          fprintf(stderr, "Fatal Error 0061.\n");
        }
        res = 0;
      }
      if(genetree_migs[gen].mignodes[mig].age <= age)
      {
        if(debug)
        {
          fprintf(stderr, "\nError: checking genetree for gen %d: ",gen);
          fprintf(stderr, "migration node %d for node %d has conflicting "
                          "age: %f, %f.", mig, node,
                          genetree_migs[gen].mignodes[mig].age, age);
        }
        else
        {
          fprintf(stderr, "Fatal Error 0062.\n");
        }
        res = 0;
      }
      for( i=0; i < genetree_migs[gen].num_migs; ++i )
      {
        if(genetree_migs[gen].living_mignodes[i] == mig)
          break;
      }
      if( i >= genetree_migs[gen].num_migs)
      {
        if(debug)
        {
          fprintf(stderr, "\nError: checking genetree for gen %d: ",gen);
          fprintf(stderr, "mignode %d doesn't appear in list of "
                          "living mignodes.", mig);
        }
        else
        {
          fprintf(stderr, "Fatal Error 0063.\n");
        }
        res = 0;
      }
      age = genetree_migs[gen].mignodes[mig].age;
      mig_band = genetree_migs[gen].mignodes[mig].migration_band;
      
      if( dataSetup.popTree->migBands[mig_band].sourcePop !=
          genetree_migs[gen].mignodes[mig].source_pop)
      {
        if(debug)
        {
          fprintf(stderr, "\nError: checking genetree for gen %d: ",gen);
          fprintf(stderr, "source population of mignode %d doesn't"
                          " match mig band %d: %d, %d.",
                  mig, mig_band,
                  genetree_migs[gen].mignodes[mig].source_pop,
                  dataSetup.popTree->migBands[mig_band].sourcePop);
        }
        else
        {
          fprintf(stderr, "Fatal Error 0064.\n");
        }
        res = 0;
      }
      if( dataSetup.popTree->migBands[mig_band].targetPop !=
              genetree_migs[gen].mignodes[mig].target_pop)
      {
        if(debug)
        {
          fprintf(stderr, "\nError: checking genetree for gen %d: ",gen);
          fprintf(stderr, "target population of mignode %d doesn't "
                          "match mig band %d: %d, %d.",
                  mig, mig_band,
                  genetree_migs[gen].mignodes[mig].target_pop,
                  dataSetup.popTree->migBands[mig_band].targetPop);
        }
        else
        {
          fprintf(stderr, "Fatal Error 0065.\n");
        }
        res = 0;
      }
      if(    age <= dataSetup.popTree->migBands[mig_band].startTime
          || age >= dataSetup.popTree->migBands[mig_band].endTime )
      {
        if(debug)
        {
          fprintf(stderr, "\nError: checking genetree for gen %d: ",gen);
          fprintf(stderr, "age of mignode %d conflicts with times of "
                          "migration band %d: %f, [%f,%f].",
              mig, mig_band, age,
              dataSetup.popTree->migBands[mig_band].startTime,
              dataSetup.popTree->migBands[mig_band].endTime);
        }
        else
        {
          fprintf(stderr, "Fatal Error 0066.\n");
        }
        res = 0;
      }
      if( !dataSetup.popTree->pops[ genetree_migs[gen].\
           mignodes[mig].target_pop ]->isAncestralTo[pop] )
      {
        if(debug)
        {
          fprintf(stderr, "\nError: checking genetree for gen %d: ",gen);
          fprintf(stderr, "target pop %d of mignode %d is not ancestral "
                          "to previous pop %d.",
                  genetree_migs[gen].mignodes[mig].target_pop,
                  mig, pop);
        }
        else
        {
          fprintf(stderr, "Fatal Error 0067.\n");
        }
        res = 0;
      }
      
      pop = genetree_migs[gen].mignodes[mig].source_pop;
      mig = findFirstMig(gen,node,age);
    }// end of while(mig>=0)

    father = getNodeFather(dataState.lociData[gen], node);
    if(father >= 0) {
      //    if(gnodes[gen][node].father >= 0) {
      /*      if(gnodes[gen][gnodes[gen][node].father].nson != 2) {
                    printf("\nError checking genetree for gen %d: ",gen);
                    printf("father %d of node %d has number of sons: %d.", gnodes[gen][node].father, node, gnodes[gen][gnodes[gen][node].father].nson);
                    res = 0;
                    } else
      */        
      if(getNodeSon(dataState.lociData[gen], father, 0) != node && getNodeSon(dataState.lociData[gen], father, 1) != node) {
        //      if(gnodes[gen][gnodes[gen][node].father].sons[0] != node  && gnodes[gen][gnodes[gen][node].father].sons[1] != node) {
    if(debug) {
      fprintf(stderr, "\nError: checking genetree for gen %d: ",gen);
      fprintf(stderr, "father %d of node %d has no matching son: %d, %d.", father, node, 
          getNodeSon(dataState.lociData[gen], father, 0), getNodeSon(dataState.lociData[gen], father, 1));
    } else {
      fprintf(stderr, "Fatal Error 0068.\n");
    }
        res = 0;
      }
      if(getNodeAge(dataState.lociData[gen], father) <= age) {      
        //      if(gnodes[gen][gnodes[gen][node].father].age <= age) {
    if(debug) {
      fprintf(stderr, "\nError: checking genetree for gen %d: ",gen);
      fprintf(stderr, "age of father %d is smaller than previous: %g, %g (diff %g).", father, getNodeAge(dataState.lociData[gen], father), age,getNodeAge(dataState.lociData[gen], father)- age);
    } else {
      fprintf(stderr, "Fatal Error 0069.\n");
    }
        res = 0;
      }
      if(!dataSetup.popTree->pops[ nodePops[gen][father] ]->isAncestralTo[pop]) {
        //      if(!dataSetup.popTree->pops[ gnodes[gen][gnodes[gen][node].father].ipop ]->isAncestralTo[pop]) {
    if(debug) {
      fprintf(stderr, "\nError: checking genetree for gen %d: ",gen);
      fprintf(stderr, "pop %d of father %d is not ancestral to previous pop %d.", nodePops[gen][father], father, pop);
    } else {
      fprintf(stderr, "Fatal Error 0070.\n");
    }
        res = 0;
      }
    } else if(node != getLocusRoot(dataState.lociData[gen])) {
    if(debug) {
      fprintf(stderr, "\nError: checking genetree for gen %d: ",gen);
      fprintf(stderr, "node %d has no father but is not root %d.", node, getLocusRoot(dataState.lociData[gen]));
    } else {
      fprintf(stderr, "Fatal Error 0071.\n");
    }
      res = 0;
    }
    
         
  }// end of for(node)

  //--- Epilogue ---
  if( num_migs != genetree_migs[gen].num_migs )
  {
    if(debug)
    {
      fprintf(stderr, "\nError: checking genetree for gen %d: ",gen);
      fprintf(stderr, "number of recorded migration nodes doesn't match:"
                      " %d, %d.",
              genetree_migs[gen].num_migs, num_migs);
    }
    else
    {
      fprintf(stderr, "Fatal Error 0072.\n");
    }
    res = 0;
  }
  for( mig_band=0; mig_band < dataSetup.popTree->numMigBands; ++mig_band )
  {
    num_migs -= genetree_stats[gen].num_migs[mig_band];
  }
  if(num_migs != 0)
  {
    if(debug)
    {
      fprintf(stderr, "\nError checking genetree for gen %d: ",gen);
      fprintf(stderr, "number of recorded migration nodes doesn't match sum "
                      "of recorded per mig band: %d, %d.",
              genetree_migs[gen].num_migs,
              genetree_migs[gen].num_migs-num_migs);
    }
    else
    {
      fprintf(stderr, "Fatal Error 0073.\n");
    }
    res = 0;
  }
  
  
  return res;
}
//-----------------------------------------------------------------------------
/* checkAll()
   invokes checkGtreeStructure on all loci and checks genetree_stats_total
   also checks data log likelihood
*/
int checkAll() {
  int gen, mig_band, pop, heredity_factor, res = 1;
  double PERCISION = 0.0000001;
  double lnLd_gen, genLnLd, dataLnLd;

  for(pop=0; pop<dataSetup.popTree->numPops; pop++) {
    genetree_stats_total_check.num_coals[pop]  = 0;
    genetree_stats_total_check.coal_stats[pop] = 0.0;
  }
  for(mig_band=0; mig_band<dataSetup.popTree->numMigBands; mig_band++) {
    genetree_stats_total_check.num_migs[mig_band]  = 0;
    genetree_stats_total_check.mig_stats[mig_band] = 0.0;
  }

  genLnLd = 0.0;
  dataLnLd = 0.0;

  for(gen=0; gen<dataSetup.numLoci; gen++) {
    if(!checkGtreeStructure(gen)) {
      fprintf(stderr, "\nError: Checking gene tree structure failed\n");
      printGenealogyAndExit(gen,0);
      return 0;
    }


    if(!checkLocusDataLikelihood(dataState.lociData[gen])) {
      fprintf(stderr, "\nError: checking recorded likelihood for gen %d!", gen);
      printGenealogyAndExit(gen,0);
      return 0;
    }

    lnLd_gen = gtreeLnLikelihood(gen);

    if(   fabs(locus_data[gen].genLogLikelihood - lnLd_gen) > PERCISION
       && fabs(1-locus_data[gen].genLogLikelihood/lnLd_gen) > PERCISION)
    {
      fprintf(stderr, "\nError: in recorded genealogy log likelihood for gen %d: (recorded = %g, actual = %g, diff = %g, 1-ratio = %g).\n",
      gen, locus_data[gen].genLogLikelihood, lnLd_gen, locus_data[gen].genLogLikelihood - lnLd_gen, 1-locus_data[gen].genLogLikelihood/lnLd_gen);
      printGenealogyAndExit(gen,0);
      return 0;
    }
    locus_data[gen].genLogLikelihood = lnLd_gen;
    dataLnLd += getLocusDataLikelihood(dataState.lociData[gen]);
    genLnLd  += lnLd_gen;

    heredity_factor = 1;
    for(pop=0; pop<dataSetup.popTree->numPops; pop++) {
      genetree_stats_total_check.num_coals[pop]  += genetree_stats[gen].num_coals[pop];
      genetree_stats_total_check.coal_stats[pop] += genetree_stats[gen].coal_stats[pop]/heredity_factor;
    }
    for(mig_band=0; mig_band<dataSetup.popTree->numMigBands; mig_band++) {
      genetree_stats_total_check.num_migs[mig_band]  += genetree_stats[gen].num_migs[mig_band];
      genetree_stats_total_check.mig_stats[mig_band] += genetree_stats[gen].mig_stats[mig_band];
    }


  }

  for(pop=0; pop<dataSetup.popTree->numPops; pop++) {
    if(   fabs(genetree_stats_total_check.coal_stats[pop] - genetree_stats_total.coal_stats[pop]) > PERCISION
       && fabs(1-genetree_stats_total_check.coal_stats[pop]/genetree_stats_total.coal_stats[pop]) > PERCISION) {
    if(debug) {
      fprintf(stderr, "\nError: checking total coal stats for pop %d: %f (saved) %f (actual) %g (difference) %g (1-ratio)",
        pop, genetree_stats_total.coal_stats[pop], genetree_stats_total_check.coal_stats[pop],
          genetree_stats_total_check.coal_stats[pop] - genetree_stats_total.coal_stats[pop],
          1 - genetree_stats_total_check.coal_stats[pop] / genetree_stats_total.coal_stats[pop]);
    } else {
      fprintf(stderr, "Fatal Error 0028.\n");
    }
      res = 0;
    }
    genetree_stats_total.coal_stats[pop] = genetree_stats_total_check.coal_stats[pop];
    if(genetree_stats_total_check.num_coals[pop] != genetree_stats_total.num_coals[pop]) {
    if(debug) {
      fprintf(stderr, "\nError: checking total num coals for pop %d: %d (saved) %d (actual)",
        pop, genetree_stats_total.num_coals[pop], genetree_stats_total_check.num_coals[pop]);
    } else {
      fprintf(stderr, "Fatal Error 0029.\n");
    }
      res = 0;
    }

    genetree_stats_total.coal_stats[pop] = genetree_stats_total_check.coal_stats[pop];
  }
  for(mig_band=0; mig_band<dataSetup.popTree->numMigBands; mig_band++) {
    if(    fabs(genetree_stats_total_check.mig_stats[mig_band] - genetree_stats_total.mig_stats[mig_band]) > PERCISION
        && fabs(1-genetree_stats_total_check.mig_stats[mig_band]/genetree_stats_total.mig_stats[mig_band]) > PERCISION) {
    if(debug) {
      fprintf(stderr, "\nError: checking total mig stats for mig band %d: %f (saved) %f (actual), %f (diff) %f (1-ratio)",
        mig_band, genetree_stats_total.mig_stats[mig_band], locus_data[gen].genetree_stats_check.mig_stats[mig_band],
          genetree_stats_total_check.mig_stats[mig_band] - genetree_stats_total.mig_stats[mig_band],
          1-genetree_stats_total_check.mig_stats[mig_band] / genetree_stats_total.mig_stats[mig_band]);
    } else {
      fprintf(stderr, "Fatal Error 0030.\n");
    }
      res = 0;
    }
    genetree_stats_total.mig_stats[mig_band] = genetree_stats_total_check.mig_stats[mig_band];
    if(genetree_stats_total_check.num_migs[mig_band] != genetree_stats_total.num_migs[mig_band]) {
    if(debug) {
      fprintf(stderr, "\nError: checking total num migs for mig band %d: %d (saved) %d (actual)",
        mig_band, genetree_stats_total.num_migs[mig_band], genetree_stats_total_check.num_migs[mig_band]);
    } else {
      fprintf(stderr, "Fatal Error 0031.\n");
    }
      res = 0;
    }

    genetree_stats_total.mig_stats[mig_band] = genetree_stats_total_check.mig_stats[mig_band];
  }
  if(fabs(dataState.dataLogLikelihood - dataLnLd) > PERCISION && fabs(1-dataState.dataLogLikelihood/dataLnLd) > PERCISION) {
    if(debug) {
      fprintf(stderr, "\nError: in recorded data log likelihood: (recorded = %g, actual = %g, diff = %g, 1-ratio = %g).\n",
          dataState.dataLogLikelihood, dataLnLd, dataState.dataLogLikelihood - dataLnLd, 1-dataState.dataLogLikelihood/dataLnLd);
    } else {
      fprintf(stderr, "Fatal Error 0032.\n");
    }
    res = 0;
  }
  dataState.dataLogLikelihood = dataLnLd;
  // check recorded total log-likelihood
  dataLnLd = (genLnLd+dataLnLd)/dataSetup.numLoci;
  if(fabs(1-dataState.logLikelihood/dataLnLd) > PERCISION) {
    if(debug) {
      fprintf(stderr, "\nError: in recorded total log likelihood: (recorded = %g, actual = %g).\n", dataState.logLikelihood, dataLnLd);
    } else {
      fprintf(stderr, "Fatal Error 0033.\n");
    }
    res = 0;
  }
  dataState.logLikelihood = dataLnLd;

  if(!res){
    printf("\n\n");
    //    exit(-1);
  }

  return res;

}
