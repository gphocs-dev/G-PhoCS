#pragma once

#include "DataLayerConstants.h"

int checkGtreeStructure(int gen);
int checkAll();

class Event;
//-----------------------------------------------------------------------------
typedef struct _CheckGtreeStructureAutoVars
{
  int gen;
  int num_lineages;
  int pop;
  int mig_band;
  int event_idx;
  int event_id;
  int node;
  int mig;
  int num_migs;
  int num_living_mig_bands;
  int father;
  int res = 1;

  Event* pCurrEvent;

  double age;
  double delta_t;
  double PRECISION = 0.0000000001;

  int living_mig_bands[MAX_MIG_BANDS];
  int pop_queue[2*NSPECIES-1];
  int pop_lins_in[2*NSPECIES-1];
} CheckGtreeStructureAutoVars;
//-----------------------------------------------------------------------------
//============================ END OF FILE ====================================
