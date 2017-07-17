#pragma once
/*============================================================================
 File: TraceLineages.h

 Trace Lineages functionality

 Version history:
 07-May-2017  evgenyidc    Extracting from patch.cpp.
 ============================================================================*/


int traceLineage (int gen, int node, int reconnect);

//----------------------------------------------------------------------------
typedef struct _TraceLineageAutoVars
{
  int gen, node, reconnect;
  int pop, event, node_id, mig_band, mig_source, proceed;
  int num_live_mig_bands, live_mig_bands[MAX_MIG_BANDS];
  double age, t;
  // these variables are needed only for reconnecting or delta_lnLd computation
  int target, num_targets, targets[2*NS-1];
  double  event_sample, rate, mig_rate, theta, heredity_factor;
} TraceLineageAutoVars;
//----------------------------------------------------------------------------

//============================ END OF FILE ====================================
