#pragma once
/*============================================================================
 File: DataLayer.h

 Class declarations of the data layer.

 Version history:
 26-Mar-2017  evgenyidc      Created. Adding Event, EVENT_CHAIN.
 ============================================================================*/
#include <vector>
#include "DataLayerConstants.h"

using namespace std;

/*-----------------------------------------------------------------------------
 * Event
 * Each event corresponds to a time band within a population where no events
 * (coalescence/migration) take place. An event is attributed with one
 * of 5 types corresponding to the event taking place at the end of
 * the interval. Events are sorted in a list according to chronology within
 * a population.
 * Actual array of events is allocated in getMem()
 *---------------------------------------------------------------------------*/

//-----------------------------------------------------------------------------
enum EventType{COAL           = 0x01,
               IN_MIG         = 0x02,
               OUT_MIG        = 0x04,
               MIG_BAND_START = 0x08,
               MIG_BAND_END   = 0x0F,
               SAMPLES_START  = 0x10,
               END_CHAIN      = 0x20,
               POP_START      = 0x40,
               POP_END        = 0x80,
               DUMMY          = 0x00};

//-----------------------------------------------------------------------------
class Event
{
public:
  Event();

public:
  int       getId()                   const;
  void      setId(int id);
  EventType getType()                 const;
  void      setType(EventType t);
  int       getNextIdx()              const;
  void      setNextIdx(int i);
  int       getPrevIdx()              const;
  void      setPrevIdx(int i);
  double    getElapsedTime()          const;
  void      setElapsedTime(double t);
  int       getNumLineages()          const;
  void      setNumLineages(int n);

public:
  void      multiplyElapsedTime(double factor);
  void      addElapsedTime(double delta);
  void      addLineages(int delta);
  int       incrementLineages();
  int       decrementLineages();
  bool      isOfType(int eventTypeMask) const;

protected:
  int        id_;
  EventType  type_;
  int        next_;
  int        prev_;
  double     elapsed_time_;	// time from last event
  int        num_lineages_;	// number of lineages before the event
};

//-----------------------------------------------------------------------------
class EventChain
{
public:
  int total_events;               // total number of events pre-allocated
                                  // to this chain

  int first_event[2*NSPECIES-1];  // pointers to first event
                                  // for every population

  int last_event[2*NSPECIES-1];   // pointers to last event
                                  // for every population

  int free_events;                // pointer to a chain of free events for use.
                                  // Always have at least one free event
  Event* events;
};

//-----------------------------------------------------------------------------

class EventChains : public vector<EventChain>
{
public:
  inline Event* getEvent(int gen, int event) const
    { return &((*this)[gen].events[event]); }
};

template<class T> class DAGsPerLocus;

extern EventChains event_chains;
//extern DAGsPerLocus<Event>* pAllDAGs;


/*-----------------------------------------------------------------------------
 * Setters/Getters section
 * */
inline int
Event::getId() const
{
  return this->id_;
}

inline void
Event::setId(int id)
{
  this->id_ = id;
}

inline EventType
Event::getType() const
{
  return this->type_;
}

inline void
Event::setType(EventType t)
{
  this->type_ = t;
}

inline int
Event::getNextIdx() const
{
  return this->next_;
}

inline void
Event::setNextIdx(int i)
{
  this->next_ = i;
}

inline int
Event::getPrevIdx() const
{
  return this->prev_;
}

inline void
Event::setPrevIdx(int i)
{
  this->prev_ = i;
}

inline double
Event::getElapsedTime() const
{
  return this->elapsed_time_;
}

inline void
Event::setElapsedTime(double t)
{
  this->elapsed_time_ = t;
}

inline int
Event::getNumLineages() const
{
  return this->num_lineages_;
}

inline void
Event::setNumLineages(int n)
{
  this->num_lineages_ = n;
}

/*-----------------------------------------------------------------------------
 * GENETREE_STATS_DELTA
 * Holds difference in stats for populations affected by changes
 * in a genetree. Not generally used for changes in theta, mig_rates.
 * For easy and fast update if changes are accepted.
 * changed_events holds the ids of all events in the interval between
 * original and new position of event. All these intervals have
 * a change in number of lineages.
 *---------------------------------------------------------------------------*/
class GENETREE_STATS_DELTA
{
public:
  void init();

public:
  int num_changed_events() const {return this->changed_events.size();}
  void clear_changed_events() {this->changed_events.clear();}
  void push_changed_event(Event* e) {this->changed_events.push_back(e);}

public:
  // id of event describing original placing of node
  int original_event;
  // id of event describing updated (new) placing of node
  int updated_event;
  // the difference in lineage number for all events affected (typically +/- 1)
  int num_lin_delta;
  // events affected by change
  vector<Event*> changed_events;
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
};

//-----------------------------------------------------------------------------
inline void
GENETREE_STATS_DELTA::init()
{
  this->num_pops_changed      =  0;
  this->num_mig_bands_changed =  0;
  this->num_lin_delta         =  0;
  this->original_event        = -1;
  this->updated_event         = -1;
  changed_events.reserve(MAX_EVENTS);
}
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
class GENETREE_STATS
{
public:
  double coal_stats[2 * NSPECIES - 1];
  double mig_stats[MAX_MIG_BANDS];
  int    num_coals[2 * NSPECIES - 1];
  int    num_migs[MAX_MIG_BANDS];
};

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
class RUBBERBAND_MIGS
{
public:
  int num_moved_events;
  int orig_events[MAX_MIGS + MAX_MIG_BANDS];
  int new_events[MAX_MIGS + MAX_MIG_BANDS];
  int pops[MAX_MIGS + MAX_MIG_BANDS];
  double new_ages[MAX_MIGS + MAX_MIG_BANDS];
};

/*-----------------------------------------------------------------------------
 * MIG_SPR_STATS
 * Holds statistics for the SPR sampling operation with migration.
 * In use in UpdateGB_MigSPR and in traceLineage.
 *---------------------------------------------------------------------------*/
class MIG_SPR_STATS
{
public:
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
};

/*-----------------------------------------------------------------------------
 * Locus_SuperStruct
 *---------------------------------------------------------------------------*/
class Locus_SuperStruct
{
public:
  GENETREE_STATS         genetree_stats_check;
  GENETREE_STATS_DELTA   genetree_stats_delta[NUM_DELTA_STATS_INSTANCES];
  MIG_SPR_STATS          mig_spr_stats;
  double                 genDeltaLogLikelihood;
  double                 genLogLikelihood;
  RUBBERBAND_MIGS        rubberband_migs;
  int                    mig_conflict_log;
};

//============================ END OF FILE ====================================
