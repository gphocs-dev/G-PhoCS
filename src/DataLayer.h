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

//-----------------------------------------------------------------------------
/* event chain
   Each event corresponds to a time band within a population where no events
   (coalescence/migration) take place. An event is attributed with one of 5 types
   corresponding to the event taking place at the end of the interval. Events are
   sorted in a list according to chronology within a population.
   Actual array of events is allocated in getMem()
*/


enum EventType{COAL           = 0x01,
               IN_MIG         = 0x02,
               OUT_MIG        = 0x04,
               MIG_BAND_START = 0x08,
               MIG_BAND_END   = 0x0F,
               SAMPLES_START  = 0x10,
               END_CHAIN      = 0x20,
               DUMMY          = 0x00};
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

class EventChain
{
public:
  int total_events;		  // total number of events pre-allocated
                                  // to this chain

  int first_event[2*NSPECIES-1];  // pointers to first event
                                  // for every population

  int last_event[2*NSPECIES-1];   // pointers to last event
                                  // for every population

  int free_events;                // pointer to a chain of free events for use.
                                  // Always have at least one free event
  Event* events;
};

extern vector<EventChain> event_chains;
//-----------------------------------------------------------------------------
//============================ END OF FILE ====================================
