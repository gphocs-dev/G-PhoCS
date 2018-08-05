#pragma once
/*============================================================================
 File: DataLayer.h

 Class declarations of the data layer.

 Version history:
 26-Mar-2017  evgenyidc      Created. Adding Event, EVENT_CHAIN.
 14-Dec-2017  evgenyidc      EventType values updated. Adding getTypeStr.
 ============================================================================*/
#include <vector>
#include <ostream>
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
               MIG_BAND_END   = 0x10,
               SAMPLES_START  = 0x20,
               END_CHAIN      = 0x40,
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
  string    getTypeStr()              const;
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

private:
  friend ostream& operator<<(ostream& os, const Event& obj);
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

extern vector<EventChain> event_chains;


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

inline string
Event::getTypeStr() const
{
  string strRes = "";
  switch(this->type_)
  {
    case COAL: strRes           = "COAL      "; break;
    case IN_MIG: strRes         = "IN_MIG    "; break;
    case OUT_MIG: strRes        = "OUT_MIG   "; break;
    case MIG_BAND_START: strRes = "MIGB_START"; break;
    case MIG_BAND_END: strRes   = "MIGB_END  "; break;
    case SAMPLES_START: strRes  = "SMPL_START"; break;
    case END_CHAIN: strRes      = "END_CHAIN "; break;
    default: strRes             = "DUMMY     "; break;
  }
  return strRes;
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

//============================ END OF FILE ====================================
