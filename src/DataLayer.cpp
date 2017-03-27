/*============================================================================
 File: DataLayer.cpp

 Classes of the Data Layer.

 Version history:
 26-Mar-2017  evgenyidc      Created. Adding Event, EVENT_CHAIN.
 ============================================================================*/

#include "DataLayer.h"

vector<EventChain> event_chains;

//=============================================================================
/*-----------------------------------------------------------------------------
 * CTOR
 * */
Event::Event():
id_(-1),
type_(DUMMY),
next_(-1),
prev_(-1),
elapsed_time_(0.0),
num_lineages_(-1)
{
}

/*-----------------------------------------------------------------------------
 * Setters/Getters section
 * */
EventType
Event::getType() const
{
  return this->type_;
}

void
Event::setType(EventType t)
{
  this->type_ = t;
}

int
Event::getId() const
{
  return this->id_;
}

void
Event::setId(int id)
{
  this->id_ = id;
}

int
Event::getNextIdx() const
{
  return this->next_;
}

void
Event::setNextIdx(int i)
{
  this->next_ = i;
}

int
Event::getPrevIdx() const
{
  return this->prev_;
}

void
Event::setPrevIdx(int i)
{
  this->prev_ = i;
}

double
Event::getElapsedTime() const
{
  return this->elapsed_time_;
}

void
Event::setElapsedTime(double t)
{
  this->elapsed_time_ = t;
}


int
Event::getNumLineages() const
{
  return this->num_lineages_;
}

void
Event::setNumLineages(int n)
{
  this->num_lineages_ = n;
}

/*-----------------------------------------------------------------------------
 * Object state modifiers
 */
void
Event::multiplyElapsedTime(double factor)
{
  this->elapsed_time_ *= factor;
}

void
Event::addElapsedTime(double delta)
{
  this->elapsed_time_ += delta;
}

void
Event::addLineages(int delta)
{
  this->num_lineages_ += delta;
}

int
Event::incrementLineages()
{
  return this->num_lineages_++;
}

int
Event::decrementLineages()
{
  return this->num_lineages_--;
}

/*-----------------------------------------------------------------------------
 * Predicates
 */
bool
Event::isOfType(int eventTypeMask) const
{
  return (this->getType() & eventTypeMask) > 0;
}

//=============================================================================
//=============================================================================

//============================ END OF FILE ====================================
