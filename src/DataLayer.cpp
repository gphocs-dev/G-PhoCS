/*============================================================================
 File: DataLayer.cpp

 Classes of the Data Layer.

 Version history:
 26-Mar-2017  evgenyidc      Created. Adding Event, EVENT_CHAIN.
 ============================================================================*/

#include "PopulationTree.h"
#include "DataLayer.h"
#include "EventsDAG.h"

EventChains  event_chains;
DAGsPerLocus<Event>* pAllDAGs;
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
	bool r = this->type_ & eventTypeMask;
  return r > 0;
}

//=============================================================================
//=============================================================================

//============================ END OF FILE ====================================
