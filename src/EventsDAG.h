#pragma once

/*============================================================================
 File: EventsDag.h

 Directed Acyclic Graph that represents both Population and Genealogy trees.

 Version history:
 29-Jun-2017  evgenyidc      Initial version reworked.
 ============================================================================*/

#include <vector>
#include <cassert>
#include <algorithm>
#include "DataLayer.h"

class Event;
using namespace std;

/*=============================================================================
 *
 * Single Event DAG node declaration
 *
 *---------------------------------------------------------------------------*/

template<class T>
class EventsDAGNode
{
public:
  EventsDAGNode():
      pPrevGenEvent_(nullptr),
      pNextGenEvent_(nullptr),
      pCoalLeftSon_(nullptr),
      pCoalRightSon_(nullptr),
      pCoalParent_(nullptr),
      pContent_(nullptr){}

  explicit EventsDAGNode(T* pContent):
      EventsDAGNode(){pContent_ = pContent;}

  T* getContent() const;
  void setContent( T* pContent);

  //@@TODO: rework the next method. Integer is bad here.
  int getType() const {pContent_->getType();}

  EventsDAGNode<T>* getPrevGenEvent() const;
  void setPrevGenEvent(EventsDAGNode<T>* pPrevGenEvent);

  EventsDAGNode<T>* getNextGenEvent() const;
  void setNextGenEvent(EventsDAGNode<T>* pNextGenEvent );

  EventsDAGNode<T>* getCoalLeftSon() const;
  void setCoalLeftSon(EventsDAGNode<T>* pCoalLeftSon);

  EventsDAGNode<T>* getCoalRightSon() const;
  void setCoalRightSon(EventsDAGNode<T>* pCoalRightSon);

  EventsDAGNode<T>* getCoalParent() const;
  void setCoalParent(EventsDAGNode<T>* pCoalParent);


protected:
  EventsDAGNode<T>* pPrevGenEvent_;
  EventsDAGNode<T>* pNextGenEvent_;
  EventsDAGNode<T>* pCoalLeftSon_;
  EventsDAGNode<T>* pCoalRightSon_;
  EventsDAGNode<T>* pCoalParent_;
  T*                pContent_;
};

/*=============================================================================
 *
 * One DAG, vector of entry points - first Event in every population
 *
 *---------------------------------------------------------------------------*/

template<class T>
class EventsDAG : public vector<EventsDAGNode<T>*>
{
public:
  EventsDAG() = default;
  virtual ~EventsDAG();
  void initDAG( int nNumOfPops );

  EventsDAGNode<T>* getFirstNode(int nPopIdx) const;
  EventsDAGNode<T>* getLastNode(int nPopIdx) const;
  EventsDAGNode<T>* getNode(int nPopIdx, int nEventIdx) const;

  bool concatenatePops(int nPrevPopIdx, int nNextPopIdx);
  bool connectCoalEvents(int nLeftPopIdx, int nLeftSonIdx,
                         int nRightPopIdx, int nRightSonIdx,
                         int nParentPopIdx, int nParentIdx);
  void addEventChainToDag(T* pEventStart, int nLen);

  void initPopulation(int nPopIdx, T* pFirstEvent);
  void appendEvent(int nPopIdx, T* pFirstEvent);
  EventsDAGNode<T>* getFirstGenEventInPop(int nPopIdx) const;
};

/*---------------------------------------------------------------------------
 * DTOR
 */
template<class T>
EventsDAG<T>::~EventsDAG()
{
  EventsDAGNode<T>* pCurrEvent = nullptr;
  EventsDAGNode<T>* pVictimEvent = nullptr;
  auto iCurr = this->begin();
  auto iEnd = this->end();
  for( ; iCurr != iEnd; ++iCurr )
  {
    pCurrEvent = *iCurr;
    while( pCurrEvent != nullptr )
    {
      pVictimEvent = pCurrEvent;
      pCurrEvent = pCurrEvent->getCoalParent();
      delete pVictimEvent;
    }
  }
  this->clear();
}
/*---------------------------------------------------------------------------
 */
template<class T> void
EventsDAG<T>::initDAG(int nNumOfPops)
{
  this->reserve(nNumOfPops);
  for( int i = 0; i < nNumOfPops; ++i )
    (*this)[i] = nullptr;
}

/*---------------------------------------------------------------------------
 */
template<class T> EventsDAGNode<T>*
EventsDAG<T>::getNode(int nPopIdx, int nEventId) const
{
  if( nPopIdx >= this->size() || nEventId < 0 || nPopIdx < 0 )
    return nullptr;
  EventsDAGNode<T>* pRes = *(this->begin() + nPopIdx);
  while(nullptr != pRes)
  {
    if( pRes->getContent()->getId() == nEventId )
      break;
    pRes = pRes->getNextGenEvent();
  }
  return pRes;
}

/*---------------------------------------------------------------------------
 */
template<class T> void
EventsDAG<T>::addEventChainToDag(T* pEventStart, int nLen)
{
  EventsDAGNode<T>* pCurrEvent = nullptr;
  pCurrEvent = new EventsDAGNode<T>(new T);
  pCurrEvent->getContent()->setType(POP_START);
  this->push_back(pCurrEvent);

  EventsDAGNode<T>* pPrevEvent = pCurrEvent;
  for( int i = 0; i < nLen; ++i )
  {
    pCurrEvent = new EventsDAGNode<T>(&(pEventStart[i]));
    pPrevEvent->setNextGenEvent(pCurrEvent);
    pCurrEvent->setPrevGenEvent(pPrevEvent);
    pPrevEvent = pCurrEvent;
  }

  pCurrEvent = new EventsDAGNode<T>(new T);
  pCurrEvent->getContent()->setType(POP_END);
  pPrevEvent->setNextGenEvent(pCurrEvent);
  pCurrEvent->setPrevGenEvent(pPrevEvent);
}

/*---------------------------------------------------------------------------
 */
template<class T> void
EventsDAG<T>::initPopulation(int nPopIdx, T* pFirstEvent)
{
  //Create POP START, POP END events
  auto* pPopStartVirtEvent = new T;
  auto* pPopStartDAGEvent = new EventsDAGNode<T>(pPopStartVirtEvent);
  (*this)[nPopIdx] = pPopStartDAGEvent;
  pPopStartDAGEvent->getContent()->setType (POP_START);
  auto* pPopEndVirtEvent = new T;
  auto* pPopEndDAGEvent = new EventsDAGNode<T>(pPopEndVirtEvent);
  pPopEndDAGEvent->getContent()->setType (POP_END);
  pPopStartDAGEvent->setNextGenEvent(pPopEndDAGEvent);
  pPopEndDAGEvent->setPrevGenEvent(pPopStartDAGEvent);
  //Insert the actual event in between
  this->appendEvent(nPopIdx, pFirstEvent);
}

/*---------------------------------------------------------------------------
 */
template<class T> void
EventsDAG<T>::appendEvent(int nPopIdx, T* pEvent)
{
  auto* pDAGEvent = new EventsDAGNode<T>(pEvent);
  EventsDAGNode<T>* pCurrDAGEvent = (*this)[nPopIdx];
  assert(pCurrDAGEvent->getContent()->getType() == POP_START);
  while( pCurrDAGEvent->getNextGenEvent()->getType() != POP_END )
    pCurrDAGEvent = pCurrDAGEvent->getNextGenEvent();

  pDAGEvent->setPrevGenEvent(pCurrDAGEvent);
  pDAGEvent->setNextGenEvent(pCurrDAGEvent->getNextGenEvent());
  pCurrDAGEvent->setNextGenEvent(pDAGEvent);
  pDAGEvent->getNextGenEvent()->setPrevGenEvent(pDAGEvent);
}
/*---------------------------------------------------------------------------
 */
template<class T> EventsDAGNode<T>*
EventsDAG<T>::getFirstGenEventInPop(int nPopIdx) const
{
  EventsDAGNode<T>* pPopStartEvent = (*this)[nPopIdx];
  assert(pPopStartEvent != nullptr);
  return pPopStartEvent->getNextGenEvent();
}

/*=============================================================================
 *
 * Repository of all DAGs
 *
 *---------------------------------------------------------------------------*/
template<class T>
class DAGsPerLocus : public vector< EventsDAG<T>* >
{
public:
  DAGsPerLocus( int nNumOfLoci, int nNumOfPops );
  virtual ~DAGsPerLocus();
  EventsDAG<T>* getDAG(int nLocusIdx) const;

protected:
  void initDAGs( int nNumOfPops );
};

/*---------------------------------------------------------------------------
 * CTOR
 */
template<class T>
DAGsPerLocus<T>::DAGsPerLocus(int nOfLoci, int nNumOfPops)
{
  this->reserve(nOfLoci);
  for( int i = 0; i < nOfLoci; ++i )
    this->push_back(new EventsDAG<T>);
  initDAGs(nNumOfPops);
}

/*---------------------------------------------------------------------------
 * DTOR
 */
template<class T>
DAGsPerLocus<T>::~DAGsPerLocus()
{
  for_each( this->begin(),
            this->end(),
            [](EventsDAG<T>* v){ delete v; } );
}

/*---------------------------------------------------------------------------
 * Initializer.
 */
template<class T>
void DAGsPerLocus<T>::initDAGs(int nNumOfPops)
{
  for_each( this->begin(),
            this->end(),
            [nNumOfPops](EventsDAG<T>* v){ v->initDAG( nNumOfPops ); } );
}

/*=============================================================================
 *
 * Single Events DAG node definition. Getters/Setters.
 *
 *---------------------------------------------------------------------------*/
template<class T> EventsDAGNode<T>*
EventsDAGNode<T>::getPrevGenEvent() const
{
  return pPrevGenEvent_;
}
template<class T>
void EventsDAGNode<T>::setPrevGenEvent(EventsDAGNode<T>* pPrevGenEvent)
{
  pPrevGenEvent_ = pPrevGenEvent;
}
template<class T> EventsDAGNode<T>*
EventsDAGNode<T>::getNextGenEvent() const
{
  return pNextGenEvent_;
}

template<class T> void
EventsDAGNode<T>::setNextGenEvent( EventsDAGNode<T>* pNextGenEvent)
{
  pNextGenEvent_ = pNextGenEvent;
}

template<class T> EventsDAGNode<T>*
EventsDAGNode<T>::getCoalLeftSon() const
{
  return pCoalLeftSon_;
}

template<class T> void
EventsDAGNode<T>::setCoalLeftSon( EventsDAGNode<T>* pCoalLeftSon )
{
  pCoalLeftSon_ = pCoalLeftSon;
}

template<class T> EventsDAGNode<T>*
EventsDAGNode<T>::getCoalRightSon() const
{
  return pCoalRightSon_;
}

template<class T> void
EventsDAGNode<T>::setCoalRightSon( EventsDAGNode<T>* pCoalRightSon)
{
  pCoalRightSon_ = pCoalRightSon;
}

template<class T> EventsDAGNode<T>*
EventsDAGNode<T>::getCoalParent() const
{
  return pCoalParent_;
}

template<class T> void
EventsDAGNode<T>::setCoalParent( EventsDAGNode<T>* pCoalParent )
{
  pCoalParent_ = pCoalParent;
}

template<class T> T*
EventsDAGNode<T>::getContent() const
{
  return pContent_;
}

template<class T> void
EventsDAGNode<T>::setContent( T* pContent )
{
  pContent_ = pContent;
}
//============================ END OF FILE ====================================
