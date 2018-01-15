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
class EventChains;
struct _POPULATION_TREE;
typedef _POPULATION_TREE PopulationTree;
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
  void createEventBefore( int nPopIdx,
                          T* pNewEvnt,
                          EventsDAGNode<T>* pNextNode );

  //EventsDAGNode<T>* getFirstNode(int nPopIdx) const;
  //EventsDAGNode<T>* getLastNode(int nPopIdx) const;
  EventsDAGNode<T>* getNode( int nPopIdx, const T* pEvent ) const;

  bool concatenatePops(int nPrevPopIdx, int nNextPopIdx);
  bool connectCoalEvents(int nLeftPopIdx, int nLeftSonIdx,
                         int nRightPopIdx, int nRightSonIdx,
                         int nParentPopIdx, int nParentIdx);

  void appendEvent(int nPopIdx, T* pFirstEvent);
  void importEventChains(int nGenIdx,
                         const EventChains* pEventChains,
                         const PopulationTree* pPopTree);
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
  {
    auto* pVirtPopStartEvent = new T();
    pVirtPopStartEvent->setType(POP_START);
    auto* pVirtPopEndEvent = new T();
    pVirtPopEndEvent->setType(POP_END);
    auto* pStartDAGNode = new EventsDAGNode<T>( pVirtPopStartEvent );
    auto* pEndDAGNode = new EventsDAGNode<T>( pVirtPopEndEvent );
    pStartDAGNode->setNextGenEvent( pEndDAGNode );
    pEndDAGNode->setPrevGenEvent( pStartDAGNode );
    this->push_back( pStartDAGNode );
  }
}

/*---------------------------------------------------------------------------
 */
template<class T> void
EventsDAG<T>::createEventBefore( int nPopIdx,
                                 T* pNewEvnt,
                                 EventsDAGNode<T>* pNextNode )
{
  assert( 0 <= nPopIdx && nPopIdx < this->size() );
  EventsDAGNode<T>* pPrevNode = pNextNode->getPrevGenEvent();
  auto* pNewDAGNode = new EventsDAGNode<T>(pNewEvnt);
  if( nullptr == pPrevNode )
    (*this)[nPopIdx] = pNewDAGNode;
  else
    pPrevNode->setNextGenEvent(   pNewDAGNode );
  pNewDAGNode->setPrevGenEvent( pPrevNode   );
  pNewDAGNode->setNextGenEvent( pNextNode   );
  pNextNode->setPrevGenEvent(   pNewDAGNode );
}

/*---------------------------------------------------------------------------
 */
template<class T> EventsDAGNode<T>*
EventsDAG<T>::getNode( int nPopIdx, const T* pEvent ) const
{
  if( nPopIdx >= this->size() || pEvent == nullptr || nPopIdx < 0 )
    return nullptr;
  EventsDAGNode<T>* pCurrNode = *(this->begin() + nPopIdx);
  while(nullptr != pCurrNode)
  {
    //Pointer comparison here
    if( pCurrNode->getContent() == pEvent )
      break;
    pCurrNode = pCurrNode->getNextGenEvent();
  }
  return pCurrNode;
}

/*---------------------------------------------------------------------------
 */
template<class T> void
EventsDAG<T>::importEventChains(int nGenIdx,
                                const EventChains* pEventChains,
                                const PopulationTree* pPopTree)
{
  int nPopIdx = -1;
  int nEventIdx = -1;
  for(nPopIdx = 0; nPopIdx < pPopTree->numPops; ++nPopIdx)
  {
    for (nEventIdx = event_chains[nGenIdx].first_event[nPopIdx];
         nEventIdx >= 0;
         nEventIdx = event_chains[nGenIdx].events[nEventIdx].getNextIdx())
    {
      T &CurrEvent = event_chains[nGenIdx].events[nEventIdx];
      if(    CurrEvent.getType() != MIG_BAND_START
          && CurrEvent.getType() != MIG_BAND_END )
      {
        appendEvent(nPopIdx, &CurrEvent);
      }
    }
  }

}

/*---------------------------------------------------------------------------
 */
template<class T> void
EventsDAG<T>::appendEvent(int nPopIdx, T* pEvent)
{
  auto* pNewDAGNode = new EventsDAGNode<T>(pEvent);
  EventsDAGNode<T>* pCurrDAGNode = (*this)[nPopIdx];
  assert( nullptr != pCurrDAGNode );
  while(    nullptr != pCurrDAGNode->getNextGenEvent()
         && (    POP_END   != pCurrDAGNode->getNextGenEvent()->getType()
              && END_CHAIN != pCurrDAGNode->getNextGenEvent()->getType() ) )
    pCurrDAGNode = pCurrDAGNode->getNextGenEvent();

  pNewDAGNode->setPrevGenEvent(pCurrDAGNode);
  pNewDAGNode->setNextGenEvent(pCurrDAGNode->getNextGenEvent());
  pCurrDAGNode->setNextGenEvent(pNewDAGNode);
  pNewDAGNode->getNextGenEvent()->setPrevGenEvent(pNewDAGNode);
}
/*---------------------------------------------------------------------------
 */
template<class T> EventsDAGNode<T>*
EventsDAG<T>::getFirstGenEventInPop(int nPopIdx) const
{
  EventsDAGNode<T>* pPopStartNode = (*this)[nPopIdx];
  assert(pPopStartNode != nullptr);
  EventsDAGNode<T>* pResult = pPopStartNode->getNextGenEvent();
  return POP_END == pResult->getType() ? nullptr : pResult;
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
  void importEventChains(int nGenIdx,
                         const EventChains* pEventChains,
                         const PopulationTree* pPopTree);
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
  this->initDAGs(nNumOfPops);
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

/*---------------------------------------------------------------------------
 * Data importer
 */
template<class T> void
DAGsPerLocus<T>::importEventChains(int nGenIdx,
                                   const EventChains* pEventChains,
                                   const PopulationTree* pPopTree)
{
  (*this)[nGenIdx]->importEventChains(nGenIdx, pEventChains,pPopTree);
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
