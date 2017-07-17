#pragma once

/*============================================================================
 File: EventsDag.h

 Directed Acyclic Graph that represents both Population and Genealogy trees.

 Version history:
 29-Jun-2017  evgenyidc      Initial version reworked.
 ============================================================================*/

#include <vector>
using namespace std;

/*-----------------------------------------------------------------------------
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

  EventsDAGNode(T* pContent):
      EventsDAGNode(){pContent_ = pContent;}

  T* getContent() const;
  void setContent( T* pContent);

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

/*-----------------------------------------------------------------------------
 *
 * EventsDAG declaration
 *
 *---------------------------------------------------------------------------*/
class Event;

template<class T>
class EventsDAG : public vector<EventsDAGNode<T>*>
{
public:
  EventsDAG(){}
  virtual ~EventsDAG();
  EventsDAGNode<T>* getNode(int nGenIdx, int nPopIdx) const;
  void addEventChainToDag(Event* pEventStart, int nLen);

};

/*---------------------------------------------------------------------------
 * DTOR
 */
template<class T>
EventsDAG<T>::~EventsDAG()
{
  EventsDAGNode<Event>* pCurrEvent = nullptr;
  EventsDAGNode<Event>* pVictimEvent = nullptr;
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
template<class T> EventsDAGNode<T>*
EventsDAG<T>::getNode(int nGenIdx, int nEventIdx) const
{
  if( nGenIdx >= this->size() )
    return nullptr;
  EventsDAGNode<T>* pRes = *(this->begin() + nGenIdx);
  for( int i = 1 ; i < nEventIdx; ++i )
  {
    if(nullptr == pRes)
      break;
    pRes = pRes->getNextGenEvent();
  }
  return pRes;
}

/*---------------------------------------------------------------------------
 */
template<class T> void
EventsDAG<T>::addEventChainToDag(Event* pEventStart, int nLen)
{
  EventsDAGNode<Event>* pStartChain = new EventsDAGNode<Event>(pEventStart);
  this->push_back(pStartChain);
  EventsDAGNode<Event>* pPrevEvent = pStartChain;
  for( int i = 0; i < nLen; ++i )
  {
    EventsDAGNode<Event>* pCurrEvent = new EventsDAGNode<Event>(&(pEventStart[i]));
    pPrevEvent->setNextGenEvent(pCurrEvent);
    pCurrEvent->setPrevGenEvent(pPrevEvent);
    pPrevEvent = pCurrEvent;
  }
}

/*-----------------------------------------------------------------------------
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
