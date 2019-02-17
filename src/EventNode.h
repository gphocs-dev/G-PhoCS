
#ifndef G_PHOCS_DAGNODE_H
#define G_PHOCS_DAGNODE_H

#include "Event.h"

/*=============================================================================
 *
 * EventNode class
 *
 * EventNode is a single node in EventsGraph.
 *
 * Contains:
 * 1. An Event object.
 * 2. Five pointers ("edges") to several other EventNode objects.
 *
 *===========================================================================*/
class EventNode {

private:

    Event       event_;            //event object
    EventNode*  pPrevGenEvent_;    //pointer to previous genealogy event
    EventNode*  pNextGenEvent_;    //pointer to next genealogy event
    EventNode*  pCoalLeftSon_;     //pointer to left son in genealogy
    EventNode*  pCoalRightSon_;    //pointer to right son in genealogy
    EventNode*  pCoalParent_;      //pointer to parent in genealogy


public:

    //constructor
    EventNode();

    //reset EventNode members
    void reset();

public:

    //getters / setters
    EventNode* getPrevGenEvent() const;
    void setPrevGenEvent(EventNode* pPrevGenEvent_);

    EventNode* getNextGenEvent() const;
    void setNextGenEvent(EventNode* pNextGenEvent_);

    EventNode* getCoalLeftSon() const;
    void setCoalLeftSon(EventNode* pCoalLeftSon_);

    EventNode* getCoalRightSon() const;
    void setCoalRightSon(EventNode* pCoalRightSon_);

    EventNode* getCoalParent() const;
    void setCoalParent(EventNode* pCoalParent_);

    Event& getEvent();
    void setEvent(const Event& event);

};


#endif //G_PHOCS_DAGNODE_H
