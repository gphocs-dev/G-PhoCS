
#ifndef G_PHOCS_DAGNODE_H
#define G_PHOCS_DAGNODE_H

#include "LocusDataConnector.h"

enum class EventType; //forward declaration

/*=============================================================================
 *
 * EventNode class
 *
 * EventNode is a single node in EventsGraph.
 *
 * Contains:
 * 1. Type of event.
 * 2. Elapsed time of event.
 * 3. Number of lineages (before the event).
 * 4. Population ID.
 * 5. Five pointers ("edges") to several other EventNode objects.
 * 6. A LocusDataConnector object.
 *===========================================================================*/
class EventNode {

private:

    EventType  type_;               //type of event
    double     elapsedTime_;        //elapsed time
    int        nLineages_;	        //number of lineages (before the event)

    int         popID_;            //population ID

    EventNode*  pPrevGenEvent_;    //pointer to previous genealogy event
    EventNode*  pNextGenEvent_;    //pointer to next genealogy event
    EventNode*  pCoalLeftSon_;     //pointer to left son in genealogy
    EventNode*  pCoalRightSon_;    //pointer to right son in genealogy
    EventNode*  pCoalParent_;      //pointer to parent in genealogy

    LocusDataConnector locusDataConnector_; //


public:

    //constructor
    EventNode();

    //reset EventNode members
    void reset();


    void incrementLineages();
    void decrementLineages();

    //add elapsed time
    void addElapsedTime(double delta);

    //return true if node is of the specified type
    bool isType(EventType type);

public:

    //getters / setters
    EventType getType() const;
    void setType(EventType type);

    double getElapsedTime() const;
    void setElapsedTime(double elapsedTime);

    int getNumLineages() const;
    void setNumLineages(int nLineages);

    int getPopID() const;
    void setPopID(int popID);

    EventNode* getPrevGenEvent() const;
    void setPrevGenEvent(EventNode* pPrevGenEvent);

    EventNode* getNextGenEvent() const;
    void setNextGenEvent(EventNode* pNextGenEvent);

    EventNode* getCoalLeftSon() const;
    void setCoalLeftSon(EventNode* pCoalLeftSon);

    EventNode* getCoalRightSon() const;
    void setCoalRightSon(EventNode* pCoalRightSon);

    EventNode* getCoalParent() const;
    void setCoalParent(EventNode* pCoalParent);

    //set left and right coal sons with same value
    void setCoalSons(EventNode* pSons);

};



/*
    Types of Event
*/
enum class EventType {
    COAL,
    IN_MIG,
    OUT_MIG,
    MIG_BAND_START, //later
    MIG_BAND_END, //later
    SAMPLES_START,
    END_CHAIN,
    POP_START,
    POP_END,
    DUMMY
};


#endif //G_PHOCS_DAGNODE_H
