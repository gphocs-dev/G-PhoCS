
#include "EventNode.h"


/*
    EventNode constructor
    Creates an Event object
    Creates a LocusDataConnector object
    Set population ID to -1
    Sets all pointers to nullptr
*/
EventNode::EventNode() :
        type_(EventType::DUMMY),
        elapsedTime_(-1),
        nLineages_(-1),
        popID_(-1),
        locusDataConnector_(),
        pPrevGenEvent_(nullptr),
        pNextGenEvent_(nullptr),
        pCoalLeftSon_(nullptr),
        pCoalRightSon_(nullptr),
        pCoalParent_(nullptr){
}

/*
    Resets EventNode content
    Sets all pointers to nullptr and resets Event object
*/
void EventNode::reset() {

    type_ = EventType::DUMMY;
    elapsedTime_ = -1;
    nLineages_ = -1;

    popID_ = -1;

    locusDataConnector_.reset();

    pPrevGenEvent_ = nullptr;
    pNextGenEvent_ = nullptr;
    pCoalLeftSon_ = nullptr;
    pCoalRightSon_ = nullptr;
    pCoalParent_ = nullptr;

}

void EventNode::addElapsedTime(double delta) {
    this->elapsedTime_ += delta;

}

/*
    returns true if node is of the specified type
*/
bool EventNode::isType(EventType type) {
    return type_ == type;
}

/******************************************************************************
 * getters / setters
******************************************************************************/

EventType EventNode::getType() const {
    return type_;
}

void EventNode::setType(EventType type) {
    type_ = type;
}

double EventNode::getElapsedTime() const {
    return elapsedTime_;
}

void EventNode::setElapsedTime(double elapsedTime) {
    elapsedTime_ = elapsedTime;
}

int EventNode::getNumLineages() const {
    return nLineages_;
}

void EventNode::setNumLineages(int nLineages) {
    nLineages_ = nLineages;
}

int EventNode::getPopID() const {
    return popID_;
}

void EventNode::setPopID(int popID) {
    popID_ = popID;
}

EventNode* EventNode::getPrevGenEvent() const {
    return pPrevGenEvent_;
}

void EventNode::setPrevGenEvent(EventNode* pPrevGenEvent) {
    pPrevGenEvent_ = pPrevGenEvent;
}

EventNode* EventNode::getNextGenEvent() const {
    return pNextGenEvent_;
}

void EventNode::setNextGenEvent(EventNode* pNextGenEvent) {
    pNextGenEvent_ = pNextGenEvent;
}

EventNode* EventNode::getCoalLeftSon() const {
    return pCoalLeftSon_;
}

void EventNode::setCoalLeftSon(EventNode* pCoalLeftSon) {
    pCoalLeftSon_ = pCoalLeftSon;
}

EventNode* EventNode::getCoalRightSon() const {
    return pCoalRightSon_;
}

void EventNode::setCoalRightSon(EventNode* pCoalRightSon) {
   pCoalRightSon_ = pCoalRightSon;
}

EventNode* EventNode::getCoalParent() const {
    return pCoalParent_;
}

void EventNode::setCoalParent(EventNode* pCoalParent) {
    pCoalParent_ = pCoalParent;
}

void EventNode::setCoalSons(EventNode *pSons) {
    pCoalRightSon_ = pSons;
    pCoalLeftSon_ = pSons;
}







