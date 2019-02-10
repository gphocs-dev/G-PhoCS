
#include "EventNode.h"


/*
    EventNode constructor
    Sets all pointers to nullptr and creates an Event object
*/
EventNode::EventNode() :
        pPrevGenEvent_(nullptr),
        pNextGenEvent_(nullptr),
        pCoalLeftSon_(nullptr),
        pCoalRightSon_(nullptr),
        pCoalParent_(nullptr),
        event_() {
}

/*
    Resets EventNode content
    Sets all pointers to nullptr and resets Event object
*/
void EventNode::reset() {

    pPrevGenEvent_ = nullptr;
    pNextGenEvent_ = nullptr;
    pCoalLeftSon_ = nullptr;
    pCoalRightSon_ = nullptr;
    pCoalParent_ = nullptr;
    event_.reset();
}

/******************************************************************************
 * getters / setters
******************************************************************************/

EventNode* EventNode::getPrevGenEvent() const {
    return pPrevGenEvent_;
}

void EventNode::setPrevGenEvent(EventNode* pPrevGenEvent_) {
    EventNode::pPrevGenEvent_ = pPrevGenEvent_;
}

EventNode* EventNode::getNextGenEvent() const {
    return pNextGenEvent_;
}

void EventNode::setNextGenEvent(EventNode* pNextGenEvent_) {
    EventNode::pNextGenEvent_ = pNextGenEvent_;
}

EventNode* EventNode::getCoalLeftSon() const {
    return pCoalLeftSon_;
}

void EventNode::setCoalLeftSon(EventNode* pCoalLeftSon_) {
    EventNode::pCoalLeftSon_ = pCoalLeftSon_;
}

EventNode* EventNode::getCoalRightSon() const {
    return pCoalRightSon_;
}

void EventNode::setCoalRightSon(EventNode* pCoalRightSon_) {
    EventNode::pCoalRightSon_ = pCoalRightSon_;
}

EventNode* EventNode::getCoalParent() const {
    return pCoalParent_;
}

void EventNode::setCoalParent(EventNode* pCoalParent_) {
    EventNode::pCoalParent_ = pCoalParent_;
}

Event& EventNode::getEvent() {
    return event_;
}

void EventNode::setEvent(const Event& event) {
    EventNode::event_ = event;
}



