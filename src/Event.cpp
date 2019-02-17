
#include "Event.h"

/*
    Event constructor
*/
Event::Event() :
    type_(EventType::DUMMY),
    elapsedTime_(0.0),
    nLineages_(-1) {
}

/*
    Resets Event content
*/
void Event::reset() {
    type_ = EventType::DUMMY;
    elapsedTime_ = 0.0;
    nLineages_ = -1;
}

/*
    Adds elapsed time
*/
void Event::addElapsedTime(double delta) {
    this->elapsedTime_ += delta;
}

/******************************************************************************
 * getters / setters
******************************************************************************/

EventType Event::getType() const {
    return type_;
}

void Event::setType(EventType type) {
    Event::type_ = type;
}

double Event::getElapsedTime() const {
    return elapsedTime_;
}

void Event::setElapsedTime(double elapsedTime) {
    Event::elapsedTime_ = elapsedTime;
}

int Event::getNumLineages() const {
    return nLineages_;
}

void Event::setNumLineages(int nLineages) {
    Event::nLineages_ = nLineages;
}






