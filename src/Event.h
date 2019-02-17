
#ifndef G_PHOCS_EVENTCLASS_H
#define G_PHOCS_EVENTCLASS_H


/*
    Types of Event
*/
enum class EventType;

/*=============================================================================
 *
 * Event class
 *
 * Event is defined by three parameters:
 * 1. Type of event.
 * 2. Elapsed time of event.
 * 3. Number of lineages (before the event).
 *
 *===========================================================================*/
class Event {

private:
    EventType  type_;           //type of event
    double     elapsedTime_;   //elapsed time
    int        nLineages_;	  //number of lineages (before the event)

public:

    //constructor
    Event();

    //reset Event content
    void reset();

    void incrementLineages();
    void decrementLineages();

    void addElapsedTime(double delta);

public:

    //getters / setters

    EventType getType() const;
    void setType(EventType type);

    double getElapsedTime() const;
    void setElapsedTime(double elapsedTime);

    int getNumLineages() const;
    void setNumLineages(int nLineages);

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


#endif //G_PHOCS_EVENTCLASS_H
