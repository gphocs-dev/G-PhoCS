
#include "EventsGraph.h"


/*
    EventsGraph constructor
    Allocates EventNode objects. Each Node contains an Event object
*/
EventsGraph::EventsGraph(DATA_SETUP* pDataSetup, int nNodes) :
        pDataSetup_(pDataSetup), nNodes_(nNodes) {


    //allocate n EventsGraph nodes consisting the events EventsGraph
    this->eventNodes_ = new EventNode[nNodes];

    //for each EventsGraph node, set its prev and next genealogy events
    for (int i = 0; i < nNodes; ++i) {

        if (i == 0){
            eventNodes_[i].setNextGenEvent(eventNodes_ + 1);
            eventNodes_[i].setPrevGenEvent(nullptr);
        }
        else if (i == nNodes-1) {
            eventNodes_[i].setNextGenEvent(nullptr);
            eventNodes_[i].setPrevGenEvent(eventNodes_ - 1);
        }
        else {
            eventNodes_[i].setNextGenEvent(eventNodes_ + 1); //point to next cell
            eventNodes_[i].setPrevGenEvent(eventNodes_ - 1); //point to prev cell
        }
    }

    //events pool points to the first EventsGraph node in eventsDAG
    // (in initialization all events in pool)
    pNodesPool_ = eventNodes_;

}

/*
    EventsGraph class destructor
*/
EventsGraph::~EventsGraph() {

    //delete array of DAGNodes
    delete eventNodes_;

}

/*
    Initializes EventsGraph
    Defines the first N events in the EventsGraph (N = nPops) as start events
    and the next N events as end events

*/
void EventsGraph::initializeDAG() {

    int nPops = pDataSetup_->popTree->numPops; //N

    //for each population define a start and end events
    for( int i = 0; i < nPops; ++i )
    {
        //set the i and the i+N cells to be start and end events, accordingly
        eventNodes_[i].getEvent().setType(EventType::POP_START);
        eventNodes_[i + nPops].getEvent().setType(EventType::POP_END);

        //next of start event points to the end event of same pop
        //prev of start event points to NULL
        eventNodes_[i].setNextGenEvent(eventNodes_ + nPops);
        eventNodes_[i].setPrevGenEvent(nullptr);

        //next of end event points to the first event of parent pop.
        //next of root's end points to NULL
        if (i == nPops-1)
            eventNodes_[i + nPops].setNextGenEvent(nullptr);
        else {
            int parent_id = pDataSetup_->popTree->popArray[i].father->id; //get parent id
            eventNodes_[i + nPops].setNextGenEvent(eventNodes_ + parent_id);
        }

        //prev of end event points to the start event of same pop
        eventNodes_[i + nPops].setPrevGenEvent(eventNodes_ - nPops);
    }

    //promote the free events pointer to the 2N cell
    pNodesPool_ = eventNodes_ + 2*nPops;
    pNodesPool_->setPrevGenEvent(nullptr);

}

/*
    Returns a free EventNode from the nodes pool
    @return: pointer to a new EventNode

*/
EventNode* EventsGraph::getNodeFromPool() {

    //if no free events
    if (pNodesPool_ == nullptr){
        //TODO: print error and exit
    }

    EventNode* pEventNode = pNodesPool_;

    //update pointer to events pool
    pNodesPool_ = pNodesPool_->getNextGenEvent();
    pNodesPool_->setPrevGenEvent(nullptr);

    pEventNode->setNextGenEvent(nullptr);

    return pEventNode;
}


/*
    Adds a new node to the free nodes pool
    @param: pointer to new a EventNode

*/
void EventsGraph::addNodeToPool(EventNode* pEventNode) {

    //reset EventNode content
    pEventNode->reset();

    //set pointers
    pNodesPool_->setPrevGenEvent(pEventNode);
    pEventNode->setNextGenEvent(pNodesPool_);
    pNodesPool_ = pEventNode;
}


/*
    Creates an event of type "dummy" before a specified event
    @param: pointer to event (EventNode) and elapsedTime

*/
void EventsGraph::createEventBefore(EventNode* pNode, double elapsedTime) {

    //get a new event from pool (default type is DUMMY)
    EventNode* pNewNode = this->getNodeFromPool();

    //set number of lineages of new event to same number as given event
    int nLineages = pNode->getEvent().getNumLineages();
    pNewNode->getEvent().setNumLineages(nLineages);

    //set elapsed time of new event
    pNewNode->getEvent().setElapsedTime(elapsedTime);

    //set elapsed time of given event
    pNode->getEvent().addElapsedTime(-elapsedTime);


    //set pointers


}

