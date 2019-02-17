
#include "EventsGraph.h"
#include "DbgErrMsgEvents.h"
#include "DataLayerConstants.h"


/*
    EventsGraph constructor
    Allocates EventNode objects. Each Node contains an Event object
*/
EventsGraph::EventsGraph(int genealogyID, int nNodes,
        DATA_SETUP* pDataSetup, DATA_STATE* pDataState,
        GENETREE_MIGS* pGenetreeMigs) :
            genealogyID_(genealogyID),
            nNodes_(nNodes),
            pDataSetup_(pDataSetup),
            pPopTree_(pDataSetup->popTree_),
            pDataState_(pDataState),
            pGenetreeMigs_(pGenetreeMigs) {


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
void EventsGraph::initializeGraph() {

    int nPops = pPopTree_->numPops; //N
    int rootPop = pPopTree_->rootPop; //root population

    //for each population define a start and end events
    for (int pop = 0; pop < nPops; ++pop) {

        //start events

        //set the i-th cell (i=pop) to be a start event
        eventNodes_[pop].getEvent().setType(EventType::POP_START);

        //next of start event points to the end event of same pop
        eventNodes_[pop].setNextGenEvent(eventNodes_ + nPops);

        //prev of start event points to NULL
        eventNodes_[pop].setPrevGenEvent(nullptr);

        //end events

        //set the (i+N)-th cell (i=pop) to be an end event
        eventNodes_[pop + nPops].getEvent().setType(EventType::POP_END);

        //set num lineages - it is important to initialize 0 incoming lineages
        eventNodes_[pop + nPops].getEvent().setNumLineages(0);

        //set elapsed time - distinguish between rootPop and rest of pops
        double t;
        if (pop == rootPop) {
            t = OLDAGE - pPopTree_->pops[rootPop]->age;
        }
        else {
            t = pPopTree_->pops[pop]->father->age - pPopTree_->pops[pop]->age;
        }
        eventNodes_[pop + nPops].getEvent().setElapsedTime(t);

        //next of end event points to the first event of parent pop.
        //next of root's end points to NULL
        if (pop == rootPop) {
            eventNodes_[pop + nPops].setNextGenEvent(nullptr);
        }
        else {
            int parent_id = pPopTree_->pops[pop]->father->id; //get parent id
            eventNodes_[pop + nPops].setNextGenEvent(eventNodes_ + parent_id);
        }

        //prev of end event points to the start event of same pop
        eventNodes_[pop + nPops].setPrevGenEvent(eventNodes_ - nPops);
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
    if (pNodesPool_ == nullptr) {
        EVENTS_FATAL_0015
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
    Returns the start event of a specified population.
    Start events allocated in the n-th cell of the events array,
    where n = population id
    @param: population id
    @return: pointer to a start EventNode

*/
EventNode* EventsGraph::getStartEvent(int pop) {
    return eventNodes_ + pop;
}


/*
    Creates an event of type "dummy" before a specified event
    @param: pointer to EventNode, elapsed time of new event,
            and (optional) event type

*/
EventNode* EventsGraph::createEventBefore(
        EventNode* pNode, double elapsed_time, EventType type) {

    //get a new event from pool (default type is DUMMY)
    EventNode* pNewNode = this->getNodeFromPool();

    //set event type
    pNewNode->getEvent().setType(type);

    //set number of lineages of new event to same number as given event
    int nLineages = pNode->getEvent().getNumLineages();
    pNewNode->getEvent().setNumLineages(nLineages);

    //set elapsed time of new event
    pNewNode->getEvent().setElapsedTime(elapsed_time);

    //decrease the elapsed time of the given event
    pNode->getEvent().addElapsedTime(-elapsed_time);

    //set pointers
    pNode->getPrevGenEvent()->setNextGenEvent(pNewNode);
    pNewNode->setPrevGenEvent(pNode->getPrevGenEvent());
    pNewNode->setNextGenEvent(pNode);
    pNode->setPrevGenEvent(pNewNode);

    return pNewNode;
}


/*
   Creates a new event in specified population at given time.
   Changes only elapsed time of subsequent event. Makes no other
   changes to events graph. numLineanges of new event is as the one
   of subsequent event.
   @param: population id, age, eventType
*/
EventNode* EventsGraph::createEvent(
        int pop, double age, EventType type) {

    ///////////////////////////////////////////////////////////////////////////////
    //why not subtract the opposite???? popTree->pops[pop]->age - age?
    ///////////////////////////////////////////////////////////////////////////////
    double delta_time = age - pPopTree_->pops[pop]->age;

    if (delta_time < 0) {
        EVENTS_FATAL_0016
        return nullptr;
    }

    if(pop != pPopTree_->rootPop &&
        age > pPopTree_->pops[pop]->father->age + 0.000001) {
        EVENTS_FATAL_0017
        return nullptr;
    }

    //find a spot for a new event
    //loop while not reaching the end event of the population
    //and while elapsed time of current event is smaller than delta time
    EventNode* pNode = this->getStartEvent(pop);
    for (;
        pNode->getEvent().getType() != EventType::POP_END &&
        pNode->getEvent().getElapsedTime() < delta_time;
        pNode = pNode->getNextGenEvent())
    {
        delta_time -= pNode->getEvent().getElapsedTime();
    }

    if (pNode->getEvent().getElapsedTime() < delta_time)
    {
        if (pNode->getEvent().getElapsedTime() < (delta_time - 0.000001)) {
            EVENTS_FATAL_0018
        }
        delta_time = pNode->getEvent().getElapsedTime();
    }

    //create the new event in the found slot, with elapsed_time = delta_time
    return this->createEventBefore(pNode, delta_time, type);

}


/*
	Constructs events graph
	Constructs everything from scratch (except for the start and end events which are defined in initaliztion)
	Typically used only for initial genetrees or for testing.
	Records number of lineages only for first events in leaf populations.
	The rest are recorded by computeGenetreeStats
*/
int EventsGraph::constructEventGraph() {


    //initialise graph with start and end events
    this->initializeGraph();

    //TODO: samples start

    // migration node events
    int gen = genealogyID_; //shorter name;
    for (int i = 0; i < genetree_migs[gen].num_migs; i++) {

        int mig = pGenetreeMigs_[gen].living_mignodes[i];

        #ifdef DEBUG_EVENT_GRAPH
                printf("living migration %d: %d in gen %d.\n", i, mig, gen);
        #endif

        //get age of migration
        int age = pGenetreeMigs_[gen].mignodes[mig].age;

        //get target population
        int popIn = pGenetreeMigs_[gen].mignodes[mig].target_pop;

        //create an incoming migration event
        EventNode* pNodeIn = this->createEvent(popIn, age, EventType::IN_MIG);
        if (!pNodeIn) {
            EVENTS_FATAL_0022
        }
        //TODO: genetree_migs[gen].mignodes[mig].target_event = event;

        //get source population
        int popOut = pGenetreeMigs_[gen].mignodes[mig].source_pop;

        //create an outgoing migration event
        EventNode* pNodeOut = this->createEvent(popOut, age, EventType::OUT_MIG);
        if (!pNodeOut) {
            EVENTS_FATAL_0023
        }
        //TODO: genetree_migs[gen].mignodes[mig].target_event = event;
    }


    // coalescent node events
    for (int node = 0; node < 2 * pDataSetup_->numSamples - 1; node++) {

        int pop = nodePops[gen][node]; //TODO: get nodePops as argument

        // if node is a leaf
        if (node < pDataSetup_->numSamples) {

            #ifdef DEBUG_EVENT_GRAPH
                printf("L%d ",pop);
            #endif
        }
        else {

            int age = getNodeAge(pDataState_->lociData[gen], node);

            //age = gnodes[gen][node].age;

            #ifdef DEBUG_EVENT_GRAPH
                printf("C%d %d %f ", pop, node, age);
            #endif

            EventNode* pNode = createEvent(pop, age, EventType::COAL);
            if (!pNode) {
                EVENTS_FATAL_0024
            }

            //TODO nodeEvents[gen][node] = event;
        }
    }

    #ifdef DEBUG_EVENT_GRAPH
        printf("\n");
    #endif

    return 0;
}



