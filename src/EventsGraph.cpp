
#include "EventsGraph.h"
#include "DbgErrMsgEvents.h"
#include "DataLayerConstants.h"


/*
    EventsGraph constructor
    Allocates EventNode objects.
    Concatenates nodes to each other.
*/
EventsGraph::EventsGraph(int genealogyID, int nNodes,
        DATA_SETUP* pDataSetup, DATA_STATE* pDataState,
        GENETREE_MIGS* pGenetreeMigs) :
            locusID_(genealogyID),
            nNodes_(nNodes),
            pDataSetup_(pDataSetup),
            pPopTree_(pDataSetup->popTree),
            pDataState_(pDataState),
            pGenetreeMigs_(pGenetreeMigs) {


    //allocate n EventsGraph nodes consisting the events EventsGraph
    this->eventsGraph_ = new EventNode[nNodes];

    //for each EventsGraph node, set its prev and next genealogy events
    for (int i = 0; i < nNodes; ++i) {

        if (i == 0){
            eventsGraph_[i].setNextGenEvent(eventsGraph_ + 1);
            eventsGraph_[i].setPrevGenEvent(nullptr);
        }
        else if (i == nNodes-1) {
            eventsGraph_[i].setNextGenEvent(nullptr);
            eventsGraph_[i].setPrevGenEvent(eventsGraph_ - 1);
        }
        else {
            eventsGraph_[i].setNextGenEvent(eventsGraph_ + 1); //point to next cell
            eventsGraph_[i].setPrevGenEvent(eventsGraph_ - 1); //point to prev cell
        }
    }

    //events pool points to the first EventsGraph node in eventsDAG
    // (in initialization all events in pool)
    pNodesPool_ = eventsGraph_;

    //number of coalescent event is constant: 2*N -1, where N is num sumples
    coalEvents_.reserve(2 * pDataSetup_->numSamples - 1);
    
}

/*
    EventsGraph class destructor
*/
EventsGraph::~EventsGraph() {
    //delete array of DAGNodes
    delete eventsGraph_;
}

/*
    Initializes EventsGraph
    In EventsGraph array, define:
     1. The first N cells to be start events (N = num pops).
     2. The next N cells to be end events.
    Total: 2N cells are reserved to start / end events.

*/
void EventsGraph::initialiseGraph() {

    int nPops = pPopTree_->numPops; //N
    int rootPop = pPopTree_->rootPop; //root population

    //for each population define a start and end events
    for (int pop = 0; pop < nPops; ++pop) {

        //start events

        //set the i-th cell (i=pop) to be a start event
        eventsGraph_[pop].setType(EventType::POP_START);

        //next of start event points to the end event of same pop
        eventsGraph_[pop].setNextGenEvent(eventsGraph_ + nPops);

        //prev of start event points to NULL
        eventsGraph_[pop].setPrevGenEvent(nullptr);

        //end events

        //set the (i+N)-th cell (i=pop) to be an end event
        eventsGraph_[pop + nPops].setType(EventType::POP_END);

        //set num lineages - it is important to initialize 0 incoming lineages
        eventsGraph_[pop + nPops].setNumLineages(0);

        //set elapsed time - distinguish between rootPop and rest of pops
        if (pop == rootPop) {
            double t = OLDAGE - pPopTree_->pops[rootPop]->age;
            eventsGraph_[pop + nPops].setElapsedTime(t);
        }
        else {
            double t = pPopTree_->pops[pop]->father->age
                        - pPopTree_->pops[pop]->age;
            eventsGraph_[pop + nPops].setElapsedTime(t);
        }

        //next of end event points to the first event of parent pop.
        //next of root's end points to NULL
        if (pop == rootPop) {
            eventsGraph_[pop + nPops].setNextGenEvent(nullptr);
        }
        else {
            int parent_id = pPopTree_->pops[pop]->father->id; //get parent id
            eventsGraph_[pop + nPops].setNextGenEvent(eventsGraph_ + parent_id);
        }

        //prev of end event points to the start event of same pop
        eventsGraph_[pop + nPops].setPrevGenEvent(eventsGraph_ + pop);
    }

    //promote the free events pointer to the 2N cell
    pNodesPool_ = eventsGraph_ + 2*nPops;
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
    return eventsGraph_ + pop;
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
    pNewNode->setType(type);

    //set number of lineages of new event to same number as given event
    int nLineages = pNode->getNumLineages();
    pNewNode->setNumLineages(nLineages);

    //set elapsed time of new event
    pNewNode->setElapsedTime(elapsed_time);

    //decrease the elapsed time of the given event
    pNode->addElapsedTime(-elapsed_time);

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
        pNode->getType() != EventType::POP_END &&
        pNode->getElapsedTime() < delta_time;
        pNode = pNode->getNextGenEvent())
    {
        delta_time -= pNode->getElapsedTime();
    }

    if (pNode->getElapsedTime() < delta_time)
    {
        if (pNode->getElapsedTime() < (delta_time - 0.000001)) {
            EVENTS_FATAL_0018
        }
    }

    //create the new event in the found slot, with elapsed_time = delta_time
    return this->createEventBefore(pNode, delta_time, type);

}


/*
	Constructs events graph
	Constructs everything from scratch
	Typically used only for initial genetrees or for testing.
	Records number of lineages only for first events in leaf populations.
	The rest are recorded by computeGenetreeStats
*/
int EventsGraph::constructEventsGraph() {

    //initialize graph with start and end events
    this->initialiseGraph();

    //create samples start events (for ancient samples)
    for (int pop = 0; pop < pPopTree_->numCurPops; pop++) {

        //create event
        EventNode* pSamplesStart =
                this->createEvent(pop,
                                  pPopTree_->pops[pop]->sampleAge,
                                  EventType::SAMPLES_START);

        //coalescent pointers of samplesStart are null
        // (this is the default, no need to code)
    }

    //create coalescent events
    int nSamples = pDataSetup_->numSamples;
    for (int node = 0; node < 2*nSamples-1; node++) {

        //get population of node
        int pop = nodePops[locusID_][node];

        // if node is not a leaf
        if (!this->isLeaf(node)) {

            //get age of node
            int age = getNodeAge(getLocusData(), node);

            #ifdef DEBUG_EVENT_GRAPH
                        printf("C%d %d %f ", pop, node, age);
            #endif

            //create a coalescent event
            EventNode* pCoal = this->createEvent(pop, age, EventType::COAL);
            if (!pCoal) {
                EVENTS_FATAL_0024
            }

            //save pointer in coalescent events array
            coalEvents_[node] = pCoal;

        } else {
            #ifdef DEBUG_EVENT_GRAPH
                        printf("L%d ",pop);
            #endif
        }
    }

    //set the coalescent pointers of the coalescent events
    //iterate on internal nodes only (and not on leaves)
    for (int node = nSamples; node < 2*nSamples-1; node++) {

        //get nodes and corresponding eventNodes of the sons
        //eventNode can be of type coalescent or samplesStart
        int nodeLeftSon = getNodeSon(this->getLocusData(), node, 0);
        int nodeRightSon = getNodeSon(this->getLocusData(), node, 1);
        EventNode* pLeftSon = this->getEventByNode(nodeLeftSon);
        EventNode* pRightSon = this->getEventByNode(nodeRightSon);

        //get node and corresponding eventNode of the father
        //if current node is the root - parent points to null
        int nodeFather;
        EventNode* pFather;
        if ( node != getLocusRoot(this->getLocusData()) ) {
            nodeFather = getNodeFather(this->getLocusData(), node);
            pFather = this->getCoalEvent(nodeFather);
        }
        else {
            pFather = nullptr;
        }

        //get eventNode of current node and set its coalescent pointers
        EventNode* pCoal = this->getCoalEvent(node);

        //set pointers of current event node
        pCoal->setCoalParent(pFather);
        pCoal->setCoalLeftSon(pLeftSon);
        pCoal->setCoalRightSon(pRightSon);

    }

    //create migration events and if needed reset coalescent pointers
    // of coalescent events
    for (int node = 0; node < 2*nSamples-1; node++) {

        //define upper and bottom pointers (the upper the closer to root)
        EventNode* pBottom = this->getEventByNode(node);
        EventNode* pUpper = pBottom->getCoalParent();

        //find migration above current node and after specified time
        int pop = nodePops[locusID_][node];
        int mig = findFirstMig(locusID_, node, pPopTree_->pops[pop]->age);

        //while there are migration events on the edge above current node
        while (mig != -1) {

            //create migration events (source and target)

            //get age of migration, and target and source populations
            double age = pGenetreeMigs_[locusID_].mignodes[mig].age;
            int target_pop = pGenetreeMigs_[locusID_].mignodes[mig].target_pop;
            int source_pop = pGenetreeMigs_[locusID_].mignodes[mig].source_pop;

            //create an incoming migration event
            EventNode* pMigIn = this->createEvent(target_pop, age, EventType::IN_MIG);
            if (!pMigIn) {
                EVENTS_FATAL_0022
            }

            //create an outgoing migration event
            EventNode* pMigOut = this->createEvent(source_pop, age, EventType::OUT_MIG);
            if (!pMigOut) {
                EVENTS_FATAL_0023
            }

            //parent of incoming migration points to outgoing migration
            //sons of incoming migration point to bottom node
            pMigIn->setCoalParent(pMigOut);
            pMigIn->setCoalSons(pBottom);

            //parent of outgoing migration points to upper node
            //sons of outgoing migration point to outgoing migration
            pMigOut->setCoalParent(pUpper);
            pMigOut->setCoalSons(pMigIn);

            //promote bottom node to current outgoing migration
            pBottom = pUpper;

            //find next migration (after time of current migration)
            mig = findFirstMig(locusID_, node, age);
        }

    }

    #ifdef DEBUG_EVENT_GRAPH
            printf("\n");
    #endif

    return 0;
}



/*
    returns true if node is a leaf
    @param: node id
    @return: boolean
*/
bool EventsGraph::isLeaf(int node) {
    if (node < pDataSetup_->numSamples)
        return true;
    return false;
}


/*
    returns samplesStart of a population
    @param: population id
    @return: pointer to a SamplesStart event
*/
EventNode* EventsGraph::getSamplesStartEvent(int pop) {

    //iterate over events and find event of type samplesStart
    EventNode* pNode = this->getStartEvent(pop)->getNextGenEvent();
    while (!pNode->isType(EventType::SAMPLES_START)) {
        pNode->getNextGenEvent();
    }
    return pNode;
}


/*
 * !!!A temp function, should be removed in future versions!!!
   Gets a node id and returns a pointer to the corresponding coalescence event
   @param: node ID
   @return: a pointer to EventNode
*/
EventNode* EventsGraph::getCoalEvent(int nodeID) {
    return coalEvents_[nodeID];
}


/*
 * !!!A temp function, should be removed in future versions!!!
   @param: node ID
   @return: pointer to eventNode
*/
EventNode* EventsGraph::getEventByNode(int nodeID) {

    //if node is a leaf - return a SamplesStart event
    if (this->isLeaf(nodeID)) {
        //get population of node
        int pop = nodePops[locusID_][nodeID];
        return this->getSamplesStartEvent(pop);
    }
    //o.w. return the corresponding coalescent event
    return this->getCoalEvent(nodeID);
}

LocusData* EventsGraph::getLocusData() {
    return pDataState_->lociData[locusID_];
}

