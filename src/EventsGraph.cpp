
#include "EventsGraph.h"
#include "DbgErrMsgEvents.h"
#include "DataLayerConstants.h"
#include <iostream>
#include <iomanip>

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

    //reserve place in vector
    //number of leaves is num_samples
    //number of coalescent events is num_samples-1,
    leafEvents_.reserve(pDataSetup_->numSamples);
    coalEvents_.reserve(pDataSetup_->numSamples - 1);
    
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




void printNumber(int x) {
    std::cout << "X:" << std::setw(6) << x << ":X\n";
}

void printStuff() {
    printNumber(528);
    printNumber(3);
    printNumber(73826);
    printNumber(37);
}


/*
    Prints EventsGraph

*/
void EventsGraph::printEventsGraph() {

    //print population tree -
    printPopulationTree(this->pDataSetup_->popTree, stderr, 1);


    printLocusGenTree(dataState.lociData[gen], stderr, nodePops[gen],nodeEvents[gen]);
    //printEventChains(stderr, gen);

    //print genealogy tree
    std::cout << "Genalogy tree:" << std::endl;;

    int nSamples = pDataSetup_->numSamples;//TODO: replace with 2*locusData->numLeaves - 1;
    for (int node = 0; node < 2*nSamples-1; node++) {
        //std::cout.precision(4);

        //std::cout << std::left;


    }

    int node, numNodes = 2*locusData->numLeaves - 1;

    fprintf(stream, "Genalogy tree:\n");
    for(node=0; node<numNodes; node++) {
        fprintf(stream, "Node %2d, age [%.10f], father (%2d), sons (%2d %2d), pop (%2d), event-id (%2d)",
                node, locusData->nodeArray[node]->age, locusData->nodeArray[node]->father,
                locusData->nodeArray[node]->leftSon, locusData->nodeArray[node]->rightSon,
                nodePops[node], nodeEvents[node]);
        if(locusData->root == node)			fprintf(stream, " - Root\n");
        else if(node<locusData->numLeaves)	fprintf(stream, " - Leaf\n");
        else								fprintf(stream, "\n");
    }

    fprintf(stream, "---------------------------------------------------------------\n");
    if(locusData->savedVersion.numChangedNodes > 0) {
        fprintf(stream, "There are %d changed nodes:",locusData->savedVersion.numChangedNodes);
        for(node=0; node<locusData->savedVersion.numChangedNodes; node++) {
            fprintf(stream, " %d",locusData->savedVersion.changedNodeIds[node]);
        }
        fprintf(stream, "\n---------------------------------------------------------------\n");
    }

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
    Creates an event of type before a specified event
    @param: pointer to EventNode, elapsed time of new event,
            and (optional) event type
    @return: pointer to the new node

*/
EventNode* EventsGraph::createEventBefore(
        EventNode* pNode, int pop, double elapsed_time, EventType type) {

    //get a new event from pool (default type is DUMMY)
    EventNode* pNewNode = this->getNodeFromPool();

    //set population ID
    pNewNode->setPopID(pop);

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
    return this->createEventBefore(pop, pNode, delta_time, type);

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
    /*
    for (int pop = 0; pop < pPopTree_->numCurPops; pop++) {

        //create event
        EventNode* pSamplesStart =
                this->createEvent(pop,
                                  pPopTree_->pops[pop]->sampleAge,
                                  EventType::SAMPLES_START);

        //coalescent pointers of samplesStart are null
        // (this is the default, no need to code)
    }*/

    //create coalescent and leaf events
    int nSamples = pDataSetup_->numSamples;
    for (int node = 0; node < 2*nSamples-1; node++) {

        //get population of node
        int pop = nodePops[locusID_][node];

        //get age of node
        int age = getNodeAge(getLocusData(), node);

        // if node is a leaf
        if (this->isLeaf(node)) {

            #ifdef DEBUG_EVENT_GRAPH
                        printf("L%d ",pop);
            #endif

            //create a coalescent event
            EventNode* pCoal = this->createEvent(pop, age, EventType::LEAF);
            if (!pCoal) {
                EVENTS_FATAL_0024
            }

            //save pointer in leaf events array
            leafEvents_[node] = pCoal;

        } else {

            #ifdef DEBUG_EVENT_GRAPH
                        printf("C%d %d %f ", pop, node, age);
            #endif

            //create a coalescent event
            EventNode* pCoal = this->createEvent(pop, age, EventType::COAL);
            if (!pCoal) {
                EVENTS_FATAL_0025
            }

            //save pointer in coalescent events array
            // (with offset of num samples-1)
            coalEvents_[node - pDataSetup_->numSamples-1] = pCoal;
        }
    }

    //set the coalescent pointers of the coalescent events
    for (int node = 0; node < 2*nSamples-1; node++) {

        //get eventNode of current node
        EventNode* pCoal = this->getEventByNode(node);

        //set parent pointer
        //get node and corresponding eventNode of the father
        //if current node is the root - parent points to null
        if ( node != getLocusRoot(this->getLocusData()) ) {

            int nodeFather = getNodeFather(this->getLocusData(), node);
            EventNode* pFather = this->getCoalEvent(nodeFather);
            pCoal->setCoalParent(pFather);
        }
        else {
            pCoal->setCoalParent(nullptr);
        }

        //set sons pointers
        //get nodes and corresponding eventNodes of the sons
        //eventNode can be of type coalescent or leaf
        //if current node is a leaf - sons point to null
        if (!this->isLeaf(node)) {

            int nodeLeftSon = getNodeSon(this->getLocusData(), node, 0);
            EventNode* pLeftSon = this->getEventByNode(nodeLeftSon);
            pCoal->setCoalLeftSon(pLeftSon);

            int nodeRightSon = getNodeSon(this->getLocusData(), node, 1);
            EventNode* pRightSon = this->getEventByNode(nodeRightSon);
            pCoal->setCoalRightSon(pRightSon);

        }
        else {
            pCoal->setCoalSons(nullptr);
        }
    }

    //create migration events and if needed reset coalescent pointers
    // of coalescent events
    for (int node = 0; node < 2*nSamples-1; node++) {

        //get event node of current node
        EventNode* pEventNode = this->getEventByNode(node);

        //define an iterator
        EventNode* pIterator = pEventNode;


        //find migration above current node and after specified time
        int mig = findFirstMig(locusID_, node, getNodeAge(getLocusData(), node));


        bool first_mig =  true;

        //while there are migration events on the edge above current node
        while (mig != -1) {

            //create migration events (source and target)

            //get age of migration, and target and source populations
            double age = pGenetreeMigs_[locusID_].mignodes[mig].age;
            int target_pop = pGenetreeMigs_[locusID_].mignodes[mig].target_pop;
            int source_pop = pGenetreeMigs_[locusID_].mignodes[mig].source_pop;

            //create an incoming migration event
            EventNode* pMigIn =
                    this->createEvent(target_pop, age, EventType::IN_MIG);
            if (!pMigIn) {
                EVENTS_FATAL_0022
            }

            //create an outgoing migration event
            EventNode* pMigOut =
                    this->createEvent(source_pop, age, EventType::OUT_MIG);
            if (!pMigOut) {
                EVENTS_FATAL_0023
            }

            //add pointer to incoming migration
            migEvents_.push_back(pMigIn);

            //parent of incoming migration points to outgoing migration
            //sons of incoming migration point to bottom node
            pMigIn->setCoalParent(pMigOut);
            pMigIn->setCoalSons(pIterator);

            //parent of outgoing migration points to upper node
            //sons of outgoing migration point to outgoing migration
            pMigOut->setCoalParent(pEventNode->getCoalParent());
            pMigOut->setCoalSons(pMigIn);

            //set parent of current node
            if (first_mig) {
                pEventNode->setCoalParent(pMigIn);
                first_mig = false;
            }

            //set one of the sons of current node's parent -
            //depends whether current node is a left or a right son
            if (pEventNode == pEventNode->getCoalParent()->getCoalLeftSon())
                pEventNode->getCoalParent()->setCoalLeftSon(pMigOut);
            else
                pEventNode->getCoalParent()->setCoalRightSon(pMigOut);

            //promote bottom node to current outgoing migration
            pIterator = pMigOut;

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
    while (!pNode->isType(EventType::POP_END)) {
        if (pNode->isType(EventType::LEAF))
            return pNode;
        pNode = pNode->getNextGenEvent();
    }
    return nullptr;
}


/*
 * !!!A temp function, should be removed in future versions!!!
   Gets a node id and returns a pointer to the corresponding coalescence event
   @param: node ID
   @return: a pointer to EventNode
*/
EventNode* EventsGraph::getCoalEvent(int nodeID) {
    int offset = nodeID - pDataSetup_->numSamples -1;
    return coalEvents_[offset];
}



/*
 * !!!A temp function, should be removed in future versions!!!
   Gets a node id and returns a pointer to the corresponding leaf event
   @param: node ID
   @return: a pointer to EventNode
*/
EventNode* EventsGraph::getLeafEvent(int nodeID) {
    return leafEvents_[nodeID];
}


/*
 * !!!A temp function, should be removed in future versions!!!
   @param: node ID
   @return: a pointer to eventNode
*/
EventNode* EventsGraph::getEventByNode(int nodeID) {

    //if node is a leaf - return a leaf event
    //o.w. return the coalescent event
    if (this->isLeaf(nodeID))
        return this->getLeafEvent(nodeID);
    else
        return this->getCoalEvent(nodeID);
}


/*
 *
   @return: a pointer to locus data of current locus
*/
LocusData* EventsGraph::getLocusData() {
    return pDataState_->lociData[locusID_];
}



