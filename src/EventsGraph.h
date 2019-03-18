
#ifndef G_PHOCS_DAG_H
#define G_PHOCS_DAG_H

#include "EventNode.h"
#include "MCMCcontrol.h"
#include "GPhoCS.h"
#include "MemoryMng.h"
#include <vector>

/*=============================================================================
 *
 * EventsGraph class
 *
 * EventsGraph is a graph whose nodes are objects of EventNode class.
 * The graph is implemented as array of pre-allocated EventNodes.
 *
 * Contains:
 * 1. Locus ID (= genealogy id).
 * 2. Total number of EventsNodes in graph.
 * 3. Array of NodeEvent objects. This is actually the graph.
 * 5. A pointer to a pool of free (unoccupied) nodes in graph.
 *    Free nodes are linked to each other, and pointer points to head.
 * 6. Vector of pointers to EventNodes of type leaf.
 * 7. Vector of pointers to EventNodes of type coalescent.
 * 8. Vector of pointers to EventNodes of type migration.
 * 9. Pointers to global structs.
 *   (DataSetup, PopulationTree, DATA_STATE, and GENETREE_MIGS structs)
 *
 *===========================================================================*/
class EventsGraph {

private:

    int locusID_;   //id of locus (genealogy id)

    int nNodes_;    //total number of nodes in graph

    EventNode* eventsGraph_;     //array of EventNodes
    EventNode* pNodesPool_;     //pointer to a pool of free nodes

    std::vector<EventNode*> leafEvents_; //vector of pointers to leaves EventNodes
    std::vector<EventNode*> coalEvents_; //vector of pointers to coalescent EventNodes
    std::vector<EventNode*> migEvents_; //vector of pointers to migration EventNodes

    DATA_SETUP*     pDataSetup_;       //pointer to DATA_SETUP struct
    PopulationTree* pPopTree_;         //pointer to PopulationTree struct
    DATA_STATE*     pDataState_;       //pointer to DATA_STATE struct
    GENETREE_MIGS*  pGenetreeMigs_;    //pointer to GENETREE_MIGS struct

public:

    //constructor
    EventsGraph(int genID, int nNodes,
                DATA_SETUP* pSetup,
                DATA_STATE* pState,
                GENETREE_MIGS* pGenetreeMigs);

    //destructor
    ~EventsGraph();

public:

    //initialize EventsGraph with start and end events
    void initialiseGraph();

    //print events graph
    void printEventsGraph();

//*****************************************************************************
    //get a free node from nodes pool
    EventNode* getNodeFromPool();

    //add a free node to nodes pool
    void addNodeToPool(EventNode* pNode);

//*****************************************************************************
    //get a pointer to pop start event of a given population
    EventNode* getStartEvent(int population);

    //get a pointer to sampleStart of a given population
    EventNode* getSamplesStartEvent(int population);

//*****************************************************************************
    //create a new event before a given event
    EventNode* createEventBefore(
            EventNode* pNode, int pop, double elapsed_time, EventType type);

    //create a new event in specified population at given time.
    EventNode* createEvent(int pop, double age,
            EventType type = EventType::DUMMY);

//*****************************************************************************


    //construct events graph
    int constructEventsGraph();


//*****************************************************************************
//temp functions which intermediate between node (in the sense of the old
// version, where node is an integer) and eventNodes

    //return true if a node is a leaf
    bool isLeaf(int node);

    //return a pointer to a leaf event by node id
    EventNode* getLeafEvent(int nodeID);

    //return a pointer to a coalescence event by node id
    EventNode* getCoalEvent(int nodeID);

    //return a pointer to event by node id
    EventNode* getEventByNode(int nodeID);

//*****************************************************************************

    LocusData* getLocusData();


//*****************************************************************************

};




#endif //G_PHOCS_DAG_H
