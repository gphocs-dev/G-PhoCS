
#ifndef G_PHOCS_DAG_H
#define G_PHOCS_DAG_H

#include "EventNode.h"
#include "MCMCcontrol.h"
#include "GPhoCS.h"
#include "MemoryMng.h"

/*=============================================================================
 *
 * EventsGraph class
 *
 * EventsGraph is a graph whose nodes are objects of EventNode class.
 * The graph is implemented as array of pre-allocated EventNodes.
 *
 * Contains:
 * 1. Genealogy ID (which is locus id)
 * 2. Array of NodeEvent objects.
 * 3. A pointer to a pool of free (unoccupied) nodes in graph.
 *    Free nodes are linked to each other, and pointer points to head.
 * 4. Total number of nodes (events) in graph.
 * 5. Pointers to global structs.
 *   (DataSetup, PopulationTree, DATA_STATE, and GENETREE_MIGS structs)
 *
 *===========================================================================*/
class EventsGraph {

private:

    int genealogyID_;            //id of genealogy (locus id)

    EventNode* eventNodes_;     //array of EventNodes
    EventNode* pNodesPool_;     //pointer to a pool of free nodes

    int nNodes_;                //total number of nodes in graph

    DATA_SETUP* pDataSetup_;        //pointer to DATA_SETUP struct
    PopulationTree* pPopTree_;      //pointer to PopulationTree struct
    DATA_STATE* pDataState_;        //pointer to DATA_STATE struct
    GENETREE_MIGS* pGenetreeMigs_;  //pointer to GENETREE_MIGS struct

public:

    //constructor
    EventsGraph(int genID, int nNodes,
                DATA_SETUP* pSetup,
                DATA_STATE* pState,
                GENETREE_MIGS* pGenetreeMigs);

    //destructor
    ~EventsGraph();

    //initialize EventsGraph with start and end events
    void initializeGraph();

    //get a free node from nodes pool
    EventNode* getNodeFromPool();

    //add a free node to nodes pool
    void addNodeToPool(EventNode* pNode);

    //get the strat event of a specified population
    EventNode* getStartEvent(int pop);

    //create a new event before a given event
    EventNode* createEventBefore(
            EventNode* pNode, double elapsed_time, EventType type);

    //create a new event in specified population at given time.
    EventNode* createEvent(int pop, double age,
            EventType type = EventType::DUMMY);

    //construct events graph
    int constructEventGraph();


};




#endif //G_PHOCS_DAG_H
