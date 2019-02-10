
#ifndef G_PHOCS_DAG_H
#define G_PHOCS_DAG_H

//#include "MCMCcontrol.h"
#include "EventNode.h"

/*=============================================================================
 *
 * EventsGraph class
 *
 * EventsGraph is a graph whose nodes are objects of EventNode class.
 * The graph is implemented as array of pre-allocated EventNodes.
 *
 * Contains:
 * 1. Array of NodeEvent objects.
 * 2. A pointer to a pool of free (unoccupied) nodes in graph.
 *    Free nodes are linked to each other, and pointer points to head.
 * 3. Total number of nodes (events) in graph.
 * 4. A pointer to the DataSetup struct.
 *
 *===========================================================================*/
class EventsGraph {

private:

    EventNode* eventNodes_;     //array of EventNodes
    EventNode* pNodesPool_;    //pointer to a pool of free nodes
    int nNodes_;               //total number of nodes in graph
    DATA_SETUP* pDataSetup_;   //pointer to dataSetup struct

public:

    //constructor
    EventsGraph(DATA_SETUP* pDataSetup, int nNodes);

    //destructor
    ~EventsGraph();

    //initialize EventsGraph with start and end events
    void initializeDAG();

    //get a free node from nodes pool
    EventNode* getNodeFromPool();

    //add a free node to nodes pool
    void addNodeToPool(EventNode* pEventNode);

    //create event before given event
    void createEventBefore(EventNode* pNode, double elapsedTime);

};




#endif //G_PHOCS_DAG_H
