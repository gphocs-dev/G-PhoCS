//
// Created by nomihadar on 3/28/19.
//

#ifndef G_PHOCS_TESTNEWDATASTRUCTURE_H
#define G_PHOCS_TESTNEWDATASTRUCTURE_H



#include "LocusGenealogy.h"
#include "LocusEmbeddedGenealogy.h"
#include "TreeNode.h"
#include "MCMCcontrol.h"
#include "GPhoCS.h"
#include "MemoryMng.h"

#include "assert.h"

#include <iostream>


void testIntervals(LocusEmbeddedGenealogy& locusEmbeddedGenealogy) {

    //get locus id and locus data
    int locusID = locusEmbeddedGenealogy.getLocusID();
    LocusData* locusData = locusEmbeddedGenealogy.getLocusData();

    //get intervals
    LocusPopIntervals& locusIntervals = locusEmbeddedGenealogy.getIntervals();

    //for each pop
    for (int pop = 0; pop < dataSetup.popTree->numPops; pop++) {

        //get first event in old structure
        int event = event_chains[locusID].first_event[pop];

        //get first (pop start) interval in the new structure
        PopInterval* pInterval = locusIntervals.getPopStart(pop)->getNext();

        //iterate both old and new structure
        for (event; event >= 0; event = event_chains[locusID].events[event].getNextIdx()) {

            //get event type
            EventType eventType = event_chains[locusID].events[event].getType();
            double t;
            switch (eventType) {

                case SAMPLES_START: {
                    assert(("type", pInterval->isType(IntervalType::SAMPLES_START)));
                }
                case COAL: {
                    assert(("type", pInterval->isType(IntervalType::COAL)));
                }
                case IN_MIG: {
                    assert(("In-mig type", pInterval->isType(IntervalType::IN_MIG)));
                }
                case OUT_MIG: {
                    assert(("Out-mig type", pInterval->isType(IntervalType::OUT_MIG)));

                    //compare num lineages
                    int eventLin = event_chains[locusID].events[event].getNumLineages();
                    int intervalLin = pInterval->getNumLineages();
                    assert(("nLineages", eventLin == intervalLin));

                    //compare elapsed time
                    double eventTime = event_chains[locusID].events[event].getElapsedTime();
                    double intervalTime = pInterval->getElapsedTime();
                    assert(("time", eventTime == intervalTime));

                    break;
                }

                case MIG_BAND_START: {

                    double eventTime = event_chains[locusID].events[event].getElapsedTime();

                    continue;
                }


                case MIG_BAND_END: {
                    continue;
                }


                case END_CHAIN: {
                    assert(("End differ", pInterval->isType(IntervalType::POP_END)));
                    break;

                }



                default: {

                    std::cout << "other event type: "  << eventType << std::endl;
                }


            }


            pInterval = pInterval->getNext();




            event_chains[locusID].events[event].getElapsedTime();



            event_chains[locusID].events[event].getPrevIdx();
            event_chains[locusID].events[event].getNextIdx();

        }

    }



}

void testLocusGenealogy(int numSamples,
        LocusEmbeddedGenealogy& locusEmbeddedGenealogy,
        GENETREE_MIGS* pGenetreeMigs) {

    //get locus id and locus data
    int locuID = locusEmbeddedGenealogy.getLocusID();
    LocusData* locusData = locusEmbeddedGenealogy.getLocusData();

    //get genealogy
    LocusGenealogy& locusGenealogy = locusEmbeddedGenealogy.getGenealogy();

    //nodePops[];

    //iterate by node id
    for (int node = 0; node < 2*numSamples-1; node++) {

        //get coal or leaf node
        TreeNode* pNode = locusGenealogy.getTreeNodeByID(node);

        //find migration above current node and after specified time
        int mig = findFirstMig(locuID, node, getNodeAge(locusData, node));;

        //while there are migration events on the edge above current node
       while (mig != -1) {

           pNode = pNode->getParent();

           int mig_node_id = pNode->getNodeId();

           //MigNode* pMigNode = locusGenealogy.getMigNode(mig_node_id);

           assert(("Mig ids differ", mig == mig_node_id));

           double age = pGenetreeMigs[locuID].mignodes[mig].age;
           mig = findFirstMig(locuID, node, age);
       }

        int parent = getNodeFather(locusData, node);
        int lSon = getNodeSon(locusData, node, 0);
        int rSon = getNodeSon(locusData, node, 1);

        int parent2 = pNode->getParent() ? pNode->getParent()->getNodeId() : -1;
        int lSon2 = pNode->getLeftSon() ? pNode->getLeftSon()->getNodeId() : -1;
        int rSon2 = pNode->getRightSon() ? pNode->getRightSon()->getNodeId() : -1;

        assert(("Parents differ", parent == parent2));
        assert(("Left sons differ", lSon == lSon2));
        assert(("Right sons differ", rSon == rSon2));

    }

}


#endif //G_PHOCS_TESTNEWDATASTRUCTURE_H
