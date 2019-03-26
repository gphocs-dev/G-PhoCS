//
// Created by nomihadar on 3/11/19.
//

#include "LocusEmbeddedGenealogy.h"
#include "DbgErrMsgIntervals.h"
#include <iostream>
#include <iomanip>

/*
	Constructs
	Constructs everything from scratch
	Typically used only for initial genetrees or for testing.
	Records number of lineages only for first events in leaf populations.
	The rest are recorded by computeGenetreeStats
*/
LocusEmbeddedGenealogy::LocusEmbeddedGenealogy(
        int locusID,
        int numIntervals,
        DATA_SETUP* pDataSetup,
        PopulationTree* pPopTree,
        DATA_STATE* pDataState,
        GENETREE_MIGS* pGenetreeMigs)

        : genealogy_(pDataSetup->numSamples, 0), //construct genealogy
          intervals_(locusID, numIntervals, pPopTree),  //construct intervals
          locusID_(locusID),
          pDataSetup_(pDataSetup),
          pPopTree_(pPopTree),
          pDataState_(pDataState),
          pGenetreeMigs_(pGenetreeMigs) {

    //create map between leaf to its pop
    // and map between pop to its leaves
    for (int node = 0; node < pDataSetup_->numSamples; node++) {
        int pop = nodePops[locusID_][node];
        leafToPop_[node] = pop;
        popToLeaves_[pop].emplace_back(node); //todo reserve?
    }

}


/*
   @return: a pointer to locus data of current locus
*/
LocusData* LocusEmbeddedGenealogy::getLocusData() {
    return pDataState_->lociData[locusID_];
}


/*
	Constructs genealogy and intervals
	Genealogy: construct branches and link to corresponding intervals,
                add mig nodes to tree
    Intervals: reset intervals, ling them to each other
                initialize with start and end intervals,
                create samples start intervals,
                create coalescent and migration intervals,
                link intervals to corresponding nodes.

	Typically used only for initial genetrees or for testing.
	Records number of lineages only for first events in leaf populations.
	The rest are recorded by computeGenetreeStats
*/
int LocusEmbeddedGenealogy::construct_genealogy_and_intervals() {

    //construct genealogy branches (edges between tree nodes)
    genealogy_.constructBranches(this->getLocusData());

    //reset intervals
    intervals_.resetIntervals();

    //link intervals to each other
    intervals_.linkIntervals();

    //add start and end intervals
    intervals_.addStartEndIntervals();

    //create samples start intervals (for ancient samples)
    for (int pop = 0; pop < pPopTree_->numCurPops; pop++) {

        //create interval
        double age = pPopTree_->pops[pop]->sampleAge;
        PopInterval* pInterval =
                intervals_.createInterval(pop, age,
                                          IntervalType::SAMPLES_START);

        if (!pInterval) {
            INTERVALS_FATAL_0024
        }
    }

    //create coalescent intervals
    //and link intervals to genealogy and vice versa
    int nSamples = pDataSetup_->numSamples;
    for (int node = 0; node < 2*nSamples-1; node++) {

        //get population and age of node
        int pop = nodePops[locusID_][node];
        double age = getNodeAge(getLocusData(), node);

        // if node is a leaf - link it to its sampleStart interval
        if (genealogy_.isLeaf(node)) {

            //get samples start interval of pop
            PopInterval* pInterval = intervals_.getSamplesStart(pop);

            //get leaf node by current node id
            LeafNode* pNode = genealogy_.getLeafNode(node);

            //samplesStart interval points to leaf node
            pInterval->setTreeNode(pNode);

            //leaf node points to samplesStart interval
            pNode->setSamplesInterval(pInterval);


        } else { //if node is not a leaf create a coal interval and link to node

            //create a coalescent interval
            PopInterval* pInterval =
                    intervals_.createInterval(pop, age, IntervalType::COAL);

            if (!pInterval) {
                INTERVALS_FATAL_0025
            }

            //get coal node by current node id
            CoalNode* pNode = genealogy_.getCoalNode(node);

            //coal interval points to coal node
            pInterval->setTreeNode(pNode);

            //coal node points to coal interval
            pNode->setCoalInterval(pInterval);
        }
    }

    //create migration intervals, add mig nodes to genealogy tree
    //and link between them
    for (int node = 0; node < 2*nSamples-1; node++) {

        //get tree node by current node id
        TreeNode* pTreeNode = genealogy_.getTreeNodeByID(node);

        //find migration above current node and after specified time
        int mig = findFirstMig(locusID_, node,
                               getNodeAge(getLocusData(), node));

        //while there are migration events on the edge above current node
        while (mig != -1) {

            //create migration events (source and target)

            //get age of migration, and target and source populations
            double age = pGenetreeMigs_[locusID_].mignodes[mig].age;
            int target_pop = pGenetreeMigs_[locusID_].mignodes[mig].target_pop;
            int source_pop = pGenetreeMigs_[locusID_].mignodes[mig].source_pop;

            //create an incoming migration interval
            PopInterval* pMigIn =
                    intervals_.createInterval(target_pop, age,
                                              IntervalType::IN_MIG);
            if (!pMigIn) {
                INTERVALS_FATAL_0022
            }

            //create an outgoing migration interval
            PopInterval* pMigOut =
                    intervals_.createInterval(source_pop, age,
                                              IntervalType::OUT_MIG);
            if (!pMigOut) {
                INTERVALS_FATAL_0023
            }

            //add a migration node to genealogy
            MigNode* pMigNode = genealogy_.addMigNode(pTreeNode);

            //mig intervals points to mig node
            pMigIn->setTreeNode(pMigNode);
            pMigOut->setTreeNode(pMigNode);

            //mig node points to incoming and outgoing intervals
            pMigNode->setInMigInterval(pMigIn);
            pMigNode->setOutMigInterval(pMigOut);

            //update tree node
            pTreeNode = pMigNode;

            //find next migration (after time of current migration)
            mig = findFirstMig(locusID_, node, age);
        }

    }

    return 0;
}


void LocusEmbeddedGenealogy::print() {

    //print population tree
    printPopulationTree(this->pDataSetup_->popTree, stderr, 1);

    //print genealogy
    std::cout << "------------------------------------------------------" << endl;
    genealogy_.printGenealogy();

    //print intervals
    std::cout << "------------------------------------------------------" << endl;
    intervals_.printIntervals();

}
