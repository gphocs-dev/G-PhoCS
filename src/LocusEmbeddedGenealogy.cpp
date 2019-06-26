
#include "LocusEmbeddedGenealogy.h"
#include "DbgErrMsgIntervals.h"
#include "PopulationTree.h"

#include <iostream>
#include <iomanip>
#include <cmath>


/*
 * LocusEmbeddedGenealogy / constructor
*/
LocusEmbeddedGenealogy::LocusEmbeddedGenealogy(
        int locusID, int numIntervals,
        DATA_SETUP* pDataSetup,
        DATA_STATE* pDataState,
        GENETREE_MIGS* pGenetreeMigs)

        : genealogy_(pDataSetup->numSamples), //construct genealogy
          intervals_(locusID, numIntervals),  //construct intervals
          locusID_(locusID),
          pDataSetup_(pDataSetup),
          pDataState_(pDataState),
          pGenetreeMigs_(pGenetreeMigs) {

}


/*
 * LocusEmbeddedGenealogy / copy-constructor
*/
LocusEmbeddedGenealogy::LocusEmbeddedGenealogy(
        const LocusEmbeddedGenealogy& other) :

        genealogy_(other.genealogy_), //copy-construct genealogy
        intervals_(other.intervals_), //copy-construct intervals
        locusID_(other.locusID_),
        pDataSetup_(other.pDataSetup_),
        pDataState_(other.pDataState_),
        pGenetreeMigs_(other.pGenetreeMigs_) {

    //copy pointers linking between genealogy and intervals

    // ******* set pointers of genealogy -> intervals *******

    //leaf nodes
    for (int i = 0; i < pDataSetup_->numSamples; i++) {
        PopInterval* pOri = other.genealogy_.getLeafNode(i)->getSamplesStart();
        PopInterval* pCopy =  intervals_.getNewPos(other.intervals_, pOri);
        genealogy_.getLeafNode(i)->setSamplesInterval(pCopy);
    }

    //coal nodes
    for (int i = pDataSetup_->numSamples; i < 2*pDataSetup_->numSamples-1; i++) {
        PopInterval* pOri = other.genealogy_.getCoalNode(i)->getCoalInterval();
        PopInterval* pCopy =  intervals_.getNewPos(other.intervals_, pOri);
        genealogy_.getCoalNode(i)->setCoalInterval(pCopy);
    }

    //mig nodes
    for (int i = 0; i < other.genealogy_.getNumMigs(); i++) {
        PopInterval* pOri = other.genealogy_.getMigNode(i)->getInMigInterval();
        PopInterval* pOri2 = other.genealogy_.getMigNode(i)->getOutMigInterval();

        PopInterval* pCopy =  intervals_.getNewPos(other.intervals_, pOri);
        PopInterval* pCopy2 =  intervals_.getNewPos(other.intervals_, pOri2);

        genealogy_.getMigNode(i)->setInMigInterval(pCopy);
        genealogy_.getMigNode(i)->setOutMigInterval(pCopy2);
    }

    // ******* set pointers of intervals -> genealogy  *******

    for (int i = 0; i < intervals_.getNumIntervals(); i++) {
        TreeNode* pOri = other.intervals_.getInterval(i)->getTreeNode();
        if (pOri) {
            TreeNode* pCopy = genealogy_.getNewPos(other.genealogy_, pOri);
            intervals_.getInterval(i)->setTreeNode(pCopy);
        }
    }

}


/*
 * getLocusData
 * @return: a pointer to locus data of current locus
*/
LocusData* LocusEmbeddedGenealogy::getLocusData() {
    return pDataState_->lociData[locusID_];
}


/*
	construct_genealogy_and_intervals
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

    //reset genealogy
    genealogy_.resetGenealogy();

    //construct genealogy branches (edges between tree nodes)
    genealogy_.constructBranches(this->getLocusData());

    //reset intervals
    intervals_.resetPopIntervals();

    //link intervals to each other
    intervals_.linkIntervals();

    //add start and end intervals
    intervals_.createStartEndIntervals();

    //create samples start intervals (for ancient samples)
    for (int pop = 0; pop < pDataSetup_->popTree->numCurPops; pop++) {

        //create interval
        double age = pDataSetup_->popTree->pops[pop]->sampleAge;
        PopInterval* pInterval =
                intervals_.createInterval(pop, age, IntervalType::SAMPLES_START);
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

            //leaf node points to samplesStart interval (but samplesStart
            // interval points to null since there are several leaves)
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
    for (int node = 0; node < 2 * nSamples - 1; node++) {

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

            //get migration band ID
            int bandId = getMigBandByPops(pDataSetup_->popTree, source_pop,
                                          target_pop)->id;

            //add a migration node to genealogy
            MigNode* pMigNode = genealogy_.addMigNode(pTreeNode, bandId);

            //set age
            pMigNode->setAge(age);

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


/*
 * computeGenetreeStats
 * compute genealogy tree statistics
*/

int LocusEmbeddedGenealogy::computeGenetreeStats() {
    return intervals_.computeGenetreeStats();
}


/*
 * recalcStats
 * recalculate statistics
*/
double LocusEmbeddedGenealogy::recalcStats(int pop) {
  return intervals_.recalcStats(pop);
}


/*
 * printEmbeddedGenealogy
 * print population tree, genealogy and intervals
*/
void LocusEmbeddedGenealogy::printEmbeddedGenealogy() {

    //print population tree
    //printPopulationTree(pDataSetup_->popTree, stderr, 1);

    //print genealogy
    std::cout << "------------------------------------------------------"
              << std::endl;
    genealogy_.printGenealogy();

    //print intervals
    std::cout << "------------------------------------------------------"
              << std::endl;
    intervals_.printIntervals();

}


/*
 * getLocusID
*/
int LocusEmbeddedGenealogy::getLocusID() {
    return locusID_;
}



/*
 * getStats
 * get a reference to statistics
*/
const GenealogyStats& LocusEmbeddedGenealogy::getStats() const {
    return intervals_.getStats();
}


/*
 * testLocusGenealogy
 * verify new genealogy data structure is consistent with the original
*/
void LocusEmbeddedGenealogy::testLocusGenealogy() {
    genealogy_.testLocusGenealogy(locusID_, this->getLocusData(),
                                  pGenetreeMigs_);
}


/*
 * testPopIntervals
 * verify new events data structure is consistent with the original
*/
void LocusEmbeddedGenealogy::testPopIntervals() {
    intervals_.testPopIntervals();
}


/*
 * testTotalStats
 * verify statistics are equal to statistics of old data structurel
*/
void LocusEmbeddedGenealogy::testGenealogyStats() {
    intervals_.testGenealogyStatistics();
}




