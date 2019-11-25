
#include "LocusEmbeddedGenealogy.h"
#include "DbgErrMsgIntervals.h"


/*
 * LocusEmbeddedGenealogy / constructor
*/
LocusEmbeddedGenealogy::LocusEmbeddedGenealogy(
        int locusID, int numIntervals,
        DATA_SETUP *pSetup,
        DATA_STATE *pState,
        GENETREE_MIGS *pGenetreeMigs)

        : locusID_(locusID),

          genealogy_(pSetup->numSamples,
                     pState->lociData[locusID]), //construct genealogy
          intervalsPro_(locusID, numIntervals),  //construct proposal intervals
          intervalsOri_(locusID, numIntervals),  //construct original intervals

          //assign pointers
          pSetup_(pSetup),
          pState_(pState),
          pGenetreeMigs_(pGenetreeMigs),

          //initialize likelihoods
          genLogLikelihood_(0), dataLogLikelihood_(0) {

}


/*
 * copy without construction
*/
void LocusEmbeddedGenealogy::copy(const LocusEmbeddedGenealogy &other) {

    genealogy_.copy(other.genealogy_); //copy genealogy
    intervalsPro_.copy(other.intervalsPro_); //copy intervals

    //copy pointers linking between genealogy and intervals

    // ******* set pointers of genealogy -> intervals *******

    PopInterval *pOri, *pCopy;

    //leaf nodes
    for (int i = 0; i < pSetup_->numSamples; i++) {
        pOri = other.genealogy_.getLeafNode(i)->getInterval();
        pCopy = intervalsPro_.getNewPos(other.intervalsPro_, pOri);
        genealogy_.getLeafNode(i)->setInterval(pCopy);
    }

    //coal nodes
    for (int i = pSetup_->numSamples; i < 2 * pSetup_->numSamples - 1; i++) {
        pOri = other.genealogy_.getCoalNode(i)->getInterval();
        pCopy = intervalsPro_.getNewPos(other.intervalsPro_, pOri);
        genealogy_.getCoalNode(i)->setInterval(pCopy);
    }

    //mig nodes
    for (int i = 0; i < other.genealogy_.getNumMigs(); i++) {
        //set in/out migration
        for (int iInterval = 0; iInterval < 2; iInterval++) {
            pOri = other.genealogy_.getMigNode(i)->getInterval(iInterval);
            pCopy = intervalsPro_.getNewPos(other.intervalsPro_, pOri);
            genealogy_.getMigNode(i)->setInterval(pCopy, iInterval);
        }
    }

    // ******* set pointers of intervals -> genealogy  *******

    for (int i = 0; i < intervalsPro_.getNumIntervals(); i++) {
        TreeNode *pOri = other.intervalsPro_.getInterval(i)->getTreeNode();
        if (pOri) {
            TreeNode *pCopy = genealogy_.getNewPos(other.genealogy_, pOri);
            intervalsPro_.getInterval(i)->setTreeNode(pCopy);
        }
    }

}


/*
 * copyIntervals
*/
void
LocusEmbeddedGenealogy::copyIntervals(bool accepted) {

    //if accepted is true, copy intervals from proposal to original
    //o.w., copy from original to proposal
    LocusPopIntervals & source = accepted? intervalsPro_ : intervalsOri_;
    LocusPopIntervals & copy = accepted? intervalsOri_ : intervalsPro_;

    //tree nodes (and only them) are copied by a shallow copy
    copy.copyIntervals(source, false);
}


/*
 * getLocusData
 * @return: a pointer to locus data of current locus
*/
LocusData *LocusEmbeddedGenealogy::getLocusData() {
    return genealogy_.getLocusData();
}


/*
 * getLogLikelihood
 * @return: locus gen log-likelihood
*/
double LocusEmbeddedGenealogy::getGenLogLikelihood() const {
    return genLogLikelihood_;
}


/*
 * getDataLogLikelihood
 * @return: locus data log-likelihood
*/
double LocusEmbeddedGenealogy::getDataLogLikelihood() const {
    return dataLogLikelihood_;
}


/*
	constructEmbeddedGenealogy
	Genealogy: construct branches and link to corresponding intervals,
                add mig nodes to tree
    Intervals: reset intervals, link them to each other
                initialize with start and end intervals,
                create samples start intervals,
                create coalescent and migration intervals,
                link intervals to corresponding nodes.

	Typically used only for initial genetrees or for testing.
	Records number of lineages only for first events in leaf populations.
	The rest are recorded by computeGenetreeStats
*/
int LocusEmbeddedGenealogy::constructEmbeddedGenealogy() {

    //reset genealogy
    genealogy_.resetGenealogy();

    //construct genealogy branches (edges between tree nodes)
    genealogy_.constructBranches();

    //reset intervals
    intervalsPro_.resetPopIntervals();

    //link intervals to each other
    intervalsPro_.linkIntervals();

    //add start and end intervals
    intervalsPro_.createStartEndIntervals();

    //create samples start intervals (for ancient samples)
    for (int pop = 0; pop < pSetup_->popTree->numCurPops; pop++) {

        //create interval
        double age = pSetup_->popTree->pops[pop]->sampleAge;
        PopInterval *pInterval =
                intervalsPro_.createInterval(pop, age,
                                          IntervalType::SAMPLES_START);
        if (!pInterval) {
            INTERVALS_FATAL_0024
        }
    }

    //create coalescent intervals
    //and link intervals to genealogy and vice versa
    int nSamples = pSetup_->numSamples;
    for (int node = 0; node < 2 * nSamples - 1; node++) {

        //get population and age of node
        int pop = nodePops[locusID_][node];
        double age = genealogy_.getNodeAgeWrap(node);

        // if node is a leaf - link it to its sampleStart interval
        if (genealogy_.isLeaf(node)) {

            //get samples start interval of pop
            PopInterval *pInterval = intervalsPro_.getSamplesStart(pop);

            //get leaf node by current node id
            LeafNode *pNode = genealogy_.getLeafNode(node);

            //leaf node points to samplesStart interval (but samplesStart
            // interval points to null since there are several leaves)
            pNode->setInterval(pInterval);


        } else { //if node is not a leaf create a coal interval and link to node

            //create a coalescent interval
            PopInterval *pInterval =
                    intervalsPro_.createInterval(pop, age, IntervalType::COAL);

            if (!pInterval) {
                INTERVALS_FATAL_0025
            }

            //get coal node by current node id
            CoalNode *pNode = genealogy_.getCoalNode(node);

            //coal interval points to coal node
            pInterval->setTreeNode(pNode);

            //coal node points to coal interval
            pNode->setInterval(pInterval);
        }
    }

    //create migration intervals, add mig nodes to genealogy tree
    //and link between them
    for (int node = 0; node < 2 * nSamples - 1; node++) {

        //get tree node by current node id
        TreeNode *pTreeNode = genealogy_.getTreeNodeByID(node);

        //find migration above current node and after specified time
        int mig = findFirstMig(locusID_, node,
                               genealogy_.getNodeAgeWrap(node));

        //while there are migration events on the edge above current node
        while (mig != -1) {

            //create migration events (source and target)

            //get age of migration, and target and source populations
            double age = pGenetreeMigs_[locusID_].mignodes[mig].age;
            int target_pop = pGenetreeMigs_[locusID_].mignodes[mig].target_pop;
            int source_pop = pGenetreeMigs_[locusID_].mignodes[mig].source_pop;

            //create an incoming migration interval
            PopInterval *pMigIn =
                    intervalsPro_.createInterval(target_pop, age,
                                              IntervalType::IN_MIG);
            if (!pMigIn) {
                INTERVALS_FATAL_0022
            }

            //create an outgoing migration interval
            PopInterval *pMigOut =
                    intervalsPro_.createInterval(source_pop, age,
                                              IntervalType::OUT_MIG);
            if (!pMigOut) {
                INTERVALS_FATAL_0023
            }

            //get migration band ID
            int bandId = getMigBandByPops(pSetup_->popTree, source_pop,
                                          target_pop)->id;

            //add a migration node to genealogy
            MigNode *pMigNode = genealogy_.addMigNode(pTreeNode, bandId);

            //set age
            pMigNode->setAge(age);

            //mig intervals points to mig node
            pMigIn->setTreeNode(pMigNode);
            pMigOut->setTreeNode(pMigNode);

            //mig node points to incoming and outgoing intervals
            pMigNode->setInterval(pMigIn, 0);
            pMigNode->setInterval(pMigOut, 1);

            //update tree node
            pTreeNode = pMigNode;

            //find next migration (after time of current migration)
            mig = findFirstMig(locusID_, node, age);
        }
    }

    //compute data likelihood
    dataLogLikelihood_ = genealogy_.getLocusDataLikelihoodWrap();

    return 0;
}


/*
 * computeGenetreeStats
 * compute genealogy tree statistics
*/

int LocusEmbeddedGenealogy::computeGenetreeStats() {
    return intervalsPro_.computeGenetreeStats();
}


/*
 * recalcStats
 * recalculate statistics
*/
void LocusEmbeddedGenealogy::recalcStats(int pop) {
    intervalsPro_.recalcStats(pop);
}


/*
 * printEmbeddedGenealogy
 * print population tree, genealogy and intervals
*/
void LocusEmbeddedGenealogy::printEmbeddedGenealogy() {

    //print population tree
    printPopulationTree(pSetup_->popTree, stderr, 1);

    //print genealogy
    std::cout << "------------------------------------------------------"
              << std::endl;
    genealogy_.printGenealogy();

    //print intervals
    std::cout << "------------------------------------------------------"
              << std::endl;
    intervalsPro_.printIntervals();

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
const GenealogyStats &LocusEmbeddedGenealogy::getStats() const {
    return intervalsPro_.getStats();
}


/*
 * testLocusGenealogy
 * verify new genealogy data structure is consistent with the original
*/
void LocusEmbeddedGenealogy::testLocusGenealogy() {
    genealogy_.testLocusGenealogy(locusID_, pGenetreeMigs_);
}


/*
 * testPopIntervals
 * verify new events data structure is consistent with the original
*/
void LocusEmbeddedGenealogy::testPopIntervals() {
    intervalsPro_.testPopIntervals();
}


/*
 * testTotalStats
 * verify statistics are equal to statistics of old data structure
*/
void LocusEmbeddedGenealogy::testGenealogyStats() {
    intervalsPro_.testGenealogyStatistics();
}


/*
 * testLogLikelihood
 * verify genealogy and data likelihood values are consistent with original
*/
void LocusEmbeddedGenealogy::testLogLikelihood() {

    assert(fabs(dataLogLikelihood_ - genealogy_.getLocusDataLikelihoodWrap()) <
           EPSILON);
    //todo: replace global locus_data
    assert(fabs(genLogLikelihood_ - locus_data[locusID_].genLogLikelihood) <
           EPSILON);



}


/*
 * testLocusEmbeddedGenealogy
 * verify
*/
void LocusEmbeddedGenealogy::testLocusEmbeddedGenealogy() {

    this->testLocusGenealogy();
    this->testPopIntervals();
    this->testGenealogyStats();
    this->testLogLikelihood();
}


/*
 *	UpdateGB_InternalNode
 *	- perturbs times of all coalescent nodes
 *	- does not change the population of the node
 *	- upper and lower bounds for new time are determined by nodes
 *	  (migration/coalescent)
 *		directly above or below that node, as well as population boundaries.
 *	- records new data log likelihood in *pointerToLnLd
 */
int LocusEmbeddedGenealogy::updateGB_InternalNode(double finetune) {

    double epsilon = 1e-15;

    if (finetune <= 0.0)
        return 0;

    //counter of num accepted of change proposals
    int accepted = 0;

    //save current values of likelihood
    double genLogLd = genLogLikelihood_;
    double dataLogLd = dataLogLikelihood_;

    //for each coal node
    int nSamples = pSetup_->numSamples;
    for (int inode = nSamples; inode < 2 * nSamples - 1; inode++) {

        //get coal node
        CoalNode *pNode = genealogy_.getCoalNode(inode);

        //get its pop
        int pop = pNode->getPop();

        //get lower time bound for new age
        double lowerBound = pSetup_->popTree->pops[pop]->age;
        lowerBound = max2(lowerBound, pNode->getLeftSon()->getAge());
        lowerBound = max2(lowerBound, pNode->getRightSon()->getAge());

        //get upper time bound for new age
        double upperBound = pop != pSetup_->popTree->rootPop
                       ? pSetup_->popTree->pops[pop]->father->age : OLDAGE;
        if (inode != genealogy_.getLocusRootWrap()) //if node is not root node
            upperBound = min2(upperBound, pNode->getParent()->getAge());

        //get current age
        double t = pNode->getAge();

        //calculate new age
        double tnew = t + finetune * rnd2normal8(locusID_);
        tnew = reflect(tnew, lowerBound, upperBound);

        //continue if difference is not significant
        if (fabs(tnew - t) < epsilon) {
            accepted++;
            continue;
        }

        //consider interval move
        double lnAcceptance = this->considerIntervalMove(pNode, tnew);

        //if proposal is accepted
        bool isAccepted =
                lnAcceptance >= 0 or rndu(locusID_) < exp(lnAcceptance);
        if (isAccepted) {

            //increase counter
            accepted++;

            //copy proposal into original
            intervalsOri_.copyIntervals(intervalsPro_, false);

            genLogLikelihood_ = intervalsPro_.computeLogLikelihood();
            dataLogLikelihood_ = genealogy_.getLocusDataLikelihoodWrap();

            resetSaved(this->getLocusData());

        } else { // reject changes and revert to saved version

            //set back node age
            pNode->setAge(t);

            //copy original into proposal
            intervalsPro_.copyIntervals(intervalsOri_, true);

            revertToSaved(this->getLocusData());
        }

        //Below code, untill END, is just for testing this function
        // should be removed later

        //revert changes
        genealogy_.adjustGenNodeAgeWrap(inode, t);
        genealogy_.computeLocusDataLikelihoodWrap(1);
        resetSaved(this->getLocusData());

        //call test function
        updateGB_InternalNode_oldDS(lowerBound, upperBound, tnew, lnAcceptance,
                                    isAccepted, inode);

#ifdef TEST_NEW_DATA_STRUCTURE
        //test all
        this->testLocusEmbeddedGenealogy();
#endif

        //END of test

    } // end of loop

    dataLogLd = (dataLogLikelihood_ - dataLogLd);
    genLogLd = (genLogLikelihood_ - genLogLd);

#ifdef ENABLE_OMP_THREADS
#pragma omp atomic
#endif

#ifdef ENABLE_OMP_THREADS
#pragma omp atomic
#endif
    pState_->dataLogLikelihood += dataLogLd;
#ifdef ENABLE_OMP_THREADS
#pragma omp atomic
#endif
    pState_->logLikelihood += (genLogLd + dataLogLd) / pSetup_->numLoci;

    return (accepted);
}


//updates old data structures (DS) by calling to old inner functions
//as well as verifies consistency of values computed
int LocusEmbeddedGenealogy::updateGB_InternalNode_oldDS(double lowerBound,
                                                        double upperBound,
                                                        double tnew,
                                                        double lnAcceptance,
                                                        bool wasAccepted,
                                                        int inode) {
    int accepted = 0;

    int pop, son;
    double t, lnacceptance, lnLd;
    double genetree_lnLd_delta;
    int mig;
    double tb[2];

    t = getNodeAge(dataState.lociData[locusID_], inode);
    pop = nodePops[locusID_][inode];

    // set lower and upper bounds for new age
    tb[0] = dataSetup.popTree->pops[pop]->age;
    if (pop != dataSetup.popTree->rootPop)
        tb[1] = dataSetup.popTree->pops[pop]->father->age;
    else
        tb[1] = OLDAGE;

    mig = findFirstMig(locusID_, inode, -1);
    if (mig >= 0)
        tb[1] = min2(tb[1], genetree_migs[locusID_].mignodes[mig].age);
    else if (inode != getLocusRoot(dataState.lociData[locusID_]))
        tb[1] = min2(tb[1], getNodeAge(dataState.lociData[locusID_],
                                       getNodeFather(
                                               dataState.lociData[locusID_],
                                               inode)));

    //assert upper bound is as calculated by new method
    assert(fabs(tb[1] - upperBound) < EPSILON);

    for (int i = 0; i < 2; i++) {
        son = getNodeSon(dataState.lociData[locusID_], inode, i);
        mig = findLastMig(locusID_, son, -1);
        if (mig >= 0)
            tb[0] = max2(tb[0], genetree_migs[locusID_].mignodes[mig].age);
        else
            tb[0] = max2(tb[0],
                         getNodeAge(dataState.lociData[locusID_], son));
    }

    //assert lower bound is as calculated by new method
    assert(fabs(tb[0] - lowerBound) < EPSILON);

    //
    assert(abs(tnew - t) >= 1e-15);

    // update node's age, and compute delta log-likelihood
    adjustGenNodeAge(dataState.lociData[locusID_], inode, tnew);
    lnLd = -getLocusDataLikelihood(dataState.lociData[locusID_]);
    lnLd += computeLocusDataLikelihood(
            dataState.lociData[locusID_], 1);

    genetree_lnLd_delta = considerEventMove(locusID_, 0,
                                            nodeEvents[locusID_][inode],
                                            pop, t, pop, tnew);
    lnacceptance = genetree_lnLd_delta + lnLd;

    //assert acceptance bound is as calculated by new method
    assert(fabs(lnacceptance - lnAcceptance) < EPSILON);

    if (wasAccepted) {
        locus_data[locusID_].genLogLikelihood += genetree_lnLd_delta;
        acceptEventChainChanges(locusID_, 0);
        resetSaved(dataState.lociData[locusID_]);
    } else {
        // reject changes and revert to saved version
        rejectEventChainChanges(locusID_, 0);
        revertToSaved(dataState.lociData[locusID_]);//->this changes node again

    }

    return accepted;
}


/*
 * considerIntervalMove
 * computes the modifications required for changing a specific genetree
   by moving an interval to a new position (in populations ancestral or
   descendant to source population). Interval can be a coalescence or migration
   event. Computes change in interval chains (by adding
   new interval and not removing the original one), and changes in statistics
   required for fast computation of genetree likelihood.
   This procedure	updates all fields of genetree_stats_delta and returns the
   delta in log-likelihood of genetree due to this step.
*/
double
LocusEmbeddedGenealogy::considerIntervalMove(TreeNode *pNode, double newAge) {

    //compute delta log-likelihood
    double lnLd = -genealogy_.getLocusDataLikelihoodWrap();
    genealogy_.adjustGenNodeAgeWrap(pNode->getNodeId(), newAge);
    lnLd += genealogy_.computeLocusDataLikelihoodWrap(1);

    //get interval pointed by tree node
    PopInterval* pInterval = pNode->getInterval();

    //get interval type
    IntervalType  type = pInterval->getType();

    //create new interval
    PopInterval *pNewInterval = intervalsPro_.createInterval(pNode->getPop(),
                                                             newAge, type);

    int deltaNLin; //delta num lineages
    PopInterval *pBottomInterval, *pTopInterval; //bottom and top intervals

    //if new age is greater than original age
    if (newAge > pNode->getAge()) {
        deltaNLin = type == IntervalType::OUT_MIG ? -1 : 1;
        pBottomInterval = pInterval;
        pTopInterval = pNewInterval;

    } else {
        deltaNLin = type == IntervalType::OUT_MIG ? 1 : -1;
        pBottomInterval = pNewInterval;
        pTopInterval = pInterval;
    }

    //compute changes in coalescence and migration statistics
    intervalsPro_.computeStatsDelta(pBottomInterval, pTopInterval, deltaNLin);

    //compute delta log-likelihood
    double delta_lnLd = this->computeLogLikelihood(true);

    //set age of tree node to new age
    pNode->setAge(newAge);

    //set pointer of new interval to current node
    pNewInterval->setTreeNode(pNode);

    //set pointer of tree node to the new interval
    int intervalIndex = type == IntervalType::OUT_MIG ? 1 : 0;
    pNode->setInterval(pNewInterval, intervalIndex);

    //detach old interval from chain and return it to pool
    intervalsPro_.returnToPool(pInterval);

    //return log acceptance
    return lnLd + delta_lnLd;

}


/*
 * computeLogLikelihood
 * Computes log-likelihood of genetree using locus statistics.
 * If computeDelta flag is on, then compute delta log-likelihood
 * between proposal and original (proposal-original)
 * @return: computed value
*/
double LocusEmbeddedGenealogy::computeLogLikelihood(bool computeDelta) {

    if (computeDelta)
        return intervalsPro_.computeLogLikelihood(&intervalsOri_);
    else
        return intervalsPro_.computeLogLikelihood();
}

void LocusEmbeddedGenealogy::updateGenLogLikelihood() {
    genLogLikelihood_ = this->computeLogLikelihood();
}






