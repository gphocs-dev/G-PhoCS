
#include "LocusEmbeddedGenealogy.h"
#include "DbgErrMsgIntervals.h"
#include "PopulationTree.h"

#include <iostream>
#include <iomanip>
#include <cmath>

/*
	Constructor
*/
LocusEmbeddedGenealogy::LocusEmbeddedGenealogy(
        int locusID,
        int numIntervals,
        DATA_SETUP *pDataSetup,
        DATA_STATE *pDataState,
        GENETREE_MIGS *pGenetreeMigs)

        : genealogy_(pDataSetup->numSamples), //construct genealogy
          intervals_(locusID, numIntervals,
                     pDataSetup->popTree),  //construct intervals

          locusID_(locusID),
          pDataSetup_(pDataSetup),
          pPopTree_(pDataSetup->popTree),
          pDataState_(pDataState),
          pGenetreeMigs_(pGenetreeMigs) {

    //create map between leaf to its pop
    // and map between pop to its leaves
    for (int node = 0; node < pDataSetup_->numSamples; node++) {
        int pop = nodePops[locusID_][node];
        leafToPop_[node] = pop;
        popToLeaves_[pop].emplace_back(node); //todo reserve?
    }

    //fill queue with pops, sorted by post order
    populationPostOrder(pPopTree_->rootPop, pop_queue_);

}


/*
   @return: a pointer to locus data of current locus
*/
LocusData *LocusEmbeddedGenealogy::getLocusData() {
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

    //reset genealogy
    genealogy_.reset();

    //construct genealogy branches (edges between tree nodes)
    genealogy_.constructBranches(this->getLocusData());

    //reset intervals
    intervals_.reset();

    //link intervals to each other
    intervals_.linkIntervals();

    //add start and end intervals
    intervals_.createStartEndIntervals();

    //create samples start intervals (for ancient samples)
    for (int pop = 0; pop < pPopTree_->numCurPops; pop++) {

        //create interval
        double age = pPopTree_->pops[pop]->sampleAge;
        PopInterval *pInterval = intervals_.createInterval(pop, age,
                                                           IntervalType::SAMPLES_START);
        if (!pInterval) {
            INTERVALS_FATAL_0024
        }
    }

    //create coalescent intervals
    //and link intervals to genealogy and vice versa
    int nSamples = pDataSetup_->numSamples;
    for (int node = 0; node < 2 * nSamples - 1; node++) {

        //get population and age of node
        int pop = nodePops[locusID_][node];
        double age = getNodeAge(getLocusData(), node);

        // if node is a leaf - link it to its sampleStart interval
        if (genealogy_.isLeaf(node)) {

            //get samples start interval of pop
            PopInterval *pInterval = intervals_.getSamplesStart(pop);

            //get leaf node by current node id
            LeafNode *pNode = genealogy_.getLeafNode(node);

            //leaf node points to samplesStart interval (but samplesStart
            // interval points to null since there are several leaves)
            pNode->setSamplesInterval(pInterval);


        } else { //if node is not a leaf create a coal interval and link to node

            //create a coalescent interval
            PopInterval *pInterval =
                    intervals_.createInterval(pop, age, IntervalType::COAL);

            if (!pInterval) {
                INTERVALS_FATAL_0025
            }

            //get coal node by current node id
            CoalNode *pNode = genealogy_.getCoalNode(node);

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
        TreeNode *pTreeNode = genealogy_.getTreeNodeByID(node);

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
            PopInterval *pMigIn =
                    intervals_.createInterval(target_pop, age,
                                              IntervalType::IN_MIG);
            if (!pMigIn) {
                INTERVALS_FATAL_0022
            }

            //create an outgoing migration interval
            PopInterval *pMigOut =
                    intervals_.createInterval(source_pop, age,
                                              IntervalType::OUT_MIG);
            if (!pMigOut) {
                INTERVALS_FATAL_0023
            }

            //add a migration node to genealogy
            MigNode *pMigNode = genealogy_.addMigNode(pTreeNode);

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


//todo: rewrite recalcStats
/* recalcStats
   Re-calculates stats for given population in given gen.
   Writes down stats in genetree_stats_check and then compares to prior stats to return the log-likelihood.
   The log-likelihood computation ignores changes in number of migs/coals!!
*/

double LocusEmbeddedGenealogy::recalcStats(int pop) {

    int id, mig_band, event;
    int live_mig_bands[MAX_MIG_BANDS];
    int num_live_mig_bands = 0;
    double t, heredity_factor = 1;
    double delta_lnLd = 0.0;

    double age;

    locus_data[locusID_].genetree_stats_check.coal_stats[pop] = 0.0;
    locus_data[locusID_].genetree_stats_check.num_coals[pop] = 0;

    //get first interval
    PopInterval *pInterval = intervals_.getFirstInterval(pop);

    //get num lineages of first interval
    int n = pInterval->getNumLineages();

    // follow intervals chain and set number of lineages per interval according
    // to previous interval also update statistics
    while (true) {

        pInterval->setNumLineages(n);

        //id = event_chains[locusID_].events[event].getId();
        MigrationBand* migBand = &dataSetup.popTree->migBands[id];

        t = pInterval->getElapsedTime();

        locus_data[locusID_].genetree_stats_check.coal_stats[pop] +=
                n * (n - 1) * t;

        age = pInterval->getAge();//todo: check
        TimeBand *timeBand = getLiveMigBands(dataSetup.popTree, pop, age);
        MigrationBand *mig;
        //for (auto pMigBand) {
        //
        //}

        for (mig_band = 0; mig_band < num_live_mig_bands; mig_band++) {
            locus_data[locusID_].genetree_stats_check.mig_stats[live_mig_bands[mig_band]] +=
                    n * t;
        }

        switch (pInterval->getType()) {
            case (IntervalType::SAMPLES_START):
                n += dataSetup.numSamplesPerPop[pop];
                break;

            case (IntervalType::COAL):
                locus_data[locusID_].genetree_stats_check.num_coals[pop]++;
                n--;
                break;

            case (IntervalType::IN_MIG):
                // figure out migration band and update its statistics
                //mig_band = genetree_migs[locusID_].mignodes[id].migration_band;
                //locus_data[locusID_].genetree_stats_check.num_migs[mig_band]++;
                n--;
                break;

            case (IntervalType::OUT_MIG):
                n++;
                break;

                /*case (MIG_BAND_START):
                    live_mig_bands[num_live_mig_bands++] = id;
                    // initialize statistics for this new migration band
                    locus_data[locusID_].genetree_stats_check.num_migs[id] = 0;
                    locus_data[locusID_].genetree_stats_check.mig_stats[id] = 0.0;
                    break;*/

                /*case (MIG_BAND_END):
                    // compare and copy stats for mig band
                    delta_lnLd -= (locus_data[locusID_].genetree_stats_check.mig_stats[id] - genetree_stats[locusID_].mig_stats[id])
                    * dataSetup.popTree->migBands[id].migRate;

                    #ifdef ENABLE_OMP_THREADS
                    #pragma omp atomic
                    #endif
                    genetree_stats_total.mig_stats[id] +=
                            locus_data[locusID_].genetree_stats_check.mig_stats[id] -
                            genetree_stats[locusID_].mig_stats[id];
    #ifdef ENABLE_OMP_THREADS
    #pragma omp atomic
    #endif
                    genetree_stats_total.num_migs[id] +=
                            locus_data[locusID_].genetree_stats_check.num_migs[id] -
                            genetree_stats[locusID_].num_migs[id];
                    genetree_stats[locusID_].mig_stats[id] = locus_data[locusID_].genetree_stats_check.mig_stats[id];
                    genetree_stats[locusID_].num_migs[id] = locus_data[locusID_].genetree_stats_check.num_migs[id];
                    // remove mig band from living list
                    for (mig_band = 0; mig_band < num_live_mig_bands; mig_band++) {
                        if (live_mig_bands[mig_band] == id)
                            break;
                    }
                    if (mig_band == num_live_mig_bands) {
                        if (debug) {
                            fprintf(stderr,
                                    "\nError: recalcStats: migration band %d not alive in population %d, gen %d.\n",
                                    event_chains[locusID_].events[event].getId(), pop, locusID_);
                            fprintf(stderr, "Live migration bands are:");
                            for (mig_band = 0; mig_band < num_live_mig_bands; mig_band++) {
                                fprintf(stderr, " %d", live_mig_bands[mig_band]);
                            }
                            fprintf(stderr, "\n");
                        } else {
                            fprintf(stderr, "Fatal Error 0025.\n");
                        }
                        printGenealogyAndExit(locusID_, -1);
                    }
                    live_mig_bands[mig_band] = live_mig_bands[--num_live_mig_bands];
                    break;
                    */
            case (IntervalType::DUMMY):
            case (IntervalType::POP_END):
                break;

            default:
                if (debug) {
                    std::cout
                            << endl
                            << "Error: recalcStats: event of unknown type "
                            << pInterval->typeToStr() << " in population "
                            << pop
                            << " gen " << locusID_ << endl;

                } else {
                    std::cout << "Fatal Error 0026." << endl;
                }
                printGenealogyAndExit(locusID_, -1);
                break;
        }// end of switch

        if (pInterval->isType(IntervalType::POP_END))
            break;


        pInterval = pInterval->getNext();

    }// end of while

    /*
    if (num_live_mig_bands != 0) {
        if (debug) {
            fprintf(stderr,
                    "\nError: recalcStats: number of live mig bands %d at end of population %d in gen %d.\n",
                    num_live_mig_bands, pop, locusID_);
        } else {
            fprintf(stderr, "Fatal Error 0027.\n");
        }
        printGenealogyAndExit(locusID_, -1);
    }

    delta_lnLd -= (locus_data[locusID_].genetree_stats_check.coal_stats[pop] -
                   genetree_stats[locusID_].coal_stats[pop]) /
                  (dataSetup.popTree->pops[pop]->theta * heredity_factor);
#ifdef ENABLE_OMP_THREADS
#pragma omp atomic
#endif
    genetree_stats_total.coal_stats[pop] +=
            (locus_data[locusID_].genetree_stats_check.coal_stats[pop] -
             genetree_stats[locusID_].coal_stats[pop]) / heredity_factor;
#ifdef ENABLE_OMP_THREADS
#pragma omp atomic
#endif
    genetree_stats_total.num_coals[pop] +=
            locus_data[locusID_].genetree_stats_check.num_coals[pop] -
            genetree_stats[locusID_].num_coals[pop];
    genetree_stats[locusID_].coal_stats[pop] = locus_data[locusID_].genetree_stats_check.coal_stats[pop];
    genetree_stats[locusID_].num_coals[pop] = locus_data[locusID_].genetree_stats_check.num_coals[pop];
    */

    return delta_lnLd;
}
/*** end of recalcStats ***/



/*	computeGenetreeStats
	Computes the statistics of a given locus.
    Assumes intervals chains are built, but number of lineages is ONLY
    set for first intervals in the leaf populations.
	Sets number of lineages for each non-leaf event by traversing the
	population tree post-order. In parallel, also records the statistics.
*/

int LocusEmbeddedGenealogy::computeGenetreeStats() {

    // go over all intervals and compute num_lineages per each interval
    // also update genetree statistics
    for (int i = 0; i < dataSetup.popTree->numPops; i++) {

        //get current pop
        int pop = pop_queue_[i];

        // if not leaf population get number of in-lineages from end-intervals
        // of son populations
        if (pop >= dataSetup.popTree->numCurPops) {

            //get num lineages of first interval of the left pop son
            int lSon = dataSetup.popTree->pops[pop]->sons[0]->id;
            int n1 = intervals_.getPopEnd(lSon)->getNumLineages();

            //get num lineages of first interval of the right pop son
            int rSon = dataSetup.popTree->pops[pop]->sons[1]->id;
            int n2 = intervals_.getPopEnd(rSon)->getNumLineages();

            //set num lineages of first interval of current pop to the sum
            //of sons' num lineages
            intervals_.getFirstInterval(pop)->setNumLineages(n1 + n2);

        } else {
            //set num lineages of first interval of current pop to 0
            intervals_.getFirstInterval(pop)->setNumLineages(0);
        }

        recalcStats(pop);

    }

    return 0;
}


/*
	print pop to leaves map
*/
void LocusEmbeddedGenealogy::printPopToLeaves() {

    std::cout << "Pop to leaves: " << std::endl;
    for (auto &x : popToLeaves_) {
        std::cout << "pop " << x.first << ", leaves: ";
        for (int leaf : x.second) {
            cout << leaf << ", ";
        }
        std::cout << std::endl;
    }
}

/*
	print leaf to map
*/
void LocusEmbeddedGenealogy::printLeafToPop() {
    std::cout << "Leaf to pop: " << std::endl;
    for (auto &x : leafToPop_) {
        std::cout << "leaf " << x.first << ", "
                  << "pop: " << x.second << std::endl;
    }
}

/*
	print population tree, genealogy and intervals
*/
void LocusEmbeddedGenealogy::printEmbeddedGenealogy() {

    //print population tree
    printPopulationTree(this->pDataSetup_->popTree, stderr, 1);

    //print pop to leaves
    std::cout << "------------------------------------------------------"
              << std::endl;
    this->printPopToLeaves();

    //print leaf to pop
    std::cout << "------------------------------------------------------"
              << std::endl;
    this->printLeafToPop();

    //print genealogy
    std::cout << "------------------------------------------------------"
              << std::endl;
    genealogy_.printGenealogy();

    //print intervals
    std::cout << "------------------------------------------------------"
              << std::endl;
    intervals_.printIntervals();

}

int LocusEmbeddedGenealogy::getLocusID() {
    return locusID_;
}


std::vector<int> &LocusEmbeddedGenealogy::getPopLeaves(int pop) {
    return popToLeaves_[pop];
}

int LocusEmbeddedGenealogy::getLeafPop(int leafId) {
    return leafToPop_[leafId];
}

void LocusEmbeddedGenealogy::testLocusGenealogy() {

    //define precision for double comparision
    double PRECISION = 0.0000000001;

    //get locus id and locus data
    LocusData* pLocusData = this->getLocusData();
    int numSamples = pDataSetup_->numSamples;

    //get all migs into a map: <node,[migsAges]>
    std::map<int, std::vector<double>> migsMap;
    for (int node = 0; node < 2 * numSamples - 1; node++) {
        int mig = findFirstMig(locusID_, node, getNodeAge(pLocusData, node));
        while (mig != -1) {
            double age = pGenetreeMigs_[locusID_].mignodes[mig].age;
            migsMap[node].emplace_back(age);
            mig = findFirstMig(locusID_, node, age);
        }
    }

    //iterate by node id
    for (int node = 0; node < 2 * numSamples - 1; node++) {

        //get coal or leaf node
        TreeNode* pNode = genealogy_.getTreeNodeByID(node);

        //get ages
        double age = getNodeAge(pLocusData, node);
        double ageNew = pNode->getAge();

        //compare ages
        assert(fabs(age - ageNew) < PRECISION);

        TreeNode* parentNew = pNode;

        //if there are migrations above
        for (double migAge : migsMap[node]) {

            parentNew = parentNew->getParent();
            double migAgeNew = parentNew->getAge();

            //compare ages
            assert(fabs(migAge - migAgeNew) < PRECISION);
        }

        //compare parents ages
        int parent = getNodeFather(pLocusData, node);
        parentNew = parentNew->getParent();
        if (parent != -1) {
            double age = getNodeAge(pLocusData, parent);
            double ageNew = parentNew->getAge();
            //compare ages
            assert(fabs(age - ageNew) < PRECISION);
        }

        //compare sons ages
        for (int son = 0; son < 2; son++) {
            int lSon = getNodeSon(pLocusData, node, son);
            TreeNode *pSonNew = pNode;

            //reverse vector
            vector<double> rev(migsMap[lSon].rbegin(), migsMap[lSon].rend());
            for (double migAge : rev) {

                pSonNew = son ? pSonNew->getRightSon() : pSonNew->getLeftSon();
                double migAgeNew = pSonNew->getAge();

                //compare ages
                assert(fabs(migAge - migAgeNew) < PRECISION);
            }

            pSonNew = son ? pSonNew->getRightSon() : pSonNew->getLeftSon();
            if (lSon != -1) {
                double age = getNodeAge(pLocusData, lSon);
                double ageNew = pSonNew->getAge();
                //compare ages
                assert(fabs(age - ageNew) < PRECISION);
            }
        }

    }

}


void LocusEmbeddedGenealogy::testPopIntervals() {


    //define precision for double comparision
    double PRECISION = 0.0000000001;

    //follow live mig bands
    std::vector<int> liveMigBand;

    //for each pop
    for (int pop = 0; pop < dataSetup.popTree->numPops; pop++) {

        //get first event in old structure
        int event = event_chains[locusID_].first_event[pop];

        //get first (pop start) interval in the new structure
        PopInterval* pInterval = intervals_.getFirstInterval(pop);

        //elapsed time of event
        double eventTime = 0;

        //iterate both old and new structures (events VS intervals)
        for (event; event >=
                    0; event = event_chains[locusID_].events[event].getNextIdx()) {

            //get event type and id (id for testing mig bands)
            EventType eventType = event_chains[locusID_].events[event].getType();

            //get elapsed time
            eventTime += event_chains[locusID_].events[event].getElapsedTime();


            switch (eventType) {
                case SAMPLES_START: {
                    //compare types
                    assert(pInterval->isType(IntervalType::SAMPLES_START));
                    break;
                }
                case COAL: {
                    //compare types
                    assert(pInterval->isType(IntervalType::COAL));
                    break;
                }
                case IN_MIG: {
                    //compare types
                    assert(pInterval->isType(IntervalType::IN_MIG));
                    break;
                }
                case OUT_MIG: {
                    //compare types
                    assert(pInterval->isType(IntervalType::OUT_MIG));
                    break;
                }
                case END_CHAIN: {
                    //compare types
                    assert(pInterval->isType(IntervalType::POP_END));
                    break;
                }
                case MIG_BAND_START: {
                    //add id of mig band to live mig bands
                    int event_id = event_chains[locusID_].events[event].getId();
                    liveMigBand.push_back(event_id);
                    break;
                }
                case MIG_BAND_END: {

                    int event_id = event_chains[locusID_].events[event].getId();

                    //remove id of mig band from live mig bands
                    liveMigBand.erase(
                            std::remove(liveMigBand.begin(), liveMigBand.end(),
                                        event_id), liveMigBand.end());
                    break;
                }

            }

            //compare live mig bands
            TimeBand* timeBand = getLiveMigBands(dataSetup.popTree, pop, pInterval->getAge());

            if (timeBand) {
                //compare num of live mig bands
                assert(liveMigBand.size() == timeBand->migBands.size());

                //for each live mig band
                for (int id : liveMigBand) {
                    //get pointer to mig band by its id
                    MigrationBand* migBand = &dataSetup.popTree->migBands[id];
                    //verify that current mig is found in the new data structure
                    assert(std::find(timeBand->migBands.begin(),
                                     timeBand->migBands.end(), migBand) != timeBand->migBands.end());
                }//end of comparing live mig band
            }


            //new structure doesn't have mig-band intervals therefor
            // if current event is a mig-band continue to next event
            // and don't compare elapsed time (n)
            if (eventType == MIG_BAND_START || eventType == MIG_BAND_END)
                continue;

            //compare num lineages
            int eventLin = event_chains[locusID_].events[event].getNumLineages();
            int intervalLin = pInterval->getNumLineages();
            assert(eventLin == intervalLin);

            //compare elapsed time
            double intervalTime = pInterval->getElapsedTime();
            assert(fabs(eventTime - intervalTime) < PRECISION);


            //get next interval
            pInterval = pInterval->getNext();

            //reset event time
            eventTime = 0;

            //reset live mig bands
            liveMigBand.clear();
        }
    }

}
