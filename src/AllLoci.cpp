
#include "AllLoci.h"

#include "cassert"

/*
    AllLoci constructor
    Allocates LocusEmbeddedGenealogy objects,
    each with a unique locus ID.
*/
AllLoci::AllLoci()
        : lociVector_(), statsTotal_(dataSetup.popTree->numPops,
                                     dataSetup.popTree->numMigBands) {

    lociVector_.reserve(dataSetup.numLoci);

    for (int locusID = 0; locusID < dataSetup.numLoci; ++locusID) {
        lociVector_.emplace_back(locusID, MAX_EVENTS, &dataSetup, &dataState, //todo: max events?
                                 genetree_migs);
    }

    //initialize mig band times
    initializeMigBandTimes(dataSetup.popTree);

}


void AllLoci::calcTotalStats() {

   statsTotal_.resetStatsTotal();

    //for each locus
    for (auto& locus : lociVector_) {

        //get statistics
        const GenealogyStats & stats = locus.getStats();

        //for each coal statistics
        for (int pop = 0; pop < dataSetup.popTree->numPops; pop++) {
            statsTotal_.coal[pop].num += stats.coal[pop].num;
            statsTotal_.coal[pop].stats += stats.coal[pop].stats;
        }

        //for each mig statistics
        for (int id = 0; id < dataSetup.popTree->numMigBands; id++) {
            statsTotal_.migs[id].num += stats.migs[id].num;
            statsTotal_.migs[id].stats += stats.migs[id].stats;
        }
    }
}


void AllLoci::testGenealogyStats() {

    //verify coal statistics are equal
    for (int pop = 0; pop < dataSetup.popTree->numPops; pop++) {
        assert(statsTotal_.coal[pop].num == genetree_stats_total.num_coals[pop]);
        assert(fabs(statsTotal_.coal[pop].stats - genetree_stats_total.coal_stats[pop]) < EPSILON);
    }

    //verify migs statistics are equal
    for (int id = 0; id < dataSetup.popTree->numMigBands; id++) {
        assert(statsTotal_.migs[id].num == genetree_stats_total.num_migs[id]);
        assert(fabs(statsTotal_.migs[id].stats - genetree_stats_total.mig_stats[id]) < EPSILON);
    }
}


void AllLoci::testLoci() {

    //for each locus
    for (auto& locus : lociVector_) {

        //construct mig bands times
        constructMigBandsTimes(dataSetup.popTree);

        //construct genealogy and intervals
        locus.construct_genealogy_and_intervals();

        //compute genealogy statistics
        locus.computeGenetreeStats();

        //test locus genealogy
        locus.testLocusGenealogy();

        //test locus intervals
        locus.testPopIntervals();

        //test locus statistics
        locus.testGenealogyStats();

        //locus.printEmbeddedGenealogy();
        //printGenealogyAndExit(locus.getLocusID(), -1);
        //break;

    }

    //calculate total statistics
    this->calcTotalStats();

    //test total statistics
    this->testGenealogyStats();
}


LocusEmbeddedGenealogy& AllLoci::getLocus(int locusID) {
    return lociVector_[locusID];
}


