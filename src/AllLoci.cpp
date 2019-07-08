
#include "AllLoci.h"

#include "cassert"

/*
    AllLoci constructor
    Allocates LocusEmbeddedGenealogy objects,
    each with a unique locus ID.
*/
AllLoci::AllLoci() : lociVector_(),
                     statsTotal_(dataSetup.popTree->numPops,
                                 dataSetup.popTree->numMigBands) {

    lociVector_.reserve(dataSetup.numLoci);

    //construct N locus embedded genealogy
    for (int locusID = 0; locusID < dataSetup.numLoci; ++locusID) {
        lociVector_.emplace_back(locusID, MAX_EVENTS, &dataSetup, &dataState,
                                 genetree_migs);
    }

    //initialize mig band times
    initializeMigBandTimes(dataSetup.popTree);

}


/*
    calcLociTotalStats
    Calculates statistics of all loci
*/
void AllLoci::calcLociTotalStats() {

   statsTotal_.resetStatsTotal();

    //for each locus
    for (auto& locus : lociVector_) {

        //get statistics
        const GenealogyStats & stats = locus.getStats();

        //for each coal statistics
        for (int pop = 0; pop < dataSetup.popTree->numPops; pop++) {
            statsTotal_.coals[pop].num += stats.coals[pop].num;
            statsTotal_.coals[pop].stats += stats.coals[pop].stats;
        }

        //for each mig statistics
        for (int id = 0; id < dataSetup.popTree->numMigBands; id++) {
            statsTotal_.migs[id].num += stats.migs[id].num;
            statsTotal_.migs[id].stats += stats.migs[id].stats;
        }
    }
}


/*
    testLociTotalStats
    Test if total statistics are consistent with old version
*/
void AllLoci::testLociTotalStats() {

    //verify coal statistics are equal
    for (int pop = 0; pop < dataSetup.popTree->numPops; pop++) {
        assert(statsTotal_.coals[pop].num == genetree_stats_total.num_coals[pop]);
        assert(fabs(statsTotal_.coals[pop].stats - genetree_stats_total.coal_stats[pop]) < EPSILON);
    }

    //verify migs statistics are equal
    for (int id = 0; id < dataSetup.popTree->numMigBands; id++) {
        assert(statsTotal_.migs[id].num == genetree_stats_total.num_migs[id]);
        assert(fabs(statsTotal_.migs[id].stats - genetree_stats_total.mig_stats[id]) < EPSILON);
    }
}


/*
    testLoci
    Test if all loci data is consistent with old version
*/
void AllLoci::testLoci() {

    //for each locus
    for (auto& locus : lociVector_) {

        //construct mig bands times
        constructMigBandsTimes(dataSetup.popTree);

        //construct genealogy and intervals
        locus.constructEmbeddedGenealogy();

        //compute genealogy statistics
        locus.computeGenetreeStats();

        //copy locus
        LocusEmbeddedGenealogy copy(locus);

        //test locus genealogy
        copy.testLocusGenealogy();

        //test locus intervals
        copy.testPopIntervals();

        //test locus statistics
        copy.testGenealogyStats();


        //locus.printEmbeddedGenealogy();
        //copy.printEmbeddedGenealogy();
        //printGenealogyAndExit(locus.getLocusID(), -1);
        //break;

    }

    //calculate total statistics
    this->calcLociTotalStats();

    //test total statistics
    this->testLociTotalStats();
}


/*
    getLocus
    Get locus by id
    @param: locus id
    @return: reference to locus embedded genealogy
*/
LocusEmbeddedGenealogy& AllLoci::getLocus(int locusID) {
    return lociVector_[locusID];
}


