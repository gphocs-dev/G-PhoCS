
#include "AllLoci.h"
#include "cassert"

/*
    AllLoci constructor
    Allocates LocusEmbeddedGenealogy objects,
    each with a unique locus ID.
*/
AllLoci::AllLoci() : lociVector_() {

    lociVector_.reserve(dataSetup.numLoci);

    for (int locusID = 0; locusID < dataSetup.numLoci; ++locusID) {
        lociVector_.emplace_back(locusID, MAX_EVENTS, &dataSetup, &dataState,
                                 genetree_migs);
    }

    //initialize mig band times
    initializeMigBandTimes(dataSetup.popTree);
}


void AllLoci::testLocus(int locusID) {

    LocusEmbeddedGenealogy locus = lociVector_[locusID];

}

void AllLoci::testLoci() {

    for (auto& locus : lociVector_) {

        //construct genealogy and intervals
        locus.construct_genealogy_and_intervals();

        //compute genealogy statistics
        locus.computeGenetreeStats();

        //construct mig bands times
        constructMigBandsTimes(dataSetup.popTree);

        locus.testLocusGenealogy();
        locus.testPopIntervals();

        //locus.print();
        //printGenealogyAndExit(locus.getLocusID(), -1);
        //break;

    }

}

LocusEmbeddedGenealogy& AllLoci::getLocus(int locusID) {
    return lociVector_[locusID];
}
