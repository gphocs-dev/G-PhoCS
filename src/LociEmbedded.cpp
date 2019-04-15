//
// Created by nomihadar on 4/1/19.
//



#include "LociEmbedded.h"
#include "assert.h"

/*
    LociEmbedded constructor
    Allocates LocusEmbeddedGenealogy objects,
    each with a unique locus ID.
*/
LociEmbedded::LociEmbedded() : lociVector_() {

    lociVector_.reserve(dataSetup.numLoci);

    for (int locusID = 0; locusID < dataSetup.numLoci; ++locusID) {
        lociVector_.emplace_back(locusID, MAX_EVENTS, &dataSetup, &dataState,
                                 genetree_migs);
    }
}


void LociEmbedded::testLocus(int locusID) {

    LocusEmbeddedGenealogy locus = lociVector_[locusID];

}

void LociEmbedded::testLoci() {

    for (auto& locus : lociVector_) {

        //construct genealogy and intervals
        locus.construct_genealogy_and_intervals();
        locus.computeGenetreeStats();
        locus.testLocusGenealogy();
        locus.testPopIntervals();
        //locus.printAll();
        //printGenealogyAndExit(locus.getLocusID(), -1);
        //break;

    }

}

LocusEmbeddedGenealogy& LociEmbedded::getLocus(int locusID) {
    return lociVector_[locusID];
}
