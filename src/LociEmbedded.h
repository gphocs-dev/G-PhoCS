//
// Created by nomihadar on 4/1/19.
//

#ifndef G_PHOCS_LOCIEMBEDDED_H
#define G_PHOCS_LOCIEMBEDDED_H

#include "LocusEmbeddedGenealogy.h"
#include <vector>


/*=============================================================================
 *
 * LociEmbedded class
 *
 * Class LociEmbedded contains all loci embedded genealogy
 *
 * Contains:
 * Vector of LocusEmbeddedGenealogy objects

 *===========================================================================*/
class LociEmbedded {

    std::vector<LocusEmbeddedGenealogy> lociVector_; //vector of loci

public:

    //contructor
    LociEmbedded();

    LocusEmbeddedGenealogy& getLocus(int locusID);

    //test locus
    void testLocus(int locusID);

    //test all loci
    void testLoci();
};


#endif //G_PHOCS_LOCIEMBEDDED_H
