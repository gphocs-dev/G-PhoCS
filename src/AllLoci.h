
#ifndef G_PHOCS_LOCIEMBEDDED_H
#define G_PHOCS_LOCIEMBEDDED_H

#include "LocusEmbeddedGenealogy.h"
#include "GenealogyStats.h"

#include <vector>


/*=============================================================================
 *
 * AllLoci class
 *
 * Class AllLoci contains all loci embedded genealogies
 *
 * Contains:
 * 1. Vector of LocusEmbeddedGenealogy objects.
 * 2. Object of total statistics.
 *===========================================================================*/
class AllLoci {

    std::vector<LocusEmbeddedGenealogy> lociVector_; //vector of loci

    GenealogyStats statsTotal_; //total statistics for all loci

public:

    //constructor
    AllLoci();

    LocusEmbeddedGenealogy& getLocus(int locusID);

    void calcTotalStats();

    void testGenealogyStats();

    //test all loci
    void testLoci();


};


#endif //G_PHOCS_LOCIEMBEDDED_H
