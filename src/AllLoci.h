
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
 * Vector of LocusEmbeddedGenealogy objects

 *===========================================================================*/
class AllLoci {

    std::vector<LocusEmbeddedGenealogy> lociVector_; //vector of loci

    GenealogyStats genealogiesStatsTotal_; //total statistics for all loci

public:

    //contructor
    AllLoci();

    LocusEmbeddedGenealogy& getLocus(int locusID);

    //test locus
    void testLocus(int locusID);

    //test all loci
    void testLoci();
};


#endif //G_PHOCS_LOCIEMBEDDED_H
