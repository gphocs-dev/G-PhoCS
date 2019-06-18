
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

    //get locus by id
    LocusEmbeddedGenealogy& getLocus(int locusID);

    //calculate statistics of all loci
    void calcLociTotalStats();

    //test if total statistics of all loci is consistent with old version
    void testLociTotalStats();

    //test if all loci data are consistent with old version
    void testLoci();

    //


};


#endif //G_PHOCS_LOCIEMBEDDED_H
