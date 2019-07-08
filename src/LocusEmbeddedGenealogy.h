//
// Created by nomihadar on 3/11/19.
//

#ifndef G_PHOCS_LOCUSALLDATA_H
#define G_PHOCS_LOCUSALLDATA_H

#include "LocusGenealogy.h"
#include "LocusPopIntervals.h"
#include "MCMCcontrol.h"
#include "GPhoCS.h"
#include "MemoryMng.h"

#include <map>


/*=============================================================================
 *
 * LocusEmbeddedGenealogy class
 *
 * Class LocusEmbeddedGenealogy combines genealogy with its corresponding
 * intervals
 *
 * Contains:
 * 1. Object of locus genealogy.
 * 2. Object of locus pop intervals.
 * 3. Locus ID.
 * 4. Pointers to old structs.
 *===========================================================================*/

class LocusEmbeddedGenealogy {

private:

    LocusGenealogy     genealogy_; //object of genealogy
    LocusPopIntervals  intervals_; //object of intervals

    int locusID_;   //id of locus (genealogy id)

    DATA_SETUP*     pDataSetup_;       //pointer to DATA_SETUP struct
    DATA_STATE*     pDataState_;       //pointer to DATA_STATE struct
    GENETREE_MIGS*  pGenetreeMigs_;    //pointer to GENETREE_MIGS struct

public:

    //constructor
    LocusEmbeddedGenealogy(int locusID, int numIntervals,
                           DATA_SETUP* pDataSetup, DATA_STATE* pDataState,
                           GENETREE_MIGS* pGenetreeMigs);

    //copy-constructor
    LocusEmbeddedGenealogy(const LocusEmbeddedGenealogy& other);

    //copy without construction
    void copy(const LocusEmbeddedGenealogy& other);

public:

    // ********************* MAIN methods *********************

    //construct both genealogy and intervals and connect between them
    int constructEmbeddedGenealogy();

    //compute genealogy tree statistics
    int computeGenetreeStats();

    //recalculate statistics
    double recalcStats(int pop);

    int updateGB_InternalNode(double finetune);

    // ********************* GET methods *********************

    //get a reference to statistics
    const GenealogyStats& getStats() const;

    //get locus ID
    int getLocusID();

    //get locus data
    LocusData* getLocusData();

    // ********************* PRINT methods *********************

    //print embedded genealogy
    void printEmbeddedGenealogy();

    // ********************* TEST methods *********************

    //test if the new genealogy data structure is consistent with the original
    void testLocusGenealogy();

    //test if the new events data structure is consistent with the original
    void testPopIntervals();

    //test if statistics are consistent with the original
    void testGenealogyStats();

};


#endif //G_PHOCS_LOCUSALLDATA_H
