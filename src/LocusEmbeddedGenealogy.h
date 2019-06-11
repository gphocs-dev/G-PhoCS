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
 * 4. Map between leaf to its pop.
 * 5. Map between pop to its leaves.
 * 6. Pointers to old structs
 *===========================================================================*/

class LocusEmbeddedGenealogy {

private:

    LocusGenealogy     genealogy_; //object of genealogy
    LocusPopIntervals  intervals_; //object of intervals

    int locusID_;   //id of locus (genealogy id)

    std::map<int,int> leafToPop_; //map between leaf to its pop
    std::map<int,std::vector<int>> popToLeaves_; //map between pop to its leaves

    DATA_SETUP*     pDataSetup_;       //pointer to DATA_SETUP struct
    PopulationTree* pPopTree_;         //pointer to PopulationTree struct
    DATA_STATE*     pDataState_;       //pointer to DATA_STATE struct
    GENETREE_MIGS*  pGenetreeMigs_;    //pointer to GENETREE_MIGS struct

public:

    //constructor
    LocusEmbeddedGenealogy(int locusID, int numIntervals,
                           DATA_SETUP* pDataSetup, DATA_STATE* pDataState,
                           GENETREE_MIGS* pGenetreeMigs);

    //construct both genealogy and intervals and connect between them
    int construct_genealogy_and_intervals();

    //get locus ID
    int getLocusID();

    //get locus data
    LocusData* getLocusData();

    //get all leaves of a population
    std::vector<int>& getPopLeaves(int pop);

    //get population of a leaf
    int getLeafPop(int leafId);

    //get a reference to statistics
    const GenealogyStats& getStats() const;

    //compute genealogy tree statistics
    int computeGenetreeStats();

    //recalculate statistics
    double recalcStats(int pop);

    // **************** print methods ****************

    //print popToLeaves
    void printPopToLeaves();

    //print leafToPop
    void printLeafToPop();

    //print all data
    void printEmbeddedGenealogy();

    // **************** test methods ****************

    //test if the new genealogy data structure is consistent with the original
    void testLocusGenealogy();

    //test if the new events data structure is consistent with the original
    void testPopIntervals();

    //test if statistics are consistent with the original
    void testGenealogyStats();


};


#endif //G_PHOCS_LOCUSALLDATA_H
