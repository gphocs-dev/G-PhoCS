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

class LocusEmbeddedGenealogy {

private:

    int locusID_;   //id of locus (genealogy id)

    LocusGenealogy     genealogy_; //object of genealogy
    LocusPopIntervals  intervals_; //object of intervals

    std::map<int,int> leafToPop_; //map between leaf to its pop
    std::map<int,std::vector<int>> popToLeaves_; //map between pop to its leaves

    int pop_queue_[2 * NSPECIES - 1]; // post-order queue of populations //todo

    DATA_SETUP*     pDataSetup_;       //pointer to DATA_SETUP struct
    PopulationTree* pPopTree_;         //pointer to PopulationTree struct
    DATA_STATE*     pDataState_;       //pointer to DATA_STATE struct
    GENETREE_MIGS*  pGenetreeMigs_;    //pointer to GENETREE_MIGS struct

public:

    //constructor
    LocusEmbeddedGenealogy(int locusID, int numIntervals,
                           DATA_SETUP *pDataSetup, DATA_STATE *pDataState,
                           GENETREE_MIGS *pGenetreeMigs);

    //construct both genealogy and intervals and connect between them
    int construct_genealogy_and_intervals();

    int computeGenetreeStats();

    double recalcStats(int pop);

    //get locus ID
    int getLocusID();

    //get locus data
    LocusData* getLocusData();

    //print popToLeaves
    void printPopToLeaves();

    //print leafToPop
    void printLeafToPop();

    //print all
    void printAll();

    //get all leaves of a population
    std::vector<int>& getPopLeaves(int pop);

    //get population of a leaf
    int getLeafPop(int leafId);

    //test genealogy new structure against the original
    void testLocusGenealogy();

    //test intervals new structure against the original
    void testPopIntervals();


};


#endif //G_PHOCS_LOCUSALLDATA_H
