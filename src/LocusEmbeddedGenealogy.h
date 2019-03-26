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

    LocusGenealogy     genealogy_; //class of genealogy
    LocusPopIntervals  intervals_; //class of intervals

    std::map<int,int> leafToPop_; //map between leaf ID to its pop
    std::map<int,std::vector<int>> popToLeaves_; //map between pop to its leaves

    DATA_SETUP*     pDataSetup_;       //pointer to DATA_SETUP struct
    PopulationTree* pPopTree_;         //pointer to PopulationTree struct
    DATA_STATE*     pDataState_;       //pointer to DATA_STATE struct
    GENETREE_MIGS*  pGenetreeMigs_;    //pointer to GENETREE_MIGS struct

public:

    //constructor
    LocusEmbeddedGenealogy(int locusID, int numIntervals, DATA_SETUP* pDataSetup,
                 PopulationTree* pPopTree, DATA_STATE* pDataState,
                 GENETREE_MIGS* pGenetreeMigs);

    //construct
    int construct_genealogy_and_intervals();

    //get locus data
    LocusData* getLocusData();

    //print
    void print();

};


#endif //G_PHOCS_LOCUSALLDATA_H
