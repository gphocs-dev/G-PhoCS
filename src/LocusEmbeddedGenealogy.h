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
 * 1. Locus ID.
 * 2. Object of locus genealogy.
 * 3. Two Objects of locus pop intervals - one for proposals and one is original.
 * 4. Log-likelihood of locus genealogy - P(gen|Model).
 * 5. Log-likelihood of locus data - P(data|gen).
 * 6. Pointers to old structs.
 *===========================================================================*/

class LocusEmbeddedGenealogy {

private:

    int locusID_;   //id of locus

    LocusGenealogy     genealogy_; //object of genealogy

    LocusPopIntervals  intervalsPro_; //object of intervals - proposal
    LocusPopIntervals  intervalsOri_; //object of intervals - original

    double genLogLikelihood_; //genealogy log-likelihood - P(gen|Model)
    double dataLogLikelihood_; //data log-likelihood - P(data|gen)

    DATA_SETUP*     pSetup_;       //pointer to DATA_SETUP struct
    DATA_STATE*     pState_;       //pointer to DATA_STATE struct
    GENETREE_MIGS*  pGenetreeMigs_;//pointer to GENETREE_MIGS struct

public:

    //constructor
    LocusEmbeddedGenealogy(int locusID, int numIntervals,
                           DATA_SETUP* pSetup, DATA_STATE* pState,
                           GENETREE_MIGS* pGenetreeMigs);

public:

    // ********************* MAIN methods *********************

    //construct both genealogy and intervals and connect between them
    int constructEmbeddedGenealogy();

    //compute genealogy tree statistics
    int computeGenetreeStats();

    //recalculate statistics
    void recalcStats(int pop);

    //
    int updateGB_InternalNode(double finetune);

    //
    double considerIntervalMove(TreeNode *pNode, double newAge);

    //compute delta log likelihood
    double computeLogLikelihood(bool computeDelta=false);


    // ********************* Copy methods *********************

    //copy without construction
    void copy(const LocusEmbeddedGenealogy& other);

    //copy intervals from proposal to original and vise versa
    void copyIntervals(bool accepted);

    // ********************* GET methods *********************

    //get a reference to statistics
    const GenealogyStats& getStats() const;

    //get log likelihood
    double getGenLogLikelihood() const;

    //get log likelihood
    double getDataLogLikelihood() const;

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

    //test if likelihood is consistent with the original
    void testLogLikelihood();

    //
    void testLocusEmbeddedGenealogy();

    //test if update internal node is consistent
    int test_updateGB_InternalNode(double lowerBound,
                                   double upperbound,
                                   double tnew,
                                   double lnacceptance,
                                   bool accepted);

};


#endif //G_PHOCS_LOCUSALLDATA_H
