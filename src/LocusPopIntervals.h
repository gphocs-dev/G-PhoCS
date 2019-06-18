//
// Created by nomihadar on 3/11/19.
//

#ifndef G_PHOCS_LOCUSPOPINTERVALS_H
#define G_PHOCS_LOCUSPOPINTERVALS_H

#include "PopInterval.h"
#include "MemoryMng.h"
#include "DataLayer.h"
#include "GenealogyStats.h"

#include <vector>

/*=============================================================================
 *
 * LocusPopIntervals class
 *
 * Class LocusPopIntervals is the structure containing the interval chains.
 * Each item in this structure will be of PopInterval class.
 * This replaces the current event chains, with the main difference being that
 * the chains of the daughter populations connect with the parent and
 * we add pop-start and pop-end intervals (which are empty - elapsed_time = 0).
 *
 * Contains:
 * 1. Array if intervals objects.
 * 2. Pointer to a pool of free intervals
 * 3. Total number of intervals
 * 4. Statistics vector of coal/migs for each pop / mig-band
 * 5. Pointer to popTree.
 *===========================================================================*/


class LocusPopIntervals {

private:

    PopInterval* intervalsArray_; //array of PopIntervals //todo: replace with vector?
    PopInterval* pIntervalsPool_; //pointer to a pool of free intervals

    int numIntervals_; //total number of intervals

    GenealogyStats stats_; //genealogy statistics of coal and migs

    int locusID_; //locus id, for error massages
    PopulationTree* pPopTree_; //pointer to PopulationTree struct

public:

    //constructor
    LocusPopIntervals(int locusID, int nIntervals);

    //copy-constructor
    LocusPopIntervals(const LocusPopIntervals& other);

    //destructor
    ~LocusPopIntervals();

    //reset intervals
    void resetPopIntervals();

    //link intervals to each other
    void linkIntervals();

    //initialize intervals array with pop-start and pop-end intervals
    void createStartEndIntervals();

    //get a free interval from intervals pool
    PopInterval* getIntervalFromPool();

    //return a free interval to intervals pool
    void returnToPool(PopInterval* pInterval);

    //create a new interval before a given interval
    PopInterval* createIntervalBefore(PopInterval* pInterval,
            int pop, double age, IntervalType type);

    //create a new interval in specified population at given time
    PopInterval* createInterval(int pop, double age, IntervalType type);

    //get a pointer to a samples start interval of a given population
    PopInterval* getSamplesStart(int pop);

    //get a pointer to a pop start interval of a given population
    PopInterval* getPopStart(int pop);

    //get a pointer to a end start interval of a given population
    PopInterval* getPopEnd(int pop);

    //get a reference to statistics
    const GenealogyStats& getStats() const;

    //compute genealogy tree statistics
    int computeGenetreeStats();

    //recalculate statistics
    double recalcStats(int pop);

    //for each pop print all intervals
    void printIntervals();

    //verify intervals are equal to old data structure
    void testPopIntervals();

    //verify statistics are equal to statistics of old data structure
    void testGenealogyStatistics();

};


#endif //G_PHOCS_LOCUSPOPINTERVALS_H
