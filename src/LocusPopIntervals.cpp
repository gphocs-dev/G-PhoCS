
#include "LocusPopIntervals.h"
#include "DbgErrMsgIntervals.h"
#include "GPhoCS.h"
#include "MemoryMng.h"
#include "TreeNode.h"

#include <iostream>


/*
    LocusPopIntervals constructor
    Allocates popIntervals objects.
    Links intervals to each other.
*/
LocusPopIntervals::LocusPopIntervals(int locusID, int nIntervals)
        : locusID_(locusID),
          numIntervals_(nIntervals),
          pPopTree_(dataSetup.popTree), //todo: get dataSetup as a pointer
          stats_(dataSetup.popTree->numPops, dataSetup.popTree->numMigBands) {

    //allocate N intervals (N = number of intervals, given as argument)
    intervalsArray_ = new PopInterval[nIntervals];

    //intervals pool points to head of intervals array
    pIntervalsPool_ = intervalsArray_;

}


/*
    LocusPopIntervals copy-constructor
    Allocates popIntervals objects.
    Links intervals to each other.
*/
LocusPopIntervals::LocusPopIntervals(const LocusPopIntervals& other) :
          numIntervals_(other.numIntervals_),
          stats_(other.stats_),
          locusID_(other.locusID_), pPopTree_(other.pPopTree_) {

    //allocate N intervals (N = number of intervals, given as argument)
    intervalsArray_ = new PopInterval[numIntervals_];

    //intervals pool points to head of intervals array
    pIntervalsPool_ = intervalsArray_;

}


/*
 * LocusPopIntervals
 * LocusPopIntervals class destructor
*/
LocusPopIntervals::~LocusPopIntervals() {
    //delete array of intervals
    delete intervalsArray_;
}


/*
 * reset
 * Links intervals to each other and reset intervals
*/
void LocusPopIntervals::resetPopIntervals() {
    //reset all intervals
    for (int i = 0; i < numIntervals_; ++i) {
        intervalsArray_[i].resetPopInterval();
    }
}


/*
    linkIntervals
    Links intervals to each other and reset intervals
*/
void LocusPopIntervals::linkIntervals() {
    //for each interval, link to prev and next intervals
    //set pointer to tree node to null
    for (int i = 0; i < numIntervals_; ++i) {
        if (i == 0) { //first interval
            intervalsArray_[i].setNext(intervalsArray_ + i + 1);
            intervalsArray_[i].setPrev(nullptr);
        } else if (i == numIntervals_ - 1) { //last interval
            intervalsArray_[i].setNext(nullptr);
            intervalsArray_[i].setPrev(intervalsArray_ + i - 1);
        } else {
            intervalsArray_[i].setNext(intervalsArray_ + i + 1); //next cell
            intervalsArray_[i].setPrev(intervalsArray_ + i - 1); //prev cell
        }
    }
}


/*
    getIntervalFromPool
    Returns a free interval from the intervals pool
    @return: pointer to a free interval
*/
PopInterval* LocusPopIntervals::getIntervalFromPool() {

    PopInterval* pInterval = pIntervalsPool_;

    //if no free intervals
    if (!pInterval) {
        INTERVALS_FATAL_0015
    }

    //update pointer to intervals pool
    pIntervalsPool_ = pIntervalsPool_->getNext();
    pIntervalsPool_->setPrev(nullptr);

    pInterval->setNext(nullptr);

    return pInterval;
}


/*
    returnToPool
    Return an interval to the interval pool
    @param: pointer to the free interval
*/
void LocusPopIntervals::returnToPool(PopInterval* pInterval) {

    //reset interval content
    pInterval->resetPopInterval();

    //set pointers
    pIntervalsPool_->setPrev(pInterval);
    pInterval->setNext(pIntervalsPool_);
    pIntervalsPool_ = pInterval;
}


/*
    createStartEndIntervals
    Initializes intervals array with pop-start and pop-end intervals
    by defining:
     1. The first N cells to be pop-start intervals (N = num pops).
     2. The next N cells to be pop-end intervals.
    Total: 2N cells are occupied with start/end intervals.
*/
void LocusPopIntervals::createStartEndIntervals() {

    int nPops = pPopTree_->numPops; //N
    int rootPop = pPopTree_->rootPop; //root population

    //for each population define a start and end intervals
    for (int pop = 0; pop < nPops; ++pop) {

        int start_index = pop;
        int end_index = pop + nPops;

        //start interval

        //set the i-th cell (i=pop) to be a pop-start interval
        intervalsArray_[start_index].setType(IntervalType::POP_START);

        //set pop
        intervalsArray_[start_index].setPopID(pop);

        //set num lineages - it is important to initialize 0 incoming lineages
        intervalsArray_[start_index].setNumLineages(0);

        //next of start interval points to the end interval of same pop
        intervalsArray_[start_index].setNext(intervalsArray_ + end_index);

        //prev of start interval points to NULL
        intervalsArray_[start_index].setPrev(nullptr);

        //set age, which is population age
        intervalsArray_[start_index].setAge(pPopTree_->pops[pop]->age);


        //end interval

        //set the (i+N)-th cell (i=pop) to be an end interval
        intervalsArray_[end_index].setType(IntervalType::POP_END);

        //set pop
        intervalsArray_[end_index].setPopID(pop);

        //set num lineages - it is important to initialize 0 incoming lineages
        intervalsArray_[end_index].setNumLineages(0);

        //set age, which is the parent population age
        // distinguish between rootPop and rest of pops
        if (pop == rootPop)
            intervalsArray_[end_index].setAge(OLDAGE);
        else
            intervalsArray_[end_index].setAge(pPopTree_->pops[pop]->father->age);


        //next of end interval points to the start interval of parent pop.
        //next of root's interval points to NULL
        if (pop == rootPop) {
            intervalsArray_[end_index].setNext(nullptr);
        }
        else {
            int parent_pop = pPopTree_->pops[pop]->father->id; //get parent pop
            intervalsArray_[end_index].setNext(intervalsArray_ + parent_pop);
        }

        //prev of end interval points to the start interval of same pop
        intervalsArray_[end_index].setPrev(intervalsArray_ + start_index);
    }

    //promote the free intervals pointer to the 2N cell
    pIntervalsPool_ = intervalsArray_ + 2*nPops;
    pIntervalsPool_->setPrev(nullptr);

}


/*
    createIntervalBefore
    Creates an interval before a specified interval
    @param: pointer to an interval, age of the new interval,
            and (optional) interval type
    @return: pointer to the new interval
*/
PopInterval *
LocusPopIntervals::createIntervalBefore(PopInterval* pInterval, int pop,
                                        double age, IntervalType type) {

    //get a new interval from pool (default type is DUMMY)
    PopInterval* pNewInterval = this->getIntervalFromPool();

    //set population ID
    pNewInterval->setPopID(pop);

    //set interval type
    pNewInterval->setType(type);

    //set number of lineages of new interval to same number as given interval
    pNewInterval->setNumLineages(pInterval->getNumLineages());

    //set age of new interval
    pNewInterval->setAge(age);

    //set pointers
    pInterval->getPrev()->setNext(pNewInterval);
    pNewInterval->setPrev(pInterval->getPrev());
    pNewInterval->setNext(pInterval);
    pInterval->setPrev(pNewInterval);

    return pNewInterval;
}


/*
   createInterval
   Creates a new interval in specified population at given time.
   Num lineages of new interval is as the one
   of subsequent interval.
   @param: population id, age, interval type
   @return: pointer to new interval
*/
PopInterval *
LocusPopIntervals::createInterval(int pop, double age, IntervalType type) {

    //check if time specified is smaller than age of target population
    if (age < pPopTree_->pops[pop]->age) {
        INTERVALS_FATAL_0016
        return nullptr;
    }

    //check if time specified is greater than age of parent population
    if(pop != pPopTree_->rootPop &&
       age > pPopTree_->pops[pop]->father->age + 0.000001) {
        INTERVALS_FATAL_0017
        return nullptr;
    }

    //find a spot for a new interval
    //loop while not reaching the end interval of the population
    //and while age of current interval is smaller than time specified
    PopInterval* pInterval = this->getPopStart(pop)->getNext();
    while (!pInterval->isType(IntervalType::POP_END) && pInterval->getAge() < age) {
           pInterval = pInterval->getNext();
    }

    //create the new interval in the found slot
    return this->createIntervalBefore(pInterval, pop, age, type);

}


/*
    getPopStart
    Returns the pop-start interval of a specified population.
    Pop-start interval is allocated in the n-th cell of the events array,
    where n = population id
    @param: population id
    @return: pointer to a pop-start interval

*/
PopInterval* LocusPopIntervals::getPopStart(int pop) {
    return intervalsArray_ + pop;
}


/*
    getPopEnd
    Returns the pop-end interval of a specified population.
    Pop-end interval is allocated in the n-th cell of the events array,
    where n = num populations + population id
    @param: population id
    @return: pointer to a pop-end interval

*/
PopInterval* LocusPopIntervals::getPopEnd(int pop) {
    return intervalsArray_ + pPopTree_->numPops + pop;
}


/*
    getSamplesStart
    returns the samplesStart interval of a population
    @param: population id
    @return: pointer to a SamplesStart interval
*/
PopInterval* LocusPopIntervals::getSamplesStart(int pop) {

    //iterate over intervals while not pop-end
    PopInterval* pInterval = this->getPopStart(pop)->getNext();
    while (!pInterval->isType(IntervalType::POP_END)) {

        //if interval is a sample start - return it
        if (pInterval->isType(IntervalType::SAMPLES_START))
            return pInterval;

        //move to next interval
        pInterval = pInterval->getNext();
    }
    return nullptr;

}


/*
    getStats
    @return: reference to statistics object
*/
const GenealogyStats& LocusPopIntervals::getStats() const {
    return stats_;
}


/*
    printIntervals
*/
void LocusPopIntervals::printIntervals() {
    std::cout << "Intervals of locus " << locusID_ << "." << std::endl;

    //for each population iterate over all intervals
    for (int pop = 0; pop < pPopTree_->numPops; ++pop) {
        PopInterval* pInterval = this->getPopStart(pop);
        //iterate intervals of current pop while not pop-end
        while (true) {
            //print interval
            pInterval->printInterval();
            if (pInterval->isType(IntervalType::POP_END))
                break;
            pInterval = pInterval->getNext();
        }
    }
}


/*	computeGenetreeStats
	Computes the statistics of a given locus.
    Assumes intervals chains are built, but number of lineages is ONLY
    set for first intervals in the leaf populations.
	Sets number of lineages for each non-leaf event by traversing the
	population tree post-order. In parallel, also records the statistics.
*/
int LocusPopIntervals::computeGenetreeStats() {
    // go over all intervals and compute num_lineages per each interval
    // also update genetree statistics
    for (int i = 0; i < pPopTree_->numPops; i++) {

        //get current pop
        int pop = pPopTree_->popsPostOrder[i];

        // if not leaf population get number of in-lineages from end-intervals
        // of son populations
        if (pop >= pPopTree_->numCurPops) {

            //get num lineages of first interval of the left pop son
            int lSon = pPopTree_->pops[pop]->sons[0]->id;
            int n1 = this->getPopEnd(lSon)->getNumLineages();

            //get num lineages of first interval of the right pop son
            int rSon = pPopTree_->pops[pop]->sons[1]->id;
            int n2 = this->getPopEnd(rSon)->getNumLineages();

            //set num lineages of current pop-start to be the sum of sons' num
            // lineages
            this->getPopStart(pop)->setNumLineages(n1 + n2);

        } else {
            //set num lineages of first interval of current pop to 0
            this->getPopStart(pop)->setNumLineages(0);
        }

        //recalculate statistics
        this->recalcStats(pop);
    }

    return 0;
}


/* recalcStats
   Re-calculates stats for a given population
*/
double LocusPopIntervals::recalcStats(int pop) {

    //define local copy of statistics

    //coal statistics (a single instance for current pop)
    GenStats coalStats;

    //create a map of mig statistics, with mig-band id as a key
    //get only mig-bands which their target pop equal current pop
    std::map<int, GenStats> migsStats; //todo: replace dynamic allocation
    for (auto pMigBand: pPopTree_->migBandsPerTarget[pop].migBands) {
        migsStats[pMigBand->id] = GenStats();
    }

    //get pop-start interval
    PopInterval* pInterval = this->getPopStart(pop);

    //get current age
    double currAge = pInterval->getAge();

    //get num lineages of first interval
    int n = pInterval->getNumLineages();

    //get live mig bands
    TimeMigBands* timeBand = getLiveMigBands(dataSetup.popTree, pop, currAge);

    // follow intervals chain and set number of lineages per interval according
    // to previous interval also update statistics
    while (true) {

        double t = min2(pInterval->getAge(), timeBand->endTime) - currAge;

        //increment coal statistics
        coalStats.stats +=  n * (n - 1) * t; //todo: more efficient calculation

        //for each live mig band update mig statistics
        for (auto pMigBand : timeBand->migBands) {
            migsStats[pMigBand->id].stats += n*t;
        }

        //update current age
        currAge += t;

        //if interval age is larger than end of time band - get next time band
        if (timeBand->endTime <= pInterval->getAge()) {
            double prevEndTime = timeBand->endTime;
            timeBand = getLiveMigBands(dataSetup.popTree, pop, currAge);

            if (!timeBand) {
                assert(pInterval->isType(IntervalType::POP_END));
                assert(prevEndTime == pInterval->getAge());
                break;
            }
            continue;
        }

        //assert
        assert(!pInterval->isType(IntervalType::POP_END));

        //switch by interval type
        switch (pInterval->getType()) {

            case (IntervalType::SAMPLES_START): {
                n += dataSetup.numSamplesPerPop[pop];
                break;
            }
            case (IntervalType::COAL): {
                n--;
                coalStats.num += 1;
                break;
            }
            case (IntervalType::IN_MIG): {
                n--;
                // figure out migration band and update its statistics
                MigNode* pMigNode = (MigNode*)pInterval->getTreeNode();
                int mig_band = pMigNode->getMigBandId();
                migsStats[mig_band].num += 1;
                break;
            }
            case (IntervalType::OUT_MIG): {
                n++;
                break;
            }
            default:
                break;

        }// end of switch

        //get next interval
        pInterval = pInterval->getNext();

        //set num lineages
        pInterval->setNumLineages(n);

    }// end of while

    double HEREDITY_FACTOR = 1;
    double deltaLnLd = 0.0;

    //for each mig statistics
    for (auto key_value : migsStats) {
        int id = key_value.first; //mig band id

        //update delta ln likelihood of mig statistics
        deltaLnLd -= (migsStats[id].stats - stats_.migs[id].stats) *
                      pPopTree_->migBands[id].migRate;

        //save mig statistics (in the pop location)
        stats_.migs[id].num = migsStats[id].num;
        stats_.migs[id].stats = migsStats[id].stats;
    }

    //update delta ln likelihood of coal statistics
    deltaLnLd -= (coalStats.stats - stats_.coals[pop].stats) /
                  (pPopTree_->pops[pop]->theta * HEREDITY_FACTOR);

    //save coal statistics (in the pop location)
    stats_.coals[pop].num = coalStats.num;
    stats_.coals[pop].stats = coalStats.stats;

    return deltaLnLd;
}


/* testPopIntervals
   test if the new events data structure is consistent with the original
*/
void LocusPopIntervals::testPopIntervals() {

    //for each pop
    for (int pop = 0; pop < pPopTree_->numPops; pop++) {

        //get first event in old structure
        int event = event_chains[locusID_].first_event[pop];

        //get pop-start interval
        PopInterval* pInterval = this->getPopStart(pop);

        //get pop age
        double eventAge = pPopTree_->pops[pop]->age;

        //assert ages equal
        assert(fabs(eventAge - pInterval->getAge()) < EPSILON);

        //vector of live mig bands (original data structure)
        std::vector<int> liveMigsOri;

        //get next interval
        pInterval = pInterval->getNext();

        //iterate both old and new structures (events VS intervals)
        for (event; event >= 0;
             event = event_chains[locusID_].events[event].getNextIdx()) {

            EventType eventType = event_chains[locusID_].events[event].getType();
            double elapsedTime = event_chains[locusID_].events[event].getElapsedTime();

            //num lineages
            //verify that nums lineages are equal
            int eventLin = event_chains[locusID_].events[event].getNumLineages();
            int intervalLin = pInterval->getNumLineages();
            assert(eventLin == intervalLin);

            //mig bands
            //if elapsed time of event is greater than epsilon, compare live mig bands
            if (elapsedTime > 2 * EPSILON) {

                double age = eventAge + elapsedTime/2;
                TimeMigBands* liveMigsNew = getLiveMigBands(pPopTree_, pop, age);

                //assert live migs band is not null
                assert(liveMigsNew != nullptr);

                //compare num of live mig bands
                assert(liveMigsOri.size() == liveMigsNew->migBands.size());

                //for each live mig band
                for (int id : liveMigsOri) {
                    //get pointer to mig band by its id
                    MigrationBand* migBand = getMigBandByID(pPopTree_, id);
                    //verify that current mig is found in the new data structure
                    assert(std::find(liveMigsNew->migBands.begin(),
                                     liveMigsNew->migBands.end(), migBand) !=
                           liveMigsNew->migBands.end());
                }

                //verify if time band contains old event
                assert(liveMigsNew->startTime <= eventAge + EPSILON);
                assert(eventAge + elapsedTime <= liveMigsNew->endTime + EPSILON);
            }

            //age
            //add elapsed time to event age
            eventAge += elapsedTime;

            //event type
            //get event type and verify that interval is of same type
            // or, in case of mig band start/end, add/remove a mig band

            switch (eventType) {
                case SAMPLES_START: {
                    assert(pInterval->isType(IntervalType::SAMPLES_START)); break;
                }
                case COAL: {
                    assert(pInterval->isType(IntervalType::COAL)); break;
                }
                case IN_MIG: {
                    assert(pInterval->isType(IntervalType::IN_MIG)); break;
                }
                case OUT_MIG: {
                    assert(pInterval->isType(IntervalType::OUT_MIG)); break;
                }
                case END_CHAIN: {
                    assert(pInterval->isType(IntervalType::POP_END)); break;
                }

                    //mig bands (don't exist in the new data structure)
                case MIG_BAND_START: {
                    //add id of mig band to live mig bands
                    int event_id = event_chains[locusID_].events[event].getId();
                    liveMigsOri.push_back(event_id);
                    continue;
                }
                case MIG_BAND_END: {
                    //remove id of mig band from live mig bands
                    int event_id = event_chains[locusID_].events[event].getId();
                    liveMigsOri.erase(
                            std::remove(liveMigsOri.begin(), liveMigsOri.end(),
                                        event_id), liveMigsOri.end());
                    continue;
                }
            }//end of switch

            //age
            //verify that ages are equal
            double intervalAge = pInterval->getAge();
            assert(fabs(eventAge - intervalAge) < EPSILON);

            //get next interval
            pInterval = pInterval->getNext();
        }

    }//end of pop loop
}


/* testGenealogyStatistics
  verify statistics are equal to statistics of old data structure
*/
void LocusPopIntervals::testGenealogyStatistics() {

    //for each pop
    for (int pop = 0; pop < pPopTree_->numPops; pop++) {

        //verify num of coal are equal
        assert(stats_.coals[pop].num == genetree_stats[locusID_].num_coals[pop]);

        //verify statistics are equal
        assert(fabs(stats_.coals[pop].stats - genetree_stats[locusID_].coal_stats[pop]) < EPSILON);
    }

    //for each mig band
    for (int id = 0; id < pPopTree_->numMigBands; id++) {

        //verify num of migs are equal
        assert(stats_.migs[id].num == genetree_stats[locusID_].num_migs[id]);

        //verify statistics are equal
        assert(fabs(stats_.migs[id].stats - genetree_stats[locusID_].mig_stats[id]) < EPSILON);
    }
}
