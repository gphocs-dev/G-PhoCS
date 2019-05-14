
#include "LocusPopIntervals.h"
#include "DbgErrMsgIntervals.h"
#include "GPhoCS.h"
#include "MemoryMng.h"

#include <iostream>

/*
    LocusPopIntervals constructor
    Allocates popIntervals objects.
    Links intervals to each other.
*/
LocusPopIntervals::LocusPopIntervals(int locusID, int nIntervals,
                                     PopulationTree* pPopTree)
        : locusID_(locusID),
          numIntervals_(nIntervals),
          pPopTree_(pPopTree) {

    //allocate N intervals (N = number of intervals, given as argument)
    intervalsArray_ = new PopInterval[nIntervals];

    //intervals pool points to head of intervals array
    pIntervalsPool_ = intervalsArray_;
    
}

/*
    LocusPopIntervals class destructor
*/
LocusPopIntervals::~LocusPopIntervals() {
    //delete array of intervals
    delete intervalsArray_;
}

/*
    Links intervals to each other and reset intervals
*/
void LocusPopIntervals::reset() {
    //reset all intervals
    for (int i = 0; i < numIntervals_; ++i) {
        intervalsArray_[i].reset();
    }
}

/*
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
    Return an interval to the interval pool
    @param: pointer to the free interval
*/
void LocusPopIntervals::returnToPool(PopInterval* pInterval) {

    //reset interval content
    pInterval->reset();

    //set pointers
    pIntervalsPool_->setPrev(pInterval);
    pInterval->setNext(pIntervalsPool_);
    pIntervalsPool_ = pInterval;
}


/*
    Initializes intervals array with pop-start and pop-end intervals
    by defining:
     1. The first N cells to be pop-start intervals (N = num pops).
     2. The next N cells to be pop-end intervals.
    Total: 2N cells are occupied to start/end intervals.
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

        //end interval

        //set the (i+N)-th cell (i=pop) to be an end interval
        intervalsArray_[end_index].setType(IntervalType::POP_END);

        //set pop
        intervalsArray_[end_index].setPopID(pop);

        //set num lineages - it is important to initialize 0 incoming lineages
        intervalsArray_[end_index].setNumLineages(0);

        //set age, distinguish between rootPop and rest of pops
        if (pop == rootPop) {
            double t = OLDAGE - pPopTree_->pops[rootPop]->age;
            intervalsArray_[end_index].setAge(OLDAGE);
        }
        else {
            double t1 = pPopTree_->pops[pop]->father->age;
            double t2 = pPopTree_->pops[pop]->age;
            intervalsArray_[end_index].setAge(t1);//todo: yes?
        }

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
    Creates an interval before a specified interval
    @param: pointer to an interval, age of the new interval,
            and (optional) interval type
    @return: pointer to the new interval
*/
PopInterval *
LocusPopIntervals::createIntervalBefore(PopInterval *pInterval, int pop,
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
    PopInterval* pInterval = this->getFirstInterval(pop);
    while (!pInterval->isType(IntervalType::POP_END) && pInterval->getAge() > age) {
           pInterval = pInterval->getNext();
    }

    if (pInterval->isType(IntervalType::SAMPLES_START)) {
        int i =2;
    }
    cout<<"&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"<<pInterval->typeToStr()<<endl;

    //if (pInterval->getElapsedTime() < delta_time) {
     //   if (pInterval->getElapsedTime() < (delta_time - 0.000001)) {
     //       INTERVALS_FATAL_0018
     //   }
     //   delta_time = pInterval->getElapsedTime();
    //}

    //create the new interval in the found slot
    return this->createIntervalBefore(pInterval, pop, age, type);

}


/*
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
    Returns the first interval (which is not the pop start)
    of a specified population.
    @param: population id
    @return: pointer of pop interval

*/
PopInterval* LocusPopIntervals::getFirstInterval(int pop) {
    return this->getPopStart(pop)->getNext();
}


/*
    returns the samplesStart interval of a population
    @param: population id
    @return: pointer to a SamplesStart interval
*/
PopInterval* LocusPopIntervals::getSamplesStart(int pop) {

    //iterate over intervals while not pop-end
    PopInterval* pInterval = this->getFirstInterval(pop);
    while (!pInterval->isType(IntervalType::POP_END)) {

        //if interval is a sample start - return it
        if (pInterval->isType(IntervalType::SAMPLES_START))
            return pInterval;

        //move to next interval
        pInterval = pInterval->getNext();
    }
    return nullptr;

}


void LocusPopIntervals::printIntervals() {

    std::cout << "Intervals of locus " << locusID_ << "." << std::endl;

    //for each population iterate over all intervals
    for (int pop = 0; pop < pPopTree_->numPops; ++pop) {

        //cout << "Intervals of pop " << pop << ":" << endl;

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
