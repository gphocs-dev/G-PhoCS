/** 
   \file PopulationTree.c 
   Functions to handle a population tree with migration bands.
	
   Contains the relevant data structures and procedure implementations
   for handling a population tree with migration bands.	
*/


#include "LocusEmbeddedGenealogy.h"
#include "PopulationTree.h"
#include "DataLayerConstants.h"

#include <stdlib.h>
#include "utils.h"
#include <math.h>
#include <algorithm>
#include "set"

/***************************************************************************************************************/
/******                                      INTERNAL CONSTANTS                                           ******/
/***************************************************************************************************************/



#define PERCISION	0.0000001

extern RandGeneratorContext RndCtx;


/***************************************************************************************************************/
/******                                     INTERNAL DATA TYPES                                           ******/
/***************************************************************************************************************/



/***************************************************************************************************************/
/******                                INTERNAL FUNCTION DECLARATIONS                                     ******/



/***************************************************************************************************************/
/******                              EXTERNAL FUNCTION IMPLEMENTATION                                     ******/
/***************************************************************************************************************/



/***********************************************************************************
 *	createPopTree
 *	- allocates basic memory for population tree (no migration bands yet)
 * 	- returns pointer to newly allocated population tree
 ***********************************************************************************/
PopulationTree* createPopTree(int numCurPops)
{
  int pop, migBand, numPops = 2*numCurPops-1;

  PopulationTree* popTree = (PopulationTree*) malloc( sizeof(PopulationTree) );
  if(popTree == nullptr)
  {
    fprintf(stderr, "\nError: Out Of Memory population tree.\n");
    exit(-1);
  }
	
  popTree->popArray = (Population*) malloc( numPops * sizeof(Population) );
  if(popTree->popArray == nullptr) {
    fprintf(stderr, "\nError: Out Of Memory population array in population tree.\n");
    exit(-1);
  }
	
  popTree->pops = (Population**) malloc( numPops * sizeof(Population*) );
  if(popTree->pops == nullptr) {
    fprintf(stderr, "\nError: Out Of Memory population array in population tree.\n");
    exit(-1);
  }
	
  popTree->isAncestralArray = (unsigned short*) malloc( numPops * numPops * sizeof(unsigned short));
  if(popTree->isAncestralArray == nullptr) {
    fprintf(stderr, "\nError: Out Of Memory boolean 2D array for isAncestrals in population tree.\n");
    exit(-1);
  }

  bzero( popTree->isAncestralArray,
         sizeof(unsigned short) * numPops * numPops );

  // allocate initial memory for migration bands
  popTree->migBands = (MigrationBand*)malloc(numPops*(numPops-1)*sizeof(MigrationBand));
  if(popTree->migBands == nullptr) {
    fprintf(stderr, "\nError: Out Of Memory migration band array in population tree.\n");
    exit(-1);
  }

  popTree->migBandIdArray = (int*)malloc(2*numPops*(numPops-1)*sizeof(int));
  if(popTree->migBandIdArray == nullptr) {
    fprintf(stderr, "\nError: Out Of Memory migration band id's array in population tree.\n");
    exit(-1);
  }

  // initialize populations
	
  popTree->numCurPops  = numCurPops;
  popTree->numPops 	   = numPops;
  popTree->numMigBands = 0;
	
  for( pop = 0; pop < numPops; ++pop )
  {
    popTree->pops[pop]                 = popTree->popArray + pop;
    popTree->pops[pop]->id             = pop;
    popTree->pops[pop]->inMigBands     = popTree->migBandIdArray
                                         + 2 * pop * (numPops-1);
    popTree->pops[pop]->outMigBands    = popTree->migBandIdArray
                                         + (2 * pop + 1) * (numPops - 1);
    popTree->pops[pop]->numInMigBands  = 0;
    popTree->pops[pop]->numOutMigBands = 0;
    popTree->pops[pop]->name[0]        = '\0';
    popTree->pops[pop]->isAncestralTo  = popTree->isAncestralArray
                                         + pop * numPops;

	  popTree->pops[pop]->age            = 0.0;
	  //MARK CHANGE
    popTree->pops[pop]->numSamples     = 0;
    popTree->pops[pop]->sampleAge      = 0.0;
    popTree->pops[pop]->theta          = 0.0;
    popTree->pops[pop]->father         = nullptr;
    popTree->pops[pop]->sons[0]        = nullptr;
    popTree->pops[pop]->sons[1]        = nullptr;
    popTree->pops[pop]->thetaPrior.alpha       = 0.0;
    popTree->pops[pop]->thetaPrior.beta        = 0.0;
    popTree->pops[pop]->thetaPrior.sampleStart = 0.0;
    popTree->pops[pop]->agePrior.alpha         = 0.0;
    popTree->pops[pop]->agePrior.beta          = 0.0;
    popTree->pops[pop]->agePrior.sampleStart   = 0.0;
	
  }

  for( migBand = 0; migBand < numPops * (numPops-1); ++migBand )
  {
    popTree->migBands[migBand].id                       = migBand;
    popTree->migBands[migBand].sourcePop                = 0;
    popTree->migBands[migBand].targetPop                = 0;
    popTree->migBands[migBand].migRate                  = 0.0;
    popTree->migBands[migBand].startTime                = 0.0;
    popTree->migBands[migBand].endTime                  = 0.0;
    popTree->migBands[migBand].migRatePrior.alpha       = 0.0;
    popTree->migBands[migBand].migRatePrior.beta        = 0.0;
    popTree->migBands[migBand].migRatePrior.sampleStart = 0.0;
  }

  return popTree;
}


/***********************************************************************************
 *	freePopTree
 *	- frees all memory allocated for population tree (in createPopTree and initMigrationBands)
 * 	- returns 0
 ***********************************************************************************/
int freePopTree(PopulationTree* popTree) {

  free(popTree->migBandIdArray);
  free(popTree->migBands);
  free(popTree->isAncestralArray);
  free(popTree->pops);
  free(popTree->popArray);
  free(popTree);

  return 0;
}



/***********************************************************************************
*	populationPostOrder
*	Computes pots-order for population subtreetree rooted at pop.
   Writes down the post order in specified array.
   Returns size of subtree.
   A recursive procedure.
***********************************************************************************/
int popPostOrder(PopulationTree *popTree, int pop, int *ordered_pops) {
    int size;

    // halting condition
    if (pop < popTree->numCurPops) {
        ordered_pops[0] = pop;
        return 1;
    }

    // compute post-order for every subtree and add root
    size = popPostOrder(popTree, popTree->pops[pop]->sons[0]->id,
                        ordered_pops);
    size += popPostOrder(popTree, popTree->pops[pop]->sons[1]->id,
                         ordered_pops + size);
    ordered_pops[size] = pop;

    return size + 1;
}

/* populationPostOrder - vector version

int AllLoci::populationPostOrder(int pop, std::vector<int>::iterator it) {
    int size;
    // halting condition
    if (pop < dataSetup.popTree->numCurPops) {
        *it = pop;
        return 1;
    }
    // compute post-order for every subtree and add root
    size = populationPostOrder(dataSetup.popTree->pops[pop]->sons[0]->id,
                               it);
    size += populationPostOrder(dataSetup.popTree->pops[pop]->sons[1]->id,
                                it + size);
    popQueue_[size] = pop;

    return size + 1;
}
*/


/***********************************************************************************
 *	printPopulationTree
 *	- prints population tree
 ***********************************************************************************/
void printPopulationTree(PopulationTree* popTree, FILE* stream, int printTauTheta)	{
  int pop, pop1, migBand;
  size_t maxNameLen=0;
  char formatStr[100];
	
  if (printTauTheta != 0 || printTauTheta == 1) //If user specified something other than boolean don't print tau & theta
    printTauTheta = 0;

  fprintf(stream, "---------------------------------------------------------------\n");
  fprintf(stream, "Current populations:\n");
  for(pop=0; pop<popTree->numPops; pop++) {
    if(strlen(popTree->pops[pop]->name) > maxNameLen)
      maxNameLen = strlen(popTree->pops[pop]->name);
  }
  sprintf(formatStr, " pop %%2d (%%%us), ", maxNameLen);
  if (printTauTheta == 1)
    sprintf(formatStr, "%stau [%%7f], theta[%%7f], ", formatStr);
  for(pop=0; pop<popTree->numPops; pop++) {
    if(pop == popTree->numCurPops) {
      fprintf(stream, "Ancestral populations:\n");
    }

    fprintf(stream, formatStr, pop, popTree->pops[pop]->name, popTree->pops[pop]->age, popTree->pops[pop]->theta);
					
    if(pop == popTree->rootPop)		fprintf(stream, "   ROOT    , ");
    else							fprintf(stream, "father (%2d), ", popTree->pops[pop]->father->id);
		
    if(pop < popTree->numCurPops)	fprintf(stream, "            , ");
    else							fprintf(stream, "sons (%2d %2d), ", popTree->pops[pop]->sons[0]->id, popTree->pops[pop]->sons[1]->id);
		
    fprintf(stream, "is ancestral to:");
    for(pop1=0; pop1<popTree->numPops; pop1++) {
      if(popTree->pops[pop]->isAncestralTo[pop1])		fprintf(stream, " %2d",pop1);
    }
    if(popTree->pops[pop]->numInMigBands > 0) {
      fprintf(stream, ". Incoming migration bands:");
      for(migBand=0; migBand<popTree->pops[pop]->numInMigBands; migBand++) {
        fprintf(stream, " %2d", migBand<popTree->pops[pop]->inMigBands[migBand]);
      }
    }
    if(popTree->pops[pop]->numOutMigBands > 0) {
      fprintf(stream, ". Outgoing migration bands:");
      for(migBand=0; migBand<popTree->pops[pop]->numOutMigBands; migBand++) {
        fprintf(stream, " %2d", migBand<popTree->pops[pop]->outMigBands[migBand]);
      }
    }
    fprintf(stream, ".\n");
  }
  fprintf(stream, "---------------------------------------------------------------\n");
  if(popTree->numMigBands > 0) {
    fprintf(stream, "%d Migration Bands:\n", popTree->numMigBands);
    for(migBand=0; migBand<popTree->numMigBands; migBand++) {
      fprintf(stream, " mig-band %2d, [%2d --> %2d], mig-rate [%7f], times [%7f - %7f].\n",
             migBand, popTree->migBands[migBand].sourcePop, popTree->migBands[migBand].targetPop,
             popTree->migBands[migBand].migRate, popTree->migBands[migBand].startTime, popTree->migBands[migBand].endTime);
    }
    fprintf(stream, "---------------------------------------------------------------\n");
  }

  return;
}
/** end of printPopulationTree **/



/***********************************************************************************
 *	getPopIdByName
 * 	- returns a population id of a population given its name (-1 if no match is found)
 * 	- used primarily to decode migration bands as specified in control file 
 * 		(called by readControlFile).
 ***********************************************************************************/
int getPopIdByName(PopulationTree* popTree, const char* name){
  int pop;

  for(pop=0; pop<popTree->numPops; pop++) {
    if(0 == strcmp(name,popTree->pops[pop]->name)) {
      return pop;
    }		
  }

  return -1;

}
/** end of getPopIdByName **/



/******************************************************************************
 *	samplePopParameters
 *	- samples population parameters according to prior average
 *	  (only thetas and taus)
 * 	- each parameter is sampled uniformly in the interval [0.9,1.1]*mean
 * 		(where mean is the prior mean for that parameter)
 * 	- makes sure a population's age does not exceed its father's
 *	- initializes all migration rates to 0.
 * 	- returns 0
 *****************************************************************************/
/* MARK: CHECK SAMPLE AGE OF SON POPULATION IF NOT ZERO AND CORRECT !!
         DONE, NEED TO CHECK !!
*/
int samplePopParameters(PopulationTree* popTree)
{
  int head, tail, migBand;
  double mean;
  Population* pop;
  Population** popQueue;
	
	
  popQueue = (Population**) malloc(popTree->numPops * sizeof(Population*));
  if( nullptr == popQueue )
  {
    fprintf(stderr,
            "\nError: Out Of Memory popQueue in samplePopParameters.\n");
    exit(-1);
  }
	
  head = tail = 0;
  popQueue[tail++] = popTree->pops[ popTree->rootPop ];
	

  // traverse population pre-order (from root down)
  // sample parameters for each population.
  // make corrections for population age, if greater than father's age	
  while( head < tail )
  {
    pop = popQueue[head++];
    //		printf("Pop %d\n", pop->id);
    // sample theta
    //		mean =  pop->thetaPrior.alpha / pop->thetaPrior.beta;
    mean =  pop->thetaPrior.sampleStart;
    pop->theta = mean * ( 0.9 + 0.2*rndu( RAND_GENERAL_SLOT ) );
    if( nullptr != pop->sons[0] )
    {
      // sample age for ancestral population
      //			mean =  pop->agePrior.alpha / pop->agePrior.beta;
      mean =  pop->agePrior.sampleStart;
      pop->age = mean * ( 0.9 + 0.2*rndu( RAND_GENERAL_SLOT ) );
      //			printf("-Tau = %f\n", pop->age);
      // if inconsistent with father's age, resample within 93% and 97% of father's age
      if(    pop->father != nullptr
          && pop->father->age < pop->age)
      {
        // first drop age to max age of sons. 
        // should be zero unless sons associated with ancient samples
        pop->age = max2( pop->sons[0]->sampleAge,
                         pop->sons[1]->sampleAge );
		    pop->age +=   (pop->father->age - pop->age)
                    * (0.93 + 0.004*rndu( RAND_GENERAL_SLOT ) );
      }
      // add sons to end of queue
      popQueue[tail++] = pop->sons[0];
      popQueue[tail++] = pop->sons[1];
			
    }
  }

  free(popQueue);
	
  // initialize migration rates to 0!
  for( migBand = 0; migBand < popTree->numMigBands; ++migBand )
  {
    popTree->migBands[migBand].migRate = 0.0;
  }
	
  computeMigrationBandTimes(popTree);
	
  return 0;
}
/** end of samplePopParameters **/



/******************************************************************************
 *	sampleMigRates
 *	- samples migration rates for all mig bands
 * 	- each rate is sampled uniformly in the interval [0.9,1.1]*mean
 * 		(where mean is the prior mean for that parameter)
 * 	- returns 0
 *****************************************************************************/
int sampleMigRates(PopulationTree* popTree)
{
  int migBand;
  double mean;
	
  // sample migration rates
  for( migBand = 0; migBand<popTree->numMigBands; ++migBand )
  {
    mean =   popTree->migBands[migBand].migRatePrior.alpha
           / popTree->migBands[migBand].migRatePrior.beta;
    popTree->migBands[migBand].migRate = mean *
                                         (0.9 +
                                          0.2 * rndu( RAND_GENERAL_SLOT ) );
    //		printf("Mig band %d, rate %f\n",migBand, popTree->migBands[migBand].migRate);
  }
	
  return 0;
	
}
/** end of sampleMigRates **/



/***********************************************************************************
 *	updateMigrationBandTimes
 *	- updates satrt and end times of given migration band according to ages of populations
 *	- returns 0 if no change was made, and 1 otherwise
 *	- if change was made, re-construct live mig-bands data structure
 ***********************************************************************************/
unsigned short updateMigrationBandTimes(PopulationTree* popTree, int migBand) {
  int sourcePop, targetPop; 
	
  unsigned short res = 0;
  double newTime;
	
  sourcePop = popTree->migBands[migBand].sourcePop;
  targetPop = popTree->migBands[migBand].targetPop;
	
  newTime = max2( popTree->pops[sourcePop]->age , popTree->pops[targetPop]->age );
  if(newTime != popTree->migBands[migBand].startTime) {
    popTree->migBands[migBand].startTime = newTime;
    res = 1;
  }

  newTime =  min2( popTree->pops[sourcePop]->father->age , popTree->pops[targetPop]->father->age );
  if(newTime != popTree->migBands[migBand].endTime) {
    popTree->migBands[migBand].endTime = newTime;
    res = 1;
  }

  //		printf("Setting migration band %d for times %f -- %f.\n",migBand, popTree->migBands[migBand].startTime, popTree->migBands[migBand].endTime);

  //if change was made, re-construct live mig-bands data structure
  //if (res) //todo call it here? and only if res?
      //constructMigBandsTimes(popTree);

  return res;
}
/** end of updateMigrationBandTimes **/



/***********************************************************************************
 *	computeMigrationBandTimes
 *	- traverses all migration bands and sets their start and end times according to 
 *		times of target and source populations.
 * 	- returns the number of migration bands with zero span
 ***********************************************************************************/
int computeMigrationBandTimes(PopulationTree* popTree) {
  int migBand, numZero = 0;

    for(migBand=0; migBand<popTree->numMigBands; migBand++) {
		
    updateMigrationBandTimes(popTree, migBand);

    if(popTree->migBands[migBand].startTime >= popTree->migBands[migBand].endTime) {
      numZero++;
      popTree->migBands[migBand].startTime = popTree->migBands[migBand].endTime = popTree->pops[popTree->migBands[migBand].targetPop]->age;
    }

  }

  return numZero;
}
/** end of computeMigrationBandTimes **/



/***************************************************************************************************************/
/******                              INTERNAL FUNCTION IMPLEMENTATION                                     ******/
/***************************************************************************************************************/



/***********************************************************************************
*	getMigBandByPops
*	- returns pointer to mig band with the given source and target population
* 	- returns pointer to migration band
***********************************************************************************/
MigrationBand* getMigBandByPops(PopulationTree* popTree, int sourcePop, int targetPop){

    //traverse all mig bands and initialize all incoming migration bands
    for(int migBand=0; migBand < popTree->numMigBands; migBand++) {

        if (popTree->migBands[migBand].sourcePop == sourcePop &&
            popTree->migBands[migBand].targetPop == targetPop)
            return &popTree->migBands[migBand];
    }
    return nullptr;
}


/***********************************************************************************
*	getMigBandById
*	- returns pointer to a mig band with the given ID
***********************************************************************************/
MigrationBand * getMigBandByID(PopulationTree* popTree, int id) {

    if (id < popTree->numMigBands)
        return &popTree->migBands[id];
    return nullptr;
}


/*******************************************************************************
 *	initializeLivingMigBands
 *	create N elements in livingMigBands vector, where N is num of pops
 *	divide migration bands into groups with same target pop
 ******************************************************************************/
void initializeMigBandTimes(PopulationTree* popTree) {

    //fill vector with N empty elements
    popTree->migBandsPerTarget.reserve(popTree->numPops);
    for (int i=0; i < popTree->numPops; i++) {
        popTree->migBandsPerTarget.emplace_back();
    }

    //for each mig band, save it in the MigBandsPerTarget element located in the
    // i'th place, where i is the id of the mig band's target pop
    for (int i=0; i < popTree->numMigBands; i++) {

        //get target pop of current mig band
        int targetPop = popTree->migBands[i].targetPop;

        //add pointer
        MigrationBand* pMigBand = &popTree->migBands[i];
        popTree->migBandsPerTarget[targetPop].pMigBands.push_back(pMigBand);

        //add mig band ID
        int id = pMigBand->id;
        popTree->migBandsPerTarget[targetPop].migBandsIDs.push_back(id);
    }
}


/*******************************************************************************
 *	constructMigBandsTimes
 *  constructs mig bands times of each pop
 *
 ******************************************************************************/
void constructMigBandsTimes(PopulationTree* popTree) {

    //for each pop reset its time bands
    for (auto& popBands : popTree->migBandsPerTarget) {
        popBands.timeMigBands.clear();
    }

    //for each target pop
    for (int pop=0; pop < popTree->numPops; pop++) {

        //get bands of current pop
        auto & popBands = popTree->migBandsPerTarget[pop];

        //create a vector of time points and fill with pop start and end,
        // and the start and end times of each mig band
        std::vector<double> timePoints;
        timePoints.reserve(2*popBands.pMigBands.size() + 2);

        //pop times
        timePoints.push_back(popTree->pops[pop]->age);
        if (pop != popTree->rootPop)
            timePoints.push_back(popTree->pops[pop]->father->age);
        else
            timePoints.push_back(OLDAGE);//todo: what age?

        //mig bands times
        for (MigrationBand* pMigBand : popBands.pMigBands) {
            timePoints.push_back(pMigBand->startTime);
            timePoints.push_back(pMigBand->endTime);
        }

        //sort time points and remove duplicates
        sort(timePoints.begin(), timePoints.end());
        timePoints.erase(unique(timePoints.begin(), timePoints.end()),
                         timePoints.end());

        //for each two successive time points
        for (auto it=timePoints.begin(); it!=timePoints.end()-1; ++it) {

            //get next time
            auto next = std::next(it,1);

            //create a time-band element which lasts between the two time points
            TimeMigBands timeMigBand{*it, *next};

            //find mig-bands which intersect with current time band
            for (MigrationBand* pMigBand : popBands.pMigBands) {

                //if time-band is contained in current mig-band
                if (pMigBand->startTime <= timeMigBand.startTime &&
                        timeMigBand.endTime <= pMigBand->endTime )

                    //add mig band
                    timeMigBand.migBands.push_back(pMigBand);
            }

            //add time-band to pop's time-bands
            popBands.timeMigBands.push_back(timeMigBand);

        }
    }
}


/*******************************************************************************
 *	getLiveMigBands
 *	for the given target pop, returns a time band containing the given age,
 *	and null if not found such.
    @param: popTree, target pop, age
    @return: pointer to a time band struct
 ******************************************************************************/
TimeMigBands *
getLiveMigBands(PopulationTree* popTree, int target_pop, double age) {

    for (auto& timeBand : popTree->migBandsPerTarget[target_pop].timeMigBands) {
        if (timeBand.startTime <= age && age < timeBand.endTime)
            return &timeBand;
    }
    return nullptr;
}


/*
 * getPopLeaves
 * get all leaves of a population
*/
std::vector<int> & getPopLeaves(PopulationTree* popTree, int pop) {
    return popTree->popToLeaves[pop];
}


/*
 * getLeafPop
 * get population of a leaf
*/
int getLeafPop(PopulationTree* popTree, int leafId) {
    return popTree->leafToPop[leafId];
}


/*
 * printPopToLeaves
 * print pop to leaves map
*/
void printPopToLeaves(PopulationTree* popTree) {

    std::cout << "Pop to leaves: " << std::endl;
    for (auto &x : popTree->popToLeaves) {
        std::cout << "pop " << x.first << ", leaves: ";
        for (int leaf : x.second) {
            std::cout << leaf << ", ";
        }
        std::cout << std::endl;
    }
}


/*
 * printLeafToPop
 * print leaf to map
*/
void printLeafToPop(PopulationTree* popTree) {
    std::cout << "Leaf to pop: " << std::endl;
    for (auto &x : popTree->leafToPop) {
        std::cout << "leaf " << x.first << ", "
             << "pop: " << x.second << std::endl;
    }
}
/***************************************************************************************************************/
/******                                        END OF FILE                                                ******/
/***************************************************************************************************************/


