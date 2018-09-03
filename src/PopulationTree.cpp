/** 
   \file PopulationTree.c 
   Functions to handle a population tree with migration bands.
	
   Contains the relevant data structures and procedure implementations
   for handling a population tree with migration bands.	
*/

#include <stdlib.h>
#include "PopulationTree.h"
#include "utils.h"
#include <math.h>
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



/****  migration band (sets) manipulation  ****/
MigrationBandSet* createMigBandSet(PopulationTree* popTree, int targetPop, double age);
int moveMigBandSource(PopulationTree* popTree, int migBand, double newAge, unsigned short startORend);



/***************************************************************************************************************/
/******                              EXTERNAL FUNCTION IMPLEMENTATION                                     ******/
/***************************************************************************************************************/



/***********************************************************************************
 *	createPopTree
 *	- allocates basic memory for population tree (no migration bands yet)
 * 	- returns pointer to newly allocated population tree
 ***********************************************************************************/
PopulationTree* createPopTree(int numCurPops)	{
  int pop, migBand, numPops = 2*numCurPops-1;
	
  PopulationTree* popTree = (PopulationTree*)malloc(sizeof(PopulationTree));
  if(popTree == NULL) {
    fprintf(stderr, "\nError: Out Of Memory population tree.\n");
    exit(-1);
  }
	
  popTree->popArray = (Population*)malloc(numPops*sizeof(Population));
  if(popTree->popArray == NULL) {
    fprintf(stderr, "\nError: Out Of Memory population array in population tree.\n");
    exit(-1);
  }
	
  popTree->pops = (Population**)malloc(numPops*sizeof(Population*));
  if(popTree->pops == NULL) {
    fprintf(stderr, "\nError: Out Of Memory population array in population tree.\n");
    exit(-1);
  }
	
  popTree->isAncestralArray = (unsigned short*)malloc(numPops*numPops*sizeof(unsigned short));
  if(popTree->isAncestralArray == NULL) {
    fprintf(stderr, "\nError: Out Of Memory boolean 2D array for isAncestrals in population tree.\n");
    exit(-1);
  }
  resetBooleanArray(popTree->isAncestralArray,numPops*numPops);

  // allocate initial memory for migration bands
  popTree->migBands = (MigrationBand*)malloc(numPops*(numPops-1)*sizeof(MigrationBand));
  if(popTree->migBands == NULL) {
    fprintf(stderr, "\nError: Out Of Memory migration band array in population tree.\n");
    exit(-1);
  }

  popTree->migBandIdArray = (int*)malloc(2*numPops*(numPops-1)*sizeof(int));
  if(popTree->migBandIdArray == NULL) {
    fprintf(stderr, "\nError: Out Of Memory migration band id's array in population tree.\n");
    exit(-1);
  }

  // initialize populations
	
  popTree->numCurPops  = numCurPops;
  popTree->numPops 	 = numPops;
  popTree->numMigBands = 0;
	
  for(pop=0; pop<numPops; pop++) {
    popTree->pops[pop] = popTree->popArray + pop;
    popTree->pops[pop]->id = pop;
    popTree->pops[pop]->inMigBands  = popTree->migBandIdArray + 2*pop*(numPops-1);
    popTree->pops[pop]->outMigBands = popTree->migBandIdArray + (2*pop+1)*(numPops-1);
    popTree->pops[pop]->numInMigBands  = 0;
    popTree->pops[pop]->numOutMigBands = 0;
    popTree->pops[pop]->name[0] = '\0';
    popTree->pops[pop]->isAncestralTo = popTree->isAncestralArray + pop*numPops;

	popTree->pops[pop]->age = 0.0;
	//MARK CHANGE
	popTree->pops[pop]->numSamples = 0;
	popTree->pops[pop]->sampleAge = 0.0;
	popTree->pops[pop]->theta = 0.0;
	popTree->pops[pop]->father = NULL;
	popTree->pops[pop]->sons[0] = NULL;
	popTree->pops[pop]->sons[1] = NULL;
	popTree->pops[pop]->thetaPrior.alpha = 0.0;
	popTree->pops[pop]->thetaPrior.beta = 0.0;
	popTree->pops[pop]->thetaPrior.sampleStart = 0.0;
	popTree->pops[pop]->agePrior.alpha = 0.0;
	popTree->pops[pop]->agePrior.beta = 0.0;
	popTree->pops[pop]->agePrior.sampleStart = 0.0;
	
  }

  for(migBand=0; migBand<numPops*(numPops-1); migBand++) {
    popTree->migBands[migBand].sourcePop = 0;
    popTree->migBands[migBand].targetPop = 0;
    popTree->migBands[migBand].migRate = 0.0;
    popTree->migBands[migBand].startTime = 0.0;
    popTree->migBands[migBand].endTime = 0.0;
    popTree->migBands[migBand].firstSet = NULL;
    popTree->migBands[migBand].lastSet = NULL;
    popTree->migBands[migBand].migRatePrior.alpha = 0.0;
    popTree->migBands[migBand].migRatePrior.beta = 0.0;
    popTree->migBands[migBand].migRatePrior.sampleStart = 0.0;
  }
  return popTree;
}



/***********************************************************************************
 *	initMigrationBands
 *	- initializes data structures for migration bands in population tree (including allocating some memory)
 * 	- sets start and end times
 * 	- for each population creates a timed-sequence of migration band sets
 * 	- returns 0
 ***********************************************************************************/
int initMigrationBands(PopulationTree* popTree) {
  int i, sourcePop, targetPop, migBand;
  MigrationBandSet *migSet;
  double startTime, endTime;
	
  // allocate memory for migration band sets and mig-band ids per pop
  popTree->migBandSetArray = (MigrationBandSet*)malloc((2*popTree->numMigBands + popTree->numPops)*sizeof(MigrationBandSet));
  if(popTree->migBandSetArray == NULL) {
    fprintf(stderr, "\nError: Out Of Memory migration band set array in population tree.\n");
    exit(-1);
  }
	
	
  // initialize stack of mig band sets
  popTree->migBandSetArray[0].prev = NULL;
  popTree->migBandSetArray[2*popTree->numMigBands + popTree->numPops - 1].next = NULL;
  for(i=0; i<2*popTree->numMigBands + popTree->numPops - 1; i++) {
    popTree->migBandSetArray[i].next = &(popTree->migBandSetArray[i+1]);
    popTree->migBandSetArray[i+1].prev = &(popTree->migBandSetArray[i]);
  }
  popTree->migBandSetStackTop = popTree->migBandSetArray;
	
  // initialize sequence of migration bands with an empty set for each population
  for(targetPop=0; targetPop<popTree->numPops; targetPop++) {
    migSet = popTree->migBandSetStackTop;
    popTree->migBandSetStackTop = migSet->next;
    migSet->next = migSet->prev = NULL;
    migSet->numMigBands = 0;
    migSet->rate = 0.0;
    migSet->age = -1.0;
    popTree->pops[targetPop]->migBandSequence = migSet;
  }
	
  // traverse all populations and initialize all incoming migration bands
  for(migBand=0; migBand<popTree->numMigBands; migBand++) {
    sourcePop = popTree->migBands[migBand].sourcePop;
    targetPop = popTree->migBands[migBand].targetPop;
    popTree->pops[sourcePop]->outMigBands[ ++popTree->pops[sourcePop]->numOutMigBands ] = migBand;
    popTree->pops[targetPop]->inMigBands [ ++popTree->pops[targetPop]->numInMigBands  ] = migBand;
		
    // start and end time are set according to source pop
    // in mig-band itself they are recorded as intersection with target pop
    startTime = popTree->pops[sourcePop]->age;
    endTime   = popTree->pops[sourcePop]->father->age;
    popTree->migBands[migBand].startTime = max2(startTime, 	popTree->pops[targetPop]->age);
    popTree->migBands[migBand].endTime   = min2(endTime, 	popTree->pops[sourcePop]->father->age);
		
    migSet = createMigBandSet(popTree, targetPop, startTime);
    popTree->migBands[migBand].firstSet = migSet;
    migSet = createMigBandSet(popTree, targetPop, endTime);
    popTree->migBands[migBand].lastSet  = migSet;
    // add index of migration band to migBandIds[] arrays of all 
    // sets between start point and end point. Also adjust rates
    for(migSet = popTree->migBands[migBand].firstSet ; migSet != popTree->migBands[migBand].lastSet; migSet = migSet->next) {
      migSet->migBandIds[ migSet->numMigBands++ ] = migBand;
      migSet->rate += popTree->migBands[migBand].migRate;
    }
  }
	
  return 0;
}
/** end of initMigrationBands **/



/***********************************************************************************
 *	freePopTree
 *	- frees all memory allocated for population tree (in createPopTree and initMigrationBands)
 * 	- returns 0
 ***********************************************************************************/
int freePopTree(PopulationTree* popTree) {
	
  free(popTree->migBandSetArray);	
  free(popTree->migBandIdArray);
  free(popTree->migBands);
  free(popTree->isAncestralArray);
  free(popTree->pops);
  free(popTree->popArray);
  free(popTree);

  return 0;
}



/***********************************************************************************
 *	printPopulationTree
 *	- prints population tree
 ***********************************************************************************/
void printPopulationTree(PopulationTree* popTree, FILE* stream, int printTauTheta)	{
  int pop, pop1, migBand;
  size_t maxNameLen=0;
  char formatStr[100];
	
  if ((!(printTauTheta == 0) || (printTauTheta == 1))) //If user specified something other than boolean don't print tau & theta
    printTauTheta = 0;

  fprintf(stream, "---------------------------------------------------------------\n");
  fprintf(stream, "Current populations:\n");
  for(pop=0; pop<popTree->numPops; pop++) {
    if(strlen(popTree->pops[pop]->name) > maxNameLen)
      maxNameLen = strlen(popTree->pops[pop]->name);
  }
  // eug:
  // Resolving the warning: 
  // format ‘%u’ expects argument of type ‘unsigned int’, 
  // but argument 3 has type ‘size_t {aka long unsigned int}’ [-Wformat=]
  // Solution - change the formatter from "u" to "zu".
  // Possible problem - corrupted reports
  sprintf(formatStr, " pop %%2d (%%%zus), ", maxNameLen);
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
  if( NULL == popQueue )
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
    mean =  pop->thetaPrior.sampleStart;
    pop->theta = mean * ( 0.9 + 0.2*rndu( RAND_GENERAL_SLOT ) );
    if( NULL != pop->sons[0] )
    {
      // sample age for ancestral population
      mean =  pop->agePrior.sampleStart;
      pop->age = mean * ( 0.9 + 0.2*rndu( RAND_GENERAL_SLOT ) );
      //			printf("-Tau = %f\n", pop->age);
      // if inconsistent with father's age, resample within 93% and 97% of father's age
      if(    pop->father != NULL
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
 *	createMigBandSet
 *	- creates a migration band set in the sequence of a given target population at a give age
 * 	- if age is equal (within PRECISION) to the age of another set in the sequence,
 * 		then returns that set. Otherwise takes a new set from the stack and copies 
 * 		all data from previous set in sequence.
 * 	- returns pointer to the new set
 ***********************************************************************************/
MigrationBandSet* createMigBandSet(PopulationTree* popTree, int targetPop, double age) {
  int i;
  MigrationBandSet *migSet, *migSetNew;
	
  // find spot for new mig-band set	
  for(migSet = popTree->pops[targetPop]->migBandSequence; 
      migSet->next != NULL && migSet->next->age < age;
      migSet = migSet->next) { ; }

  if(migSet->next != NULL && migSet->next->age - age < PERCISION) {
    // no need to create new set - use old one
    return migSet->next;
  }
  // take new set from stack and copy all data (except for age)
  migSetNew = popTree->migBandSetStackTop;
  popTree->migBandSetStackTop = migSetNew->next;
  migSetNew->next = migSet->next;
  migSetNew->prev = migSet;
  migSetNew->numMigBands = migSet->numMigBands;
  migSetNew->rate = migSet->rate;
  migSetNew->age = age;
  for(i=0; i<migSetNew->numMigBands; i++) {
    migSetNew->migBandIds[i] = migSet->migBandIds[i];
  }
  migSet->next = migSetNew;
  if(migSetNew->next != NULL) {
    migSetNew->next->prev = migSetNew;
  }
  return migSetNew;
}
/** end of createMigBandSet **/



/***********************************************************************************
 *	moveMigBandSource
 *	- moves start/end time of source population for specific migration band
 * 	- if startORend == 0, moves end time, otherwise (startORend == 1), moves start time.
 * 	- returns 0
 ***********************************************************************************/
int moveMigBandSource(PopulationTree* popTree, int migBand, double newAge, unsigned short startORend) {
  int i, targetPop;
  unsigned short addORremove;
  MigrationBandSet *migSet, *migSetBottom, *migSetTop;
  double oldAge;
	
  targetPop = popTree->migBands[migBand].targetPop;
  migSet = createMigBandSet(popTree, targetPop, newAge);
  // consider all 4 possible scenarios
  if(startORend) {
    oldAge = popTree->migBands[migBand].startTime;
    if(oldAge < newAge) {
      // moving start time up - need to remove mig band from sets
      addORremove = 0;
      migSetTop = migSet;
      migSetBottom = popTree->migBands[migBand].firstSet;
    } else {
      // moving start time down - need to add mig band from sets
      addORremove = 1;
      migSetBottom = migSet;
      migSetTop = popTree->migBands[migBand].firstSet;
    }
  } else {
    oldAge = popTree->migBands[migBand].endTime;
    if(oldAge < newAge) {
      // moving end time up - need to add mig band from sets
      addORremove = 1;
      migSetTop = migSet;
      migSetBottom = popTree->migBands[migBand].lastSet->next;
    } else {
      // moving end time down - need to remove mig band from sets
      addORremove = 0;
      migSetBottom = migSet;
      migSetTop = popTree->migBands[migBand].lastSet->next;
    }
		
  }

  if(addORremove) {
    // add migBand to every set between bottom and top (excluding)
    for( migSet=migSetBottom; migSet->next != migSetTop; migSet = migSet->next) {
      migSet->migBandIds[ migSet->numMigBands++ ] = migBand;
      migSet->rate += popTree->migBands[migBand].migRate;
    }
  } else {
    // remove migBand from every set between bottom and top (excluding)
    for( migSet=migSetBottom; migSet->next != migSetTop; migSet = migSet->next) {
      for(i=0; i<migSet->numMigBands; i++) {
        if(migSet->migBandIds[i] == migBand)	break;
      }
      if(i >= migSet->numMigBands) {
        fprintf(stderr, "\nError: moveMigBandSource: moving migration band %d %s from %g to %g.\n",
               migBand, (startORend) ? ("start time") : ("end time"), oldAge, newAge);
        fprintf(stderr, "Could not find migration band in set at time %g\n", migSet->age);
        printPopulationTree(popTree, stderr, 1);
        exit(-1);
      }
      migSet->migBandIds[i] = migSet->migBandIds[--migSet->numMigBands];
      migSet->rate -= popTree->migBands[migBand].migRate;
    }
  }

  return 0;
}
/** end of moveMigBandSource **/



/***************************************************************************************************************/
/******                                        END OF FILE                                                ******/
/***************************************************************************************************************/
