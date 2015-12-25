/** 
   \file GenealogyInPopTree.c 
   Geneaology data structures, functions for generating proposed changes and tracking changes in genealogy, event manipulation, printing genealogies.
	
   Contains the relevant data structures and procedure implementations
   for holding detailed information about a genealogy embedded in a population tree.
   This information is stored in form of a chain of intervals.
   Each interval resides in a specific population.
   Transition between one interval to the next is defined either by a change
   in number of lineages or by a change in migration bands.
	
   NO MIGRATION YET
   NO CHANGE TO STRUCTURE OF TREE YET

*/


#include "utils.h"
#include "PopulationTree.h"
#include "GenealogyInPopTree.h"


/***************************************************************************************************************/
/******                                      INTERNAL CONSTANTS                                           ******/
/***************************************************************************************************************/




/***************************************************************************************************************/
/******                                              DATA TYPES                                           ******/
/***************************************************************************************************************/

/** EventType
    Defining different types of intervals.
*/
typedef enum EVENT_TYPE {
  COAL, 				/**< coalescent event */
  IN_MIG, 				/**< incoming migration event	(forward-time view) */
  OUT_MIG, 				/**< outgoing migration event	(forward-time view) */
  MIG_BAND_START,		/**< indication that a migration band is activated at this point */
  MIG_BAND_END,			/**< indication that a migration band is deactivated at this point */
  POP_START,	 		/**< bottom event for every population */
  POP_END,				/**< top event for every population */
  LEAF,					/**< event for leaf of genealogy (leaves are kept below POP_START events) */
  DUMMY					/**< dummy event - assigned to events which are modified by a proposed genealogy change */
} EventType;



/**  Event
     Data type which corresponds to a specific interval

	- Each interval points to next and previous intervals in linear chain
		- mostly next pointer is used
		- POP_END event points up to POP_START event of father population
		- POP_START events do not point to prev
	- Each COAL, IN/OUT_MIG interval indicates the GenBranch on which it resides
	- Each COAL, IN/OUT_MIG interval points to father in genealogy which is
		next event along GenBranch on which it resides
	- Each COAL interval points to two children in genealogy
		(2 previous events along 2 GenBranches below coalescent node)
	- Each IN/OUT_MIG interval points to a single child in genealogy
		(previous event along GenBranch on which it resides)
*/
typedef struct EVENT Event;

/** Event

    Data tpe which corresponds to a specific interval
    @note For COAL events: id of branch above event 
	@note For IN_MIG events: id of migration band 
    @note For MIG_BAND_START/MIG_BAND_END the id of the migration band
*/
struct EVENT {							
  EventType type;			    /**< type of interval */
  Event *next;                  /**< pointer to next interval in chain */
  Event *prev;		            /**< pointers to  previous interval in chain */
  Event* father;			    /**< pointer to father interval in genealogy (for COAL IN/OUT_MIG) */
  Event* sons[2];			    /**< pointers to child(ren) interval in genealogy (for COAL IN/OUT_MIG) */
  double age;				    /**< end time of interval */
  int popId;				    /**< id of population where interval resides */
  MigrationBandSet* migBands;	/**< pointer to migration band set for this interval */
  int numLineages;		        /**< number of lineages crossing interval */
  int numLineages_bak;	        /**< number of lineages before proposed change (kept in negative format - see changeNumLineages()) */
  int auxId;				    /**< auxiliary id: */
};								



/**	GenBranch
    Data type which corresponds to a branch in the locus genealogy
    @note For each branch we hold information on the node below it:
    - Event for that node (COAL or LEAF)
    - LikelihoodNode id for the purpose of data likelihood calculations
*/
typedef struct GEN_BRANCH{
  Event*	bottomEvent;		/**< pointer to event below branch */
  int		likelihoodNodeId;	/**< id of likelihood node */
} GenBranch;



/**	GenStats
    Data type which holds statistics for locus genealogy for quick computation of likelihood of genealogy (Pr[G_i | PopTree])
*/
typedef struct GEN_STATS{
  double* 	popCoalStats;			/**< array of coalescent stats - one for each population */
  double* 	migBandStats;			/**< array of migration  stats - one for each migration band */
  int*		popNumCoals;			/**< array holding number of coalescent events in each population */
  int*		migBandNumMigs;			/**< array holding number of migration  events in each migration band */
} GenStats;



/**	GenChanges
    Data type which holds changes made to genealogy due to proposed modification to accommodate quick acceptance/rejection of change
*/
typedef struct GEN_CHANGES{
  // changes in stats
  GenStats	newStats;				    /**< genealogy statistics after change */
  int 		numAffectedPops;		    /**< length of affectedPops[] array */
  int*		affectedPops;			    /**< array of id's of pops affected by change */
  unsigned short* isAffectedPop;		/**< boolean array indicating affected populations */
  int 		numAffectedMigBands;	    /**< length of affectedMigBands[] array */
  int*		affectedMigBands;		    /**< array of id's of migration bands affected by change */
  unsigned short* isAffectedMigBand;	/**< boolean array indicating affected migration bands */
	
  // general consequential changes
  int			numLinChangedEvents;	/**< length of linChangedEvents[] array */
  Event**		linChangedEvents;		/**< an array of pointers to events whose number of lineages changed */
	
  // moved events
  int			numMovedEvents;			/**< length of movedEventsOld[] and movedEventsOld[] arrays */
  Event**		movedEventsOld;			/**< array of pointers to events before move */
  Event**		movedEventsNew;			/**< array of pointers to events after move */
	 
  // replaced events (migration nodes and migration start/end signals) 
  int			numAddedEvents;			/**< length of addedEvents[] array */
  Event**		addedEvents;			/**< array of pointers to newly added events */
  int			numRemovedEvents;		/**< length of removedEvents[] array */
  Event**		removedEvents;			/**< array of pointers to events marked for removal */
	
  // specific data kept for SPR operations (Subtree Prune and Regraft)
  Event*		prunedCoalEvent;		/**< pruned coalescent event (at top of pruned edge) */
  Event*		newCoalEvent;			/**< point of regraft */
  Event*		eventBelowRegraft;		/**< event just below regraft event on target edge */
  unsigned short prunedChildOfEvent;	/**< indicating whether pruning was done on sides (0 or 1 of prunedCoalEvent) */
	
} GenChanges;



/***********************************************************************************
 *	LocusGenealogy
 *	- Data type which holds genealogy for a given locus
 *	- Holds array of events (only used for allocation):
 *		* Total number of events per locus is the following:
 *			2*NUM_POPS + 2*NUM_MIG_BANDS + NUM_COALS + 2*maxMigrations + SPARE  ,where SPARE = 4 + 2*maxMigrations
 *		* New events are taken from the stack, which should never be left empty !
 *	- Holds LocusData data structure for locus (see LocusDataLikelihood.h)
 *	- Holds LocusGenStats data structure.
 ***********************************************************************************/
struct LOCUS_GENEALOGY{
  // general locus-specific information
  int 		locusId;					// id (or number) of locus
  int 		numSamples;					// number of samples (leaves in genealogy)
  int 		numBranches;				// number of branches in genealogy = 2*numSamples-1
  int 		numMigNodes;				// number of migration nodes in genealogy (dynamic - up to maxMigrations)
  double		heredityFactor;				// heredity factor for locus
  PopulationTree* popTree;				// pointer to the population tree
  LocusData* 	locusData;					// pointer to LocusData for maintaining data and computing data likelihood
  GenStats 	genStats;					// GenStats structure for computing genealogy likelihood
  GenChanges	genChanges;					// GenChanges structure for maintaining a proposed change to genealogy 
  GenBranch*	genBranches;				// array of genealogy branches
  int			rootBranchId;				// id of rot branch
  Event**		popStartEvents;				// array of pointers to all POP_START
  Event*		eventStackTop;				// pointer to top event in free event stack
  Event*		eventArray;					// an array of event consisting of the total pool used for this locus
};



/***************************************************************************************************************/
/******                                  GLOBAL DATA STRUCTURES                                           ******/
/***************************************************************************************************************/

int maxMigrations;		// global maximum number of migration events (set by setMaxMigrations)

/***************************************************************************************************************/
/******                                INTERNAL FUNCTION DECLARATIONS                                     ******/
/***************************************************************************************************************/



/****  generating proposed changes in genealogy  ****/
int		pruneGenSubtree(LocusGenealogy* locusGen, int branchId);
Event*	findLineageInInterval(PopulationTree* popTree, Event* event, int targetNum);
int		migrateAndCoalesceLineage(LocusGenealogy* locusGen, int branchId);
int 	proposeEventMove(LocusGenealogy* locusGen, Event* event, double newAge);
void	changeLocusMutationRate(LocusGenealogy* locusGen, double newRate);


/****  tracking changes in genealogy  ****/
int 	addAffectedPop(LocusGenealogy* locusGen, int popId);
int 	addAffectedMigBand(LocusGenealogy* locusGen, int migBand);
int 	addAffectedMigBandSet(LocusGenealogy* locusGen, MigrationBandSet* migBandSet);
int 	changeNumLineages(LocusGenealogy* locusGen, Event* event, int deltaNumLins);
int		computeNewGenStats(LocusGenealogy* locusGen);
double 	computeDeltaLogLikelihood(LocusGenealogy* locusGen);
int		acceptGenChanges(LocusGenealogy* locusGen);
int		resetGenChanges(LocusGenealogy* locusGen);

/****  event manipulation  ****/
int 	initEvents(LocusGenealogy* locusGen, int numEvents);
int 	recycleEvent(LocusGenealogy* locusGen, Event* event);
Event* 	createEventBefore(LocusGenealogy* locusGen, Event* event, double age);
Event* 	createEventAtAge(LocusGenealogy* locusGen, Event* refEvent, double age);


/****  printing  ****/
void 	printEvent(Event* event);
void 	printGenealogy(LocusGenealogy* locusGen);


/***************************************************************************************************************/
/******                              EXTERNAL FUNCTION IMPLEMENTATION                                     ******/
/***************************************************************************************************************/



/***********************************************************************************
 *	initLocusGenealogy
 *	- allocates all memory resources for locus (other than ones for data)
 * 	- returns 0
 ***********************************************************************************/
int initLocusGenealogy(LocusGenealogy* locusGen, PopulationTree* popTree, int id, int numSamples) {
	
  int pop, migBand,
    numEdges = 2*numSamples-1,
    numPops  = locusGen->popTree->numPops,
    numMigBands = locusGen->popTree->numMigBands,
		
    // max number of events is 2*NUM_POPS + 2*NUM_MIG_BANDS + NUM_LEAFS + NUM_COALS + 2*maxMigrations + SPARE,
    // where SPARE = 4 + 2*maxMigrations
    numEvents = 4 + 2*numPops + 2*numMigBands + numEdges + 4*maxMigrations;
	
  int* intArray;
  double* doubleArray;
  unsigned short* boolArray;
  Event** eventPtrArray;

  initEvents(locusGen,numEvents);
	
	
  locusGen->genBranches = (GenBranch*)malloc(numEdges*sizeof(GenBranch));
  if(locusGen->genBranches == NULL) {
    fprintf(stderr, "\nError: Out Of Memory genBranches for locus %d.\n",locusGen->locusId);
    exit(-1);
  }

	
  // allocate array of event pointers for the following purposes:
  // 1. pop start event
  // 2. events whose numLineages have changed
  // 3. moved events (old and updated versions)
  // 4. removed and added events (migration nodes)
  eventPtrArray = (Event**)malloc((numPops + numEvents + 6*maxMigrations + 2)*sizeof(Event*));
  if(eventPtrArray == NULL) {
    fprintf(stderr, "\nError: Out Of Memory Event pointer array for locus %d.\n",locusGen->locusId);
    exit(-1);
  }
	
  locusGen->popStartEvents = eventPtrArray;
  locusGen->genChanges.linChangedEvents 	= eventPtrArray+numPops;
  locusGen->genChanges.movedEventsOld 	= eventPtrArray+numPops+numEvents;
  locusGen->genChanges.movedEventsNew 	= eventPtrArray+numPops+numEvents+maxMigrations+1;
  locusGen->genChanges.addedEvents	 	= eventPtrArray+numPops+numEvents+2*maxMigrations+2;
  locusGen->genChanges.removedEvents	 	= eventPtrArray+numPops+numEvents+4*maxMigrations+2;
	
	
  // allocate array of double for locus genealogy stats (current and changed)
  doubleArray = (double*)malloc(2*(numPops+numMigBands)*sizeof(double));
  if(doubleArray == NULL) {
    fprintf(stderr, "\nError: Out Of Memory array of doubles for locus %d.\n",locusGen->locusId);
    exit(-1);
  }
	
  locusGen->genStats->popCoalStats = doubleArray;
  locusGen->genStats->migBandStats = doubleArray+numPops;
  locusGen->genChanges.newStats->popCoalStats = doubleArray+numPops+numMigBands;
  locusGen->genChanges.newStats->migBandStats = doubleArray+2*numPops+numMigBands;
	
  // allocate array of integers for coal and mig counts and ids of affected pops and migbands
  intArray = (int*)malloc(3*(numPops+numMigBands)*sizeof(int));
  if(intArray == NULL) {
    fprintf(stderr, "\nError: Out Of Memory array of integers for locus %d.\n",locusGen->locusId);
    exit(-1);
  }
	
  locusGen->genStats->popNumCoals = intArray;
  locusGen->genStats->migBandNumMigs = intArray+numPops;
  locusGen->genChanges.newStats->popNumCoals = intArray+numPops+numMigBands;
  locusGen->genChanges.newStats->migBandNumMigs = intArray+2*numPops+numMigBands;
  locusGen->genChanges.affectedPops = intArray+2*(numPops+numMigBands);
  locusGen->genChanges.affectedMigBands = intArray+3*numPops+2*numMigBands;
	
  // allocate array of booleans for indicators of affected pops and migbands
  boolArray = (unsigned short*)malloc((numPops+numMigBands)*sizeof(unsigned short));
  if(boolArray == NULL) {
    fprintf(stderr, "\nError: Out Of Memory array of booleans for locus %d.\n",locusGen->locusId);
    exit(-1);
  }
	
  locusGen->genChanges.isAffectedPop = boolArray;
  locusGen->genChanges.isAffectedMigBand = boolArray+numPops;
	
  // set and initialize various fields
	
  locusGen->numBranches = numEdges;
  locusGen->numSamples  = numSamples;
  locusGen->locusId     = id;
  locusGen->popTree = popTree;
  locusGen->numMigNodes = 0;
  locusGen->heredityFactor = 1.0;
	
  locusGen->genChanges.numAffectedPops = 0;
  locusGen->genChanges.numAffectedMigBands = 0;
  locusGen->genChanges.numLinChangedEvents = 0;
  locusGen->genChanges.numMovedEvents = 0;
  locusGen->genChanges.numAddedEvents = 0;
  locusGen->genChanges.numRemovedEvents = 0;
  locusGen->genChanges.prunedCoalEvent = NULL;
  locusGen->genChanges.newCoalEvent = NULL;

  for(pop=0; pop<popTree->numPops, pop++) {
    locusGen->genStats.popNumCoals[pop] = 0;
    locusGen->genStats.popCoalStats[pop] = 0.0;
    locusGen->genChanges.newStats.popNumCoals[pop] = 0;
    locusGen->genChanges.newStats.popCoalStats[pop] = 0.0;
    locusGen->genChanges.isAffectedPop[pop] = 0;
  }
  for(migBand=0; migBand<popTree->numMigBands, migBand++) {
    locusGen->genStats.migBandNumMigs[migBand] = 0;
    locusGen->genStats.migBandStats[migBand] = 0.0;
    locusGen->genChanges.newStats.migBandNumMigs[migBand] = 0;
    locusGen->genChanges.newStats.migBandStats[migBand] = 0.0;
    locusGen->genChanges.isAffectedMigBand[migBand] = 0;
  }

  return 0;
}
/** end of initLocusGenealogy **/



/***********************************************************************************
 *	initLocusData
 *	- allocates all memory resources for locus data and initializes them
 * 	- receives pattern info to copy into locus data
 * 	- returns 0
 ***********************************************************************************/
int initLocusData(LocusGenealogy* locusGen, int numSeqPatterns, char** patternArray, int* patternCounts) {
	
  // allocate memory for locus data likelihood structure
  locusGen->locusData = createLocusData (locusGen->numSamples, numSeqPatterns, patternArray, patternCounts);
  if(locusGen->locusData == NULL) {
    printf(" for locus %d.\n",locusGen->locusId);
    exit(-1);
  }

  return 0;
}
/** end of initLocusData **/



/***********************************************************************************
 *	freeLocusMemeory
 *	- frees all memory resources and data structures for locus
 * 	- returns 0
 ***********************************************************************************/
int freeLocusMemeory(LocusGenealogy* locusGen) {
	
  free(locusGen->eventArray);
  free(locusGen->genBranches);
  free(locusGen->popStartEvents);
  free(locusGen->genStats->popCoalStats);
  free(locusGen->genStats->popNumCoals);
  free(locusGen->genChanges.isAffectedPop);
  freeLocusData(locusGen->locusData);
	
  return 0;
}
/** end of freeLocusMemory **/



/***********************************************************************************
 *	setMaxNumMigrations
 *	- sets global maximum for number of migrations in all loci
 * 	- should be called once before initializing any of the loci, and only then.
 ***********************************************************************************/
void setMaxNumMigrations(int maxNumMigrations) {
	
  maxMigrations = maxNumMigrations;	
  return;
}
/** end of setMaxNumMigrations **/



/***********************************************************************************
 *	addLeaf
 *	- adds a leaf node in genealogy in indicated population
 * 	- used when initializing a genealogy
 * 	- initializes branch for this leaf
 * 	- returns 0
 ***********************************************************************************/
int addLeaf(LocusGenealogy* locusGen, int leafId, int popId) {
  // does not use createEvent, because events are below POP_START
  Event* newEvent = locusGen->eventStackTop;
  if(newEvent == NULL) {
    fprintf(stderr, "\nError: Empty event pool in locus %d.\n", locusGen->locusId);
    printGenealogy(locusGen);
    exit(-1);
  }
  locusGen->eventStackTop = newEvent->next;
	
  newEvent->type = LEAF;
  newEvent->next = locusGen->popStartEvents[popId];
  locusGen->genBranches[leafId].bottomEvent = newEvent;
  locusGen->genBranches[leafId].likelihoodNodeId = leafId;
	
  return 0;
}
/** end of addLeaf **/



/***********************************************************************************
 *	generateRandomGenealogy
 *	- generates a random genealogy according to population tree parameters
 * 	- simulated lineages one-by-one from leaves up in the population tree until they 
 * 		coalesce with some other lineage.
 * 	- calls migrateAndCoalesceLineage in each iteration
 * 	- returns 0
 ***********************************************************************************/
int generateRandomGenealogy(LocusGenealogy* locusGen) {
  int leafId = 0;
  int targetBranchId, newNodeId;
  Event	*topEvent, *event;
  Event**	 eventArray;
  int numEvents, ev;
	
  // find top-most event for initialization
  for(topEvent = locusGen->popStartEvents[ locusGen->popTree->rootPop ];
      topEvent->next != NULL;
      topEvent=topEvent->next) {	;  }
	
  // freely coalesce first lineage along population tree (with possible migrations
  locusGen->numBranches = 1;		// signal for migrateAndCoalesceLineage()
  migrateAndCoalesceLineage(locusGen, leafId);
  computeNewGenStats(locusGen);
  // if there are any migration events along first edge, chain them
  if(0 < (numEvents = locusGen->genChanges.numAddedEvents)) {
    eventArray 	= locusGen->genChanges.addedEvents;
    eventArray[0]->sons[0] = locusGen->genBranches[leafId].bottomEvent;
    locusGen->genBranches[leafId].bottomEvent->father = eventArray[0];
    for(ev=0; ev<numEvents-1; ev++) {
      eventArray[ev]->father = eventArray[ev+1];
      eventArray[ev+1]->sons[0] = eventArray[ev];
      if(ev%2 == 0) {
        event->type = IN_MIG;
        locusGen->genStats.migBandNumMigs[ event->auxId ]++;					
      } else {
        event->type = OUT_MIG;					
      }
    }
    topEvent->sons[0] = eventArray[numEvents-1];
    eventArray[numEvents-1]->father = topEvent;		
    locusGen->genChanges.numAddedEvents = 0;
  } else{
    // attach new coal event to bottom event of branch
    topEvent->sons[0] = locusGen->genBranches[leafId].bottomEvent;
    locusGen->genBranches[leafId].bottomEvent->father = topEvent;		
  }
	
  // atatch rest of leaves one by one
  for(leafId=1; leafId<locusGen->numSamples; leafId++) {
    migrateAndCoalesceLineage(locusGen, leafId);
    for(event=locusGen->genChanges.eventBelowRegraft; 
        event->type == OUT_MIG;
        event = event->sons[0]->sons[0])   { ; }
		
    targetBranchId = event->auxId;
    newNodeId = attachLeaf(	locusGen->locusData, 
                            locusGen->genBranches[leafId].likelihoodNodeId, 
                            locusGen->genBranches[targetBranchId].likelihoodNodeId, 
                            locusGen->genChanges.newCoalEvent->age);
					
    computeNewGenStats(locusGen);
    performRegraft(locusGen, leafId, leafId+locusGen->numSamples-1);
    locusGen->genBranches[ leafId+locusGen->numSamples-1 ].likelihoodNodeId = newNodeId;
    if(topEvent->sons[0] == locusGen->genBranches[ leafId+locusGen->numSamples-1 ].bottomEvent) {
      locusGen->rootBranchId = leafId+locusGen->numSamples-1;
    }
  }
	
  return 0;
}
/** end of generateRandomGenealogy **/



/***********************************************************************************
 *	changeLocusMutationRate
 *	- changes the locus-specific mutation rate to a new value.
 ***********************************************************************************/
void changeLocusMutationRate(LocusGenealogy* locusGen, double newRate)	{
  setLocusMutationRate(locusGen->locusData, newRate);
  return;
}
/** end of generateRandomGenealogy **/



/***********************************************************************************
 *	perturbNodeAges
 *	- perturb ages of all nodes (COAL, MIG) without modifying structure of genealogy.
 *	- input finetune argument defines the size of interval in which age is perturbed
 *	- returns the proportion of successful proposals.
 *	Note: we consider also migration events along root branch.
 ***********************************************************************************/
double	perturbNodeAges(LocusGenealogy* locusGen, double finetune) {
  int acceptCount, totalCount;
  double logAcceptanceRatio;
  double ageLB, ageUB, newAge;
  int branchId;
  unsigned short isMigEvent;
  Event* event;
	
  for(branchId=0; branchId<locusGen->numBranches; branchId++) {
    event=locusGen->genBranches[branchId].bottomEvent;
    if(event->type == LEAF) {
      // for leaf edges, perturb all migration events
      event = event->father;
      if(event->type == COAL)		continue;
    }

    // consider first COAL event at base of branch
    // and then traverse through all migration events
    while(1) {
      // slightly different operation for MIG o COAL events
      isMigEvent = (event->type == IN_MIG);
      totalCount++;
      // define bounds for perturbation according to father and children
      // and migration band times or population times
      if(isMigEvent) {
        ageLB = max2(event->sons[0]->age, locusGen->popTree->migBands[event->auxId].startTime);
        ageUB = min2(event->father->father->age,  locusGen->popTree->migBands[event->auxId].endTime);
      } else {
        ageLB = max2(event->sons[0]->age, event->sons[1]->age);
        ageLB = max2(ageLB, locusGen->popTree->pops[event->popId]->age);
        if(event->popId == locusGen->popTree->rootPop) {
          ageUB = event->father->age;
        } else {
          ageUB = min2(event->father->age, locusGen->popTree->pops[event->popId]->father->age);			
        }
      }
			
      // sample new age
      newAge = event->age + finetune*rnd2normal8();
      newAge = reflect(newAge, ageLB, ageUB);
			
      // move event, and if migration event, change also outgoing event
      proposeEventMove(locusGen, event, newAge);
      if(isMigEvent) {
        proposeEventMove(locusGen, event->father, newAge);
      }
			
      computeNewGenStats(locusGen);
      logAcceptanceRatio = computeDeltaLogLikelihood(locusGen);

      if(!isMigEvent) {
        adjustGenNodeAge(locusGen->locusData, locusGen->genBranches[branchId].likelihoodNodeId, newAge);
        logAcceptanceRatio += computeLocusDataLikelihood(locusGen->locusData, 1);
      }

      if(logAcceptanceRatio >= 0 || rndu() < exp(logAcceptanceRatio)) {
        // accept genealogy change
        acceptCount++;
        acceptGenChanges(locusGen);				
        resetSaved(locusGen->locusData);
      } else {
        // reject genealogy change
        resetGenChanges(locusGen);
        revertToSaved(locusGen->locusData);
      }
			
      // move on to next event
      if(isMigEvent) {
        event = event->father->father;
      } else {
        event = event->father;
      }
			
      // stop when next event is not migration
      if(event->type != IN_MIG)		break;
    }// end of while(1)

    if(event->type != COAL) {
      fprintf(stderr, "\nError: perturbNodeAges: locus %d, branch %d\n", locusGen->locusId, branchId);
      fprintf(stderr, "When traversing up along branch, reached the following event:\n");
      printEvent(event);
      printGenealogy(locusGen);
      exit(-1);
    }

  }// end for(branchId)

  return ((double)acceptCount)/totalCount;
}
/** end of perturbNodeAges **/

	
		
/***********************************************************************************
 *	pruneAndReCoalesceSubtrees
 *	- traverses all branches (excluding one above root), prunes the subtree below edge
 *		and re-coalesces the lineage from the root of the pruned subtree back into 
 *		remaining genealogy, according to model parameters (pop sizes and mig rates).
 *	- returns the proportion of successful proposals.
 ***********************************************************************************/
double	pruneAndReCoalesceSubtrees(LocusGenealogy* locusGen) {
  int acceptCount, totalCount;
  double logAcceptanceRatio;
  int branchId, targetBranchId;
  int res;
  Event* event;
	
  acceptCount = totalCount = 0;
	
  for(branchId=0; branchId<locusGen->numBranches; branchId++) {
    if(branchId == locusGen->rootBranchId)
      continue;
    totalCount++;
    pruneGenSubtree(locusGen, branchId);
    res = migrateAndCoalesceLineage(locusGen, branchId);
		
    if(res==0) {
      // if re-coalescence was successful, implement changes in locusGen->locusData
      //and compute delta log-likelihood of data
			
      // traverse down all migration events in target branch
      for(event=locusGen->genChanges.eventBelowRegraft; 
          event->type == OUT_MIG;
          event = event->sons[0]->sons[0])   { ; }
			
      targetBranchId = event->auxId;
      res = executeGenSPR(	locusGen->locusData, 
                            locusGen->genBranches[branchId].likelihoodNodeId,
                            locusGen->genBranches[targetBranchId].likelihoodNodeId,
                            locusGen->genChanges.newCoalEvent->age);
      logAcceptanceRatio = computeLocusDataLikelihood(locusGen->locusData, 1);

      if(logAcceptanceRatio >= 0 || rndu() < exp(logAcceptanceRatio)) {
        // accept genealogy change
        acceptCount++;
		
        event = locusGen->genChanges.newCoalEvent;
        computeNewGenStats(locusGen);
        acceptGenChanges(locusGen);
        // no need to give id for branch above regrafting point
        performRegraft(locusGen, branchId, -1);	
        resetSaved(locusGen->locusData);
        // check if root was switched
        if(res == 1) {
          // branch is regrafted above old root
          locusGen->rootBranchId = event->auxId;
        } else if(res == 2) {
          // branch was just below root and regrafted in another place
          // move up on genealogy and record branch id's until reaching root
          for( ; event->father != NULL; event = event->father) {
            if(event->type == COAL) {
              locusGen->rootBranchId = event->auxId;							
            }
          }
        }
        continue;
      }
    }
		
    // reject genealogy change
		
    resetGenChanges(locusGen);
    revertToSaved(locusGen->locusData);
		
  }// end of for(branchId)
		
  return ((double)acceptCount)/totalCount;
}
/** end of pruneAndReCoalesceSubtrees **/



/***************************************************************************************************************/
/******                              INTERNAL FUNCTION IMPLEMENTATION                                     ******/
/***************************************************************************************************************/



/***********************************************************************************
 *	migrateAndCoalesceLineage
 *	- propose a genalogy modification which re-coalesces a pruned subtree.
 *	- records in genChanges the following things:
 * 		* affected pops and mig bands
 * 		* events whose number of lineages is affected
 * 		* removed migration events
 * 		* removed coalescent event
 *	- returns -1 if too many migration events are sampled along the way
 *	- returns 0, otherwise
 *	- Note: called after pruneGenSubtree() or in order to generate random genealogy
 * 		from scratch
 ***********************************************************************************/
int migrateAndCoalesceLineage(LocusGenealogy* locusGen, int branchId) {
  Event* event;
  double age, rate, migRate, theta, eventSampler;
  int targetNum;
  int i, sourcePop, migBand, numMigBands, *migBandArray;
  int numMigs = locusGen->numMigNodes - locusGen->genChanges.numRemovedEvents/2;
	
  event = locusGen->genBranches[branchId].bottomEvent;
  if(event->type == LEAF)
    event = event->next;
	
  age = event->age;
  theta = locusGen->popTree->pops[ event->popId ]->theta/locusGen->heredityFactor;
  migRate = event->migBands->rate;
  while(event->next != NULL) {
    event = event->next;
    // sample time of next event
    rate = migRate + 2/theta * event->numLineages;
    // rate can be zero if simulating first lineage
    if(rate > 0.0)
      age += rndexp( 1/rate );
   		
    // if top-most event, add slack to allow for coalescence
    if(event->next == NULL) {
      event->age = age + ROOT_SLACK;
    }
   		
    // if event is outside of interval, move on to next interval
    if(rate <= 0.0 || age >= event->age) {
      // no event sampled in this interval
      changeNumLineages(locusGen,event,1);
      age = event->age;
      switch(event->type) {
      case(POP_END):
        // move on to POP_START event and reset rates
        event 	= event->next;
        if(event == NULL)
          break;
        theta 	= locusGen->popTree->pops[ event->popId ]->theta/locusGen->heredityFactor;
        migRate = event->migBands->rate;
        addAffectedPop(locusGen, event->popId);
        addAffectedMigBandSet(locusGen, event->migBands);
        break;
      case(MIG_BAND_END):
        migRate = event->next->migBands->rate;
        break;
      case(MIG_BAND_START):
        addAffectedMigBand(locusGen, event->auxId);
        migRate = event->next->migBands->rate;
        break;
      default:
        break;
      }
      continue;
    }
   		
    // sample type of event (COAL or MIG)
    eventSampler = rate * rndu();
   		
    if(eventSampler >= migRate) {
      // coalescence event -> create event, sample lineage
      // for coalescence and finish up
      event = createEventBefore(locusGen, event, age);
      changeNumLineages(locusGen,event,1);
      locusGen->genChanges.newCoalEvent = event;
      /***********************************
       * randomly sample a lineage which crosses event.
       * eventSampler-migRate is a uniform RV in the
       * interval [0, 2/theta * numLins].
       * So the following simulates a uniform RV in
       * [0, numLins].
       ***********************************/
      targetNum = (int)(eventSampler-migRate) * (theta/2);
      event = findLineageInInterval(locusGen->popTree, event, targetNum);
      if(event == NULL) {
        fprintf(stderr, "\nError: migrateAndCoalesceLineage: locus %d, branch %d\n", locusGen->locusId, branchId);
        fprintf(stderr, "Couldn't find %d lineages crossing event.\n", targetNum);
        printEvent(event);
        printGenealogy(locusGen);
        exit(-1);
      }
      locusGen->genChanges.eventBelowRegraft = event;
      return 0; 	
    }

    // sample migration event
    if(numMigs >= maxMigrations) {
      // migration events exceed quota
      return -1;
    }
		
    // choose migration band according to eventSampler (< migRate)
    numMigBands = event->migBands->numMigBands;
    migBandArray = event->migBands->migBandIds;
    for(i=0; eventSampler > 0 && i < numMigBands; i++) {
      migBand = migBandArray[i];
      eventSampler -= locusGen->popTree->migBands[migBand].migRate;
    }
    if(eventSampler > 0) {
      fprintf(stderr, "\nError: migrateAndCoalesceLineage locus %d, branch %d.\n", locusGen->locusId, branchId);
      fprintf(stderr, "Could not find migration band when sampling migration event (eventSampler=%g, migRate=%g)\n",
             eventSampler, migRate);
      printEvent(event);
      printGenealogy(locusGen);
      exit(-1);
    }
		
    // create IN_MIG event
    numMigs++;
    event = createEventBefore(locusGen, event, age);
    event->auxId = migBand;
    changeNumLineages(locusGen,event,1);
    locusGen->genChanges.addedEvents[ locusGen->genChanges.numAddedEvents++ ] = event;
    //create OUT_MIG event, and move to it (continue up from that point)
    sourcePop = locusGen->popTree->migBands[migBand].sourcePop;
    event = createEventAtAge(locusGen, locusGen->popStartEvents[sourcePop], age);
    locusGen->genChanges.addedEvents[ locusGen->genChanges.numAddedEvents++ ] = event;
    addAffectedPop(locusGen, sourcePop);
    addAffectedMigBandSet(locusGen, event->migBands);
		
  }// end of while(event!= NULL)
	
  // this part should be reachable only if simulating a single lineage
  // in an empty genealogy (when generating a random genealogy).
  if(locusGen->numBranches > 1) {
    fprintf(stderr, "\nError: migrateAndCoalesceLineage: locus %d, branch %d\n", locusGen->locusId, branchId);
    fprintf(stderr, "Unknown reason for leaving while(1) loop.\n");
    printGenealogy(locusGen);
    exit(-1);
  }
}
/** end of migrateAndCoalesceLineage **/



/***********************************************************************************
 *	findLineageInInterval
 *	- finds a lineage which crosses a certain given interval
 * 	- interval is indicated by event, and targetNum is an integer between 0 and
 * 		event->numLineages-1 indicating which lineage to choose
 * 	- returns pointer to last event residing on lineage (branch) below the time
 * 		of the reference event
 * 	- returns NULL if couldn't find enough candidates
 ***********************************************************************************/
Event* findLineageInInterval(PopulationTree* popTree, Event* event, int targetNum) {
  double eventAge = event->age;
  unsigned short* isAncestralToEventPop = popTree->pops[ event->popId ]->isAncestralTo;
	
  // traverse events above event, and mark COAL or IN_MIG events
  // whose children lie below event
  for( ; event != NULL; event = event->next) {
    // this is a discarded event on pruned branch
    if(event->type == DUMMY)	continue;
		
    if(event->sons[0] != NULL && event->sons[0]->age < eventAge && isAncestralToEventPop[ event->sons[0]->popId ]) {
      if(--targetNum < 0)		return event->sons[0];
      if(event->sons[1] != NULL && event->sons[1]->age < eventAge && isAncestralToEventPop[ event->sons[1]->popId ]) {
        if(--targetNum < 0)		return event->sons[1];
      }
    }
		
  }// end of for

  // if didn't find targetNum valid candidates, return NULL
  return NULL;
}
/** end of findLineageInInterval **/



/***********************************************************************************
 *	pruneGenSubtree
 *	- propose a genalogy modification which prunes the subtree below a given branch.
 *	- records in genChanges the following things:
 * 		* affected pops and mig bands
 * 		* events whose number of lineages is affected
 * 		* removed migration events
 * 		* removed coalescent event
 *	- returns 0
 *	- Note: original events are not removed
 ***********************************************************************************/
int pruneGenSubtree(LocusGenealogy* locusGen, int branchId) {
  Event *event, *fatherEvent, *childEvent;	

  event = childEvent = locusGen->genBranches[branchId].bottomEvent;
  fatherEvent = event->father;
  if(event->type == LEAF)
    event = event->next;
	
  addAffectedPop(locusGen, event->popId);
  addAffectedMigBandSet(locusGen,event->migBands);

  while(1) {
    while(event != fatherEvent) {
      event = event->next;
      changeNumLineages(locusGen,event,-1);
      switch(event->type) {
      case(POP_END):
        event = event->next;
        addAffectedPop(locusGen, event->popId);
        addAffectedMigBandSet(locusGen,event->migBands);
        break;
      case(MIG_BAND_START):
        addAffectedMigBand(locusGen,event->auxId);
        break;
      default:
        break;
      }
    }
    // reached father - check if MIG or COAL
    switch(event->type) {
    case(COAL):
      // reached top
      event->type = DUMMY;	// for the purpose of lineage counting later on
      locusGen->genChanges.prunedCoalEvent = event;
      locusGen->genChanges.prunedChildOfEvent = (event->sons[0] != childEvent);
      return 0;
    case(IN_MIG):
      // mark migration events for removal and move on to appropriate OUT_MIG event
      event->type = DUMMY;	// for the purpose of lineage counting later on
      locusGen->genChanges.removedEvents[ locusGen->genChanges.numRemovedEvents++ ] = event;
      event = event->father;
      locusGen->genChanges.removedEvents[ locusGen->genChanges.numRemovedEvents++ ] = event;
      fatherEvent = event->father;
      childEvent = event;
      addAffectedPop(locusGen, event->popId);
      addAffectedMigBandSet(locusGen,event->migBands);
      break;
    default:
      fprintf(stderr, "\nError: pruneGenSubtree: locus %d, branch %d\n", locusGen->locusId, branchId);
      fprintf(stderr, "Illegal type of father event.\n");
      printEvent(event);
      printGenealogy(locusGen);
      exit(-1);
    }
		
  }
	
  fprintf(stderr, "\nError: pruneGenSubtree: locus %d, branch %d\n", locusGen->locusId, branchId);
  fprintf(stderr, "Unknown reason for leaving while(1) loop.\n");
  printGenealogy(locusGen);
  exit(-1);
}
/** end of pruneGenSubtree **/



/***********************************************************************************
 *	proposeEventMove
 *	- propose a genalogy modification which moves an event (coalescent or migration)
 *		to a different age (but within same population).
 *	- records in genChanges the moved event (old and new copies), affected pops and mig bands
 *		as well as all events whose number of lineages is affected
 *	- returns 0
 *	- Note: original event is not removed
 *	- Note: for moving a migration event, you actually need to move two events together
 ***********************************************************************************/
int proposeEventMove(LocusGenealogy* locusGen, Event* event, double newAge) {

  Event *bottomEvent, *topEvent, *newEvent;
  int deltaNumLineages;	
	
  newEvent = createEventAtAge(locusGen, event, newAge);
  newEvent->type = event->type;
  newEvent->father;
	
  locusGen->genChanges.movedEventsOld[ locusGen->genChanges.numMovedEvents ] = event;
  locusGen->genChanges.movedEventsNew[ locusGen->genChanges.numMovedEvents ] = newEvent;
  locusGen->genChanges.numMovedEvents++;
	
  if(newAge > event->age) {
    bottomEvent = event->next;
    topEvent = newEvent;
    deltaNumLineages = (event->type == OUT_MIG) ? (-1) : (1);
  } else {
    bottomEvent = newEvent->next;
    topEvent = event;
    deltaNumLineages = (event->type == OUT_MIG) ? (1) : (-1);
  }
	
  event = bottomEvent;
  addAffectedPop(locusGen, event->popId);
  addAffectedMigBandSet(locusGen, event->migBands);
  while(event!=NULL) {
    changeNumLineages(locusGen, event, deltaNumLineages);
    if(event == topEvent)
      break;
		
    if(event->type == MIG_BAND_START) {
      addAffectedMigBand(locusGen, event->auxId);
    }
    event = event->next;
  }
	
  return 0;	
	
}
/** end of proposeEventMove **/



/***********************************************************************************
 *	computeNewGenStats
 *	- computes genealogy stats for modified genealogy (after proposed change)
 * 	- does not update number of coals or migs (for SPR changes)
 * 	- uses changes tracked in locusGen->genChanges to do this
 * 	- returns 0
 ***********************************************************************************/
int computeNewGenStats(LocusGenealogy* locusGen) {
  Event* event;
  Event** eventArray;
  int ev, numEvents, i, numMigBands, *migBandIds;
  int numLinDiff;
  double deltaStats;
	
  // consider number-of-lineages changes and update coal and mig stats accordingly
  numEvents 	= locusGen->genChanges.numLinChangedEvents;
  eventArray 	= locusGen->genChanges.linChangedEvents;
  for(ev=0; ev<numEvents; ev++) {
    event = eventArray[ev];
    numLinDiff = event->numLineages + event->numLineages_bak;
    if(numLinDiff != 0) {
      // migration stats for interval equal numLins * elapsedTime
      // so difference in stats is elapsedTime * deltaNumLins
      deltaStats = numLinDiff*(event->age - event->prev->age);
      migBandIds = event->migBands->migBandIds;
      numMigBands = event->migBands->numMigBands;
      for(i=0; i<numMigBands; i++) {
        locusGen->genChanges.newStats.migBandStats[ migBandIds[i] ] += deltaStats;				
      }

      // coalescence stats for interval equal numLins * (numLins-1) * elapsedTime
      // so difference in stats is elapsedTime * ( deltaNumLins^2 - deltaNumLins + 2*deltaNumLins*numLins_old ) ==
      // deltaStats * ( deltaNumLins - 1 + 2*deltaNumLins*numLins_old )
      locusGen->genChanges.newStats.popCoalStats[ event->popId ] += deltaStats * (numLinDiff - 1 - 2*event->numLineages_bak);
    }
  }// end of for(ev) - linChangedEvents
	
  return 0;
}
	



/***********************************************************************************
 *	computeDeltaLogLikelihood
 *	- computes delta in log likelihood due to changes in coal and mig stats
 * 	- difference in log-likelihood does not consider changes in number of migration
 * 		and coalescent events (for SPR operations).
 ***********************************************************************************/
double computeDeltaLogLikelihood(LocusGenealogy* locusGen) {
  int i, pop, migBand;
  double rate, deltaLogLikelihood = 0.0;
	
  for(i=0; i<locusGen->genChanges.numAffectedPops; i++) {
    pop = locusGen->genChanges.affectedPops[i];
    rate = 1/(locusGen->popTree->pops[pop]->theta * locusGen->heredityFactor);
    deltaLogLikelihood += rate * (locusGen->genStats.popCoalStats[pop] - locusGen->genChanges.newStats.popCoalStats[pop]);
  }
		
  for(i=0; i<locusGen->genChanges.numAffectedMigBands; i++) {
    migBand = locusGen->genChanges.affectedMigBands[i];
    rate = locusGen->popTree->migBands[migBand].migRate;
    deltaLogLikelihood += rate * (locusGen->genStats.migBandStats[migBand] - locusGen->genChanges.newStats.migBandStats[migBand]);
  }
	
  return deltaLogLikelihood;
}
/** end of computeDeltaLogLikelihood **/



/***********************************************************************************
 *	acceptGenChanges
 *	- "hard wires" changes proposed for genealogy and discards old version
 *	- copies new stats (coal/mig)
 *	- resets numLineages_bak of all events hose numLineages has been changed
 *	- removes all old moved events and "activates" new mved events
 * 	- makes sure to reset all relevant counters so that changes are reset
 *	- does not hard wire specific changes for SPR modifications (see performRegraft)
 * 	- this procedure recycles (old moved) events
 * 	- returns 0
 ***********************************************************************************/
int	acceptGenChanges(LocusGenealogy* locusGen) {
  Event *event, *event1;
  Event **eventArray, **eventArray1;
  int ev, numEvents;
  int i, pop, migBand, numMigBands, *migBandIds;
	
  // copy affected coal/mig stats
  for(i=0; i<locusGen->genChanges.numAffectedPops; i++) {
    pop = locusGen->genChanges.affectedPops[i];
    locusGen->genChanges.isAffectedPop[pop] = 0;
    locusGen->genStats.popCoalStats[pop] = locusGen->genChanges.newStats.popCoalStats[pop];
  }
  locusGen->genChanges.numAffectedPops = 0;
	
  for(i=0; i<locusGen->genChanges.numAffectedMigBands; i++) {
    migBand = locusGen->genChanges.affectedMigBands[i];
    locusGen->genChanges.isAffectedMigBand[migBand] = 0;
    locusGen->genStats.migBandStats[migBand] = locusGen->genChanges.newStats.migBandStats[migBand];
  }
  locusGen->genChanges.numAffectedMigBands = 0;

  // reset numLineages_bak for all events in linChangedEvents[] array
  numEvents 	= locusGen->genChanges.numLinChangedEvents;
  eventArray 	= locusGen->genChanges.linChangedEvents;
  for(ev=0; ev<numEvents; ev++) {
    eventArray[ev]->numLineages_bak = eventArray[ev]->numLineages;
  }
  locusGen->genChanges.numLinChangedEvents = 0;
	
  // replace old moved events with new versions	
  if(0 < (numEvents = locusGen->genChanges.numMovedEvents)) {
    // copy father/children pointers from old to new event
    eventArray 	= locusGen->genChanges.movedEventsOld;
    eventArray1	= locusGen->genChanges.movedEventsNew;
    for(ev=0; ev<numEvents; ev++) {
      event 	= eventArray[i];
      event1	= eventArray1[i];
      event1->type = event->type;
      event1->auxId = event->auxId;
      if(event->type == COAL) {
        locusGen->genBranches[ event->auxId ].bottomEvent = event1;
      }
      i = (event->father->sons[0] != event);
      event->father->sons[i] = event1;
      event1->sons[0] = event->sons[0];
      event1->sons[1] = event->sons[1];
      event->sons[0]->father - event1;
      event->sons[1]->father - event1;
      recycleEvent(locusGen, event);
    }
    locusGen->genChanges.numMovedEvents = 0;		
		
  }
	
  return 0;
}
/** end of acceptGenChanges **/



/***********************************************************************************
 *	resetGenChanges
 *	- resets all recorded changes in genealogy
 *	- removes all proposed added events (including newCoalEvent for SPR)
 * 	- returns 0
 ***********************************************************************************/
int	resetGenChanges(LocusGenealogy* locusGen) {
  Event*	event;
  Event**	eventArray;
  int ev, numEvents;
  int i, pop, migBand, numMigBands, *migBandIds;
	
  // reset affected coal/mig stats
  resetBooleanArray(locusGen->genChanges.isAffectedPop, locusGen->popTree->numPops);
  locusGen->genChanges.numAffectedPops = 0;
  resetBooleanArray(locusGen->genChanges.isAffectedMigBand, locusGen->popTree->numMigBands);
  locusGen->genChanges.numAffectedMigBands = 0;

  // reset numLineages to -numLineages_bak for all events in linChangedEvents[] array
  // THIS IS PRETTY TIME CONSUMING FOR REJECTION. CONSIDER CHANGING !!!
  numEvents 	= locusGen->genChanges.numLinChangedEvents;
  eventArray 	= locusGen->genChanges.linChangedEvents;
  for(ev=0; ev<numEvents; ev++) {
    eventArray[ev]->numLineages = eventArray[ev]->numLineages_bak = -eventArray[ev]->numLineages_bak;
  }
  locusGen->genChanges.numLinChangedEvents = 0;
	
  // remove new moved events	
  if(0 < (numEvents = locusGen->genChanges.numMovedEvents)) {
    eventArray	= locusGen->genChanges.movedEventsNew;
    for(ev=0; ev<numEvents; ev++) {
      recycleEvent(locusGen, eventArray[ev]);
    }
    locusGen->genChanges.numMovedEvents = 0;		
  }
	
  // remove added events	
  if(0 < (numEvents = locusGen->genChanges.numAddedEvents)) {
    eventArray	= locusGen->genChanges.addedEvents;
    for(ev=0; ev<numEvents; ev++) {
      recycleEvent(locusGen, eventArray[ev]);
    }
    locusGen->genChanges.numAddedEvents = 0;		
  }
	
  // reassign types to removed migration nodes and pruned coal nodes	
  if(0 < (numEvents = locusGen->genChanges.numRemovedEvents)) {
    // copy father/children pointers from old to new event
    eventArray	= locusGen->genChanges.removedEvents;
    for(ev=0; ev<numEvents; ev++) {
      if(eventArray[ev]->type == DUMMY)	eventArray[ev]->type = IN_MIG;
    }
    locusGen->genChanges.numRemovedEvents = 0;		
  }
	
  if(NULL != locusGen->genChanges.prunedCoalEvent) {
    locusGen->genChanges.prunedCoalEvent->type = COAL;
  }	
	
  // remove new SPR COAL event
  if(NULL != locusGen->genChanges.newCoalEvent) {
    recycleEvent(locusGen,locusGen->genChanges.newCoalEvent);
  }
	
  return 0;
}
/** end of resetGenChanges **/




/***********************************************************************************
 *	performRegraft
 *	- "hard wires" changes proposed by an SPR proposal or generating a genealogy from scratch
 *	- branchId is the id of the branch being regrafted
 *	- newBranchId is the id of the branch above regrafted branch
 *		(only given when generating a new genealogy)
 * 	- recycles all old migration events
 *	- chains new migration events
 *	- adjusts locusGen->genStats.migBandNumMigs[] array accordingly
 * 	- performs the (prune and) regrafting itself
 *	- adjusts locusGen->genStats.popNumCoals[] array, if necessary
 * 	- returns 0
 ***********************************************************************************/
int	performRegraft(LocusGenealogy* locusGen, int branchId, int newBranchId) {
	
  int 	ev, numEvents;
  Event**	eventArray;
  Event*	event;

  // remove migration events marked for removal and record change in migBandNumMigs[]
  if(0 < (numEvents = locusGen->genChanges.numRemovedEvents)) {
    eventArray 	= locusGen->genChanges.removedEvents;
    for(ev=0; ev<numEvents; ev++) {
      event = eventArray[ev];
      if(event->type == IN_MIG) {
        locusGen->genStats.migBandNumMigs[ event->auxId ]--;					
      }
      recycleEvent(locusGen, event);
    }
    locusGen->genChanges.numRemovedEvents = 0;		
  }
	
  // add newly added migration events and record change in migBandNumMigs[]
  if(0 < (numEvents = locusGen->genChanges.numAddedEvents)) {
    eventArray 	= locusGen->genChanges.addedEvents;
    eventArray[0]->sons[0] = locusGen->genBranches[branchId].bottomEvent;
    locusGen->genBranches[branchId].bottomEvent->father = eventArray[0];
    for(ev=0; ev<numEvents-1; ev++) {
      eventArray[ev]->father = eventArray[ev+1];
      eventArray[ev+1]->sons[0] = eventArray[ev];
      if(ev%2 == 0) {
        event->type = IN_MIG;
        locusGen->genStats.migBandNumMigs[ event->auxId ]++;					
      } else {
        event->type = OUT_MIG;					
      }
    }
    locusGen->genChanges.newCoalEvent->sons[0] = eventArray[numEvents-1];
    eventArray[numEvents-1]->father = locusGen->genChanges.newCoalEvent;		
    locusGen->genChanges.numAddedEvents = 0;
  } else{
    // attach new coal event to bottom event of branch
    locusGen->genChanges.newCoalEvent->sons[0] = locusGen->genBranches[branchId].bottomEvent;
    locusGen->genBranches[branchId].bottomEvent->father = locusGen->genChanges.newCoalEvent;		
  }
	
  // performs actual topology change
  if(newBranchId < 0) {
    // bypass old (pruned) event and remove it
    event = locusGen->genChanges.prunedCoalEvent;
    ev = (event->father->sons[0] != event);
    event->father->sons[ev] = event->sons[ 1-locusGen->genChanges.prunedChildOfEvent ];
    event->father->sons[ev]->father = event->father;
    // set id of new branch as old one
    branchId = event->auxId;
    recycleEvent(locusGen, event);
  }

  event = locusGen->genChanges.newCoalEvent;
  event->type 	= COAL;
  event->auxId 	= newBranchId;
  locusGen->genBranches[newBranchId].bottomEvent = event;
  event->sons[1]	= locusGen->genChanges.eventBelowRegraft;
  event->father	= locusGen->genChanges.eventBelowRegraft->father;
  event->sons[1]->father = event;
  ev = (event->father->sons[0] != event->sons[1]);
  event->father->sons[ev] = event;
  locusGen->genChanges.newCoalEvent = NULL;
	
  // update number of coalescent events, if changed
  if(locusGen->genChanges.prunedCoalEvent == NULL) {
    // this is for when genealogy is generated
    locusGen->genStats.popNumCoals[ locusGen->genChanges.newCoalEvent->popId ]++;
  } else if(locusGen->genChanges.prunedCoalEvent->popId != locusGen->genChanges.newCoalEvent->popId)	{
    locusGen->genStats.popNumCoals[ locusGen->genChanges.prunedCoalEvent->popId ]--;
    locusGen->genStats.popNumCoals[ locusGen->genChanges.newCoalEvent->popId ]++;
  }
	
	
  return 0;
}
/** end of performRegraft **/



/***********************************************************************************
 *	addAffectedPop
 *	- adds a population id to the list of pops affected by change
 * 	- checks to see if not already marked as affected
 * 	- returns 0
 ***********************************************************************************/
int addAffectedPop(LocusGenealogy* locusGen, int popId){
	
  if(!locusGen->genChanges.isAffectedPop[popId]) {
    locusGen->genChanges.isAffectedPop[popId] = 1;
    locusGen->genChanges.affectedPops[ locusGen->genChanges.numAffectedPops++ ] = popId;
    locusGen->genChanges.newStats.popCoalStats[popId] = locusGen->genStats.popCoalStats[popId];
    locusGen->genChanges.newStats.popNumCoals[popId]  = locusGen->genStats.popNumCoals[popId];
  }
	
  return 0;
}
/** end of addAffectedPop **/



/***********************************************************************************
 *	addAffectedMigBand
 *	- adds a migration band id to the list of mig-bands affected by change
 * 	- checks to see if not already marked as affected
 * 	- returns 0
 ***********************************************************************************/
int addAffectedMigBand(LocusGenealogy* locusGen, int migBand){
	
  if(!locusGen->genChanges.isAffectedMigBand[migBand]) {
    locusGen->genChanges.isAffectedMigBand[migBand] = 1;
    locusGen->genChanges.affectedMigBands[ locusGen->genChanges.numAffectedMigBands++ ] = migBand;
    locusGen->genChanges.newStats.migBandStats[migBand] = locusGen->genStats.migBandStats[migBand];
    locusGen->genChanges.newStats.migBandNumMigs[migBand]  = locusGen->genStats.migBandNumMigs[migBand];
  }
	
  return 0;
}
/** end of addAffectedMigBand **/



/***********************************************************************************
 *	addAffectedMigBandSet
 *	- adds a set of migration band id's to the list of mig-bands affected by change
 * 	- calls addAffectedMigBand
 * 	- returns 0
 ***********************************************************************************/
int addAffectedMigBandSet(LocusGenealogy* locusGen, MigrationBandSet* migBandSet){
  int i;
	
  for(i=0; i<migBandSet->numMigBands; i++) {
    addAffectedMigBand(locusGen, migBandSet->migBandIds[i]);
  }
	
  return 0;
}
/** end of addAffectedMigBandSet **/



/***********************************************************************************
 *	changeNumLineages
 *	- changes number of lineages which cross interval indicated by event
 *	- if number of lineages has not been changed before, records old number in 
 * 		event->numLineages_bak (in negative form).
 *	- returns 0
 ***********************************************************************************/
int changeNumLineages(LocusGenealogy* locusGen, Event* event, int deltaNumLins){
  if(event->numLineages_bak > 0) {
    event->numLineages_bak = -event->numLineages_bak;
  }
  event->numLineages += deltaNumLins;
}
/** end of changeNumLineages **/



/***********************************************************************************
 *	initEvents
 *	- initializes events for genealogy
 *	- allocates memory for all events, sets up event stack, and creates start and end
 *		events for every population.
 *	- returns 0
 *	NEED TO ADD MIGRATION BAND SHIT !!
 ***********************************************************************************/
int initEvents(LocusGenealogy* locusGen, int numEvents){
  int i, pop;
  Event* event;

  locusGen->eventArray = event = (Event*)malloc(numEvents*sizeof(Event));
	
  if(event == NULL) {
    fprintf(stderr, "\nError: Out Of Memory Event Array for locus %d.\n",locusGen->locusId);
    exit(-1);
  }
	
  // set up POP_START events
  for(pop=0; pop<locusGen->popTree->numPops; pop++) {
    locusGen->popStartEvents[pop] = event;
    event->type = POP_START;
    event->popId = pop;
    event->age   = locusGen->popTree->pops[pop]->age;
    event->numLineages = event->numLineages_bak = 0;
    event->auxId = -1;
    event->migBands = NULL;
    event->sons[0] = event->sons[1] = event->father = NULL;
    event->prev = NULL;
    event++;		
  }
  // set up POP_END events
  for(pop=0; pop<locusGen->popTree->numPops; pop++) {
    event->type = POP_END;
    event->popId = pop;
    event->age   = locusGen->popTree->pops[pop]->father->age;
    event->numLineages = event->numLineages_bak = 0;
    event->auxId = -1;
    event->migBands = NULL;
    event->sons[0] = event->sons[1] = event->father = NULL;
    event->prev = locusGen->popStartEvents[pop];
    event->prev->next = event;
    event->next = locusGen->popStartEvents[ locusGen->popTree->pops[pop]->father->id ];
    event++;		
  }

  // put remaining events in stack
  locusGen->eventStackTop = event;
  for(i=2*pop<locusGen->popTree->numPops; i<numEvents; i++) {
    event->type = DUMMY;
    event->popId = -1;
    event->age   = -1;
    event->numLineages = event->numLineages_bak = 0;
    event->auxId = -1;
    event->migBands = NULL;
    event->sons[0] = event->sons[1] = event->father = NULL;
    event->prev = NULL;
    if(i<numEvents-1) {
      event->next = (++event);
    } else {
      event->next = NULL;
    }
  }
	
  return 0;
	
}
/** end of initEvents **/



/***********************************************************************************
 *	recycleEvent
 *	- removes event from chain and moves it to top of stack
 *	- returns 0
 ***********************************************************************************/
int recycleEvent(LocusGenealogy* locusGen, Event* event){

  event->next->prev = event->prev;
  event->prev->next = event->next;

  event->next = locusGen->eventStackTop;
  locusGen->eventStackTop = event;

  return 0;
}
/** end of recycleEvent **/



/***********************************************************************************
 *	createEventBefore
 *	- creates new event immediately before specified event
 *	- does not set type, father, sons, auxId
 *	- takes new event from top of stack
 *	- returns pointer to new event
 ***********************************************************************************/
Event* createEventBefore(LocusGenealogy* locusGen, Event* event, double age) {

  Event* newEvent = locusGen->eventStackTop;
  if(newEvent == NULL) {
    fprintf(stderr, "\nError: Empty event pool in locus %d.\n", locusGen->locusId);
    printGenealogy(locusGen);
    exit(-1);
  }
  locusGen->eventStackTop = newEvent->next;
	
  newEvent->next = event;
  newEvent->prev = event->prev;
  event->prev = newEvent;
  event->prev->next = event;
  newEvent->type = DUMMY;
  newEvent->numLineages = newEvent->numLineages_bak = event->numLineages;
  newEvent->popId = event->popId;
  newEvent->migBands = event->migBands;
  newEvent->age = age;

  return newEvent;
}
/** end of createEventBefore **/



/***********************************************************************************
 *	createEventAtAge
 *	- creates new event at given age in given pop
 * 	- starts looking from refEvent.
 *	- does not set type, father, sons, auxId
 *	- calls createEventBefore()
 *	- returns pointer to new event
 ***********************************************************************************/
Event* createEventAtAge(LocusGenealogy* locusGen, Event* refEvent, double age) {

  Event* newEvent;
	
  if(refEvent->age <= age) {
    // search up
    for( ; refEvent->type != POP_END && refEvent->age <= age; refEvent = refEvent->next) {
      // do nothing
    }
    // if age of POP_END is below age in root pop, extend POP_END age
    if(refEvent->age <= age) {
      if(refEvent->popId == locusGen->popTree->rootPop) {
        refEvent->age = age + ROOT_SLACK;
      } else {
        fprintf(stderr, "\nError: createEventAtAge: trying to create event in locus %d, pop %d at age %g (age of father pop is %g).\n",
               locusGen->locusId, refEvent->popId, age, refEvent->age);
        printEvent(refEvent);
        printGenealogy(locusGen);
        exit(-1);
      }
    }
  } else {
    // search down
    for( ; refEvent->type != POP_START && refEvent->prev->age > age; refEvent = refEvent->prev) {
      // do nothing
    }
    // if age of POP_END is below age in root pop, extend POP_END age
    if(refEvent->type == POP_START) {
      fprintf(stderr, "\nError: createEventAtAge: trying to create event in locus %d, pop %d at age %g (age of pop is %g).\n",
             locusGen->locusId, refEvent->popId, age, refEvent->age);
      printEvent(refEvent);
      printGenealogy(locusGen);
      exit(-1);
    }
  }
	
	
  newEvent = createEventBefore(locusGen, refEvent, age);

  return newEvent;
}
/** end of createEventAtAge **/



/***********************************************************************************
 *	printEvent
 *	- prints event
 ***********************************************************************************/
void printEvent(Event* event) {

  return;
}
/** end of printEvent **/



/***********************************************************************************
 *	printEvent
 *	- prints entire locus genalogy
 ***********************************************************************************/
void printGenealogy(LocusGenealogy* locusGen) {

  return;
}
/** end of printGenealogy **/

	


/***************************************************************************************************************/
/******                                        END OF FILE                                                ******/
/***************************************************************************************************************/
