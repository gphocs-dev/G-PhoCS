#ifndef POPULATION_TREE_H
#define POPULATION_TREE_H

/** 
	\file PopulationTree.h 
    Data Structures to implement handling a population tree with migration bands
	
	Contains the relevant data structures and procedure implementations
	for handling a population tree with migration bands.	
*/


#include "utils.h"
/***************************************************************************************************************/
/******                                      EXTERNAL CONSTANTS                                           ******/
/***************************************************************************************************************/



/***************************************************************************************************************/
/******                                              DATA TYPES                                           ******/
/***************************************************************************************************************/



/***********************************************************************************
*	GammaPrior
*	- Holds parameters for gamma prior
***********************************************************************************/
typedef struct GAMMA_PRIOR {
	double	sampleStart;	// mean for sample start (sampling 10% below and above)
	double	alpha;		// alpha for gamma prior
	double	beta;		// beta  for gamma prior
} GammaPrior;



/***********************************************************************************
*	MigrationBandSet
*	- Holds relevant info for set of migration band active at a certain time interval
* 		at their target population
***********************************************************************************/
typedef struct MIGRATION_BAND_SET MigrationBandSet;

struct MIGRATION_BAND_SET {
	int		numMigBands;		// length of migBandIds[] array
	int*	migBandIds;			// an array of migration band id's
	double	rate;
	double	age;				// start time of band set
	MigrationBandSet* next;		// next migration band set (back in time)
	MigrationBandSet* prev;		// previous migration band set (next in time)
};



/***********************************************************************************
*	MigrationBand
*	- Holds relevant info for migration band
***********************************************************************************/
typedef struct MIGRATION_BAND {
	int		sourcePop;			// id of source population
	int		targetPop;			// id of target population
	double	migRate;			// migration rate for band
//	double	upperBound;			// upper bound for uniform prior
	GammaPrior	migRatePrior;	// parameters for gamma-prior of migration rate - NOT IN USE !!
	double	startTime;			// start time for migration band
	double	endTime;			// end time for migration band
	MigrationBandSet* firstSet;	// set of migration bands right after this one starts
	MigrationBandSet* lastSet;	// set of migration bands right before this one ends
} MigrationBand;



/***********************************************************************************
*	Population
*	- Holds relevant info for population
***********************************************************************************/
typedef struct POPULATION Population;

struct POPULATION {
	char	name[STRING_LENGTH];		// name of population
	int 	id;					// population identifier
	int     numSamples;         // number of samples for current population (leaf in tree).
	double	age;				// start time for population (for ancestral population)
	double	sampleAge;			// age of samples in current population - used for extinct populations
	unsigned short updateSampleAge;			// set to 1 iff algorithm is to update the sample age
	double	theta;				// theta parameter for population
	GammaPrior	thetaPrior;		// parameters for gamma-prior of theta
	GammaPrior	agePrior;		// parameters for gamma-prior of age (for ancestral population)
	Population*	father;			// father population
	Population*	sons[2];		// two child populations (for ancestral population)
	unsigned short*	isAncestralTo;	// a boolean array indicating all descendant populations

	// migration bands
	int 	numInMigBands;		// length of inMigBands[] array
	int* 	inMigBands;			// array of in migration band id's
	int 	numOutMigBands;		// length of outMigBands[] array
	int* 	outMigBands;		// array of out migration band id's
	MigrationBandSet* migBandSequence;	// pointer to first migration band set in the sequence
										// of (incoming) sets active along the population
};



/***********************************************************************************
*	PopulationTree
*	- Holds relevant info for population
***********************************************************************************/
typedef struct POPULATION_TREE {
	int 	numCurPops;			// number of current populations
	int 	numPops;			// number of populations in tree ( = 2*numCurPops-1 )
	int		numMigBands;		// number of migration bands in tree
	int		rootPop;			// id of root population
	Population**	pops;		// an array of pointers to populations
	MigrationBand*	migBands;	// an array of migration bands
	MigrationBandSet* 	migBandSetStackTop;	// pointer to top of MigrationBandSet stack

	Population*			popArray;			// pointer to allocated memory for all populations
	MigrationBandSet* 	migBandSetArray;	// pointer to allocated memory for MigrationBandSets
	unsigned short*		isAncestralArray;	// pointer to allocated memory for isAncestralTo[] arrays
	int*				migBandIdArray;		// pointer to allocated memory for in/out migband arrays for pops
} PopulationTree;



/***************************************************************************************************************/
/******                               EXTERNAL FUNCTION DECLARATIONS                                      ******/
/***************************************************************************************************************/



/***********************************************************************************
*	createPopTree
*	- allocates basic memory for population tree (no migration bands yet)
* 	- returns pointer to newly allocated population tree
***********************************************************************************/
PopulationTree*	createPopTree(int numPops);



/***********************************************************************************
*	initMigrationBands
*	- initializes data structures for migration bands in population tree (including allocating some memory)
* 	- sets start and end times
* 	- for each population creates a timed-sequence of migration band sets
* 	- returns 0
***********************************************************************************/
int initMigrationBands(PopulationTree* popTree);



/***********************************************************************************
*	freePopTree
*	- frees all memory allocated for population tree
* 	- returns 0
***********************************************************************************/
int	freePopTree(PopulationTree* popTree);


	
/***********************************************************************************
*	printPopulationTree
*	- prints population tree
***********************************************************************************/
void printPopulationTree(PopulationTree* popTree, FILE* stream, int printTauTheta);



/***********************************************************************************
 *	getPopIdByName
 * 	- returns a population id of a population given its name (-1 if no match is found)
 * 	- used primarily to decode migration bands as specified in control file 
 * 		(called by readControlFile).
 ***********************************************************************************/
int getPopIdByName(PopulationTree* popTree, const char* name);



/***********************************************************************************
*	samplePopParameters
*	- samples population parameters according to prior average (only thetas and taus)
* 	- each parameter is sampled uniformly in the interval [0.9,1.1]*mean
* 		(where mean is the prior mean for that parameter)
* 	- makes sure a population's age does not exceed its father's
*	- initializes all migration rates to 0.
* 	- returns 0
***********************************************************************************/
int samplePopParameters(PopulationTree* popTree);



/***********************************************************************************
*	sampleMigRates
*	- samples migration rates for all mig bands
* 	- each rate is sampled uniformly in the interval [0.9,1.1]*mean
* 		(where mean is the prior mean for that parameter)
* 	- returns 0
***********************************************************************************/
int sampleMigRates(PopulationTree* popTree);



/***********************************************************************************
*	updateMigrationBandTimes
*	- updates start and end times of given migration band according to ages of populations
*	- returns 0 if no change was made, and 1 otherwise
***********************************************************************************/
unsigned short updateMigrationBandTimes(PopulationTree* popTree, int migBand);



/***********************************************************************************
*	computeMigrationBandTimes
*	- traverses all migration bands and sets their start and end times according to 
*		times of target and source populations.
* 	- returns the number of migration bands with zero span
***********************************************************************************/
int computeMigrationBandTimes(PopulationTree* popTree);

	
	
/***************************************************************************************************************/
/******                                        END OF FILE                                                ******/
/***************************************************************************************************************/
#endif
