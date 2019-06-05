#ifndef POPULATION_TREE_H
#define POPULATION_TREE_H

/** 
	\file PopulationTree.h 
    Data Structures to implement handling a population tree with migration bands
	
	Contains the relevant data structures and procedure implementations
	for handling a population tree with migration bands.	
*/


#include "utils.h"
#include <map>
#include <vector>
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
*	MigrationBand
*	- Holds relevant info for migration band
***********************************************************************************/
typedef struct MIGRATION_BAND {
    int     id;                 //id of migration band
	int		sourcePop;			// id of source population
	int		targetPop;			// id of target population
	double	migRate;			// migration rate for band
//	double	upperBound;			// upper bound for uniform prior
	GammaPrior	migRatePrior;	// parameters for gamma-prior of migration rate - NOT IN USE !!
	double	startTime;			// start time for migration band
	double	endTime;			// end time for migration band

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
	#define	LEFT	0 			// index of left son population in sons array
	#define	RIGHT	1 			// index of right son population in sons array
	unsigned short*	isAncestralTo;	// a boolean array indicating all descendant populations

	// migration bands
	int 	numInMigBands;		// length of inMigBands[] array
	int* 	inMigBands;			// array of in migration band id's
	int 	numOutMigBands;		// length of outMigBands[] array
	int* 	outMigBands;		// array of out migration band id's

};


/***********************************************************************************
*	PopulationTree
*	- Holds relevant info for population
***********************************************************************************/
typedef struct TIME_MIG_BANDS TimeMigBands;
typedef struct MIG_BANDS_PER_TARGET_POP MigBandsPerTarget;

typedef struct _POPULATION_TREE
{
	int 	numCurPops;			// number of current populations
	int 	numPops;			// number of populations in tree ( = 2*numCurPops-1 )
	int		numMigBands;		// number of migration bands in tree
	int		rootPop;			// id of root population
	Population**	pops;		// an array of pointers to populations
	MigrationBand*	migBands;	// an array of migration bands

	Population*			popArray;			// pointer to allocated memory for all populations
	unsigned short*		isAncestralArray;	// pointer to allocated memory for isAncestralTo[] arrays
	int*				migBandIdArray;		// pointer to allocated memory for in/out migband arrays for pops

    // migBandsPerTarget is a vector of size num-pops,
    //contains MigBandsPerTarget structs.
    std::vector<MigBandsPerTarget> migBandsPerTarget;

} PopulationTree;


/***********************************************************************************
*	TIME_MIG_BANDS
*	Time Mig band is defined by a start time and end time, and it contains
 *	a vector of migration bands active within that period
 *	(migration bands share same target pop)
***********************************************************************************/
typedef struct TIME_MIG_BANDS
{
    double startTime;  //start time of band
    double endTime;    //end time of band
    std::vector<MigrationBand*> migBands; //pointers to mig bands

    TIME_MIG_BANDS() : startTime(-1), endTime(-1) {}; //constructor
    TIME_MIG_BANDS(double s, double e) : startTime(s), endTime(e) {};//constructor

} TimeMigBands;


/***********************************************************************************
*	MIG_BANDS_PER_TARGET_POP
*	Struct contains the two following vectors which are associated with each pop:
    1. Vector of mig bands which pop is their target pop.
    2. Vector of time bands, where each time band holds a vector of mig bands
        active in that time band.
***********************************************************************************/
typedef struct MIG_BANDS_PER_TARGET_POP
{
    std::vector<MigrationBand*> migBands;
    std::vector<TimeMigBands> timeMigBands;
} MigBandsPerTarget;



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
*	- if change was made, re-construct live mig-bands data structure
***********************************************************************************/
unsigned short updateMigrationBandTimes(PopulationTree* popTree, int migBand);


/***********************************************************************************
*	computeMigrationBandTimes
*	- traverses all migration bands and sets their start and end times according to 
*		times of target and source populations.
* 	- returns the number of migration bands with zero span
***********************************************************************************/
int computeMigrationBandTimes(PopulationTree* popTree);


/***********************************************************************************
*	getMigBandByPops
*	- returns pointer to a mig band with the given source and target populations
***********************************************************************************/
MigrationBand *
getMigBandByPops(PopulationTree* popTree, int sourcePop, int targetPop);


/***********************************************************************************
*	getMigBandById
*	- returns pointer to a mig band with the given ID
***********************************************************************************/
MigrationBand * getMigBandByID(PopulationTree* popTree, int id);


/*******************************************************************************
 *	initializeLivingMigBands
 *	create N elements in livingMigBands vector, where N is num of pops
 *	divide migration bands into groups with same target pop
 ******************************************************************************/
void initializeMigBandTimes(PopulationTree* popTree);


/*******************************************************************************
 *	constructMigBandsTimes
 *
 *
 ******************************************************************************/
void constructMigBandsTimes(PopulationTree* popTree);


/*******************************************************************************
 *	getLiveMigBands
 *	for the given target pop, returns a time band containing the given age,
 *	and null if not found such.
    @param: popTree, target pop, age
    @return: pointer to a time band struct
 ******************************************************************************/
TimeMigBands *
getLiveMigBands(PopulationTree* popTree, int target_pop, double age);



/***************************************************************************************************************/
/******                                        END OF FILE                                                ******/
/***************************************************************************************************************/
#endif
