#ifndef LOCUS_DATA_LIKELIHOOD_H
#define LOCUS_DATA_LIKELIHOOD_H

/** 
	\file LocusGenealogy.h
    DataStrucutres for locus-specific genealogy.
	
	Contains the relevant data structures and function declarations
	for a locus-specific genealogy.
	
*/

#include "PopulationTree.h"


/***************************************************************************************************************/
/******                                              DATA TYPES                                           ******/
/***************************************************************************************************************/



/***********************************************************************************
*	LocusGenealogy
*	- Data type which holds genealogy for a given locus
*	- implemented in LocusGenealogy.c
***********************************************************************************/
typedef struct LOCUS_GENEALOGY LocusGenealogy;



/***************************************************************************************************************/
/******                               EXTERNAL FUNCTION DECLARATIONS                                      ******/
/***************************************************************************************************************/



/***********************************************************************************************/
/******                          INITIALIZATION FUNCTIONS                                 ******/
/***********************************************************************************************/



/***********************************************************************************
*	initLocusGenealogy
*	- allocates all memory resources for locus (other than ones for data)
* 	- returns 0
***********************************************************************************/
int		initLocusGenealogy(LocusGenealogy* locusGen, PopulationTree* popTree, int id, int numSamples);



/***********************************************************************************
*	initLocusData
*	- allocates all memory resources for locus data and initializes them
* 	- receives pattern info to copy into locus data
* 	- returns 0
***********************************************************************************/
int		initLocusData(LocusGenealogy* locusGen, int numSeqPatterns, char** patternArray, int* patternCounts);



/***********************************************************************************
*	freeLocusMemeory
*	- frees all memory resources and data structures for locus
* 	- returns 0
***********************************************************************************/
int		freeLocusMemeory(LocusGenealogy* locusGen);



/***********************************************************************************
*	setMaxNumMigrations
*	- sets global maximum for number of migrations in all loci
* 	- should be called once before initializing any of the loci, and only then.
***********************************************************************************/
void	setMaxNumMigrations(int maxNumMigrations);



/***********************************************************************************************/
/******                           GENALOGY MODIFICATIONS                                  ******/
/***********************************************************************************************/



/***********************************************************************************
*	addLeaf
*	- adds a leaf node in genealogy in indicated population
* 	- used when initializing a genealogy
* 	- initializes branch for this leaf
* 	- returns 0
***********************************************************************************/
int		addLeaf(LocusGenealogy* locusGen, int leafId, int popId);



/***********************************************************************************
*	generateRandomGenealogy
*	- generates a random genealogy according to population tree parameters
* 	- simulated lineages one-by-one from leaves up in the population tree until they 
* 		coalesce with some other lineage.
* 	- calls migrateAndCoalesceLineage in each iteration
* 	- returns 0
***********************************************************************************/
int		generateRandomGenealogy(LocusGenealogy* locusGen);



/***********************************************************************************
*	changeLocusMutationRate
*	- changes the locus-specific mutation rate to a new value.
***********************************************************************************/
void	changeLocusMutationRate(LocusGenealogy* locusGen, double newRate);



/***********************************************************************************
*	perturbNodeAges
*	- perturb ages of all nodes (COAL, MIG) without modifying structure of genealogy.
*	- input finetune argument defines the size of interval in which age is perturbed
*	- returns the proportion of successful proposals.
*	Note: we consider also migration events along root branch.
***********************************************************************************/
double	perturbNodeAges(LocusGenealogy* locusGen, double finetune);

	
		
/***********************************************************************************
*	pruneAndReCoalesceSubtrees
*	- traverses all branches (excluding one above root), prunes the subtree below edge
*		and re-coalesces the lineage from the root of the pruned subtree back into 
*		remaining genealogy, according to model parameters (pop sizes and mig rates).
*	- returns the proportion of successful proposals.
***********************************************************************************/
double	pruneAndReCoalesceSubtrees(LocusGenealogy* locusGen);



/***************************************************************************************************************/
/******                                        END OF FILE                                                ******/
/***************************************************************************************************************/
#endif
