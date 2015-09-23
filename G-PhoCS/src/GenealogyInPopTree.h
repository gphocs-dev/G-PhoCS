#ifndef GENEALOGY_IN_POP_TREE_H
#define GENEALOGY_IN_POP_TREE_H

/** 
	\file GenealogyInPopTree.h
    Data Structures dealing with information about a genealogy embeded in a population tree.
	
	Contains the relevant declarations for data structures and procedures dealing with
	holding detailed information about a genealogy embedded in a population tree.
	This information is stored in form of a chain of intervals.
	Each interval resides in a specific population.
	Transition between one interval to the next is defined either by a change
	in number of lineages or by a change in migration bands.
	
	NO MIGRATION YET
	NO CHANGE TO STRUCTURE OF TREE YET
	
*/




/***************************************************************************************************************/
/******                                              DATA TYPES                                           ******/
/***************************************************************************************************************/



/**	LocusGenealogy
    Data type which holds genealogy for a given locus
    @note Implemented in LocusGenealogy.c
*/
typedef struct LOCUS_GENEALOGY LocusGenealogy;



/***************************************************************************************************************/
/******                               EXTERNAL FUNCTION DECLARATIONS                                      ******/
/***************************************************************************************************************/



/***********************************************************************************************/
/******                          INITIALIZATION FUNCTIONS                                 ******/
/***********************************************************************************************/



/**	initLocusGenealogy
    Allocates all memory resources for locus (other than ones for data)
    @pre setMaxNumMigrations called 
    @param locusGen Locus Genealogy to allocate memory for 
    @param popTree Population tree populated from control file so we know how much memory to allocate
    @param id Number used to identify this new locus
    @param numSamples Number of samples for this locus
    @return 0
*/
int		initLocusGenealogy(LocusGenealogy* locusGen, PopulationTree* popTree, int id, int numSamples);



/**	initLocusData
    Allocates all memory resources for locus data and initializes them
    @param locusGen Locus Genealogy that has already had its genealogy initialized with initLocusGenealogy()
    @param numSeqPatterns Number of patterns in patternArray
    @param patternArray Array of sequence patterns
    @param patternCounts Number of occurances for each pattern
    @return 0
*/
int		initLocusData(LocusGenealogy* locusGen, int numSeqPatterns, char** patternArray, int* patternCounts);



/** freeLocusMemeory
    Frees all memory resources and data structures for locus
    @return 0
*/
int		freeLocusMemeory(LocusGenealogy* locusGen);



/**	setMaxNumMigrations
    Sets Global Maximum for number of migrations in all loci
    @param maxNumMigrations Global maximum for number of migrations in all loci
    @note Should be called once before initializing any of the loci, and only then.
*/
void	setMaxNumMigrations(int maxNumMigrations);



/***********************************************************************************************/
/******                           GENALOGY MODIFICATIONS                                  ******/
/***********************************************************************************************/



/**	addLeaf
    Adds a leaf node in genealogy in indicated population
    @param locusGen Locus Genealogy that has been initialized
    @param leafId Number to identify this new leaf by
    @param popId Population ID for event after this one
 	@note Used when initializing a genealogy
    @note Initializes branch for this leaf
    @return 0
*/
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



/**	changeLocusMutationRate
    Changes the locus-specific mutation rate to a new value.
    @param locusGen Locus to have its mutation rate modified
    @param newRate New mutation rate for locus
*/
void	changeLocusMutationRate(LocusGenealogy* locusGen, double newRate);



/**	perturbNodeAges
    Perturb ages of all nodes (COAL, MIG) without modifying structure of genealogy.
    @param locusGen Locus to containing nodes whos ages will be perturbed
    @param finetune Finetune parameter used to adjust how much the ages are perturbed
    @note we consider also migration events along root branch.
    @return The proportion of successful proposals.
*/
double	perturbNodeAges(LocusGenealogy* locusGen, double finetune);

	
		
/**	pruneAndReCoalesceSubtrees
    Prunes subtree and re-coalesces.
    Traverses all branches (excluding one above root), prunes the subtree below edge
    and re-coalesces the lineage from the root of the pruned subtree back into 
	remaining genealogy, according to model parameters (pop sizes and mig rates).
    @return The proportion of successful proposals.
*/
double	pruneAndReCoalesceSubtrees(LocusGenealogy* locusGen);



/***************************************************************************************************************/
/******                                        END OF FILE                                                ******/
/***************************************************************************************************************/
#endif
