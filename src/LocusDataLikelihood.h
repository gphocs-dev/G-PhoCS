#ifndef LOCUS_DATA_LIKELIHOOD_H
#define LOCUS_DATA_LIKELIHOOD_H
/** 
	\file LocusDataLikelihood.h
    Compute likelihood of data given locus-specific genealogy.

	This header file describes the interface functions for computing
	likelihood of data given locus-specific genealogy.
	When changes are made to genealogy, old version is being saved, in case changes 
	are rejected.
	Procedures for reading sequence data and doing initial analysis are taken from
	PAML and original MCMCcoal (NEED TO CHANGE THIS !!)
	
*/



#include "GenericTree.h"

#define OPT1
#define OPT2_not

/***************************************************************************************************************/
/******                                              DATA TYPES                                           ******/
/***************************************************************************************************************/



/***********************************************************************************
 *	LikelihoodNode
 *	- Data type which holds data for computing conditional probability
 *		of all base assignments to a certain genealogy (coalescent) node.
 ***********************************************************************************/
typedef struct LIKELIHOOD_NODE {
    int father;              // father of node in genealogy (-1 for root)
    int leftSon, rightSon;        // sons of node in genealogy
    double age;              // age of node
    double *conditionalProbs;      // array conditional probabilities for base assignment at node (array of length CODE_SIZE * numPatterns)
} LikelihoodNode;


/***********************************************************************************
 *	PreviousVersion
 *	- Data type which saves the unmodified version of the genealogy
 *		in case a certain genealogy change is rejected by sampler.
 ***********************************************************************************/
typedef struct PREVIOUS_VERSION {
    double dataLogLikelihood;      // original log-likelihood of data given genealogy
    int root;              // old root (if root was changed by SPR)
    unsigned short copyAll;      // this flag is turned on when all nodes are copied to saved
    unsigned short *recalcConditionals;  // array of booleans indicating for which nodes we need to recalc conditionals
    int numChangedNodes;        // number of nodes affected by proposed change
    int *changedNodeIds;        // array of ids (or indices) of nodes changed
    int numChangedConditionals;      // number of nodes whose conditional probabilities were changed
    int *changedCondIds;        // array ids (or indices) of nodes whose conditionals were changed
    LikelihoodNode **savedNodes;    // an array of pointers to previous versions.
} PreviousVersion;


/***********************************************************************************
 *	LocusSeqData
 *	- Data type which holds a summary of the sequence data for a given locus
 ***********************************************************************************/
typedef struct LOCUS_SEQ_DATA {
    int numPatterns;      // number of column patterns in alignment
    int numLivePatterns;    // number of patterns relevant for likelihood computation
    int *patternCount;      // array (of length numPatterns) of counts for each pattern
    int *numPhases;        // number of phases per pattern
    int *patternList;      // list of relevant patterns for likelihood computations (established when computing likelihood)
} LocusSeqData;


/***********************************************************************************
 *	LocusData
 *	- Data structure which holds likelihood of a locus and also saves
 *		the relevant information for quick recomputation of this likelihood
 *		given changes in the locus genealogy.
 * 	- typedef is done in LocusDataLikelihood.h
 ***********************************************************************************/
struct LOCUS_LIKELIHOOD {
    unsigned short hetMode;      // mode for computing likelihood of het alignment columns (0, 1, or 2)
    int numLeaves;          // number of leaves in genealogy
    double dataLogLikelihood;    // log-likelihood of data, given genealogy
    double mutationRate;      // relative locus-specific mutation rate
    int root;            // index of root in nodeArray[]
    LikelihoodNode **nodeArray;    // array of pointers to nodes
    LocusSeqData seqData;      // holds sequence data for locus
    PreviousVersion savedVersion;  // notes on changes proposed to genealogy

    // pointers for allocated memory
    double *doubleArray_m;
    int *intArray_m;
    LikelihoodNode *nodeArray_m;
};

/***********************************************************************************
*	LocusData
*	- Data structure which holds likelihood of a locus and also saves
*		the relevant information for quick re-computation of this likelihood
*		given changes in the locus genealogy.
***********************************************************************************/
typedef struct LOCUS_LIKELIHOOD LocusData;



/***************************************************************************************************************/
/******                               EXTERNAL FUNCTION DECLARATION                                       ******/
/***************************************************************************************************************/



/***********************************************************************************
*	createLocusData
*	- creates a new LocusData structure and allocates memory for all structures other than data.
* 	- receives as input the number of samples in the locus (genealogy leaves)
*		and the mode in which to compute likelihood of het columns in alignment:
 *		0 indicates using only first phase for each heterozygote
*		1 indicates using average likelihood of all phasings for each heterozygote
*		2 indicates using phasing with maximum likelihood for each heterozygote
* 	- returns a pointer to the structure.
***********************************************************************************/
LocusData *createLocusData(int numLeaves, unsigned short hetMode);


/***********************************************************************************
*	initializeLocusData
* 	- initializes all data structures for locus data
*	- computes leaf likelihoods for all phased patterns
* 	- if patternCounts != NULL, sets all pattern counts (otherwise, set them to 0).
*	- returns 0 if all OK, and -1 otherwise
***********************************************************************************/
int initializeLocusData(LocusData *locusData, char **patternArray, int numPatterns, int *numPhases, int *patternCounts);


/***********************************************************************************
*	freeLocusData
*	- frees all allocated memory for LocusData
* 	- returns 0
***********************************************************************************/
int freeLocusData(LocusData *locusData);


/***********************************************************************************
*	attachLeaf - UNUSED
*	- attaches a leaf to existing sub-genealogy
* 	- used when generating a starting genealogy
* 	- leafId indicates the id of the leaf being attached
* 	- target indicates the node above which leaf is attached (at a given age)
* 	- id of attachment node is set to numLeaves-1 more than id of attached leaf
* 	- returns the id of the attachment node
***********************************************************************************/
int attachLeaf_UNUSED(LocusData *locusData, int leafId, int target, double age);


/***********************************************************************************
*	setLocusMutationRate
*	- sets the mutation rate for locus to given value
***********************************************************************************/
void setLocusMutationRate(LocusData *locusData, double newRate);


/***********************************************************************************
*	getLocusMutationRate
*	- returns the mutation rate for locus
***********************************************************************************/
double getLocusMutationRate(LocusData *locusData);


/***********************************************************************************
*	computeAllConditionals
*	- computes all conditional likelihoods at all nodes of tree under all patterns
*	- calls recursive procedure computeConditionalJC to recompute conditional probabilities
*		and override all old versions.
*	- also initializes all pattern counts to zero
*	- returns 0 
***********************************************************************************/
int computeAllConditionals(LocusData *locusData);


/***********************************************************************************
*	computeLocusDataLikelihood
*	- computes log-likelihood of data at a given locus, given its genealogy
*	- calls recursive procedure computeConditionalJC to recompute conditional probabilities
*		and makes sure to save old versions.
*	- if useOldConditionals == 1, uses previously computed conditionals, when possible.
*	- otherwise, recomputes everything from scratch
*	- returns the log-likelihood
***********************************************************************************/
double computeLocusDataLikelihood(LocusData *locusData, unsigned short useOldConditionals);


/***********************************************************************************
*	computePatternLogLikelihood
*	- receives a list of patterns and their counts and computes their log likelihood 
*	- uses pre-computed conditionals and doesn't modify anything in the data structure.
*	- returns the log likelihood
***********************************************************************************/
double computePatternLogLikelihood(LocusData *locusData, int numPatterns, int *patternIds, int *patternCounts);


/***********************************************************************************
*	!!!!! FOR DEBUGGING !!!!!
***********************************************************************************/
double computeLocusDataLikelihood_deb(LocusData *locusData, unsigned short useOldConditionals);


/***********************************************************************************
*	addSitePatterns
*	- adds a set of site patterns to live pattern set
*	- receives ids of patterns to add and their respective counts
*	- if revertToSaved == 1, then conditionals do not need to be computed (since they already exist),
*		and new likelihood is the saved one.
*	- returns delta in log likelihood of this step
***********************************************************************************/
double addSitePatterns(LocusData *locusData, int numPatterns, int *patternIds, int *patternCounts,
                       unsigned short revertToSaved);


/***********************************************************************************
*	reduceSitePatterns
*	- reduces the counts of a set of site patterns in likelihood computation
*	- receives ids of patterns to reduce and the respective counts
*	- if revertToSaved == 1, uses saved likelihood to compute new likelihood
*	- returns delta in log likelihood of this step
***********************************************************************************/
double reduceSitePatterns(LocusData *locusData, int numPatterns, int *patternIds, int *patternCounts,
                          unsigned short revertToSaved);


/***********************************************************************************
*	checkLocusDataLikelihood
*	- re-computes log-likelihood of data at a given locus, given its genealogy
*	- if inconsistency is found in computed log likelihood, checks inconsistencies in 
*		all recorded conditional probabilities
*	- returns 1 if all is OK, and 0 if inconsistencies were found
***********************************************************************************/
int checkLocusDataLikelihood(LocusData *locusData);


/***********************************************************************************
*	revertToSaved
*	- reverts locus data structure (genealogy and conditional likelihoods) to saved version
*	- returns 0
***********************************************************************************/
int revertToSaved(LocusData *locusData);


/***********************************************************************************
*	resetSaved
*	- resets saved data and maintains updates
*	- returns 0
***********************************************************************************/
int resetSaved(LocusData *locusData);


/***********************************************************************************
*	adjustGenNodeAge
*	- adjusts nodes age and saves old version
*	- nodeId is id of node and age is new age
*	- returns 0
***********************************************************************************/
int adjustGenNodeAge(LocusData *locusData, int nodeId, double age);


/***********************************************************************************
*	scaleAllNodeAges
*	- scales all node ages in genealogy with a given multiplicative factor
*	- returns the delta in log-likelihood of suggested step
*	- saves all original ages and conditionals
***********************************************************************************/
double scaleAllNodeAges(LocusData *locusData, double factor);


/***********************************************************************************
*	executeGenSPR
*	- executes an SPR operation on locus genealogy
*	- subtree_root is node id for root of pruned subtree,targetBranch is id of node below
*		branch where subtree should be regrafted, and age indicates time of regrafting.
*	- returns 0 if root node remains the same
*	- otherwise, returns 1 if subtree is a child of the root AFTER regrafting
*	- otherwise, returns 2 if subtree was a child of the root BEFORE regrafting
***********************************************************************************/
int executeGenSPR(LocusData *locusData, int subtreeRoot, int targetBranch, double age);


/***********************************************************************************
*	copyGenericTreeToLocus
* 	- copies a generic tree into likelihood tree
*	- assumes label1 holds age of node
*	- returns 0
***********************************************************************************/
int copyGenericTreeToLocus(LocusData *locusData, GenericBinaryTree *genericTree);


/***********************************************************************************
*	printLocusGenTree
* 	- prints the genealogy tree (including supplied population and event id)
*	- prints each node in a separate line, according to id order
*	- if there are nodes which are considered for changes, print out their id's
***********************************************************************************/
void printLocusGenTree(LocusData *locusData, FILE *stream, int *nodePops, int *nodeEvents);


/***********************************************************************************
*	printLocusDataStats
* 	- prints stats on alignment patterns
*	- stats are outputted in one line (with newline) to stdout in the following order:
*	- num hom patterns, num het patterns with 2 phases, , num het patterns with 4 phases...
*	- for each of the above stats, prints consecutively the number of distinct patterns 
*		and the number of columns corresponding to it in the alignment
*	- maxLogPhases is an upper bound on the number of phased hets per pattern
***********************************************************************************/
void printLocusDataStats(LocusData *locusData, int maxLogPhases);


/***********************************************************************************
*	printLocusDataPatterns
* 	- prints the alignment of the locus, pattern by pattern to output file.
*	- patterns are printed in columns, according to leaf id.
*	- below each pattern is its multiplicity in the alignment
*	- het patterns are represented by all phasings 
*		(mult is written below first phasing)
***********************************************************************************/
void printLocusDataPatterns(LocusData *locusData, FILE *outFile);


/***********************************************************************************
 *	computePairwiseLCAs
 *	- procedure for computing a 2D matrix with the ids of the LCAs (least
 *    common ancestors) of all pairs of leaves
 *  - returns 1 if successful, 0 otherwise
 ***********************************************************************************/
int computePairwiseLCAs(LocusData *locusData, int **lcaMatrix, int *leafArray_aux);


/***********************************************************************************
 *	getSortedAges
 *	- procedure for computing a sorted list of node ages (from most recent to root)
 *  - assumes array given as input has 2x space for all internal nodes (also for auxiliary space for sorting)
 *  - calls recursive procedure getSortedAges_rec on root
 *  - returns 1 if successful, 0 otherwise
 ***********************************************************************************/
int getSortedAges(LocusData *locusData, double *ageArray);



/*  a sequence of "get functions" for various attributes - AVOID USING !!!   */



/***********************************************************************************
*	getLocusDataLikelihood
*	- returns the log likelihood of the locus as last recorded (no new computations)
***********************************************************************************/
double getLocusDataLikelihood(LocusData *locusData);


/***********************************************************************************
*	getLocusRoot
*	- returns the root node id
***********************************************************************************/
int getLocusRoot(LocusData *locusData);


/***********************************************************************************
*	getNodeAge
*	- returns the age of a node
***********************************************************************************/
double getNodeAge(LocusData *locusData, int nodeId);


/***********************************************************************************
*	getNodeFather
*	- returns the id of the father of a node
***********************************************************************************/
int getNodeFather(LocusData *locusData, int nodeId);


/***********************************************************************************
*	getNodeSon
*	- returns the id of a son of a node (son = 0 -> left, son = 1 -> right)
***********************************************************************************/
int getNodeSon(LocusData *locusData, int nodeId, unsigned short son);



/***************************************************************************************************************/
/******                                        END OF FILE                                                ******/
/***************************************************************************************************************/
#endif
