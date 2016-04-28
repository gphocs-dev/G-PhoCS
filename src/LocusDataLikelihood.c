/** 
   \file LocusDataLikelihood.c 
   DataStructures to hold all locus sequence data & functions to compute likelihood given specific genealogy.	

   Contains the relevant data structures and procedure implementations
   for holding all locus sequence data and computing likelihood of data 
   given a specific genealogy.
*/



#include "LocusDataLikelihood.h"
#include "utils.h"
#include <stdlib.h>
#include <math.h>
#include <string.h>


/***************************************************************************************************************/
/******                                       INTERNAL CONSANTS                                           ******/
/***************************************************************************************************************/



#define CODE_SIZE	4



/***************************************************************************************************************/
/******                                         INTERNAL DATA TYPES                                       ******/
/***************************************************************************************************************/

double global_prob_val;

/***********************************************************************************
 *	LocusSeqData
 *	- Data type which holds a summary of the sequence data for a given locus
 ***********************************************************************************/
typedef struct LOCUS_SEQ_DATA {
  int numPatterns;			// number of column patterns in alignment
  int numLivePatterns;		// number of patterns relevant for likelihood computation
  int* patternCount;			// array (of length numPatterns) of counts for each pattern
  int* numPhases;				// number of phases per pattern
  int* patternList;			// list of relevant patterns for likelihood computations (established when computing likelihood)
} LocusSeqData;



/***********************************************************************************
 *	LikelihoodNode
 *	- Data type which holds data for computing conditional probability
 *		of all base assignments to a certain genealogy (coalescent) node.
 ***********************************************************************************/
typedef struct LIKELIHOOD_NODE {
  int father;							// father of node in genealogy (-1 for root)
  int leftSon, rightSon;				// sons of node in genealogy
  double age;							// age of node
  double* conditionalProbs;				// array conditional probabilities for base assignment at node (array of length CODE_SIZE * numPatterns)
} LikelihoodNode;



/***********************************************************************************
 *	PreviousVersion
 *	- Data type which saves the unmodified version of the genealogy
 *		in case a certain genealogy change is rejected by sampler.
 ***********************************************************************************/
typedef struct PREVIOUS_VERSION {
  double dataLogLikelihood;			// original log-likelihood of data given genealogy
  int root;							// old root (if root was changed by SPR)
  unsigned short	copyAll;			// this flag is turned on when all nodes are copied to saved
  unsigned short* recalcConditionals;	// array of booleans indicating for which nodes we need to recalc conditionals
  int numChangedNodes;				// number of nodes affected by proposed change
  int* changedNodeIds;				// array of ids (or indices) of nodes changed
  int numChangedConditionals;			// number of nodes whose conditional probabilities were changed
  int* changedCondIds;				// array ids (or indices) of nodes whose conditionals were changed
  LikelihoodNode** savedNodes;		// an array of pointers to previous versions.
} PreviousVersion;



/***********************************************************************************
 *	LocusData
 *	- Data structure which holds likelihood of a locus and also saves
 *		the relevant information for quick recomputation of this likelihood
 *		given changes in the locus genealogy.
 * 	- typedef is done in LocusDataLikelihood.h
 ***********************************************************************************/
struct LOCUS_LIKELIHOOD {
  unsigned short hetMode;		// mode for computing likelihood of het alignment columns (0, 1, or 2)
  int numLeaves;				// number of leaves in genealogy
  double dataLogLikelihood;		// log-likelihood of data, given genealogy
  double mutationRate;			// relative locus-specific mutation rate
  int root;						// index of root in nodeArray[]
  LikelihoodNode** nodeArray;	// array of pointers to nodes
  LocusSeqData seqData;			// holds sequence data for locus
  PreviousVersion savedVersion;	// notes on changes proposed to genealogy
	
  // pointers for allocated memory
  double* doubleArray_m;
  int* intArray_m;
  LikelihoodNode* nodeArray_m;
};


/***************************************************************************************************************/
/******                                INTERNAL FUNCTION DECLARATIONS                                     ******/
/***************************************************************************************************************/



int computeConditionalJC_new (LocusData* locusData, int nodeId, int numPatterns, int* patternIds, unsigned short overideOld);
int computeConditionalJC (LocusData* locusData, int nodeId, int numPatterns, int* patternIds, unsigned short overideOld);
double	computeEdgeConditionalJC(double edgeLength);
int	copyNodeToSaved(LocusData* locusData, int nodeId, unsigned short recalcConditionals);
int	copyNodeConditionals(LocusData* locusData, int nodeId);
int computeLeafConditionals(LocusData* locusData, char* patternString);
void computeSubtreeConditionals (double* sonConditionals, double* parentConditionals, double* edgeConditionals);
void computeSubtreeConditionals_new (double* sonConditionals, double* parentConditionals, double* edgeSubstProb);
int computePairwiseLCAs_rec (LocusData* locusData, int nodeId, int** lcaMatrix, int* leafArray, int arrayOffset, int* numLeaves_out);
int getSortedAges_rec (LocusData* locusData, int nodeId, double* sortedAges, double* sortedAges_aux, int arrayOffset, int* numInternalNodes_out);



/***************************************************************************************************************/
/******                              EXTERNAL FUNCTION IMPLEMENTATION                                     ******/
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
LocusData* createLocusData (int numLeaves, unsigned short hetMode) {
  LocusData* locusData;
  int* intArray;
  int node, numNodes = 2*numLeaves-1;
	
//  printf("Creating locus with %d leaves and %d nodes.\n", numLeaves, numNodes);
  
  locusData = (LocusData*)malloc(sizeof(LocusData));
  if(locusData == NULL) {
    fprintf(stderr, "\nError: Out Of Memory when allocating space for locusData ");
    return NULL;
  }

  intArray = (int*)malloc(2*numNodes*sizeof(int));
  if(intArray == NULL) {
    fprintf(stderr, "\nError: Out Of Memory when allocating space for locusData array of integers (for changed node ids) ");
    return NULL;
  }
  locusData->savedVersion.changedNodeIds = intArray;
  locusData->savedVersion.changedCondIds = intArray+numNodes;

 	
  locusData->savedVersion.recalcConditionals = (unsigned short*)malloc(numNodes*sizeof(unsigned short));
  if(locusData->savedVersion.recalcConditionals == NULL) {
    fprintf(stderr, "\nError: Out Of Memory when allocating space for locusData recalcConditionals boolean array ");
    return NULL;
  }

	
  // every vertex has two copies - for genealogy changes
  locusData->nodeArray_m = (LikelihoodNode*)malloc(2*numNodes*sizeof(LikelihoodNode));
  if(locusData->nodeArray_m == NULL) {
    fprintf(stderr, "\nError: Out Of Memory when allocating space for locusData array of nodes ");
    return NULL;
  }

	
  locusData->nodeArray = (LikelihoodNode**)malloc((numNodes*sizeof(LikelihoodNode*)));
  if(locusData->nodeArray == NULL) {
    fprintf(stderr, "\nError: Out Of Memory when allocating space for locusData array of node pointers ");
    return NULL;
  }

	
  locusData->savedVersion.savedNodes = (LikelihoodNode**)malloc((numNodes*sizeof(LikelihoodNode*)));
  if(locusData->savedVersion.savedNodes == NULL) {
    fprintf(stderr, "\nError: Out Of Memory when allocating space for locusData array of saved node pointers ");
    return NULL;
  }

	
  //set up entries of locusData
  locusData->hetMode = hetMode;
  locusData->numLeaves = numLeaves;
  locusData->mutationRate = 1.0;
  locusData->root = -1;
  locusData->doubleArray_m = NULL;
  locusData->dataLogLikelihood = 0.0;
  locusData->savedVersion.dataLogLikelihood = 0.0;

  locusData->seqData.numPatterns = 0;

  // initialize node data structures (other than conditional array
  for(node=0; node<2*numNodes; node++) {
    locusData->nodeArray_m[node].conditionalProbs = NULL;
    locusData->nodeArray_m[node].age = 0.0;
    locusData->nodeArray_m[node].father = -1;
    locusData->nodeArray_m[node].leftSon = -1;
    locusData->nodeArray_m[node].rightSon = -1;
  }
	
  // initialize node pointers
  for(node=0; node<numNodes; node++) {
    locusData->nodeArray[node] = &locusData->nodeArray_m[2*node];			
    locusData->savedVersion.savedNodes[node] = &locusData->nodeArray_m[2*node + 1];			
  }
	
	

  // reset saved version for locus data
  resetSaved(locusData);

  global_prob_val = ((1-exp(-4*0.0000001/3.0)) / 4.0 );

  return locusData;
}
/** end of createLocusData **/



/***********************************************************************************
 *	initializeLocusData
 * 	- initializes all data structures for locus data
 *	- computes leaf likelihoods for all phased patterns
 * 	- if patternCounts != NULL, sets all pattern counts (otherwise, set them to 0).
 *	- returns 0 if all OK, and -1 otherwise
 ***********************************************************************************/
int initializeLocusData(LocusData* locusData, char** patternArray, int numPatterns, int* numPhases, int* patternCounts)	{
	
  int node, patt, unphasedPatt;
	
  // auxiliary arrays
  char *patternString;
	
  if(locusData == NULL)		return -1;
	
  //	printf("Initializing locus data likelihood with %d patterns.\n",numPatterns);
	
  // allocate seqData memory (pattern frequencies, conditional arrays, and numPhases)
  locusData->doubleArray_m = (double*)malloc(2*(2*locusData->numLeaves-1)*CODE_SIZE*numPatterns*sizeof(double));
  if(locusData->doubleArray_m == NULL) {
    fprintf(stderr, "\nError: Out Of Memory when allocating space for locusData array of doubles (for conditional probabilities of genealogy nodes) in initializeLocusData().\n");
    return -1;
  }
	
  locusData->intArray_m = (int*)malloc(numPatterns*3*sizeof(int));
  if(locusData->intArray_m == NULL) {
    fprintf(stderr, "\nError: Out Of Memory when alloating space for locusData->intArray_m in initializeLocusData().\n");
    return -1;
  }
  locusData->seqData.numPhases = locusData->intArray_m;
  locusData->seqData.patternList = locusData->intArray_m + numPatterns;
  locusData->seqData.patternCount = locusData->intArray_m + 2*numPatterns;

  for(node=0; node < 2*locusData->numLeaves-1; node++) {
    locusData->nodeArray[node]->conditionalProbs = locusData->doubleArray_m + (2*node)*numPatterns*CODE_SIZE;
    locusData->savedVersion.savedNodes[node]->conditionalProbs = locusData->doubleArray_m + (2*node+1)*numPatterns*CODE_SIZE;
  }

  // initialize leaf conditionals for hom patterns
  locusData->seqData.numPatterns = 0;
  locusData->seqData.numLivePatterns = 0;
  unphasedPatt = 0;
  for(patt=0; patt<numPatterns; patt++) {
    locusData->seqData.numPhases[patt] = numPhases[patt];
    patternString = patternArray[patt];
    if(patternCounts != NULL && numPhases[patt]> 0) {
      locusData->seqData.patternCount[patt] = patternCounts[unphasedPatt];
      unphasedPatt++;
    } else {
      locusData->seqData.patternCount[patt] = 0;
    }
    if(0 > computeLeafConditionals(locusData, patternString)) {
      fprintf(stderr, "Error: Error while computing conditionals for Hom pattern #%d: ",patt+1);
      for(node=0; node<locusData->numLeaves; node++) {
        fprintf(stderr, "%c",patternString[node]);
      }
      fprintf(stderr, "\n");
      return -1;
    }
  }
  if(patternCounts != NULL) {
    locusData->seqData.numLivePatterns = numPatterns;
  } else {
    locusData->seqData.numLivePatterns = 0;
  }		
	
  //	printLocusDataPatterns(locusData,stdout);
	
  return 0;
}	
/** end of initializeLocusData **/



/***********************************************************************************
 *	freeLocusData
 *	- frees all allocated memory for LocusData
 * 	- returns 0
 ***********************************************************************************/
int freeLocusData (LocusData* locusData) {
	
  if(locusData->doubleArray_m != NULL) free(locusData->doubleArray_m);
  if(locusData->intArray_m != NULL) free(locusData->intArray_m);
  free(locusData->nodeArray);
  free(locusData->nodeArray_m);
  free(locusData->savedVersion.savedNodes);
  free(locusData->savedVersion.recalcConditionals);
  free(locusData->savedVersion.changedNodeIds);
  free(locusData);
	
  return 0;
}
/** end of freeLocusData **/



/***********************************************************************************
 *	attachLeaf - UNUSED
 *	- attaches a leaf to existing sub-genealogy
 * 	- used when generating a starting genealogy
 * 	- leafId indicates the id of the leaf being attached
 * 	- target indicates the node above which leaf is attached (at a given age)
 * 	- id of attachment node is set to numLeaves-1 more than id of attached leaf
 * 	- returns the id of the attachment node
 ***********************************************************************************/
/* MARK: ADD ARGUMENT FOR LEAF AGE IF WE WANT TO USE ANCIENT SAMPLES FOR LEAVES
*/
int attachLeaf_UNUSED (LocusData* locusData, int leafId, int target, double age)  {
  int newnodeId = leafId + locusData->numLeaves-1;
  LikelihoodNode* node = locusData->nodeArray[newnodeId];

  node->age = age;
  node->father = locusData->nodeArray[target]->father;
  node->leftSon = leafId;
  node->rightSon = target;
  locusData->nodeArray[leafId]->father = newnodeId;
  locusData->nodeArray[target]->father = newnodeId;
  if(node->father < 0) {
    locusData->root = newnodeId;
  } else if(locusData->nodeArray[node->father]->leftSon == target) {
    locusData->nodeArray[node->father]->leftSon = newnodeId;
  } else {
    locusData->nodeArray[node->father]->rightSon = newnodeId;			
  }
	
  return newnodeId;
}
/** end of attachLeaf **/



/***********************************************************************************
 *	setLocusMutationRate
 *	- sets the mutation rate for locus to given value
 ***********************************************************************************/
void setLocusMutationRate (LocusData* locusData, double newRate)	{
  locusData->mutationRate = newRate;
  return;
}
/** end of setLocusMutationRate **/



/***********************************************************************************
 *	getLocusMutationRate
 *	- returns the mutation rate for locus
 ***********************************************************************************/
double getLocusMutationRate (LocusData* locusData)	{
  return locusData->mutationRate;
}
/** end of getLocusMutationRate **/


/***********************************************************************************
 *	computeAllConditionals
 *	- computes all conditional likelihoods at all nodes of tree under all patterns
 *	- calls recursive procedure computeConditionalJC to recompute conditional probabilities
 *		and override all old versions.
 *	- also initializes all pattern counts to zero
 *	- returns 0 
 ***********************************************************************************/
int computeAllConditionals (LocusData* locusData)  {
  int  pattId, phase, numLivePatterns;

  // set to compute likelihood under all patterns 
  for(pattId=0, numLivePatterns=0; pattId<locusData->seqData.numPatterns; pattId++) {
    locusData->seqData.patternCount[pattId] = 0;
    for(phase=0; phase<locusData->seqData.numPhases[pattId]; phase++) {
      locusData->seqData.patternList[numLivePatterns] = pattId+phase;
      numLivePatterns++;
    }
  }
#ifdef OPT1	
  computeConditionalJC_new(locusData, locusData->root, numLivePatterns, locusData->seqData.patternList,/*overrideOld=*/ 1); 
#else
  computeConditionalJC(locusData, locusData->root, numLivePatterns, locusData->seqData.patternList,/*overrideOld=*/ 1); 
#endif
  return 0;
}
/** end of computeAllConditionals **/



/***********************************************************************************
 *	computeLocusDataLikelihood
 *	- computes log-likelihood of data at a given locus, given its genealogy
 *	- calls recursive procedure computeConditionalJC to recompute conditional probabilities
 *		and makes sure to save old versions.
 *	- if useOldConditionals == 1, uses previously computed conditionals, when possible.
 *	- otherwise, recomputes everything from scratch
 *	- returns the log-likelihood
 ***********************************************************************************/
double computeLocusDataLikelihood (LocusData* locusData, unsigned short useOldConditionals)  {
  int res, node;
  int  patt, pattId, phase, numLivePatterns, conditional, numConditionals;
  double prob;
	
  if(locusData->seqData.numLivePatterns == 0) return 0.0;
	
  if(!useOldConditionals) {
    for(node = locusData->numLeaves; node < 2*locusData->numLeaves-1; node++) {
      copyNodeConditionals(locusData,node);
    }
  }
	
  //	printf("saving old likelihood %g.\n",locusData->dataLogLikelihood);
  locusData->savedVersion.dataLogLikelihood = locusData->dataLogLikelihood;

  for(pattId=0, numLivePatterns=0; pattId<locusData->seqData.numPatterns; pattId++) {
    if(locusData->seqData.patternCount[pattId] > 0) {
      for(phase=0; phase<locusData->seqData.numPhases[pattId]; phase++) {
        locusData->seqData.patternList[numLivePatterns] = pattId+phase;
        numLivePatterns++;
      }
    }
  }
  if(numLivePatterns != locusData->seqData.numLivePatterns) {
    fprintf(stderr, "Error: there should be %d live patterns and there are %d.\n",locusData->seqData.numLivePatterns, numLivePatterns);
    exit(-1);
  } else {
    //		printf("%d live patterns.\n",numLivePatterns);
  }
	
#ifdef OPT1	
res = computeConditionalJC_new(locusData, locusData->root, numLivePatterns, locusData->seqData.patternList, !useOldConditionals);
#else
res = computeConditionalJC(locusData, locusData->root, numLivePatterns, locusData->seqData.patternList, !useOldConditionals);
#endif
	

  if(!res)	return locusData->dataLogLikelihood;
	
  locusData->dataLogLikelihood = 0.0; 
	
  //	printf("Locus likelihood computation:\n");
	
  // sum over root conditionals assuming uniform distribution at root
  for(patt=0; patt<numLivePatterns; patt+=locusData->seqData.numPhases[pattId]) {
    pattId = locusData->seqData.patternList[patt];
    prob = 0.0;
    numConditionals = CODE_SIZE*locusData->seqData.numPhases[pattId];
    for(conditional=0; conditional<numConditionals; conditional++) {
      prob += locusData->nodeArray[ locusData->root ]->conditionalProbs[pattId*CODE_SIZE+conditional];
    }
    locusData->dataLogLikelihood += log(prob/numConditionals) * locusData->seqData.patternCount[pattId];
  }
	
  //	printf("new likelihood is %g.\n",locusData->dataLogLikelihood);
  return locusData->dataLogLikelihood;
}
/** end of computeLocusDataLikelihood **/



/***********************************************************************************
 *	computePatternLogLikelihood
 *	- receives a list of patterns and their counts and computes their log likelihood 
 *	- uses pre-computed conditionals and doesn't modify anything in the data structure.
 *	- returns the log likelihood
 ***********************************************************************************/
double computePatternLogLikelihood (LocusData* locusData, int numPatterns, int* patternIds, int* patternCounts)  {

  int patt, pattId, numConditionals, conditional;
  double prob, logLikelihood;
	
  // sum over root conditionals assuming uniform distribution at root
  logLikelihood = 0.0;
  for(patt=0; patt<numPatterns; patt++) {
    pattId = patternIds[patt];
    prob = 0.0;
    numConditionals = CODE_SIZE*locusData->seqData.numPhases[pattId];
    for(conditional=0; conditional<numConditionals; conditional++) {
      prob += locusData->nodeArray[ locusData->root ]->conditionalProbs[pattId*CODE_SIZE+conditional];
    }
    logLikelihood += log(prob/numConditionals) * patternCounts[patt];
  }
	
  return logLikelihood;


}
/** end of computePatternLogLikelihood **/


/***********************************************************************************
 *	!!!!!FOR DEBUGGING !!!!!
 ***********************************************************************************/
double computeLocusDataLikelihood_deb (LocusData* locusData, unsigned short useOldConditionals)  {
  int res, node;
  int  patt, pattId, phase, numLivePatterns, conditional, numConditionals;
  double prob;
	
  if(locusData->seqData.numLivePatterns == 0) return 0.0;
	
  if(!useOldConditionals) {
    for(node = locusData->numLeaves; node < 2*locusData->numLeaves-1; node++) {
      copyNodeConditionals(locusData,node);
    }
  }
	
  //	printf("saving old likelihood %g.\n",locusData->dataLogLikelihood);
  locusData->savedVersion.dataLogLikelihood = locusData->dataLogLikelihood;

  for(pattId=0, numLivePatterns=0; pattId<locusData->seqData.numPatterns; pattId++) {
    if(locusData->seqData.patternCount[pattId] > 0) {
      for(phase=0; phase<locusData->seqData.numPhases[pattId]; phase++) {
        locusData->seqData.patternList[numLivePatterns] = pattId+phase;
        numLivePatterns++;
      }
    }
  }
  if(numLivePatterns != locusData->seqData.numLivePatterns) {
    fprintf(stderr, "Error: there should be %d live patterns and there are %d.\n",locusData->seqData.numLivePatterns, numLivePatterns);
    return -1;
  }
	
#ifdef OPT1	
  res = computeConditionalJC_new(locusData, locusData->root, numLivePatterns, locusData->seqData.patternList, !useOldConditionals);
#else
  res = computeConditionalJC(locusData, locusData->root, numLivePatterns, locusData->seqData.patternList, !useOldConditionals);
#endif
	

  if(!res)	return locusData->dataLogLikelihood;
	
  locusData->dataLogLikelihood = 0.0; 
	
  //	printf("Locus likelihood computation:\n");
	
  // sum over root conditionals assuming uniform distribution at root
  for(patt=0; patt<numLivePatterns; patt+=locusData->seqData.numPhases[pattId]) {
    pattId = locusData->seqData.patternList[patt];
    prob = 0.0;
    numConditionals = CODE_SIZE*locusData->seqData.numPhases[pattId];
    printf("pattern %d accumulative conditional:",pattId+1);
    for(conditional=0; conditional<numConditionals; conditional++) {
      prob += locusData->nodeArray[ locusData->root ]->conditionalProbs[pattId*CODE_SIZE+conditional];
      printf(" %g",locusData->nodeArray[ locusData->root ]->conditionalProbs[pattId*CODE_SIZE+conditional]);
    }
    printf("\n");
    locusData->dataLogLikelihood += log(prob/numConditionals) * (double)locusData->seqData.patternCount[pattId];
  }
	
  //	printf("new likelihood is %g.\n",locusData->dataLogLikelihood);
  return locusData->dataLogLikelihood;
}
/** end of computeLocusDataLikelihood **/



/***********************************************************************************
 *	addSitePatterns
 *	- adds a set of site patterns to live pattern set
 *	- receives ids of patterns to add and their respecitve counts
 *	- if revertToSaved == 1, then conditionals do not need to be computed (since they already exist),
 *		and new likelihood is the saved one.
 *	- returns delta in log likelihood of this step
 ***********************************************************************************/
double addSitePatterns (LocusData* locusData, int numPatterns, int* patternIds, int* patternCounts, unsigned short revertToSaved)  {
  int patt, pattId, phase;
  int numNewPatterns;
  int* newPatterns = locusData->seqData.patternList;		// use this pre-allocated space to list
  double deltaLogLikelihood;
	

  if(numPatterns == 0) {
    locusData->savedVersion.dataLogLikelihood = locusData->dataLogLikelihood;
    return 0;
  }
	
  // compute new counts of affected patterns
  // and all phased versions of newly introduced patterns
  numNewPatterns = 0;
  for(patt=0; patt<numPatterns; patt++) {
    pattId = patternIds[patt];
    if(locusData->seqData.patternCount[pattId] == 0) {
      for(phase=0; phase<locusData->seqData.numPhases[pattId]; phase++) {
        newPatterns[numNewPatterns] = pattId+phase;
        numNewPatterns++;
      }
    }
    locusData->seqData.patternCount[pattId] += patternCounts[patt];
  }
	
  locusData->seqData.numLivePatterns += numNewPatterns;
	
  if(revertToSaved) {
    // revert to saved version
    deltaLogLikelihood = locusData->savedVersion.dataLogLikelihood - locusData->dataLogLikelihood;
    locusData->dataLogLikelihood = locusData->savedVersion.dataLogLikelihood;
    return deltaLogLikelihood;
  }
	
  // save current log likelihood (in case we need to revert)
  locusData->savedVersion.dataLogLikelihood = locusData->dataLogLikelihood;

  if(numNewPatterns > 0) {		
#ifdef OPT1	
    computeConditionalJC_new(locusData, locusData->root, numNewPatterns, newPatterns, /*overide old conditionals*/ 1);
#else
    computeConditionalJC(locusData, locusData->root, numNewPatterns, newPatterns, /*overide old conditionals*/ 1);
#endif
  }
	
	
  deltaLogLikelihood = computePatternLogLikelihood(locusData, numPatterns, patternIds, patternCounts);

  locusData->dataLogLikelihood += deltaLogLikelihood;
	
  return deltaLogLikelihood;
	
	
}
/** end of addSitePatterns **/



/***********************************************************************************
 *	reduceSitePatterns
 *	- reduces the counts of a set of site patterns in likelihood computation
 *	- receives ids of patterns to reduce and the respective counts
 *	- if revertToSaved == 1, uses saved likelihood to compute new likelihood
 *	- returns delta in log likelihood of this step
 ***********************************************************************************/
double reduceSitePatterns (LocusData* locusData, int numPatterns, int* patternIds, int* patternCounts, unsigned short revertToSaved)  {
  int patt, pattId;
  int numRemovedPatterns;
  double deltaLogLikelihood;
	

  if(numPatterns == 0) {
    locusData->savedVersion.dataLogLikelihood = locusData->dataLogLikelihood;
    return 0;
  }
	
  // compute new counts of affected patterns
  // and all phased versions of newly introduced patterns
  numRemovedPatterns = 0;
  for(patt=0; patt<numPatterns; patt++) {
    pattId = patternIds[patt];
    locusData->seqData.patternCount[pattId] -= patternCounts[patt];
    if(locusData->seqData.patternCount[pattId] == 0) {
      numRemovedPatterns+=locusData->seqData.numPhases[pattId];
    } else if(locusData->seqData.patternCount[pattId] < 0) {
      fprintf(stderr, "Error: Error in removing site patterns from likelihood computation. Pattern %d has negative count %d.\n",
             pattId+1, locusData->seqData.patternCount[pattId]);
      exit(-1);
    }
  }
	
  locusData->seqData.numLivePatterns -= numRemovedPatterns;
	
  if(revertToSaved) {
    // revert to saved version
    deltaLogLikelihood = locusData->savedVersion.dataLogLikelihood - locusData->dataLogLikelihood;
    locusData->dataLogLikelihood = locusData->savedVersion.dataLogLikelihood;
    return deltaLogLikelihood;
  }
	
  // save current log likelihood (in case we need to revert)
  locusData->savedVersion.dataLogLikelihood = locusData->dataLogLikelihood;
	
  deltaLogLikelihood = - computePatternLogLikelihood(locusData, numPatterns, patternIds, patternCounts);
	
  locusData->dataLogLikelihood += deltaLogLikelihood;
	
  return deltaLogLikelihood;
	
	
}
/** end of reduceSitePatterns **/



/***********************************************************************************
 *	checkLocusDataLikelihood
 *	- re-computes log-likelihood of data at a given locus, given its genealogy
 *	- if inconsistency is found in computed log likelihood, checks inconsistencies in 
 *		all recorded conditional probabilities
 *	- returns 1 if all is OK, and 0 if inconsistencies were found
 ***********************************************************************************/
int checkLocusDataLikelihood (LocusData* locusData) {
  int node, patt, conditional, numConditionals;
  double *savedConds, *newConds;
	
  computeLocusDataLikelihood (locusData,/*do not use old conditionals*/ 0);
	
  // if likelihood is compatible, no need for further checks
  if(locusData->dataLogLikelihood == locusData->savedVersion.dataLogLikelihood || fabs(1- locusData->dataLogLikelihood/locusData->savedVersion.dataLogLikelihood) < 0.000000001) {
    resetSaved(locusData);
    return 1;
  }

	
  printf("\nInconsistent locus log-likelihood (saved %g, recomputed %g, diff %g).\n", 
         locusData->savedVersion.dataLogLikelihood,locusData->dataLogLikelihood,locusData->savedVersion.dataLogLikelihood-locusData->dataLogLikelihood);
  printf("Checking conditionals...\n");

  for(node=0; node<2*locusData->numLeaves-1; node++) {
    savedConds = locusData->savedVersion.savedNodes[node]->conditionalProbs;
    newConds   = locusData->nodeArray[node]->conditionalProbs;
    for(	patt=0; 
            patt<locusData->seqData.numPatterns; 
            patt+=locusData->seqData.numPhases[patt], 
              savedConds+=numConditionals, 
              newConds+=numConditionals) {
			
      numConditionals = CODE_SIZE*locusData->seqData.numPhases[patt];
      if(locusData->seqData.patternCount[patt] == 0)		continue;
      for(conditional=0; conditional<numConditionals; conditional++) {
        if( newConds[conditional] != savedConds[conditional]) {
          printf("Inconsistent conditionals in node %d, patt %d, phased base %d (saved %g, recomputed %g).\n", 
                 node, patt, conditional, savedConds[conditional], newConds[conditional]);
        }
      }
    }
  }


  resetSaved(locusData);
  return 0;
}	
/** end of checkLocusDataLikelihood **/



/***********************************************************************************
 *	revertToSaved
 *	- reverts locus data structure (genealogy and conditional likelihoods) to saved version
 *	- returns 0
 ***********************************************************************************/
int revertToSaved(LocusData* locusData) {
  int nodeId, i;
  LikelihoodNode* tmp_node;
  LikelihoodNode** tmp_nodeArray;
  double* conditionalPointer;
	
	
  //	printf("Reverting to saved version of locus\n");
	
  // copy old likelihood
  locusData->dataLogLikelihood = locusData->savedVersion.dataLogLikelihood;
	
  // replace root, if necessary
  if(locusData->savedVersion.root >= 0) {
    locusData->root = locusData->savedVersion.root;
    locusData->savedVersion.root = -1;
  }

  // if all were copied, switch all back
  if(locusData->savedVersion.copyAll) {
    //		printf("\ncopy all");
    tmp_nodeArray = locusData->nodeArray;
    locusData->nodeArray = locusData->savedVersion.savedNodes;
    locusData->savedVersion.savedNodes = tmp_nodeArray;
    resetSaved(locusData);		
    return 0;
  }
		
  //	printf("Reverting to saved version of locus (2)\n");
  if(locusData->savedVersion.numChangedConditionals == 0 && locusData->savedVersion.numChangedNodes == 0) {
    //		printf("\nNo changed conditionals");
    return 0;
  }
  //	printf("Reverting to saved version of locus (3)\n");
	
  // replace all changed nodes with saved versions	
  for(i=0; i<locusData->savedVersion.numChangedNodes; i++) {
    nodeId = locusData->savedVersion.changedNodeIds[i];
    //		printf("Reverting node %d\n",nodeId);
    // switch entire pointer to node likelihood struct
    tmp_node = locusData->nodeArray[nodeId];
    locusData->nodeArray[nodeId] = locusData->savedVersion.savedNodes[nodeId];
    locusData->savedVersion.savedNodes[nodeId] = tmp_node;
    if(locusData->savedVersion.recalcConditionals[nodeId]) {
      // if conditionals were updated, do nothing
      locusData->savedVersion.recalcConditionals[nodeId] = 0;
    } else {
      // if conditionals were not updated, copy back pointer to conditional array
      conditionalPointer = locusData->nodeArray[nodeId]->conditionalProbs;
      locusData->nodeArray[nodeId]->conditionalProbs = locusData->savedVersion.savedNodes[nodeId]->conditionalProbs;
      locusData->savedVersion.savedNodes[nodeId]->conditionalProbs = conditionalPointer;
    }
  }
	
  // for nodes which have not been changed, but whose conditional probabilities have been recomputed
  // switch pointers to conditional array
  for(i=0; i<locusData->savedVersion.numChangedConditionals; i++) {
    nodeId = locusData->savedVersion.changedCondIds[i];
    if(locusData->savedVersion.recalcConditionals[nodeId]) {
      // if conditionals were not already copied in previous loop, copy them now
      conditionalPointer = locusData->nodeArray[nodeId]->conditionalProbs;
      locusData->nodeArray[nodeId]->conditionalProbs = locusData->savedVersion.savedNodes[nodeId]->conditionalProbs;
      locusData->savedVersion.savedNodes[nodeId]->conditionalProbs = conditionalPointer;
      locusData->savedVersion.recalcConditionals[nodeId] = 0;
    }			
  }

	
  locusData->savedVersion.numChangedNodes = 0;
  locusData->savedVersion.numChangedConditionals = 0;


  return 0;
}
/** end of revertToSaved **/



/***********************************************************************************
 *	resetSaved
 *	- resets saved data and maintains updates
 *	- locus is number of locus
 *	- returns 0
 ***********************************************************************************/
int resetSaved(LocusData* locusData) {
	
  locusData->savedVersion.copyAll = 0;
  locusData->savedVersion.numChangedNodes = 0;
  locusData->savedVersion.numChangedConditionals = 0;
  locusData->savedVersion.root = -1;
  locusData->savedVersion.dataLogLikelihood = locusData->dataLogLikelihood;
	
  // copy zeros into recalcConditionals[] array for internal nodes
  resetBooleanArray(locusData->savedVersion.recalcConditionals, 2*locusData->numLeaves-1);
	
  return 0;
}
/** end of resetSaved **/



/***********************************************************************************
 *	adjustGenNodeAge
 *	- adjusts nodes age and saves old version
 *	- nodeId is id of node and age is new age
 *	- returns 0
 ***********************************************************************************/
int	adjustGenNodeAge(LocusData* locusData, int nodeId, double age)		{
  copyNodeToSaved(locusData, nodeId, 1);		// conditionals should be recomputed for such node
  
  locusData->nodeArray[nodeId]->age = age;
  
  return 0;
}
/** end of adjustGenNodeAge **/



/***********************************************************************************
 *	scaleAllNodeAges
 *	- scales all node ages in genealogy with a given multiplicative factor
 *	- returns the delta in log-likelihood of suggested step
 *	- saves all original ages and conditionals
 ***********************************************************************************/
/* MARK: NOTE THAT RESCALING OF AGES IS DONE ALSO FOR LEAVES ASSOCIATED WITH
         ANCIENT SAMPLES
*/
double	scaleAllNodeAges(LocusData* locusData, double factor){
  int nodeId;
  double oldLnLd = locusData->dataLogLikelihood;
	
  // mark all nodes as copied
  locusData->savedVersion.copyAll = 1;
	

  // save all fathers of leaves (for quick switch back to saved)
  for(nodeId = 0; nodeId < locusData->numLeaves; nodeId++) {
    locusData->savedVersion.savedNodes[nodeId]->father = locusData->nodeArray[nodeId]->father;
  }
	
  // traverse all internal nodes and change age by factor
//  for(nodeId = locusData->numLeaves; nodeId<2*locusData->numLeaves-1; nodeId++) {
  for(nodeId = 0; nodeId<2*locusData->numLeaves-1; nodeId++) {
    adjustGenNodeAge(locusData, nodeId, factor*locusData->nodeArray[nodeId]->age);
  }
	
  computeLocusDataLikelihood(locusData, /*reuse conditionals*/ 1);
	
  return locusData->dataLogLikelihood - oldLnLd;
}
/** end of scaleAllNodeAges **/



/***********************************************************************************
 *	executeGenSPR
 *	- executes an SPR operation on locus genealogy
 *	- subtree_root is node id for root of pruned subtree,targetBranch is id of node below
 *		branch where subtree should be regrafted, and age indicates time of regrafting.
 *	- returns 0 if root node remains the same
 *	- otherwise, returns 1 if subtree is a child of the root AFTER regrafting
 *	- otherwise, returns 2 if subtree was a child of the root BEFORE regrafting
 ***********************************************************************************/
int	executeGenSPR(LocusData* locusData, int subtreeRoot, int targetBranch, double age)		{

  int father, grandpa, sibling, targetFather;

  targetFather = locusData->nodeArray[targetBranch]->father;
  father  = locusData->nodeArray[subtreeRoot]->father;
  grandpa = locusData->nodeArray[father]->father;
  sibling = locusData->nodeArray[father]->leftSon + locusData->nodeArray[father]->rightSon - subtreeRoot;
	
  //	printf("\nProposing SPR of subtree rooted at node %d, father %d (age %f), sibling %d, target branch %d, new age %f.\n",
  //				subtreeRoot, father, locusData->nodeArray[father]->age, sibling, targetBranch, age);
	
  adjustGenNodeAge(locusData, father, age);
	
  // if target branch is above sibling or father nodes, no topological change is required
  if(targetBranch == sibling || targetBranch == father)
    return 0;

  // prune subtree

  // mark sibling as changed (due to change in father pointer) - do not save conditionals
  copyNodeToSaved(locusData, sibling, 0);
  locusData->nodeArray[sibling]->father = grandpa;
  if(grandpa >= 0) {
    // mark grandpa as changed - conditionals should be changed
    // root is replaced for sibling later, if father was root
    copyNodeToSaved(locusData, grandpa, 1);
    if(locusData->nodeArray[grandpa]->leftSon == father) {
      locusData->nodeArray[grandpa]->leftSon  = sibling;
    } else {
      locusData->nodeArray[grandpa]->rightSon = sibling;			
    }
  }
	
  // regraft subtree
  locusData->nodeArray[father]->father	= targetFather;
  locusData->nodeArray[father]->leftSon 	= subtreeRoot;
  locusData->nodeArray[father]->rightSon	= targetBranch;

  /********************************************************************
   * we have to make sure not to save a node twice (in two "identities")
   * possibilities for duplicate identities include:
   * - target       == grandpa
   * - targetFather == grandpa
   * - targetFather == sibling
   ********************************************************************/
	
  // mark targetBranch as changed (due to change in father pointer) - do not save conditionals
  if(targetBranch != grandpa) {	
    copyNodeToSaved(locusData, targetBranch, 0);
  }
  locusData->nodeArray[targetBranch]->father = father;

  // if regrafting above root, set new root to 'father' and return 1
  if(targetFather < 0) {
    locusData->savedVersion.root = targetBranch;
    locusData->root = father;
    return 1;
  }

  // mark target's father as changed - conditionals should be changed
  if(targetFather == sibling) {
    copyNodeConditionals(locusData, targetFather);
  } else if(targetFather != grandpa) {
    copyNodeToSaved(locusData, targetFather, 1);
  }
	
  if(locusData->nodeArray[targetFather]->leftSon == targetBranch) {
    locusData->nodeArray[targetFather]->leftSon  = father;
  } else {
    locusData->nodeArray[targetFather]->rightSon = father;			
  }

  // if father was root, set new root to 'sibling' and return 2
  if(grandpa < 0) {
    locusData->savedVersion.root = father;
    locusData->root = sibling;
    return 2;		
  }
	
  return 0;	
}
/** end of executeGenSPR **/



/***********************************************************************************
 *	copyGenericTreeToLocus
 * 	- copies a generic tree into likelihood tree
 *	- assumes label1 holds age of node
 *	- returns 0
 ***********************************************************************************/
int copyGenericTreeToLocus(LocusData* locusData, GenericBinaryTree* genericTree) {
  int node, numNodes = 2*locusData->numLeaves - 1;

  locusData->root = genericTree->rootId;
  for(node=0; node<numNodes; node++) {
    locusData->nodeArray[node]->father   = genericTree->father[node];
    locusData->nodeArray[node]->leftSon  = genericTree->leftSon[node];
    locusData->nodeArray[node]->rightSon = genericTree->rightSon[node];
    locusData->nodeArray[node]->age      = genericTree->label1[node];
  }

  return 0;	
}
/** end of copyGenericTreeToLocus **/




/***********************************************************************************
 *	printLocusGenTree
 * 	- prints the genealogy tree (including supplied population and event id)
 *	- prints each node in a separate line, according to id order
 *	- if there are nodes which are considered for changes, print out their id's
 ***********************************************************************************/
void printLocusGenTree(LocusData* locusData, FILE* stream, int* nodePops, int* nodeEvents)	{
	
  int node, numNodes = 2*locusData->numLeaves - 1;

  fprintf(stream, "Genalogy tree:\n");
  for(node=0; node<numNodes; node++) {
    fprintf(stream, "Node %2d, age [%.10f], father (%2d), sons (%2d %2d), pop (%2d), event-id (%2d)",
           node, locusData->nodeArray[node]->age, locusData->nodeArray[node]->father,
           locusData->nodeArray[node]->leftSon, locusData->nodeArray[node]->rightSon,
           nodePops[node], nodeEvents[node]);
    if(locusData->root == node)			fprintf(stream, " - Root\n");
    else if(node<locusData->numLeaves)	fprintf(stream, " - Leaf\n");
    else								fprintf(stream, "\n");
  }
	
  fprintf(stream, "---------------------------------------------------------------\n");
  if(locusData->savedVersion.numChangedNodes > 0) {
    fprintf(stream, "There are %d changed nodes:",locusData->savedVersion.numChangedNodes);
    for(node=0; node<locusData->savedVersion.numChangedNodes; node++) {
      fprintf(stream, " %d",locusData->savedVersion.changedNodeIds[node]);
    }
    fprintf(stream, "\n---------------------------------------------------------------\n");
  }
			
  return;
}
/** end of printLocusGenTree **/




/***********************************************************************************
 *	printLocusDataStats
 * 	- prints stats on alignment patterns
 *	- stats are outputted in one line (with newline) to stdout in the following order:
 *	- num total patterns, num hom patterns, num het patterns with 2 phases, , num het patterns with 4 phases...
 *	- for each of the above stats, prints consecutively the number of distinct patterns 
 *		and the number of columns corresponding to it in the alignment
 *	- maxStat is an upper bound on the stats collected
 ***********************************************************************************/
void printLocusDataStats(LocusData* locusData, int maxStat)	{
  int patt, numHetsPerPatt, numPhases;
  int *numColArray, *numPattArray = (int*)malloc(2*(maxStat+1)*sizeof(int));
	
  if(numPattArray == NULL) {
    fprintf(stderr, "Error: Out Of Memory het stats array in printLocusDataStats().\n");
    exit(-1);
  }
  numColArray = numPattArray + maxStat + 1;
	
  // initialize counts
  for(numHetsPerPatt=0; numHetsPerPatt<=maxStat; numHetsPerPatt++) {
    numColArray [numHetsPerPatt] = 0;
    numPattArray[numHetsPerPatt] = 0;
  }
	
  // compute stats for all hets
  for(patt=0; patt<locusData->seqData.numPatterns; patt+= locusData->seqData.numPhases[patt]) {
    // set numHetsPerPatt = log_2( locusData->seqData.numPhases[patt] )
    numPhases = locusData->seqData.numPhases[patt];
    numHetsPerPatt = 0;
    while( numPhases > 1) {
      numPhases /= 2;
      numHetsPerPatt++;
    }
	
    numPattArray[numHetsPerPatt]++;
    numColArray [numHetsPerPatt] += locusData->seqData.patternCount[patt];
  }
	
  printf("%d",locusData->seqData.numPatterns);
  for(numHetsPerPatt=0; numHetsPerPatt<= maxStat; numHetsPerPatt++) {
    printf("\t%d\t%d",numPattArray[numHetsPerPatt],numColArray [numHetsPerPatt]);
  }
	
  printf("\n");
	
  free(numPattArray);
  return;
}



/***********************************************************************************
 *	printLocusDataPatterns
 * 	- prints the alignment of the locus, pattern by pattern to output file.
 *	- patterns are printed in columns, according to leaf id.
 *	- below each pattern is its multiplicity in the alignment
 *	- het patterns are represented by all phasings 
 *		(mult is written below first phasing)
 ***********************************************************************************/
void printLocusDataPatterns(LocusData* locusData, FILE* outFile)	{
  static const char baseSymbols[] = "TCAGYKWSMRN";

  int leaf, patt, base, ambigSize, firstBase, secondBase, phase;
  double sumConds;
  double* leafConditionals;
  char ch;
	
	
  fprintf(outFile,"\n%d phased patterns:",locusData->seqData.numPatterns);
	
  for(leaf=0; leaf<locusData->numLeaves; leaf++) {
    leafConditionals = locusData->nodeArray[leaf]->conditionalProbs;
    fprintf(outFile,"\n%5d",leaf+1);
    for(patt=0; patt<locusData->seqData.numPatterns; patt++, leafConditionals += CODE_SIZE) {
      ambigSize = 0;
      firstBase = secondBase = -1;
      sumConds = 0.0;
			
      for(base=0; base<CODE_SIZE; base++) {
        sumConds += leafConditionals[base];
        if(leafConditionals[base] > 0.0) {
          ambigSize++;
          if(firstBase < 0)		firstBase = base;
          else if(secondBase < 0)	secondBase = base;
        }
      }
      if(ambigSize<1 || ambigSize>4 || firstBase<0 || firstBase>3 || secondBase>3 || (ambigSize>1 && secondBase<=firstBase) || (sumConds != 1.0 && sumConds != 4.0)) {
        fprintf(stderr, "\nError: Fatal Error in parsing leaf conditionals for leaf %d, pattern %d.\n",leaf,patt);
      }
			
			
      if(ambigSize == 1) {
        ch = baseSymbols[firstBase];
      } else if(ambigSize == 2) {
        if(firstBase == 0) {
          ch = baseSymbols[3+secondBase];
        } else if(firstBase == 1) {
          ch = baseSymbols[5+secondBase];
        } else {
          ch = 'R';
        }
      } else if(ambigSize == 4) {
        ch = 'N';
      } else {
        ch = 'X';
      }				
      fprintf(outFile,"%5c", ch);
    }//end of for(patt)
  }// end of for(leaf)
	
  fprintf(outFile,"\ncount");
  for(patt=0; patt<locusData->seqData.numPatterns; patt+=locusData->seqData.numPhases[patt]) {
    fprintf(outFile,"%5d",locusData->seqData.patternCount[patt]);
    for(phase=1; phase < locusData->seqData.numPhases[patt]; phase++) {
      fprintf(outFile,"     ");
    }
  }
  if(patt != locusData->seqData.numPatterns) {
    fprintf(stderr, "\nError: in locus data total number of phased patterns is %d, but recorded to be %d.\n", patt, locusData->seqData.numPatterns);
  }
	
  fprintf(outFile,"\n");
		
					
  return;
}
/** end of printLocusDataPatterns **/



/***********************************************************************************
 *	computePairwiseLCAs
 *	- procedure for computing a 2D matrix with the ids of the LCAs (least
 *    common ancestors) of all pairs of leaves
 *  - calls recursive procedure computePairwiseLCAs_rec on root
 *  - returns 1 if successful, 0 otherwise
 ***********************************************************************************/
int computePairwiseLCAs (LocusData* locusData, int** lcaMatrix, int* leafArray_aux){
  int numLeaves;
  int res = computePairwiseLCAs_rec (locusData, locusData->root, lcaMatrix, leafArray_aux, 0, &numLeaves);
  
  if(!res || numLeaves != locusData->numLeaves) {
    return 0;
  }
  return 1;
}
/** end of computePairwiseLCAs **/



/***********************************************************************************
 *	getSortedAges
 *	- procedure for computing a sorted list of node ages (from most recent to root)
 *  - assumes array given as input has 2x space for all internal nodes (also for auxiliary space for sorting)
 *  - calls recursive procedure getSortedAges_rec on root
 *  - returns 1 if successful, 0 otherwise
 ***********************************************************************************/
int getSortedAges (LocusData* locusData, double* ageArray){
  int numInternalNodes;
  int res = getSortedAges_rec(locusData, locusData->root, ageArray, ageArray+locusData->numLeaves-1, 0, &numInternalNodes);
	
  if(!res || numInternalNodes != locusData->numLeaves-1) {
    return 0;
  }
  return 1;
}
/** end of getSortedAges **/


/***********************************************************************************
 *	getLocusDataLikelihood
 *	- returns the log likelihood of the locus as last recorded (no new computations are performed)
 ***********************************************************************************/
double getLocusDataLikelihood (LocusData* locusData){
  return locusData->dataLogLikelihood;
}
/** end of getLocusDataLikelihood **/



/***********************************************************************************
 *	getLocusRoot
 *	- returns the root node id
 ***********************************************************************************/
int getLocusRoot (LocusData* locusData)	{
  return locusData->root;
}
/** end of getLocusRoot **/



/***********************************************************************************
 *	getNodeAge
 *	- returns the age of a node
 ***********************************************************************************/
double getNodeAge (LocusData* locusData, int nodeId)	{
  return locusData->nodeArray[nodeId]->age;
}
/** end of getNodeAge **/



/***********************************************************************************
 *	getNodeFather
 *	- returns the id of the father of a node
 ***********************************************************************************/
int getNodeFather (LocusData* locusData, int nodeId)	{
  return locusData->nodeArray[nodeId]->father;
}
/** end of getNodeFather **/



/***********************************************************************************
 *	getNodeSon
 *	- returns the id of a son of a node (son = 0 -> left, son = 1 -> right)
 ***********************************************************************************/
int getNodeSon (LocusData* locusData, int nodeId, unsigned short son)	{
  if(son)		return locusData->nodeArray[nodeId]->rightSon;
  else		return locusData->nodeArray[nodeId]->leftSon;
}
/** end of getNodeSon **/



/***************************************************************************************************************/
/******                              INTERNAL FUNCTION IMPLEMENTATION                                     ******/
/***************************************************************************************************************/



/***********************************************************************************
 *	computeLeafConditionals
 *	- computes conditional probabilities for leaves given a sequence pattern in a char array
 *	- advances locusData->numPattern accrodingly
 *	- patterns should be entered one by one
 * 	- note that this is done before sampling begins (when data structure is initialized)
 * 		and these conditionals remain untouched after this point.
 * 	- we allow only nucleotide, het-ambiguity symbols or 'N' characters in a pattern (11 possibilities total)
 *	- saves same values in "saved conditionals array" - these do not change
 *	- returns 0, if all OK and -1, if bad pattern
 ***********************************************************************************/
int computeLeafConditionals(LocusData* locusData, char* patternString)	{
  int leaf, base;
  double *conditionals, *conditionalsForSaved;

  for(leaf=0; leaf<locusData->numLeaves; leaf++) {
    conditionals = locusData->nodeArray[leaf]->conditionalProbs + locusData->seqData.numPatterns * CODE_SIZE;
    conditionalsForSaved = locusData->savedVersion.savedNodes[leaf]->conditionalProbs + locusData->seqData.numPatterns * CODE_SIZE;
    for(base=0; base<CODE_SIZE; base++) {
      conditionals[base] = 0.0;
    }
/*    if(patternString[leaf] != 'T' && patternString[leaf] != 'C') {
		printf("found character %c.\n",patternString[leaf]);
		patternString[leaf] = 'N';
	}
*/
  switch(patternString[leaf]) {
    case('T'):
      conditionals[0] = 1.0;
      break;
    case('C'):
      conditionals[1] = 1.0;
      break;
    case('A'):
//		printf("AAAAA - Error!!\n");
      conditionals[2] = 1.0;
      break;
    case('G'):
//		printf("GGGGG - Error!!\n");
      conditionals[3] = 1.0;
      break;
      /*			case('Y'):
                    conditionals[0] = 0.5;
                    conditionals[1] = 0.5;
                    break;
                    case('K'):
                    conditionals[0] = 0.5;
                    conditionals[2] = 0.5;
                    break;
                    case('W'):
                    conditionals[0] = 0.5;
                    conditionals[3] = 0.5;
                    break;
                    case('S'):
                    conditionals[1] = 0.5;
                    conditionals[2] = 0.5;
                    break;
                    case('M'):
                    conditionals[1] = 0.5;
                    conditionals[3] = 0.5;
                    break;
                    case('R'):
                    conditionals[2] = 0.5;
                    conditionals[3] = 0.5;
                    break;
      */
    case('N'):
//		printf("NNNNN - Error!!\n");
      conditionals[0] = 1.0;
      conditionals[1] = 1.0;
      conditionals[2] = 1.0;
      conditionals[3] = 1.0;
      break;
    default:
      fprintf(stderr, "\nError: Unexpected character '%c' for leaf %d in pattern.\n",leaf, patternString[leaf]);
      return -1;
    }// end of switch
    // copy conditionals to saved
    for(base=0; base<CODE_SIZE; base++) {
      conditionalsForSaved[base] = conditionals[base];
    }
		
  }// end of for(leaf)


  // advance number of patterns in locus
  locusData->seqData.numPatterns++;
  return 0;
}
/** end of computeLeafConditionals **/



/***********************************************************************************
 *	computeConditionalJC
 *	- RECURSIVE PROCEDURE
 *	- computes conditional probabilities for a subtree of the genealogy at a specified locus
 *	- nodeId indicates the root of the subtree
 *	- conditional probabilities are written in conditionalProbs[] array of node
 *	- returns 1 if conditionals had to be recomputed, and 0 if old ones were used
 *	- recomputations are needed if this node has been modified or if recomputations were
 *		made in one of its subtrees
 *	- if overideOld == 1, then does not save old version
 ***********************************************************************************/
int computeConditionalJC (LocusData* locusData, int nodeId, int numPatterns, int* patternIds, unsigned short overideOld)		{
  int res;
  int patt, pattId, base;
  double edgeLength;
  LikelihoodNode *node, *leftSon, *rightSon;
  double leftEdgeConditionalProb[2];
  double rightEdgeConditionalProb[2];

	
  // leaf case - no need for any computation
  if(nodeId < locusData->numLeaves)
    return 0;
	
  node = locusData->nodeArray[nodeId];
  // internal node - compute conditionals for both children
  res = computeConditionalJC(locusData, node->leftSon, numPatterns, patternIds, overideOld);
  res = computeConditionalJC(locusData, node->rightSon, numPatterns, patternIds, overideOld) || res;
	
  //	printf("Computing conditionals for node %d:\n",nodeId);

  if(!overideOld && !res && !locusData->savedVersion.recalcConditionals[nodeId])
    return 0;
	
  // save old conditional probabilities (if haven't already been saved)
  if(!overideOld) {
    copyNodeConditionals(locusData,nodeId);
  }

	
  leftSon = locusData->nodeArray[ node->leftSon ];
  rightSon = locusData->nodeArray[ node->rightSon ];

  edgeLength = locusData->mutationRate * (node->age - leftSon->age);
  leftEdgeConditionalProb[1] = computeEdgeConditionalJC(edgeLength);
  leftEdgeConditionalProb[0] = 1 - 3.0*leftEdgeConditionalProb[1];

  edgeLength = locusData->mutationRate * (node->age - rightSon->age);
  rightEdgeConditionalProb[1] = computeEdgeConditionalJC(edgeLength);
  rightEdgeConditionalProb[0] = 1 - 3.0*rightEdgeConditionalProb[1];
	

  for (patt=0; patt < numPatterns; patt++) {
    pattId = patternIds[patt];
    // initialize conditionals
    for(base=0; base<CODE_SIZE; base++)  {
      node->conditionalProbs[CODE_SIZE*pattId + base] = 1.0;
    }
    //		printf("edge (%d,%d)", nodeId,node->leftSon);
    computeSubtreeConditionals(&(leftSon->conditionalProbs[CODE_SIZE*pattId]),&(node->conditionalProbs[CODE_SIZE*pattId]),leftEdgeConditionalProb);
    //		printf(", edge (%d,%d)", nodeId,node->rightSon);
    computeSubtreeConditionals(&(rightSon->conditionalProbs[CODE_SIZE*pattId]),&(node->conditionalProbs[CODE_SIZE*pattId]),rightEdgeConditionalProb);
    //		printf(".\n");
  }
               
  return 1;
}
/** end of computeConditionalJC **/



/***********************************************************************************
 *	computeSubtreeConditionals
 *	- computes conditional probabilities for a subtree rooted at some edge
 *	- assumes conditionals at bottom of edge are given in sonConditionals CODE_SIZE-long array
 *	- uses son conditionals to compute parentConditionals (CODE_SIZE-long array)
 *	- does this through the use of 2-long array edgeConditionals
 *	- multiplies the values in parentConditionals with contribution from son
 ***********************************************************************************/
void computeSubtreeConditionals (double* sonConditionals, double* parentConditionals, double* edgeConditionals)		{
  int sonState, sonBase, fatherBase;
  double prob;

  for(sonBase=0; sonBase<CODE_SIZE; sonBase++)  {
    //		printf(" S%d=%3lf",sonBase,sonConditionals[sonBase]);
  }
  // first determine the state of the son (if a nucleotide leaf, or missing data)
  // state -1 means non-base and non-N states are observed
  // state i=0..3 means base i has prob=1, and all others have prob=0
  // state 4 means that all previous bases have prob=1 (at least 2)
  sonState = -1;
  for(sonBase=0; sonBase<CODE_SIZE; sonBase++)  {
    prob = sonConditionals[sonBase];
    if(prob == 0.0) {
      if(sonState < 4) {
        continue;
      } else {
        sonState = -1;
        break;
      }
    } else if(prob == 1.0) {
      if(sonState == -1) {
        sonState = sonBase;
      } else if(sonState == 4) {
        continue;
      } else if(sonState == 0 && sonBase == 1) {  
        sonState = 4;
      } else {
        sonState = -1;
        break;
      }
    } else {
      sonState = -1;
      break;
    }
  }// end of for(sonBase)
		
  if(sonState == 4) {	// son is missing data
    return;
    //		printf(" son is 'N',");
  } else if(sonState >= 0) { // son is specific base
    //		printf(" son is base %d,",sonState+1);
    for(fatherBase=0; fatherBase<CODE_SIZE; fatherBase++)  {
      parentConditionals[fatherBase] *= edgeConditionals[fatherBase != sonState];
      //			printf(" %3lf",parentConditionals[fatherBase]);
    }
  } else {
    for(fatherBase=0; fatherBase<CODE_SIZE; fatherBase++)  {
      prob = 0.0;
      for(sonBase=0; sonBase<CODE_SIZE; sonBase++)  {
        prob += edgeConditionals[fatherBase != sonBase] * sonConditionals[sonBase];
      }
      parentConditionals[fatherBase] *= prob;
      //			printf(" %3lf",parentConditionals[fatherBase]);
    }
  }
		
  return;
}
/** end of computeSubtreeConditionals **/



/***********************************************************************************
 *	computeConditionalJC_new
 *	- RECURSIVE PROCEDURE
 *	-> SAME AS ORIGINAL LOGIC BUT CALLS computeSubtreeConditionals_new
 *	- computes conditional probabilities for a subtree of the genealogy at a specified locus
 *	- nodeId indicates the root of the subtree
 *	- conditional probabilities are written in conditionalProbs[] array of node
 *	- returns 1 if conditionals had to be recomputed, and 0 if old ones were used
 *	- recomputations are needed if this node has been modified or if recomputations were
 *		made in one of its subtrees
 *	- if overideOld == 1, then does not save old version
 ***********************************************************************************/
int computeConditionalJC_new (LocusData* locusData, int nodeId, int numPatterns, int* patternIds, unsigned short overideOld)		{
  int res;
  int patt, pattId, base;
  double edgeLength;
  LikelihoodNode *node, *leftSon, *rightSon;
  double leftEdgeConditionalProb[2];
  double rightEdgeConditionalProb[2];

	
  // leaf case - no need for any computation, but propagate up, if leaf age has changed
  if(nodeId < locusData->numLeaves) {
    if(locusData->savedVersion.recalcConditionals[nodeId]) {
      //printf("leaf %d updated. Will have to recompute parent.\n",nodeId);
      return 100;
    }
    return locusData->savedVersion.recalcConditionals[nodeId];
  }
  node = locusData->nodeArray[nodeId];
  // internal node - compute conditionals for both children
  res = computeConditionalJC_new(locusData, node->leftSon, numPatterns, patternIds, overideOld);
  res = computeConditionalJC_new(locusData, node->rightSon, numPatterns, patternIds, overideOld) + res;

  //	printf("Computing conditionals for node %d:\n",nodeId);

  if(!overideOld && !res && !locusData->savedVersion.recalcConditionals[nodeId])
    return 0;
	
  // save old conditional probabilities (if haven't already been saved)
  if(!overideOld) {
    copyNodeConditionals(locusData,nodeId);
  }

	
  leftSon = locusData->nodeArray[ node->leftSon ];
  rightSon = locusData->nodeArray[ node->rightSon ];


  edgeLength = locusData->mutationRate * (node->age - leftSon->age);
  leftEdgeConditionalProb[0] = computeEdgeConditionalJC(edgeLength);
  leftEdgeConditionalProb[1] = 1 - 4.0*leftEdgeConditionalProb[0];

  edgeLength = locusData->mutationRate * (node->age - rightSon->age);
  rightEdgeConditionalProb[0] = computeEdgeConditionalJC(edgeLength);
  rightEdgeConditionalProb[1] = 1 - 4.0*rightEdgeConditionalProb[0];
	
  if(res>10) {
    //printf("Node %d parent of %d,%d, one of which changed. Ages %g, %g, %g. Conditional probs: %g, %g.\n",
    //         nodeId, node->leftSon, node->rightSon,node->age, leftSon->age, rightSon->age,leftEdgeConditionalProb[0],rightEdgeConditionalProb[0]);
  }

  for (patt=0; patt < numPatterns; patt++) {
    pattId = patternIds[patt];
#ifdef OPT2
    if(nodeId != locusData->root && locusData->seqData.numBases[pattId] == 1) {
	  node->conditionalProbs[CODE_SIZE*pattId] = 
	    leftSon->conditionalProbs[CODE_SIZE*pattId]*(leftEdgeConditionalProb[0]+leftEdgeConditionalProb[1])*
	    rightSon->conditionalProbs[CODE_SIZE*pattId]*(rightEdgeConditionalProb[0]+rightEdgeConditionalProb[1]);
		
      for(base=1; base<CODE_SIZE; base++)  {
        node->conditionalProbs[CODE_SIZE*pattId + base] = 0.0;
      }
      continue;
	}
#endif
    // initialize conditionals
    for(base=0; base<CODE_SIZE; base++)  {
      node->conditionalProbs[CODE_SIZE*pattId + base] = 1.0;
    }
    //		printf("edge (%d,%d)", nodeId,node->leftSon);
    computeSubtreeConditionals_new(&(leftSon->conditionalProbs[CODE_SIZE*pattId]),&(node->conditionalProbs[CODE_SIZE*pattId]),leftEdgeConditionalProb);
    //		printf(", edge (%d,%d)", nodeId,node->rightSon);
    computeSubtreeConditionals_new(&(rightSon->conditionalProbs[CODE_SIZE*pattId]),&(node->conditionalProbs[CODE_SIZE*pattId]),rightEdgeConditionalProb);
    //		printf(".\n");
  }
               
  return 1;
}
/** end of computeConditionalJC_new **/



/***********************************************************************************
 *	computeSubtreeConditionals_new
 *	-> SAME AS ORIGINAL LOGIC BUT WITHOUT THE DOUBLE LOOP
 *	- computes conditional probabilities for a subtree rooted at some edge
 *	- assumes conditionals at bottom of edge are given in sonConditionals CODE_SIZE-long array
 *	- uses son conditionals to compute parentConditionals (CODE_SIZE-long array)
 *	- does this through the use of 2-long array edgeSubstProb of p=edge probability of (non-identity) transition of  and 1-4p
 *	- multiplies the values in parentConditionals with contribution from son
 ***********************************************************************************/
void computeSubtreeConditionals_new (double* sonConditionals, double* parentConditionals, double* edgeSubstProb)		{
  int base;
  double probSum, probSumTimesSubst;

  probSum = 0.0;
  for(base=0; base<CODE_SIZE; base++)  {
//		printf(" S%d=%3lf",base,sonConditionals[base]);
	probSum += sonConditionals[base];
  }// end of for(sonBase)
  
  if(probSum >= CODE_SIZE) {
//	printf(" son is 'N',");
	  return;
  }

  probSumTimesSubst = probSum * edgeSubstProb[0];
  
  for(base=0; base<CODE_SIZE; base++)  {
      parentConditionals[base] *= (probSumTimesSubst + sonConditionals[base]*edgeSubstProb[1]);
      //			printf(" %3lf",parentConditionals[fatherBase]);
  }
		
  return;
}
/** end of computeSubtreeConditionals_new **/


/***********************************************************************************
 *	computePairwiseLCAs_rec
 *	- recursive procedure for computing a 2D matrix with the ids of the LCAs (least
 *    common ancestors) of all pairs of leaves
 *  - starts with a designated internal node (nodeID)
 *  - recursively computes the set of right / left leaves under that node (leafArray, arrayOffset,  numLeaves_out)
 *  - returns 1 if successful, 0 otherwise
 ***********************************************************************************/
int computePairwiseLCAs_rec (LocusData* locusData, int nodeId, int** lcaMatrix, int* leafArray, int arrayOffset, int* numLeaves_out){
  LikelihoodNode *node;
  int res, numLeftLeaves, numRightLeaves, l, r;
	
  if(nodeId < 0) {
	return 0;
  }

  // leaf case - just put node id in leaf array
  if(nodeId < locusData->numLeaves) {
	leafArray[arrayOffset] = nodeId;
	*numLeaves_out   = 1;
    return 1;
  }
  // internal node
	
  node = locusData->nodeArray[nodeId];
  
  // recursively call for left and right sons and aggregate children list
  numLeftLeaves = numRightLeaves = 0;
  res = computePairwiseLCAs_rec (locusData, node->leftSon,  lcaMatrix, leafArray, arrayOffset, &numLeftLeaves);
  if(!res || arrayOffset+numLeftLeaves >= locusData->numLeaves) {
	return 0;
  }
  res = computePairwiseLCAs_rec (locusData, node->rightSon, lcaMatrix, leafArray, arrayOffset + numLeftLeaves, &numRightLeaves);
  if(!res) {
	return 0;
  }

  // fill in matrix elements with node id
  for(l=arrayOffset; l<arrayOffset+numLeftLeaves; l++) {
	for(r=arrayOffset+numLeftLeaves; r<arrayOffset+numLeftLeaves+numRightLeaves; r++) {
	  lcaMatrix [leafArray[l]] [leafArray[r]] = nodeId;
	  lcaMatrix [leafArray[r]] [leafArray[l]] = nodeId;
	}
  }
  
  *numLeaves_out = numLeftLeaves + numRightLeaves;
  return 1;
}
/** end of computePairwiseLCAs_rec **/


/***********************************************************************************
 *	getSortedAges_rec
 *	- recursive procedure for computing sorted list of node ages
 *  - starts with a designated internal node (nodeID)
 *  - recursively computes the set of right / left leaves under that node (leafArray, arrayOffset,  numLeaves_out)
 *  - assumes array is doubled for auxiliary space
 *  - returns 1 if successful, 0 otherwise
 ***********************************************************************************/
int getSortedAges_rec (LocusData* locusData, int nodeId, double* sortedAges, double* sortedAges_aux, int arrayOffset, int* numInternalNodes_out){
  LikelihoodNode *node;
  int res, numLeftInternalNodes, numRightInternalNodes, l, r, i;
	
  if(nodeId < 0) {
	return 0;
  }

  // leaf case - just put node id in leaf array
  if(nodeId < locusData->numLeaves) {
	*numInternalNodes_out   = 0;
    return 1;
  }
  // internal node  -- need space to add entry
  if(arrayOffset >= locusData->numLeaves-1) {
	return 0;
  }
	
  node = locusData->nodeArray[nodeId];
  
  // recursively call for left and right sons and aggregate children list
  numLeftInternalNodes = numRightInternalNodes = 0;
  res = getSortedAges_rec(locusData, node->leftSon,  sortedAges, sortedAges_aux, arrayOffset, &numLeftInternalNodes);
  if(!res) {
	return 0;
  }

  res = getSortedAges_rec(locusData, node->rightSon, sortedAges, sortedAges_aux, arrayOffset + numLeftInternalNodes, &numRightInternalNodes);
  if(!res) {
	return 0;
  }

  sortedAges[arrayOffset + numLeftInternalNodes + numRightInternalNodes] = locusData->nodeArray[nodeId]->age;
  *numInternalNodes_out = numLeftInternalNodes + numRightInternalNodes + 1;

  if (numLeftInternalNodes == 0 || numRightInternalNodes == 0){
	  return 1;
  }

  // move to auxilliary array
  for(i=0; i<=numLeftInternalNodes+numRightInternalNodes; i++) {
	  sortedAges_aux[i] = sortedAges[arrayOffset+i];
  }
  // merge sorted lists
  l = 0;
  r = numLeftInternalNodes;
  i = arrayOffset;
  while(i<=arrayOffset+numLeftInternalNodes+numRightInternalNodes && l<numLeftInternalNodes && r<=numLeftInternalNodes+numRightInternalNodes) {
	  if(sortedAges_aux[l] < sortedAges_aux[r]) {
		  sortedAges[i] = sortedAges_aux[l];
		  i++;
		  l++;
	  } else {
		  sortedAges[i] = sortedAges_aux[r];
		  i++;
		  r++;
	  }
  }
  if(l>=numLeftInternalNodes) {
	  while(i<=arrayOffset+numLeftInternalNodes+numRightInternalNodes && r<=numLeftInternalNodes+numRightInternalNodes) {
		  sortedAges[i] = sortedAges_aux[r];
		  i++;
		  r++;
	  }
  } else {
	  while(i<=arrayOffset+numLeftInternalNodes+numRightInternalNodes && l<numLeftInternalNodes) {
		  sortedAges[i] = sortedAges_aux[l];
		  i++;
		  l++;
	  }
  }
  if(i != arrayOffset+numLeftInternalNodes+numRightInternalNodes+1 || l != numLeftInternalNodes || r!=numLeftInternalNodes+numRightInternalNodes+1) {
	  return 0;
  }
  return 1;
}
/** end of getSortedAges_rec **/


/***********************************************************************************
 *	computeEdgeConditionalJC
 *	- computes probability of seeing a specific transition along a given egde
 *	- edgeLength is the length of the edge
 ***********************************************************************************/
double computeEdgeConditionalJC (double edgeLength)		{
	
	
//	return global_prob_val;
//	edgeLength = 0.0000001;
	
  if (edgeLength <  -0.0001) {
	  if(debug) {
    	fprintf(stderr, "\nWarning: Negative edge length %g in computeEdgeConditionalsJC", edgeLength);
	  }
  }
 	
  if (edgeLength < 1e-100) {
    return 0.0;
  }
	
  return ((1-exp(-4*edgeLength/3.0)) / 4.0 );  
}
/** end of computeEdgeConditionalJC **/



/***********************************************************************************
 *	copyNodeToSaved
 *	- copies data from current version of node to saved version
 *	- marks node as changed
 *	- nodeId is id of node being switched
 *	- if recalcConditionals is 1, then calls copyNodeConditionals to switch pointers
 *		to conditional probability array
 *	- IT IS IMPORTANT THAT THIS PROCEDURE IS CALLED NO MORE THAN ONCE FOR EVERY NODE
 *		IN EACH SAMPLING ITERATION 
 *	- returns 0
 ***********************************************************************************/
int copyNodeToSaved(LocusData* locusData, int nodeId, unsigned short recalcConditionals) {

  if(recalcConditionals)		copyNodeConditionals(locusData,nodeId);

  //	printf("Copying node %d, to saved version (%d node copied).\n",nodeId,locusData->savedVersion.numChangedNodes+1);
	
  locusData->savedVersion.changedNodeIds[ locusData->savedVersion.numChangedNodes++ ] = nodeId;				
  locusData->savedVersion.savedNodes[nodeId]->age		 = locusData->nodeArray[nodeId]->age;
  locusData->savedVersion.savedNodes[nodeId]->father	 = locusData->nodeArray[nodeId]->father;
  locusData->savedVersion.savedNodes[nodeId]->leftSon	 = locusData->nodeArray[nodeId]->leftSon;
  locusData->savedVersion.savedNodes[nodeId]->rightSon = locusData->nodeArray[nodeId]->rightSon;
  return 0;
}
/** end of copyNodeToSaved **/



/***********************************************************************************
 *	copyNodeConditionals
 *	- switches pointers to arrays of conditional probabilities between current and saved
 *		versions of node
 *	- marks node as one whose conditionals were changed (recalculated)
 *	- nodeId is id of node being switched
 *	- returns 1, if node was already marked, and 0  otherwise
 ***********************************************************************************/
int copyNodeConditionals(LocusData* locusData, int nodeId) {
  double* conditionalPointer;

  // if node is already recorded as changed, do nothing
  if(locusData->seqData.numPatterns <= 0 || locusData->savedVersion.recalcConditionals[nodeId])
    return 1;
		
  locusData->savedVersion.changedCondIds[locusData->savedVersion.numChangedConditionals++] = nodeId;
  locusData->savedVersion.recalcConditionals[nodeId] = 1;

  // switch updated and saved versions of node
  conditionalPointer = locusData->nodeArray[nodeId]->conditionalProbs;
  locusData->nodeArray[nodeId]->conditionalProbs = locusData->savedVersion.savedNodes[nodeId]->conditionalProbs;
  locusData->savedVersion.savedNodes[nodeId]->conditionalProbs = conditionalPointer;
	
	
  return 0;
}
/** end of copyNodeConditionals **/



/***************************************************************************************************************/
/******                                        END OF FILE                                                ******/
/***************************************************************************************************************/
