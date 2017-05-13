#ifndef ALIGNMENT_PROCESSOR_H
#define ALIGNMENT_PROCESSOR_H
/** 
	\file AlignmentProcessor.h 
	Defines structs: LocusProfile, AlignmentData, AlignmentHets
*/

#include <stdio.h>
#include "utils.h"
/******************************************************************************************************/
/******                                      CONSTANTS                                           ******/
/******************************************************************************************************/



/***************************************************************************************************************/
/******                                              DATA TYPES                                           ******/
/***************************************************************************************************************/



/** LocusProfile
   Holds info about site profiles observed in locus alignment
*/
typedef struct LOCUS_PROFILE_STRUCT{
  int		numPatterns;	/**< number of distinct patterns observed in locus alignment */
  char 	name[NAME_LENGTH];	/**< id of the genealogy used to model this locus */
  int		genealogyId;	/**< genealogy id, for sampling representative genealogies */
  int*	patternIds;		    /**< list of patterns Ids */
  int*	patternCounts;	    /**< number of occurrences per observed pattern */
}LocusProfile;


/***************************************************************************************************************/
/******                                  GLOBAL DATA STRUCTURES                                           ******/
/***************************************************************************************************************/



/** AlignmentData
  Holds info about site patterns observed in a series of locus alignments
*/
struct ALIGNMENT_DATA_STRUCT{
  int		numSamples;		    /**< number of maximal samples to read per alignment */
  char**	sampleNames;	    /**< names of all samples */
  unsigned short*	isDiploid;	/**< binary array indicating the diploid samples */
  int		numPatterns;	    /**< number of distinct site patterns observed */
  char**	patternArray;	    /**< list of all site patterns observed */
  int		numLoci;		    /**< number of loci analyzed (or to be analized) */
  LocusProfile*	locusProfiles;	/**< a list of locus profiles */
};
extern struct ALIGNMENT_DATA_STRUCT AlignmentData;


/** AlignmentHets
  Holds info about site patterns observed in a series of locus alignments
*/
struct PHASED_PATTERNS_STRUCT{
  int		numHaploids;		/**< number of haploids representing samples */
  int		numPhasedPatterns;	/**< number of site patterns with het genotype */
  int*	numPhases;				/**< number of phases per pattern (only indicated in first phase) */
  char**	patternArray;		/**< list of all phased site patterns observed */
  int		numLoci;			/**< number of loci analyzed (or to be analized) */
  LocusProfile*	locusProfiles;	/**< a list of locus profiles per phased pattern */
};

extern struct PHASED_PATTERNS_STRUCT PhasedPatterns;


/******************************************************************************************************/
/******                                FUNCTION DECLARATIONS                                     ******/
/******************************************************************************************************/



/**	initAlignmentData
    Initializes alignment data structures and initial space
    @param numLoci Number of Loci to analyze
    @param numSamples Total number of samples to analyze
    @param initSeqLength Tentative maximum  sequence length of Alignment (treated dynamically)
    @param initNumPatterns Tentative maximum number of patterns in Alignment (treated dynamically)
    @param sampleNames Array of char strings containing names of all the samples defined in the control file
    @return 0 if successful (-1 if allocation problems)
*/
int	initAlignmentData(int numLoci, int numSamples, int initSeqLength, int initNumPatterns, char** sampleNames);



/**	finalizeAlignmentData
    Finalizes data structure after all preprocessing is done (mostly compact memory and frees unnecessary memory usage)
    @return 0
*/
int	finalizeAlignmentData();



/**	freeAlignmentData
    Frees all memory allocated for alignement data
	@return 0
*/
int	freeAlignmentData();



/**	printPatterns
    Prints array containing patterns to stdout
    @param patternArray Array containing patterns to print
    @param patternCounts Array containing number of occurances for each pattern
    @param numPatterns Number of patterns containing in patternArray
    @param numSamples Number of samples per pattern
*/
void	printSitePatterns(char** patternArray, int* patternCounts, int numPatterns, int numSamples);



/**	printLocusProfiles
    Prints profiles of all loci
*/
void	printLocusProfiles();


/** printAlignmentError
    Prints content of error message encountered durring alignment 
*/
void printAlignmentError();


/**	readSeqFile
    Reads series of locus alignment from file, initializes Alignment data structures, performs initial processing of all alignments into site patterns
    @param seqFileName Path to file containing sequences to read in
    @param numSample Number of samples to read according to the control file
    @param sampleNames List of sample names from the control file, used to associate sequences with samples
    @param numLociToRead Maximum number of loci to read if defined in control file (ignored if -1)
    @warning has to be called before all other processing procedures can be called
    @note if numLociToRead is positive and smaller than number of loci in file, reads only the first numLociToRead loci
    @return 0 if successful (-1 otherwise)
*/
int	readSeqFile(const char* seqFileName, int numSample, char** sampleNames, int numLociToRead);



/**	readSeqs
    Reads sequences of a single locus from file and writes them in preallocated internal space
    @param seqFile File descriptor containing sequences to read from
    @param numSeqs Number of sequences in this locus to read from file
    @param seqLength Length of each sequence to read from file
    @param seqArray Preallocated space used to save sequences read from file
    @param locus Number defining current locus that is being read form file
    @note Assumes each sequence is preceded by a name in a predefined list
    @note Orders sequences according to their names, using AlignmentData.sampleNames as a reference from names to indices
    @note Reads sequence, while ignoring any white spaces (except newline)
    @note Capitalizes all bases and checks to see if they are legitimate nucleotides or ambiguities:T,C,A,G ; U,Y,R,M,K,S,W,H,B,V,D ; N,?,-
    @return 0 if OK (-1 if file in bad format)
*/
int	readSeqs(FILE* seqFile, int numSeqs, int seqLength, char** seqArray, int locus);


/**	processLocusAlignment
    Reads through alignment columns and identifies site patterns (according to JC symmetries)
    @param seqArray Single locus alignment, read in from file, these are the sequences to process
    @param seqLength Length of the alignment
    @param locusProfile allocated LocusProfile to populate with alignment data
    @note Adds new site patterns to patternArray
    @note Records locus pattern profile in locusProfile
    @return 0 if successful (-1 otherwise)
*/
int	processLocusAlignment(char** seqArray, int seqLength, LocusProfile* locusProfile);



/**	processHetPatterns
    Creates a phased version of het patterns
    @param patternArray UnPhased Het patterns to process
    @param patternCounts For each pattern in patternArray, the number of times it occured
    @param numPatterns Number of patterns in patternArray
    @param breakSymmetries If == 1, takes into consideration symmetry breaking schemes for hets
    @param phasedPatternArray_ptr Phased version of hets saved here (array of length maxNumPhasedPatterns pre-allocated), more space allocated if needed and maxNumPhasedPatterns is updated accordingly
    @param numPhasesArray_ptr For each pattern the number of phases associated with it (pointer to array - may be reallocated by function)
    @param maxNumPhasedPatterns Maximum number of phased patterns (may be updated by function)
    @note Assumes first element in array (phasedPatternArray) points to entire space needed.
    @return Num phased patterns if successful (-1 otherwise)
*/
int	processHetPatterns(char** patternArray, int* patternCounts, int numPatterns, unsigned short breakSymmetries, char*** phasedPatternArray_ptr, int** numPhasesArray_ptr, int* maxNumPhasedPatterns);

  

/**	countInformativePatterns
    Returns number of informative sites in subset of samples.
    @param includeSample Binary array of size AlignmentData.numSamples, that selects subset of samples in which to look for informative sites, a 1 indicates to include the site, a 0 means not to include the site
    @return Number of informative sites in subset of samples
*/
int	countInformativePatterns(unsigned short* includeSample);



/**	fourGameteTest
    Performs the 4-gamete test on all loci
    @note Prints out potential violations (some printed cases might not be actual violations - need to check by eye)
    @return 0 if all is OK, and -1 otherwise
*/
int	fourGameteTest();



/**	getPatternTypes
    Sorts patterns into types and prints how many columns are observed of each type
    @note A type consists of number of occurrences of each base (in decreasing order)
    @return 0 if all ok (-1 otherwise)
*/
int	getPhasedPatternTypes();



/***************************************************************************************************************/
/******                                        END OF FILE                                                ******/
/***************************************************************************************************************/

#endif
