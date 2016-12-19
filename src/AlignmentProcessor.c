/** 
   \file AlignmentProcessor.c
   Functions to read in and process sequence file, handle haploids, find patterns, handle ambiguity characters, perform tests
*/



/******************************************************************************************************/
/******                                      INCLUDES                                            ******/
/******************************************************************************************************/
#include "utils.h"
#include "AlignmentProcessor.h"

/***************************************************************************************************************/
/******                                              DATA TYPES                                           ******/
/***************************************************************************************************************/


/** Different types of bases */
typedef enum{
  NUCLEOTIDE,		/**< nucleotide base (A,C,T,G) */
  COMPLETE_AMBIG,	/**< complete ambiguity character (N,?,-) */
  PARTIAL_AMBIG,	/**< partial ambiguity character (Y,R,M,K,S,W,H,B,V,D) */
  NO_BASE			/**< illegal base type */
} BASE_TYPE;



/***************************************************************************************************************/
/******                                  GLOBAL DATA STRUCTURES                                           ******/
/***************************************************************************************************************/



/** GlobalSpace
    Holds pointers to global space used by various data processing procedures
*/
struct GLOBAL_SPACE_STRUCT{
  char 	errorMessage[1000];  	/**< error message for alignment preprocessing - OVERFLOW IS NEVER CHECKED FOR !!! */
  char* 	errorMessageEnd;	/**< end of error message for alignment preprocessing (for appending messages) */
  int		maxSeqLength;		/**< maximum sequence length considered */
  int*	intArray;			    /**< global integer array for various purposes */
  char*	seqSpace;			    /**< space for all sequences */
  int		maxNumPatterns;		/**< maximum number of patterns that can be held in patternSpace */
  //	char*	patternSpace;  	/**< space for all patterns */
  char*	fourColumns;		    /**< work space for 4 alignment columns */
}AlignmentGlobal;



/** baseTransformation
    Holds all 24 permutations of base identities, and their projections on ambiguity characters
    @note Used for computing site patterns under JC
*/
short baseTransformations[24][15];


/** cannonizedBaseSymbols
    Standardized order of base characters (including ambiguities)
*/
static const char cannonizedBaseSymbols[] = "TCAGYWKMSRVDBHN";
//static const char cannonizedBaseSymbols[] = "TCAGYKWSMRVDBHN";   BAD ORDER !!!



/******************************************************************************************************/
/******                                INTERNAL FUNCTION DECLARATION                             ******/
/******************************************************************************************************/




/** cannonizeJCpattern
    Canonizes an alignment column according to JC symmetry
    @param column String containing bases from column of alignment (no termination character required)
    @param pattern Pre-allocated string into which to write pattern
    @param numSeqs Number of sequences in alignment (length of column)
    @warning makes use of global variables 'baseTransformations' and 'cannonizedBaseSymbols'
    @return 0 if all OK, and -1 otherwise
*/
int cannonizeJCpattern(const char* column, char* pattern, int numSeqs);

/** Return the type of base (nucleotide, ambiguity, ...) for a given base (i.e. A, T, R, ...)
    @param base Char that is a valid base or ambiguity
    @return A base type telling if it is nucleotide, ambituity, etc.
*/
BASE_TYPE getBaseType(char base);

/** Finds alignment column in given pattern array.
    @param column String of bases for column (no terminating char)
    @param patternArray List of patterns to search through for column (2D Array of size numPatterns * numSeqs)
    @param numSeqs Number of sequences represented in column (length of column string)
    @param numPatterns Number of patterns in patternArray
    @returns column id of patternArray identical to column, if exists, otherwise -1
*/
int findPattern(const char* column, char** patternArray, int numSeqs, int numPatterns);

/** Initializes global baseTransformation array
    Initializes global 2D array for base transformations where every row of the array is a permutation of A,C,G,T (24 total).
    The row describes the permutation and its impact on the ambiguity character
*/
void initializeBaseTransformations();

/** Increases (doubles) the size of globally allocated array for patterns 
    @return 0 if OK, -1 if error (error will be saved in AlignmentGlobal.errorMssageEnd)
    @note called only when adding pattern to array (in processLocusAlignment())
*/
int increasePatternArraySize();

/** Computes symmetry breaking scheme for hets.
    Goes through the site patterns containing het genotypes and determines 
    a het breaking scheme such that each diploid sequence is arbitrarily 
    phased at no more than a SINGLE COLUMN in each locus.
    @parm patternArray Het patterns to compute breaks for
    @param patternCounts The number of occurances for each pattern
    @param numPatterns Number of total patterns in patternArray
    @param numSamples Number of samples in each pattern in patternArray
    @param symmetryBreaks Computed symmetry breaks written out here (2D array pre-allocated size numPatterns * numSamples) 
    @note symmetry breaks indicate hets to arbitrary phase per pattern
**/
int	computeHetSymmetryBreaks(char** patternArray, int* paternCounts, int numPatterns, int numSamples, unsigned short** symmetryBreaks);

/** Computes all phases of a column (given as haploid samples)
    @param column String containing bases to make up a column of the pattern array with each sample a haploid
    @param numHaploids Number of samples in the column (length of column string)
    @param perturbSample Binary array of indicating which haploid pairs to pertub (length numHaploids).  If entry is 1, then perterb that haploid with the previous one.
    @param phasedColumns Phased patterns written here (pre-allocated)
    @warning writes on top of perturbSample to go over all states of perturbations, but returns it to original state
    @return Total number of phases if OK, -1 if error (error stored in AlignmentGlobal.errorMessageEnd)
*/
int getAllPhases(char* column, int numHaploids, unsigned short* perturbSample, char** phasedColumns);

/** Translates ambiguity character to a sequence of bases in given string pointer 
    @param ch Ambiguity character
    @param outp Sequence of bases that represent the ambiguity character made up only of (A,C,G,T)
    @note if character is not a nucleotide symbol or 2-wise ambiguity symbol, writes 2 'N's
    @return Type of input character (ambiguity, nucleotide, etc.)
 */
BASE_TYPE translateAmbiguity(char ch, char* outp);

/** Translates base pair to ambiguity character
    @param String of 2 characters representing a base pair (i.e. AT) to translate into an ambiguity
    @return Ambiguity character (i.e. R) or if error return 'X'
*/
char translateToAmbiguity(char* basePair);

/** Perform 4 gamete test between two potential sites 
    @param patternId1 ID number of first pattern in AlignmentData.patternArray
    @param patternID2 ID number of second pattern in AlignmentData.patternArray
    @return 0 if passed, 1 if failed without doubt, 2 if passes due to greedy phasing
*/
int	twoSiteFourGameteTest(int patternId1, int patternId2);



/******************************************************************************************************/
/******                                FUNCTION IMPLEMENTATION                                     ******/
/******************************************************************************************************/



/***********************************************************************************
 *	initAlignmentData
 * 	- initializes alignment data structures and initial space
 *	- receives number of samples and number of loci to analyze
 *	- also receives TENTATIVE sequence length and total number of patterns (which are treated dynamically)
 *	- returns 0 if successful (-1 if allocation problems)
 ***********************************************************************************/
int	initAlignmentData(int numLoci, int numSamples, int initSeqLength, int initNumPatterns, char** sampleNames) {
  int sample;
	
  AlignmentData.numSamples = numSamples;
  AlignmentData.numLoci = numLoci;
  AlignmentData.numPatterns = 0;
	
  AlignmentData.sampleNames = NULL;
  AlignmentData.isDiploid = NULL;
  AlignmentData.patternArray = NULL;
  AlignmentData.locusProfiles = NULL;
	
  AlignmentGlobal.errorMessageEnd = AlignmentGlobal.errorMessage;
  AlignmentGlobal.maxSeqLength = initSeqLength;
  AlignmentGlobal.maxNumPatterns = initNumPatterns;
	
  AlignmentGlobal.seqSpace = NULL;
  //	AlignmentGlobal.patternSpace = NULL;
  AlignmentGlobal.intArray = NULL;
  AlignmentGlobal.fourColumns = NULL;
	
  PhasedPatterns.numLoci = 0;
  PhasedPatterns.numHaploids = 0;
  PhasedPatterns.numPhasedPatterns = 0;
  PhasedPatterns.numPhases = NULL;
  PhasedPatterns.patternArray = NULL;
  PhasedPatterns.locusProfiles = NULL;

  // allocate space
	
  AlignmentData.isDiploid = (unsigned short*)malloc(numSamples*sizeof(unsigned short));
  if(AlignmentData.isDiploid == NULL) {
    AlignmentGlobal.errorMessageEnd += 
      sprintf(AlignmentGlobal.errorMessageEnd,"Out Of Memory AlignmentData.isDiploid at initAlignmentData.\n");
    return -1;
  }			
	
  if(sampleNames != NULL) {
    AlignmentData.sampleNames = (char**)malloc(numSamples*sizeof(char*));
    if(AlignmentData.sampleNames == NULL) {
      AlignmentGlobal.errorMessageEnd += 
        sprintf(AlignmentGlobal.errorMessageEnd,"Out Of Memory AlignmentData.sampleNames at initAlignmentData.\n");
      return -1;
    }			
	
    AlignmentData.sampleNames[0] = (char*)malloc(numSamples*NAME_LENGTH*sizeof(char));
    if(AlignmentData.sampleNames[0] == NULL) {
      AlignmentGlobal.errorMessageEnd += 
        sprintf(AlignmentGlobal.errorMessageEnd,"Out Of Memory AlignmentData.sampleNames[0] at initAlignmentData.\n");
      return -1;
    }
    // initialize seqNames array
    for(sample = 0; sample<numSamples; sample++) {
      AlignmentData.sampleNames[sample] = AlignmentData.sampleNames[0] + sample*NAME_LENGTH;
      strncpy(AlignmentData.sampleNames[sample],sampleNames[sample], NAME_LENGTH-1);
      //			printf("Sample %d is %s.\n", sample+1, AlignmentData.sampleNames[sample]);
	  if(sampleNames[sample] == NULL || sampleNames[sample][0] == '\0') {
		  if(sample == 0) {
      		AlignmentGlobal.errorMessageEnd += 
        		sprintf(AlignmentGlobal.errorMessageEnd,"First sample must have a name specified.\n");
			return -1;
		  }
		  AlignmentData.isDiploid[sample-1] = AlignmentData.isDiploid[sample] = 1;
	  } else {
		  AlignmentData.isDiploid[sample] = 0;
	  }
		  
    }
  }



  AlignmentData.locusProfiles = (LocusProfile*)malloc(numLoci*sizeof(LocusProfile));
  if(AlignmentData.locusProfiles == NULL) {
    AlignmentGlobal.errorMessageEnd += 
      sprintf(AlignmentGlobal.errorMessageEnd,"Out Of Memory AlignmentData.locusProfiles at initAlignmentData.\n");
    return -1;
  }			

  AlignmentData.patternArray = (char**)malloc(initNumPatterns*sizeof(char*));
  if(AlignmentData.patternArray == NULL) {
    AlignmentGlobal.errorMessageEnd += 
      sprintf(AlignmentGlobal.errorMessageEnd,"Out Of Memory AlignmentData.patternArray at initAlignmentData.\n");
    return -1;
  }			

  AlignmentData.patternArray[0] = (char*)malloc(initNumPatterns*numSamples*sizeof(char));
  if(AlignmentData.patternArray[0] == NULL) {
    AlignmentGlobal.errorMessageEnd += 
      sprintf(AlignmentGlobal.errorMessageEnd,"Out Of Memory AlignmentData.patternArray[0] at initAlignmentData.\n");
    return -1;
  }			

  AlignmentGlobal.seqSpace = (char*)malloc(initSeqLength*numSamples*sizeof(char));
  if(AlignmentGlobal.seqSpace == NULL) {
    AlignmentGlobal.errorMessageEnd += 
      sprintf(AlignmentGlobal.errorMessageEnd,"Out Of Memory AlignmentGlobal.seqSpace at initAlignmentData.\n");
    return -1;
  }			

  AlignmentGlobal.fourColumns = (char*)malloc(4*(numSamples+1)*sizeof(char));
  if(AlignmentGlobal.fourColumns == NULL) {
    AlignmentGlobal.errorMessageEnd += 
      sprintf(AlignmentGlobal.errorMessageEnd,"Out Of Memory AlignmentGlobal.fourColumns at initAlignmentData.\n");
    return -1;
  }			
	
  AlignmentGlobal.intArray = (int*)malloc(initSeqLength*2*sizeof(int));
  if(AlignmentGlobal.intArray == NULL) {
    AlignmentGlobal.errorMessageEnd += 
      sprintf(AlignmentGlobal.errorMessageEnd,"Out Of Memory AlignmentGlobal.intArray at initAlignmentData.\n");
    return -1;
  }

  for(sample = 0; sample<numLoci; sample++) {
    AlignmentData.locusProfiles[sample].patternIds = NULL;
    AlignmentData.locusProfiles[sample].patternCounts = NULL;
  }
		
  // compute all base transformations
  initializeBaseTransformations();
	
  return 0;

}
/* end of initAlignmentData */



/***********************************************************************************
 *	finalizeAlignmentData
 * 	- finalizes data structure after all preprocessing is done
 *	- mostly compact memory and frees unnecessary memory usage
 *	- returns 0
 ***********************************************************************************/
int	finalizeAlignmentData() {
	
  // TBA !!!
  return 0;
}
/* end of finalizeAlignmentData */



/***********************************************************************************
 *	freeAlignmentData
 * 	- frees all memory allocated for alignment data
 *	- returns 0
 ***********************************************************************************/
int	freeAlignmentData() {
  int locus;
	
  // frees locus Profiles arrays
  if(AlignmentData.locusProfiles != NULL) {
    for(locus = 0; locus<AlignmentData.numLoci; locus++) {
      if(AlignmentData.locusProfiles[locus].patternIds != NULL) {
        free(AlignmentData.locusProfiles[locus].patternIds);
      }
    }
    free(AlignmentData.locusProfiles);
  }
  AlignmentData.locusProfiles = NULL;
	
  if(AlignmentData.isDiploid != NULL) {
    free(AlignmentData.isDiploid);
  }			
  AlignmentData.isDiploid = NULL;
	
  if(AlignmentData.sampleNames != NULL) {
    if(AlignmentData.sampleNames[0] != NULL) {
      free(AlignmentData.sampleNames[0]);
    }
    free(AlignmentData.sampleNames);
  }			
  AlignmentData.sampleNames = NULL;
	
	
  if(AlignmentData.patternArray != NULL) {
    if(AlignmentData.patternArray[0] != NULL) {
      free(AlignmentData.patternArray[0]);
    }			
    AlignmentData.patternArray[0] = NULL;
    free(AlignmentData.patternArray);
  }			
  AlignmentData.patternArray = NULL;

  if(PhasedPatterns.numPhases != NULL) {
    free(PhasedPatterns.numPhases);
  }			
  PhasedPatterns.numPhases = NULL;

  if(PhasedPatterns.numPhases != NULL) {
    free(PhasedPatterns.numPhases);
  }			
  PhasedPatterns.numPhases = NULL;

  if(PhasedPatterns.patternArray != NULL) {
    if(PhasedPatterns.patternArray[0] != NULL) {
      free(PhasedPatterns.patternArray[0]);
    }			
    PhasedPatterns.patternArray[0] = NULL;
    free(PhasedPatterns.patternArray);
  }			
  PhasedPatterns.patternArray = NULL;

  if(PhasedPatterns.locusProfiles != NULL) {
    for(locus = 0; locus<PhasedPatterns.numLoci; locus++) {
      if(PhasedPatterns.locusProfiles[locus].patternIds != NULL) {
        free(PhasedPatterns.locusProfiles[locus].patternIds);
      }
    }
    free(PhasedPatterns.locusProfiles);
  }
  PhasedPatterns.locusProfiles = NULL;
	
  if(AlignmentGlobal.seqSpace != NULL) {
    free(AlignmentGlobal.seqSpace);
  }			
  AlignmentGlobal.seqSpace = NULL;
	
  if(AlignmentGlobal.fourColumns != NULL) {
    free(AlignmentGlobal.fourColumns);
  }			
  AlignmentGlobal.fourColumns = NULL;
	
  if(AlignmentGlobal.intArray != NULL) {
    free(AlignmentGlobal.intArray);
  }
  AlignmentGlobal.intArray = NULL;

  return 0;
}
/* end of freeAlignmentData */



/***********************************************************************************
*	printSitePatterns
* 	- prints pattern array
***********************************************************************************/
void	printSitePatterns(char** patternArray, int* patternCounts, int numPatterns, int numSamples) {
	int sample;
	int patt;
	
	printf("%d patterns:\n",numPatterns);
	printf("\n");
	for(sample=0; sample<numSamples; sample++) {
		printf("%3d ",sample+1);
		for(patt=0; patt<numPatterns; patt++) {
			printf("   %c ",patternArray[patt][sample]);
		}
		printf("\n");		
	}
	printf("   ");
	for(patt=0; patt<numPatterns; patt++) {
		printf("%5d",patternCounts[patt]);
	}
	printf("\n");
	return;
}
/* end of printSitePatterns */



/***********************************************************************************
 *	printLocusProfiles
 * 	- prints profiles of all loci
 ***********************************************************************************/
void	printLocusProfiles() {
  int locus;
  int patt;
	
  printf("%d loci:\n",AlignmentData.numLoci);

  for(locus=0; locus<AlignmentData.numLoci; locus++) {
    printf("%5d, %2d patterns: ",locus+1,AlignmentData.locusProfiles[locus].numPatterns);
    for(patt=0; patt<AlignmentData.locusProfiles[locus].numPatterns; patt++) {
      printf(" %d (X %d) ",AlignmentData.locusProfiles[locus].patternIds[patt]+1, AlignmentData.locusProfiles[locus].patternCounts[patt]);
    }
    printf("\n");		
  }
  printf("\n");
  return;
}
/* end of printLocusProfiles */



/***********************************************************************************
 *	readSeqFile
 * 	- reads series of locus alignment from file
 *	- initializes Alignment data structures
 *	- has to be called before all other processing procedures can be called
 *	- performs initial processing of all alignments into site patterns
 *	- receives the list of sample names, to associate sequences with samples
 *	- if numLociToRead is positive and smaller than number of loci in file, 
 *		reads only the first numLociToRead loci
 *	- returns 0 if successful (-1 otherwise)
 ***********************************************************************************/
int	readSeqFile(const char* seqFileName, int numSamples, char** sampleNames, int numLociToRead) {
  
  int INIT_SEQ_LENGTH = 1000;
  int INIT_NUM_PATTERNS = 1000;
  
  FILE *fseq=fopen(seqFileName,"r");
  char *token;
  char** seqArray;	// array of pointers to sequences in locus alignment
  int numLoci, sample;
  int locus, seqLength, numLocusSamples;
  int res;
  
  int* sampleSeqInFile = NULL;
  
  char line[STRING_LENGTH];
  
  if(fseq == NULL) {
    fprintf(stderr, "Error: Could not find sequence file '%s' in readSeqFile().\n",seqFileName);
    return -1;
  }
  
  seqArray = (char**)malloc(numSamples*sizeof(char*));
  if(seqArray == NULL) {
    fprintf(stderr, "Error: Out Of Memory allocating seqArray in readSeqFile().\n");
	fclose(fseq);
    return -1;
  }
  
  sampleSeqInFile = (int*)malloc(numSamples*sizeof(int));
  if(sampleSeqInFile == NULL) {
    fprintf(stderr, "Error: Out Of Memory allocating sampleSeqInFile in readSeqFile().\n");
    free(seqArray);
    fclose(fseq);
    return -1;
  }
  
  for(sample=0; sample<numSamples; sample++) {
	  if(sampleNames[sample] == NULL || sampleNames[sample][0] == '\0') {
		  // for second sample reserved for diploid genomes
		  sampleSeqInFile[sample] = 1;
	  } else {
		  sampleSeqInFile[sample] = 0;
	  }
  }
		  
  
  //Discard empty lines at beginning of file
  while(fgets(line, STRING_LENGTH, fseq) == NULL && !feof(fseq));

  //Get first string of non-empty line
  token = strtokCS(line, parseFileDelims);
  if(token == NULL) {
    fprintf(stderr, "\nError: Unexpected End of File when trying to read number of loci from seq file.\n");
    free(seqArray);
    free(sampleSeqInFile);
    fclose(fseq);
    return -1;
 }
  //See if the string is the expected integer
  res = sscanf(token, "%d", &numLoci);
  if (res != 1) {
    fprintf(stderr, "\nError: Expected number of loci when reading sequence file, got %s\n", token);
    free(seqArray);
    free(sampleSeqInFile);
    fclose(fseq);
    return -1;
  }
  if(numLoci <= 0) {
    fprintf(stderr, "\nError: At least one locus must be specified in the sequence file.\n");
    free(seqArray);
    free(sampleSeqInFile);
    fclose(fseq);
    return -1;
  }

  if(numLociToRead >0 && numLociToRead < numLoci) {
    printf("Reading sequence data...  %d first loci (out of %d) in file, as specified by user.\n", numLociToRead, numLoci);
    numLoci = numLociToRead;
  } else {
    printf("Reading sequence data...  %d loci, as specified in sequence file %s.\n", numLoci, seqFileName);
  }
  
  if (verbose)
    printf("Initializing alignment data structure for %d samples and %d loci...", numSamples, numLoci);
  res = initAlignmentData(numLoci, numSamples, INIT_SEQ_LENGTH, INIT_NUM_PATTERNS, sampleNames);
  if(res<0) {
    printAlignmentError();
    freeAlignmentData();
    free(seqArray);
    free(sampleSeqInFile);
    fclose(fseq);
    return -1;
  }
  
  printf("Reading loci (.=100 loci): ");
  for(locus=0; locus<numLoci; locus++) {
    
    token = NULL;
    while(!feof(fseq)) {
      if (fgets(line, STRING_LENGTH, fseq) == NULL) {
       token = NULL;
       break;
      }
      token = strtokCS(line, parseFileDelims); 
      if (token != NULL)
        break;
    }

    if(token == NULL) {
      	fprintf(stderr, "\nError: Sequence file says to use %d loci, but the sequence file only contains %d loci.\n", numLoci, locus);
    	printAlignmentError();
    	freeAlignmentData();
    	free(seqArray);
    	free(sampleSeqInFile);
    	fclose(fseq);
      	return -1;
    } else if (strlen(token) >= (NAME_LENGTH -2))
      fprintf(stderr, "\nWarning: Locus names can only be %d characters long. Truncated loci name for loci #%d.\n", NAME_LENGTH-1, locus+1);
    strncpy(AlignmentData.locusProfiles[locus].name, token, NAME_LENGTH);

    token = strtokCS(NULL, parseFileDelims);
    if(token == NULL) {
		fprintf(stderr, "\nError: Unexpected end of file when trying to read number of samples for locus %d.\n", locus+1);
    	freeAlignmentData();
    	free(seqArray);
		free(sampleSeqInFile);
    	fclose(fseq);
 		return -1;
    }
    res = sscanf(token, "%d", &numLocusSamples);
    if (res != 1) {
 	    fprintf(stderr, "\nError: Expected number of locus samples, got '%s'.\n", token);
    	freeAlignmentData();
    	free(seqArray);
		free(sampleSeqInFile);
    	fclose(fseq);
 		return -1;
    }

    if(numLocusSamples <= 0) {
      	fprintf(stderr, "\nError: Every Locus must have one or more samples, locus #%d did not.\n", locus+1);
    	freeAlignmentData();
    	free(seqArray);
		free(sampleSeqInFile);
    	fclose(fseq);
 		return -1;
    }
    
    token = strtokCS(NULL, parseFileDelims);
    if(token == NULL) {
     	fprintf(stderr, "\nError: Unexpected end of file when trying to read sequence length for locus %d.\n", locus+1);
    	freeAlignmentData();
    	free(seqArray);
		free(sampleSeqInFile);
    	fclose(fseq);
 		return -1;
    }
    res = sscanf(token, "%d", &seqLength);
    if (res != 1) {
      	fprintf(stderr, "\nError: Expected sequence length, got '%s'.\n", token);
    	freeAlignmentData();
    	free(seqArray);
		free(sampleSeqInFile);
    	fclose(fseq);
 		return -1;
    }
    /*		if(numLocusSamples > numSamples) {
            AlignmentGlobal.errorMessageEnd += 
            sprintf(AlignmentGlobal.errorMessageEnd,"Number of samples %d in locus %d exceeds maximum of %d.\n",numLocusSamples,locus+1,numSamples);
            printf("\n");
            printAlignmentError();
            freeAlignmentData();
            free(seqArray);
            return -1;
            }
    */		

	// read sequences from file into memory (seqArray)
    res = readSeqs(fseq, numLocusSamples, seqLength, seqArray, locus);
    if(res<0) {
      AlignmentGlobal.errorMessageEnd += 
        sprintf(AlignmentGlobal.errorMessageEnd,"Error occurred while reading sequences for locus %d.\n",locus+1);
      printAlignmentError();
      freeAlignmentData();
	  free(sampleSeqInFile);
      free(seqArray);
      fclose(fseq);
      return -1;
    } else {
      for(sample = 0; sample<AlignmentData.numSamples; sample++) {
        if(sampleNames[sample] != NULL && sampleNames[sample][0] != '\0' && seqArray[sample] != NULL)	
          sampleSeqInFile[sample] = 1;
      }
    }
    
    // transform sequence info to patterns, and collapse patterns for JC
    res = processLocusAlignment(seqArray, seqLength, &(AlignmentData.locusProfiles[locus]));
    if(res<0) {
      AlignmentGlobal.errorMessageEnd += 
        sprintf(AlignmentGlobal.errorMessageEnd,"Error occurred while processing alignment of locus %d to site patterns.\n",locus+1);
      printAlignmentError();
      freeAlignmentData();
	  free(sampleSeqInFile);
      free(seqArray);
      fclose(fseq);
      return -1;
    }
    //		printf("+");
    //		fflush(stdout);
    
    if((locus+1)%100 == 0) {
      printf(".");
      if((locus+1)%1000 == 0) {
        printf(" ");
        if((locus+1)%10000 == 0) {
          printf("\n");
        }
      }
    }
    fflush(stdout);
  } // end of for(locus)

  // print samples that are not represented in seq file
  for(sample = 0; sample<AlignmentData.numSamples; sample++) {
    if(sampleSeqInFile[sample] == 0) {  
      AlignmentGlobal.errorMessageEnd += 
        sprintf(AlignmentGlobal.errorMessageEnd,"Sample name '%s' was defined in the control file, but no samples for this name exist in the sequence file.\n",sampleNames[sample]);
      printAlignmentError();
	  fclose(fseq);
	  free(seqArray);
	  free(sampleSeqInFile);
      return -1;
    }
  }

  
  fclose(fseq);
  free(seqArray);
  free(sampleSeqInFile);

  if (verbose)  
    printf(" Done.\n");
  else
    printf("\n");
	
  return 0;
    
}
/* end of readSeqFile */



/***********************************************************************************
 *	readSeqs
 * 	- reads sequences from file and writes them in preallocated internal space
 *	- points to this space with seqArray output argument
 * 	- assumes each sequence is preceded by a name in a predefined list
 *	- orders sequences according to their names, using AlignmentData.sampleNames as a reference from names to indices
 * 	- reads sequence, while ignoring any white spaces
 * 	- capitalizes all bases and checks to see if they are legitimate nucleotides or ambiguities:
 * 		T,C,A,G ; U,Y,R,M,K,S,W,H,B,V,D ; N
 * 	- returns 0 if all is fine (-1 if file in bad format)
 ***********************************************************************************/
int	readSeqs(FILE* seqFile, int numSeqs, int seqLength, char** seqArray, int locus) {
	
  // globally saved values
  int numTotalSamples = AlignmentData.numSamples;

  int res, seq, site, seqIndex;
  char sampleName[NAME_LENGTH];
  
  char* seqSpace;
  char ch;
	
  // reallocate space for sequences, if necessary - MOVE TO HIGHER FUNCTION !!!
  if(seqLength > AlignmentGlobal.maxSeqLength) {
    //		printf("Increasing max sequence length to %d.\n",seqLength);
    AlignmentGlobal.maxSeqLength = seqLength;
    free(AlignmentGlobal.seqSpace);
    AlignmentGlobal.seqSpace = (char*)malloc(seqLength*numTotalSamples*sizeof(char));
    if(AlignmentGlobal.seqSpace == NULL) {
      AlignmentGlobal.errorMessageEnd += 
        sprintf(AlignmentGlobal.errorMessageEnd,"Out Of Memory reallocating AlignmentGlobal.seqSpace with seq length %d.\n",seqLength);
      return -1;
    }			
    free(AlignmentGlobal.intArray);
    AlignmentGlobal.intArray = (int*)malloc(2*seqLength*sizeof(int));
    if(AlignmentGlobal.intArray == NULL) {
      AlignmentGlobal.errorMessageEnd += 
        sprintf(AlignmentGlobal.errorMessageEnd,"Out Of Memory reallocating AlignmentGlobal.intArray with seq length %d.\n",seqLength);
      return -1;
    }			
		
  }
  seqSpace = AlignmentGlobal.seqSpace;
	
  // initialize pointers to sequences
  // NULL will indicate a sequence to a specific sample was not read ('N')
  for(seq = 0; seq<numTotalSamples; seq++) {
    seqArray[seq] = NULL;
  }// end of for(seq)
	
  for(seq = 0; seq<numSeqs; seq++) {
    // read sample name, and discard
    res = readStringFromFile(seqFile, NAME_LENGTH, sampleName);
    if (res < 1) {
      AlignmentGlobal.errorMessageEnd += 
        sprintf(AlignmentGlobal.errorMessageEnd,"Encountered unexpected EOF at seq %d.\n",seq+1);
      return -1;
    } else if (res > (NAME_LENGTH-1)) {
      fprintf(stderr, "\nWarning: sample names can only be %d characters long, name of sample %d was truncated to %s.\n", (NAME_LENGTH-1), seq+1, sampleName);
    }
    
    if(AlignmentData.sampleNames == NULL) {
		seqIndex = seq;
	} else {
		for(seqIndex = 0; seqIndex<numTotalSamples; seqIndex++) {
          if(0 == strcmp(sampleName, AlignmentData.sampleNames[seqIndex]))	{
            break;
          }
        }
	}
       
	if(seqIndex >= numTotalSamples) {
        // if no match found, skip sample
        flushLine(seqFile);
        continue;
      }
      
    // initialize pointer to sequence 
    seqArray[seqIndex] = seqSpace + seqIndex*seqLength;

    // skip white spaces before sequence
    for(ch=fgetc(seqFile); ch != EOF && (isspace(ch)> 0); ch = fgetc(seqFile)) { ; }
//    fputc(ch, seqFile);

// read actual sequence (seqLength characters)
    for(site=0; site<seqLength; site++) {
      if(ch == '\n') {
        AlignmentGlobal.errorMessageEnd += 
          sprintf(AlignmentGlobal.errorMessageEnd,"Sequence for sample %s contained only %d bases instead of the expected %d bases as defined in the sequence file.\n", sampleName, site, seqLength);
        return -1;
      }
			
      if(ch == EOF) {
        AlignmentGlobal.errorMessageEnd += 
          sprintf(AlignmentGlobal.errorMessageEnd,"Unexpected EOF while reading sequence %d.\n", seq+1);
        return -1;
      }

      if (isspace(ch) != 0) {
        AlignmentGlobal.errorMessageEnd += 
          sprintf(AlignmentGlobal.errorMessageEnd,"Whitespace found in site %d for sample %s. No whitespaces (tab, space, etc.) permitted inside sequences.\n", site+1, sampleName);
        return -1;
      }
            
      ch = (char)toupper(ch);
/*      if(!AlignmentData.isDiploid[seqIndex] && PARTIAL_AMBIG == getBaseType(ch)) {
        AlignmentData.isDiploid[seqIndex] = 1;
        //				printf("Sample %d (%s) found to be diploid.\n", seqIndex+1, AlignmentData.sampleNames[seqIndex]);
        AlignmentData.numHaploids ++;
      } else 
*/	  if(NO_BASE == getBaseType(ch)) {
        AlignmentGlobal.errorMessageEnd += 
          sprintf(AlignmentGlobal.errorMessageEnd,"Illegal base type '%c' found in site %d of sample %s.\n",ch,site+1,sampleName);
        return -1;
      } else if(PARTIAL_AMBIG == getBaseType(ch) && !AlignmentData.isDiploid[seqIndex]) {
        AlignmentGlobal.errorMessageEnd += 
          sprintf(AlignmentGlobal.errorMessageEnd,"Ambiguity character '%c' found in site %d of haploid sample %s.\n",ch,site+1,sampleName);
        return -1;
	  }
		  

      seqArray[seqIndex][site] = ch;
      ch = fgetc(seqFile);
   }// end of for(site)
		
    // discard of remainder of line
	if(!isspace(ch)) {
        AlignmentGlobal.errorMessageEnd += 
			sprintf(AlignmentGlobal.errorMessageEnd,"Sequence for sample %s might be too long than specified (%d bases). Found character %c at position %d.\n",
					sampleName, seqLength, ch,seqLength+1);
        return -1;
	}
		
    if(ch != '\n')
		flushLine(seqFile);
		
  }// end of for(seq)

  return 0;
}
/* end of readSeqs */



/***********************************************************************************
 *	processLocusAlignment
 * 	- reads through alignment columns and identifies site patterns (according to JC symmetries)
 *	- adds new site patterns to patternArray
 *	- records locus pattern profile in locusProfile
 * 	- returns 0 if successful (-1 otherwise)
 ***********************************************************************************/
int	processLocusAlignment(char** seqArray, int seqLength, LocusProfile* locusProfile)	{
						
  int numSamples = AlignmentData.numSamples;
  char** patternArray = AlignmentData.patternArray;
  char *column, *pattern;
  int  *patternIds, *patternCounts;
  int seq, site, patt, pattId, numPatterns;
  unsigned short notAllNs;	// flag for informative columns
  int res;

  // use global space
  column = AlignmentGlobal.fourColumns;
  column[numSamples]= '\0';
  pattern = column + numSamples+1;
  pattern[numSamples] = '\0';
  patternIds = AlignmentGlobal.intArray;
  patternCounts = patternIds + seqLength;
	
  numPatterns = 0;	
  for(site=0; site<seqLength; site++) {
    // extract column (and augment with missing data)
    notAllNs = 0;
    for(seq=0; seq<numSamples; seq++) {
      if(seqArray[seq] == NULL) {
        column[seq] = 'N';
      } else {
        column[seq] = seqArray[seq][site];
        if(column[seq] != 'N') 	{
          notAllNs = 1;
        }
      }
    }

    // if column has all-Ns, do not consider in analysis
    if(!notAllNs) {
      //			printf(" %d",site+1);
      continue;
    }
		
    res = cannonizeJCpattern(column, pattern, numSamples);
    if(res < 0) {
      AlignmentGlobal.errorMessageEnd += 
        sprintf(AlignmentGlobal.errorMessageEnd,"Error while canonizing  site %d under JC model.\n",site+1);
      return -1;
    }
    if(res == 1) {
      //			printf("%5d %s\n",site+1,pattern);
    }
    // see if site pattern appears in patternArray
    pattId = findPattern(pattern, patternArray, numSamples, AlignmentData.numPatterns);
    if(pattId >= 0) {
      for(patt=0; patt<numPatterns; patt++) {
        if(patternIds[patt] == pattId) {
          patternCounts[patt]++;
          //					printf("Found old pattern %d, new to locus.\n",pattId+1);
          break;
        }
      }
      if(patt >= numPatterns) {
        //				printf("Found old pattern %d, (%d in locus).\n",pattId+1,patt+1);
        patternIds[numPatterns] = pattId;
        patternCounts[numPatterns] = 1;
        numPatterns++;
      }
    } else {
      // pattern not previously observed in any alignment
      //			printf("Pattern %d %s from column %s.\n",AlignmentData.numPatterns+1,pattern,column);
      //			fflush(stdout);
      patternIds[numPatterns] = AlignmentData.numPatterns;
      patternCounts[numPatterns] = 1;
      numPatterns++;
      if(AlignmentData.numPatterns >= AlignmentGlobal.maxNumPatterns) {
        if(0 > increasePatternArraySize()) {
          return -1;
        }
        //				printf("-max patterns %d-",AlignmentGlobal.maxNumPatterns);
        patternArray = AlignmentData.patternArray;
      }
      patternArray[AlignmentData.numPatterns] = AlignmentData.patternArray[0] + AlignmentData.numPatterns*numSamples;
      //			memcpy((void*)patternArray[AlignmentData.numPatterns],(void*)pattern,numSamples*sizeof(char));
      strncpy(patternArray[AlignmentData.numPatterns],pattern,numSamples);
      AlignmentData.numPatterns++;
    }

  } // end of for(site)
	
	
  if(locusProfile == NULL) {
    AlignmentGlobal.errorMessageEnd += 
      sprintf(AlignmentGlobal.errorMessageEnd,"Supplied locus profile data structure is NULL.\n");
    return -1;
  }
  locusProfile->numPatterns = numPatterns;
  if(numPatterns == 0) {
    //		printf("no patterns??\n");
    locusProfile->patternCounts = locusProfile->patternIds = NULL;
    fflush(stdout);
  } else {
    locusProfile->patternIds  = (int*)malloc(2*numPatterns*sizeof(int));
    if(locusProfile->patternIds == NULL) {
      AlignmentGlobal.errorMessageEnd += 
        sprintf(AlignmentGlobal.errorMessageEnd,"Out Of Memory int array for locusProfile");
      return -1;
    }
    locusProfile->patternCounts = locusProfile->patternIds + numPatterns;
    for(patt = 0; patt<numPatterns; patt++) {
      locusProfile->patternIds[patt] = patternIds[patt];
      locusProfile->patternCounts[patt] = patternCounts[patt];
    }
  }
	
  return 0;
}
/* end of processLocusAlignment */



/***********************************************************************************
 *	processHetPatterns
 * 	- processes all het patterns and creates phased versions of these
 *	- receives a list of patterns to process
 *	- if breakSymmetries == 1, takes into consideration symmetry breaking schemes for hets
 *	- writes phased versions of these in phasedPatternArray (array of length maxNumPhasedPatterns)
 *	- assumes first element in array (phasedPatternArray) points to entire space needed
 *	- if needs more space, allocates it and updates maxNumPhasedPatterns.
 * 	- returns num phased patterns if successful (-1 otherwise)
 ***********************************************************************************/
int	processHetPatterns(char** patternArray, int* patternCounts, int numPatterns, unsigned short breakSymmetries, char*** phasedPatternArray_ptr, int** numPhasesArray_ptr, int* maxNumPhasedPatterns)	{
	
  int numPhasedPatterns, phaseCount;
  int numSamples = AlignmentData.numSamples;
	
  int res, sample, patt, newPatt;
  
  char ** phasedPatternArray = *phasedPatternArray_ptr;
  
  int * numPhasesArray = *numPhasesArray_ptr;
	
  BASE_TYPE baseType;
	
  char*	haploidColumn;
  unsigned short*	perturbHaploid = NULL;
	
  unsigned short** symmetryBreaks = NULL;
	
  // allocate auxilliary memory for symmetry breaking
  symmetryBreaks = (unsigned short**)malloc(numPatterns*sizeof(unsigned short*));
  if(symmetryBreaks == NULL) {
    AlignmentGlobal.errorMessageEnd += 
      sprintf(AlignmentGlobal.errorMessageEnd,"Out Of Memory allocating symmetryBreaks in processHetPatterns().\n");
    return -1;
  }

  perturbHaploid = (unsigned short*)malloc((numPatterns+2)*numSamples*sizeof(unsigned short));
  if(perturbHaploid == NULL) {
    AlignmentGlobal.errorMessageEnd += 
      sprintf(AlignmentGlobal.errorMessageEnd,"Out Of Memory allocating perturbHaploid in processHetPatterns().\n");
    free(symmetryBreaks);
    return -1;
  }
  for(patt=0; patt<numPatterns; patt++) {
    symmetryBreaks[patt] = perturbHaploid + (patt+2)*numSamples;
  }
	
  haploidColumn = (char*)malloc(numSamples*sizeof(char));
  if(haploidColumn == NULL) {
    AlignmentGlobal.errorMessageEnd += 
      sprintf(AlignmentGlobal.errorMessageEnd,"Out Of Memory allocating haploidColumn in processHetPatterns().\n");
    free(perturbHaploid);
    free(symmetryBreaks);
    return -1;
  }

  //	printf("Computing het symmetry breaks...\n");
  // compute symmetry breaking scheme
  res = computeHetSymmetryBreaks(patternArray, patternCounts, numPatterns, numSamples, symmetryBreaks);
  if(res < 0) {
    AlignmentGlobal.errorMessageEnd += 
      sprintf(AlignmentGlobal.errorMessageEnd,"Error while breaking symmetries.\n");
    free(haploidColumn);
    free(perturbHaploid);
    free(symmetryBreaks);
    return -1;
  }

  //	printf("Done.\n");

	
  // compute number of phased patterns
  numPhasedPatterns = 0;
  for(patt=0; patt<numPatterns; patt++) {
    phaseCount = 1;
    for(sample=0; sample<numSamples; sample++) {
      if(getBaseType(patternArray[patt][sample]) == PARTIAL_AMBIG && (!breakSymmetries || !symmetryBreaks[patt][sample])) {
        phaseCount *=2;
      }			
    }
    numPhasedPatterns += phaseCount;
  }

	
  // allocate additional memory, if needed
  if(numPhasedPatterns > *maxNumPhasedPatterns) {
	  if(debug) {
    	printf("reallocating phased pattern array to size %d.\n", numPhasedPatterns);
	  }
    free(phasedPatternArray[0]);
    free(phasedPatternArray);
   free(numPhasesArray);
    numPhasesArray = (int*)malloc(numPhasedPatterns*sizeof(int));
    if(numPhasesArray == NULL) {
      AlignmentGlobal.errorMessageEnd += 
        sprintf(AlignmentGlobal.errorMessageEnd,"Error while reallocating numPhasesArray.\n");
      free(haploidColumn);
      free(perturbHaploid);
      free(symmetryBreaks);
      return -1;
    }
    phasedPatternArray = (char**)malloc(numPhasedPatterns*sizeof(char*));
    if(phasedPatternArray == NULL) {
      AlignmentGlobal.errorMessageEnd += 
        sprintf(AlignmentGlobal.errorMessageEnd,"Error while reallocating phasedPatternArray.\n");
      free(haploidColumn);
      free(perturbHaploid);
      free(symmetryBreaks);
      return -1;
    }
    phasedPatternArray[0] = (char*)malloc(numPhasedPatterns*numSamples*sizeof(char));
    if(phasedPatternArray[0] == NULL) {
      AlignmentGlobal.errorMessageEnd += 
        sprintf(AlignmentGlobal.errorMessageEnd,"Error while reallocating space for phasedPatternArray.\n");
      free(haploidColumn);
      free(symmetryBreaks);
      return -1;
    }
    *phasedPatternArray_ptr = phasedPatternArray;
	*numPhasesArray_ptr = numPhasesArray;
    *maxNumPhasedPatterns = numPhasedPatterns;
  }
			
	
  for(newPatt=0; newPatt<numPhasedPatterns; newPatt++) {
    phasedPatternArray[newPatt] = phasedPatternArray[0] + numSamples*newPatt;
    numPhasesArray[newPatt] = 0;
  }
	

  //	printf("Generating phased versions of all patterns (with symmetry breaks).\n");
  for(patt=0, newPatt=0; patt<numPatterns;patt++) {
    //		printf("pattern %d --> %d\n",patt+1,newPatt+1);
    for(sample=0; sample<numSamples; sample++) {
      if(!AlignmentData.isDiploid[sample]) {
        perturbHaploid[sample] = 0;
        haploidColumn[sample] = patternArray[patt][sample];
      } else {
        perturbHaploid[sample] = 0;
        baseType = translateAmbiguity(patternArray[patt][sample], &haploidColumn[sample]);
        if(baseType == PARTIAL_AMBIG && (!breakSymmetries || !symmetryBreaks[patt][sample])) {
          perturbHaploid[sample+1] = 1;
        } else {
          perturbHaploid[sample+1] = 0;
        }		
        // skip next space reserved for second haploid in diploid sequence
        sample++;
      }
    }
    numPhasesArray[newPatt] = getAllPhases(haploidColumn, numSamples,perturbHaploid,&phasedPatternArray[newPatt]);
    if(numPhasesArray[newPatt] < 1) {
      AlignmentGlobal.errorMessageEnd += 
        sprintf(AlignmentGlobal.errorMessageEnd,"Error while computing phases of pattern %d in processHetPatterns().\n", patt+1);
      free(haploidColumn);
      free(perturbHaploid);
      free(symmetryBreaks);
      return -1;
    }
			
    newPatt += numPhasesArray[newPatt];
  }// end of for(patt)
	
  //	printf("Done.\n");
	
	
  free(haploidColumn);
  free(perturbHaploid);
  free(symmetryBreaks);
						
  return numPhasedPatterns;
}
/* end of processHetPatterns */



/***********************************************************************************
 *	countInformativePatterns()
 * 	- returns number of informative sites in subset of samples indicated by input binary array
 *	- returns 0 if all ok (-1 otherwise)
 ***********************************************************************************/
int	countInformativePatterns(unsigned short* includeSample) {
			
  int numInformative, patt, sample;
  char base;
	
  numInformative = 0;
  for(patt=0; patt<AlignmentData.numPatterns; patt++) {
    base = 'N';
    for(sample=0; sample<AlignmentData.numSamples; sample++) {
      if(includeSample[sample] && AlignmentData.patternArray[patt][sample] != 'N') {
        if(base == 'N') {
          base = AlignmentData.patternArray[patt][sample];
        } else {
          if(base != AlignmentData.patternArray[patt][sample]) {
            numInformative++;
            break;
          }
        }
      }
    }// end of for(sample)
  }// end of for(patt)
  return numInformative;
}
/* end of countInformativePatterns */



/***********************************************************************************
 *	getPatternTypes()
 * 	- sorts patterns into types
 *	- prints how many columns are observed of each type
 *	- a type consists of number of occurrences of each base (in decreasing order)
 *	- returns 0 if all ok (-1 otherwise)
 ***********************************************************************************/
int	getPhasedPatternTypes() {
			
  int numPatterns 		= AlignmentData.numPatterns;
  int numLoci 			= AlignmentData.numLoci;
  int numSamples 			= PhasedPatterns.numHaploids;
  ;
	
  char** 	patternArray 	= PhasedPatterns.patternArray;
	
  BASE_TYPE baseType;
	
  int i, j, patt, pattId, locus, sample, base, type, numTypes;
	
  int** patternTypes;	// 2D array of all distinct pattern types
	
  int* typeCounts;	// the number of columns for each type
	
  int* pattern2type;	// for each pattern indicates its type
	
	
  patternTypes = (int**)malloc(numPatterns*sizeof(int*));
  if(patternTypes == NULL) {
    AlignmentGlobal.errorMessageEnd += 
      sprintf(AlignmentGlobal.errorMessageEnd,"Out Of Memory patternTypes array in getPatternTypes()");
    return -1;
  }
	
  patternTypes[0] = (int*)malloc(numPatterns*6*sizeof(int));
  if(patternTypes[0] == NULL) {
    AlignmentGlobal.errorMessageEnd += 
      sprintf(AlignmentGlobal.errorMessageEnd,"Out Of Memory patternTypes[0] array in getPatternTypes()");
    free(patternTypes);
    return -1;
  }
		
  typeCounts = patternTypes[0] + numPatterns*4;
  pattern2type = patternTypes[0] + numPatterns*5;
	
  numTypes = 0;
  // assign to each pattern a pattern type
  for(patt=0, pattId=0; patt<numPatterns; patt++, pattId+=PhasedPatterns.numPhases[pattId] ) {
    for(base=0; base<4; base++) {
      patternTypes[numTypes][base] = 0;
    }
    for(sample=0; sample<numSamples; sample++) {
      baseType = getBaseType(patternArray[pattId][sample]);
      if(baseType == COMPLETE_AMBIG) continue;
      if(baseType != NUCLEOTIDE) {
        AlignmentGlobal.errorMessageEnd += 
          sprintf(AlignmentGlobal.errorMessageEnd,"Illegal base found: %c.\n",patternArray[pattId][sample]);
        free(patternTypes[0]);
        free(patternTypes);
        return -1;
      }
				
      base = strchr(cannonizedBaseSymbols,patternArray[pattId][sample])-cannonizedBaseSymbols; // base in [0..3]
      patternTypes[numTypes][base]++;
    }
		
    // sort base counts
    for(i=1; i<4; i++) {
      for(j=0; j<i; j++) {
        if(patternTypes[numTypes][i] > patternTypes[numTypes][j]) {
          patternTypes[numTypes][i] = patternTypes[numTypes][i] + patternTypes[numTypes][j];
          patternTypes[numTypes][j] = patternTypes[numTypes][i] - patternTypes[numTypes][j];
          patternTypes[numTypes][i] = patternTypes[numTypes][i] - patternTypes[numTypes][j];
        }
      }
    }
		
    // search for type in previously observed types
    for(type=0; type<numTypes; type++) {
      for(base=0; base<4; base++) {
        if(patternTypes[numTypes][base] != patternTypes[type][base]) break;
      }
			
      if(base == 4) break;
    }
		
    pattern2type[patt] = type;

    if(type == numTypes) {
      //			printf("Type %d found at pattern %d:",numTypes+1,patt+1);
      //			for(base=0; base<4; base++) {
      //				printf(" %d",patternTypes[numTypes][base]);
      //			}
      //			printf("\n");
      typeCounts[numTypes] = 0;
      patternTypes[numTypes+1] = patternTypes[numTypes] + 4;
      numTypes++;
    }
			
  }// end of for(patt)
	
  // compute type counts;
  for(locus=0; locus<numLoci; locus++) {
    for(patt=0; patt<AlignmentData.locusProfiles[locus].numPatterns; patt++) {
      pattId = AlignmentData.locusProfiles[locus].patternIds[patt];
      typeCounts[ pattern2type[pattId] ] += AlignmentData.locusProfiles[locus].patternCounts[patt];
    }
  }// end of for(locus)
	
  // print statistics
  for(type=0; type<numTypes; type++) {
    printf("( %2d ) [ %8d ] <",type+1,typeCounts[type]);
    for(base=0; base<4; base++) {
      printf(" %2d",patternTypes[type][base]);
    }
    printf(" >\n");
  }// end of for(type)
	
  free(patternTypes[0]);
  free(patternTypes);
  return 0;
}
/* end of getPatternTypes */



/***********************************************************************************
 *	fourGameteTest()
 * 	- performs the 4-gamete test on all loci
 *	- prints out potential violations (some printed cases might not be actual violations - need to check by eye)
 *	- returns 0 if all is OK, and -1 otherwise
 ***********************************************************************************/
int	fourGameteTest() {
			
  int numPatterns 		= AlignmentData.numPatterns;
  int numLoci 			= AlignmentData.numLoci;
  int numSamples 			= AlignmentData.numSamples;

  int res;
	
  int isViolated, numViolated;
	
  char** 	patternArray 	= AlignmentData.patternArray;
	
  BASE_TYPE baseType;
	
  int counts[4];
  char ambig[2];
	
  // 0 - non-informative or singleton, 1 - potential conflict, 2 - tri-allelic sites with more than a singleton for the 3rd allele  
  unsigned short* patternStatus;	
	
  int patt, base, i, j, sample;
	
  int locus, patt1, patt2, pattId1, pattId2;
	
  patternStatus = (unsigned short*)malloc(numPatterns*sizeof(unsigned short));
  if(patternStatus == NULL) {
    AlignmentGlobal.errorMessageEnd += 
      sprintf(AlignmentGlobal.errorMessageEnd,"Out Of Memory patternStatus array in fourGameteTest()");
    return -1;
  }
		
	
  // perform base counts
  printf("Classifying sites.\n");
  for(patt=0; patt<numPatterns; patt++) {
    for(base=0; base<4; base++) {
      counts[base] = 0;
    }
    for(sample=0; sample<numSamples; sample++) {
      baseType = getBaseType(patternArray[patt][sample]);
      if(baseType == COMPLETE_AMBIG) continue;
			
      translateAmbiguity(patternArray[patt][sample], ambig);
      for(i=0; i<2; i++) {
        base = strchr(cannonizedBaseSymbols,ambig[i])-cannonizedBaseSymbols; // base in [0..3]
        counts[base]++;
      }
    }
		
    // sort base counts
    for(i=1; i<4; i++) {
      for(j=0; j<i; j++) {
        if(counts[i] > counts[j]) {
          counts[i] = counts[i] + counts[j];
          counts[j] = counts[i] - counts[j];
          counts[i] = counts[i] - counts[j];
        }
      }
    }
		
    if(counts[1] < 2) {
      patternStatus[patt] = 0;
    } else if(counts[2] > 1) {
      patternStatus[patt] = 2;
    } else {
      patternStatus[patt] = 1;
    }
			
		
  }// end of for(patt)
	
  // check for 4-gamete conflicts
  numViolated = 0;
  for(locus=0; locus<numLoci; locus++) {
    //		printf("Locus %d.\n",locus+1); 
    isViolated = 0;
    for(patt1=1; patt1<AlignmentData.locusProfiles[locus].numPatterns; patt1++) {
      pattId1 = AlignmentData.locusProfiles[locus].patternIds[patt1];
      if(patternStatus[pattId1] == 0) {
        continue;
      }
      for(patt2=0; patt2<patt1; patt2++) {
        pattId2 = AlignmentData.locusProfiles[locus].patternIds[patt2];
        if(patternStatus[pattId2] == 0) {
          res = 0;
        } else if(patternStatus[pattId1] == 2 || patternStatus[pattId2] == 2) {
          res = 3;
        } else {
          //					printf("Locus %d, patterns %d and %d.\n",locus+1,pattId1+1,pattId2+1);
          res = twoSiteFourGameteTest(pattId1,pattId2);
          //					printf("res=%d.\n",res);
        }
				
        if(res>0) {
          if(!isViolated) {
            isViolated = 1;
            numViolated++;
          }
          printf("potential conflict at locus %5d %25s, patterns ",locus+1, AlignmentData.locusProfiles[locus].name);
          for(sample=0; sample<numSamples; sample++) {
            printf("%c",patternArray[pattId1][sample]);
          }
          printf(" and ");
          for(sample=0; sample<numSamples; sample++) {
            printf("%c",patternArray[pattId2][sample]);
          }
          printf("- %d\n",res);
        }
      }// end of for(patt2)
    }// end of for(patt1)
  }// end of for(locus)
	
  free(patternStatus);
	
  printf("Total of %d violated loci.\n",numViolated);
  return 0;
}
/* end of fourGameteTest */



/***************************************************************************************************************/
/******                                INTERNAL FUNCTION IMPLEMENTATION                                   ******/
/***************************************************************************************************************/



void printAlignmentError() {
  printf("\n");
  fprintf(stderr, "Error while processing alignments:\n%s",AlignmentGlobal.errorMessage);
  AlignmentGlobal.errorMessageEnd = AlignmentGlobal.errorMessage;
  return;
}



/***********************************************************************************
 *	getBaseType
 * 	- returns the base type of a given character (see BASE_TYPE type)
 ***********************************************************************************/
BASE_TYPE getBaseType(char base)	{
	
  char* ptr = strchr(cannonizedBaseSymbols,base);
	
  if(ptr == NULL)		return NO_BASE;
  if(ptr<cannonizedBaseSymbols+4) 	return NUCLEOTIDE;
  if(ptr<cannonizedBaseSymbols+14)	return PARTIAL_AMBIG;
	
  return COMPLETE_AMBIG;
}
/** end of getBaseType **/



/***********************************************************************************
 *	findPattern
 * 	- finds alignment column in given pattern array
 *	- column is given as a string (no terminating char)
 *	- pattern array is a numPatterns X numSeqs  2D array
 *	- returns column id of patternArray identical to column, if exists
 *	- returns -1 of no such pattern exists
 *	- (better to use when pattern and patternArrray are canonized some way
 ***********************************************************************************/
int findPattern(const char* column, char** patternArray, int numSeqs, int numPatterns)	{
  int seq, patt;

  for(patt=0; patt<numPatterns; patt++) {
    for(seq=0; seq<numSeqs; seq++) {
      if(patternArray[patt][seq] != 
         column[seq])	break;
    }
    // pattern match is found - break out of for(patt) loop
    if(seq >= numSeqs) {
      return patt;
    }
  }
		
  // at this point, no match was found
  return -1;
}
/* end of findPattern */




/***********************************************************************************
 *	initializeBaseTransformations
 * 	- initializes global 2D array for base transformations
 *	- every row of the array is a permutation of A,C,G,T (24 total)
 *	- the row describes the permutation and its impact on the ambiguity characters
 ***********************************************************************************/
void initializeBaseTransformations() {
  int perm, base1,base2,baseMap1,baseMap2,ambig,ambigMap;
	
  if (verbose)	
    printf("Initializing base transformations....");
	
  // set up all 24 permutations
  baseTransformations[0][0] = 0;  baseTransformations[0][1] = 1;   baseTransformations[0][2] = 2;   baseTransformations[0][3] = 3; 
  baseTransformations[1][0] = 0;  baseTransformations[1][1] = 1;   baseTransformations[1][2] = 3;   baseTransformations[1][3] = 2; 
  baseTransformations[2][0] = 0;  baseTransformations[2][1] = 2;   baseTransformations[2][2] = 1;   baseTransformations[2][3] = 3; 
  baseTransformations[3][0] = 0;  baseTransformations[3][1] = 2;   baseTransformations[3][2] = 3;   baseTransformations[3][3] = 1; 
  baseTransformations[4][0] = 0;  baseTransformations[4][1] = 3;   baseTransformations[4][2] = 2;   baseTransformations[4][3] = 1; 
  baseTransformations[5][0] = 0;  baseTransformations[5][1] = 3;   baseTransformations[5][2] = 1;   baseTransformations[5][3] = 2; 
  baseTransformations[6][0] = 1;  baseTransformations[6][1] = 0;   baseTransformations[6][2] = 2;   baseTransformations[6][3] = 3; 
  baseTransformations[7][0] = 1;  baseTransformations[7][1] = 0;   baseTransformations[7][2] = 3;   baseTransformations[7][3] = 2; 
  baseTransformations[8][0] = 1;  baseTransformations[8][1] = 2;   baseTransformations[8][2] = 0;   baseTransformations[8][3] = 3; 
  baseTransformations[9][0] = 1;  baseTransformations[9][1] = 2;   baseTransformations[9][2] = 3;   baseTransformations[9][3] = 0; 
  baseTransformations[10][0] = 1; baseTransformations[10][1] = 3;  baseTransformations[10][2] = 2;  baseTransformations[10][3] = 0; 
  baseTransformations[11][0] = 1; baseTransformations[11][1] = 3;  baseTransformations[11][2] = 0;  baseTransformations[11][3] = 2; 
  baseTransformations[12][0] = 2; baseTransformations[12][1] = 1;  baseTransformations[12][2] = 0;  baseTransformations[12][3] = 3; 
  baseTransformations[13][0] = 2; baseTransformations[13][1] = 1;  baseTransformations[13][2] = 3;  baseTransformations[13][3] = 0; 
  baseTransformations[14][0] = 2; baseTransformations[14][1] = 0;  baseTransformations[14][2] = 1;  baseTransformations[14][3] = 3; 
  baseTransformations[15][0] = 2; baseTransformations[15][1] = 0;  baseTransformations[15][2] = 3;  baseTransformations[15][3] = 1; 
  baseTransformations[16][0] = 2; baseTransformations[16][1] = 3;  baseTransformations[16][2] = 0;  baseTransformations[16][3] = 1; 
  baseTransformations[17][0] = 2; baseTransformations[17][1] = 3;  baseTransformations[17][2] = 1;  baseTransformations[17][3] = 0; 
  baseTransformations[18][0] = 3; baseTransformations[18][1] = 1;  baseTransformations[18][2] = 2;  baseTransformations[18][3] = 0; 
  baseTransformations[19][0] = 3; baseTransformations[19][1] = 1;  baseTransformations[19][2] = 0;  baseTransformations[19][3] = 2; 
  baseTransformations[20][0] = 3; baseTransformations[20][1] = 0;  baseTransformations[20][2] = 2;  baseTransformations[20][3] = 1; 
  baseTransformations[21][0] = 3; baseTransformations[21][1] = 0;  baseTransformations[21][2] = 1;  baseTransformations[21][3] = 2; 
  baseTransformations[22][0] = 3; baseTransformations[22][1] = 2;  baseTransformations[22][2] = 0;  baseTransformations[22][3] = 1; 
  baseTransformations[23][0] = 3; baseTransformations[23][1] = 2;  baseTransformations[23][2] = 1;  baseTransformations[23][3] = 0;
	
  // set all ambiguities according to basic permuations
  for(perm=0; perm<24; perm++) {
    //4-way ambigs: N is N
    baseTransformations[perm][14] = 14;
		
    for(base1=0; base1<4; base1++) {
      // 3-way ambigs are determined directly by the basic permuation
      baseTransformations[perm][base1+10] = baseTransformations[perm][base1]+10;
			
      // 2-way ambigs are determined by configurations of base-pairs. Use clever trick for this
      for(base2=base1+1; base2<4; base2++) {
        ambig = 2*base1 + base2 + 3;
        if(ambig==10) ambig = 9;

        if(baseTransformations[perm][base1] < baseTransformations[perm][base2]) {
          baseMap1 = baseTransformations[perm][base1];
          baseMap2 = baseTransformations[perm][base2];
        } else {
          baseMap1 = baseTransformations[perm][base2];
          baseMap2 = baseTransformations[perm][base1];
        }

        ambigMap = 2*baseMap1 + baseMap2 + 3;
        if(ambigMap==10) ambigMap = 9;

        baseTransformations[perm][ambig] = ambigMap;
      }
    }// end of for(base1)
		
    /*		printf("Perm %2d:",perm);
            for(base1=0;base1<15;base1++) {
			printf(" %2d",baseTransformations[perm][base1]);
            }
            printf("\n");
    */		

		
  }// end of for(perm)
  if (verbose)
    printf("Done.\n");
  return;
}


/**  Canonizes an alignment column according to JC symmetry  */
int cannonizeJCpattern(const char* column, char* pattern, int numSeqs)	{
	
  int isLivingTrans[24];
	
  int trans, base, map, seq;
	
  int countTs = 0, countNonTs = 0;
	
  char* ptr;
	
	
  // initialize
  for(trans=0; trans<24; trans++) {
    isLivingTrans[trans] = 1;
  }
	
  for(seq=0; seq<numSeqs; seq++) {
    ptr = strchr(cannonizedBaseSymbols,column[seq]);
    if(ptr == NULL) {
      AlignmentGlobal.errorMessageEnd += 
        sprintf(AlignmentGlobal.errorMessageEnd,"Illegal base symbol %c.\n",column[seq]);
      return -1;
    }
    base = ptr - cannonizedBaseSymbols;
		
    // go through all living transformations and choose lowest living mapping for base
    map = 100;
    for(trans=0; trans<24; trans++) {
      if(isLivingTrans[trans] && baseTransformations[trans][base] < map) map = baseTransformations[trans][base];
    }
    if(map > 14) {
      AlignmentGlobal.errorMessageEnd += 
        sprintf(AlignmentGlobal.errorMessageEnd,"No valid mapping found for %d character of column %s (map = %d).\n",seq+1,column,map);
      return -1;
    }
		
    if(map == 0) {
      countTs++;
    } else if(map < 14) { 
      countNonTs++;
    }
    // remove all transformations not agreeing with map
    for(trans=0; trans<24; trans++) {
      if(isLivingTrans[trans] && baseTransformations[trans][base] != map) 	isLivingTrans[trans] = 0;
    }
		
    // set pattern character according to map
    pattern[seq] = cannonizedBaseSymbols[map];

		
  }// end of for(seq)

  if(countTs > 1 && countNonTs > 1) {
    return 1;
  }
  return 0;

}
/* end of cannonizeJCpattern */



/***********************************************************************************
 *	increasePatternArraySize
 * 	- increases global array allocated for site patterns
 *	- doubles it by 2
 *	- called only when adding pattern to array (in processLocusAlignment())
 ***********************************************************************************/
int increasePatternArraySize(){
  char* origPointer = AlignmentData.patternArray[0];
  int patt;
	
  AlignmentGlobal.maxNumPatterns *= 2;
	
  //	printf("increasing pattern array to size %d.\n",AlignmentGlobal.maxNumPatterns);
	
  AlignmentData.patternArray[0] = realloc(AlignmentData.patternArray[0], AlignmentGlobal.maxNumPatterns*AlignmentData.numSamples*sizeof(char));
  if(AlignmentData.patternArray[0] == NULL) {
    AlignmentGlobal.errorMessageEnd += 
      sprintf(AlignmentGlobal.errorMessageEnd,"Out Of Memory reallocating AlignmentData.patternArray[0] to %d patterns.\n",AlignmentGlobal.maxNumPatterns);
    return -1;
  }
	
  AlignmentData.patternArray = realloc(AlignmentData.patternArray, AlignmentGlobal.maxNumPatterns*sizeof(char*));
  if(AlignmentData.patternArray == NULL) {
    AlignmentGlobal.errorMessageEnd += 
      sprintf(AlignmentGlobal.errorMessageEnd,"Out Of Memory reallocating AlignmentData.patternArray to %d patterns.\n",AlignmentGlobal.maxNumPatterns);
    return -1;
  }

  if(AlignmentData.patternArray[0] != origPointer) {
    //		printf("reallocating in different place, num patterns = %d.\n",AlignmentData.numPatterns);
    for(patt = 0; patt<AlignmentData.numPatterns; patt++) {
      AlignmentData.patternArray[patt] = AlignmentData.patternArray[0] + patt*AlignmentData.numSamples;
    }
  }
	
  return 0;
}
/* end of increasePatternArraySize */



/***********************************************************************************
 *	computeHetSymmetryBreaks()
 * 	- goes through the site patterns containing het genotypes and determines 
 *		a het breaking scheme such that each diploid sequence is arbitrarily 
 *		phased at no more than a SINGLE COLUMN in each locus.
 *	- writes output in symmetryBreaks (numPatterns X numSamples) 2D array (indicating hets to arbitrarily phase per pattern):
 *	- returns 0 if all ok (-1 otherwise)
 ***********************************************************************************/
int	computeHetSymmetryBreaks(char** patternArray, int* patternCounts, int numPatterns, int numSamples, unsigned short** symmetryBreaks)	{
  int numLivePatterns;
  int patt, patt1, pattId1; // UNUSED, pattId;
  int chosenPatt, sample, sample1;
  int numBreaks; // UNUSED  breakCount
  double maxScore;

  int* livePatterns;		// array length numLivePatterns (the ids of all het patetrns elligible for arbitrary phasing)
  int* liveIndex;			// array length numPatterns (the index in livePatterns in which a live pattern resides)
  int* numLiveHets;		// array length numHetPatterns (number of hets that can be broken in pattern)
  double* pattScores;		// scores on the basis of which patterns are selected for arbitrary phasing
	
  int** liveHets;			// 2D array (numPatterns X numDiploids) listing live diploid samples per het pattern
	
	
  // allocate auxiliary memory
  pattScores = (double*)malloc(numPatterns*sizeof(double));
  if(pattScores == NULL) {
    AlignmentGlobal.errorMessageEnd += 
      sprintf(AlignmentGlobal.errorMessageEnd,"Out Of Memory allocating pattScores in computeHetSymmetryBreaks().\n");
    return -1;
  }
	
  liveHets = (int**)malloc(numPatterns*sizeof(int*));
  if(liveHets == NULL) {
    AlignmentGlobal.errorMessageEnd += 
      sprintf(AlignmentGlobal.errorMessageEnd,"Out Of Memory allocating liveHets in computeHetSymmetryBreaks().\n");
    free(pattScores);
    return -1;
  }
	
  liveHets[0] = (int*)malloc(numPatterns*(3+numSamples)*sizeof(int));
  if(liveHets[0] == NULL) {
    AlignmentGlobal.errorMessageEnd += 
      sprintf(AlignmentGlobal.errorMessageEnd,"Out Of Memory allocating liveHets[0] in computeHetSymmetryBreaks().\n");
    free(pattScores);
    free(liveHets);
    return -1;
  }
  livePatterns = liveHets[0] + numPatterns*numSamples;
  liveIndex = livePatterns + numPatterns;
  numLiveHets = liveIndex + numPatterns;

  // initialize 2D arrays
  for(patt = 0; patt<numPatterns; patt++) {
    liveHets[patt] = liveHets[0] + patt*numSamples;
    numLiveHets[patt] = 0;
    for(sample=0; sample<numSamples; sample++) {
      symmetryBreaks[patt][sample] = 0;
    }
  }
	
  // determine the het patterns
  // compute pattern score by 2^{num_hets}

  numLivePatterns = 0;
  numBreaks = 0;
  maxScore = -1.0;
  chosenPatt = -1;
  for(patt=0; patt<numPatterns; patt++) {
    //		printf("Pattern %d of %d.\n",patt+1,numPatterns);
    liveIndex[patt] = -1;
    numLiveHets[patt] = 0;
    pattScores[patt] = -1.0;
    if(patternCounts[patt] > 1) {
      //			printf("count is %d, so no symmetry breaks.\n", patternCounts[patt]);
      continue;
    }
    for(sample=0; sample<numSamples; sample++) {
			if(getBaseType(patternArray[patt][sample]) == PARTIAL_AMBIG) {
        liveHets[patt][ numLiveHets[patt] ] = sample;
        numLiveHets[patt]++;
        pattScores[patt] *= 2;
        numBreaks++;
        if(numLiveHets[patt] <= 1) {
          //					printf("het samples:",patt+1);
          livePatterns[numLivePatterns] = patt;
          liveIndex[patt] = numLivePatterns;
          numLivePatterns++;
          pattScores[patt] = 2.0;
        }
        //				printf(" %d",sample+1);
      }
    }// end of for(sample)
    //		printf(" score %lf.\n", pattScores[patt]);
    if(maxScore < pattScores[patt]) {
      chosenPatt = patt;
      maxScore = pattScores[patt];
    }
  }// end of for(patt)
	
  //	printf("%d het patterns found and %d potential breaks.\n",numLivePatterns,numBreaks);
	
	
  /*	printf("%d live patterns:",numLivePatterns);
        for(patt=0;patt<numLivePatterns; patt++) {
		printf(" %d (score = %.1lf, hets:",livePatterns[patt]+1,pattScores[livePatterns[patt]]);
		for(sample=0; sample<numLiveHets[livePatterns[patt]]; sample++) {
        printf(" %d",liveHets[livePatterns[patt]][sample]+1);
		}
		printf(")");
        }
        printf("\n");
  */


  //	printf("Finding symmetry breaks...\n");
// UNUSED  breakCount = 0;
  // break symmetries until maxScore is negative
  while(maxScore > 0.0) {
    // choose sample for arbitrary phasing
    numLiveHets[chosenPatt]--;
    sample = liveHets[chosenPatt][ numLiveHets[chosenPatt] ];
    symmetryBreaks[chosenPatt][sample] = 1;
    //		printf("Breaking symmetry of pattern %d in sample %d, %d hets, score %.1lf.\n",chosenPatt+1,sample+1,numLiveHets[chosenPatt]+1, maxScore);
		
    // adjust score of chosen patt
    if(numLiveHets[chosenPatt] <= 0) {
      pattScores[chosenPatt] = -1.0;
    } else {
      pattScores[chosenPatt] /= 2;
    }
		
    // go over all other patterns to search for new choice
    maxScore = pattScores[chosenPatt];
// UNUSED    pattId = chosenPatt;
    for(patt1=0; patt1<numLivePatterns; ) {
      pattId1 = livePatterns[patt1];
      //			printf(" live pattern %d (id=%d)\n",patt1+1,pattId1+1);
      if(liveIndex[pattId1] < 0 || liveIndex[pattId1] >= numLivePatterns || livePatterns[ liveIndex[pattId1] ] != pattId1) {
        AlignmentGlobal.errorMessageEnd += 
          sprintf(AlignmentGlobal.errorMessageEnd,"Error in keeping live index %d of pattern %d .\n", liveIndex[pattId1]+1, pattId1+1);
        free(pattScores);
        free(liveHets[0]);
        free(liveHets);
        return -1;
      }
      for(sample1=0; sample1<numLiveHets[pattId1]; sample1++) {
        if(liveHets[pattId1][sample1] == sample) {
          numLiveHets[pattId1]--;
          liveHets[pattId1][sample1] = liveHets[pattId1][ numLiveHets[pattId1] ];
          //					printf("patterns %d and %d share sample %d.\n",pattId1+1,pattId+1,sample);
          break;
        }
      }// end of for(sample1)
      if(numLiveHets[pattId1] > 0) {
        patt1++;
      } else {
        // mark pattern as ineligible
				
        //				printf("pattern %d became ineligible, switching with pattern %d (live index=%d).\n",pattId1+1, livePatterns[numLivePatterns-1]+1, numLivePatterns);
        numLivePatterns--;
        livePatterns[ liveIndex[pattId1] ] = livePatterns[numLivePatterns];
        liveIndex[ livePatterns[ liveIndex[pattId1] ] ] = liveIndex[pattId1];
        liveIndex[pattId1] = -1;
        pattScores[pattId1] = -1.0;
      }
			
      if(maxScore < pattScores[pattId1]) {
        maxScore = pattScores[pattId1];
        chosenPatt = pattId1;
      }
    }// end of for(patt1)
				
		
			
	
    /*
      printf("%d live patterns:",numLivePatterns);
      for(patt=0;patt<numLivePatterns; patt++) {
      printf(" %d (score = %.1lf, hets:",livePatterns[patt]+1,pattScores[livePatterns[patt]]);
      for(sample=0; sample<numLiveHets[livePatterns[patt]]; sample++) {
      printf(" %d",liveHets[livePatterns[patt]][sample]+1);
      }
      printf(")");
      }
      printf("\n");
      **/
  }// end of while(maxScore > 0.0)
	
	
	
  // free auxilliary allocated memory
  free(pattScores);
  free(liveHets[0]);
  free(liveHets);
	
  return 0;
}
/* end of computeHetSymmetryBreaks */



/***********************************************************************************
 *	computeHetSymmetryBreaks_allLoci()
 *	- this version is not currently used. It considers together the alignment of all loci
 * 	- goes through the site patterns containing het genotypes and determines 
 *		a het breaking scheme such that each diploid sequence is arbitrarily 
 *		phased at no more than a SINGLE COLUMN in each locus.
 *	- writes output in symmetryBreaks (numPatterns X numSamples) 2D array (indicating hets to arbitrarily phase per pattern):
 *	- also updates AlginemntData.isDiploid
 *	- returns 0 if all ok (-1 otherwise)
 ***********************************************************************************/
int	computeHetSymmetryBreaks_allLoci(char** patternArray, int numPatterns, int numSamples, unsigned short** symmetryBreaks)	{
  int numLoci = AlignmentData.numLoci;
  int numLivePatterns;
  int patt, patt1, pattId, pattId1;
  int chosenPatt, sample, sample1, locus;
  int breakCount, numPhases, numBreaks;
  double maxScore;

  int* livePatterns;		// array length numLivePatterns (the ids of all het patetrns elligible for arbitrary phasing)
  int* liveIndex;			// array length numPatterns (the index in livePatterns in which a live pattern resides)
  int* numLiveHets;		// array length numHetPatterns (number of hets that can be broken in pattern)
  int* numLociInPatt;		// array length numHetPatterns (number of loci in which het pattern appears)
  double* pattScores;		// scores on the basis of which patterns are selected for arbitrary phasing
	
  int** liveHets;			// 2D array (numPatterns X numDiploids) listing live diploid samples per het pattern
  int** pattern2locus;	// 2D array (numPatterns X numLoci) listing loci containing this site pattern
  unsigned short** patternCoexist;	// 2D array (numPatterns X numPatterns) indicating whether 2 patterns exist in same locus
	
	
  // allocate auxiliary memory
  patternCoexist = (unsigned short**)malloc(numPatterns*sizeof(unsigned short*));
  if(patternCoexist == NULL) {
    AlignmentGlobal.errorMessageEnd += 
      sprintf(AlignmentGlobal.errorMessageEnd,"Out Of Memory allocating patternCoexist in computeHetSymmetryBreaks().\n");
    return -1;
  }			
	
  patternCoexist[0] = (unsigned short*)malloc(numPatterns*numPatterns*sizeof(unsigned short));
  if(patternCoexist[0] == NULL) {
    AlignmentGlobal.errorMessageEnd += 
      sprintf(AlignmentGlobal.errorMessageEnd,"Out Of Memory allocating patternCoexist[0] in computeHetSymmetryBreaks().\n");
    free(patternCoexist);
    return -1;
  }
	
  pattScores = (double*)malloc(numPatterns*sizeof(double));
  if(pattScores == NULL) {
    AlignmentGlobal.errorMessageEnd += 
      sprintf(AlignmentGlobal.errorMessageEnd,"Out Of Memory allocating pattScores in computeHetSymmetryBreaks().\n");
    free(patternCoexist[0]);
    free(patternCoexist);
    return -1;
  }
	
  liveHets = (int**)malloc(2*numPatterns*sizeof(int*));
  if(liveHets == NULL) {
    AlignmentGlobal.errorMessageEnd += 
      sprintf(AlignmentGlobal.errorMessageEnd,"Out Of Memory allocating liveHets in computeHetSymmetryBreaks().\n");
    free(patternCoexist[0]);
    free(patternCoexist);
    free(pattScores);
    return -1;
  }
  pattern2locus = liveHets + numPatterns;
	
  liveHets[0] = (int*)malloc(numPatterns*(4+numSamples+numLoci)*sizeof(int));
  if(liveHets[0] == NULL) {
    AlignmentGlobal.errorMessageEnd += 
      sprintf(AlignmentGlobal.errorMessageEnd,"Out Of Memory allocating liveHets[0] in computeHetSymmetryBreaks().\n");
    free(patternCoexist[0]);
    free(patternCoexist);
    free(pattScores);
    free(liveHets);
    return -1;
  }
  pattern2locus[0] = liveHets[0] + numPatterns*numSamples;
  livePatterns = pattern2locus[0] + numPatterns*numLoci;
  liveIndex = livePatterns + numPatterns;
  numLiveHets = liveIndex + numPatterns;
  numLociInPatt = numLiveHets + numPatterns;

  // initialize 2D arrays
  for(patt = 0; patt<numPatterns; patt++) {
    patternCoexist[patt] = patternCoexist[0] + patt*numPatterns;
    liveHets[patt] = liveHets[0] + patt*numSamples;
    pattern2locus[patt] = pattern2locus[0] + patt*numLoci;
		
    numLiveHets[patt] = 0;
    numLociInPatt[patt] = 0;
    for(sample=0; sample<numSamples; sample++) {
      symmetryBreaks[patt][sample] = 0;
    }
  }
	
  // determine the diploid samples and the het patterns
  // initialize pattern score by 2^{num_hets}

  printf("Determining het patterns...\n");

  numLivePatterns = 0;
  numPhases = 0;
  numBreaks = 0;
  for(patt=0; patt<numPatterns; patt++) {
    liveIndex[patt] = -1;
    numLiveHets[patt] = 0;
    numLociInPatt[patt] = 0;
    pattScores[patt] = -1.0;
    for(patt1=0; patt1<numPatterns; patt1++) {
      patternCoexist[patt][patt1] = 0;
    }
    for(sample=0; sample<numSamples; sample++) {
      if(getBaseType(AlignmentData.patternArray[patt][sample]) == PARTIAL_AMBIG) {
        liveHets[patt][ numLiveHets[patt] ] = sample;
        numLiveHets[patt]++;
        pattScores[patt] *= 2;
        numBreaks++;
        if(numLiveHets[patt] <= 1) {
          //					printf("New het pattern found: %d, with het samples:",patt+1);
          livePatterns[numLivePatterns] = patt;
          liveIndex[patt] = numLivePatterns;
          numLivePatterns++;
          pattScores[patt] = 2.0;
        }
        //				printf(" %d",sample+1);
      }
    }// end of for(sample)
    if(pattScores[patt] > 0) {
      numPhases += (int)pattScores[patt];
    } else {
      numPhases ++;
    }
  }// end of for(patt)
	
  //	printf("%d het patterns found, with total number of phases %d, and %d potential breaks.\n",numLivePatterns,numPhases,numBreaks);
	
  //	printf("Determining candidates for symmetry breaking (appearing in < 2 columns per locus), and pattern coexistence...\n");
  // check het patterns eligible for symmetry breaking
  // consider only patterns that appear once in each locus
  // also compute patterns which coexist in the same locus
  for(locus=0; locus<numLoci; locus++) {
    for(patt=0; patt<AlignmentData.locusProfiles[locus].numPatterns; patt++) {
      pattId = AlignmentData.locusProfiles[locus].patternIds[patt];

			
      if(pattScores[pattId] < 0.0) {
        //				printf("Pattern %d in locus %d already ineligible.\n",pattId+1,locus+1);
        continue;
      }
			
      if(AlignmentData.locusProfiles[locus].patternCounts[patt] <= 1) {
        //				printf("Pattern %d found with %d occurrences in locus %d.\n",pattId+1,AlignmentData.locusProfiles[locus].patternCounts[patt],locus+1);
        pattern2locus[pattId][ numLociInPatt[pattId] ] = locus;
        numLociInPatt[pattId]++;
			
        for(patt1=0; patt1<patt; patt1++) {
          pattId1 = AlignmentData.locusProfiles[locus].patternIds[patt1];
          if(pattScores[pattId] > 0) {
            patternCoexist[pattId][pattId1] = patternCoexist[pattId1][pattId] = 1;
          }
        }// end of fo(patt1)
      } else {
        // mark pattern as ineligible
        pattScores[pattId] = -1.0;
      }
    }// end of for(patt)
  }// end of for(locus)
	
  printf("Done.\n");
	
  /*	printf("%d live patterns:",numLivePatterns);
        for(patt=0;patt<numLivePatterns; patt++) {
		printf(" %d (score = %.1lf, hets:",livePatterns[patt]+1,pattScores[livePatterns[patt]]);
		for(sample=0; sample<numLiveHets[livePatterns[patt]]; sample++) {
        printf(" %d",liveHets[livePatterns[patt]][sample]+1);
		}
		printf(")");
        }
        printf("\n");
  */
  // set pattern score (of eligible het patterns) to 2^{num_hets} / num_loci
  // compute coexisting patterns
  // and find pattern with maximal score


  printf("Adjusting scores and finding initial eligible patterns.\n");

  maxScore = -1.0;
  chosenPatt = -1;
  for(patt=0; patt<numLivePatterns; patt++) {
    pattId = livePatterns[patt];
    if(pattScores[pattId] > 0.0) {
      pattScores[pattId] /= numLociInPatt[pattId];
    } else {
      if(liveIndex[pattId] < 0 || liveIndex[pattId] >= numLivePatterns || livePatterns[ liveIndex[pattId] ] != pattId) {
        AlignmentGlobal.errorMessageEnd += 
          sprintf(AlignmentGlobal.errorMessageEnd,"Error in keeping indices for live patterns.\n");
        free(patternCoexist[0]);
        free(patternCoexist);
        free(pattScores);
        free(liveHets[0]);
        free(liveHets);
        return -1;
      }
			
      numLivePatterns--;
      /*				printf("Pattern %d (live id %d) found inelligible due to %d occurrences in locus %d (score %lf).\n",
                        pattId+1,liveIndex[pattId]+1, AlignmentData.locusProfiles[locus].patternCounts[patt],locus+1,pattScores[pattId]);
                        printf("Switching with pattern %d (live id %d).\n",livePatterns[numLivePatterns]+1, numLivePatterns+1);
      */			
      livePatterns[ liveIndex[pattId] ] = livePatterns[numLivePatterns];
      liveIndex[ livePatterns[ liveIndex[pattId] ] ] = liveIndex[pattId];
      liveIndex[pattId] = -1;
      numLociInPatt[pattId] = 0;
      numBreaks -= numLiveHets[pattId];
      numLiveHets[pattId] = 0;
    }
    if(maxScore < pattScores[pattId]) {
      chosenPatt = pattId;
      maxScore = pattScores[pattId];
    }
  }// end of for(patt)

  printf("%d eligible patterns found, with %d total breaks.\n",numLivePatterns,numBreaks);

  printf("Finding symmetry breaks...\n");
  breakCount = 0;
  // break symmetries until maxScore is negative
  //	maxScore = -1.0;
  while(maxScore > 0.0) {
    // choose sample for arbitrary phasing
    numLiveHets[chosenPatt]--;
    sample = liveHets[chosenPatt][ numLiveHets[chosenPatt] ];
    symmetryBreaks[chosenPatt][sample] = 1;
    //		printf("Breaking symmetry of pattern %d in sample %d, %d hets, score %.1lf.\n",chosenPatt+1,sample+1,numLiveHets[chosenPatt]+1, maxScore);
		
    // adjust score of chosen patt
    if(numLiveHets[chosenPatt] <= 0) {
      pattScores[chosenPatt] = -1.0;
    } else {
      pattScores[chosenPatt] /= 2;
    }
		
    // go over all other patterns to search for new choice
    // remove sample from liveHets of all patterns coexisting with this one in same locus
    maxScore = pattScores[chosenPatt];
    pattId = chosenPatt;
    for(patt1=0; patt1<numLivePatterns; ) {
      pattId1 = livePatterns[patt1];
      //			printf(" live pattern %d (id=%d)\n",patt1+1,pattId1+1);
      if(liveIndex[pattId1] < 0 || liveIndex[pattId1] >= numLivePatterns || livePatterns[ liveIndex[pattId1] ] != pattId1) {
        AlignmentGlobal.errorMessageEnd += 
          sprintf(AlignmentGlobal.errorMessageEnd,"Error in keeping live index %d of pattern %d .\n", liveIndex[pattId1]+1, pattId1+1);
        free(patternCoexist[0]);
        free(patternCoexist);
        free(pattScores);
        free(liveHets[0]);
        free(liveHets);
        return -1;
      }
      if(patternCoexist[pattId][pattId1]) {
        for(sample1=0; sample1<numLiveHets[pattId1]; sample1++) {
          if(liveHets[pattId1][sample1] == sample) {
            numLiveHets[pattId1]--;
            liveHets[pattId1][sample1] = liveHets[pattId1][ numLiveHets[pattId1] ];
            //						printf("pattern %d coexists with pattern %d, and they share sample %d.\n",pattId1+1,pattId+1,sample);
            break;
          }
        }// end of for(sample1)
      }
      if(numLiveHets[pattId1] > 0) {
        patt1++;
      } else {
        // mark pattern as ineligible
				
        //				printf("pattern %d became ineligible, switching with pattern %d (live index=%d).\n",pattId1+1, livePatterns[numLivePatterns-1]+1, numLivePatterns);
        numLivePatterns--;
        livePatterns[ liveIndex[pattId1] ] = livePatterns[numLivePatterns];
        liveIndex[ livePatterns[ liveIndex[pattId1] ] ] = liveIndex[pattId1];
        liveIndex[pattId1] = -1;
        pattScores[pattId1] = -1.0;
        numLociInPatt[pattId1] = 0;
      }
			
      if(maxScore < pattScores[pattId1]) {
        maxScore = pattScores[pattId1];
        chosenPatt = pattId1;
      }
    }// end of for(pattId1)
				
		
    // eliminate patterns that coexist with this one
			
	
    /*
      printf("%d live patterns:",numLivePatterns);
      for(patt=0;patt<numLivePatterns; patt++) {
      printf(" %d (score = %.1lf, hets:",livePatterns[patt]+1,pattScores[livePatterns[patt]]);
      for(sample=0; sample<numLiveHets[livePatterns[patt]]; sample++) {
      printf(" %d",liveHets[livePatterns[patt]][sample]+1);
      }
      printf(")");
      }
      printf("\n");
      **/
    breakCount++;
    if(breakCount%100 == 0) {
      printf(".");
      if(breakCount%1000 == 0) {
        printf(" ");
        if(breakCount%10000 == 0) {
          printf("\n");
        }
      }
      fflush(stdout);
    }
  }// end of while(maxScore > 0.0)
	
	
	
  // free auxiliary allocated memory
  free(patternCoexist[0]);
  free(patternCoexist);
  free(pattScores);
  free(liveHets[0]);
  free(liveHets);
	
  return 0;
}
/* end of computeHetSymmetryBreaks_allLoci */




/***********************************************************************************
 *	getAllPhases
 * 	- returns all phased versions of a given alignment column
 *	- the input column is given such that each sample is a haploid
 *	- perturbSample is a binary array of length numHaploids indicating which
 *		haploid pairs to perturb. If a given entry is 1, then perturbs that haploid with the previous one.
 *	- writes on top of perturbSample to go over all states of perturbations, but returns it to original state
 *	- writes phased patterns onto phasedColumns array (pre-allocated)
 *	- returns total number of phases
 ***********************************************************************************/
int getAllPhases(char* column, int numHaploids, unsigned short* perturbSample, char** phasedColumns)	{
  int haploid, numPhases, phase;
  unsigned short flipPhase;
	
  // copy initial phase
  for(haploid=0, numPhases = 1; haploid<numHaploids; haploid++) {
    phasedColumns[0][haploid] = column[haploid];
    if(perturbSample[haploid] > 0) numPhases *=2;
  }
	
  // non-phased pairs remain 0, phased pairs alternates between 1 and 2
  for(phase=1; phase<numPhases; phase++) {

    // compute next phasing
    flipPhase = 1;
    for(haploid=0; haploid<numHaploids; haploid++) {
      phasedColumns[phase][haploid] = phasedColumns[phase-1][haploid];
      if(flipPhase && perturbSample[haploid] > 0) {
        if(haploid <= 0 || perturbSample[haploid-1] != 0) {
          AlignmentGlobal.errorMessageEnd += 
            sprintf(AlignmentGlobal.errorMessageEnd,"Illegal setting for perturbSample array in haploid pair (%d,%d).\n",haploid,haploid+1);
          return -1;
        }
        phasedColumns[phase][haploid] = phasedColumns[phase-1][haploid-1];
        phasedColumns[phase][haploid-1] = phasedColumns[phase-1][haploid];
        if(perturbSample[haploid] == 1) {
          perturbSample[haploid] = 2;
          // last phase to flip
          flipPhase = 0;
        } else {
          perturbSample[haploid] = 1;
        }
      }
    }// end of for(haploid)

    // if all phases were flipped back to 1, break out of loop
    if(flipPhase) {
      AlignmentGlobal.errorMessageEnd += 
        sprintf(AlignmentGlobal.errorMessageEnd,"All phases were flipped in getAllPhases in phase %d.\n",phase+1);
      return -1;
    }
  }// end of for(numPhases)
	
  return numPhases;

}
/* end of getAllPhases */



/***********************************************************************************
 *	translateAmbiguity
 * 	- translates ambiguity character to a pair of bases in given string pointer
 *	- if character is not a nucleotide symbol or 2-wise ambiguity symbol, writes 2 'N's
 *	- returns base type of character
 ***********************************************************************************/
BASE_TYPE translateAmbiguity(char ch, char* outp) {

  BASE_TYPE baseType = getBaseType(ch);

  switch(ch) {
  case('Y'):
    outp[0] = 'T';
    outp[1] = 'C';
    break;
  case('K'):
    outp[0] = 'T';
    outp[1] = 'G';
    break;
  case('W'):
    outp[0] = 'T';
    outp[1] = 'A';
    break;
  case('S'):
    outp[0] = 'C';
    outp[1] = 'G';
    break;
  case('M'):
    outp[0] = 'A';
    outp[1] = 'C';
    break;
  case('R'):
    outp[0] = 'A';
    outp[1] = 'G';
    break;
  case('T'):
  case('C'):
  case('A'):
  case('G'):
    outp[0] = outp[1] = ch;
    break;
  default:
    outp[0] = outp[1] = 'N';
    break;
  }
	
  return baseType;
}
/* end of translateAmbiguity */



/***********************************************************************************
 *	translateToAmbiguity
 * 	- translates a base pair into an ambiguity character
 *	- if characters are not a nucleotide symbol or 2-wise ambiguity symbol, writes 2 'N's
 *	- returns ambiguity character
 ***********************************************************************************/
char translateToAmbiguity(char* basePair) {
  char ch;

  if(basePair[0] == basePair[1]) {
    ch = basePair[0];
  } else {
    if(basePair[0] > basePair[1]) {
      ch = basePair[0];
      basePair[0] = basePair[1];
      basePair[1] = ch;
    }
    switch(basePair[0]) {
    case('A'):
      switch(basePair[1]) {
      case('C'):
        ch = 'M';
        break;
      case('G'):
        ch = 'R';
        break;
      case('T'):
        ch = 'W';
        break;
      default:
        ch = 'X';
        break;
      }
      break;
    case('C'):
      switch(basePair[1]) {
      case('G'):
        ch = 'S';
        break;
      case('T'):
        ch = 'Y';
        break;
      default:
        ch = 'X';
        break;
      }
      break;
    case('G'):
      switch(basePair[1]) {
      case('T'):
        ch = 'K';
        break;
      default:
        ch = 'X';
        break;
      }
      break;
    default:
      ch = 'X';
      break;
    }
  }
	
  return ch;
}
/* end of translateToAmbiguity */



/***********************************************************************************
 *	twoSiteFourGameteTest()
 * 	- performs the 4-gamete test between two potential sites
 *	- returns 0, if passes test
 *	- returns 1, if fails test without doubt
 *	- returns 2, if passes test due to greedy phasing
 ***********************************************************************************/
int	twoSiteFourGameteTest(int patternId1, int patternId2) {
  int numSamples = AlignmentData.numSamples;

  int res;
  int i, j, config, sample, numConfigs, maxI;
  BASE_TYPE baseType1;
  BASE_TYPE baseType2;
  char* pattern1 = AlignmentData.patternArray[patternId1];
  char* pattern2 = AlignmentData.patternArray[patternId2];	
  char ambig1[2];
  char ambig2[2];
  int existsPairConfig[4]; // existence of pair configurations
  int optionalConfigs[4];
  int numNew1, numNew2;
	
  for(config=0; config<4; config++) {
    existsPairConfig[config] = 0;
  }
	
  numConfigs = 0;

  // first pass - do only non-ambiguous pairs
  for(sample=0; sample<numSamples; sample++) {
    baseType1 = getBaseType(pattern1[sample]);
    baseType2 = getBaseType(pattern2[sample]);
		
    if(baseType1 == COMPLETE_AMBIG || baseType2 == COMPLETE_AMBIG) {
      continue;
    }
    translateAmbiguity(pattern1[sample], ambig1);
    translateAmbiguity(pattern2[sample], ambig2);
		
    if(baseType1 == NUCLEOTIDE && baseType2 == NUCLEOTIDE) {
      // check only first character
      maxI = 0;
    } 
    else if(baseType1 == NUCLEOTIDE || baseType2 == NUCLEOTIDE) {
      maxI = 1;
    }
    // both are ambiguity characters
    else {
      continue;
    } 
		
    for(i=0; i<=maxI; i++) {
      if(ambig1[i] == ambig2[i]) {
        if(ambig1[i] == 'T') {
          config = 0;
        } else {
          config = 3;
        }
      } else {
        if(ambig1[i] == 'T') {
          config = 1;
        } else {
          config = 2;
        }
      }				
      if(existsPairConfig[config] == 0) {	
        //				printf(" config %d",config);
        existsPairConfig[config] = 1;
        numConfigs++;
      }
    }
  }// end of for(sample)
	
	
  if(numConfigs == 4) {
    return 1;
  }
	
  if(numConfigs == 3) {
    res = 1; // if failed, then without a doubt
  } else {
    res = 2;
  }
  //	printf("\n");
  // second pass - do only ambiguous pairs
  for(sample=0; sample<numSamples; sample++) {
    baseType1 = getBaseType(pattern1[sample]);
    baseType2 = getBaseType(pattern2[sample]);
		
    if(!(baseType1 == PARTIAL_AMBIG && baseType2 == PARTIAL_AMBIG)) {
      continue;
    }
    translateAmbiguity(pattern1[sample], ambig1);
    translateAmbiguity(pattern2[sample], ambig2);
		
    // figure out all pair configurations
    for(i=0; i<2; i++) {
      for(j=0; j<2; j++) {
        if(ambig1[i] == ambig2[j]) {
          if(ambig1[i] == 'T') {
            config = 0;
          } else {
            config = 3;
          }
        } else {
          if(ambig1[i] == 'T') {
            config = 1;
          } else {
            config = 2;
          }
        }
        optionalConfigs[2*i+j] = config;
      }
    }
    //		printf("ambiguous pair %d, with configs (%d,%d) or (%d,%d).\n",sample,optionalConfigs[0],optionalConfigs[3],optionalConfigs[1],optionalConfigs[2]);
		
    numNew1 = (existsPairConfig[ optionalConfigs[0] ] == 0) + (existsPairConfig[ optionalConfigs[3] ] == 0);
    numNew2 = (existsPairConfig[ optionalConfigs[1] ] == 0) + (existsPairConfig[ optionalConfigs[2] ] == 0);
		
    if(numNew2 < numNew1) {
      //			printf("choosing pair (%d,%d).\n",optionalConfigs[1],optionalConfigs[2]);
      optionalConfigs[0] = optionalConfigs[1];
      optionalConfigs[1] = optionalConfigs[2];
    } else {
      //			printf("choosing pair (%d,%d).\n",optionalConfigs[0],optionalConfigs[3]);
      //optionalConfigs[0] = optionalConfigs[0];
      optionalConfigs[1] = optionalConfigs[3];
    }
			
    for(i=0; i<2; i++) {		
      if(existsPairConfig[ optionalConfigs[i] ] == 0) {	
        existsPairConfig[ optionalConfigs[i] ] = 1;
        numConfigs++;
      }
    }
		
  }// end of for(sample)
	
	
	
	
  if (numConfigs < 4) {
    return 0;
  } else {
    return res;
  }
}
/* end of twoSiteFourGameteTest */



/***************************************************************************************************************/
/******                                        END OF FILE                                                ******/
/***************************************************************************************************************/
