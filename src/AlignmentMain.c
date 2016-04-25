/** 
   \file AlignmentMain.c
    Functions for performing 4 gamete test and analyzing patterns in samples.
   
*/

/******************************************************************************************************/
/******                                      INCLUDES                                            ******/
/******************************************************************************************************/

#include "AlignmentProcessor.h"
#include "utils.h"


/***************************************************************************************************************/
/******                                  GLOBAL DATA STRUCTURES                                           ******/
/***************************************************************************************************************/



/******************************************************************************************************/
/******                                INTERNAL FUNCTION DECLARATION                             ******/
/******************************************************************************************************/



/******************************************************************************************************/
/******                                         M A I N                                          ******/
/******************************************************************************************************/

/** Prints to stdout, usage text for programs Haploid2Diploid, AnalyzePatterns, 4GameteTest
    @param Name of the program (i.e. Haploid2Diploid)
*/
void usage(char* progName) {
  printf("Usage: %s <seq-file> <out-file> <num-samples> [-l num-loci] <diploid-list>\n",progName);
  printf("       - seq-file: sequence file to read\n");
  printf("       - out-file: output file in which to write diploid alignments\n");
  printf("       - num-samples: number of haploid samples in each alignment\n");
  printf("       - diploid-list: list of haploid indices (1-based) indicating diploid pairs (i indicates a pair (i,i+1))\n");
  printf("       - -l num-loci: max number of loci to transform (optional)\n");

  return;
}



/** Main function for executable 'AnalyzePatterns '
    @param argc Number of command line arguments
    @param argv Array of command line arguments
    @return return value of program
    @note This executable is not built as part of GPhoCS
*/
int main_analyze_patterns(int argc, char*argv[]) {
	
  char sampleFile[STRING_LENGTH];
  char seqFile[STRING_LENGTH];
  char** seqNames;
  int numSamples, sample, numInformative, numLoci = -1;
  int res;
  FILE *fsample;
	
  if(argc<= 2) {
    printf("Usage: AlignmentProcessor <samples-file> <seq-file> [num-loci]\n");
    return 1;
  }
	
  strcpy(sampleFile,argv[1]);
  strcpy(seqFile,argv[2]);
  if(argc == 4) {
    res=sscanf(argv[3],"%d",&numLoci);
    if(res < 1) {
      fprintf(stderr, "Error: 3rd option (%s) cannot act as number of loci.\n",argv[3]);
      return 1;
    }
  }
  fsample=(FILE*)fopen(sampleFile,"r");
  if(fsample == NULL) {
    fprintf(stderr, "Error: Could not open sample file %s.\n", sampleFile);
    return 1;
  }
  res = fscanf(fsample,"%d",&numSamples);
  seqNames = (char**)malloc(numSamples*sizeof(char*));
  if(seqNames == NULL) {
    fprintf(stderr, "Error: Out Of Memory seqNames array in main().\n");
    fclose(fsample);
    return 1;
  }
  seqNames[0] = (char*)malloc(NAME_LENGTH*numSamples*sizeof(char));
  if(seqNames == NULL) {
    fprintf(stderr, "Error: Out Of Memory seqNames space in main().\n");
    free(seqNames);
    fclose(fsample);
    return 1;
  }

  for(sample=0; sample<numSamples; sample++) {
    seqNames[sample] = seqNames[0] + sample*NAME_LENGTH;
    res=fscanf(fsample,"%s",seqNames[sample]);
    if(res < 1) {
      fprintf(stderr, "Error: Could not read sample number %d from file.\n",sample+1);
      free(seqNames[0]);
      free(seqNames);
      fclose(fsample);
      return 1;
    }
  }
	
  fclose(fsample);

  //	readSeqFile(seqFile,numSamples,seqNames,numLoci);
  readSeqFile(seqFile,numSamples,NULL,numLoci);

  free(seqNames[0]);
  free(seqNames);
	
  printf("Total %d patterns.\n",AlignmentData.numPatterns);
	
	
  res = processHetPatterns(/*break the symmetries*/ 1);
  if(res<0) {
    fprintf(stderr, "Error: problem occurred while processing het patterns.\n");
    printAlignmentError();
    freeAlignmentData();
    return -1;
  }
	
  res = getPhasedPatternTypes();
  if(res<0) {
    fprintf(stderr, "Error: problem occurred while computing pattern types.\n");
    printAlignmentError();
    freeAlignmentData();
    return -1;
  }
	
  freeAlignmentData();
  return 0;

}
// end of main 



/** Main function for executable '4GameteTest'
    @param argc Number of command line arguments
    @param argv Array of command line arguments
    @return return value of program (1=error 0=OK)
    @note This executable is not bult as part of GPhoCS
*/
int main_4gam_test(int argc, char*argv[]) {
	
  char sampleFile[STRING_LENGTH];
  char seqFile[STRING_LENGTH];
  char** seqNames;
  int numSamples, sample, numInformative, numLoci = -1;
  int res;
  FILE *fsample;
	
  if(argc<= 2) {
    printf("Usage: AlignmentProcessor <samples-file> <seq-file> [num-loci]\n");
    return 1;
  }
	
  strcpy(sampleFile,argv[1]);
  strcpy(seqFile,argv[2]);
  if(argc == 4) {
    res=sscanf(argv[3],"%d",&numLoci);
    if(res < 1) {
      fprintf(stderr, "Error: 3rd option (%s) cannot act as number of loci.\n",argv[3]);
      return 1;
    }
  }
  fsample=(FILE*)fopen(sampleFile,"r");
  if(fsample == NULL) {
    fprintf(stderr, "Error: Could not open sample file %s.\n", sampleFile);
    return 1;
  }
  res = fscanf(fsample,"%d",&numSamples);
  seqNames = (char**)malloc(numSamples*sizeof(char*));
  if(seqNames == NULL) {
    fprintf(stderr, "Error: Out Of Memory seqNames array in main().\n");
    fclose(fsample);
    return 1;
  }
  seqNames[0] = (char*)malloc(NAME_LENGTH*numSamples*sizeof(char));
  if(seqNames == NULL) {
    fprintf(stderr, "Error: Out Of Memory seqNames space in main().\n");
    free(seqNames);
    fclose(fsample);
    return 1;
  }

  for(sample=0; sample<numSamples; sample++) {
    seqNames[sample] = seqNames[0] + sample*NAME_LENGTH;
    res=fscanf(fsample,"%s",seqNames[sample]);
    if(res < 1) {
      fprintf(stderr, "Error: Could not read sample number %d from file.\n",sample+1);
      free(seqNames[0]);
      free(seqNames);
      fclose(fsample);
      return 1;
    }
  }
	
  fclose(fsample);

  readSeqFile(seqFile,numSamples,seqNames,numLoci);
  //	readSeqFile(seqFile,numSamples,NULL,numLoci);

  free(seqNames[0]);
  free(seqNames);
	
  printf("Performing 4 gamete test.\n");
  res = fourGameteTest();
  if(res<0) {
    fprintf(stderr, "Error: problem occurred while performing 4 gamete test.\n");
    printAlignmentError();
    freeAlignmentData();
    return -1;
  }
  printf("Done.\n");
	
  freeAlignmentData();
  return 0;

}
// end of main 



/***************************************************************************************************************/
/******                                        END OF FILE                                                ******/
/***************************************************************************************************************/
