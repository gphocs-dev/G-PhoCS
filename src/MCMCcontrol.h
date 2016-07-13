#ifndef MCMC_CONTROL_FILE_H
#define MCMC_CONTROL_FILE_H
/** 
	\file MCMCcontrol.h
    Read and processes and holds control information for MCMC

	This header file describes the interface functions for reading and 
	processing a control file.
	
*/

# include "PopulationTree.h"

/***************************************************************************************************************/
/******                                              DATA TYPES                                           ******/
/***************************************************************************************************************/




/*********
 * UpdateStats used to hold information on MCMC update steps
 *		(finetunes, acceptance rates, etc.)
 *********/
typedef struct UPDATE_STAST{
  double coalTime;
  double SPR;
  double migTime;
  double theta;
  double migRate;
  double *taus;
  double locusRate;
  double admix;
  double mixing;
} UpdateStats;



/***************************************************************************************************************/
/******                                  GLOBAL DATA STRUCTURES                                           ******/
/***************************************************************************************************************/



/*********
 * i/o setup
 *********/
struct IO_SETUP {
	char seedFileName[NAME_LENGTH];			// name of random seed file
	char seqFileName[NAME_LENGTH];			// name of sequence file
	// char debugFileName[NAME_LENGTH];		// name of debug trace file
	char admixFileName[NAME_LENGTH];		// name of admixture trace file
	char rateFileName[NAME_LENGTH];			// name of locus-rate file
	char traceFileName[NAME_LENGTH];		// name of trace file (for MCMC trace output)
	char flatStatsFileName[NAME_LENGTH];	// name of coalescent stats file (for MCMC model evaluation)
	char cladeStatsFileName[NAME_LENGTH];	// name of clade coalescent stats file
	int samplesPerLog;						// number of samples for which to generate a log summary in stdout
	int logsPerLine;						// number of sample logs per log line
	
	FILE*	traceFile;						// trace file
	FILE*	debugFile;						// debugging file
	FILE*	admixFile;						// admixture stats file
	FILE*	flatStatsFile;					// coalescent stats file
	FILE*	cladeStatsFile;					// clade stats file
	FILE**	nodeStatsFile;					// coalescent stats files for nodes (one per pop)
}  ioSetup;


/*********
 * mcmc setup
 *********/
struct MCMC_SETUP {
	// general info
	int numParameters;				// number of model parameters
	unsigned short useData;			// flag which indicated whether to use sequence data or to sample from prior 
	int randomSeed;					// random seed used
	
  // sampling info
	int numSamples;					// number of sampling iterations to perform
	int burnin;						// number of iterations for burnin
	int sampleSkip;					// number of samples to skip between each recorded sample
	int startMig;						// number of generations to skip before starting to sample migrations
	int genetreeSamples;				// number of gene tree updates per each population parameter update
	
	unsigned short allowAdmixture;	// flag which is turned on when allowing admixed samples in model
	unsigned short mutRateMode;	// flag which is turned on when constant mutation rates are assumed across loci
	double varRatesAlpha;			// alpha for a Dirichlet distribution of variable rates across loci
	int genRateRef;					// reference genealogy for updates in rate
	
	// finetune parameters
	UpdateStats finetunes;
	unsigned short doMixing;			// flag which is turned on to allow usage of mising procedure (default is 1)
	int findFinetunes;					//if == 1, dynamically search for finetunes
	int findFinetunesSamplesPerStep;	//if using find-finetunes, this is the number of samples to take before adjusting finetune values
	int findFinetunesNumSteps;  		//if using find-finetunes, this is the number of steps before settling in
	
	double* printFactors;			// array of factors in which to output parameters (allocated in readControlFile)
//  char traceFileTitle[500];
}  mcmcSetup;



/*********
 * data setup
 *********/
struct DATA_SETUP {
	int numLoci;					// number of loci in data
	int numSamples;					// number of total samples
	int maxSamples;					// maximum number of samples for allocation purposes
	int numPopPartitions;			// number of partitions to break each pop into for stats
	
	int* numSamplesPerPop;			// number of samples per population
	char**		sampleNames;		// array of sample names - ordered according to population order
	PopulationTree* popTree;		// population tree	
} dataSetup;



/*********
 * admixed samples - for admixture
 *********/
struct ADMIXED_SAMPLES {
	int number;		// number of admixed samples
	int* samples;	// list of admixed samples
	int** popPairs;	// list of population pairs (one per sample)
	int* index;		// index for each admixed sample (-1 for non-admixed)
} admixed_samples;

/***************************************************************************************************************/
/******                               EXTERNAL FUNCTION DECLARATION                                       ******/
/***************************************************************************************************************/



/***********************************************************************************
 *	initGeneralInfo
 * 	- initializes control and I/O settings to default settings
 *	- returns 0
 ***********************************************************************************/
int initGeneralInfo();



/***********************************************************************************
 *	readControlFile
 * 	- reads control file and initializes control and I/O settings
 *	- returns 0, if all OK, and -1 otherwise.
 ***********************************************************************************/
int readControlFile(char* controlFileName);



/***********************************************************************************
 *	readSecondaryControlFile
 * 	- reads secondary control file with only general info and mig bands
 *	- returns 0, if all OK, and -1 otherwise.
 ***********************************************************************************/
int readSecondaryControlFile(char* controlFileName);



/***********************************************************************************
 *	checkSettings
 * 	- tests validity and completeness of settings collected in control file(s)
 *	- returns number of errors found
 ***********************************************************************************/
int checkSettings();



/***********************************************************************************
 *	printPriorSettings
 * 	- prints the prior settings onto standard output
 *	- returns 0
 ***********************************************************************************/
int printPriorSettings();



/***********************************************************************************
 *	finalizeNumParameters
 * 	- determines number of parameters in the model, and finalizes printFactor array
 *	- returns 0
 ***********************************************************************************/
int finalizeNumParameters();



/***************************************************************************************************************/
/******                                        END OF FILE                                                ******/
/***************************************************************************************************************/
#endif
