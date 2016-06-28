/** 
   \file GPhoCS.c 
    The MAIN file containing code for the G-PhoCS program, implementation of the MCMC sampling procedure
 */


#include "utils.h"
#include "MCMCcontrol.h"
#include "AlignmentProcessor.h"
#include "GenericTree.h"
#include "PopulationTree.h"
#include <getopt.h>
#include "LocusDataLikelihood.h"				// NEXTGEN: switch to LocusGenealogy.h !!!
#include <unistd.h>
#include <stdlib.h>

/******************************************************************************************************/
/******                                      CONSTANTS                                           ******/
/******************************************************************************************************/

#define LOG_STEPS_NOT
#define CHECKALL_NOT
#define CHECKCLADE

#define NUM_TYPES						5
#define TARGET_ACCEPTANCE_PERCENT       35
#define TARGET_ACCEPTANCE_RANGE		    5
#define FINETUNE_RESOLUTION             0.0000001
#define MAX_FINETUNE					10
#define ACCEPTANCE_FUDGE 				2

int typeCount[NUM_TYPES];

/***************************************************************************************************************/
/******                                              DATA TYPES                                           ******/
/***************************************************************************************************************/



/***************************************************************************************************************/
/******                                  GLOBAL DATA STRUCTURES                                           ******/
/***************************************************************************************************************/



/*********
 * data setup
 *********/
struct DATA_STATE {
	double	logLikelihood;			// average log-likelihood per genealogy of data given pop tree - {ln[P(X|Z)]+ln[P(Z|M)]}/numLoci
	double	dataLogLikelihood;		// log likelihood (not averaged) of data given all genealogies - ln[P(X|Z,M,T)]
	double	genealogyLogLikelihood;	// log likelihood of all genealogies given model & parameters - ln[P(Z|M,T)]
	double	rateVar;				// the actual variance in locus-specific mutation rate

	LocusData**	lociData;			// array of LocusData data structures (of length numLoci). (allocated in processAlignments)
} dataState;



/*********
 * misc stats
 *********/
struct MISC_STATS{
	int rubberband_mig_conflicts;	// number of rubber band conflicts with migration nodes
	int spr_zero_targets;			// number of times an SPR event encounters zero target edges
	//	int small_interval;				// very small interval for moving coalescent event or migration event
	int not_enough_migs;			// number of times not enough pre-allocated space for migration nodes
	//	double spr_lnld_disc;			// the size of the smallest discrepancy in log-likelihood computation
} misc_stats;

static struct option long_options[] =
{
		{"help",     no_argument,       0, 'h'},
		{"verbose",  no_argument,       0, 'v'},
		{0, 0, 0, 0}
};

/******************************************************************************************************/
/******                                FUNCTION DECLARATIONS                                     ******/
/******************************************************************************************************/


void printUsage(char *programName);
int processAlignments();
int readRateFile(const char* fileName);
int initLociWithoutData();
void printParamVals (double paramVals[], int startParam, int endParam, FILE* out);
int recordTypes ();
int recordParamVals (double paramVals[]);
int performMCMC();
void printGenealogyAndExit(int gen, int errStatus);
int freeAllMemory();


/** sampling functions **/
int UpdateGB_InternalNode(double finetune);
int UpdateGB_MigrationNode(double finetune);
int UpdateGB_MigSPR();
int UpdateTheta(double finetune);
int UpdateMigRates (double finetune);
void UpdateTau(double *finetunes, int *accepted);
void UpdateSampleAge(double *finetunes, int *accepted);
int UpdateLocusRate(double finetune);
int UpdateAdmixCoeffs(double finetune);
int mixing(double finetune);

/** patch for intermediate G-PhoCS version **/
#include "patch.c"

/******************************************************************************************************/
/******                              FUNCTION IMPLEMENTATION                                     ******/
/******************************************************************************************************/



void printUsage(char *programName) {
	printf("Usage: %s <control-file-name> [secondary-control-file-name] [options].\n",programName);
	printf("-v, --verbose     Print more information at the beginning of the program\n");
	printf("-h, --help\n");
	printf("See manual for more help.\n");
}



/***********************************************************************************
 *	main function
 ***********************************************************************************/
int main (int argc, char*argv[]) {

	int res, c, option_index;

	printf("**********************************************************************\n");
	printf("*                G-Phocs ver 1.2.3        Sept 2015                  *\n");
	printf("**********************************************************************\n\n");
	starttime();
	if(argc<=1) {
		printUsage(argv[0]);
		exit(-1);
	}

	debug = 0;

	while (1) {
		// getopt_long stores the option index here. 
		option_index=0;
		c = getopt_long (argc, argv, "hv",
				long_options, &option_index);

		// Detect the end of the options. 
		if (c == -1)
			break;

		switch (c) {
		case 'v':
			verbose = 1;
			break;

		case 'h':
			printUsage(argv[0]);
			exit(-1);
			break;

		case '?':
			// getopt_long already printed an error message.
			break;

		default:
			abort();
		}
	}

	if(argv[optind] == NULL) {
		printUsage(argv[0]);
		exit(-1);
	}

	printf("Reading control settings from file %s...\n",argv[optind]);
	initGeneralInfo();
	res = readControlFile(argv[optind]);
	if(res != 0) {
		exit(-1);
	}

	// secondary control file
	if(argv[optind+1] != NULL) { 
		printf("Reading control settings from secondary file %s...\n",argv[optind+1]);
		res = readSecondaryControlFile(argv[optind+1]);
		if(res != 0) {
			exit(-1);
		}
	}
	if(dataSetup.popTree->numCurPops > NSPECIES) {
		printf("Error: defined too many populations (%d), maximum allowed is %d.\n", dataSetup.popTree->numCurPops, NSPECIES);
		printf("Please set NSPECIES constant at top of patch.c source file to at least %d, recomplie, and re-run.\n", dataSetup.popTree->numCurPops);
		exit(-1);
	}
	if(dataSetup.popTree->numMigBands > MAX_MIG_BANDS) {
		printf("Error: defined too many migration bands (%d), maximum allowed is %d.\n", dataSetup.popTree->numMigBands, MAX_MIG_BANDS);
		printf("Please set MAX_MIG_BANDS constant at top of patch.c source file to at least %d, recomplie, and re-run.\n", dataSetup.popTree->numMigBands);
		exit(-1);
	}
	printf("Done.\n");

	res = checkSettings();
	finalizeNumParameters();

	if(res > 0) {
		fprintf(stderr, "Found %d errors when processing control settings.\n",res);
		exit(-1);
	}

	if(mcmcSetup.randomSeed < 0) {
		mcmcSetup.randomSeed=abs(2*(int)time(NULL)+1);
	}
	setSeed(mcmcSetup.randomSeed);
	if(verbose) {
		printPriorSettings();
		printf("\nRandom seed set to %d\n", mcmcSetup.randomSeed);
	}

	if(mcmcSetup.useData) {
		res = processAlignments();
		if(res < 0) {
			exit(-1);
		}
	} else {
		res = initLociWithoutData();
		if(res < 0) {
			exit(-1);
		}
	}

	if(dataSetup.numSamples > NS) {
		printf("Error: defined too many samples (%d), maximum allowed is %d.\n", dataSetup.numSamples, NS);
		printf("Please set NS constant at top of patch.c source file to at least %d, recomplie, and re-run.\n", dataSetup.numSamples);
		exit(-1);
	}

	GetMem();
	printf("\n");

	performMCMC();
	// MAYBE PERFORM SOME SUMMARIES HERE BEFORE CLOSING ??? !!!
	//	printf("Summarizing statistics, time reset.");
	//	fprintf(fout,"\nSummary of MCMC results:\n");
	//	DescriptiveStatistics(fout, mcmcout, 50, 20, 0);


	freeAlignmentData();
	freeAllMemory();
	exit(0);
}
/** end of main **/



/***********************************************************************************
 *	processAlignments
 * 	- processes alignment data from file.
 *	- deals with het genotypes by summing over all phases
 *	- initializes data structures for gene trees using data
 *	- returns 0
 ***********************************************************************************/
int processAlignments(){
	int res, gen;
	int patt, numPhasedPatterns, maxNumPatterns; // UNUSED numPatterns;
	int* numPhasesArray;
	char** patternArray;
	char** phasedPatternArray;
	int totalNumPatterns, totalPhasedPattern;

	res = readSeqFile(ioSetup.seqFileName,dataSetup.numSamples,dataSetup.sampleNames,dataSetup.numLoci);
	if(res < 0) {
		//fprintf(stderr, "Error: Problem occurred while reading sequence file.\n");
		//printAlignmentError();
		//freeAlignmentData();
		return -1;
	}
	dataSetup.numLoci = AlignmentData.numLoci;

	if(verbose)
		printf("Found %d patterns in %d loci over %d samples.\n", AlignmentData.numPatterns, dataSetup.numLoci, dataSetup.numSamples);

	maxNumPatterns = 4*AlignmentData.numPatterns;

	patternArray = (char**)malloc(AlignmentData.numPatterns*sizeof(char*));
	if(patternArray == NULL) {
		fprintf(stderr, "Error: Out Of Memory when trying to allocate patternArray in processAlignments.\n");
		freeAlignmentData();
		return -1;
	}
	phasedPatternArray = (char**)malloc(maxNumPatterns*sizeof(char*));
	if(phasedPatternArray == NULL) {
		fprintf(stderr, "Error: Out Of Memory when trying to allocate phasedPatternArray in processAlignments.\n");
		freeAlignmentData();
		free(patternArray);
		return -1;
	}

	phasedPatternArray[0] = (char*)malloc(dataSetup.numSamples*maxNumPatterns*sizeof(char));
	if(phasedPatternArray[0] == NULL) {
		fprintf(stderr, "Error: Out Of Memory when trying to allocate phasedPatternArray space in processAlignments.\n");
		freeAlignmentData();
		free(phasedPatternArray);
		free(patternArray);
		return -1;
	}

	numPhasesArray = (int*)malloc(maxNumPatterns*sizeof(int));
	if(numPhasesArray == NULL) {
		fprintf(stderr, "Error: Out Of Memory when trying to allocate numPhasesArray in processAlignments.\n");
		freeAlignmentData();
		free(phasedPatternArray[0]);
		free(phasedPatternArray);
		free(patternArray);
		return -1;
	}

	// initialize gene trees
	if(verbose)
		printf("Initializing %d genealogies with %d leaves...\n",dataSetup.numLoci, dataSetup.numSamples);
	dataState.lociData = (LocusData**)malloc(dataSetup.numLoci*sizeof(LocusData*));
	if(dataState.lociData == NULL) {
		fprintf(stderr, "Error: Out Of Memory when trying to allocate lociData array.\n");
		freeAlignmentData();
		free(patternArray);
		free(phasedPatternArray[0]);
		free(phasedPatternArray);
		free(numPhasesArray);
		return -1;
	}
	totalNumPatterns = totalPhasedPattern = 0;
	for(gen=0; gen<dataSetup.numLoci; gen++) {
		//    		printf("\n gen %d.\n",gen+1);
		dataState.lociData[gen] = createLocusData(dataSetup.numSamples, 1);
		if(dataState.lociData[gen] == NULL) {
			fprintf(stderr, "Error: Out Of Memory when creating genealogy %d.\n",gen+1);
			freeAlignmentData();
			free(patternArray);
			free(phasedPatternArray[0]);
			free(phasedPatternArray);
			free(numPhasesArray);
			return -1;
		}

		// UNUSED    numPatterns = 0;
		for(patt=0; patt<AlignmentData.locusProfiles[gen].numPatterns; patt++) {
			patternArray[patt] = AlignmentData.patternArray[ AlignmentData.locusProfiles[gen].patternIds[patt] ];
		}
		numPhasedPatterns = processHetPatterns(patternArray, AlignmentData.locusProfiles[gen].patternCounts, AlignmentData.locusProfiles[gen].numPatterns,
				/*break symmetries*/ 1, &phasedPatternArray, &numPhasesArray, &maxNumPatterns);

		if(numPhasedPatterns < 0) {
			fprintf(stderr, "Error: Number of phased patterns is negative, unable to process het patterns for genealogy %d.\n",gen+1);
			printAlignmentError();
			freeAlignmentData();
			free(patternArray);
			free(phasedPatternArray[0]);
			free(phasedPatternArray);
			free(numPhasesArray);
			return -1;
		}

		totalNumPatterns += AlignmentData.locusProfiles[gen].numPatterns;
		totalPhasedPattern += numPhasedPatterns;

		res = initializeLocusData(dataState.lociData[gen], phasedPatternArray, numPhasedPatterns,
				numPhasesArray, AlignmentData.locusProfiles[gen].patternCounts);

		if(res < 0) {
			fprintf(stderr, "Error: Unable to initialize locus data, which is necessary to initialize genealogy %d.\n",gen+1);
			freeAlignmentData();
			free(patternArray);
			free(phasedPatternArray[0]);
			free(phasedPatternArray);
			free(numPhasesArray);
			return -1;
		}

	}// end of for(gen)


	free(patternArray);
	free(phasedPatternArray[0]);
	free(phasedPatternArray);
	free(numPhasesArray);
	if(verbose)
		printf("Done. Total of %d patterns (%lf average per locus) transformed to %d phased patterns (%lf average per locus).\n",
				totalNumPatterns, ((double)totalNumPatterns)/dataSetup.numLoci, totalPhasedPattern, ((double)totalPhasedPattern)/dataSetup.numLoci);

	return 0;

}
/** end of processAlignments **/



/***********************************************************************************
 *	initLociWithoutData
 * 	- initializes locus info when no sequence data is available
 *	- done for sampling according to priors
 *	- returns 0 if OK, -1 if error
 ***********************************************************************************/
int initLociWithoutData () {
	int locus;

	if(dataSetup.numLoci <= 0) {
		fprintf(stderr, "Error: when using no sequence data, a positive number of loci should be explicitly specified in control file.\n");
		return (-1);
	}
	printf("Initializing locus data without sequences for %d loci.\n",dataSetup.numLoci);

	dataState.lociData = (LocusData**)malloc(dataSetup.numLoci*sizeof(LocusData*));
	if(dataState.lociData == NULL) {
		fprintf(stderr, "\n Error: Out Of Memory array of locus data pointers.\n");
		return (-1);
	}

	// create empty loci without sequence data
	for(locus=0; locus<dataSetup.numLoci; locus++) {
		dataState.lociData[locus] = createLocusData(dataSetup.numSamples, 0);
		if(dataState.lociData[locus] == NULL) {
			fprintf(stderr, "\n Error: Unable to create locus data formatting for locus %d.\n",locus+1);
			return (-1);
		}
	}
	return 0;
}
/** end of initLociWithoutData **/



/***********************************************************************************
 *	readRateFile
 * 	- reads locus mutation rates from file - and normalizes (so that average is 1)
 * 	- returns 0, if no errors, and -1 otherwise.
 ***********************************************************************************/
int  readRateFile(const char* fileName)  {

	FILE *frate=fopen(fileName,"r");
	int res, locus; // UNUSED numZero;
	double rateSum,tmp;
	double* rates = (double*)malloc(dataSetup.numLoci*sizeof(double));


	if(frate == NULL) {
		fprintf(stderr, "Error: Could not find/read rate file %s.\n",fileName);
		exit(-1);
	}

	if(rates == NULL) {
		fprintf(stderr, "Error: Out Of Memory: Could not allocate rates array in readRateFile() function.\n");
		exit(-1);
	}

	printf("Reading locus rates from file %s... ",fileName);
	// UNUSED  numZero=0;
	rateSum = 0.0;
	for(locus=0; locus<dataSetup.numLoci; locus++) {
		res = fscanf(frate,"%lf",&rates[locus]);
		if(res != 1) {
			fprintf(stderr, "Error: Cannot read rate for locus %d.\n",locus+1);
			return -1;
		}

		rateSum += rates[locus];
		if(rates[locus] <= 0.0) {
			fprintf(stderr, "Error: Locus %d has non-positive (%g) rate.\n", locus+1, rates[locus]);
			return -1;
		}
	}// end of for(locus)

	if(locus<dataSetup.numLoci) {
		fprintf(stderr, "Error: Number of loci read from rate file (%d), was less than the number of loci (%d) specified in the sequence file.\n", locus, dataSetup.numLoci);
		free(rates);
		fclose(frate);
		return -1;
	}

	if(fscanf(frate,"%lf",&tmp) == 1) {
		fprintf(stderr, "Error: Rate file contains more than the %d loci specified in the sequence file.\n",dataSetup.numLoci);
		free(rates);
		fclose(frate);
		return -1;
	}

	fclose(frate);

	rateSum /= dataSetup.numLoci;

	if(verbose)
		printf("Average rate %f. Normalizing for average 1... ", rateSum);

	dataState.rateVar = 0.0;
	for(locus=0; locus<dataSetup.numLoci; locus++) {
		setLocusMutationRate(dataState.lociData[locus],rates[locus]/rateSum);
		dataState.rateVar += (rates[locus]/rateSum-1)*(rates[locus]/rateSum-1);
	}// end of for(locus)
	dataState.rateVar /= dataSetup.numLoci;

	if(verbose)
		printf("Rate variance set at %f.\n", dataState.rateVar);
	free(rates);
	return 0;
}
/** end of readRateFile **/



/***********************************************************************************
 *	freeAllMemory
 * 	- in charge of freeing all memory allocated by sampler
 *	- returns 0
 ***********************************************************************************/
int	freeAllMemory() {
	int gen,i;

	//Closing files
	if(ioSetup.debugFile != NULL)	fclose(ioSetup.debugFile);
	if(ioSetup.traceFile != NULL)	fclose(ioSetup.traceFile);
	if(ioSetup.coalStatsFile != NULL)fclose(ioSetup.coalStatsFile);
	if(ioSetup.nodeStatsFile != NULL) {
		for(i=0; i<3*dataSetup.popTree->numPops; i++) {
			fclose(ioSetup.nodeStatsFile[i]);
		}
		free(ioSetup.nodeStatsFile);
	}
	if(ioSetup.admixFile != NULL)	fclose(ioSetup.admixFile);
	//Freeing print factors array
	free(mcmcSetup.printFactors);
	// NEXTGEN - change this to locusGenealogy !!
	//Freeing population tree
	freePopTree(dataSetup.popTree);
	//Freeing loci data array
	for(gen=0; gen<dataSetup.numLoci; gen++) {
		freeLocusData(dataState.lociData[gen]);
	}
	//Freeing loci data
	free(dataState.lociData);
	//Freeing sample arrays
	free(dataSetup.numSamplesPerPop);
	free(dataSetup.sampleNames[0]);
	free(dataSetup.sampleNames);

	//Freeing genetree migs
	free(genetree_migs);

	//Freeing alignment data
	freeAlignmentData();

	// NEXTGEN - NEED TO REMOVE THIS PART !!!
	//Freeing event chains
	free(event_chains[0].events);
	free(event_chains);
	free(genetree_stats);
	free(rubberband_migs);
	free(nodePops[0]);
	free(nodePops);
	//Done

	return 0;
}
/** end of freeAllMemory **/



/***********************************************************************************
 *	printGenealogyAndExit
 * 	- prints current information on a (faulty) genealogy (including population tree), and frees all memory
 *	- if errStatus != 0, exits with that status. Otherwise don't exit
 ***********************************************************************************/
void printGenealogyAndExit(int gen, int errStatus) {

	if(debug) {
		printPopulationTree(dataSetup.popTree, stderr, 1);
		printLocusGenTree(dataState.lociData[gen], stderr, nodePops[gen], nodeEvents[gen]);
		printEventChains(stderr, gen);
	}
	freeAllMemory();

	if(errStatus != 0)		exit(errStatus);

	return;
}
/** end of printGenealogyAndExit **/



/***********************************************************************************
 *	recordTypes
 * 	- records and prints the number of genealogies of each type as defined in code
 *	- returns 0;
 ***********************************************************************************/
int recordTypes ()	{
	int type, gen; // UNUSED targetPop;
	int numMigs;


	// UNUSED  targetPop = dataSetup.popTree->migBands[0].targetPop;

	for(type=0; type<NUM_TYPES; type++) {
		typeCount[type] = 0;
	}

	for(gen=0; gen<dataSetup.numLoci; gen++) {
		numMigs = genetree_stats[gen].num_migs[0];
		// lineages at end of pop
		/*		numAncestralLineages = event_chains[gen].events[targetPop].num_lineages;
            if(numAncestralLineages == 2 && numMigs == 0) {
            typeCount[0]++;
            } else if(numAncestralLineages == 1 && numMigs == 0) {
            typeCount[1]++;
            } else if(numAncestralLineages == 1 && numMigs == 1) {
            typeCount[2]++;
            } else if(numAncestralLineages == 0 && numMigs == 1) {
            typeCount[3]++;
            } else if(numAncestralLineages == 0 && numMigs == 2) {
            typeCount[4]++;
            } else {
            printf("\n Illegal gen type for gen %d. num-migs = %d, num-ancestral-lineages = %d.\n",gen,numMigs, numAncestralLineages);
            printGenealogyAndExit(gen,-1);
            }
		 */
		if(numMigs < NUM_TYPES-1) {
			typeCount[numMigs]++;
		} else {
			typeCount[NUM_TYPES-1]++;
		}
	}// end of for(gen)

	//    for(type=0; type<NUM_TYPES; type++) { //We removed the output of type_* to the trace file
	//  fprintf(ioSetup.traceFile," %04.3f",((double)typeCount[type])/dataSetup.numLoci);
	// }

	return 0;
}
/** end of recordTypes **/


/***********************************************************************************
 *	printParamVals
 * 	- prints parameter values to out file (without newline)
 *	- uses print factors to factor parameters
 ***********************************************************************************/
void printParamVals (double* paramVals, int startParam, int endParam, FILE* out)	{
	int i;
	for(i=startParam; i<endParam; i++) {
		fprintf(out, "%8.5f\t",paramVals[i]*mcmcSetup.printFactors[i]);
	}
}
/** end of printParamVals **/



/***********************************************************************************
 *	recordAdmixtureCounts
 * 	- records admixture coefficients
 *	- returns 0
 ***********************************************************************************/
int recordAdmixtureCounts() {
	int sample, locus;
	for(sample=0; sample<admixed_samples.number; sample++) {
		admixture_status.admixtureCounts[sample] = 0.0;
		for(locus=0; locus<dataSetup.numLoci; locus++) {
			if(nodePops[locus][ admixed_samples.samples[sample] ] == admixed_samples.popPairs[sample][1]) {
				admixture_status.admixtureCounts[sample] += 1.0;
			}
		}
		for(locus=0; locus<admixture_status.numSampledLoci; locus++) {
			if(nodePops[ admixture_status.sampledLoci[locus] ][ admixed_samples.samples[sample] ] == admixed_samples.popPairs[sample][1]) {
				admixture_status.sampleLocusAdmixRate[sample][locus] += 1.0;
			}
		}
		//		admixture_status.admixtureCoefficients[sample] /= dataSetup.numLoci;
	}

	return 0;

}
/** end of recordAdmixtureCounts **/


/***********************************************************************************
 *	recordParamVals
 * 	- previously collectx()
 *	- records values of the model parameters in a given double array
 *	- returns 0
 ***********************************************************************************/
int recordParamVals (double* paramVals)		{

	int pop, migBand, ind, sample;

	ind = 0;

	for(pop=0; pop<dataSetup.popTree->numPops; pop++) {
		paramVals[ind++] = dataSetup.popTree->pops[pop]->theta;
	}
	for(pop=dataSetup.popTree->numCurPops; pop<dataSetup.popTree->numPops; pop++) {
		paramVals[ind++] = dataSetup.popTree->pops[pop]->age;
	}
	for(migBand=0; migBand<dataSetup.popTree->numMigBands; migBand++) {
		paramVals[ind++] = dataSetup.popTree->migBands[migBand].migRate;
	}

	// record ages of ancient populations
	for(pop=0; pop<dataSetup.popTree->numCurPops;pop++) {
		if(dataSetup.popTree->pops[pop]->updateSampleAge || dataSetup.popTree->pops[pop]->sampleAge > 0.0) {
			paramVals[ind++] = dataSetup.popTree->pops[pop]->sampleAge;
		}
	}

	//  recordAdmixtureCounts();
	for(sample=0; sample<admixed_samples.number; sample++) {
		paramVals[ind++] = admixture_status.admixtureCoefficients[sample];
	}

	// if variable mutation rates, record the relative rate of first gen
	// if variable mutation rates, record the variance in locus-specific mutation rate
	if(mcmcSetup.mutRateMode == 1) {
		//		paramVals[ind++] = getLocusMutationRate (dataState.lociData[ mcmcSetup.locusRef ]);
		paramVals[ind++] = sqrt(dataState.rateVar);
	}
	return 0;
}
/** end of recordParamVals **/



/***********************************************************************************
 *	getLogPrior
 * 	- compute log of prior distribution [ WITHOUT LOCUS-SPECIFIC MUT RATES ]
 ***********************************************************************************/
/*** computes gamma density function  ***/
double getLogGammaDist(double alpha, double beta, double val) {
	double logP = 0;
	if(alpha!=1) logP -= lgamma(alpha);
	return logP + alpha*log(beta) + (alpha-1)*log(val) - beta*val; 
}
/*** AUXILIARY FUNCTION  ***/

double getLogPrior ()		{

	int pop, migBand;

	double logPrior = 0.0;

	for(pop=0; pop<dataSetup.popTree->numPops; pop++) {
		logPrior +=  getLogGammaDist(
				dataSetup.popTree->pops[pop]->thetaPrior.alpha,
				dataSetup.popTree->pops[pop]->thetaPrior.beta,
				dataSetup.popTree->pops[pop]->theta  );
	}
	for(pop=dataSetup.popTree->numCurPops; pop<dataSetup.popTree->numPops; pop++) {
		logPrior +=  getLogGammaDist(
				dataSetup.popTree->pops[pop]->agePrior.alpha,
				dataSetup.popTree->pops[pop]->agePrior.beta,
				dataSetup.popTree->pops[pop]->age  );
	}
	for(migBand=0; migBand<dataSetup.popTree->numMigBands; migBand++) {
		logPrior +=  getLogGammaDist(
				dataSetup.popTree->migBands[migBand].migRatePrior.alpha,
				dataSetup.popTree->migBands[migBand].migRatePrior.beta,
				dataSetup.popTree->migBands[migBand].migRate  );
	}

	return logPrior;
}
/** end of getLogPrior **/


/***********************************************************************************
 *	printCoalStats
 * 	- prints coalescent statistics to file
 *  - this includes stats needed to compute null likelihood as well as stats on coal times
 *	- returns 0
 ***********************************************************************************/
int printCoalStats (int iteration)		{

	double logPrior = 0;

	// print header
	if(iteration<0) {
		fprintf(ioSetup.coalStatsFile, "iter\tcoalStatFlat\tnumCoalFlat\tmigStatFlat\tnumMigFlat\tlogPrior\tlogLikelihood\tdataLogLikelihood\tgenealogyLogLikelihood");
		fprintf(ioSetup.coalStatsFile, "\n");

		return 0;
	}// end of if(iteration<0)

	logPrior = getLogPrior();


	fprintf(ioSetup.coalStatsFile, "%7d", iteration);
	fprintf(ioSetup.coalStatsFile, "\t%9f\t%9d\t%9f\t%9d\t%9f\t%9f\t%9f\t%9f",
			genetree_stats_flat.coal_stats_flat,
			genetree_stats_flat.num_coals_total,
			genetree_stats_flat.mig_stats_flat,
			genetree_stats_flat.num_migs_total,
			logPrior,
			dataState.logLikelihood*dataSetup.numLoci,
			dataState.dataLogLikelihood,
			dataState.genealogyLogLikelihood
	);

	fprintf(ioSetup.coalStatsFile, "\n");
	fflush(ioSetup.coalStatsFile);

	return 0;
}
/** end of printCoalStats **/


/***********************************************************************************
 *	initializeAdmixtureStructures
 * 	- initializes data structure for tracing admixture
 * - returns 0
 ***********************************************************************************/
int initializeAdmixtureStructures()	{
	int sample, locus;

	admixture_status.numSampledLoci = dataSetup.numLoci;


	admixed_samples.index = (int*)malloc(admixed_samples.number * sizeof(int));

	for(sample=0; sample<dataSetup.numSamples; sample++) {
		admixed_samples.index[sample] = -1;
	}
	for(sample=0; sample<admixed_samples.number; sample++) {
		admixed_samples.index[ admixed_samples.samples[sample] ] = sample;
		//		printf(" %d\n", admixed_samples.samples[sample]+1);
	}
	//	printf("\n");

	admixture_status.sampledLoci = (int*)malloc(admixture_status.numSampledLoci * sizeof(int));
	admixture_status.admixtureCoefficients = (double*)malloc(admixed_samples.number * sizeof(double));
	admixture_status.sampleLocusAdmixRate = (double**)malloc(admixed_samples.number * sizeof(double*));
	admixture_status.sampleLocusAdmixRate[0] = (double*)malloc(admixture_status.numSampledLoci * admixed_samples.number * sizeof(double));
	admixture_status.admixtureCounts = (int*)malloc(admixed_samples.number * sizeof(int));

	mcmcSetup.numParameters += admixed_samples.number;
	//	printf("Reallocating mcmcSetup.printFactors from %d parameters to %d.\n", mcmcSetup.numParameters-admixed_samples.number, mcmcSetup.numParameters);
	mcmcSetup.printFactors = (double*) realloc(mcmcSetup.printFactors, mcmcSetup.numParameters*sizeof(double));
	for(sample=mcmcSetup.numParameters-admixed_samples.number; sample<mcmcSetup.numParameters; sample++) {
		mcmcSetup.printFactors[sample] = 1.0;
	}
	for(sample=0; sample<admixed_samples.number; sample++) {
		admixture_status.admixtureCoefficients[sample] = 0.5;
		admixture_status.sampleLocusAdmixRate[sample] = admixture_status.sampleLocusAdmixRate[0] + sample * admixture_status.numSampledLoci;
		for(locus=0; locus<admixture_status.numSampledLoci; locus++) {
			admixture_status.sampleLocusAdmixRate[sample][locus] = 0.0;
		}
	}


	for(locus=0; locus<admixture_status.numSampledLoci; locus++) {
		admixture_status.sampledLoci[locus] = locus;
		//		admixture_status.sampledLoci[locus] = (int)(rndu()*dataSetup.numLoci);
	}
	return 0;

}
/** end of initializeAdmixtureStructures **/



/***********************************************************************************
 *	initializeMCMC
 * 	- initializes the state of the MCMC.
 * - returns total number of coalescent nodes in initialized trees
 ***********************************************************************************/
int initializeMCMC()	{

	int gen;
	int totalCoals;
	GenericBinaryTree* tree;	// tree for generating random genealogies
	double totalMutationRate = 0.0;
	double mutationRate = 0.0;

	// initialize population parameters (according to prior)
	samplePopParameters(dataSetup.popTree);

	// initialize gen-specific mutation rates
	// if variable, sample uniformly in [0.8 , 1.2], and normalize
	// so that average rate is 1.0
	if(mcmcSetup.mutRateMode == 0) {
		for(gen=0; gen<dataSetup.numLoci; gen++) {
			setLocusMutationRate(dataState.lociData[gen],1);
		}
		dataState.rateVar = 0.0;
		if(verbose)
			printf("Constant mutation rate for all loci.\n");
	} else if(mcmcSetup.mutRateMode == 2) {
		if(0 != readRateFile(ioSetup.rateFileName)) {
			fprintf(stderr, "Error: Unable to reading rate file '%s'. Aborting !!\n", ioSetup.rateFileName);
			return(-1);
		}
	} else {
		totalMutationRate = 0.0;
		for(gen=0; gen<dataSetup.numLoci; gen++) {
			mutationRate = 0.8 + 0.4*rndu();
			totalMutationRate += mutationRate;
			setLocusMutationRate(dataState.lociData[gen], mutationRate);
		}
		// normalize
		totalMutationRate /= dataSetup.numLoci;
		dataState.rateVar = 0.0;
		for(gen=0; gen<dataSetup.numLoci; gen++) {
			mutationRate = getLocusMutationRate(dataState.lociData[gen]) / totalMutationRate;
			setLocusMutationRate(dataState.lociData[gen], mutationRate);
			dataState.rateVar += (mutationRate-1)*(mutationRate-1);
		}
		dataState.rateVar /= dataSetup.numLoci;
		mcmcSetup.genRateRef = 0;
	}



	// initialize genealogies by random sampling (according to pop parameters
	if(verbose)
		printf("Initialize %d genealogies... ",dataSetup.numLoci);
	totalCoals=0;
	dataState.logLikelihood = 0.0;
	dataState.dataLogLikelihood = 0.0;
	tree = createGenericTree(dataSetup.numSamples);
	if(tree == NULL) {
		fprintf(stderr, "\nError:Out Of Memory generic tree at beginning of MCMC.\n");
		exit(-1);
	}
	//	printSitePatterns();
	//	printLocusProfiles();
	for(gen=0; gen<dataSetup.numLoci; gen++) {
		//		printf(" %d",gen+1);
		//		fflush(stdout);
		totalCoals += dataSetup.numSamples-1;
		GetRandomGtree(tree, gen);
		copyGenericTreeToLocus(dataState.lociData[gen],tree);
		constructEventChain(gen);
		computeGenetreeStats(gen);
		genLogLikelihood[gen] = gtreeLnLikelihood(gen);
		dataState.logLikelihood += genLogLikelihood[gen];
		//		printLocusDataPatterns(dataState.lociData[gen],stdout);
		//		printf("likelihood computation for gen %d:\n",gen+1);
		//		printf("-----------------------------------\n",gen+1);
		dataState.dataLogLikelihood += computeLocusDataLikelihood (dataState.lociData[gen], /* reset values */ 0);
		resetSaved(dataState.lociData[gen]);
	}
	if(verbose)
		printf("Done.\n");
	freeGenericTree(tree);
	computeTotalStats();
	dataState.logLikelihood = (dataState.logLikelihood + dataState.dataLogLikelihood) /dataSetup.numLoci;

	return totalCoals;

}/** end of initializeMCMC **/

void printTraceFileHeader() {

	int pop, migBand, sample;

	fprintf(ioSetup.traceFile, "Sample");
	for (pop = 0; pop < dataSetup.popTree->numPops; pop++) {
		fprintf(ioSetup.traceFile, "\ttheta_%s",
				dataSetup.popTree->pops[pop]->name);
	}
	for (pop = dataSetup.popTree->numCurPops; pop < dataSetup.popTree->numPops;
			pop++) {
		fprintf(ioSetup.traceFile, "\ttau_%s",
				dataSetup.popTree->pops[pop]->name);
	}
	for (migBand = 0; migBand < dataSetup.popTree->numMigBands; migBand++) {
		fprintf(ioSetup.traceFile, "\tm_%s->%s",
				dataSetup.popTree->pops[dataSetup.popTree->migBands[migBand].sourcePop]->name,
				dataSetup.popTree->pops[dataSetup.popTree->migBands[migBand].targetPop]->name);
	}
	for (pop = 0; pop < dataSetup.popTree->numCurPops; pop++) {
		if (dataSetup.popTree->pops[pop]->updateSampleAge
				|| dataSetup.popTree->pops[pop]->sampleAge > 0.0) {
			fprintf(ioSetup.traceFile, "\ttau_%s",
					dataSetup.popTree->pops[pop]->name);
		}
	}
	for (sample = 0; sample < admixed_samples.number; sample++) {
		fprintf(ioSetup.traceFile, "\tA%d[%s]",
				admixed_samples.samples[sample],
				dataSetup.popTree->pops[admixed_samples.popPairs[sample][1]]->name);
	}
	if (mcmcSetup.mutRateMode == 1)
		fprintf(ioSetup.traceFile, "\tVariance-Mut");

	fprintf(ioSetup.traceFile, "\tData-ld-ln\tFull-ld-ln\tGene-ld-ln\n");
}


/* shortcut to checkall(), print error message and exit
 */
void checkAllAndExit(char errorMessage[], int iteration) {
	if (!checkAll()) {
		fprintf(stderr, errorMessage, iteration);
		exit(-1);
	}
}

/***********************************************************************************
 *	performMCMC
 * 	- main procedure in program.
 *	- performs the MCMC random walk
 ***********************************************************************************/
int performMCMC()	{


	int totalCoals;
	int  gen, iteration, pop, sample;
	double *paramVals, *paramMeans, *doubleArray;
	UpdateStats acceptanceCounts, acceptancePercents, finetuneMaxes, finetuneMins;

	int acceptCount;
	int *acceptCountArray = malloc(sizeof(int) *dataSetup.popTree->numPops);
	int i, j, logCount, totalNumMigNodes;
	int numSamplesPerLog, logsPerLine;

	unsigned short findingFinetunes = 0; // set to 1 while dynamically searching for finetunes
	unsigned short recordCoalStats = (0 != strcmp(ioSetup.nodeStatsFileName,"NONE")); // set to 1 for recording coal stats

	char timeString[STRING_LENGTH];
	char fileName[NAME_LENGTH];


	ioSetup.traceFile = fopen(ioSetup.traceFileName,"w");
	if(ioSetup.traceFile == NULL) {
		fprintf(stderr, "Error: Could not open trace file %s.\n", ioSetup.traceFileName);
		return(-1);
	}

	if(recordCoalStats) {
		ioSetup.coalStatsFile = fopen(ioSetup.nodeStatsFileName,"w");
		if(ioSetup.coalStatsFile == NULL) {
			fprintf(stderr, "Error: Could not open coal stats file %s.\n", fileName);
			return(-1);
		}

		printCoalStats(-1);
	}

#ifdef LOG_STEPS
	ioSetup.debugFile = fopen("G-PhoCS-debug.txt","w");
#endif
	if(admixed_samples.number > 0) {
		initializeAdmixtureStructures();
		strncpy(ioSetup.admixFileName, "admixture-trace.out", NAME_LENGTH);
		//	ioSetup.debugFile = fopen(ioSetup.debugFileName,"w");
	}

	printTraceFileHeader();

	printf("Starting MCMC: %d burnin, %d running, sampled every %d iteration(s).\n", mcmcSetup.burnin, mcmcSetup.numSamples, mcmcSetup.sampleSkip);
	//    printf("Updating genealogies %d times between parameter updates, and starting to sample migration after %d iterations\n", mcmcSetup.genetreeSamples, mcmcSetup.startMig);

	totalCoals = initializeMCMC();
	if(totalCoals <= 0) {
		printf("Error while initializing MCMC.\n");
		return -1;
	}
	// allocate and initialize parameter value arrays
	printf("There are %d parameters in the model.\n", mcmcSetup.numParameters);
	doubleArray = (double*)malloc((2*mcmcSetup.numParameters + 4*dataSetup.popTree->numPops)*sizeof(double));
	if(doubleArray == NULL) {
		fprintf(stderr, "\nError: Out Of Memory while allocating double array in performMCMC.\n");
		exit(-1);
	}
	paramVals = doubleArray;
	paramMeans = paramVals + mcmcSetup.numParameters;
	acceptanceCounts.taus = paramMeans + mcmcSetup.numParameters;
	acceptancePercents.taus = acceptanceCounts.taus + dataSetup.popTree->numPops;
	finetuneMaxes.taus = acceptancePercents.taus + dataSetup.popTree->numPops;
	finetuneMins.taus = finetuneMaxes.taus + dataSetup.popTree->numPops;

	recordParamVals(paramVals);


	if(verbose) {
		printf("Initial parameters: ");
		printParamVals(paramVals, 0, mcmcSetup.numParameters, stdout);
		printf(", data log-likelihood=%g\n\n", dataState.logLikelihood);
	}
	// title for log
	printf("Samples   CoalTimes MigTimes  SPRs      Thetas    MigRates ");
	if(admixed_samples.number > 0) {
		printf("    AdmxCoefs ");
	}
	for(pop=0;pop<dataSetup.popTree->numPops;pop++) {
		if(pop>=dataSetup.popTree->numCurPops || dataSetup.popTree->pops[pop]->updateSampleAge ) {
			printf("TAU_%2d    ", pop);
		}
	}

	printf("RbberBnd  MutRates  Mixing    | DATA-ln-ld |  TIME\n");
	printf("-----------------------------------------------------------------------------------------------------------------------------------------------------\n");

	fflush(stdout);

	// reset statistics and counts before main loop
	logCount = 1;
	for(i=0; i<mcmcSetup.numParameters; i++) {
		paramMeans[i] = paramVals[i];
	}

	totalNumMigNodes = 0;

	misc_stats.rubberband_mig_conflicts = 0;
	misc_stats.spr_zero_targets = 0;
	//	misc_stats.small_interval = 0;
	misc_stats.not_enough_migs = 0;
	//	misc_stats.spr_lnld_disc = 0.0;


	finetuneMaxes.coalTime = MAX_FINETUNE;
	finetuneMaxes.migTime = MAX_FINETUNE;
	finetuneMaxes.theta = MAX_FINETUNE;
	finetuneMaxes.migRate = MAX_FINETUNE;
	finetuneMaxes.locusRate = MAX_FINETUNE;
	finetuneMaxes.admix = MAX_FINETUNE;
	finetuneMaxes.mixing = MAX_FINETUNE;

	finetuneMins.coalTime = 0.0;
	finetuneMins.migTime = 0.0;
	finetuneMins.theta = 0.0;
	finetuneMins.migRate = 0.0;
	finetuneMins.locusRate = 0.0;
	finetuneMins.admix = 0.0;
	finetuneMins.mixing = 0.0;
	acceptanceCounts.coalTime		= 0;
	acceptanceCounts.migTime		= 0;
	acceptanceCounts.SPR			= 0;
	acceptanceCounts.theta		= 0;
	acceptanceCounts.migRate		= 0;
	acceptanceCounts.locusRate	= 0;
	acceptanceCounts.admix		= 0;
	acceptanceCounts.mixing		= 0;
	for(pop=0;pop<dataSetup.popTree->numPops;pop++) {
		acceptanceCounts.taus[pop]		= 0;
	}


	// initialize finetunes for dynamic search
	if(!mcmcSetup.findFinetunes) {
		numSamplesPerLog = ioSetup.samplesPerLog;
		logsPerLine = ioSetup.logsPerLine;
	} else {
		findingFinetunes = 1;
		numSamplesPerLog = mcmcSetup.findFinetunesSamplesPerStep;
		logsPerLine = 1;
		printf("   ---  Dynamically finding finetune settings for the first %d samples, updating finetunes every %d samples  ---- \n",
				numSamplesPerLog*mcmcSetup.findFinetunesNumSteps, numSamplesPerLog);
		printf("-----------------------------------------------------------------------------------------------------------------------------------------------------\n");

		if(mcmcSetup.finetunes.coalTime < 0) mcmcSetup.finetunes.coalTime = 1.0;
		if(mcmcSetup.finetunes.migTime < 0) mcmcSetup.finetunes.migTime = 1.0;
		if(mcmcSetup.finetunes.theta < 0) mcmcSetup.finetunes.theta = 1.0;
		if(mcmcSetup.finetunes.migRate < 0) mcmcSetup.finetunes.migRate = 1.0;
		if(mcmcSetup.finetunes.locusRate < 0) mcmcSetup.finetunes.locusRate = 1.0;
		if(mcmcSetup.finetunes.admix < 0) mcmcSetup.finetunes.admix = 1.0;
		if(mcmcSetup.finetunes.mixing < 0) mcmcSetup.finetunes.mixing = 1.0;
		for(pop=0; pop<dataSetup.popTree->numPops; pop++) {
			if(mcmcSetup.finetunes.taus[pop] < 0) mcmcSetup.finetunes.taus[pop] = 1.0;
			finetuneMaxes.taus[pop] = MAX_FINETUNE;
			finetuneMins.taus[pop] = 0.0;
		}
	}





#ifdef CHECKALL
	checkAllAndExit(stderr, "\nError:  --   Aborting before starting MCMC.\n\n");
#endif


	for(iteration = -mcmcSetup.burnin;
			iteration < mcmcSetup.numSamples;
			iteration++) 	{

		for(j=0; j<mcmcSetup.genetreeSamples; j++) {


			// update COALESCENCE NODE ages
			acceptCount = UpdateGB_InternalNode(mcmcSetup.finetunes.coalTime);
			acceptanceCounts.coalTime += acceptCount;

#ifdef CHECKALL
			checkAllAndExit(stderr, "\nERROR:  --  Aborting after UpdateGB_InternalNode in MCMC iteration %d.\n\n",iteration);
#endif

			// update MIGRATION NODE ages
			acceptCount = UpdateGB_MigrationNode(mcmcSetup.finetunes.migTime);
			acceptanceCounts.migTime += acceptCount;
			// count number of events for acceptance ratio
			for(i=0; i<dataSetup.popTree->numMigBands; i++) {
				totalNumMigNodes += genetree_stats_total.num_migs[i];
			}

#ifdef CHECKALL
			checkAllAndExit(stderr, "\nError:  --  Aborting after UpdateGB_MigrationNode in MCMC iteration %d.\n\n",iteration);
#endif

			// update GENEALOGY TOPOLOGY (including migration events)
			acceptCount = UpdateGB_MigSPR();
			acceptanceCounts.SPR += acceptCount;

#ifdef CHECKALL
			checkAllAndExit(stderr, "\nError:  --  Aborting after UpdateGB_MigSPR in MCMC iteration %d.\n\n",iteration);
#endif

			// update individual LOCUS MUTATION rates
			if(mcmcSetup.mutRateMode == 1) {
				acceptCount = UpdateLocusRate(mcmcSetup.finetunes.locusRate);
				acceptanceCounts.locusRate += acceptCount;

#ifdef CHECKALL
				checkAllAndExit(stderr, "\nError:  --  Aborting after UpdateLocusRate in MCMC iteration %d.\n\n",iteration);
#endif
			}
		}// end of for(j)


		// update THETAs
		acceptCount = UpdateTheta(mcmcSetup.finetunes.theta);
		acceptanceCounts.theta += acceptCount;

#ifdef CHECKALL
		checkAllAndExit(stderr, "\nError:  --  Aborting after UpdateTheta in MCMC iteration %d.\n\n",iteration);
#endif

		// update MIGRATION RATEs
		if(iteration > mcmcSetup.startMig) {
			acceptCount = UpdateMigRates(mcmcSetup.finetunes.migRate);
			acceptanceCounts.migRate += acceptCount;

#ifdef CHECKALL
			checkAllAndExit(stderr, "\nError:  --  Aborting after UpdateMigRates in MCMC iteration %d.\n\n",iteration);
#endif
		}

		// update TAUs

		UpdateTau(mcmcSetup.finetunes.taus, acceptCountArray);
		for(pop=0;pop<dataSetup.popTree->numPops;pop++) {
			acceptanceCounts.taus[pop] += acceptCountArray[pop];
		}

#ifdef CHECKALL
		checkAllAndExit(stderr, "\nError:  --  Aborting after UpdateTau in MCMC iteration %d.\n\n",iteration);
#endif

		UpdateSampleAge(mcmcSetup.finetunes.taus, acceptCountArray);
		for(pop=0;pop<dataSetup.popTree->numPops;pop++) {
			acceptanceCounts.taus[pop] += acceptCountArray[pop];
		}

#ifdef CHECKALL
		checkAllAndExit(stderr, "\nError:  --  Aborting after UpdateSampleAge in MCMC iteration %d.\n\n",iteration);
#endif


		// update admixture coefficients
		acceptCount = UpdateAdmixCoeffs(mcmcSetup.finetunes.admix);
		acceptanceCounts.admix += acceptCount;

#ifdef CHECKALL
		checkAllAndExit("\nError:  --  Aborting after AdmixtureCoefficients in MCMC iteration %d.\n\n",iteration);
#endif


		// mix all population parameters by multiplying with a factor
		// NO MIXING
		if(mcmcSetup.doMixing) {
			acceptCount = mixing(mcmcSetup.finetunes.mixing);
			acceptanceCounts.mixing += acceptCount;
		}



#ifdef CHECKALL
		checkAllAndExit("\n  --  Aborting after mixing() in iteration %d.\n\n",iteration);
#endif

		// synchronize events due to possible inconsistencies caused by rescaling of ages (mixing and rubber band).
		for(gen=0; gen<dataSetup.numLoci; gen++) {
			if (!synchronizeEvents(gen)) {
				printf("\n  --  Aborting due to problems found when synchronizing data structures for locus #%d after MCMC iteration %d.\n\n",iteration, gen+1);
				printGenealogyAndExit(gen,-1);
			}
		}

#ifdef CHECKALL
		checkAllAndExit("\n  --  Aborting after MCMC iteration %d.\n\n", iteration);
#endif


		// record parameters, means, and print to trace, if appropriate
		recordParamVals(paramVals);
		for(i=0; i<mcmcSetup.numParameters; i++) {
			paramMeans[i] = paramMeans[i]*((double)logCount/(logCount+1)) + paramVals[i]/(logCount+1);
		}

		// start sampling migrations
		if(iteration == mcmcSetup.startMig) {
			sampleMigRates(dataSetup.popTree);
			// adjust likelihoods to newly sampled migration rates
			for(gen=0; gen<dataSetup.numLoci; gen++) {
				dataState.logLikelihood -= genLogLikelihood[gen]/dataSetup.numLoci;
				genLogLikelihood[gen] = gtreeLnLikelihood(gen);
				dataState.logLikelihood += genLogLikelihood[gen]/dataSetup.numLoci;
			}
		}

		dataState.genealogyLogLikelihood = dataState.logLikelihood*dataSetup.numLoci - dataState.dataLogLikelihood;

		if(iteration >= 0 && iteration%(mcmcSetup.sampleSkip+1) == 0) {

			fprintf(ioSetup.traceFile, "%d\t", iteration);
			printParamVals(paramVals,0,mcmcSetup.numParameters, ioSetup.traceFile);
			fprintf(ioSetup.traceFile,"%.9f\t%.9f\t%.9f\n", dataState.logLikelihood, dataState.dataLogLikelihood, dataState.genealogyLogLikelihood);
			fflush(ioSetup.traceFile);



			if(recordCoalStats) {
				computeFlatStats();
//				computeNodeStats();
				computeCladeStats();
#ifdef CHECKCLADE
				test_validateRootCladeVsFlatStats();
#endif
				printCoalStats(iteration);
			}

			if(admixed_samples.number > 0 && iteration % 1000 == 0) {
				ioSetup.admixFile = fopen(ioSetup.admixFileName,"w");
				if(ioSetup.admixFile == NULL) {
					fprintf(stderr, "Error: Could not open admixture trace file %s.\n", ioSetup.admixFileName);
					return(-1);
				}
				fprintf(ioSetup.admixFile, "%d", iteration);
				for(sample=0; sample<admixed_samples.number; sample++) {
					for(gen=0; gen<admixture_status.numSampledLoci; gen++) {
						//				fprintf(ioSetup.admixFile,"\t%d", nodePops[ admixture_status.sampledLoci[gen] ][ admixed_samples.samples[sample] ]);
						fprintf(ioSetup.admixFile,"\t%lf", admixture_status.sampleLocusAdmixRate[sample][ admixture_status.sampledLoci[gen] ] / (iteration+1));
					}
				}
				fprintf(ioSetup.admixFile,"\n");
				fclose(ioSetup.admixFile);
				ioSetup.admixFile = NULL;
			}
		}

		logCount++;

		// print log
		if((iteration+1)%numSamplesPerLog == 0) {
			// print the 8 acceptance ratios
			if (!checkAll()) {
				fprintf(stderr, "\nError:  --  Aborting when logging after MCMC iteration %d, due to data structure inconsistency.\n\n",iteration);
				exit(-1);
			}

			acceptancePercents.coalTime = acceptanceCounts.coalTime*100.0/(((double)logCount)*totalCoals*mcmcSetup.genetreeSamples);
			acceptancePercents.migTime = acceptanceCounts.migTime*100.0/(totalNumMigNodes+0.000001);
			acceptancePercents.SPR = acceptanceCounts.SPR*100.0/(((double)logCount)*2*totalCoals*mcmcSetup.genetreeSamples);
			acceptancePercents.theta = acceptanceCounts.theta*100.0/(((double)logCount)*dataSetup.popTree->numPops);
			acceptancePercents.migRate = acceptanceCounts.migRate*100.0/(((double)logCount)*dataSetup.popTree->numMigBands+0.000001);
			for(pop=0;pop<dataSetup.popTree->numPops;pop++) {
				acceptancePercents.taus[pop] = acceptanceCounts.taus[pop]*100.0/(double)logCount;
			}
			acceptancePercents.locusRate = acceptanceCounts.locusRate*100.0/(((double)logCount)*(dataSetup.numLoci-1)*mcmcSetup.genetreeSamples);
			acceptancePercents.mixing = acceptanceCounts.mixing*100.0/logCount;
			if(admixed_samples.number > 0) {
				acceptancePercents.admix = acceptanceCounts.admix*100.0/(((double)logCount)*admixed_samples.number);
			} else {
				acceptancePercents.admix = 0.0;
			}


			printf("\r%7d   %5.1f%%    %5.1f%%    %5.1f%%    %5.1f%%    %5.1f%%    ",
					(iteration+1),
					acceptancePercents.coalTime,
					acceptancePercents.migTime,
					acceptancePercents.SPR,
					acceptancePercents.theta,
					acceptancePercents.migRate);
			if(admixed_samples.number > 0) {
				printf("%5.1f%%    ", acceptancePercents.admix);
			}

			for(pop=0;pop<dataSetup.popTree->numPops;pop++) {
				if(pop>=dataSetup.popTree->numCurPops || dataSetup.popTree->pops[pop]->updateSampleAge) {
					printf("%5.1f%%    ", acceptancePercents.taus[pop]);
				}
			}

			printf("%6.1f%%    %5.1f%%    %5.1f%%    ",
					misc_stats.rubberband_mig_conflicts*100.0/(logCount*(dataSetup.popTree->numPops - dataSetup.popTree->numCurPops)),
					acceptancePercents.locusRate,
					acceptancePercents.mixing);

			// print parameter means
			//printParamVals(paramMeans,dataSetup.popTree->numPops,mcmcSetup.numParameters,stdout);

			// print data log likelihood
			printf("|%12.6f|",dataState.logLikelihood);


			printf(" %s", printtime(timeString));
			if((iteration+1)%(numSamplesPerLog*logsPerLine) == 0) {
				printf("\n");
			}
			fflush(stdout);


			// dynamically adjust finetunes, if applicable
			if(findingFinetunes) {
				if(acceptancePercents.coalTime > TARGET_ACCEPTANCE_PERCENT+TARGET_ACCEPTANCE_RANGE) {
					finetuneMins.coalTime = mcmcSetup.finetunes.coalTime; //raise the minimal finetune value
					if(finetuneMaxes.coalTime - finetuneMins.coalTime < FINETUNE_RESOLUTION) { // recompute maximal finetune
						//				finetuneMaxes.coalTime = MAX_FINETUNE;
						if(finetuneMaxes.coalTime >= MAX_FINETUNE) {
							finetuneMaxes.coalTime = finetuneMins.coalTime = MAX_FINETUNE;
						} else {
							finetuneMaxes.coalTime *= 2.0;
						}
					}
				} else if(acceptancePercents.coalTime < TARGET_ACCEPTANCE_PERCENT-TARGET_ACCEPTANCE_RANGE) {
					finetuneMaxes.coalTime = mcmcSetup.finetunes.coalTime; //lower the maximal finetune value
					if(finetuneMaxes.coalTime - finetuneMins.coalTime < FINETUNE_RESOLUTION) { // recompute minimal finetune
						//				finetuneMins.coalTime = 0.0;
						finetuneMins.coalTime /= 2.0;
					}
				}
				mcmcSetup.finetunes.coalTime = 0.5*(finetuneMaxes.coalTime + finetuneMins.coalTime); //finetune set to midpoint of interval

				if(acceptancePercents.migTime > TARGET_ACCEPTANCE_PERCENT+TARGET_ACCEPTANCE_RANGE) {
					finetuneMins.migTime = mcmcSetup.finetunes.migTime; //raise the minimal finetune value
					if(finetuneMaxes.migTime - finetuneMins.migTime < FINETUNE_RESOLUTION) { // recompute maximal finetune
						//				finetuneMaxes.migTime = MAX_FINETUNE;
						if(finetuneMaxes.migTime >= MAX_FINETUNE) {
							finetuneMaxes.migTime = finetuneMins.migTime = MAX_FINETUNE;
						} else {
							finetuneMaxes.migTime *= 2.0;
						}
					}
				} else if(acceptancePercents.migTime < TARGET_ACCEPTANCE_PERCENT-TARGET_ACCEPTANCE_RANGE) {
					finetuneMaxes.migTime = mcmcSetup.finetunes.migTime; //lower the maximal finetune value
					if(finetuneMaxes.migTime - finetuneMins.migTime < FINETUNE_RESOLUTION) { // recompute minimal finetune
						//				finetuneMins.migTime = 0.0;
						finetuneMins.migTime /= 2.0;
					}
				}
				mcmcSetup.finetunes.migTime = 0.5*(finetuneMaxes.migTime + finetuneMins.migTime); //finetune set to midpoint of interval

				if(acceptancePercents.theta > TARGET_ACCEPTANCE_PERCENT+TARGET_ACCEPTANCE_RANGE) {
					finetuneMins.theta = mcmcSetup.finetunes.theta; //raise the minimal finetune value
					if(finetuneMaxes.theta - finetuneMins.theta < FINETUNE_RESOLUTION) { // recompute maximal finetune
						//				finetuneMaxes.theta = MAX_FINETUNE;
						if(finetuneMaxes.theta >= MAX_FINETUNE) {
							finetuneMaxes.theta = finetuneMins.theta = MAX_FINETUNE;
						} else {
							finetuneMaxes.theta *= 2.0;
						}
					}
				} else if(acceptancePercents.theta < TARGET_ACCEPTANCE_PERCENT-TARGET_ACCEPTANCE_RANGE) {
					finetuneMaxes.theta = mcmcSetup.finetunes.theta; //lower the maximal finetune value
					if(finetuneMaxes.theta - finetuneMins.theta < FINETUNE_RESOLUTION) { // recompute minimal finetune
						//				finetuneMins.theta = 0.0;
						finetuneMins.theta /= 2.0;
					}
				}
				mcmcSetup.finetunes.theta = 0.5*(finetuneMaxes.theta + finetuneMins.theta); //finetune set to midpoint of interval

				if(acceptancePercents.migRate > TARGET_ACCEPTANCE_PERCENT+TARGET_ACCEPTANCE_RANGE) {
					finetuneMins.migRate = mcmcSetup.finetunes.migRate; //raise the minimal finetune value
					if(finetuneMaxes.migRate - finetuneMins.migRate < FINETUNE_RESOLUTION) { // recompute maximal finetune
						//				finetuneMaxes.migRate = MAX_FINETUNE;
						if(finetuneMaxes.migRate >= MAX_FINETUNE) {
							finetuneMaxes.migRate = finetuneMins.migRate = MAX_FINETUNE;
						} else {
							finetuneMaxes.migRate *= 2.0;
						}
					}
				} else if(acceptancePercents.migRate < TARGET_ACCEPTANCE_PERCENT-TARGET_ACCEPTANCE_RANGE) {
					finetuneMaxes.migRate = mcmcSetup.finetunes.migRate; //lower the maximal finetune value
					if(finetuneMaxes.migRate - finetuneMins.migRate < FINETUNE_RESOLUTION) { // recompute minimal finetune
						//				finetuneMins.migRate = 0.0;
						finetuneMins.migRate /= 2.0;
					}
				}
				mcmcSetup.finetunes.migRate = 0.5*(finetuneMaxes.migRate + finetuneMins.migRate); //finetune set to midpoint of interval

				if(acceptancePercents.admix > TARGET_ACCEPTANCE_PERCENT+TARGET_ACCEPTANCE_RANGE) {
					finetuneMins.admix = mcmcSetup.finetunes.admix; //raise the minimal finetune value
					if(finetuneMaxes.admix - finetuneMins.admix < FINETUNE_RESOLUTION) { // recompute maximal finetune
						//				finetuneMaxes.migRate = MAX_FINETUNE;
						if(finetuneMaxes.admix >= MAX_FINETUNE) {
							finetuneMaxes.admix = finetuneMins.admix = MAX_FINETUNE;
						} else {
							finetuneMaxes.admix *= 2.0;
						}
					}
				} else if(acceptancePercents.admix < TARGET_ACCEPTANCE_PERCENT-TARGET_ACCEPTANCE_RANGE) {
					finetuneMaxes.admix = mcmcSetup.finetunes.admix; //lower the maximal finetune value
					if(finetuneMaxes.admix - finetuneMins.admix < FINETUNE_RESOLUTION) { // recompute minimal finetune
						//				finetuneMins.admix = 0.0;
						finetuneMins.admix /= 2.0;
					}
				}
				mcmcSetup.finetunes.admix = 0.5*(finetuneMaxes.admix + finetuneMins.admix); //finetune set to midpoint of interval

				if(acceptancePercents.locusRate > TARGET_ACCEPTANCE_PERCENT+TARGET_ACCEPTANCE_RANGE) {
					finetuneMins.locusRate = mcmcSetup.finetunes.locusRate; //raise the minimal finetune value
					if(finetuneMaxes.locusRate - finetuneMins.locusRate < FINETUNE_RESOLUTION) { // recompute maximal finetune
						//				finetuneMaxes.locusRate = MAX_FINETUNE;
						if(finetuneMaxes.locusRate >= MAX_FINETUNE) {
							finetuneMaxes.locusRate = finetuneMins.locusRate = MAX_FINETUNE;
						} else {
							finetuneMaxes.locusRate *= 2.0;
						}
					}
				} else if(acceptancePercents.locusRate < TARGET_ACCEPTANCE_PERCENT-TARGET_ACCEPTANCE_RANGE) {
					finetuneMaxes.locusRate = mcmcSetup.finetunes.locusRate; //lower the maximal finetune value
					if(finetuneMaxes.locusRate - finetuneMins.locusRate < FINETUNE_RESOLUTION) { // recompute minimal finetune
						//				finetuneMins.locusRate = 0.0;
						finetuneMins.locusRate /= 2.0;
					}
				}
				mcmcSetup.finetunes.locusRate = 0.5*(finetuneMaxes.locusRate + finetuneMins.locusRate); //finetune set to midpoint of interval

				if(acceptancePercents.mixing > TARGET_ACCEPTANCE_PERCENT+TARGET_ACCEPTANCE_RANGE) {
					finetuneMins.mixing = mcmcSetup.finetunes.mixing; //raise the minimal finetune value
					if(finetuneMaxes.mixing - finetuneMins.mixing < FINETUNE_RESOLUTION) { // recompute maximal finetune
						//				finetuneMaxes.mixing = MAX_FINETUNE;
						if(finetuneMaxes.mixing >= MAX_FINETUNE) {
							finetuneMaxes.mixing = finetuneMins.mixing = MAX_FINETUNE;
						} else {
							finetuneMaxes.mixing *= 2.0;
						}
					}
				} else if(acceptancePercents.mixing < TARGET_ACCEPTANCE_PERCENT-TARGET_ACCEPTANCE_RANGE) {
					finetuneMaxes.mixing = mcmcSetup.finetunes.mixing; //lower the maximal finetune value
					if(finetuneMaxes.mixing - finetuneMins.mixing < FINETUNE_RESOLUTION) { // recompute minimal finetune
						//				finetuneMins.mixing = 0.0;
						finetuneMins.mixing /= 2.0;
					}
				}
				mcmcSetup.finetunes.mixing = 0.5*(finetuneMaxes.mixing + finetuneMins.mixing); //finetune set to midpoint of interval

				for(pop=dataSetup.popTree->numCurPops;pop<dataSetup.popTree->numPops;pop++) {
					if(acceptancePercents.taus[pop] > TARGET_ACCEPTANCE_PERCENT+TARGET_ACCEPTANCE_RANGE) {
						finetuneMins.taus[pop] = mcmcSetup.finetunes.taus[pop]; //raise the minimal finetune value
						if(finetuneMaxes.taus[pop] - finetuneMins.taus[pop] < FINETUNE_RESOLUTION) { // recompute maximal finetune
							//					finetuneMaxes.taus[pop] = MAX_FINETUNE;
							if(finetuneMaxes.taus[pop] >= MAX_FINETUNE) {
								finetuneMaxes.taus[pop] = finetuneMins.taus[pop] = MAX_FINETUNE;
							} else {
								finetuneMaxes.taus[pop] *= 2.0;
							}
						}
					} else if(acceptancePercents.taus[pop] < TARGET_ACCEPTANCE_PERCENT-TARGET_ACCEPTANCE_RANGE) {
						finetuneMaxes.taus[pop] = mcmcSetup.finetunes.taus[pop]; //lower the maximal finetune value
						if(finetuneMaxes.taus[pop] - finetuneMins.taus[pop] < FINETUNE_RESOLUTION) { // recompute minimal finetune
							//					finetuneMins.taus[pop] = 0.0;
							finetuneMins.taus[pop] /= 2.0;
						}
					}
					mcmcSetup.finetunes.taus[pop] = 0.5*(finetuneMaxes.taus[pop] + finetuneMins.taus[pop]); //finetune set to midpoint of interval
				}
				printf("          %-9.7lf %-9.7lf           %-9.7lf %-9.7lf %-9.7lf ",
						mcmcSetup.finetunes.coalTime,
						mcmcSetup.finetunes.migTime,
						mcmcSetup.finetunes.theta,
						mcmcSetup.finetunes.migRate,
						mcmcSetup.finetunes.admix);
				for(pop=dataSetup.popTree->numCurPops; pop<dataSetup.popTree->numPops; pop++) {
					printf("%-9.7lf ", mcmcSetup.finetunes.taus[pop]);
				}
				printf("          %-9.7lf %-9.7lf \n",
						mcmcSetup.finetunes.locusRate,
						mcmcSetup.finetunes.mixing);
			}// end of if(findingFinetunes)

			// reset log
			logCount = 1;
			for(i=0; i<mcmcSetup.numParameters; i++) {
				paramMeans[i] = paramVals[i];
			}
			totalNumMigNodes = 0;
			misc_stats.rubberband_mig_conflicts = 0;
			misc_stats.spr_zero_targets = 0;
			//			misc_stats.small_interval = 0;
			misc_stats.not_enough_migs = 0;
			//			misc_stats.spr_lnld_disc = 0.0;

			acceptanceCounts.coalTime		= 0;
			acceptanceCounts.migTime		= 0;
			acceptanceCounts.SPR			= 0;
			acceptanceCounts.theta		= 0;
			acceptanceCounts.migRate		= 0;
			for(pop=0; pop<dataSetup.popTree->numPops;pop++) {
				acceptanceCounts.taus[pop]			= 0;
			}
			acceptanceCounts.locusRate	= 0;
			acceptanceCounts.mixing		= 0;
			acceptanceCounts.admix		= 0;

			if(findingFinetunes && iteration+1 >= mcmcSetup.findFinetunesSamplesPerStep *mcmcSetup.findFinetunesNumSteps) {

				findingFinetunes = 0;
				numSamplesPerLog = ioSetup.samplesPerLog;
				logsPerLine = ioSetup.logsPerLine;
				printf("\n");
				printf("-------------------------------------  finetunes  ------------------------------------\n");
				printf("          %8lf  %8lf            %8lf  %8lf  ",
						mcmcSetup.finetunes.coalTime,
						mcmcSetup.finetunes.migTime,
						mcmcSetup.finetunes.theta,
						mcmcSetup.finetunes.migRate);
				for(pop=0;pop<dataSetup.popTree->numPops;pop++) {
					printf("%8lf  ", mcmcSetup.finetunes.taus[pop]);
				}
				printf("          %8lf  %8lf  \n",
						mcmcSetup.finetunes.locusRate,
						mcmcSetup.finetunes.mixing);
				printf("--------------------------------------------------------------------------------------\n");
			}

		}// print log

	}// end of main loop - for(iteration)


	free(doubleArray);
	free(acceptCountArray);
	printf("\nMCMC finished. Time used: %s\n", printtime(timeString));

	return 0;
}
/** end of performMCMC **/




/******************************************************************************************************/
/******                                 SAMPLING FUNCTIONS                                       ******/
/******************************************************************************************************/



/***********************************************************************************
 *	UpdateGB_InternalNode
 *	- perturbs times of all coalescent nodes in all gene trees
 * 	- does not change the population of the node
 *	- upper and lower bounds for new time are determined by nodes (migration/coalescent)
 *		directly above or below that node, as well as population boundaries.
 *	- records new data log likelihood in *pointerToLnLd
 ***********************************************************************************/
int UpdateGB_InternalNode (double finetune) {
	int  accepted=0, gen, i, inode, son;
	int  pop;
	double lnacceptance, lnLd, t,tnew;
	double tb[2];
	int mig;
	double genetree_lnLd_delta;

	if(finetune <= 0.0) {
		return 0;
	}

	for(gen=0; gen<dataSetup.numLoci; gen++) {
		for(inode=dataSetup.numSamples; inode<2*dataSetup.numSamples-1; inode++) {

			t = getNodeAge(dataState.lociData[gen],inode);
			pop = nodePops[gen][inode];

			// set lower and upper bounds for new age
			tb[0] = dataSetup.popTree->pops[pop]->age;
			if(pop != dataSetup.popTree->rootPop) {
				tb[1] = dataSetup.popTree->pops[pop]->father->age;
			} else {
				tb[1] = OLDAGE;
			}

			mig = findFirstMig(gen, inode,-1);
			if(mig >= 0) {
				tb[1] = min2(tb[1],genetree_migs[gen].mignodes[mig].age);
			} else if(inode != getLocusRoot(dataState.lociData[gen])) {
				tb[1] = min2(tb[1], getNodeAge(dataState.lociData[gen], getNodeFather(dataState.lociData[gen],inode)));
			}
			for(i=0;i<2;i++) {
				son = getNodeSon(dataState.lociData[gen],inode,i);
				mig = findLastMig(gen, son,-1);
				if(mig >= 0) {
					tb[0] = max2(tb[0],genetree_migs[gen].mignodes[mig].age);
				} else {
					tb[0] = max2(tb[0],getNodeAge(dataState.lociData[gen],son));
				}
			}
			tnew = t + finetune*rnd2normal8();
			tnew = reflect(tnew, tb[0], tb[1]);
			/**/
			if(fabs(tnew-t) < 1e-15)	{
				accepted++;
				continue;
			}
			/**/
#ifdef LOG_STEPS
			fprintf(ioSetup.debugFile, "  gen %d, internal node %d, proposing age shift: %g-->%g, ",gen, inode, t, tnew);
#endif
			// update node's age, and compute delta log-likelihood
			adjustGenNodeAge(dataState.lociData[gen], inode, tnew);
			lnLd = -getLocusDataLikelihood(dataState.lociData[gen]);
			lnLd += computeLocusDataLikelihood(dataState.lociData[gen], /*reuse old conditionals*/ 1);

			genetree_lnLd_delta = considerEventMove(gen, 0, nodeEvents[gen][inode], pop, t, pop, tnew);
			lnacceptance = genetree_lnLd_delta + lnLd;
			//			printf("done.\n");

#ifdef LOG_STEPS
			fprintf(ioSetup.debugFile, "lnacceptance = %g, ",lnacceptance);
#endif

			if(lnacceptance>=0 || rndu()<exp(lnacceptance)) {
#ifdef LOG_STEPS
				fprintf(ioSetup.debugFile, "accepting.\n");
#endif
				accepted++;
				genLogLikelihood[gen] += genetree_lnLd_delta;
				dataState.dataLogLikelihood += lnLd;
				dataState.logLikelihood += (genetree_lnLd_delta + lnLd)/dataSetup.numLoci;
				acceptEventChainChanges(gen, 0);
				resetSaved(dataState.lociData[gen]);
			}
			else {
				// reject changes and revert to saved version
#ifdef LOG_STEPS
				fprintf(ioSetup.debugFile, "rejecting.\n");
#endif
				rejectEventChainChanges(gen, 0);
				revertToSaved(dataState.lociData[gen]);
			}
			/**
         if (!checkAll()) {
         printf("\n  --  Aborting after UpdateGB_InternalNode(%d,%d).\n",gen,inode);
         exit(-1);
         }
			 **/
		}// end of for(inode)
	}// end of for(gen)

	return(accepted);
}
/** end of UpdateGB_InternalNode **/



/***********************************************************************************
 *	UpdateGB_MigrationNode
 *	- perturbs times of all migration nodes in all gene trees
 * 	- essentailly follows the same lines as UpdateGB_InternalNode
 *	- likelihood of data given genetree is not altered by this step
 ***********************************************************************************/
int UpdateGB_MigrationNode(double finetune) {

	int gen, m, mignode, mig_below, mig_above, node_below;
	int  pop_source, pop_target, event_source, event_target;
	int accepted = 0;
	double lnacceptance, t,tnew;
	double t_bounds[2];
	double genetree_lnLd_delta;

	int father;

	if(finetune <= 0.0) {
		return 0;
	}

	for(gen=0; gen<dataSetup.numLoci; gen++) {

		for(m=0; m<genetree_migs[gen].num_migs; m++) {
			mignode = genetree_migs[gen].living_mignodes[m];

			t = genetree_migs[gen].mignodes[mignode].age;
			pop_source = genetree_migs[gen].mignodes[mignode].source_pop;
			pop_target = genetree_migs[gen].mignodes[mignode].target_pop;
			event_source = genetree_migs[gen].mignodes[mignode].source_event;
			event_target = genetree_migs[gen].mignodes[mignode].target_event;
			node_below = genetree_migs[gen].mignodes[mignode].gtree_branch;


			// determine upper and lower bounds for new time
			// start up with start and end times of migration band
			// then bound according to events right below or above the migration event
			t_bounds[0] = dataSetup.popTree->migBands[ genetree_migs[gen].mignodes[mignode].migration_band ].startTime;
			t_bounds[1] = dataSetup.popTree->migBands[ genetree_migs[gen].mignodes[mignode].migration_band ].endTime;

			mig_below = findLastMig(gen, node_below,t);
			mig_above = findFirstMig(gen,node_below,t);
			if(mig_below >= 0) {
				t_bounds[0] = max2(t_bounds[0] , genetree_migs[gen].mignodes[mig_below].age);
			} else {
				t_bounds[0] = max2(t_bounds[0] ,getNodeAge(dataState.lociData[gen],node_below));
			}

			if(mig_above >= 0) {
				t_bounds[1] = min2(t_bounds[1] , genetree_migs[gen].mignodes[mig_above].age);
			} else {
				father = getNodeFather(dataState.lociData[gen],node_below);
				if(father < 0) {
					//					printf("ERROR UpdateGB_MigrationNode: migration event %d in gen %d is on edge above root.\n",mignode,gen);
					// genealogy root can actually be in population below root (under migration scenarios)
					//printGenealogyAndExit(gen,-1);
					//					printf("\n Migration event %d above genealogy root at gen %d.\n",mignode,gen);
					t_bounds[1] = min2(t_bounds[1],OLDAGE);
				} else {
					t_bounds[1] = min2(t_bounds[1] , getNodeAge(dataState.lociData[gen],father));
				}
			}

			/*
        if(t_bounds[1]-t_bounds[0] < 0.00000001) {
        //				printf("\n interval of size %g for mignode %d in gen %d.\n",t_bounds[1]-t_bounds[0], mignode, gen);
        misc_stats.small_interval++;
        continue;
        }
			 */
			// Note: migration node cannot move to another population
			// because it is restricted to specific band

			tnew = t + finetune*rnd2normal8();
			tnew = reflect(tnew, t_bounds[0], t_bounds[1]);
			if(fabs(tnew-t) < 1e-15)	{
				accepted++;
				continue;
			}

#ifdef LOG_STEPS
			fprintf(ioSetup.debugFile, "  gen %d, migration node %d, proposing age shift: %g-->%g, ",gen, mignode, t, tnew);
#endif

			//         printEventChains(gen);

			genetree_lnLd_delta = considerEventMove(gen, 0, event_source, pop_source, t, pop_source, tnew);
			//			printEventChains(gen);
			genetree_lnLd_delta += considerEventMove(gen, 1, event_target, pop_target, t, pop_target, tnew);
			lnacceptance = genetree_lnLd_delta;


#ifdef LOG_STEPS
			fprintf(ioSetup.debugFile, "lnacceptance = %g, ",lnacceptance);
#endif
			if(lnacceptance>=0 || rndu()<exp(lnacceptance)) {
#ifdef LOG_STEPS
				fprintf(ioSetup.debugFile, "accepting.\n");
#endif
				accepted++;
				genLogLikelihood[gen] += genetree_lnLd_delta;
				dataState.logLikelihood += genetree_lnLd_delta / dataSetup.numLoci;
				acceptEventChainChanges(gen, 0);
				acceptEventChainChanges(gen, 1);
				genetree_migs[gen].mignodes[mignode].age = tnew;
			}
			else {
#ifdef LOG_STEPS
				fprintf(ioSetup.debugFile, "rejecting.\n");
#endif
				rejectEventChainChanges(gen, 0);
				rejectEventChainChanges(gen, 1);
			}
		}// end of for(mignode)
	}// end of for(gen)

	return(accepted);
}
/** end of UpdateGB_MigrationNode **/



/***********************************************************************************
 *	UpdateGB_MigSPR
 *	- traverses all edges in all gene trees, detaches them and lets them
 * 	- randomly recoalesce with the remaining tree (with migrations).
 ***********************************************************************************/
int UpdateGB_MigSPR ()	{

	int res, accepted=0;
	int i, mig, mig_band, gen, node, event;
	int pop, father, sibling, target, father_pop_old;
	double lnacceptance, lnLd, t_new, heredity_factor; // UNUSED, t_old;
	// UNUSED  unsigned short didAccept;
	unsigned short altPop, admixSwitch;	// flag for admixed samples
	int admixIndex=-1, oldPop=-1, newPop=-1;
	//	double	genetree_lnLd, genetree_lnLd_new;

	for(gen=0; gen<dataSetup.numLoci; gen++) {
		heredity_factor = 1.0;
		//		printLocusGenTree(dataState.lociData[gen], nodePops[gen], nodeEvents[gen]);;
		//		printEventChains(gen);

		//		genetree_lnLd = gtreeLnLikelihood(gen);

		for(node=0; node<2*dataSetup.numSamples-1; node++) {
			//			printLocusGenTree(dataState.lociData[gen], nodePops[gen], nodeEvents[gen]);;
			//			printEventChains(gen);

			// tree root is NOT skipped - the root can be below root pop in the presence of migration
			// in this case, we might want to introduce/remove migration events from the edge above the root
			if(node == getLocusRoot(dataState.lociData[gen]))		continue;

			father = getNodeFather(dataState.lociData[gen], node);
			// record for later
			// UNUSED      t_old = getNodeAge(dataState.lociData[gen], father);
			father_pop_old = nodePops[gen][father];
			sibling = getNodeSon(dataState.lociData[gen], father, 0) + getNodeSon(dataState.lociData[gen], father, 1) - node;

			// trace original lineage and collect delta-stats from pruned version
			// of genetree to complete (original) version of genetree.
			// in genetree_stats_delta[0].

#ifdef LOG_STEPS
			fprintf(ioSetup.debugFile, "  gen %d, node %d, detaching father %d, pop %d, ",gen, node, father, father_pop_old);
#endif
			traceLineage(gen,node,0);

			// trace new lineage from node until reconnected
			// this also re-grafts edge onto tree.
			// collect delta-stats from pruned version to new version in genetree_stats_delta[1].
			// if res < 0, then could not reconnect (due to too many migrations)

			// if node corresponds to admixed sample, resample population assignment
			altPop = 0;
			admixSwitch = 0;
			if(node<dataSetup.numSamples && admixed_samples.number > 0 && admixed_samples.index[node] >= 0) {
				admixIndex = admixed_samples.index[node];
				oldPop = nodePops[gen][node];
				// consider alternative population w.p. admixture_status.admixtureCoefficients[node]
				if(rndu() < admixture_status.admixtureCoefficients[admixIndex]) {
					altPop = 1;
				}
				newPop = admixed_samples.popPairs[admixIndex][altPop];
				admixSwitch = (oldPop != newPop);
				nodePops[gen][node] = newPop;
				if(admixSwitch) {
					//			fprintf(ioSetup.debugFile, "Proposing to switch pop of sample %d in gen %d, from pop %d to pop %d.",node+1, gen+1, oldPop, newPop);
				}
			}

			res = traceLineage(gen,node,1);

#ifdef LOG_STEPS
			fprintf(ioSetup.debugFile, "proposing reattachment at edge above node %d, pop %d, age %g. ",
					mig_spr_stats.target, mig_spr_stats.father_pop_new, getNodeAge(dataState.lociData[gen], father));
#endif
			lnLd = -getLocusDataLikelihood(dataState.lociData[gen]);
			lnLd += computeLocusDataLikelihood(dataState.lociData[gen], /*reuse old conidtionals*/ 1);
			lnacceptance = lnLd;

#ifdef LOG_STEPS
			fprintf(ioSetup.debugFile, "lnacceptance = %g, ",lnacceptance);
#endif

			if(res >= 0 &&(lnacceptance>=0 || rndu()<exp(lnacceptance))) {
#ifdef LOG_STEPS
				fprintf(ioSetup.debugFile, "accepting.\n");
#endif
				accepted++;
				// UNUSED        didAccept = 1;
				genLogLikelihood[gen] += (mig_spr_stats.genetree_delta_lnLd[1] - mig_spr_stats.genetree_delta_lnLd[0]);
				dataState.dataLogLikelihood += lnLd;
				dataState.logLikelihood += (lnLd - mig_spr_stats.genetree_delta_lnLd[0] + mig_spr_stats.genetree_delta_lnLd[1])/dataSetup.numLoci;

				if(admixSwitch) {
					//			fprintf(ioSetup.debugFile, " Accepting with delta log likelihood %lf.\n",lnacceptance);
					genLogLikelihood[gen] +=   log(1/admixture_status.admixtureCoefficients[admixIndex] - 1) * (1 - 2*altPop);
					dataState.logLikelihood += log(1/admixture_status.admixtureCoefficients[admixIndex] - 1) * (1 - 2*altPop) / dataSetup.numLoci;
				}

				// change pointers of migration nodes to genetree edges
				// note that id of father might have changed (because of root)
				target = mig_spr_stats.target;
				t_new = getNodeAge(dataState.lociData[gen], father);
				for(i=0; i<genetree_migs[gen].num_migs; i++) {
					mig = genetree_migs[gen].living_mignodes[i];
					if(genetree_migs[gen].mignodes[mig].gtree_branch == father) {
						genetree_migs[gen].mignodes[mig].gtree_branch = sibling;
						/*						printf("\nSwitching node id's for mignode %d in gen %d from %d to %d due to MIG_SPR on node %d.\n",
                                    mig,gen,father_old,sib, node);
                                    printf("\nAttaching node %d to target edge %d at time %g.\n", node, target, t_new);
                                    printGtree(1);		*/
					}
					// this is to transfer "sibling" migrations into "father"
					// in case where target == father.
					if(target == father)		target = sibling;
					if(genetree_migs[gen].mignodes[mig].gtree_branch == target && genetree_migs[gen].mignodes[mig].age >= t_new) {
						genetree_migs[gen].mignodes[mig].gtree_branch = father;
						/*						printf("\nSwitching node id's for mignode %d in gen %d from %d to %d due to MIG_SPR on node %d.\n",
                                    mig,gen,target,nodes[node].father, node);
                                    printf("\nAttaching node %d to target edge %d at time %g.\n", node, target, t_new);
                                    printGtree(1);		*/
					}
				}

				// DO WE NEED TO DO ANYTHING IN CASE OF ROOT SWAP ??? !!
				/*
          event = nodes[tree.root].event_id;
          i = event_chains[gen].events[event].node_id;
          if(i != tree.root) {
          //printf("\n\n Root Swap.\n");
          event_chains[gen].events[event].node_id = tree.root;
          event = nodes[i].event_id;
          event_chains[gen].events[event].node_id = i;
          misc_stats.root_swaps++;
          }
				 */
				// remove old coalescent event and configure new one
				removeEvent(gen, mig_spr_stats.father_event_old);
				event_chains[gen].events[mig_spr_stats.father_event_new].type = COAL;
				event_chains[gen].events[mig_spr_stats.father_event_new].node_id = father;
				nodeEvents[gen][father] = mig_spr_stats.father_event_new;

				// update number of coalescences for new and old father populations (if necessary)
				if(mig_spr_stats.father_pop_new != father_pop_old) {
					nodePops[gen][father] = mig_spr_stats.father_pop_new;
					genetree_stats[gen].num_coals[father_pop_old] --;
					genetree_stats_total.num_coals[father_pop_old] --;
					genetree_stats[gen].num_coals[mig_spr_stats.father_pop_new] ++;
					genetree_stats_total.num_coals[mig_spr_stats.father_pop_new] ++;
				}


				// configure new migration nodes/events
				// remove old nodes/events
				replaceMigNodes(gen, node);

				// add lineage to all events in new path to new father
				for(i=0; i<genetree_stats_delta[1].num_changed_events; i++) {
					event = genetree_stats_delta[1].changed_events[i];
					event_chains[gen].events[event].num_lineages++;
				}

				// apply changes to genetree stats
				// changes in number of coals and migs are applied separately
				for(mig_band=0; mig_band<dataSetup.popTree->numMigBands; mig_band++) {
					genetree_stats[gen].mig_stats[mig_band] += (genetree_stats_delta[1].mig_stats_delta[mig_band] - genetree_stats_delta[0].mig_stats_delta[mig_band]);
					genetree_stats_total.mig_stats[mig_band]  += (genetree_stats_delta[1].mig_stats_delta[mig_band] - genetree_stats_delta[0].mig_stats_delta[mig_band]);
				}

				for(pop=0; pop<dataSetup.popTree->numPops; pop++) {
					genetree_stats[gen].coal_stats[pop] += genetree_stats_delta[1].coal_stats_delta[pop]-genetree_stats_delta[0].coal_stats_delta[pop];
					genetree_stats_total.coal_stats[pop]  += (genetree_stats_delta[1].coal_stats_delta[pop]-genetree_stats_delta[0].coal_stats_delta[pop])/heredity_factor;
				}


				resetSaved(dataState.lociData[gen]);
				/*
          genetree_lnLd_new = gtreeLnLikelihood(gen);
          misc_stats.spr_lnld_disc = max2(misc_stats.spr_lnld_disc , 
          fabs((genetree_lnLd_new - genetree_lnLd) -(mig_spr_stats.genetree_delta_lnLd[1] - mig_spr_stats.genetree_delta_lnLd[0])));
          genetree_lnLd = genetree_lnLd_new;
				 */
			}
			else {
				// UNUSED        didAccept = 0;
#ifdef LOG_STEPS
				fprintf(ioSetup.debugFile, "rejecting.\n");
#endif
				// remove all added migration events
				if(res >= 0) {
					removeEvent(gen, mig_spr_stats.father_event_new);
				}
				for(i=0; i<mig_spr_stats.num_new_migs; i++) {
					removeEvent(gen, mig_spr_stats.new_migs_in[i]);
					removeEvent(gen, mig_spr_stats.new_migs_out[i]);
				}
				// return reduced lineage to all events of original edge
				for(i=0; i<genetree_stats_delta[0].num_changed_events; i++) {
					event = genetree_stats_delta[0].changed_events[i];
					event_chains[gen].events[event].num_lineages++;
				}
				// change back population assignment (due to admixture)
				if(admixSwitch) {
					nodePops[gen][node] = oldPop;
					//			fprintf(ioSetup.debugFile, " Rejecting with delta log likelihood %lf.\n",lnacceptance);
				}
				revertToSaved(dataState.lociData[gen]);
			}
			/**
         if(!checkLocusDataLikelihood(dataState.lociData[gen])) {
        	 printf("\nError checking recorded likelihood for gen %d!", gen);
        	 printf("\n  --  Aborting after UpdateGB_MigSPR(%d,%d).\n",gen,node);
        	 printf("\nProposing SPR of subtree rooted at node %d, father %d (age %f-->%f), sibling %d --> target %dn and %s.\n", 
        	 node, father, t_old, t_new, sibling, target, (didAccept)?("accepting"):("rejecting"));
         	exit(-1);
         	res = 0;
         }

          if (!checkAll()) {
          	printf("\n  --  Aborting after UpdateGB_MigSPR(%d,%d).\n",gen,node);
          	printf("\nProposing SPR of subtree rooted at node %d, father %d (age %f-->%f), sibling %d --> target %d.\n", 
         	node, father, t_old, t_new, sibling, target);
         	exit(-1);
          }
			 **/
		} // end for(node)
	} // end for(gen)

	return(accepted);
}
/** end of UpdateGB_MigSPR **/



/***********************************************************************************
 *	UpdateAdmixCoeffs
 *	- perturbs admixture coefficients for all samples
 *	- this step does not affect data likelihood
 *	- this step doe not affect any other recorded statistic as well
 ***********************************************************************************/
int UpdateAdmixCoeffs (double finetune)	{
	int sample, node, gen, accepted=0;
	double coeffOld, coeffNew, logCoeffRatio, logCompCoeffRatio, lnacceptance;

	double deltaLogLikelihood;

	if(finetune <= 0.0) {
		return 0;
	}

	// first record the admixture counts for likelihood computation
	recordAdmixtureCounts();

	for(sample=0; sample<admixed_samples.number; sample++) {
		// record previous coefficient
		coeffOld = admixture_status.admixtureCoefficients[sample];
		// perform MULTIPLICATIVE UPDATE
		coeffNew = coeffOld + finetune*rnd2normal8();
		coeffNew = reflect(coeffNew, 0, 1);

		logCoeffRatio = log(coeffNew/coeffOld);
		logCompCoeffRatio = log((1-coeffNew)/(1-coeffOld));

		deltaLogLikelihood = admixture_status.admixtureCounts[sample] * logCoeffRatio + (dataSetup.numLoci - admixture_status.admixtureCounts[sample]) * logCompCoeffRatio;

		// acceptance ratio according to delta log likelihood
		lnacceptance = deltaLogLikelihood;


#ifdef LOG_STEPS
		fprintf(ioSetup.debugFile, "lnacceptance = %g, ",lnacceptance);
#endif
		if(lnacceptance>=0 || rndu()<exp(lnacceptance)) {
#ifdef LOG_STEPS
			fprintf(ioSetup.debugFile, "accepting.\n");
#endif
			accepted++;
			node = admixed_samples.samples[sample];
			for(gen=0; gen<dataSetup.numLoci; gen++) {
				if(nodePops[gen][node] == admixed_samples.popPairs[sample][0]) {
					genLogLikelihood[gen] += logCompCoeffRatio;
				} else {
					genLogLikelihood[gen] += logCoeffRatio;
				}
			}
			dataState.logLikelihood += deltaLogLikelihood/dataSetup.numLoci;
			admixture_status.admixtureCoefficients[sample] = coeffNew;
		} else {
#ifdef LOG_STEPS
			fprintf(ioSetup.debugFile, "rejecting.\n");
#endif
		}

	}// end of for(pop)

	return(accepted);
}
/** end of UpdateAdmixCoeffs **/



/***********************************************************************************
 *	UpdateTheta
 *	- perturbs all population thetas
 *	- this step does not affect data likelihood
 *	- this step does not affect any other recorded statistic as well
 ***********************************************************************************/
int UpdateTheta (double finetune)	{
	int pop, gen, accepted=0;
	double thetaold, thetanew, c, lnc, lnacceptance;

	double deltaLogLikelihood;

	if(finetune <= 0.0) {
		return 0;
	}

	for(pop=0; pop<dataSetup.popTree->numPops; pop++) {
		// record previous theta
		thetaold = dataSetup.popTree->pops[pop]->theta;
		// perform MULTIPLICATIVE UPDATE
		lnc = finetune*rnd2normal8();
		c = exp(lnc);
		thetanew = thetaold * c;

		#ifdef LOG_STEPS
				fprintf(ioSetup.debugFile, "  pop %d, proposing theta shift: %g-->%g, ",pop, thetaold,thetanew);
		#endif

		// acceptance ratio according to proposal prior
		lnacceptance = lnc + lnc * (dataSetup.popTree->pops[pop]->thetaPrior.alpha - 1) - (thetanew-thetaold) * dataSetup.popTree->pops[pop]->thetaPrior.beta;

		// delta in log-likelihood of all genealogies
		deltaLogLikelihood =
				-( 			lnc				* genetree_stats_total.num_coals[pop]	+
						(1/thetanew - 1/thetaold) 	* genetree_stats_total.coal_stats[pop] 	);

		lnacceptance += deltaLogLikelihood;

		#ifdef LOG_STEPS
				fprintf(ioSetup.debugFile, "lnacceptance = %g, ",lnacceptance);
		#endif

		if (lnacceptance>=0 || rndu()<exp(lnacceptance)) {
			#ifdef LOG_STEPS
						fprintf(ioSetup.debugFile, "accepting.\n");
			#endif
			accepted++;
			for(gen=0; gen<dataSetup.numLoci; gen++) {
				genLogLikelihood[gen] -= ( lnc * genetree_stats[gen].num_coals[pop] + (1/thetanew - 1/thetaold) * genetree_stats[gen].coal_stats[pop] );
			}
			dataState.logLikelihood += deltaLogLikelihood/dataSetup.numLoci;
			dataSetup.popTree->pops[pop]->theta = thetanew;
		} else {
		#ifdef LOG_STEPS
					fprintf(ioSetup.debugFile, "rejecting.\n");
		#endif
		}

	}// end of for(pop)

	return(accepted);
}
/** end of UpdateTheta **/



/***********************************************************************************
 *	UpdateMigRates
 *	- perturbs all migration rates (for all bands)
 *	- this step does not affect data likelihood
 *	- this step does not affect any other recorded statistic as well
 ***********************************************************************************/
int UpdateMigRates (double finetune) {

	int mig_band, gen, accepted=0;
	double old_rate, new_rate, c,lnc, lnacceptance, deltaLogLikelihood;

	if(finetune <= 0.0) {
		return 0;
	}

	for(mig_band=0; mig_band<dataSetup.popTree->numMigBands; mig_band++) {


		// record old rate and propose new migration rate
		old_rate = dataSetup.popTree->migBands[mig_band].migRate;
		// ADDITIVE UPDATE !!!
		//		new_rate = fabs(old_rate + finetune*rnd2normal8());
		//		c = new_rate/old_rate;
		//		lnc = log(new_rate/old_rate);

		// perform MULTIPLICATIVE UPDATE
		lnc = finetune*rnd2normal8();
		c = exp(lnc);
		new_rate = old_rate * c;
#ifdef LOG_STEPS
		fprintf(ioSetup.debugFile, "  mig-band %d, proposing rate shift: %g-->%g, ",mig_band, old_rate, new_rate);
#endif
		//		if(new_rate < 0.0000001 || new_rate > MAX_MIG_RATE) {
#ifdef LOG_STEPS
		//			fprintf(ioSetup.debugFile, "rejecting - out of bounds.\n");
#endif
		//			continue;
		//		}

		//   	printf("\nUpdating migration rate for band %d (%d-->%d): %g-->%g.",
		//  				mig_band, dataSetup.popTree->migBands[mig_band].sourcePop, dataSetup.popTree->migBands[mig_band].targetPop, old_rate, new_rate);

		// GAMA PRIOR
		if(new_rate < 0.00001)		continue;
		lnacceptance = lnc + lnc * 	(dataSetup.popTree->migBands[mig_band].migRatePrior.alpha - 1) -
				(new_rate-old_rate) *dataSetup.popTree->migBands[mig_band].migRatePrior.beta;

		// UNIFORM PRIOR
		// if new rate is outside boundary, reject
		//		if(new_rate < 0.00001 || new_rate > dataSetup.popTree->migBands[mig_band].upperBound)	continue;
		//		lnacceptance = lnc;

		// diff in log-likelihood
		deltaLogLikelihood =
				( 			lnc			* genetree_stats_total.num_migs[mig_band]	-
						(new_rate - old_rate) 	* genetree_stats_total.mig_stats[mig_band]   );

		lnacceptance += deltaLogLikelihood;

#ifdef LOG_STEPS
		fprintf(ioSetup.debugFile, "lnacceptance = %g, ",lnacceptance);
#endif

		if(lnacceptance>=0 || rndu()<exp(lnacceptance)) {
#ifdef LOG_STEPS
			fprintf(ioSetup.debugFile, "accepting.\n");
#endif
			accepted++;
			for(gen=0; gen<dataSetup.numLoci; gen++) {
				genLogLikelihood[gen] += ( lnc * genetree_stats[gen].num_migs[mig_band] - (new_rate - old_rate) * genetree_stats[gen].mig_stats[mig_band] );
			}
			dataSetup.popTree->migBands[mig_band].migRate = new_rate;
			dataState.logLikelihood += deltaLogLikelihood/dataSetup.numLoci;
		} else {
#ifdef LOG_STEPS
			fprintf(ioSetup.debugFile, "rejecting.\n");
#endif
		}
	}// end of for(mig_band)
	return(accepted);
}
/** end of UpdateMigRates **/




/***********************************************************************************
 *	UpdateTau
 *	- perturbs ancestral population ages
 *	- performs rubber band operation on ancestral population and two daughter populations
 *	- checks for conflicts - due to migration bands and/or migration nodes
 *	- this step affects data likelihood
 *	- this step also affects all the recorded statistics
 ***********************************************************************************/
void UpdateTau (double *finetunes, int *accepted)	{

	int k, ancestralPop, inode, gen;
	int ntj[2], ntj_gen[2];
	double tauold, taunew, taub[2], taufactor[2];
	double lnacceptance=0; //, lnLd;

	int i, mig, mig1, pop=-1, mig_conflict, event;
	int num_affected_mig_bands, affected_mig_bands[MAX_MIG_BANDS], start_or_end[MAX_MIG_BANDS];
	double new_band_ages[MAX_MIG_BANDS];
	int mig_band, ntj_gen1[2];
	double new_age, age;
	double dataDeltaLnLd, genDeltaLnLd;

	int sourcePop, targetPop, fatherNode, sons[2];
	unsigned short inORout=-1;			// for potentially conflicting migration events
	unsigned short isRoot = 0;
	unsigned short res;
	//UNUSED unsigned short didAccept = 0;


	for(ancestralPop=dataSetup.popTree->numCurPops; ancestralPop<dataSetup.popTree->numPops; ancestralPop++) {
		accepted[ancestralPop] = 0;
		isRoot = (ancestralPop == dataSetup.popTree->rootPop);

		tauold = dataSetup.popTree->pops[ancestralPop]->age;
		dataDeltaLnLd = 0.0;
		genDeltaLnLd  = 0.0;
		sons[0] = dataSetup.popTree->pops[ancestralPop]->sons[0]->id;
		sons[1] = dataSetup.popTree->pops[ancestralPop]->sons[1]->id;
		taub[0] = max2(dataSetup.popTree->pops[sons[0]]->age , dataSetup.popTree->pops[sons[1]]->age);
		// MARK CHANGE
		taub[0] = max2(taub[0], dataSetup.popTree->pops[sons[0]]->sampleAge);
		taub[0] = max2(taub[0], dataSetup.popTree->pops[sons[1]]->sampleAge);
		if(isRoot) {
			taub[1] = OLDAGE;
		} else {
			taub[1] = dataSetup.popTree->pops[ancestralPop]->father->age;
		}


		// modify bounds according to migration bands - make sure all migration bands stay alive
		for(mig_band=0; mig_band<dataSetup.popTree->numMigBands; mig_band++) {
			sourcePop = dataSetup.popTree->migBands[mig_band].sourcePop;
			targetPop = dataSetup.popTree->migBands[mig_band].targetPop;
			if(sourcePop == ancestralPop || targetPop == ancestralPop) {
				taub[1] = min2(taub[1], dataSetup.popTree->migBands[mig_band].endTime);
			} else if( sourcePop == sons[0] || sourcePop == sons[1] || targetPop == sons[0] || targetPop == sons[1] ) {
				taub[0] = max2(taub[0], dataSetup.popTree->migBands[mig_band].startTime);
			}
		}

		//sample new time
		taunew = tauold + finetunes[ancestralPop]*rnd2normal8();
		taunew = reflect(taunew,taub[0],taub[1]);
		// set new age for now to compute migration band times, but restore later !!!
		dataSetup.popTree->pops[ancestralPop]->age = taunew;
#ifdef LOG_STEPS
		fprintf(ioSetup.debugFile, "  ancestral pop %d, proposing age shift: %g-->%g, ",ancestralPop, tauold, taunew);
#endif      	
		// set rubberband factors
		for(k=0; k<2; k++)
			taufactor[k]=(taunew-taub[k])/(tauold-taub[k]);

		// events in root population are inversely scaled (for some reason)
		if(isRoot) {
			taufactor[1] = taufactor[0];
		}

		// compute new times for mig bands and find affected migration bands
		//		printf("\nAffected migration bands, when adjusting ancestral pop %d split %f --> %f (taub[0] = %g, taub[1] = %g, taufactor[0] = %g, taufactor[1] = %g:\n",
		//					ancestralPop, tauold, taunew,taub[0], taub[1],taufactor[0], taufactor[1]);
		num_affected_mig_bands = 0;
		for(mig_band=0; mig_band<dataSetup.popTree->numMigBands; mig_band++) {
			sourcePop = dataSetup.popTree->migBands[mig_band].sourcePop;
			targetPop = dataSetup.popTree->migBands[mig_band].targetPop;
			res = updateMigrationBandTimes(dataSetup.popTree,mig_band);

			if((sourcePop == sons[0] && targetPop == sons[1]) || (sourcePop == sons[1] && targetPop == sons[0])) {
				// migration bands between son populations are actually not affected,
				// but this makes future conditions simpler

			} else if(targetPop == ancestralPop) {
				// mig bands entering rubber-banded populations are not affected
				// by the standard rubber band. We factor the times artificially
				// so that after rubber-band they will be in the right spot
				if(dataSetup.popTree->migBands[mig_band].endTime < taub[1]) {
					//					printf("    mig band %d, type 1a.\n",mig_band);
					affected_mig_bands[num_affected_mig_bands] = mig_band;
					start_or_end[num_affected_mig_bands] = 0;		//indicate that end time has changed
					new_band_ages[num_affected_mig_bands] = taub[1] + ( dataSetup.popTree->migBands[mig_band].endTime - taub[1] ) / taufactor[1] ;
					num_affected_mig_bands++;
				}
				// do same with start time (if not bounded before and after by age of ancestralPop)
				if(dataSetup.popTree->migBands[mig_band].startTime < taub[1] && dataSetup.popTree->pops[sourcePop]->age > min2(tauold,taunew)) {
					//					printf("    mig band %d, type 1b.\n",mig_band);
					affected_mig_bands[num_affected_mig_bands] = mig_band;
					start_or_end[num_affected_mig_bands] = 1;		//indicate that start time has changed
					new_band_ages[num_affected_mig_bands] = taub[1] + ( dataSetup.popTree->migBands[mig_band].startTime - taub[1] ) / taufactor[1] ;
					num_affected_mig_bands++;
				}
			} else if(targetPop  == sons[0] || targetPop  == sons[1]) {
				// same idea as previous condition, but with lower rubberband
				if(dataSetup.popTree->migBands[mig_band].startTime > taub[0]) {
					//					printf("    mig band %d, type 2a.\n",mig_band);
					affected_mig_bands[num_affected_mig_bands] = mig_band;
					start_or_end[num_affected_mig_bands] = 1;		//indicate that start time has changed
					new_band_ages[num_affected_mig_bands] = taub[0] + ( dataSetup.popTree->migBands[mig_band].startTime - taub[0] ) / taufactor[0] ;
					num_affected_mig_bands++;
				}
				// do same with start time (if not bounded before and after by age of ancestralPop)
				if(dataSetup.popTree->migBands[mig_band].endTime > taub[0] && dataSetup.popTree->pops[sourcePop]->father->age < max2(tauold,taunew)) {
					affected_mig_bands[num_affected_mig_bands] = mig_band;
					start_or_end[num_affected_mig_bands] = 0;		//indicate that end time has changed
					new_band_ages[num_affected_mig_bands] = taub[0] + ( dataSetup.popTree->migBands[mig_band].endTime - taub[0] ) / taufactor[0] ;
					num_affected_mig_bands++;
					//					printf("    mig band %d, type 2b. New end time = %g, so end time is artificially set to %g\n",
					//						   mig_band, dataSetup.popTree->migBands[mig_band].endTime, new_band_ages[num_affected_mig_bands-1]);
				}
			} else if(res && sourcePop == ancestralPop) {
				// start time of this migration band changes with proposed change in tau
				//				printf("    mig band %d, type 3.\n",mig_band);
				affected_mig_bands[num_affected_mig_bands] = mig_band;
				start_or_end[num_affected_mig_bands] = 1;		//indicate that start time has changed
				new_band_ages[num_affected_mig_bands] = dataSetup.popTree->migBands[mig_band].startTime;
				num_affected_mig_bands++;
			} else if(res && (sourcePop  == sons[0] || sourcePop  == sons[1])) {
				//				printf("    mig band %d, type 4.\n",mig_band);
				// end time of this migration band changes with proposed change in tau
				affected_mig_bands[num_affected_mig_bands] = mig_band;
				start_or_end[num_affected_mig_bands] = 0;		//indicate that end time has changed
				new_band_ages[num_affected_mig_bands] = dataSetup.popTree->migBands[mig_band].endTime;
				num_affected_mig_bands++;
			}
			//			else {
			//				printf("not affected.\n");
			//			}




			/*
      //			printf("mig band %d (%d-->%d) - ",mig_band, sourcePop,targetPop);
      if(	sourcePop == ancestralPop && res) {
      // start time of this migration band changes with proposed change in tau
      //				printf("above split.\n",mig_band);
      affected_mig_bands[num_affected_mig_bands] = mig_band;
      start_or_end[num_affected_mig_bands] = 1;		//indicate that start time has changed
      new_band_ages[num_affected_mig_bands] = max2(taunew, dataSetup.popTree->pops[targetPop]->age);
      num_affected_mig_bands++;
      } else if( 	(sourcePop  == sons[0] || sourcePop  == sons[1]) && 
      (targetPop  != sons[0] && targetPop  != sons[1]) && dataSetup.popTree->pops[targetPop]->father->age > min2(taunew,tauold)) {
      //				printf("below split.\n",mig_band);
      // end time of this migration band changes with proposed change in tau
      affected_mig_bands[num_affected_mig_bands] = mig_band;
      start_or_end[num_affected_mig_bands] = 0;		//indicate that end time has changed
      new_band_ages[num_affected_mig_bands] = min2(dataSetup.popTree->pops[targetPop]->father->age , taunew);
      num_affected_mig_bands++;
      } 
      //			else {
      //				printf("not affected.\n");
      //			}
			 */
		} // end of for(mig_band)

		// restoring old time for various computations - DIRTY !!!!
		dataSetup.popTree->pops[ancestralPop]->age = tauold;


		//      printf("\nConsidering UpdateTau for split of pop %d: %g-->%g. Num affected migration bands is %d.",ancestralPop,tauold,taunew, num_affected_mig_bands);

#ifdef DEBUG_RUBBERBAND
		printf("Performing rubber band on split %d: times %g --> %g. Upper/lower bounds - %f / %f, factors: %f / %f.\n",
				ancestralPop, tauold,taunew,taub[0],taub[1],taufactor[0],taufactor[1]);
#endif				

		lnacceptance = 	log(taunew/tauold) * (dataSetup.popTree->pops[ancestralPop]->agePrior.alpha - 1) -
				(taunew-tauold)   *  dataSetup.popTree->pops[ancestralPop]->agePrior.beta;

		dataDeltaLnLd = 0.0;
		genDeltaLnLd  = 0.0;

		// initialize -  no migration conflicts, and number of moved nodes
		mig_conflict = 0;
		ntj[0]=ntj[1]=0;
		// implement rubberband on all gen genealogies
		for(gen=0; gen<dataSetup.numLoci; gen++) {

#ifdef CHECK_OPERATIONS
			ntj_gen[0] = ntj_gen[1] = 0;
			for(inode=dataSetup.numSamples; inode<2*dataSetup.numSamples-1; inode++) {
				t=getNodeAge(dataState.lociData[gen], inode);
				if (t>=taub[0] && t<taub[1] &&
						(nodePops[gen][inode] == ancestralPop || nodePops[gen][inode] == sons[0] || nodePops[gen][inode] == sons[1])) {
					k = (t>=tauold && !isRoot); /* k=0: below; 1: above */
					//					k = (t>=tauold); /* k=0: below; 1: above */
					ntj_gen[k]++;
				}
			}
#endif

			// deal with rubber-banded migration nodes, and their representation in
			// non-rubberbanded populations


			rubberband_migs[gen].num_moved_events = 0;
			event = -1;
			new_age = 0.0;
			ntj_gen1[0] = ntj_gen1[1] = 0;
			for(i=0; i<genetree_migs[gen].num_migs; i++) {
				mig = genetree_migs[gen].living_mignodes[i];
				mig_band  = genetree_migs[gen].mignodes[mig].migration_band;
				sourcePop = genetree_migs[gen].mignodes[mig].source_pop;
				targetPop = genetree_migs[gen].mignodes[mig].target_pop;
				age       = genetree_migs[gen].mignodes[mig].age;

				if(age < taub[0] || age > taub[1])		continue;

				// we assume here that there are no in/out migrations from root population.
				if((sourcePop == sons[0] && targetPop == sons[1]) || (sourcePop == sons[1] && targetPop == sons[0])) {
					// migration bands between son populations are actually not affected,
					// but this makes future conditions simpler
					ntj_gen1[0]++;
				} else if(sourcePop == ancestralPop) {
					inORout = 1;			// indicating out migration
					event = genetree_migs[gen].mignodes[mig].target_event;
					pop = targetPop;
					new_age = taub[1] + taufactor[1]*(age-taub[1]);
					// rubberBand only counts migrations coming into pops
					ntj_gen1[1]++;
				} else if(targetPop == ancestralPop) {
					inORout = 0;			// indicating in migration
					event = genetree_migs[gen].mignodes[mig].source_event;
					pop = sourcePop;
					new_age = taub[1] + taufactor[1]*(age-taub[1]);
					ntj_gen1[1]++;
				} else if(	(sourcePop == sons[0] || sourcePop == sons[1]) && genetree_migs[gen].mignodes[mig].age > taub[0] ) {
					inORout = 1;			// indicating out migration
					event = genetree_migs[gen].mignodes[mig].target_event;
					pop = targetPop;
					new_age = taub[0] + taufactor[0]*(age-taub[0]);
					// rubberBand only counts migrations coming into pops
					ntj_gen1[0]++;
				} else if(	(targetPop == sons[0] || targetPop == sons[1]) && genetree_migs[gen].mignodes[mig].age > taub[0] ) {
					inORout = 0;			// indicating in migration
					event = genetree_migs[gen].mignodes[mig].source_event;
					pop = sourcePop;
					new_age = taub[0] + taufactor[0]*(age-taub[0]);
					ntj_gen1[0]++;
				}


				if(event >= 0) {
					inode = genetree_migs[gen].mignodes[mig].gtree_branch;
					// check for conflicts
					if(new_age >= dataSetup.popTree->migBands[mig_band].endTime) {
						mig_conflict = 1;
						//						fprintf(ioSetup.debugMiscFile, "Mig conflict of type 1, mig-band %d, end time %g. ",mig_band, dataSetup.popTree->migBands[mig_band].endTime);
					} else if(new_age <= dataSetup.popTree->migBands[mig_band].startTime) {
						mig_conflict = 1;
						//						fprintf(ioSetup.debugMiscFile, "Mig conflict of type 4, mig-band %d, start time %g. ",mig_band, dataSetup.popTree->migBands[mig_band].startTime);
					} else if(inORout == 0 && new_age > age) {
						// an incoming migration event can conflict with event directly above it
						fatherNode = getNodeFather(dataState.lociData[gen], inode);
						mig1 = findFirstMig(gen,inode,genetree_migs[gen].mignodes[mig].age);
						if(mig1 >=0 &&
								genetree_migs[gen].mignodes[mig1].source_pop != ancestralPop &&
								genetree_migs[gen].mignodes[mig1].source_pop != sons[0] &&
								genetree_migs[gen].mignodes[mig1].source_pop != sons[1] &&
								new_age >= genetree_migs[gen].mignodes[mig1].age) {
							mig_conflict = 1;
							//							fprintf(ioSetup.debugMiscFile, "Mig conflict of type 2, mig-node %d, age %g." ,mig1, genetree_migs[gen].mignodes[mig1].age);
						} else if(fatherNode >= 0 && new_age >= getNodeAge(dataState.lociData[gen], fatherNode)) {
							mig_conflict = 1;
							//							fprintf(ioSetup.debugMiscFile, "Mig conflict of type 3, father %d, pop %d, age %g. ",fatherNode, nodePops[gen][fatherNode], getNodeAge(dataState.lociData[gen], fatherNode));
						}
					} else if(inORout == 1 && new_age < age) {
						// outgoing migration events can conflict with event directly below it
						mig1 = findLastMig(gen,inode,genetree_migs[gen].mignodes[mig].age);
						if(mig1 >=0 &&
								genetree_migs[gen].mignodes[mig1].target_pop != ancestralPop &&
								genetree_migs[gen].mignodes[mig1].target_pop != sons[0] &&
								genetree_migs[gen].mignodes[mig1].target_pop != sons[1] &&
								new_age <= genetree_migs[gen].mignodes[mig1].age) {
							mig_conflict = 1;
							//							fprintf(ioSetup.debugMiscFile, "Mig conflict of type 5, mig-node %d, age %g." ,mig1, genetree_migs[gen].mignodes[mig1].age);
						} else if(new_age <=getNodeAge(dataState.lociData[gen], inode)) {
							mig_conflict = 1;
							//							fprintf(ioSetup.debugMiscFile, "Mig conflict of type 6, node %d, pop %d, age %g, ",inode, nodePops[gen][inode], getNodeAge(dataState.lociData[gen], inode));
						}
					}
					if(mig_conflict)	break;

					//					printf("\n- Locus %d: moving migration event %d from age %g to age %g.",
					//						    gen, event, genetree_migs[gen].mignodes[mig].age, new_age);
					rubberband_migs[gen].orig_events[rubberband_migs[gen].num_moved_events] = event;
					rubberband_migs[gen].pops[rubberband_migs[gen].num_moved_events] = pop;
					rubberband_migs[gen].new_ages[rubberband_migs[gen].num_moved_events] = new_age;
					rubberband_migs[gen].num_moved_events++;
					event = -1;
				}
			} // end for(mignode)

			if(mig_conflict)  {
				//				printf("Mignodes requiring special updates:\n");
				//				for(i=0; i<rubberband_migs[gen].num_moved_events; i++) {
				//					printf("mignode %d, event %d, pop %d, new age = %f.\n", event_chains[gen].events[rubberband_migs[gen].orig_events[i]].node_id,
				//							rubberband_migs[gen].orig_events[i], rubberband_migs[gen].pops[i], rubberband_migs[gen].new_ges[i]);
				//				printEventChains(gen);

				//				revertToSaved(dataState.lociData[gen]);
				rubberband_migs[gen].num_moved_events = 0;
				break;
			}

			// create new events for affected migration bands
			for(i=0; i<num_affected_mig_bands; i++) {
				mig_band = affected_mig_bands[i];
				targetPop = dataSetup.popTree->migBands[mig_band].targetPop;
				for(event = event_chains[gen].first_event[targetPop]; event >= 0; event = event_chains[gen].events[event].next) {
					if(event_chains[gen].events[event].node_id == mig_band &&
							((start_or_end[i] && event_chains[gen].events[event].type == MIG_BAND_START) ||
									event_chains[gen].events[event].type == MIG_BAND_END  ) )     break;
				}
				if(event<0) {
					if(debug) {
						fprintf(stderr, "\nError: UpdateTau: couldn't find event for migration band %d in gen %d.\n", mig_band,gen);
					} else {
						fprintf(stderr, "Fatal Error 0074.\n");
					}
					printGenealogyAndExit(gen, -1);
				}
				rubberband_migs[gen].orig_events[rubberband_migs[gen].num_moved_events] = event;
				rubberband_migs[gen].pops[rubberband_migs[gen].num_moved_events] = targetPop;
				rubberband_migs[gen].new_ages[rubberband_migs[gen].num_moved_events] = new_band_ages[i];
				rubberband_migs[gen].num_moved_events++;
			}

			// compute residual effects of rubber-band (before actual rubber-band
			genDeltaLogLikelihood[gen] = rubberBandRipple(gen, 1 /*do changes*/);

			// compute rubber band
			if(isRoot) {
				genDeltaLogLikelihood[gen]  += rubberBand(gen, ancestralPop, taub[0], tauold, taufactor[1], 0 /*don't change chain*/, &ntj_gen1[1]);
			} else {
				genDeltaLogLikelihood[gen]  += rubberBand(gen, ancestralPop, taub[1], tauold, taufactor[1], 0 /*don't change chain*/, &ntj_gen1[1]);
			}
			genDeltaLogLikelihood[gen] += rubberBand(gen, sons[0], taub[0], tauold, taufactor[0], 0 /*don't change chain*/, &ntj_gen1[0]);
			genDeltaLogLikelihood[gen] += rubberBand(gen, sons[1], taub[0], tauold, taufactor[0], 0 /*don't change chain*/, &ntj_gen1[0]);

			genDeltaLnLd += genDeltaLogLikelihood[gen];

#ifdef CHECK_OPERATIONS
			if(	(!isRoot && (ntj_gen[0] != ntj_gen1[0] || ntj_gen[1] != ntj_gen1[1])) ||
					( isRoot && (ntj_gen[0] != ntj_gen1[0]+ ntj_gen1[1] || ntj_gen[1] != 0))){
				fprintf(stderr, "Error: UpdateTau has incorrect computation of modified nodes for gen %d, split %d: below %d, %d  ; above %d, %d.\n",
						gen, ancestralPop, ntj_gen[0] , ntj_gen1[0] , ntj_gen[1] , ntj_gen1[1]);
				printGenealogyAndExit(gen, -1);
			}
			ntj_gen1[0] = ntj_gen[0];
			ntj_gen1[1] = ntj_gen[1];
#endif         
			ntj[0]+=ntj_gen1[0];
			ntj[1]+=ntj_gen1[1];

			if(ntj_gen1[0]+ntj_gen1[1]) {
				dataDeltaLnLd -= getLocusDataLikelihood(dataState.lociData[gen]);
				dataDeltaLnLd += computeLocusDataLikelihood(dataState.lociData[gen], /*reuse old conditionals*/ 1);
			}
		}// end for(gen) - genealogy updates by rubberband

		lnacceptance += dataDeltaLnLd + genDeltaLnLd + ntj[0]*log(taufactor[0]) + ntj[1]*log(taufactor[1]);
		//		lnacceptance += totalDeltaLnLd;

		//		fprintf(ioSetup.debugMiscFile, "Ancestral pop %d %g --> %g. Mig node age %g --> %g.\n",ancestralPop,tauold,taunew,age,new_age);


#ifdef LOG_STEPS
		if(mig_conflict) {
			fprintf(ioSetup.debugFile, "migration conflict at gen %d, ", gen);
		} else {
			fprintf(ioSetup.debugFile, "lnacceptance = %g, ",lnacceptance);
		}
#endif	
		//No migraton conflict   Positive acceptance    accept some even if acceptance below 0 based on random probability
		if(!mig_conflict && (lnacceptance>=0 || rndu()<exp(lnacceptance))) {

#ifdef LOG_STEPS
			fprintf(ioSetup.debugFile, "accepting.\n");
#endif
			// UNUSED      didAccept = 1;
			accepted[ancestralPop]++;
			dataState.dataLogLikelihood += dataDeltaLnLd;
			dataState.logLikelihood += (dataDeltaLnLd + genDeltaLnLd) / dataSetup.numLoci;
			for(gen=0; gen<dataSetup.numLoci; gen++) {
				genLogLikelihood[gen] += genDeltaLogLikelihood[gen];
				// change gene trees, event chains, and likelihoods
				if(isRoot) {
					rubberBand(gen, ancestralPop, taub[0], tauold, taufactor[1], 1 /*change chain*/, &ntj_gen1[1]);
				} else {
					rubberBand(gen, ancestralPop, taub[1], tauold, taufactor[1], 1 /*change chain*/, &ntj_gen1[1]);
				}
				rubberBand(gen, sons[0], taub[0], tauold, taufactor[0], 1 /*change chain*/, &ntj_gen[0]);
				rubberBand(gen, sons[1], taub[0], tauold, taufactor[0], 1 /*change chain*/, &ntj_gen[0]);

				// accept genealogy changes
				resetSaved(dataState.lociData[gen]);

				// remove original added events for migrations and migration bands
				for(i=0; i<rubberband_migs[gen].num_moved_events; i++) {
					// set pointers from mignodes to new events
					mig = event_chains[gen].events[rubberband_migs[gen].new_events[i]].node_id;
					if(event_chains[gen].events[rubberband_migs[gen].new_events[i]].type == IN_MIG) {
						genetree_migs[gen].mignodes[mig].target_event = rubberband_migs[gen].new_events[i];
						// adjust ages of mignodes for migrations out of rubberband
						// the ones coming in are adjusted in rubberBand.
						genetree_migs[gen].mignodes[mig].age = rubberband_migs[gen].new_ages[i];
					} else if(event_chains[gen].events[rubberband_migs[gen].new_events[i]].type == OUT_MIG) {
						genetree_migs[gen].mignodes[mig].source_event = rubberband_migs[gen].new_events[i];
					}
					removeEvent(gen,rubberband_migs[gen].orig_events[i]);
				}
				rubberband_migs[gen].num_moved_events = 0;
			}// end of for(gen) - implement genealogy changes

			// commit to new split time [ no need to change mig-band times - these were changed already ]
			dataSetup.popTree->pops[ancestralPop]->age = taunew;

			if(isRoot) {
				adjustRootEvents();
			}
		}
		else {
#ifdef LOG_STEPS
			fprintf(ioSetup.debugFile, "rejecting.\n");
#endif
			// UNUSED       didAccept = 0;
			// restore old times of migration bands
			computeMigrationBandTimes(dataSetup.popTree);
			if(mig_conflict) {
				//				printf("(migration conflict at gen %d)\n",gen);
				misc_stats.rubberband_mig_conflicts++;
			}
			else {
				gen = dataSetup.numLoci-1;
			}
			// start from gen before last and redo changes
			for( ; gen >= 0; --gen) {
				// redo changes in events for migrations and mig bands.
				revertToSaved(dataState.lociData[gen]);
				rubberBandRipple(gen, 0 /*redo changes*/);
			}
		}
		/**
           if (!checkAll()) {
           printf("\n  --  Aborting after UpdateTau for ancestral pop %d, accepted = %d.\n",ancestralPop, didAccept);
           if(didAccept) {
           printf("Updating from time %g to time %g.\n",tauold, taunew);
           }
           if(mig_conflict) printf("Found migration conflict.\n");
           exit(-1);
           }
           //		else {
           //printf("\n  UpdateTau for ancestral pop %d OK.",ancestralPop);
           // }

	for(gen=0; gen<dataSetup.numLoci; gen++) {
		if (!synchronizeEvents(gen)) {
			printf("\n  --  Aborting due to problems found when synchronizing data structures for locus #%d after UpdateTAU(pop=%d) %g-->%g (accepted=%d).\n\n",
				   gen+1, ancestralPop, tauold, taunew, didAccept);
			printGenealogyAndExit(gen,-1);
		}
	}
		 **/
	}// end of for(ancestralPop)
}
/** end of UpdateTau **/



/***********************************************************************************
 *	UpdateSampleAge
 *	- perturbs non-zero sample ages of sampled populations
 *	- performs rubber band operation on part of population above/below
 *	- checks for conflicts - due to migration bands and/or migration nodes
 *	- this step affects data likelihood
 *	- this step also affects all the recorded statistics (similar to updateTau)
 ***********************************************************************************/
void UpdateSampleAge (double *finetunes, int *accepted)	{

	int k, pop, inode, gen;
	int ntj[2], ntj_gen[2];
	double tauold, taunew, taub[2], taufactor[2];
	double lnacceptance=0; //, lnLd;

	int i, mig, mig1, migPop=-1, mig_conflict, event;
	int num_affected_mig_bands, affected_mig_bands[MAX_MIG_BANDS], start_or_end[MAX_MIG_BANDS];
	double new_band_ages[MAX_MIG_BANDS];
	int mig_band, ntj_gen1[2];
	double new_age, age;
	double dataDeltaLnLd, genDeltaLnLd;

	int sourcePop, targetPop, fatherNode;
	unsigned short inORout=-1;			// for potentially conflicting migration events
	// UNUSED    unsigned short didAccept = 0;


	for(pop=0; pop<dataSetup.popTree->numCurPops; pop++) {
		accepted[pop] = 0;
		if(!dataSetup.popTree->pops[pop]->updateSampleAge) continue;

		tauold = dataSetup.popTree->pops[pop]->sampleAge;
		dataDeltaLnLd = 0.0;
		genDeltaLnLd  = 0.0;
		taub[0] = 0.0;
		taub[1] = dataSetup.popTree->pops[pop]->father->age;


		//sample new time
		taunew = tauold + finetunes[pop]*rnd2normal8();
		taunew = reflect(taunew,taub[0],taub[1]);
		// set new age for now to compute migration band times, but restore later !!!
		dataSetup.popTree->pops[pop]->sampleAge = taunew;
#ifdef LOG_STEPS
		fprintf(ioSetup.debugFile, "  pop %d, proposing sample age shift: %g-->%g, ",pop, tauold, taunew);
#endif      	
		// set rubberband factors
		for(k=0; k<2; k++)
			taufactor[k]=(taunew-taub[k])/(tauold-taub[k]);


		num_affected_mig_bands = 0;
		for(mig_band=0; mig_band<dataSetup.popTree->numMigBands; mig_band++) {
			sourcePop = dataSetup.popTree->migBands[mig_band].sourcePop;
			targetPop = dataSetup.popTree->migBands[mig_band].targetPop;

			if(targetPop == pop) {
				// mig bands entering rubber-banded populations are not affected
				// by the standard rubber band. We factor the times artificially
				// so that after rubber-band they will be in the right spot
				if(dataSetup.popTree->migBands[mig_band].endTime < taub[1] && dataSetup.popTree->migBands[mig_band].endTime > taub[0]) {
					//					printf("    mig band %d, type 1a.\n",mig_band);
					affected_mig_bands[num_affected_mig_bands] = mig_band;
					start_or_end[num_affected_mig_bands] = 0;		//indicate that end time has changed
					age = dataSetup.popTree->migBands[mig_band].endTime;
					new_band_ages[num_affected_mig_bands] =  taub[age>taunew] + (age-taub[age>taunew])/taufactor[age>taunew];
					num_affected_mig_bands++;
				}
				// do same with start time
				if(dataSetup.popTree->migBands[mig_band].startTime < taub[1] && dataSetup.popTree->migBands[mig_band].startTime > taub[0]) {
					//					printf("    mig band %d, type 1b.\n",mig_band);
					affected_mig_bands[num_affected_mig_bands] = mig_band;
					start_or_end[num_affected_mig_bands] = 1;		//indicate that start time has changed
					age = dataSetup.popTree->migBands[mig_band].startTime;
					new_band_ages[num_affected_mig_bands] =  taub[age>taunew] + (age-taub[age>taunew])/taufactor[age>taunew];
					num_affected_mig_bands++;
				}
			}
			//else {
			//	 printf("not affected.\n");
			//}
		} // end of for(mig_band)


		// restoring old time for various computations - DIRTY !!!!
		dataSetup.popTree->pops[pop]->sampleAge = tauold;



#ifdef DEBUG_RUBBERBAND
		printf("Performing rubber band on pop %d for sample age: times %g --> %g. Upper/lower bounds - %f / %f, factors: %f / %f.\n",
				pop, tauold,taunew,taub[0],taub[1],taufactor[0],taufactor[1]);
#endif				


		//SAMPLEAGE: do we want to have a prior Gamma distribution associated with sample age?
		lnacceptance = 	log(taunew/tauold) * (dataSetup.popTree->pops[pop]->agePrior.alpha - 1) -
				(taunew-tauold)   *  dataSetup.popTree->pops[pop]->agePrior.beta;

		dataDeltaLnLd = 0.0;
		genDeltaLnLd  = 0.0;

		// initialize -  no migration conflicts, and number of moved nodes
		mig_conflict = 0;
		ntj[0]=ntj[1]=0;
		// implement rubberband on all gen genealogies
		for(gen=0; gen<dataSetup.numLoci; gen++) {

#ifdef CHECK_OPERATIONS
			ntj_gen[0] = ntj_gen[1] = 0;
			for(inode=dataSetup.numSamples; inode<2*dataSetup.numSamples-1; inode++) {
				t=getNodeAge(dataState.lociData[gen], inode);
				if (t>=taub[0] && t<taub[1] &&
						(nodePops[gen][inode] == pop)) {
					k = (t>=tauold); /* k=0: below; 1: above */
					//					k = (t>=tauold); /* k=0: below; 1: above */
					ntj_gen[k]++;
				}
			}
#endif

			// deal with rubber-banded migration nodes, and their representation in
			// non-rubberbanded populations


			// printPopulationTree(dataSetup.popTree, stderr, 1);
			// printLocusGenTree(dataState.lociData[gen], stderr, nodePops[gen], nodeEvents[gen]);
			// printEventChains(stderr, gen);

			rubberband_migs[gen].num_moved_events = 0;
			event = -1;
			new_age = 0.0;
			ntj_gen1[0] = ntj_gen1[1] = 0;
			for(i=0; i<genetree_migs[gen].num_migs; i++) {
				mig = genetree_migs[gen].living_mignodes[i];
				mig_band  = genetree_migs[gen].mignodes[mig].migration_band;
				sourcePop = genetree_migs[gen].mignodes[mig].source_pop;
				targetPop = genetree_migs[gen].mignodes[mig].target_pop;
				age       = genetree_migs[gen].mignodes[mig].age;

				if(age < taub[0] || age > taub[1])		continue;

				// we assume here that there are no in/out migrations from root population.
				if(sourcePop == pop) {
					inORout = 1;			// indicating out migration
					event = genetree_migs[gen].mignodes[mig].target_event;
					migPop = targetPop;
					new_age = taub[age>tauold] + taufactor[age>tauold]*(age-taub[age>tauold]);
					// rubberBand only counts migrations coming into pops
					ntj_gen1[age>tauold]++;
				} else if(targetPop == pop) {
					inORout = 0;			// indicating in migration
					event = genetree_migs[gen].mignodes[mig].source_event;
					migPop = sourcePop;
					new_age = taub[age>tauold] + taufactor[age>tauold]*(age-taub[age>tauold]);
					ntj_gen1[age>tauold]++;
				}


				if(event >= 0) {
					inode = genetree_migs[gen].mignodes[mig].gtree_branch;
					// check for conflicts
					if(new_age >= dataSetup.popTree->migBands[mig_band].endTime) {
						mig_conflict = 1;
					} else if(new_age <= dataSetup.popTree->migBands[mig_band].startTime) {
						mig_conflict = 1;
					} else if(inORout == 0 && new_age > age) {
						// an incoming migration event can conflict with event directly above it
						fatherNode = getNodeFather(dataState.lociData[gen], inode);
						mig1 = findFirstMig(gen,inode,genetree_migs[gen].mignodes[mig].age);
						if(mig1 >=0 &&
								genetree_migs[gen].mignodes[mig1].source_pop != pop && new_age >= genetree_migs[gen].mignodes[mig1].age) {
							mig_conflict = 1;
						} else if(fatherNode >= 0 && new_age >= getNodeAge(dataState.lociData[gen], fatherNode)) {
							mig_conflict = 1;
						}
					} else if(inORout == 1 && new_age < age) {
						// outgoing migration events can conflict with event directly below it
						mig1 = findLastMig(gen,inode,genetree_migs[gen].mignodes[mig].age);
						if(mig1 >=0 &&
								genetree_migs[gen].mignodes[mig1].target_pop != pop && new_age <= genetree_migs[gen].mignodes[mig1].age) {
							mig_conflict = 1;
						} else if(new_age <=getNodeAge(dataState.lociData[gen], inode)) {
							mig_conflict = 1;
						}
					}
					if(mig_conflict)	break;

					rubberband_migs[gen].orig_events[rubberband_migs[gen].num_moved_events] = event;
					rubberband_migs[gen].pops[rubberband_migs[gen].num_moved_events] = migPop;
					rubberband_migs[gen].new_ages[rubberband_migs[gen].num_moved_events] = new_age;
					rubberband_migs[gen].num_moved_events++;
					event = -1;
				}
			} // end for(mignode)

			if(mig_conflict)  {
				rubberband_migs[gen].num_moved_events = 0;
				break;
			}

			// create new events for affected migration bands
			for(i=0; i<num_affected_mig_bands; i++) {
				mig_band = affected_mig_bands[i];
				targetPop = dataSetup.popTree->migBands[mig_band].targetPop;
				for(event = event_chains[gen].first_event[targetPop]; event >= 0; event = event_chains[gen].events[event].next) {
					if(event_chains[gen].events[event].node_id == mig_band &&
							((start_or_end[i] && event_chains[gen].events[event].type == MIG_BAND_START) ||
									event_chains[gen].events[event].type == MIG_BAND_END  ) )     break;
				}
				if(event<0) {
					if(debug) {
						fprintf(stderr, "\nError: UpdateSampleAge: couldn't find event for migration band %d in gen %d.\n", mig_band,gen);
					} else {
						fprintf(stderr, "Fatal Error 0174.\n");
					}
					printGenealogyAndExit(gen, -1);
				}
				rubberband_migs[gen].orig_events[rubberband_migs[gen].num_moved_events] = event;
				rubberband_migs[gen].pops[rubberband_migs[gen].num_moved_events] = targetPop;
				rubberband_migs[gen].new_ages[rubberband_migs[gen].num_moved_events] = new_band_ages[i];
				rubberband_migs[gen].num_moved_events++;
			}

			// compute residual effects of rubber-band (before actual rubber-band
			genDeltaLogLikelihood[gen] = rubberBandRipple(gen, 1 /*do changes*/);

			genDeltaLogLikelihood[gen]  += rubberBand(gen, pop, taub[1], tauold, taufactor[1], 0 /*don't change chain*/, &ntj_gen1[1]);
			genDeltaLogLikelihood[gen]  += rubberBand(gen, pop, taub[0], tauold, taufactor[0], 0 /*don't change chain*/, &ntj_gen1[0]);

			genDeltaLnLd += genDeltaLogLikelihood[gen];

#ifdef CHECK_OPERATIONS
			if(	(ntj_gen[0] != ntj_gen1[0] || ntj_gen[1] != ntj_gen1[1]) ||
					(ntj_gen[0] != ntj_gen1[0]+ ntj_gen1[1] || ntj_gen[1] != 0)){
				fprintf(stderr, "Error: UpdateSampleAge has incorrect computation of modified nodes for gen %d, split %d: below %d, %d  ; above %d, %d.\n",
						gen, pop, ntj_gen[0] , ntj_gen1[0] , ntj_gen[1] , ntj_gen1[1]);
				printGenealogyAndExit(gen, -1);
			}
			ntj_gen1[0] = ntj_gen[0];
			ntj_gen1[1] = ntj_gen[1];
#endif         
			ntj[0]+=ntj_gen1[0];
			ntj[1]+=ntj_gen1[1];

			dataDeltaLnLd -= getLocusDataLikelihood(dataState.lociData[gen]);
			dataDeltaLnLd += computeLocusDataLikelihood(dataState.lociData[gen], /*reuse old conditionals*/ 1);
		}// end for(gen) - genealogy updates by rubberband

		lnacceptance += dataDeltaLnLd + genDeltaLnLd + ntj[0]*log(taufactor[0]) + ntj[1]*log(taufactor[1]);


#ifdef LOG_STEPS
		if(mig_conflict) {
			fprintf(ioSetup.debugFile, "migration conflict at gen %d, ", gen);
		} else {
			fprintf(ioSetup.debugFile, "lnacceptance = %g, ",lnacceptance);
		}
#endif	
		//No migraton conflict   Positive acceptance    accept some even if acceptance below 0 based on random probability
		if(!mig_conflict && (lnacceptance>=0 || rndu()<exp(lnacceptance))) {

#ifdef LOG_STEPS
			fprintf(ioSetup.debugFile, "accepting.\n");
#endif
			// UNUSED        didAccept = 1;
			accepted[pop]++;
			dataState.dataLogLikelihood += dataDeltaLnLd;
			dataState.logLikelihood += (dataDeltaLnLd + genDeltaLnLd) / dataSetup.numLoci;
			for(gen=0; gen<dataSetup.numLoci; gen++) {
				genLogLikelihood[gen] += genDeltaLogLikelihood[gen];
				// change gene trees, event chains, and likelihoods
				rubberBand(gen, pop, taub[1], tauold, taufactor[1], 1 /*change chain*/, &ntj_gen1[1]);
				rubberBand(gen, pop, taub[0], tauold, taufactor[0], 1 /*change chain*/, &ntj_gen[0]);

				// accept genealogy changes
				resetSaved(dataState.lociData[gen]);

				// remove original added events for migrations and migration bands
				for(i=0; i<rubberband_migs[gen].num_moved_events; i++) {
					// set pointers from mignodes to new events
					mig = event_chains[gen].events[rubberband_migs[gen].new_events[i]].node_id;
					if(event_chains[gen].events[rubberband_migs[gen].new_events[i]].type == IN_MIG) {
						genetree_migs[gen].mignodes[mig].target_event = rubberband_migs[gen].new_events[i];
						// adjust ages of mignodes for migrations out of rubberband
						// the ones coming in are adjusted in rubberBand.
						// THIS IS TO ENSURE CORRECTLY ADDRESSING MIGRATIONS BETWEEN TWO CHILDREN POPULATIONS
						// AFFECTED BY THE RUBBER BAND
						genetree_migs[gen].mignodes[mig].age = rubberband_migs[gen].new_ages[i];
					} else if(event_chains[gen].events[rubberband_migs[gen].new_events[i]].type == OUT_MIG) {
						genetree_migs[gen].mignodes[mig].source_event = rubberband_migs[gen].new_events[i];
					}
					removeEvent(gen,rubberband_migs[gen].orig_events[i]);
				}
				rubberband_migs[gen].num_moved_events = 0;
			}// end of for(gen) - implement genealogy changes

			// commit to new sample age [ no need to change mig-band times - these were changed already ]
			dataSetup.popTree->pops[pop]->sampleAge = taunew;

		}
		else {
#ifdef LOG_STEPS
			fprintf(ioSetup.debugFile, "rejecting.\n");
#endif
			if(mig_conflict) {
				//				printf("(migration conflict at gen %d)\n",gen);
				misc_stats.rubberband_mig_conflicts++;
			}
			else {
				gen = dataSetup.numLoci-1;
			}
			// start from gen before last and redo changes
			for( ; gen >= 0; --gen) {
				// redo changes in events for migrations and mig bands.
				revertToSaved(dataState.lociData[gen]);
				rubberBandRipple(gen, 0 /*redo changes*/);
			}
		}
		/**
    if (!checkAll()) {
      printf("\n  --  Aborting after UpdateSampleAge for pop %d, accepted = %d.\n",pop, didAccept);
      if(didAccept) {
         printf("Updating from time %g to time %g.\n",tauold, taunew);
      }
      if(mig_conflict) printf("Found migration conflict.\n");
      exit(-1);
   }
   // else {
   //    printf("\n  UpdateSampleAge for ancestral pop %d OK.",pop);
   // }

	for(gen=0; gen<dataSetup.numLoci; gen++) {
		if (!synchronizeEvents(gen)) {
			printf("\n  --  Aborting due to problems found when synchronizing data structures for locus #%d after UpdateSampleAge(pop=%d) %g-->%g (accepted=%d).\n\n",
				   gen+1, pop, tauold, taunew, didAccept);
			printGenealogyAndExit(gen,-1);
		}
	}
		 **/
	}// end of for(pop)
}
/** end of UpdateSampleAge **/



/***********************************************************************************
 *	UpdateLocusRate
 *	- perturbs locus-specific mutation rates
 *	- for now does not estimate heredity multipliers !!!
 *	- reference genealogy is the first one (arbitrarily)
 *	- for each locus other than referencw, changes rate and changes rate of reference
 *		accordingly (to maintain an average rate of 1)
 *	- this step affects data likelihood
 *	- this step does not affect any of the recorded statistics
 ***********************************************************************************/
int UpdateLocusRate (double finetune)	{
	int  accepted = 0, gen, genRateRef = mcmcSetup.genRateRef;
	double lnacceptance, lnLd;
	double rold, rnew, rrefold, rrefnew;
	//	double hold, hnew;
	//	double delta_lnLd, factor;
	//	int pop;

	if(finetune <= 0.0) {
		return 0;
	}

	for(gen=0; gen<dataSetup.numLoci; gen++) {
		if(gen == genRateRef)		continue;
		rrefold = getLocusMutationRate(dataState.lociData[genRateRef]);
		rold = getLocusMutationRate(dataState.lociData[gen]);
		rnew = rold + finetune*rnd2normal8();
		rnew = reflect(rnew, 0, rold+rrefold);
		setLocusMutationRate(dataState.lociData[gen], rnew);
		rrefnew = rrefold + rold - rnew;
		setLocusMutationRate(dataState.lociData[genRateRef], rrefnew);
#ifdef LOG_STEPS
		fprintf(ioSetup.debugFile, "  gen %d, proposing locus-rate shift: %g-->%g, (ref rate %g-->%g),  ",gen, rold, rnew, rrefold, rrefnew);
#endif

		// prior acceptance rate
		// should correspond to Dirichlet with all alphas set to varRatesAlpha
		// maybe this needs to be factored by the number of loci??
		// Does not matter when varRatesAlpha = 1.0 !
		// CHECK THIS !!!
		lnacceptance = (mcmcSetup.varRatesAlpha - 1) * log((rnew*rrefnew)/(rold*rrefold));

		// compute delta in log likelihood of gen and reference gen
		lnLd  = - ( getLocusDataLikelihood(dataState.lociData[gen]) + getLocusDataLikelihood(dataState.lociData[genRateRef]) );
		lnLd += computeLocusDataLikelihood(dataState.lociData[gen], /*recompute from scratch*/ 0);
		lnLd += computeLocusDataLikelihood(dataState.lociData[genRateRef], /*recompute from scratch*/ 0);

		lnacceptance += lnLd;

#ifdef LOG_STEPS
		fprintf(ioSetup.debugFile, "lnacceptance = %g, ",lnacceptance);
#endif
		if(lnacceptance>=0 || rndu()<exp(lnacceptance)) {
#ifdef LOG_STEPS
			fprintf(ioSetup.debugFile, "accepting.\n");
#endif
			accepted ++;
			dataState.dataLogLikelihood += lnLd;
			dataState.logLikelihood += lnLd/dataSetup.numLoci;
			resetSaved(dataState.lociData[gen]);
			resetSaved(dataState.lociData[genRateRef]);
			dataState.rateVar += ( rnew*rnew + rrefnew*rrefnew - rold*rold - rrefold*rrefold ) / dataSetup.numLoci;
		}
		else {
#ifdef LOG_STEPS
			fprintf(ioSetup.debugFile, "rejecting.\n");
#endif
			setLocusMutationRate(dataState.lociData[gen], rold);
			setLocusMutationRate(dataState.lociData[genRateRef], rrefold);
			revertToSaved(dataState.lociData[gen]);
			revertToSaved(dataState.lociData[genRateRef]);
		}
	}// end of for(gen)
	return(accepted);
}
/** end of UpdateLocusRate **/



/***********************************************************************************
 *	mixing
 *	- scales all population parameters by  a uniform constant
 *	- migration rates are scaled in the other direction from taus/thetas
 *	- this step doesn't really change the genealogy statistics, but just the data likelihood
 ***********************************************************************************/
/* MARK: NOTE THAT THIS STEP CHANGES THE AGE OF ANCIENT POPULATION SAMPLE
 */
int mixing (double finetune)	{
	double xold, xnew, c, lnc, lnacceptance, dataDeltaLnLd, genDeltaLnLd;

	int gen, i, mig, mig_band, pop, num_events;

	unsigned short rejectIssue = 0;	// a flag which indicates if found any issue that results in a-priori rejection


	if(finetune <= 0.0) {
		return 0;
	}

	// sample a multiplicative factor
	lnc = finetune*rnd2normal8();
	c = exp(lnc);
#ifdef LOG_STEPS
	fprintf(ioSetup.debugFile, "  mixing factor %g (log=%g), ",c, lnc);
#endif

	// compute number of coalescent and migration events whose age is scaled
	num_events=0;
	for(pop=0; pop<dataSetup.popTree->numPops; pop++) {
		num_events += genetree_stats_total.num_coals[pop];
	}
	for(mig_band=0; mig_band<dataSetup.popTree->numMigBands; mig_band++) {
		num_events += genetree_stats_total.num_migs[mig_band];
	}

	// proposal ratio - note that migration rates are scaled in the other direction
	lnacceptance = lnc * (2*dataSetup.popTree->numPops - dataSetup.popTree->numCurPops - dataSetup.popTree->numMigBands + num_events);
	//	lnacceptance = lnc;

	dataDeltaLnLd = 0.0;
	genDeltaLnLd = 0.0;
	for(pop=0; pop<dataSetup.popTree->numPops; pop++) {
		xold = dataSetup.popTree->pops[pop]->theta;
		dataSetup.popTree->pops[pop]->theta = xnew = xold*c;
		lnacceptance += lnc * (dataSetup.popTree->pops[pop]->thetaPrior.alpha - 1) - (xnew-xold) *dataSetup.popTree->pops[pop]->thetaPrior.beta;
		// change in genetree likelihoods is not in the stats
		// because times and rates are scaled together.
		// the difference is only in coalescence/migration densities.
		// this actually cancels out with proposal ratio.
		genDeltaLnLd -= lnc * genetree_stats_total.num_coals[pop];
		if(pop<dataSetup.popTree->numCurPops && dataSetup.popTree->pops[pop]->sampleAge > 0.0) {
			dataSetup.popTree->pops[pop]->sampleAge *= c;
		}
	}
	for(pop=dataSetup.popTree->numCurPops; pop<dataSetup.popTree->numPops; pop++) {
		// consider also current populations with ancient samples
		xold = dataSetup.popTree->pops[pop]->age;
		dataSetup.popTree->pops[pop]->age = xnew = xold*c;
		lnacceptance += lnc * (dataSetup.popTree->pops[pop]->agePrior.alpha - 1) - (xnew-xold) *dataSetup.popTree->pops[pop]->agePrior.beta;
	}
	for(mig_band=0; mig_band<dataSetup.popTree->numMigBands; mig_band++) {
		// migration rates are scaled inversely
		xold = dataSetup.popTree->migBands[mig_band].migRate;
		dataSetup.popTree->migBands[mig_band].migRate = xnew = xold/c;
		// see if migration rate got out of bounds
		//		if(xnew < 0.00001 || xnew > dataSetup.popTree->migBands[mig_band].upperBound) {
		//		if(xold > 0 && xnew < 0.0000001 || xnew > MAX_MIG_RATE) {
#ifdef LOG_STEPS
		//			fprintf(ioSetup.debugFile, "mig rate out of bound for mig-band %d, ",mig_band);
#endif
		//			rejectIssue = 1;
		//		}
		// GAMMA PRIOR
		lnacceptance += -lnc * (dataSetup.popTree->migBands[mig_band].migRatePrior.alpha - 1) - (xnew-xold) *dataSetup.popTree->migBands[mig_band].migRatePrior.beta;
		dataSetup.popTree->migBands[mig_band].startTime *= c;
		dataSetup.popTree->migBands[mig_band].endTime *= c;
		// change in genetree likelihoods is not in the stats
		// because times and rates are scaled together.
		// the difference is only in coalescence/migration densities.
		// this actually cancels out with proposal ratio.
		genDeltaLnLd -= lnc * genetree_stats_total.num_migs[mig_band];
	}

	if(!rejectIssue) {
		// adjust all gen genealogies
		for(gen=0; gen<dataSetup.numLoci; gen++) {
			// scale age of nodes and compute delta likelihood
			dataDeltaLnLd += scaleAllNodeAges(dataState.lociData[gen], c);
		}

		lnacceptance += (dataDeltaLnLd + genDeltaLnLd);

#ifdef LOG_STEPS
		fprintf(ioSetup.debugFile, "lnacceptance = %g, ",lnacceptance);
#endif
		if(lnacceptance>=0 || rndu()<exp(lnacceptance)) { /* accept */
#ifdef LOG_STEPS
			fprintf(ioSetup.debugFile, "accepting.\n");
#endif
			for(gen=0; gen<dataSetup.numLoci; gen++) {
				resetSaved(dataState.lociData[gen]);
				for(i=0; i<genetree_migs[gen].num_migs; i++) {
					mig = genetree_migs[gen].living_mignodes[i];
					genetree_migs[gen].mignodes[mig].age *= c;
				}

				genLogLikelihood[gen] -=  lnc * (dataSetup.numSamples - 1 + genetree_migs[gen].num_migs);

				// update all statistics by the constant
				for(pop=0; pop<dataSetup.popTree->numPops; pop++) {
					genetree_stats[gen].coal_stats[pop] *= c;
				}
				for(mig_band=0;mig_band<dataSetup.popTree->numMigBands; mig_band++) {
					genetree_stats[gen].mig_stats[mig_band] *= c;
				}

				// update elapsed times of all valid events
				for(i=0; i<event_chains[gen].total_events; i++) {
					if(event_chains[gen].events[i].elapsed_time > 0) event_chains[gen].events[i].elapsed_time *= c;
				}
			}// end of for(gen)

			// update total statistics by the constant
			for(pop=0; pop<dataSetup.popTree->numPops; pop++) {
				genetree_stats_total.coal_stats[pop] *= c;
			}
			for(mig_band=0;mig_band<dataSetup.popTree->numMigBands; mig_band++) {
				genetree_stats_total.mig_stats[mig_band] *= c;
			}

			dataState.dataLogLikelihood += dataDeltaLnLd;
			dataState.logLikelihood += (dataDeltaLnLd + genDeltaLnLd)/dataSetup.numLoci;
			adjustRootEvents();

			return 1;
		}
	}

	// reject
#ifdef LOG_STEPS
	fprintf(ioSetup.debugFile, "rejecting.\n");
#endif
	if(!rejectIssue) {
		for(gen=0; gen<dataSetup.numLoci; gen++) {
			revertToSaved(dataState.lociData[gen]);
		}
	}

	// revert to old parameters and genealogies
	for(pop=0; pop<dataSetup.popTree->numPops; pop++) {
		dataSetup.popTree->pops[pop]->theta /= c;
	}
	//  for(pop=dataSetup.popTree->numCurPops; pop<dataSetup.popTree->numPops; pop++) {
	for(pop=0; pop<dataSetup.popTree->numPops; pop++) {
		dataSetup.popTree->pops[pop]->age /= c;
		if(pop<dataSetup.popTree->numCurPops && dataSetup.popTree->pops[pop]->sampleAge > 0.0) {
			dataSetup.popTree->pops[pop]->sampleAge /= c;
		}
	}
	for(mig_band=0; mig_band<dataSetup.popTree->numMigBands; mig_band++) {
		// migration rates are scaled inversely
		dataSetup.popTree->migBands[mig_band].migRate   *= c;
		dataSetup.popTree->migBands[mig_band].startTime /= c;
		dataSetup.popTree->migBands[mig_band].endTime   /= c;
	}

	return 0;
}
/** end of mixing **/



/***************************************************************************************************************/
/******                                        END OF FILE                                                ******/
/***************************************************************************************************************/
