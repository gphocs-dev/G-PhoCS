/** 
	\file MCMCcontrol.c
    Read and processes and holds control information for MCMC

	Contains the relevant data structures and procedure implementations for
	reading and processing a control file.
	
*/

#include "MCMCcontrol.h"


/***************************************************************************************************************/
/******                                  GLOBAL INTERNAL DATA STRUCTURES                                  ******/
/***************************************************************************************************************/



/*********
 * globalSetup
 *********/
struct GLOBAL_SETUP {
	double alpha;
	double beta;
	double printFactor;
	double migAlpha;
	double migBeta;
	double migFactor;
	double finetuneTaus;
}  globalSetup;




/***************************************************************************************************************/
/******                                INTERNAL FUNCTION DECLARATIONS                                     ******/
/***************************************************************************************************************/



char* getNextToken(FILE* file, char* space);
int expectNextToken(FILE* file, const char* expectedToken, char* tokenSpace);
int countTokens(FILE* file, char* countToken, char* endToken, char* tokenSpace);

int readGeneralInfo(FILE* fctl);
int readCurrentPops(FILE* fctl);
int readAncestralPops(FILE* fctl);
int readMigrationBands(FILE* fctl);

int readSampleLine(char* sampleLine, int pop);
int parseSampleNames();



/***************************************************************************************************************/
/******                              EXTERNAL FUNCTION IMPLEMENTATION                                     ******/
/***************************************************************************************************************/



/***********************************************************************************
 *	initGeneralInfo
 * 	- initializes control and I/O settings to default settings
 *	- returns 0
 ***********************************************************************************/
int initGeneralInfo() {
	
	globalSetup.printFactor = 1.0;
	globalSetup.alpha = -1.0;
	globalSetup.beta	= -1.0;
	globalSetup.migFactor = 1.0;
	globalSetup.migAlpha = -1.0;
	globalSetup.migBeta	= -1.0;
	globalSetup.finetuneTaus = -1.0;
	
	strcpy(ioSetup.seqFileName,"NONE");
	strcpy(ioSetup.rateFileName,"NONE");
	strcpy(ioSetup.nodeStatsFileName, "NONE");
	strcpy(ioSetup.combStatsFileName, "NONE");
	strcpy(ioSetup.traceFileName, "mcmc-trace.out");

	ioSetup.samplesPerLog 	= 100;
	ioSetup.logsPerLine 	= 100;

	mcmcSetup.randomSeed = -1;
	mcmcSetup.useData  = 0;
	mcmcSetup.startMig = 0;
	mcmcSetup.mutRateMode = 0;
	mcmcSetup.allowAdmixture = 0;
	mcmcSetup.numSamples = 10000;
	mcmcSetup.doMixing   = 1;
	mcmcSetup.burnin     = 0;
	mcmcSetup.sampleSkip = 0;
	mcmcSetup.findFinetunes = 0;;
	mcmcSetup.findFinetunesSamplesPerStep = 100;
	mcmcSetup.findFinetunesNumSteps = 100;
	mcmcSetup.genetreeSamples = 1;
	mcmcSetup.finetunes.coalTime = -1.0;
	mcmcSetup.finetunes.migTime = -1.0;
	mcmcSetup.finetunes.theta = -1.0;
	mcmcSetup.finetunes.migRate = -1.0;
	mcmcSetup.finetunes.taus = NULL;
	mcmcSetup.finetunes.locusRate =  -1.0;
	mcmcSetup.finetunes.mixing = -1.0;

	mcmcSetup.finetunes.admix = -1.0;

	dataSetup.numLoci = -1;
	dataSetup.numPopPartitions = 0;
	return 0;
} 
/** end of initGeneralInfo **/



/***********************************************************************************
 *	readControlFile
 * 	- reads control file and initializes control and I/O settings
 *	- returns 0, if all OK, and -1 otherwise.
 ***********************************************************************************/
int readControlFile(char* controlFileName) {
	
	int numErrors = 0;
	FILE *fctl=(FILE*)fopen(controlFileName,"r");

	if(fctl == NULL) {
		fprintf(stderr, "Error: Could not open control file '%s'.\n", controlFileName);
		return -1;
	}

	numErrors += readGeneralInfo(fctl);
	
	if(numErrors > 0) {
		fprintf(stderr, "Found %d errors when parsing GENERAL-INFO in control file %s.\n", numErrors, controlFileName);
		fclose(fctl);
		return -1;
	}
	

	numErrors += readCurrentPops(fctl);
	
	numErrors += parseSampleNames();
	
	if(numErrors > 0) {
		fprintf(stderr, "Found %d errors when parsing CURRENT-POPS in control file %s.\n", numErrors, controlFileName);
		fclose(fctl);
		return -1;
	}
	numErrors += readAncestralPops(fctl);
	
	if(numErrors > 0) {
		fprintf(stderr, "Found %d errors when parsing ANCESTRAL-POPS in control file %s.\n", numErrors, controlFileName);
		fclose(fctl);
		return -1;
	}

	numErrors += readMigrationBands(fctl);
	
	if(numErrors > 0) {
		fprintf(stderr, "Found %d errors when parsing MIG-BANDS in control file %s.\n", numErrors, controlFileName);
		fclose(fctl);
		return -1;
	}

	
	fclose(fctl);
	return 0;
}
/** end of readControlFile **/



/***********************************************************************************
 *	readSecondaryControlFile
 * 	- reads secondary control file with only general info and mig bands
 *	- returns 0, if all OK, and -1 otherwise.
 ***********************************************************************************/
int readSecondaryControlFile(char* controlFileName) {
	
	int numErrors = 0;
	FILE *fctl=(FILE*)fopen(controlFileName,"r");

	if(fctl == NULL) {
		fprintf(stderr, "Error: Could not open secondary control file '%s'.\n", controlFileName);
		return -1;
	}
	

	// read general info
	numErrors += readGeneralInfo(fctl);
	
	if(numErrors > 0) {
		fprintf(stderr, "Found %d errors when parsing GENERAL-INFO in secondary control file %s.\n", numErrors, controlFileName);
		fclose(fctl);
		return -1;
	}
	

	numErrors += readMigrationBands(fctl);
	
	if(numErrors > 0) {
		fprintf(stderr, "Found %d errors when parsing MIG-BANDS in secondary control file %s.\n", numErrors, controlFileName);
		fclose(fctl);
		return -1;
	}

	fclose(fctl);
	return 0;
}
/** end of readSecondaryControlFile **/



/***********************************************************************************
 *	checkSettings
 * 	- tests validity and completeness of settings collected in control file(s)
 *	- returns number of errors found
 ***********************************************************************************/
int checkSettings() {
	int numErrors = 0;
	int pop, migBand;

	if(0 != strcmp("NONE",ioSetup.seqFileName)) {
		mcmcSetup.useData = 1;
	} else {    //Require user to use a sequence file.
		fprintf(stderr, "Error: Sequences file ('seq-file') is not defined in the control file.\n");
		numErrors++;  
	}

	if(0 == strcmp("NONE",ioSetup.nodeStatsFileName) && dataSetup.numPopPartitions > 0) {
		fprintf(stderr, "Error: number of population partitions is %d, but no stats file name was specified.\n", dataSetup.numPopPartitions);
		numErrors++;  
	}
	
	if(!mcmcSetup.findFinetunes) {	
		// check that all finetunes were specified
		// do not check for taus (those might be specified individually
		// do not check for admixture - might need only if there are admixed individuals
		if(mcmcSetup.finetunes.coalTime < 0.0) {
			fprintf(stderr, "Error: positive finetune for coal-time should be specified.\n");
			numErrors++;
		}
		if(mcmcSetup.finetunes.migTime < 0.0) {
			fprintf(stderr, "Error: positive finetune for mig-time should be specified.\n");
			numErrors++;
		}
		if(mcmcSetup.finetunes.theta < 0.0) {
			fprintf(stderr, "Error: positive finetune for theta should be specified.\n");
			numErrors++;
		}
		if(mcmcSetup.finetunes.migRate < 0.0) {
			fprintf(stderr, "Error: positive finetune for mig-rate should be specified.\n");
			numErrors++;
		}
		if(mcmcSetup.mutRateMode == 1 && mcmcSetup.finetunes.locusRate < 0.0) {
			fprintf(stderr, "Error: positive finetune for locus-rate should be specified.\n");
			numErrors++;  
		}
		if(mcmcSetup.finetunes.mixing < 0.0) {
			fprintf(stderr, "Error: positive finetune for mixing should be specified.\n");
			 numErrors++;
		}
	}
	
	if(ioSetup.samplesPerLog <= 0) {
		fprintf(stderr, "Warning: samples-per-log must be 1 or greater, adjusting to 100.\n");
		ioSetup.samplesPerLog = 100;    
	}	
	if(ioSetup.logsPerLine <= 0) {
		fprintf(stderr, "Warning: logs-per-line must be 1 or greater, adjusting to 100.\n");
		ioSetup.logsPerLine = 100;
	}
	
	// set start point of tau for sampling for prior mean, if not pre-set
	for(pop=dataSetup.popTree->numCurPops; pop<dataSetup.popTree->numPops; pop++) {
		if(dataSetup.popTree->pops[pop]->agePrior.sampleStart <= 0) {
			dataSetup.popTree->pops[pop]->agePrior.sampleStart =  dataSetup.popTree->pops[pop]->agePrior.alpha / dataSetup.popTree->pops[pop]->agePrior.beta;
		}
	}
	
	// check population specific settings
	for(pop=0; pop<dataSetup.popTree->numPops; pop++) {
		if(dataSetup.popTree->pops[pop]->thetaPrior.alpha < 0) {
			fprintf(stderr, "Error: gamma prior alpha parameter not set for theta of pop %s (%d).\n", dataSetup.popTree->pops[pop]->name,pop+1);
			numErrors++;
		}
		if(dataSetup.popTree->pops[pop]->thetaPrior.beta < 0) {
			fprintf(stderr, "Error: gamma prior beta argument not set for theta of pop %s (%d).\n", dataSetup.popTree->pops[pop]->name,pop+1);
			numErrors++;
		}
		// set start point of theta for sampling for prior mean
		dataSetup.popTree->pops[pop]->thetaPrior.sampleStart =  dataSetup.popTree->pops[pop]->thetaPrior.alpha / dataSetup.popTree->pops[pop]->thetaPrior.beta;

		// for ancestral pops
		if(pop >= dataSetup.popTree->numCurPops) {
			if(dataSetup.popTree->pops[pop]->agePrior.alpha < 0) {
				fprintf(stderr, "Error: gamma prior alpha parameter not set for tau of ancestral pop %s (%d).\n", dataSetup.popTree->pops[pop]->name,pop+1);
				numErrors++;
			}
			if(dataSetup.popTree->pops[pop]->agePrior.beta < 0) {
				fprintf(stderr, "Error: gamma prior beta argument not set for tau of ancestral pop %s (%d).\n", dataSetup.popTree->pops[pop]->name,pop+1);
				numErrors++;
			}
			if( !mcmcSetup.findFinetunes && mcmcSetup.finetunes.taus[pop] < 0) {
				fprintf(stderr, "Error: finetune not set for update of tau of ancestral pop %s (%d).\n", dataSetup.popTree->pops[pop]->name,pop+1);
				numErrors++;
			}

			
			// check to see that population tau priors  and start times are lower than parent's
			if(dataSetup.popTree->rootPop != pop && 
				( dataSetup.popTree->pops[pop]->father->agePrior.alpha/dataSetup.popTree->pops[pop]->father->agePrior.beta < 
					dataSetup.popTree->pops[pop]->agePrior.alpha/dataSetup.popTree->pops[pop]->agePrior.beta) ) {
						fprintf(stderr, "\nError:Conflicting priors for ancestral population ages found for pop %s, and parent pop %s.\n",
              					dataSetup.popTree->pops[pop]->name, dataSetup.popTree->pops[pop]->father->name);
						numErrors++;
			}
			if(dataSetup.popTree->rootPop != pop && 
				(dataSetup.popTree->pops[pop]->father->agePrior.sampleStart < dataSetup.popTree->pops[pop]->agePrior.sampleStart) ) {
						fprintf(stderr, "\nError:Conflicting initalization settings for ancestral population ages found for pop %s, and parent pop %s.\n",
              					dataSetup.popTree->pops[pop]->name, dataSetup.popTree->pops[pop]->father->name);
						numErrors++;
			}
		}// end of if(pop>=numCurPops)
		else {
			if( ( dataSetup.popTree->pops[pop]->father->agePrior.alpha/dataSetup.popTree->pops[pop]->father->agePrior.beta < 
					dataSetup.popTree->pops[pop]->sampleAge ) ) {
						fprintf(stderr, "\nError:Conflicting prior for ancestral population age for parent pop %s and sample age for pop %s.\n",
              					dataSetup.popTree->pops[pop]->father->name , dataSetup.popTree->pops[pop]->name );
						numErrors++;
			}
			if( (dataSetup.popTree->pops[pop]->father->agePrior.sampleStart < dataSetup.popTree->pops[pop]->sampleAge) ) {
						fprintf(stderr, "\nError:Conflicting initialization for ancestral population age for parent pop %s and sample age for pop %s (%g,%g).\n",
              					dataSetup.popTree->pops[pop]->father->name , dataSetup.popTree->pops[pop]->name,
								dataSetup.popTree->pops[pop]->father->agePrior.sampleStart,dataSetup.popTree->pops[pop]->sampleAge
   							);
						numErrors++;
			}
		}
	}// end of for(pop)

	for(migBand=0; migBand<dataSetup.popTree->numMigBands; migBand++) {
		// check that priors are set
		if(dataSetup.popTree->migBands[migBand].migRatePrior.alpha < 0) {
			fprintf(stderr, "Error: gamma prior alpha argument not set for mig-rate of mig-band (#%d).\n", migBand+1);
			numErrors++;
		}
		if(dataSetup.popTree->migBands[migBand].migRatePrior.beta < 0) {
			fprintf(stderr, "Error: gamma prior beta argument not set for mig-rate of mig-band (#%d).\n", migBand+1);
			numErrors++;
		}
		
		// set start point of migration rate for sampling for prior mean
		dataSetup.popTree->migBands[migBand].migRatePrior.sampleStart =  
				dataSetup.popTree->migBands[migBand].migRatePrior.alpha / dataSetup.popTree->migBands[migBand].migRatePrior.beta;
	}// end of for(migBand)
  return numErrors;
}
/** end of checkSettings **/



/******************************************************************************
 *	printPriorSettings
 * 	- prints the prior settings onto standard output
 *	- returns 0
 *****************************************************************************/
int printPriorSettings()
{
  int pop, migBand;
  char string[1000];
  double alpha,beta;
	
  // print out tree and population parameters
  printPopulationTree(dataSetup.popTree, stdout, 0);
  printf("\n");
  printf("---------------------------------------------------------------\n");
  

  printf("\nGamma prior: mean +- SE for theta's ,tau's and m's\n");
  for(pop=0; pop<dataSetup.popTree->numPops; pop++)
  {
    sprintf(string, "theta_%s: ", dataSetup.popTree->pops[pop]->name);
    alpha = dataSetup.popTree->pops[pop]->thetaPrior.alpha;
    beta  = dataSetup.popTree->pops[pop]->thetaPrior.beta;
    printf("%-15s %9.5f +- %9.5f\n", string, alpha/beta, sqrt(alpha)/beta);
  }
  printf("---------------------------------------------------------------\n");

  for(pop=dataSetup.popTree->numCurPops; pop<dataSetup.popTree->numPops; pop++)
  {
    sprintf(string, "tau_%s: ", dataSetup.popTree->pops[pop]->name);
    alpha = dataSetup.popTree->pops[pop]->agePrior.alpha;
    beta  = dataSetup.popTree->pops[pop]->agePrior.beta;
    printf("%-15s %9.5f +- %9.5f\n", string, alpha/beta, sqrt(alpha)/beta);
  }
  printf("---------------------------------------------------------------\n");
	
  if(dataSetup.popTree->numMigBands > 0)
  {
    for(migBand=0; migBand<dataSetup.popTree->numMigBands; migBand++)
    {
      alpha = dataSetup.popTree->migBands[migBand].migRatePrior.alpha;
      beta = dataSetup.popTree->migBands[migBand].migRatePrior.beta;
      printf("m_%s->%s:  %9.5f +- %9.5f\n",
             dataSetup.popTree->pops[ \
               dataSetup.popTree->migBands[migBand].sourcePop ]->name,
             dataSetup.popTree->pops[ \
               dataSetup.popTree->migBands[migBand].targetPop ]->name,
             alpha/beta, sqrt(alpha)/beta);
    }
    printf("---------------------------------------------------------------\n");
  }
  return 0;
}
/** end of printPriorSettings **/



/***********************************************************************************
 *	finalizeNumParameters
 * 	- determines number of parameters in the model, and finalizes printFactor array
 *	- returns 0
 ***********************************************************************************/
/* MARK: ADD PARAMETERS FOR AGE OF ANCIENT POPULATIONS
         DONE. NEED TO CHECK !!
*/
int finalizeNumParameters() {
	int numAncientPops = 0;
	int param;
	
	for(param=0; param<dataSetup.popTree->numCurPops; param++) {
		if(dataSetup.popTree->pops[param]->updateSampleAge || dataSetup.popTree->pops[param]->sampleAge > 0.0)  numAncientPops++;
	}

	mcmcSetup.numParameters = 
		2*dataSetup.popTree->numPops - dataSetup.popTree->numCurPops + 
		dataSetup.popTree->numMigBands + 
		numAncientPops +
		admixed_samples.number + 
		(mcmcSetup.mutRateMode == 1);
		
	
	
	
	mcmcSetup.printFactors = (double*)realloc(mcmcSetup.printFactors,mcmcSetup.numParameters*sizeof(double));
	if(mcmcSetup.printFactors == NULL) {
		fprintf( stderr, "Error: oom reallocating mcmcSetup.printFactors when finalizing parameters.\n");
		exit(-1);
	}
	
	for(param=2*dataSetup.popTree->numPops - dataSetup.popTree->numCurPops + dataSetup.popTree->numMigBands; 
		param<2*dataSetup.popTree->numPops - dataSetup.popTree->numCurPops + dataSetup.popTree->numMigBands+numAncientPops;
		param++) {
			mcmcSetup.printFactors[param] = globalSetup.printFactor;
	}
	
	//set print factors to 1 for all non-scaled parameters
	for(param=2*dataSetup.popTree->numPops - dataSetup.popTree->numCurPops + numAncientPops + dataSetup.popTree->numMigBands; param<mcmcSetup.numParameters; param++) {
		mcmcSetup.printFactors[param] = 1.0;
	}
	return 0;
}
/** end of finalizeNumParameters **/



/***************************************************************************************************************/
/******                              INTERNAL FUNCTION IMPLEMENTATION                                     ******/
/***************************************************************************************************************/



/***********************************************************************************
 *	getNextToken
 * 	- reads next string in file.
 *	- ignores lines starting with '#'
 *	- string is put in 'space', which should be long enough
 *	- a pointer to space is returned
 ***********************************************************************************/
char* getNextToken(FILE* file, char* space) {
  int res;
  res =	fscanf(file,"%s",space);
  if(0 > res)
    return space;
  while(!feof(file) && space[0] == '#')   {
    flushLine(file);
    res = fscanf(file,"%s",space);
    if(0 > res)
      break;
  }
  if(feof(file))
    strncpy(space, "EOF", 3);
	
  return space;
}
/** end of getNextToken **/



/*****************************************************************************
 *	expectNextToken
 * 	- reads file until reach a pre-specified token
 *	- uses given tokenSpace to read tokens
 *	- outputs warnings if reads other tokens before reaching expected token
 *  - returns -1 if reached EOF in process (0 otherwise)
 *	- a pointer to space is returned
 ****************************************************************************/
int expectNextToken(FILE* file, const char* expectedToken, char* tokenSpace) 
{
  int numErrors = 0;
  if( NULL == expectedToken )
    return 1;
  if( NULL == tokenSpace )
    return 1;

  // read general info
  while( 0 != strcmp( expectedToken,
                      getNextToken( file,
                                    tokenSpace ) ) ) 
  {
    if( feof( file ) ) 
      break;
    
    ++numErrors;
    fprintf( stderr, 
             "Error: unexpected token '%s' before '%s'. "
             "Will Ignore this.\n",tokenSpace, expectedToken );
  }
  return numErrors;
}

/** end of expectNextToken **/



/***********************************************************************************
 *	countTokens
 * 	- counts number of occurrences of given token in file until end token is observed
 *	- if end token is NULL, reads through end of file
 *	- rewinds file pointer to starting point
 *	- receives pointer for space
 *	- returns count
 ***********************************************************************************/
int countTokens(FILE* file, char* countToken, char* endToken, char* tokenSpace) {
	int count = 0;
	fpos_t filePos;	

	//record current position of file
	fgetpos(file, &filePos);
	while(!feof(file)) {
		getNextToken(file, tokenSpace);
		if(endToken == NULL || 0 == strcmp(endToken,tokenSpace)) {
			break;
		}
		if(0 == strcmp(countToken,tokenSpace)) {
			count++;
		}
	} // end of while
	
	//rewind file to starting point
	fsetpos(file, &filePos); 

	return count;
}
/** end of countTokens **/



/***********************************************************************************
 *	readGeneralInfo
 * 	- reads GENERAL-INFO part of control file and initializes control and I/O settings
 *	- returns number of errors found in parsing
 ***********************************************************************************/
int readGeneralInfo(FILE* fctl) {
	
	int numErrors = 0;
	char token[STRING_LENGTH];
	char* token2;
	char line[STRING_LENGTH];
	
	numErrors += expectNextToken(fctl,"GENERAL-INFO-START",token);
	if(numErrors > 0) {
		return numErrors;
	}


	while (1) {
		// read token
		getNextToken(fctl,token);
	    if (0 == strcmp("GENERAL-INFO-END", token))
			break;  
   
		// read next token - value of some argument
		if(NULL == fgets(line,STRING_LENGTH-1,fctl)) {
			fprintf(stderr,"Error: unexpected end of file or other error inside GENERAL-INFO module.\n");
			return numErrors+1;
		}
		token2 = strtokCS(line,parseFileDelims);
		if(token2 == NULL) {
			fprintf(stderr,"Error: unable to read value for %s.\n",token);
			numErrors++;
			continue;
		}

		if(0 == strcmp("seq-file",token)) {
			strncpy(ioSetup.seqFileName, token2,  NAME_LENGTH-1);
		} else if(0 == strcmp("trace-file",token)) {
			strncpy(ioSetup.traceFileName, token2, NAME_LENGTH-1);
		} else if(0 == strcmp("coal-stats-file",token)) {
			strncpy(ioSetup.nodeStatsFileName, token2, NAME_LENGTH-1);
		} else if(0 == strcmp("comb-stats-file",token)) {
			strncpy(ioSetup.combStatsFileName, token2, NAME_LENGTH-1);
		} else if(0 == strcmp("num-pop-partitions",token)) {
			if (sscanf(token2, "%d", &dataSetup.numPopPartitions) != 1 || dataSetup.numPopPartitions <= 0) {
				fprintf(stderr,"Error: value for num-pop-partitions should be positive integer, got %s.\n", token2);
				numErrors++;
			}
		} else if(0 == strcmp("num-loci",token)) {
			if (sscanf(token2, "%d", &dataSetup.numLoci) != 1 || dataSetup.numLoci <= 0) {
				fprintf(stderr,"Error: value for num-loci should be positive integer, got %s.\n", token2);
				numErrors++;
			}
		} else if(0 == strcmp("random-seed",token)) {
			if (sscanf(token2, "%d", &mcmcSetup.randomSeed) != 1) {
				fprintf(stderr,"Error: value for random-seed should be integer, got %s.\n", token2);
				numErrors++;
			}
		} else if (0 == strcmp("burn-in",token)) {
			if (sscanf(token2, "%d", &mcmcSetup.burnin) != 1 || mcmcSetup.burnin < 0) {
				fprintf(stderr,"Error: value for burnin should be non-negative integer, got %s.\n", token2);
				numErrors++;
			}
		} else if(0 == strcmp("mcmc-iterations",token)) {
			if (sscanf(token2, "%d", &mcmcSetup.numSamples) != 1 || mcmcSetup.numSamples <= 0) {
				fprintf(stderr,"Error: value for mcmc-iterations should be positive integer, got %s.\n", token2);
				numErrors++;
			}
		} else if(0 == strcmp("mcmc-sample-skip",token)) {
			if (sscanf(token2, "%d", &mcmcSetup.sampleSkip) != 1 || mcmcSetup.sampleSkip < 0) {
				fprintf(stderr,"Error: value for mcmc-sample-skip should be non-negative integer, got %s.\n", token2);
				numErrors++;
			}
		} else if(0 == strcmp("start-mig",token)) {
			if (sscanf(token2, "%d", &mcmcSetup.startMig) != 1 || mcmcSetup.startMig < 0) {
				fprintf(stderr,"Error: value for start-mig should be non-negative integer, got %s.\n", token2);
				numErrors++;
			}
		} else if(0 == strcmp("no-mixing",token)) {
			mcmcSetup.doMixing=0;
		} else if(0 == strcmp("iterations-per-log",token)) {
			if (sscanf(token2, "%d", &ioSetup.samplesPerLog) != 1) {
				fprintf(stderr,"Error: value for iterations-per-log should be integer, got %s.\n", token2);
				numErrors++;
			}
		} else if(0 == strcmp("logs-per-line",token)) {
			if (sscanf(token2, "%d", &ioSetup.logsPerLine) != 1) {
				fprintf(stderr,"Error: value for logs-per-line should be integer, got %s.\n", token2);
				numErrors++;
			}
		} else if(0 == strcmp("tau-theta-print",token)) {
			if (sscanf(token2, "%lf", &globalSetup.printFactor) != 1) {
				fprintf(stderr,"Error: value for tau-theta-print should be floating point number, got %s.\n", token2);
				numErrors++;
			}
		} else if(0 == strcmp("tau-theta-alpha",token)) {
			if (sscanf(token2, "%lf", &globalSetup.alpha) != 1) {
				fprintf(stderr,"Error: value for tau-theta-alpha should be floating point number, got %s.\n", token2);
				numErrors++;
			}
		} else if(0 == strcmp("tau-theta-beta",token)) {
			if (sscanf(token2, "%lf", &globalSetup.beta) != 1) {
				fprintf(stderr,"Error: value for tau-theta-beta should be floating point number, got %s.\n", token2);
				numErrors++;
			}
		} else if(0 == strcmp("mig-rate-print",token)) {
			if (sscanf(token2, "%lf", &globalSetup.migFactor) != 1) {
				fprintf(stderr,"Error: value for mig-rate-print should be floating point number, got %s.\n", token2);
				numErrors++;
			}
		} else if(0 == strcmp("mig-rate-alpha",token)) {
			if (sscanf(token2, "%lf", &globalSetup.migAlpha) != 1) {
				fprintf(stderr,"Error: value for mig-rate-alpha should be floating point number, got %s.\n", token2);
				numErrors++;
			}
		} else if(0 == strcmp("mig-rate-beta",token)) {
			if (sscanf(token2, "%lf", &globalSetup.migBeta) != 1) {
				fprintf(stderr,"Error: value for mig-rate-beta should be floating point number, got %s.\n", token2);
				numErrors++;
			}
/*	    } else if(0 == strcmp("admixture",token)) {
			if(0 == strcmp("TRUE",token2)) {
				mcmcSetup.allowAdmixture = 1;
			} else if(0 == strcmp("FALSE",token2)) {
				mcmcSetup.allowAdmixture = 0;
			} else {
				fprintf(stderr, "Error: value of admixture should be TRUE or FALSE, got %s.\n",token2);
				numErrors++;
			}		
*/	    } else if(0 == strcmp("locus-mut-rate",token)) {
			if(0 == strcmp("CONST",token2)) {
				mcmcSetup.mutRateMode = 0;
			} else if(0 == strcmp("FIXED",token2)) {
				token2 = strtokCS(NULL, parseFileDelims);
				if(token2 == NULL) {
					fprintf(stderr,"Error: unable to read filename for fixed locus mutation rates.\n");
					numErrors++;
					continue;
				}
				strncpy(ioSetup.rateFileName, token2, NAME_LENGTH-1);
				mcmcSetup.mutRateMode = 2;
			} else if(0 == strcmp("VAR",token2)) {
				token2 = strtokCS(NULL, parseFileDelims);
				if(token2 == NULL ||(sscanf(token2, "%lf", &mcmcSetup.varRatesAlpha) != 1)) {
					fprintf(stderr, "Error: unable to read floating point alpha parameter for Dirichlet prior of mutation rate variation in locus-mut-rate.\n");
					numErrors++;
				}
				mcmcSetup.mutRateMode = 1;
			} else {
				fprintf(stderr, "Error: value of const-rate should be CONST, FIXED, or VAR, got %s.\n",token2);
				numErrors++;
			}		
		} else if(0 == strcmp("finetune-coal-time",token)) {
			if (sscanf(token2, "%lf", &mcmcSetup.finetunes.coalTime) != 1) {
				fprintf(stderr,"Error: value for finetune-coal-time should be floating point number, got %s.\n", token2);
				numErrors++;
			}
		} else if(0 == strcmp("finetune-mig-time",token)) {
			if (sscanf(token2, "%lf", &mcmcSetup.finetunes.migTime) != 1) {
				fprintf(stderr,"Error: value for finetune-mig-time should be floating point number, got %s.\n", token2);
				numErrors++;
			}
		} else if(0 == strcmp("finetune-theta",token)) {
			if (sscanf(token2, "%lf", &mcmcSetup.finetunes.theta) != 1) {
				fprintf(stderr,"Error: value for finetune-theta should be floating point number, got %s.\n", token2);
				numErrors++;
			}
		} else if(0 == strcmp("finetune-mig-rate",token)) {
			if (sscanf(token2, "%lf", &mcmcSetup.finetunes.migRate) != 1) {
				fprintf(stderr,"Error: value for finetune-mig-rate should be floating point number, got %s.\n", token2);
				numErrors++;
			}
		} else if(0 == strcmp("finetune-tau",token)) {
			if (sscanf(token2, "%lf", &globalSetup.finetuneTaus) != 1) {
				fprintf(stderr,"Error: value for finetune-tau should be floating point number, got %s.\n", token2);
				numErrors++;
			}
		} else if(0 == strcmp("finetune-locus-rate",token)) {
			if (sscanf(token2, "%lf", &mcmcSetup.finetunes.locusRate) != 1) {
				fprintf(stderr,"Error: value for finetune-locus-rate should be floating point number, got %s.\n", token2);
				numErrors++;
			}
		} else if(0 == strcmp("finetune-mixing",token)) {
			if (sscanf(token2, "%lf", &mcmcSetup.finetunes.mixing) != 1) {
				fprintf(stderr,"Error: value for finetune-mixing should be floating point number, got %s.\n", token2);
				numErrors++;
			}
/*		} else if(0 == strcmp("finetune-admix",token)) {
			if (sscanf(token2, "%lf", &mcmcSetup.finetunes.admix) != 1) {
				fprintf(stderr,"Error: value for finetune-admix should be floating point number, got %s.\n", token2);
				numErrors++;
			}
*/		} else if (0 == strcmp("find-finetunes",token)) {
			if (0 == strcmp(token2, "TRUE")) {
				mcmcSetup.findFinetunes = 1;
			} else if(0 != strcmp(token2, "FALSE")) {
				fprintf(stderr,"Error: value of find-finetunes should be TRUE or FALSE, got '%s'.\n", token2);
				numErrors++;
			}
		} else if(0 == strcmp("find-finetunes-num-steps",token)) {
			if (sscanf(token2, "%d", &mcmcSetup.findFinetunesNumSteps) != 1 || mcmcSetup.findFinetunesNumSteps <= 0) {
				fprintf(stderr,"Error: value for find-finetunes-num-steps should be positive integer, got %s.\n", token2);
				numErrors++;
			}
		} else if(0 == strcmp("find-finetunes-samples-per-step",token)) {
			if (sscanf(token2, "%d", &mcmcSetup.findFinetunesSamplesPerStep) != 1 || mcmcSetup.findFinetunesSamplesPerStep <= 0) {
				fprintf(stderr,"Error: value for find-finetunes-samples-per-step should be positive integer, got %s.\n", token2);
				numErrors++;
			}
		} else {
			fprintf(stderr, "Error: argument '%s' is not accepted in GENERAL-INFO module.\n",token);
			numErrors++;
		}
	}// end of while(0 != strcmp("GENERAL-INFO-END",token))

	return numErrors;
}
/** end of readGeneralInfo **/



/***********************************************************************************
 *	readCurrentPops
 * 	- reads CURRENT-POPS part of control file
 *	- returns number of errors found
 ***********************************************************************************/
/* MARK: ADD OPTION TO READ IN POP AGE
         DONE. NEED TO CHECK !!
*/
int readCurrentPops(FILE* fctl) {
	int curPops = 0, pop;
	int numErrors = 0;
	char token[STRING_LENGTH];
	char* token2;
	char line[STRING_LENGTH];
	
	
	numErrors += expectNextToken(fctl,"CURRENT-POPS-START",token);
	if(feof(fctl)){
		fprintf(stderr, "Error: unexpected end of file before CURRENT-POPS-START.\n");
		return numErrors+1;
	}

	curPops = countTokens(fctl, "POP-START", "CURRENT-POPS-END", token);
	
	if(curPops <=0) {
		fprintf(stderr, "Error: could not find any POP items in CURRENT-POPS module.\n");
		return numErrors+1;
	}
	
	dataSetup.popTree = createPopTree(curPops);

	// allocate temporary array for print factors, might expand later on, with migrations, etc.
	mcmcSetup.printFactors = (double*)malloc((3*dataSetup.popTree->numCurPops - 2)*sizeof(double));
	if(mcmcSetup.printFactors == NULL) {
		fprintf(stderr, "\n Error: oom allocating printFactors array.\n");
		exit(-1);
	}
	
	dataSetup.numSamplesPerPop = (int*)malloc(dataSetup.popTree->numCurPops*sizeof(int));
	if(dataSetup.numSamplesPerPop == NULL) {
		fprintf(stderr, "\nError: oom allocating maxSamplesPerPop array.\n");
		exit(-1);
	}
		
	// start up sample name array
	dataSetup.maxSamples = 0;
	dataSetup.numSamples = 0;
	for(pop=0; pop<dataSetup.popTree->numCurPops; pop++) {
		numErrors += expectNextToken(fctl,"POP-START",token);
		if(feof(fctl)) {
			fprintf(stderr, "Error: unexpected end of file before POP-START.\n");
			return numErrors+1;
		}
		
		// initialize population arguments
		dataSetup.numSamplesPerPop[pop] = 0;
		dataSetup.popTree->pops[pop]->thetaPrior.alpha = globalSetup.alpha;
		dataSetup.popTree->pops[pop]->thetaPrior.beta  = globalSetup.beta;
		dataSetup.popTree->pops[pop]->age = 0.0;
		dataSetup.popTree->pops[pop]->isAncestralTo[pop] = 1;
		mcmcSetup.printFactors[pop] = globalSetup.printFactor;
		
		while(1) {
			// read token
			getNextToken(fctl,token);
	    	if (0 == strcmp("POP-END", token))
				break;  
   
			// read next token - value of some argument
			if(NULL == fgets(line,STRING_LENGTH-1,fctl)) {
				fprintf(stderr,"Error: unexpected end of file or other error inside POP module.\n");
				return numErrors+1;
			}
			
			// read sample name list
			if (0 == strcmp("samples", token)) {
				numErrors += readSampleLine(line,pop);
			} else {
				token2 = strtokCS(line, parseFileDelims);
				if(token2 == NULL) {
					fprintf(stderr,"Error: unable to read value for %s.\n",token);
					numErrors++;
					continue;
				}
				if(0 == strcmp("name",token)) {
					strncpy(dataSetup.popTree->pops[pop]->name, token2,  NAME_LENGTH-1);
				} else if(0 == strcmp("theta-print",token)) {
					if (sscanf(token2, "%lf", &mcmcSetup.printFactors[pop]) != 1) {
						fprintf(stderr,"Error: value for POP theta-print should be floating point number, got %s.\n", token2);
						numErrors++;
					}
				} else if(0 == strcmp("theta-alpha",token)) {
					if (sscanf(token2, "%lf", &dataSetup.popTree->pops[pop]->thetaPrior.alpha) != 1) {
						fprintf(stderr,"Error: value for POP theta-alpha should be floating point number, got %s.\n", token2);
						numErrors++;
					}
				} else if(0 == strcmp("theta-beta",token)) {
					if (sscanf(token2, "%lf", &dataSetup.popTree->pops[pop]->thetaPrior.beta) != 1) {
						fprintf(stderr,"Error: value for POP theta-beta should be floating point number, got %s.\n", token2);
						numErrors++;
					}
				} else if(0 == strcmp("age",token)) {
					if (sscanf(token2, "%lf", &dataSetup.popTree->pops[pop]->sampleAge) != 1) {
						fprintf(stderr,"Error: value for POP age should be floating point number, got %s.\n", token2);
						numErrors++;
					}
					token2 = strtokCS(NULL, parseFileDelims);
					if(strlen(token2)>1 || (token2[0] != 'f' && token2[0] != 'e')) {
						fprintf(stderr,"Error: POP age can be set to fixed (f) or estimated (e), not %s.\n", token2);
						numErrors++;
					}
					if(token2[0] == 'f') {
						dataSetup.popTree->pops[pop]->updateSampleAge = 0;
						if(dataSetup.popTree->pops[pop]->sampleAge != 0.0) {
							mcmcSetup.doMixing = 0;
						}
					} else {
						dataSetup.popTree->pops[pop]->updateSampleAge = 1;
					}					
				} else {
					fprintf(stderr, "Error: argument '%s' is not accepted in POP module of CURRENT-POPS.\n",token);
					numErrors++;
				}
			}
		}// end of while(1)

		// check that name is set
		if(dataSetup.popTree->pops[pop]->name[0] == '\0') {
			fprintf(stderr, "Error: no name is given for current pop %d.\n", pop+1);
			numErrors++;
		}
	}// end of for(pop)
    
	numErrors += expectNextToken(fctl,"CURRENT-POPS-END",token);
	if(feof(fctl)) {
		fprintf(stderr, "Error: unexpected end of file before CURRENT-POPS-END.\n");
		numErrors++;
	}

	return numErrors;
}
/** end of readCurrentPops **/



/***********************************************************************************
 *	readAncestralPops
 * 	- reads ANCESTRAL-POPS part of control file
 *	- returns number of errors found
 ***********************************************************************************/
int readAncestralPops(FILE* fctl) {
	int numErrors = 0;
	char token[STRING_LENGTH];
	char* token2;
	char line[STRING_LENGTH];

	int pop, pop1, pop2, son;

	
	numErrors += expectNextToken(fctl,"ANCESTRAL-POPS-START",token);
	if(feof(fctl)){
		fprintf(stderr, "Error: unexpected end of file before ANCESTRAL-POPS-START.\n");
		return numErrors+1;
	}

	// create array for per-population TAU split time finetunes, and initialize with global finetune
	mcmcSetup.finetunes.taus = (double*)malloc(sizeof(double) * dataSetup.popTree->numPops);
	for(pop=0; pop<dataSetup.popTree->numPops;pop++) {
		mcmcSetup.finetunes.taus[pop] = globalSetup.finetuneTaus;
  	}
  	
				
	for(pop=dataSetup.popTree->numCurPops; pop<dataSetup.popTree->numPops; pop++) {
		numErrors += expectNextToken(fctl,"POP-START",token);
		if(feof(fctl)) {
			fprintf(stderr, "Error: could not find module POP for ancestral population %d (expecting %d populations total).\n",pop+1, dataSetup.popTree->numPops);
			return numErrors+1;
		}
		
		// initialize population arguments
		dataSetup.popTree->pops[pop]->thetaPrior.alpha = globalSetup.alpha;
		dataSetup.popTree->pops[pop]->thetaPrior.beta  = globalSetup.beta;
		dataSetup.popTree->pops[pop]->agePrior.alpha   = globalSetup.alpha;
		dataSetup.popTree->pops[pop]->agePrior.beta    = globalSetup.beta;
		dataSetup.popTree->pops[pop]->agePrior.sampleStart = -1;
		dataSetup.popTree->pops[pop]->isAncestralTo[pop] = 1;
		mcmcSetup.printFactors[pop] = globalSetup.printFactor;	// theta
		mcmcSetup.printFactors[pop + dataSetup.popTree->numCurPops-1] = globalSetup.printFactor; //tau
		
		while(1) {
			// read token
			getNextToken(fctl,token);
	    	if (0 == strcmp("POP-END", token))
				break;  
   
			// read next token - value of some argument
			if(NULL == fgets(line,STRING_LENGTH-1,fctl)) {
				fprintf(stderr,"Error: unexpected end of file or other error inside POP module.\n");
				return numErrors+1;
			}
			
			token2 = strtokCS(line, parseFileDelims);
			if(token2 == NULL) {
				fprintf(stderr,"Error: unable to read value for %s.\n",token);
				numErrors++;
				continue;
			}
			
			if(0 == strcmp("name",token)) {
				strncpy(dataSetup.popTree->pops[pop]->name, token2,  NAME_LENGTH-1);
			} else if(0 == strcmp("children",token)) {
				for(son=0; son<2; son++) {
					pop1 = getPopIdByName(dataSetup.popTree, token2);
					if(pop1<0) {
						fprintf(stderr, "Error: pop child name '%s' unrecognized for ancestral pop %d.\n", token2, pop+1);
						numErrors++;
					}
					
					dataSetup.popTree->pops[pop]->sons[son] = dataSetup.popTree->pops[pop1];
					if(dataSetup.popTree->pops[pop1]->father != NULL) {
						fprintf(stderr, "Error: population %s already has a parent defined already (%d in addition to %d).\n",
							dataSetup.popTree->pops[pop1]->name, dataSetup.popTree->pops[pop1]->father->id+1,pop+1);
						numErrors++;
					} else {
						dataSetup.popTree->pops[pop1]->father = dataSetup.popTree->pops[pop];
					}
					// copy ancestral array of child to parent
					for(pop2=0; pop2<dataSetup.popTree->numPops; pop2++) {
						if(dataSetup.popTree->pops[pop1]->isAncestralTo[pop2]) {
							dataSetup.popTree->pops[pop]->isAncestralTo[pop2] = 1;
						}
					}

					token2 = strtokCS(NULL, parseFileDelims);
					if(son<1 && token2 == NULL) {
						fprintf(stderr,"Error: second child of ancestral pop %d is missing.\n",pop+1);
						numErrors++;
						break;
					}
				}// end of for(son)
			} else if(0 == strcmp("theta-print",token)) {
				if (sscanf(token2, "%lf", &mcmcSetup.printFactors[pop]) != 1) {
					fprintf(stderr,"Error: value for POP theta-print should be floating point number, got %s.\n", token2);
					numErrors++;
				}
			} else if(0 == strcmp("theta-alpha",token)) {
				if (sscanf(token2, "%lf", &dataSetup.popTree->pops[pop]->thetaPrior.alpha) != 1) {
					fprintf(stderr,"Error: value for POP theta-alpha should be floating point number, got %s.\n", token2);
					numErrors++;
				}
			} else if(0 == strcmp("theta-beta",token)) {
				if (sscanf(token2, "%lf", &dataSetup.popTree->pops[pop]->thetaPrior.beta) != 1) {
					fprintf(stderr,"Error: value for POP theta-beta should be floating point number, got %s.\n", token2);
					numErrors++;
				}
			} else if(0 == strcmp("tau-print",token)) {
				if (sscanf(token2, "%lf", &mcmcSetup.printFactors[dataSetup.popTree->numCurPops+pop-1]) != 1) {
					fprintf(stderr,"Error: value for POP tau-print should be floating point number, got %s.\n", token2);
					numErrors++;
				}
			} else if(0 == strcmp("tau-alpha",token)) {
				if (sscanf(token2, "%lf", &dataSetup.popTree->pops[pop]->agePrior.alpha) != 1) {
					fprintf(stderr,"Error: value for POP tau-alpha should be floating point number, got %s.\n", token2);
					numErrors++;
				}
			} else if(0 == strcmp("tau-beta",token)) {
				if (sscanf(token2, "%lf", &dataSetup.popTree->pops[pop]->agePrior.beta) != 1) {
					fprintf(stderr,"Error: value for POP tau-beta should be floating point number, got %s.\n", token2);
					numErrors++;
				}
			} else if(0 == strcmp("tau-initial",token)) {
				if (sscanf(token2, "%lf", &dataSetup.popTree->pops[pop]->agePrior.sampleStart) != 1) {
					fprintf(stderr,"Error: value for POP tau-initial should be floating point number, got %s.\n", token2);
					numErrors++;
				}
			} else if(0 == strcmp("finetune-tau",token)) {
				if (sscanf(token2, "%lf", &mcmcSetup.finetunes.taus[pop]) != 1) {
					fprintf(stderr,"Error: value for POP finetune-tau should be floating point number, got %s.\n", token2);
					numErrors++;
				}
				if (globalSetup.finetuneTaus > 0.0 && verbose) {
					fprintf(stderr, "Warning: overriding global tau finetune %lf for ancestral population %d, with population specific finetune %lf.\n", 
										globalSetup.finetuneTaus, pop+1, mcmcSetup.finetunes.taus[pop]);
				}

			} else {
				fprintf(stderr, "Error: argument '%s' is not accepted in POP module of ANCESTRAL-POPS.\n",token);
				numErrors++;
			}
		}// end of while(1)

		// check that name is set
		if(dataSetup.popTree->pops[pop]->name[0] == '\0') {
			fprintf(stderr, "Error: no name is given for ancestral pop %d.\n", pop+1);
			numErrors++;
		}
		
		// check that both children are set
		for(son=0; son<2; son++) {
			if(dataSetup.popTree->pops[pop]->sons[son] == NULL) {
				fprintf(stderr, "Error: son #%d is not set for ancestral pop %d.\n", son+1, pop+1);
				numErrors++;
			}
		}
		
	}// end of for(pop)
	
    if(pop != dataSetup.popTree->numPops) {
		fprintf(stderr, "Error: with %d currnet pops, expecting to see a total of %d pops, but reading %d.\n", dataSetup.popTree->numCurPops, dataSetup.popTree->numPops, pop);
		numErrors++;
	}
	// set root pop
	dataSetup.popTree->rootPop = dataSetup.popTree->numPops-1;
		
	numErrors += expectNextToken(fctl,"ANCESTRAL-POPS-END",token);
	if(feof(fctl)) {
		fprintf(stderr, "Error: unexpected end of file before ANCESTRAL-POPS-END.\n");
		numErrors++;
	}
	

	return numErrors;
}
/** end of readAncestralPops **/



/***********************************************************************************
 *	readMigrationBands
 * 	- reads MIG-BANDS part of control file
 *	- returns number of errors found
 ***********************************************************************************/
int readMigrationBands(FILE* fctl) {
	int numErrors = 0;
	char token[STRING_LENGTH];
	char* token2;
	char line[STRING_LENGTH];

	int sourcePop, targetPop, migBand;

	
	numErrors += expectNextToken(fctl,"MIG-BANDS-START",token);
	if(feof(fctl)){
		// do not consider eof at this point to be an error
		return numErrors;
	}

	// determine number of migration bands
	dataSetup.popTree->numMigBands = countTokens(fctl, "BAND-START", "MIG-BANDS-END", token);

	//reallocate space for print factors
	mcmcSetup.printFactors = (double*)realloc(mcmcSetup.printFactors, (3*dataSetup.popTree->numCurPops - 2 + dataSetup.popTree->numMigBands)*sizeof(double));
	if(mcmcSetup.printFactors == NULL) {
		fprintf( stderr, "Error: oom reallocating mcmcSetup.printFactors for mig bands.\n");
		exit(-1);
	}

			
	for(migBand=0; migBand<dataSetup.popTree->numMigBands; migBand++) {
		numErrors += expectNextToken(fctl,"BAND-START",token);
		if(feof(fctl)) {
			fprintf(stderr, "Error: unexpected end of file before satrt of mig band %d.\n", migBand);
			return numErrors;
		}
		
		// initialize mig band arguments
		dataSetup.popTree->migBands[migBand].sourcePop = -1;
		dataSetup.popTree->migBands[migBand].targetPop = -1;
		dataSetup.popTree->migBands[migBand].migRatePrior.alpha = globalSetup.migAlpha;
		dataSetup.popTree->migBands[migBand].migRatePrior.beta  = globalSetup.migBeta;
		mcmcSetup.printFactors[migBand + 2*dataSetup.popTree->numPops - dataSetup.popTree->numCurPops] = globalSetup.migFactor;
       
		while(1) {
			// read token
			getNextToken(fctl,token);
	    	if (0 == strcmp("BAND-END", token))
				break;  
   
			// read next token - value of some argument
			if(NULL == fgets(line,STRING_LENGTH-1,fctl)) {
				fprintf(stderr,"Error: unexpected end of file or other error inside BAND module.\n");
				return numErrors+1;
			}
			
			token2 = strtokCS(line, parseFileDelims);
			if(token2 == NULL) {
				fprintf(stderr,"Error: unable to read value for %s.\n",token);
				numErrors++;
				continue;
			}
			
			if(0 == strcmp("source",token)) {
				sourcePop = getPopIdByName(dataSetup.popTree, token2);
				if(sourcePop<0) {
					fprintf(stderr, "Error: invalid name '%s' for source pop of mig-band %d.\n", token2, migBand+1);
					numErrors++;
				} else {
					dataSetup.popTree->migBands[migBand].sourcePop = sourcePop;
				}
			} else if(0 == strcmp("target",token)) {
				targetPop = getPopIdByName(dataSetup.popTree, token2);
				if(targetPop<0) {
					fprintf(stderr, "Error: invalid name '%s' for target pop of mig-band %d.\n", token2, migBand+1);
					numErrors++;
				} else {
					dataSetup.popTree->migBands[migBand].targetPop = targetPop;
				}
			} else if(0 == strcmp("mig-rate-print",token)) {
				if (sscanf(token2, "%lf", &mcmcSetup.printFactors[3*dataSetup.popTree->numCurPops-2+migBand]) != 1) {
					fprintf(stderr,"Error: value for mig-rate-print should be floating point number, got %s.\n", token2);
					numErrors++;
				}
			} else if(0 == strcmp("mig-rate-alpha",token)) {
				if (sscanf(token2, "%lf", &dataSetup.popTree->migBands[migBand].migRatePrior.alpha) != 1) {
					fprintf(stderr,"Error: value for mig-rate-alpha should be floating point number, got %s.\n", token2);
					numErrors++;
				}
			} else if(0 == strcmp("mig-rate-beta",token)) {
				if (sscanf(token2, "%lf", &dataSetup.popTree->migBands[migBand].migRatePrior.beta) != 1) {
					fprintf(stderr,"Error: value for mig-rate-beta should be floating point number, got %s.\n", token2);
					numErrors++;
				}
			} else {
				fprintf(stderr, "Error: argument '%s' is not accepted in BAND module.\n",token);
				numErrors++;
			}
		}// end of while(1)

		// check that source and target pops are set and not ancestral to one another
		if(dataSetup.popTree->migBands[migBand].sourcePop == -1) {
			fprintf(stderr, "Error: source population for migration band %d was not defined.\n", migBand+1);
			return numErrors+1;
		}
		if(dataSetup.popTree->migBands[migBand].targetPop == -1) {
			fprintf(stderr, "Error: target population for migration band %d was not defined.\n", migBand+1);
			return numErrors+1;
		}
		if(dataSetup.popTree->pops[ dataSetup.popTree->migBands[migBand].sourcePop ]->isAncestralTo[ dataSetup.popTree->migBands[migBand].targetPop ]) {
			fprintf(stderr, "Error: source pop for migration band %d is an ancestor of its target pop.\n\t\tMigration bands can only be placed between two populations which may have co-occured.\n",migBand+1);
			return numErrors+1;
		}      
		if(dataSetup.popTree->pops[ dataSetup.popTree->migBands[migBand].targetPop ]->isAncestralTo[ dataSetup.popTree->migBands[migBand].sourcePop ]) {
			fprintf(stderr, "Error: target pop for migration band %d is an ancestor of its source pop.\n\t\tMigration bands can only be placed between two populations which may have co-occured.\n",migBand+1);
         printPopulationTree(dataSetup.popTree, stdout, 0);
			return numErrors+1;
		}
	}// end of for(migBand)
	
	numErrors += expectNextToken(fctl,"MIG-BANDS-END",token);
	if(feof(fctl)) {
		fprintf(stderr, "Error: unexpected end of file before MIG-BANDS-END.\n");
		numErrors++;
	}
	
	// read until end of file
	getNextToken(fctl,token);
	while(!feof(fctl)) {
		numErrors++;
		fprintf(stderr, "Error: ignoring token '%s' after MIG-BANDS.\n",token);
		getNextToken(fctl,token);
	}

	return numErrors;
}
/** end of readMigrationBands **/



/***********************************************************************************
 *	readSampleLine
 * 	- reads line of sample names in current population module
 *	- sampleLine: string for line
 *	- pop: index of pop
 *	- numExistingSamples
 *	- updates saved data in dataSetup structure
 *	- returns number of errors found
 ***********************************************************************************/
int readSampleLine(char* sampleLine, int pop) {
	char saveLine[STRING_LENGTH];
	char* str;
	unsigned short nameORformat;
	
	int numErrors = 0;
	
	int sample, sampleIndex = dataSetup.numSamples;
	
	if(dataSetup.numSamplesPerPop[pop] > 0) {
		fprintf(stderr, "Error: observed two sample lists for pop %d.\n",pop+1);
		return numErrors+1;
	}
	
	// copy line for second parsing round
	strncpy(saveLine,sampleLine, STRING_LENGTH-1);
	
	// determine number of samples - one per haploid, two per diploid
	str = strtokCS(sampleLine, parseFileDelims);
	nameORformat = 1;
	while(str != NULL) {
		if(nameORformat == 0) {
			if(strlen(str) != 1 || (str[0] != 'h' && str[0] != 'd')) {
				fprintf(stderr, "Error: faulty format %s for sample pop %d. Expected h or d.\n",str,pop+1);
				numErrors++;
			}
			dataSetup.numSamplesPerPop[pop]++;
			if(str[0] == 'd') {
				dataSetup.numSamplesPerPop[pop]++;
			}
		}
		nameORformat = !nameORformat;						
		str = strtokCS(NULL, parseFileDelims);
	}// end of while(str != NULL)
				
	if(nameORformat != 1) {
		fprintf(stderr, "Error: uneven terms in sample line for pop %d.\n",pop+1);
		numErrors++;
	}
	if(dataSetup.numSamplesPerPop[pop] < 1) {
		fprintf(stderr,"Error: no samples provided for pop %d.\n", pop+1);
		numErrors++;
	}
				
	// if need to increase allocation for sample data structures
	if(dataSetup.numSamples+dataSetup.numSamplesPerPop[pop] > dataSetup.maxSamples) {
		dataSetup.maxSamples += 15+dataSetup.numSamplesPerPop[pop];
		if(dataSetup.sampleNames == NULL) {
			dataSetup.sampleNames    = (char**)malloc(dataSetup.maxSamples*sizeof(char*));
			dataSetup.sampleNames[0] = (char*)malloc(dataSetup.maxSamples*NAME_LENGTH*sizeof(char));
		} else {
			dataSetup.sampleNames    = (char**)realloc(dataSetup.sampleNames, dataSetup.maxSamples*sizeof(char*));
			dataSetup.sampleNames[0] = (char*)realloc(dataSetup.sampleNames[0], dataSetup.maxSamples*NAME_LENGTH*sizeof(char));
		}	
		if(dataSetup.sampleNames == NULL || dataSetup.sampleNames[0] == NULL) {
			fprintf(stderr, "Error: oom allocating sampleNames pointers.\n");
			exit(-1);
		}
		for(sample=1; sample<dataSetup.maxSamples; sample++) {
			dataSetup.sampleNames[sample] = dataSetup.sampleNames[sample-1] + NAME_LENGTH;
		}
	}
				
	// read sample names
	str = strtokCS(saveLine, parseFileDelims);
	nameORformat = 1;
	while(str != NULL) {
		if(nameORformat == 1) {
			strncpy(dataSetup.sampleNames[sampleIndex], str, NAME_LENGTH-1);
			sampleIndex++;
		} else if(str[0] == 'd') {
			dataSetup.sampleNames[sampleIndex][0] = '\0';
			sampleIndex++; // save two samples for diploid
		}
		nameORformat = !nameORformat;						
		str = strtokCS(NULL, parseFileDelims);
	}// end of while(str != NULL)
	
	if(sampleIndex - dataSetup.numSamples != dataSetup.numSamplesPerPop[pop]) {
		fprintf(stderr,"Error: expecting to see %d samples in pop %d, but was able to read only %d.\n", 
							dataSetup.numSamplesPerPop[pop], pop+1, sampleIndex - dataSetup.numSamples);
		numErrors++;
	}
	
	dataSetup.numSamples = sampleIndex;

	return numErrors;
}
/** end of readSampleLine **/



/***********************************************************************************
 *	parseSampleNames
 *	- looks for samples that appear in several populations
 *	- if admixture modeling is enabled, then allows two appearances. Otherwise, only one.
 *	- also setsp up admixture data structures if needed
 *	- returns number of errors found
 ***********************************************************************************/
int parseSampleNames() {
	int pop, sample, sample1, sampleMatch, skip, numOccur;
	int numErrors = 0;
	int* samples2pops = NULL;
	
	admixed_samples.number = 0;
	// allocate memory in admixed_samples data structures
	if(mcmcSetup.allowAdmixture) {
		samples2pops = (int*)malloc(dataSetup.numSamples*sizeof(int));
		admixed_samples.popPairs = (int**)malloc(dataSetup.numSamples*sizeof(int*));
		admixed_samples.samples = (int*)malloc(dataSetup.numSamples*sizeof(int));
		admixed_samples.popPairs[0] = (int*)malloc(2*dataSetup.numSamples*sizeof(int));
		admixed_samples.index = (int*)malloc(dataSetup.numSamples*sizeof(int));

		// initialize samples2pops array
		sample1 = 0;
		sample = 0;
		for(pop=0; pop<dataSetup.popTree->numCurPops; pop++) {
			sample1+= dataSetup.numSamplesPerPop[pop];
			for( ; sample<sample1; sample++) {
				samples2pops[sample] = pop;
			}// end of for(sample)
		}// end of for(pop)
	}
	
	
	
	for(sample=0; sample<dataSetup.numSamples; sample++) {
		if(dataSetup.sampleNames[sample][0] == '\0') {
			continue;
		}
		numOccur = 1;
		sampleMatch = -1;
		for(sample1=sample+1; sample1<dataSetup.numSamples; sample1++) {
			if(0 == strcmp(dataSetup.sampleNames[sample] , dataSetup.sampleNames[sample1])) {
				if(numOccur == 1) {
					sampleMatch = sample1;
				}
				numOccur++;
				break;
			}
		}// end of for(sample1)
		
		if(numOccur > 2) {
			fprintf(stderr, "Error: sample %s appears more than twice in control file.\n", dataSetup.sampleNames[sample]);
			numErrors++;
			continue;
		}
		if(!mcmcSetup.allowAdmixture) {
			if(numOccur > 1) {
				fprintf(stderr, "Error: admixture is not allowed, so sample %s cannot appear more than once in control file.\n",dataSetup.sampleNames[sample]);
				numErrors++;
			}
			continue;
		}
		
		
		if(sampleMatch >= 0) {
			skip = 1;
			admixed_samples.samples [admixed_samples.number] = sample;
			admixed_samples.popPairs[admixed_samples.number][0] = samples2pops[sample];
			admixed_samples.popPairs[admixed_samples.number][1] = samples2pops[sampleMatch];
			admixed_samples.number++;
			if(admixed_samples.number < dataSetup.numSamples) {
				admixed_samples.popPairs[admixed_samples.number] = admixed_samples.popPairs[admixed_samples.number-1] + 2;
			}
			// for diploid sample, define both haploids as admixed
			if(dataSetup.sampleNames[sample+1][0] == '\0') {
				skip++;
				admixed_samples.samples [admixed_samples.number] = sample+1;
				admixed_samples.popPairs[admixed_samples.number][0] = samples2pops[sample];
				admixed_samples.popPairs[admixed_samples.number][1] = samples2pops[sampleMatch];
				admixed_samples.number++;
				if(admixed_samples.number < dataSetup.numSamples) {
					admixed_samples.popPairs[admixed_samples.number] = admixed_samples.popPairs[admixed_samples.number-1] + 2;
				}
				if(sampleMatch>=dataSetup.numSamples || dataSetup.sampleNames[sampleMatch+1][0] != '\0') {
					fprintf(stderr, "Error: sample %s is admixed, but first copy is defined as diploid, and second as haploid.\n", dataSetup.sampleNames[sample]);
					numErrors++;
					continue;
				}
			} else {
				if(sampleMatch<dataSetup.numSamples && dataSetup.sampleNames[sampleMatch+1][0] == '\0') {
					fprintf(stderr, "Error: sample %s is admixed, but first copy is defined as haploid, and second as diploid.\n", dataSetup.sampleNames[sample]);
					numErrors++;
					continue;
				}
			}
			// erase later mention of sample and move back all samples in list
			dataSetup.numSamplesPerPop[ samples2pops[sample] ] -= skip;
			dataSetup.numSamples -= skip;
			for(sample1=sampleMatch; sample1<dataSetup.numSamples; sample1++) {
				samples2pops[sample1] = samples2pops[sample1+skip];
				strcpy(dataSetup.sampleNames[sample1],dataSetup.sampleNames[sample1+skip]); 
			}// end of for(sample1)
		}// end of if(sampleMatch >= 0)
	}// end of for(sample)

	if(mcmcSetup.allowAdmixture) {
		free(samples2pops);
	}

	return numErrors;
}
/** end of parseSampleNames **/



/***************************************************************************************************************/
/******                                        END OF FILE                                                ******/
/***************************************************************************************************************/
