/*
 * GPhoCS.h
 *
 *  Created on: Feb 4, 2017
 *      Author: ron
 */

#ifndef SRC_GPHOCS_H_
#define SRC_GPHOCS_H_

#include "LocusDataLikelihood.h"
#include <getopt.h>

// --- CONSTANTS --------------------------------------------------------------

#define LOG_STEPS_NOT
#define CHECKALL_NOT
//#define CHECKALL

#define NUM_TYPES                    5
#define TARGET_ACCEPTANCE_PERCENT    35
#define TARGET_ACCEPTANCE_RANGE      5
#define FINETUNE_RESOLUTION          0.0000001
#define MAX_FINETUNE                 10
#define ACCEPTANCE_FUDGE             2

#define GPHOCS_VERSION_NUM          "1.3"
#define GPHOCS_VERSION_DATE         "Jan. 2017"

int typeCount[NUM_TYPES];

// --- GLOBAL DATA STRUCTURES -------------------------------------------------

// Data setup. "Singleton"
struct DATA_STATE {
		              	  	  	  	  	   // average log-likelihood per genealogy of data given pop tree
		double logLikelihood;			   //  (ln[P(X|Z)]+ln[P(Z|M)])/numLoci

										   // log likelihood (not averaged) of data given all genealogies:
		double dataLogLikelihood;          //   ln[P(X|Z,M,T)]

										   // log likelihood of all genealogies given model & parameters -
		double genealogyLogLikelihood;	   //   ln[P(Z|M,T)]

		double rateVar;                    // the actual variance in locus-specific
                                           // mutation rate
		LocusData** lociData;              // array of LocusData data structures
                                           // (of length numLoci).
                                           // (allocated in processAlignments)
} dataState;

// Miscellaneous statistics. "Singleton"
struct MISC_STATS {
		int rubberband_mig_conflicts;      // number of rubber band conflicts
                                           // with migration nodes

		int spr_zero_targets;              // number of times an SPR event
                                           // encounters zero target edges

		int not_enough_migs;               // number of times not enough pre-
                                           // allocated space for migration nodes

// int small_interval;                     // very small interval for moving
// coalescent event or migration event

// double spr_lnld_disc;                   // the size of the smallest discrepancy
//in log-likelihood computation
} misc_stats;

static struct option long_options[] = { { "help", no_argument, 0, 'h' },
                                        { "verbose", no_argument, 0, 'v' },
                                        { "nthreads", no_argument, 0, 'n' },
                                        { 0, 0, 0, 0 } };

// --- FUNCTION DECLARATIONS --------------------------------------------------

void printUsage(char *programName);
int processAlignments();
int readRateFile(const char* fileName);
int initLociWithoutData();
void printParamVals(double paramVals[], int startParam, int endParam, FILE* o);
int recordTypes();
int recordParamVals(double paramVals[]);
int performMCMC();
void printGenealogyAndExit(int gen, int errStatus);
int freeAllMemory();

// Sampling functions
int UpdateGB_InternalNode(double finetune);  // step 1: update coalescent times
int UpdateGB_MigrationNode(double finetune); // step 2: update migration times
int UpdateGB_MigSPR();                       // step 3: update genealogy struct
int UpdateTheta(double finetune);            // step 4: No to MT
int UpdateMigRates(double finetune);         // step 5: No to MT,
                                             //         update migration bands

void UpdateTau(double *finetunes,            // step 6: update tau
               int *accepted);               //         More difficult to MT,
                                             //         most time consuming

void UpdateSampleAge(double *finetunes,      // similar to update Tau.
                     int *accepted);         // Modifies the time

int UpdateLocusRate(double finetune);
int UpdateAdmixCoeffs(double finetune);      // Shall be skipped ("ledaleg")
int mixing(double finetune);

void allocateAllMemory() ; // TODO - tidy header file
int combStatsActivated();
double getLogPrior() ;



#endif /* SRC_GPHOCS_H_ */
