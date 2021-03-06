#ifndef UTILS_H   
#define UTILS_H
/**
   \file utils.h 
   Utility definitions (Miscellaneous).
	
	Contains some utility functions and definitions.
	
*/


#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <float.h>
#include <time.h>


#define STRING_LENGTH 	500
#define NAME_LENGTH		200
#define min2(a,b)  ((a)<(b)?(a):(b))
#define max2(a,b)  ((a)>(b)?(a):(b))
#define swap2(a,b) ({a += b ; b = a-b ; a -= b ;})
#define rndexp(nLocusIdx, mean) (-(mean)*log(rndu(nLocusIdx)))

#define	ROOT_SLACK	0.001		// default length of top-most interval

double reflect(double x, double a, double b);

//---- Random Generator related stuff -----------------------------------------
#define RAND_GENERAL_SLOT (RndCtx.nOfSlots-1)
typedef struct _RandGeneratorContext
{
  int nOfSlots;

  unsigned int* rndu_z;
  unsigned int* rndu_w;
  unsigned int* rndu_x;
  unsigned int* rndu_y;

  double* m2s2_kernel;
  double* m2N_kernel;
  double* s2N_kernel;

  double* rndgamma2_b;
  double* rndgamma2_h;
  double* rndgamma2_ss;

  double* rndgamma1_a;
  double* rndgamma1_p;
  double* rndgamma1_uf;
  double* rndgamma1_ss;
  double* rndgamma1_d;
} RandGeneratorContext;

void initRandomGenerator( int nNumLoci, unsigned int seed );
double rndnormal( int nLocusIdx );
double rnd2normal8( int nLocusIdx );
double rndu( int nLocusIdx );
double rndgamma( int nLocusIdx, double s );
double rndgamma1( int nLocusIdx, double s );
double rndgamma2( int nLocusIdx, double s );
//-----------------------------------------------------------------------------

double PointChi2 (double prob, double v);
char* printtime(char timestr[]);
void starttime(void);
void resetBooleanArray(unsigned short* booleanArray, int arrayLength);
void turnOnBooleanArray(unsigned short* booleanArray, int arrayLength);
double logSumOfExponents(double* logArray, int arrayLength, unsigned short* indsToSum);


void mergeArrays(double array[], int numEntries, int borderPoint, double tmpArray[]);
double mergeSort(double array[], int numEntries, double tmpArray[]);

void flushLine(FILE* readFile);
int readStringFromFile(FILE *stream, int bufferLength, char *destination);
char * strtokCS( char * str, const char * delimiters);

extern int debug;
extern int verbose;
extern char parseFileDelims[];



typedef struct {
	time_t start;
	time_t accumulated;
} TIMER_METHOD;
typedef struct {
	TIMER_METHOD UpdateGB_InternalNode;
	TIMER_METHOD UpdateGB_MigrationNode;
	TIMER_METHOD UpdateGB_MigSPR;
	TIMER_METHOD UpdateTau;
	TIMER_METHOD UpdateMigRates;
	TIMER_METHOD mixing;
	TIMER_METHOD UpdateTheta;
	TIMER_METHOD UpdateSampleAge;
	TIMER_METHOD UpdateLocusRate;
	TIMER_METHOD UpdateAdmixCoeffs;
	TIMER_METHOD MCMCIterations;
} TIMERS;

enum METHOD_NAME
{
	T_UpdateGB_InternalNode,
	T_UpdateGB_MigrationNode,
	T_UpdateGB_MigSPR,
	T_UpdateTau,
	T_UpdateMigRates,
	T_mixing,
	T_UpdateTheta,
	T_UpdateSampleAge,
	T_UpdateLocusRate,
	T_UpdateAdmixCoeffs,
	T_MCMCIterations
} ;


void setEndTimeMethod(enum METHOD_NAME method);
void printMethodTimes();
void setStartTimeMethod(enum METHOD_NAME method);
char *printtime_i(int t , char timestr[]);
#endif
