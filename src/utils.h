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
#define rndexp(mean) (-(mean)*log(rndu()))

void test(const char* message);

double rndu (void);
double reflect(double x, double a, double b);
void setSeed (unsigned int seed);	
double rndnormal(void);
double rnd2normal8(void);
double rndgamma(double s);
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

#endif
