#include <stdio.h>


#ifndef SRC_COMBSTATSHEADER_H_
#define SRC_COMBSTATSHEADER_H_


void printCombStats(int iteration, FILE* file);

void printCombStatsHeader(FILE* file);
void printSpecificCombHeader(int clade, FILE* file);
void printSpecificPopHeader(int pop, FILE* file);
void printSpecificMigHeader(int mig_band, FILE* file);
void printCombStats(int iteration, FILE* file);
void printSpecificCombStats(int clade, FILE* file);
void printSpecificPopStats(int pop, FILE* file);
void printSpecificMigBandStats(int mig_band, FILE* file);


#endif /* SRC_COMBSTATSHEADER_H_ */
