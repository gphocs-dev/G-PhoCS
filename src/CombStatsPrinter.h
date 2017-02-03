#include <stdlib.h>


#ifndef SRC_COMBSTATSHEADER_H_
#define SRC_COMBSTATSHEADER_H_

/** clade stats printing functions **/

void printCladeStats(int iteration, FILE* file);


void printCladeStatsHeader(FILE* file);
void printSpecificCladeHeader(int clade, FILE* file);
void printSpecificPopHeader(int pop, FILE* file);
void printSpecificMigHeader(int mig_band, FILE* file);
void printCladeStats(int iteration, FILE* file);
void printSpecificCladeStats(int clade, FILE* file);
void printSpecificPopStats(int pop, FILE* file);
void printSpecificMigBandStats(int mig_band, FILE* file);


#endif /* SRC_COMBSTATSHEADER_H_ */
