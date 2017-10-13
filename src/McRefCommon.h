#ifndef SRC_MCREFCOMMON_H_
#define SRC_MCREFCOMMON_H_

#include "DataLayer.h"
#include "LocusDataLikelihood.h"

#define TRUE 1
#define FALSE 0

#define DEBUG_COMB_STATS FALSE
#define DEBUG_TAU_BOUNDS TRUE
#define SET_COMB_AGE_TO_ZERO FALSE

double calculateCoalStats(double *elapsed_times, int *num_lineages, int size);

int isLeafPop(int pop);

int areChildrenLeaves(int pop);

int isAncestralTo(int father, int son);

int getSon(int pop, int SON);

double getPopAge(int pop);

int getSourcePop(int mig);

int getTargetPop(int mig);

bool hasNextEvent(EventChain chain, int event);

char *getPopName(int popId);

bool areAlmostEqual(double eventAge, double combAge);

double relativeDistance(double dbl1, double dbl2);

double requiredRelativePrecision();

int lca(int pop1, int pop2);

int getPopFather(int popId);

bool isLeafNode(LikelihoodNode *node);

int getNodePop(int nodeId, int gen);

LikelihoodNode *getNode(int nodeId, int gen);

int numNodes();


/*=== Debug functions ===*/
const char *getEventTypeName(int eventType);

char *getMigName(int mig);

char *concat(const char *s1, const char *s2); // NOTE - this should not be used in production. the memory isn't released

#endif