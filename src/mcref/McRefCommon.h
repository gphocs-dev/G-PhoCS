#ifndef SRC_MCREFCOMMON_H_
#define SRC_MCREFCOMMON_H_

#include "../DataLayer.h"
#include "../LocusDataLikelihood.h"

#define TRUE 1
#define FALSE 0

#define DEBUG_COMB_STATS FALSE
#define DEBUG_TAU_BOUNDS FALSE
#define SET_COMB_AGE_TO_ZERO FALSE

double calculateCoalStats(double *elapsed_times, int *num_lineages, int size);

int isLeafPop(int pop);

int isLeafNode(int node, int gen);

int areChildrenLeaves(int pop);

int isAncestralTo(int ancestor, int descendant);

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

int getFather(int popId);

int getNodePop(int nodeId, int gen);

LikelihoodNode *getNode(int nodeId, int gen);

extern int **lca_pops; // "cache" of LCA of all pops. for example, lca_pops[1][2] is the lca of pops 1 & 2. Filled once by computeLca().
void computeLcas();
int migLcaPop(int nodeId, int gen, int defaultLcaPop);
int getMigNodeAbove(int nodeId, int gen, double requiredAge);

#endif