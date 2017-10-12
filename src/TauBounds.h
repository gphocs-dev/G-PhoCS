#include <stdio.h>


#ifndef G_PHOCS_TAUBOUNDS_H
#define G_PHOCS_TAUBOUNDS_H

#include "LocusDataLikelihood.h"

void initializeTauBounds();

void initializeLcaPops(int gen);

void updateTauBoundsOfDescendants(int pop, double bound);

void allocateTauBoundsMem();

void calculateTauBounds1();

void calculateLociTauBounds(int nodeId, int gen);

void calculateTauBounds2();

void printTauBoundsHeader(FILE *file);

void printTauBounds(int iteration, FILE *file);

//TODO - move to util

int lca(int pop1, int pop2);

int getPopFather(int popId);

bool isLeafNode(LikelihoodNode *node);

int getLeafNodePop(int nodeId, int gen);

LikelihoodNode *getNode(int nodeId, int gen);

int numNodes();

#endif //G_PHOCS_TAUBOUNDS_H
