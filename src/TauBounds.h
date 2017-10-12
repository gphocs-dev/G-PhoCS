#include <stdio.h>


#ifndef G_PHOCS_TAUBOUNDS_H
#define G_PHOCS_TAUBOUNDS_H

#include "LocusDataLikelihood.h"

void initializeTauBounds();

void allocateTauBoundsMem();

void calculateTauBounds1();

void calculateLociTauBounds(int nodeId, int gen);

void calculateTauBounds2();

void printTauBoundsHeader(FILE *file);

void printTauBounds(int iteration, FILE *file);




//TODO - move to util
bool isLeafNode(LikelihoodNode *node);
LikelihoodNode *getNode(int nodeId, int gen);


#endif //G_PHOCS_TAUBOUNDS_H
