#include <stdio.h>


#ifndef G_PHOCS_TAUBOUNDS_H
#define G_PHOCS_TAUBOUNDS_H

#include "LocusDataLikelihood.h"

void initializeBounds();

void propagateBoundsAcrossPopTree();

void allocateTauBoundsMem();

void calculateTauBounds();

void calculateTauUpperBounds1();

void calculateTauUpperBounds2();

void calculateTauLowerBounds();

int calculateLocusTauUpperBounds(int nodeId, int gen);

int migLcaPop(int nodeId, int gen, int defaultLcaPop);

double updateLowerBoundsOfDescendants(int pop);

void updateUpperBoundsOfDescendants(int pop, double bound);

void printTauBoundsHeader(FILE *file);

void printTauBounds(int iteration, FILE *file);

void assertBoundsAreMonotonousAscending();

void assertTausAreSmallerThanBounds();

void runTauBoundsAssertions();

void computeLcas();

#endif //G_PHOCS_TAUBOUNDS_H
