#include <stdio.h>


#ifndef G_PHOCS_TAUBOUNDS_H
#define G_PHOCS_TAUBOUNDS_H

#include "../LocusDataLikelihood.h"

extern double *tau_ubounds;
extern double *tau_lbounds;

void initializeBounds();

void propagateBoundsAcrossPopTree();

void allocateTauBoundsMem();

void calculateTauBounds();

void calculateUpperAndLowerBounds();

int calculateLocusTauUpperBounds(int nodeId, int gen);

double updateLowerBoundsOfDescendants(int pop);

void updateUpperBoundsOfDescendants(int pop, double bound);

void printTauBoundsHeader(FILE *file);

void printTauBounds(int iteration, FILE *file);

void assertBoundsAreMonotonousAscending();

void assertTausAreBetweenBounds();

void runTauBoundsAssertions();

#endif //G_PHOCS_TAUBOUNDS_H
