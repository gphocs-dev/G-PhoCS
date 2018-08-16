#include <stdio.h>


#ifndef G_PHOCS_TAUBOUNDS_H
#define G_PHOCS_TAUBOUNDS_H

#include "LocusDataLikelihood.h"

void initializeBounds();

void updateTauBoundsOfDescendants(int pop, double bound);

void propagateBoundsDownPopTree();

void allocateTauBoundsMem();

void calculateTauBounds();

void calculateTauUpperBounds1();
void calculateTauUpperBounds2();

void calculateTauLowerBounds();

void propagateLowerBoundUpwards(int pop, double bound);
void propagateUpperBoundDownwards(int pop, double bound);

int calculateLocusTauBounds(int nodeId, int gen);

void printTauBoundsHeader(FILE *file);

void printTauBounds(int iteration, FILE *file);

void assertBoundsAreMonotonousAscending();

void assertTausAreSmallerThanBounds();

void runTauBoundsAssertions();

void computeLcas();

#endif //G_PHOCS_TAUBOUNDS_H
