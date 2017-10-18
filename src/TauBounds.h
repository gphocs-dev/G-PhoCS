#include <stdio.h>


#ifndef G_PHOCS_TAUBOUNDS_H
#define G_PHOCS_TAUBOUNDS_H

#include "LocusDataLikelihood.h"

void initializeBounds();

void updateTauBoundsOfDescendants(int pop, double bound);

void propagateBoundsDownPopTree();

void allocateTauBoundsMem();

void calculateTauBounds();

int calculateLociTauBounds(int nodeId, int gen);

void printTauBoundsHeader(FILE *file);

void printTauBounds(int iteration, FILE *file);

void assertBoundsAreMonotonousAscending();

void assertTausAreSmallerThanBounds();

void runTauBoundsAssertions();

#endif //G_PHOCS_TAUBOUNDS_H
