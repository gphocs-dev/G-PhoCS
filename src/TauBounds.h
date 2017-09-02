#include <stdio.h>

#ifndef G_PHOCS_TAUBOUNDS_H
#define G_PHOCS_TAUBOUNDS_H

void initializeTauBounds();
void allocateTauBoundsMem();
void calculateTauBounds();
void printTauBoundsHeader(FILE *file);
void printTauBounds(int iteration, FILE *file);
#endif //G_PHOCS_TAUBOUNDS_H
