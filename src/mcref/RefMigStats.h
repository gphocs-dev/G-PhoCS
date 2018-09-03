#ifndef G_PHOCS_REFMIGSTATS_H
#define G_PHOCS_REFMIGSTATS_H

#include "../LocusDataLikelihood.h"

void calculateReferenceMigrationStats();

void initRefMigStats();

int calculateReferenceGenMigStats(int nodeId, double parentNodeAge, int gen);

double migBandIntersection(int mb, double bottomAge, double topAge);

// ============================

void allocateRefMigStatsMem();

void printRefMigStatsHeader(FILE *file);

void printReferenceMigrationStats(int iteration, FILE *file);

#endif //G_PHOCS_REFMIGSTATS_H
