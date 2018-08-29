#ifndef G_PHOCS_REFMIGSTATS_H
#define G_PHOCS_REFMIGSTATS_H

#include "../LocusDataLikelihood.h"

void calculateReferenceMigrationStats();

void initRefMigStats();

int calculateReferenceGenMigStats(int nodeId, int gen);

void addChildEdgesToMigStats(const LikelihoodNode *currentNode, int gen, int lcaPop);

double migBandIntersection(int mb, double fromAge, double toAge);

// ============================

void allocateRefMigStatsMem();

void printRefMigStatsHeader(FILE *file);

void printReferenceMigrationStats(int iteration, FILE *file);

#endif //G_PHOCS_REFMIGSTATS_H
