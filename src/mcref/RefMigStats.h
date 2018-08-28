#ifndef G_PHOCS_REFMIGSTATS_H
#define G_PHOCS_REFMIGSTATS_H

void calculateReferenceMigrationStats();

void initRefMigStats();

int calculateReferenceGenMigStats(int nodeId, int gen);

// ============================

void allocateRefMigStatsMem();

void printRefMigStatsHeader(FILE *file);

void printReferenceMigrationStats(int iteration, FILE *file);

#endif //G_PHOCS_REFMIGSTATS_H
