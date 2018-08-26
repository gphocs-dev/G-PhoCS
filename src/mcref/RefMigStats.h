#ifndef G_PHOCS_REFMIGSTATS_H
#define G_PHOCS_REFMIGSTATS_H

void allocateRefMigStatsMem();

void calculateReferenceMigrationStats();

void printRefMigStatsHeader(FILE *file);

void printReferenceMigrationStats(int iteration, FILE *file);

#endif //G_PHOCS_REFMIGSTATS_H
