#ifndef G_PHOCS_REFMIGSTATS_H
#define G_PHOCS_REFMIGSTATS_H

void allocateRefMigStatsMem();

void printRefMigStatsHeader(FILE *file);

void calculateReferenceMigrationStats();

void initRefMigStats() ;

void printReferenceMigrationStats(int iteration, FILE *file);

#endif //G_PHOCS_REFMIGSTATS_H
