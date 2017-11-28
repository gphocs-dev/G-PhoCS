#ifndef G_PHOCS_HYPOTHESISPRINTER_H
#define G_PHOCS_HYPOTHESISPRINTER_H

#include <stdio.h>


void printHypStatsHeader(FILE *hypStatsFile);

void printOnePopHeader(int pop, FILE *file);

void printOneMigBandHeader(int migband, FILE *file);

void printHypStats(int iteration, FILE *file);

void printOnePopStats(int pop, FILE *file);

void printOneMigBandStats(int migband, FILE *file);


#endif
