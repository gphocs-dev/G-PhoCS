/*
 * CombComb.h
 *
 *  Created on: Feb 2, 2017
 *      Author: ron
 */

#ifndef SRC_COMBSTATS_H_
#define SRC_COMBSTATS_H_




// clade_stats calculation functions
void computeCombStats();
void initCombStats();
void initSpecificCombStats(int clade);
void computeCombNumMigs();
void computeCombMigStats();
void computeCombNumCoals();
void computeCombNumCoals_rec(int pop);
void computeCombCoalStats();
void computeCombCoalStats_rec(int clade, int gen);
void fillupLeafCombStats(int clade, int gen);
void appendPopToComb(int clade, int gen, int startingPoint);
void fillupCombStats(int clade, int gen);
void addChildenIntoCombStats(int clade, int gen);
void mergeChildern(int clade, int gen);
void addChildrenCombStats(int clade, int gen);
void addCurrentPopIntoCombStats(int clade, int gen);
double getCoalStats(double* elapsed_times, int* num_lineages, int size);

int	isLeafPopulation(int pop);
int printLine();


#endif /* SRC_COMBSTATS_H_ */
