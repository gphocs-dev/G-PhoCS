#ifndef SRC_MCREFCOMMON_H_
#define SRC_MCREFCOMMON_H_

#define TRUE 1
#define FALSE 0

double calculateCoalStats(double* elapsed_times, int* num_lineages, int size);

int	isLeaf(int pop);
int areChildrenLeaves(int pop);
int isAncestralTo(int father, int son);


const char* getEventTypeName(int eventType);

int getSon(int pop, int SON);
int getSourcePop(int mig);
int getTargetPop(int mig);
char* getPopName(int pop);
char* getMigName(int mig);
char* concat(const char *s1, const char *s2); // TODO - THIS SHOULD NOT BE USED  IN PRODUCTION SINCE THE MEMORY ISN'T RELEASED
#endif /* SRC_MCREFCOMMON_H_ */
