#include "McRefCommon.h"



double calculateCoalStats(double* elapsed_times, int* num_lineages, int size){
	int n;
	double t;
	double result = 0.0;
	for( int i = 0 ; i < size ; i++) {
		n = num_lineages[i];
		t = elapsed_times[i];
		result += n*(n-1)*t;
	}
	return result;
}
