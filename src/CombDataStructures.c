
struct COMB_STATS{
	double coal_stats_total, mig_stats_total;
	int num_coals_total, num_migs_total;

	double* sorted_ages; 	// temporary array to hold sorted node ages of a specific genealogy.
	double* elapsed_times;	// temporary array to hold elapsed_time between adjacent events of a specific genealogy.
	int* num_lineages; 		// temporary array to hold num of lineages of a specific genealogy.
	int* event_types; 		// temporary array to hold event types

	int num_events; 		// temporary number of events in clade.
	// used as an "array length" variable for sortedAges, event_types & num_lineages arrays
	double debug_total_error;

} *comb_stats;
