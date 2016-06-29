/** 
   \file patch.c 
   Old MCMCcoal code remaining in the GPhoCS implementation, global memory allocation, correctness checks.

 */


#define MAX_MIG_BANDS	40			// max migration bands in the population tree
#define MAX_MIGS		10			// max migration events per genealogy
#define NSPECIES		20			// max # of species
#define NS				200			// max # of sequences
#define OLDAGE			999			// upper bound on age (can be extended...)
#define MAX_EVENTS   (NS + 2*NSPECIES + 3*MAX_MIGS)
#define NUM_DELTA_STATS_INSTANCES 2

#define DEBUG_NODE_CHANGE_NOT
#define DEBUG_RUBBERBAND_NOT



/***************************************************************************************************************/
/******                                  GLOBAL DATA STRUCTURES                                           ******/
/***************************************************************************************************************/


double* genLogLikelihood;
double* genDeltaLogLikelihood;			// for updateTau
double averageMigTimes[MAX_MIG_BANDS];

/* node surrogates.
 */
int**	nodePops;			// a 2D array (numLoci X numNodes) for populations per genealogy node.
int**	nodeEvents;			// a 2D array (numLoci X numNodes) for event ids per genealogy node.

struct ADMIXTURE_STATUS {
	int 	numSampledLoci;				// number of loci for which to show admixture status
	int*	sampledLoci;				// array of sampled loci
	double** sampleLocusAdmixRate;		// array of admixture rates (up to the curent point of sampling) for each locus and each admixed sample
	int*	admixtureCounts;			// number of loci in which sample is in alternative population
	double* admixtureCoefficients;		// estimated coefficients for all samples
} admixture_status;

/* genetree_migs is a struct which contains information about
   migration events in a specific genealogy.
   Array of structs (of length numLoci is allocated in GetMem().
 */
struct GENETREE_MIGS {
	int num_migs;							// number of migration event in genetree
	int living_mignodes[MAX_MIGS];	// array of indices for mignodes for dynamic managing of migration nodes
	struct MIGNODE{
		int gtree_branch;						// id of gene tree node below the migration node
		int migration_band;					// id of migration band relevant to this node
		int target_pop, source_pop;		// source and target populations of migration event (backwards view)
		int target_event, source_event;	// ids of event corresponding to this migration node.
		double age;								// time of event
	} mignodes[MAX_MIGS];
}* genetree_migs;


/* event chain
   Each event corresponds to a time band within a population where no events
   (coalescence/migration) take place. An event is attributed with one of 5 types
   corresponding to the event taking place at the end of the interval. Events are
   sorted in a list according to chronology within a population.
   Actual array of events is allocated in getMem()
 */
struct EVENT_CHAIN{
	int total_events;					// total number of events pre-allocated to this chain
	int first_event[2*NSPECIES-1];		// pointers to first event for every population
	int last_event[2*NSPECIES-1];		// pointers to last event for every population
	int free_events;					// pointer to a chain of free events for use. Always have at least one free event
	struct EVENT{
		enum event_type {COAL, IN_MIG, OUT_MIG, MIG_BAND_START, MIG_BAND_END, SAMPLES_START, END_CHAIN, DUMMY} type;
		int node_id, next, prev;
		double elapsed_time;			// time from last event
		int num_lineages;				// number of lineages before the event
	}* events;
} *event_chains;

char* getEventTypeName(enum event_type type){ //TODO - is this duplicate??
	switch (type)
	{
		case COAL: return "COAL";
		case IN_MIG: return "IN_MIG";
		case OUT_MIG: return "OUT_MIG";
		case MIG_BAND_START: return "MIG_BAND_START";
		case MIG_BAND_END: return "MIG_BAND_END";
		case SAMPLES_START: return "SAMPLES_START";
		case END_CHAIN: return "END_CHAIN";
		case DUMMY: return "DUMMY";
	}
	return "ERROR - BAD event_type !";
}


/*	genetree stats
	holds relevant statistics for fast computation of probability of tree
	given model parameters (split times, pop sizes and migration rates). 
	Holds statistics for every population is species tree.
	-- NEED TO ADD DOCUMENTATION FOR THIS --  !!!
	genetree_stats[gen] holds relevant information for genetree of genealogy 'gen'. 
    array allocated in GetMem().
	genetree_stats_total holds sum of statistics for all loci. coal_stats here
    considers also all gen-specific heredity factors (but not thetas).
 */
struct GENETREE_STATS{
	double coal_stats[2*NSPECIES-1], mig_stats[MAX_MIG_BANDS];
	int num_coals[2*NSPECIES-1], num_migs[MAX_MIG_BANDS];
} *genetree_stats, *genetree_stats_total_partitioned, genetree_stats_total, genetree_stats_check, genetree_stats_total_check;						



/*	genetree stats flat
	holds relevant statistics for fast computation of probability of tree
	given a null model (single population).
    array allocated in GetMem().
 */
struct GENETREE_STATS_FLAT{
	double coal_stats_flat, mig_stats_flat;
	int num_coals_total, num_migs_total;
	double* sortedAgesArray;
} genetree_stats_flat;						


/*	clade_stats - holds statistics of all clades of tree.
 *  clade_stats[pop] holds statistics of the clade under (and including) pop.
 *  arrays allocated in GetMem().
 */
struct CLADE_STATS{
	double coal_stats_total, mig_stats_total;
	int num_coals_total, num_migs_total;

	double* sorted_ages; 	// temporary array to hold sorted node ages of a specific genealogy.
	double* elapsed_times;	// temporary array to hold elapsed_time between adjacent events of a specific genealogy.
	int* num_lineages; 		// temporary array to hold num of lineages of a specific genealogy.
	int* event_types; 		// temporary array to hold event types

	int num_events; 		// temporary number of events in clade.
	// used as an "array length" variable for sortedAges, event_types & num_lineages arrays

} *clade_stats;








/*	genetree node stats
	holds relevant statistics for coalescent distribution across populations
    array allocated in GetMem().
 */
struct GENETREE_NODE_STATS{
	/*** memory allocation ***/
	double*   doubleArray;
	double**  doublePtrArray;
	double*** doublePtrPtrArray;
	int*      intArray;
	int**     intPtrArray;

	/*** auxilliary matrices for computation ***/
	int**     lcaMatrix;            // matrix of LCAs for all pairs of leaves
	int*      leafArray;            // array of leaves for LCA computation
	int*      firstNodesInPop;      // array with the first coal node in each pop
	double*   firstNodesAges;       // array with the age of the first coal node in each pop
	double*   nodeAges;              // array with the ages of all internal nodes

	/*** main stat matrices ***/
	double*** probCoalMatrix;       // Prob[sample pair coalesce in pop]
	double*** probFirstCoalMatrix;  // Prob[sample pair is the first coalescence in pop]
	double*** coalTimeMatrix;       // mean coal time for coalescence of pair (given they coalesce in pop)
} genetree_node_stats;						


/*	genetree_stats_delta
	Holds difference in stats for populations affected by changes in a genetree.
	Not generally used for changes in theta, mig_rates.
	For easy and fast update if changes are accepted. 
	changed_events holds the ids of all events in the interval between original and 
	new position of event. All these intervals have a change in number of lineages.
 */
struct GENETREE_STATS_DELTA{
	int original_event;							// id of event describing original placing of node
	int updated_event;							// id of event describing updated (new) placing of node
	int num_lin_delta;							// the difference in lineage number for all events affected (typically +1 or -1)
	int num_changed_events;						// number of events affected by change
	int changed_events[MAX_EVENTS];			// an array of ids of events affected by change
	int num_pops_changed;						// number of population affected by change
	int pops_changed [2*NSPECIES-1];			// an array of populations affected by change
	int num_mig_bands_changed;					// number of migration bands affected by change
	int mig_bands_changed[MAX_MIG_BANDS];	// an array of migration bands affected by change
	double coal_stats_delta[2*NSPECIES-1];	// difference in coalescence statistics per population affected
	double mig_stats_delta[MAX_MIG_BANDS];	// difference in migration statistics per migration band affected
} genetree_stats_delta[NUM_DELTA_STATS_INSTANCES];

/*	mig_spr_stats
	holds statistics for the SPR sampling operation with migration.
	In use in UpdateGB_MigSPR and in traceLineage.
 */
struct MIG_SPR_STATS{
	int father_event_old, father_event_new;
	int father_pop_new;
	int target;
	int num_old_migs, num_new_migs;
	int old_migs[MAX_MIGS], new_migs_in[MAX_MIGS], new_migs_out[MAX_MIGS], new_migs_bands[MAX_MIGS];
	double new_migs_ages[MAX_MIGS];
	double genetree_delta_lnLd[2];
} mig_spr_stats;

/* rubberband_migs
   structure for holding data on migration events out of rubber-banded populations
   which are affected by rubber-band operation (per gen).
   num_moved_events - total number of affected events. (in/out migrations and start/end of migration bands).
   orig_events      - original copies of events.
   new_events       - new copies of events.
   pops				  - population in which each event resides.
   new_ages			  - age of each new event.

   rubberband_migs is an array of size numLoci allocated in getMem().
 */

struct RUBBERBAND_MIGS{

	int num_moved_events, orig_events[MAX_MIGS+MAX_MIG_BANDS], new_events[MAX_MIGS+MAX_MIG_BANDS], pops[MAX_MIGS+MAX_MIG_BANDS];
	double new_ages[MAX_MIGS+MAX_MIG_BANDS];

} *rubberband_migs;



/******************************************************************************************************/
/******                                FUNCTION DECLARATIONS                                     ******/
/******************************************************************************************************/



int findLastMig(int gen, int node_id, double time);
int findFirstMig(int gen, int node_id, double time);
int Coalescence1Pop (PopulationTree* popTree, GenericBinaryTree* tree, int gen, int pop, int* livingLineages);
int removeEvent(int gen, int event);

// auxiliary functions
int getLineagesAtInterval(int gen, int start_event, int pop, int exc_node, int* out_array);
int getEdgesForTimePop(int gen, double time, int pop, int exc_node, int* out_array);
int populationPostOrder(int pop, int* ordered_pops);

// event chain functions
int checkAll();
int checkGtreeStructure(int gen);
int synchronizeEvents(int gen);
int printEventChains(FILE* stream, int gen);
int createEvent(int gen, int pop, double age);
int createEventBefore(int gen, int pop, int event, double elapsed_time);
int constructEventChain(int gen);
void computeFlatStats();
int computeNodeStats();
int computeTotalStats();




double recalcStats(int gen, int pop);
int recalcStats_partitioned(int gen, int pop);
int computeGenetreeStats(int gen);
double gtreeLnLikelihood(int gen);
double computeDeltaLnLd(int gen, int instance);
int computeMigStatsDelta(int instance, double bottom_age, int bottom_pop, double top_age, int num_lins_delta);
int computeCoalStatsDelta(int instance, int gen, int bottom_event, int bottom_pop, int top_event, int num_lins_delta);
double considerEventMove(int gen, int instance, int event_id, int source_pop, double original_age, int target_pop, double new_age);
int acceptEventChainChanges(int gen, int instance);
int rejectEventChainChanges(int gen, int instance);
double rubberBandRipple(int gen, int do_or_redo);
double rubberBand(int gen, int pop, double static_point, double moving_point, double factor, unsigned short postORpre, int* out_num_events);
int replaceMigNodes(int gen, int node);
int traceLineage(int gen, int node, int reconnect);



/******************************************************************************************************/
/******                               FUNCTION IMPLEMENTATION                                    ******/
/******************************************************************************************************/


/***********************************************************************************
 *	GetMem
 * 	- allocates all auxiliary memory (not allocated elsewhere)
 *	- returns 0
 ***********************************************************************************/
int GetMem (void)		{
	int gen, count, i, maxNodes = 2*dataSetup.numSamples-1;
	int leaf1, leaf2;
	int pop;

	genLogLikelihood = (double*)malloc(2*dataSetup.numLoci*sizeof(double));
	if(genLogLikelihood == NULL) {
		fprintf(stderr, "\nError: Out Of Memory genLogLikelihood.\n");
		exit(-1);
	}
	genDeltaLogLikelihood = genLogLikelihood + dataSetup.numLoci;

	// get mem for auxiliary node arrays
	nodePops = (int**)malloc(2*dataSetup.numLoci*sizeof(int*));
	if(nodePops == NULL) {
		fprintf(stderr, "\nError: Out Of Memory nodePop array.\n");
		exit(-1);
	}
	nodeEvents = nodePops + dataSetup.numLoci;
	nodePops[0] = (int*)malloc(2*dataSetup.numLoci*maxNodes*sizeof(int));
	if(nodePops[0] == NULL) {
		fprintf(stderr, "\nError: Out Of Memory nodePop space.\n");
		exit(-1);
	}
	for(gen=0; gen<dataSetup.numLoci; gen++) {
		nodePops[gen] = nodePops[0] + gen*maxNodes;
		nodeEvents[gen] = nodePops[gen] + dataSetup.numLoci*maxNodes;
	}

	// get memory for genetree_migs and initialize
	// get memory for event_chains and stats and initialize
	genetree_migs=(struct GENETREE_MIGS*)malloc(dataSetup.numLoci*sizeof(struct GENETREE_MIGS));
	if(genetree_migs == NULL) {
		fprintf(stderr, "\nError: Out Of Memory genetree migs.\n");
		exit(-1);
	}
	genetree_stats=(struct GENETREE_STATS*)malloc(dataSetup.numLoci*sizeof(struct GENETREE_STATS));
	if(genetree_stats == NULL) {
		fprintf(stderr, "\nError: Out Of Memory genetree stats.\n");
		exit(-1);
	}
	genetree_stats_total_partitioned=(struct GENETREE_STATS*)malloc(dataSetup.numPopPartitions*sizeof(struct GENETREE_STATS));
	if(genetree_stats_total_partitioned == NULL) {
		fprintf(stderr, "\nError: Out Of Memory genetree stats total partitioned.\n");
		exit(-1);
	}
	rubberband_migs=(struct RUBBERBAND_MIGS*)malloc(dataSetup.numLoci*sizeof(struct RUBBERBAND_MIGS));
	if(rubberband_migs == NULL) {
		fprintf(stderr, "\nError: Out Of Memory rubberband migs.\n");
		exit(-1);
	}
	event_chains=(struct EVENT_CHAIN*)malloc(dataSetup.numLoci*sizeof(struct EVENT_CHAIN));
	if(event_chains == NULL) {
		fprintf(stderr, "\nError: Out Of Memory event chains.\n");
		exit(-1);
	}
	// compute max number of events for all loci
	// max number of events should cover all possible coalescences (+1 auxiliary),
	// migrations (X2), added migrations (X2) and population endings, and migration bands (start+end+changed event).
	int max_events = 2*dataSetup.numSamples+ 4*MAX_MIGS + 3*dataSetup.popTree->numMigBands + dataSetup.popTree->numPops + 10;
	for(gen=0,count=0; gen<dataSetup.numLoci; gen++) {
		event_chains[gen].total_events = max_events;
		count += event_chains[gen].total_events;
	}
	event_chains[0].events=(struct EVENT*)malloc(count*sizeof(struct EVENT));
	if(event_chains[0].events == NULL) {
		fprintf(stderr, "\nError: Out Of Memory event space.\n");
		exit(-1);
	}
	for(gen=0; gen<dataSetup.numLoci; gen++) {
		if(gen > 0) {
			event_chains[gen].events = event_chains[gen-1].events + event_chains[gen-1].total_events;
		}
		genetree_migs[gen].num_migs = 0;
		//initialize mignodes
		for(i=0; i<MAX_MIGS; i++) {
			genetree_migs[gen].mignodes[i].age=0;
			genetree_migs[gen].mignodes[i].target_pop=-1;
			genetree_migs[gen].mignodes[i].source_pop=-1;
			genetree_migs[gen].mignodes[i].gtree_branch=-1;
			genetree_migs[gen].mignodes[i].source_event=-1;
			genetree_migs[gen].mignodes[i].target_event=-1;
			genetree_migs[gen].mignodes[i].migration_band=-1;
		}
	}


	genetree_stats_flat.sortedAgesArray = (double*)malloc(maxNodes*sizeof(double));
	if(genetree_stats_flat.sortedAgesArray == NULL) {
		fprintf(stderr, "\nError: Out Of Memory event chains.\n");
		exit(-1);
	}
	clade_stats=(struct CLADE_STATS*)malloc(dataSetup.popTree->numPops*sizeof(struct CLADE_STATS));
	for (int clade = 0; clade < dataSetup.popTree->numPops; clade++) {
		clade_stats[clade].sorted_ages = (double*)malloc(max_events*sizeof(double));
		clade_stats[clade].elapsed_times = (double*)malloc(max_events*sizeof(double));
		clade_stats[clade].num_lineages = (int*)malloc(max_events*sizeof(int));
		clade_stats[clade].event_types = (int*)malloc(max_events*sizeof(int));
		clade_stats[clade].coal_stats_total = 0.0;

		clade_stats[clade].mig_stats_total = 0.0;


		if((clade_stats == NULL) ||
				(!clade_stats[clade].sorted_ages) || (!clade_stats[clade].elapsed_times) ||
				(!clade_stats[clade].num_lineages)|| (!clade_stats[clade].event_types)) {
			fprintf(stderr, "\nError: Out Of Memory clade_stats\n");
			exit(-1);
		}
	}

	// initialize genetree_stats_delta
	for(i=0;i < NUM_DELTA_STATS_INSTANCES; i++) {
		genetree_stats_delta[i].num_pops_changed = 0;
		genetree_stats_delta[i].num_mig_bands_changed = 0;
		genetree_stats_delta[i].num_changed_events = 0;
		genetree_stats_delta[i].num_lin_delta = 0;
		genetree_stats_delta[i].original_event = -1;
		genetree_stats_delta[i].updated_event = -1;
	}

	// initialize genetree_node_stats

	count = 3*dataSetup.numSamples*dataSetup.numSamples*dataSetup.popTree->numPops + maxNodes + dataSetup.popTree->numPops;
	genetree_node_stats.doubleArray = (double*)malloc(count*sizeof(double));
	if(genetree_node_stats.doubleArray == NULL) {
		fprintf(stderr, "\nError: Out Of Memory genetree node stats - doubles.\n");
		exit(-1);
	}

	count = 3*dataSetup.numSamples*dataSetup.numSamples;
	genetree_node_stats.doublePtrArray = (double**)malloc(count*sizeof(double*));
	if(genetree_node_stats.doublePtrArray == NULL) {
		fprintf(stderr, "\nError: Out Of Memory genetree node stats - double pointers.\n");
		exit(-1);
	}

	count = 3*dataSetup.numSamples;
	genetree_node_stats.doublePtrPtrArray = (double***)malloc(count*sizeof(double*));
	if(genetree_node_stats.doublePtrPtrArray == NULL) {
		fprintf(stderr, "\nError: Out Of Memory genetree node stats - double pointers^2.\n");
		exit(-1);
	}

	count = dataSetup.numSamples*dataSetup.numSamples + dataSetup.numSamples + dataSetup.popTree->numPops;
	genetree_node_stats.intArray = (int*)malloc(count*sizeof(int));
	if(genetree_node_stats.intArray == NULL) {
		fprintf(stderr, "\nError: Out Of Memory genetree node stats - ints.\n");
		exit(-1);
	}

	count = dataSetup.numSamples;
	genetree_node_stats.intPtrArray = (int**)malloc(count*sizeof(int*));
	if(genetree_node_stats.intPtrArray == NULL) {
		fprintf(stderr, "\nError: Out Of Memory genetree node stats - int pointers.\n");
		exit(-1);
	}

	genetree_node_stats.lcaMatrix           =  genetree_node_stats.intPtrArray;
	genetree_node_stats.probCoalMatrix      =  genetree_node_stats.doublePtrPtrArray;
	genetree_node_stats.probFirstCoalMatrix =  genetree_node_stats.doublePtrPtrArray +    dataSetup.numSamples;
	genetree_node_stats.coalTimeMatrix      =  genetree_node_stats.doublePtrPtrArray +  2*dataSetup.numSamples;

	for(leaf1=0; leaf1<dataSetup.numSamples; leaf1++) {
		genetree_node_stats.lcaMatrix[leaf1]           = genetree_node_stats.intArray + leaf1*dataSetup.numSamples;
		genetree_node_stats.probCoalMatrix[leaf1]      = genetree_node_stats.doublePtrArray + 3*leaf1*dataSetup.numSamples;
		genetree_node_stats.probFirstCoalMatrix[leaf1] = genetree_node_stats.probCoalMatrix[leaf1] + dataSetup.numSamples;
		genetree_node_stats.coalTimeMatrix[leaf1]      = genetree_node_stats.probCoalMatrix[leaf1] + 2*dataSetup.numSamples;
		for(leaf2=0; leaf2<dataSetup.numSamples; leaf2++) {
			genetree_node_stats.probCoalMatrix[leaf1][leaf2]      = genetree_node_stats.doubleArray + 3*dataSetup.popTree->numPops*(leaf1*dataSetup.numSamples+leaf2);
			genetree_node_stats.probFirstCoalMatrix[leaf1][leaf2] = genetree_node_stats.probCoalMatrix[leaf1][leaf2] + dataSetup.popTree->numPops;
			genetree_node_stats.coalTimeMatrix[leaf1][leaf2]      =  genetree_node_stats.probCoalMatrix[leaf1][leaf2] + 2*dataSetup.popTree->numPops;

		}
	}
	genetree_node_stats.leafArray       = genetree_node_stats.intArray + dataSetup.numSamples*dataSetup.numSamples;
	genetree_node_stats.firstNodesInPop = genetree_node_stats.leafArray + dataSetup.numSamples;
	genetree_node_stats.firstNodesAges  = genetree_node_stats.doubleArray + 3*dataSetup.popTree->numPops*dataSetup.numSamples*dataSetup.numSamples;
	genetree_node_stats.nodeAges        = genetree_node_stats.firstNodesAges + dataSetup.popTree->numPops;

	return 0;
}
/** end of GetMem **/



/***********************************************************************************
 *	FreeMem
 * 	- Frees all auxiliary memory allocated in GetMem
 *	- returns 0
 ***********************************************************************************/
int FreeMem (void)		{
	free(genLogLikelihood);
	free(nodePops[0]);
	free(nodePops);
	free(genetree_migs);
	free(genetree_stats);
	free(genetree_stats_total_partitioned);
	free(rubberband_migs);
	free(event_chains[0].events);
	free(event_chains);
	free(admixed_samples.samples);
	free(admixed_samples.index);
	free(admixed_samples.popPairs[0]);
	free(admixed_samples.popPairs);
	free(admixture_status.sampledLoci);
	free(admixture_status.admixtureCounts);
	free(admixture_status.sampleLocusAdmixRate[0]);
	free(admixture_status.sampleLocusAdmixRate);
	free(admixture_status.admixtureCoefficients);
	free(genetree_stats_flat.sortedAgesArray);
	free(genetree_node_stats.doubleArray);
	free(genetree_node_stats.doublePtrArray);
	free(genetree_node_stats.doublePtrPtrArray);
	free(genetree_node_stats.intArray);
	free(genetree_node_stats.intPtrArray);

	return 0;
}
/** end of FreeMem **/



/***********************************************************************************
 *	GetRandomGtree
 * 	- generates genealogy according to model parameters (without migration)
 *	- uses subroutine Coalescence1Pop and global int nodeNumber_global
 *	- returns 0
 ***********************************************************************************/
int nextAvailableNodeId;


int GetRandomGtree(GenericBinaryTree* tree, int gen)		{
	int  pop;
	int* livingLineages = (int*) malloc(dataSetup.numSamples*sizeof(int));

	if(livingLineages == NULL) {
		fprintf(stderr, "\nError: Out Of Memory living lineages at gen %d.\n",gen);
		exit(-1);
	}
	// convert numSamplesPerPop to cumulative form
	for(pop=1; pop<dataSetup.popTree->numCurPops; pop++) {
		dataSetup.numSamplesPerPop[pop] += dataSetup.numSamplesPerPop[pop-1];
	}
	nextAvailableNodeId = dataSetup.numSamples;
	Coalescence1Pop(dataSetup.popTree, tree, gen, dataSetup.popTree->rootPop, livingLineages);

	tree->rootId = nextAvailableNodeId - 1;

	// convert numSamplesPerPop back
	for(pop=dataSetup.popTree->numCurPops-1; pop>0; pop--) {
		dataSetup.numSamplesPerPop[pop] -= dataSetup.numSamplesPerPop[pop-1];
	}
	free(livingLineages);

	return 0;
}
/** end of GetRandomGtree **/


/***********************************************************************************
 *	Coalescence1Pop
 * 	- simulates coalescence along a single population using population parameters in
 *		given population tree (no migration).
 *	- recursive subroutine used only by GetRandomGtree().
 *	- returns number of lineages coming out of population
 ***********************************************************************************/
/* MARK: FIX INITIAL TREE BUILDING PROCEDURE - USE POP AGE FOR NODE AGE
         DONE. NEED TO CHECK !!
 */
int Coalescence1Pop (PopulationTree* popTree, GenericBinaryTree* tree, int gen, int pop, int* livingLineages)		{

	int numLineages, node1, node2, choice;
	double t, T;

	// if current population initialize number of incoming lineages,
	// otherwise, run on both subtrees first
	if(pop < popTree->numCurPops) {
		if(pop>0) {
			node1 = dataSetup.numSamplesPerPop[pop-1];
		} else {
			node1 = 0;
		}
		numLineages = dataSetup.numSamplesPerPop[pop] - node1;
		for(node2=0; node2<numLineages; node2++) {
			livingLineages[node2] = node1+node2;
			nodePops[gen][node1+node2] = pop;
			nodeEvents[gen][node1+node2] = -1;
			tree->leftSon[node1+node2] = -1;
			tree->rightSon[node1+node2] = -1;
			tree->father[node1+node2] = -1;
			tree->label1[node1+node2] = popTree->pops[pop]->sampleAge;
		}

	} else {
		numLineages  = Coalescence1Pop(popTree, tree, gen, popTree->pops[pop]->sons[0]->id, livingLineages);
		numLineages += Coalescence1Pop(popTree, tree, gen, popTree->pops[pop]->sons[1]->id, livingLineages+numLineages);
	}

	T = popTree->pops[pop]->age;
	if(pop<popTree->numCurPops) {
		T = popTree->pops[pop]->sampleAge;
	}
	for( ; numLineages>1; numLineages--, nextAvailableNodeId++) {
		t = rndexp(popTree->pops[pop]->theta/(numLineages*(numLineages-1.)));
		T += t;
		if(pop != popTree->rootPop && T > popTree->pops[pop]->father->age)
			break;

		// choose first living node and remove from list
		choice = (int)(numLineages*rndu());
		node1 = livingLineages[choice];
		livingLineages[choice] = livingLineages[numLineages-1];

		// choose second living node and replace with new one
		choice = (int)((numLineages-1)*rndu());
		node2 = livingLineages[choice];
		livingLineages[choice] = nextAvailableNodeId;

		tree->rightSon[nextAvailableNodeId] = node1;
		tree->leftSon[nextAvailableNodeId]  = node2;
		tree->father[nextAvailableNodeId] = -1;
		tree->label1[nextAvailableNodeId]  = T;
		tree->father[node1] = nextAvailableNodeId;
		tree->father[node2] = nextAvailableNodeId;
		nodePops[gen][nextAvailableNodeId] = pop;

	}

	return numLineages;
}
/** end of Coalescence1Pop **/




/******************************************************************************************************/
/******                                 MIGRATION FUNCTIONS                                      ******/
/******************************************************************************************************/



/* returns the id of the last migration event on the edge above given node
   and before specified time. If specified time is negative, then ignores it.
   Returns -1 if no relevant migration events exist.
 */
int findLastMig(int gen, int node_id, double age) {
	int i,mig,last_mig = -1;
	for(i=0; i<genetree_migs[gen].num_migs; i++) {
		mig = genetree_migs[gen].living_mignodes[i];
		if(genetree_migs[gen].mignodes[mig].gtree_branch != node_id)
			continue;
		if((age < 0 || genetree_migs[gen].mignodes[mig].age < age) &&
				(last_mig < 0 || genetree_migs[gen].mignodes[mig].age > genetree_migs[gen].mignodes[last_mig].age)) {
			last_mig = mig;
		}
	}

	return last_mig;
}
/** findLastMig **/


/* returns the id of the first migration event on the edge above given node
   and after specified time. If specified time is negative, then returns first 
   event on specified edge.
   Returns -1 if no relevant migration events exist.
 */
int findFirstMig(int gen, int node_id, double age) {
	int i,mig,first_mig = -1;
	for(i=0; i<genetree_migs[gen].num_migs; i++) {
		mig = genetree_migs[gen].living_mignodes[i];
		if(genetree_migs[gen].mignodes[mig].gtree_branch != node_id)
			continue;
		if(genetree_migs[gen].mignodes[mig].age > age &&
				(first_mig < 0 || genetree_migs[gen].mignodes[mig].age < genetree_migs[gen].mignodes[first_mig].age)) {
			first_mig = mig;
		}
	}

	return first_mig;
}
/** findFirstMig **/



/*	getLineagesAtInterval - UNUSED
	returns the set of nodes for which the lineage above them lives in a given interval.
	event - indicates the interval
	pop   - the population in which event lives
	exc_node - a node to exclude in result
	out_array - array in which to put all node ids

	returns 0.
 */
int getLineagesAtInterval_UNUSED(int gen, int start_event, int pop, int exc_node, int* out_array) {

	int i,id, node, event;

	int num_targets = event_chains[gen].events[start_event].num_lineages;
	int exc_nodes[2*NS-1];

	for(i=0; i<2*dataSetup.numSamples-1; i++) {
		exc_nodes[i] = 0;
	}
	exc_nodes[exc_node] = 1;

	event = start_event;

	while(num_targets > 0) {
		id = event_chains[gen].events[event].node_id;
		switch(event_chains[gen].events[event].type) {
			case(COAL):
																																node = getNodeSon(dataState.lociData[gen], id, 0);
			//node = nodes[id].sons[0];
			if(!exc_nodes[node]) {
				out_array[--num_targets] = node;
			}
			node = getNodeSon(dataState.lociData[gen], id, 1);
			//node = nodes[id].sons[1];
			if(!exc_nodes[node]) {
				out_array[--num_targets] = node;
			}
			exc_nodes[id] = 1;
			break;
			case(IN_MIG):
																																node = genetree_migs[gen].mignodes[id].gtree_branch;
			if(!exc_nodes[node]) {
				out_array[--num_targets] = node;
			}
			break;
			case(OUT_MIG):
    																															exc_nodes[ genetree_migs[gen].mignodes[id].gtree_branch ] = 1;
			break;
			case(MIG_BAND_START):
			case(SAMPLES_START):
			// DIDN'T ADD OPERATION FOR THIS CASE BECAUSE FUNCTION WENT OUT OF USE
			case(MIG_BAND_END):
			case(END_CHAIN):
			case(DUMMY):
			break;
		}

		event = event_chains[gen].events[event].next;

		if(event < 0) {
			if(dataSetup.popTree->pops[pop]->father == NULL) {
				if(debug) {
					fprintf(stderr, "\nError: getLineagesAtInterval reached top event and could not find %d lineages for event %d (excluding %d).\n",
							event_chains[gen].events[event].num_lineages, start_event, exc_node);
					fprintf(stderr, "Found the following lineages: ");
					for(i=num_targets; i<event_chains[gen].events[event].num_lineages; i++) {
						fprintf(stderr, "%d ",out_array[i]);
					}
					fprintf(stderr, ".\n");

				} else {
					fprintf(stderr, "Fatal Error 0001.\n");
				}
				printGenealogyAndExit(gen,-1);

			}
			pop = dataSetup.popTree->pops[pop]->father->id;
			event = event_chains[gen].first_event[pop];
		}
	}

	if(num_targets < 0) {
		if(debug) {
			fprintf(stderr, "\nError: getLineagesAtInterval found more than %d lineages for event %d (excluding %d).\n",
					event_chains[gen].events[event].num_lineages, start_event, exc_node);
			fprintf(stderr, "Found the following lineages: ");
			for(i=num_targets; i<event_chains[gen].events[event].num_lineages; i++) {
				fprintf(stderr, "%d ",out_array[i]);
			}
			fprintf(stderr, ".\n");
		} else {
			fprintf(stderr, "Fatal Error 0002.\n");
		}

		printGenealogyAndExit(gen,-1);
	}

	return 0;

}
/** getLineagesAtInterval_UNUSED **/



/* returns a list of edges living at time 'time' in a population 'pop'.
   Excludes edge above exc_node.
   The list is constructed in given array and the number of edges is the returned value.
   the life span of an edge includes its child node and does not include its parent node.
 */

int getEdgesForTimePop(int gen, double time, int pop, int exc_node, int* out_array) {
	int node, mig, pop1, num_edges = 0;
	int node_father;

	// sanity check
	if(dataSetup.popTree->pops[pop]->age > time+0.0000001) {
		if(debug) {
			fprintf(stderr, "\nError: getEdgesForTimePop: population %d of age %g greater than requested age %g.\n",pop, dataSetup.popTree->pops[pop]->age, time);
		} else {
			fprintf(stderr, "Fatal Error 0003.\n");
		}
		return 0;
	}



	for(node=0; node<2*dataSetup.numSamples-1; node++) {
		node_father = getNodeFather(dataState.lociData[gen], node);
		// check time span of edge
		if(	node == exc_node ||
				getNodeAge(dataState.lociData[gen],node) > time ||
				(node_father >= 0 && getNodeAge(dataState.lociData[gen], node_father) <= time)) {
			//		if(node == exc_node || nodes[node].age > time || (node != tree.root && nodes[nodes[node].father].age <= time)) {
			continue;
		}

		// if pop in question is root and edge lives at time in question, then it must live in root.
		if(pop == dataSetup.popTree->rootPop) {
			out_array[num_edges++] = node;
			continue;
		}

		// find the population below pop such that edge has no migration events above that point
		mig = findLastMig(gen, node, time);
		pop1 = (mig >= 0) ? (genetree_migs[gen].mignodes[mig].source_pop) : (nodePops[gen][node]);

		// if pop is ancestral to pop1, then edge lives in pop at time in question.
		if(dataSetup.popTree->pops[pop]->isAncestralTo[pop1]) {
			out_array[num_edges++] = node;
		}

	}

	return num_edges;
}
/** getEdgesForTimePop **/




/*************************************  EVENT CHAIN PROCEDURES   ****************************************/

/*	rubberBand
	Computes the effect of a rubber band expansion/contraction of a specified region 
	within a population in  agiven gen. 
	static_point - the age of the hold-point of the rubber band
	moving_point - the original age of the moving point of the rubber band
    (this distinction is necessary to compute new node's age).
	factor - the factor of expansion/contraction
	postORpre == 1 if called after acceptance (and changes should be actually implemented
    Note that if postORpre == 0, then event chains are not modified, but the rubberband influence
    is computed on them

	change is to be implemented within event chains and stat structures
    0 if only delta log-likelihood is to be computed
	out_mum_changes -  output - number of migration and coalescence events affected by change

	function returns delta in log likelihood of modification
 */

double rubberBand(int gen, int pop, double static_point, double moving_point, double factor, unsigned short postORpre, int* out_num_events) {

	int i, event, mig_band, node_id, living_mig_bands[MAX_MIG_BANDS], num_mig_bands;
	int num_lins, count_events, flag;
	double age, delta_time, heredity_factor = 1;
	double mig_rate, mig_stats_delta, coal_stats_delta, lnLd_delta;

	double 	factor_minus_one = factor - 1.0;

	double	start_time = min2(static_point, moving_point),
			end_time   = max2(static_point, moving_point);

	double age1;

	// special consideration for root.
	// CONSIDER CHANGING THIS !!!
	if(pop == dataSetup.popTree->rootPop) {
		start_time = moving_point;
		end_time = OLDAGE;
	}

	//initialization
	num_mig_bands = 0;
	count_events = 0;
	// coal_stats_delta is accumulative, whereas mig_stats_delta is not
	coal_stats_delta = 0.0;
	lnLd_delta = 0.0;
	mig_rate = 0.0;
	event = event_chains[gen].first_event[pop];
	age = dataSetup.popTree->pops[pop]->age;
	flag = (age >= start_time);		// flag = 1 when statistics are to be changed

#ifdef DEBUG_RUBBERBAND
	printf("\nG(%d) P(%d) || ",gen+1,pop);
#endif

	while(age<end_time) {
		delta_time = min2(event_chains[gen].events[event].elapsed_time, end_time - age);
		// already advance age (before further changes are made)
		age += delta_time;


		// if first time exceeding start_time, raise flag
		if(!flag && age > start_time) {
			flag = 1;
			delta_time = age-start_time;
		}

		// update affected statistics for interval
		if(flag) {
			// delta_time signifies the difference between new and old elapsed times of the interval
			delta_time *= factor_minus_one;
			num_lins = event_chains[gen].events[event].num_lineages;

			mig_stats_delta   = delta_time * num_lins;
			//coal_stats_delta is accumulative
			coal_stats_delta += mig_stats_delta * (num_lins-1) ;

			//take coal_stats_delta into account after loop is over
			lnLd_delta -= mig_stats_delta * mig_rate;

#ifdef DEBUG_RUBBERBAND
			printf("E%d a%3g t%d n%d t%3g-%3g | ",
					event, age, event_chains[gen].events[event].type,event_chains[gen].events[event].node_id,
					event_chains[gen].events[event].elapsed_time, event_chains[gen].events[event].elapsed_time+delta_time);
			/*
		printf("E%d n%d l%d t%3g-%3g l%3f-%3g | ",
             event, event_chains[gen].events[event].node_id, num_lins,event_chains[gen].events[event].elapsed_time,
             event_chains[gen].events[event].elapsed_time+delta_time,
             -event_chains[gen].events[event].elapsed_time*num_lins*(num_lins-1)/(dataSetup.popTree->pops[pop]->theta * heredity_factor),
             -(event_chains[gen].events[event].elapsed_time*num_lins + mig_stats_delta)*(num_lins-1)/ (dataSetup.popTree->pops[pop]->theta * heredity_factor));
			 */
#endif			
			// if required to implement changes, change elapsed time
			// and update migration stats (coal stats are updated after the loop)
			if(postORpre) {
				event_chains[gen].events[event].elapsed_time += delta_time;
				for(mig_band=0; mig_band<num_mig_bands; mig_band++) {
					genetree_stats[gen].mig_stats[living_mig_bands[mig_band]] += mig_stats_delta;
					genetree_stats_total.mig_stats[living_mig_bands[mig_band]] += mig_stats_delta;
				}
			}
		}// end if(flag)

		// finish updates, and don't modify node at this point
		if(age >= end_time && event_chains[gen].events[event].type != SAMPLES_START) {
			break;
		}

		// actions according to event type
		node_id = event_chains[gen].events[event].node_id;
		switch(event_chains[gen].events[event].type) {
			// modify event counts and rescale node ages, if necessary
			case(COAL):
					if(flag) {
						count_events++;
						// if done after acceptance, no need to do anything
						if(!postORpre) {
							age1 = getNodeAge(dataState.lociData[gen], node_id);
							age1 += (age1 - static_point)*factor_minus_one;
							// fprintf(stderr,"node %d, gen %d, age %g-->%g, static point %g \n", node_id, gen, getNodeAge(dataState.lociData[gen], node_id), age1,static_point);
#ifdef DEBUG_RUBBERBAND
							printf("age %g-->%g | ", getNodeAge(dataState.lociData[gen], node_id), age1);
#endif
							adjustGenNodeAge(dataState.lociData[gen], node_id, age1);
						}
					}// end of if(flag)
			break;
			case(SAMPLES_START):
					if(flag && dataSetup.popTree->pops[pop]->sampleAge > 0) {
						// if done after acceptance or if update is below, no need to do anything
						if(static_point < moving_point && !postORpre) {
							age1 = dataSetup.popTree->pops[pop]->sampleAge;
							age1 += (age1 - static_point)*factor_minus_one;
							//printf("-->sample age %g-->%g \n", dataSetup.popTree->pops[pop]->sampleAge, age1);
#ifdef DEBUG_RUBBERBAND
							printf("sample age %g-->%g | ", dataSetup.popTree->pops[pop]->sampleAge, age1);
#endif
							for(i=0; i<dataSetup.numSamples; i++) {
								if(nodePops[gen][i] == pop) {
									//printf("-->adjusting age of node %d in gen %4d to %g.\n",i,gen, age1);
									adjustGenNodeAge(dataState.lociData[gen], i, age1);
								}
							}
						}
					}
			break;
			// case(OUT_MIG):
			case(IN_MIG):
					if(flag) {
						// change ages of mignodes for in-migs
						if(postORpre) {
							// fprintf(stderr,"switching mig node %d in gen %d from age %g to age ",node_id,gen,genetree_migs[gen].mignodes[node_id].age);
							genetree_migs[gen].mignodes[node_id].age += (genetree_migs[gen].mignodes[node_id].age - static_point)*factor_minus_one;
							// fprintf(stderr,"%g, static point %g\n",genetree_migs[gen].mignodes[node_id].age, static_point);
						}
					}
			break;
			// modify living migration bands and accumulative migration rate
			case(MIG_BAND_START):
					mig_rate += dataSetup.popTree->migBands[node_id].migRate;
			living_mig_bands[num_mig_bands++] = node_id;
			break;
			case(MIG_BAND_END):
					mig_rate -= dataSetup.popTree->migBands[node_id].migRate;
			for(i=0; i<num_mig_bands; i++) {
				if(node_id == living_mig_bands[i])
					break;
			}
			if(i == num_mig_bands) {
				if(debug) {
					fprintf(stderr, "\nError: collectStats: migration band %d ended without starting in gen %d.\n", node_id, gen);
				} else {
					fprintf(stderr, "Fatal Error 0004.\n");
				}
				printGenealogyAndExit(gen,-1);
			}
			living_mig_bands[i] = living_mig_bands[--num_mig_bands];
			break;
			case(END_CHAIN):
#ifdef DEBUG_RUBBERBAND
													if(age < end_time) printf("\nrubber band for pop %d, gen %d ended at end-chain (%g time to go).\n",pop,gen, end_time- age);
#endif
			age = end_time;
			break;
			default:
				break;
		}// end of switch

		event = event_chains[gen].events[event].next;
	}// end of while

	// after loop is done, consider coalescence stats
	if(postORpre) {
		genetree_stats[gen].coal_stats[pop] += coal_stats_delta;
		genetree_stats_total.coal_stats[pop]  += coal_stats_delta/heredity_factor;
	}

#ifdef DEBUG_RUBBERBAND
	if(lnLd_delta != 0.0) {
		printf("After rubber band, lnLd considering only migration is %g.\n",lnLd_delta);
	}
#endif

	lnLd_delta -= coal_stats_delta / (dataSetup.popTree->pops[pop]->theta * heredity_factor);

#ifdef DEBUG_RUBBERBAND
	printf("\n");
#endif
	*out_num_events += count_events;

	return lnLd_delta;
}




/* rubberBandRipple
   Considers changes made by rubber band operation on "distant" populations.
   Uses list of migration nodes pre-collected in rubberband_migs strucure.
   do_or_redo flag is 0 if calculation done before accepting/rejecting,
   and 1 if needs to be redone (after rejecting).
   Returns delta in log-likelihood.

 */

double rubberBandRipple(int gen, int do_or_redo) {
	int i, pop, new_event, orig_event, affected_pops[2*NSPECIES-1];
	double delta_lnLd = 0.0;

	if(rubberband_migs[gen].num_moved_events == 0) return 0.0;

	//	printf("\nRipple");

	for(pop=0; pop<dataSetup.popTree->numPops; pop++) {
		affected_pops[pop] = 0;
	}

	for(i=0; i<rubberband_migs[gen].num_moved_events; i++) {
		pop = rubberband_migs[gen].pops[i];
		orig_event = rubberband_migs[gen].orig_events[i];
		affected_pops[pop] = 1;
		if(do_or_redo) {
			new_event = rubberband_migs[gen].new_events[i] = createEvent(gen, pop, rubberband_migs[gen].new_ages[i]);
			if( new_event < 0 ) {
				if(debug) {
					fprintf(stderr, "Error: problem creating event in rubber band ripple.\n");
				} else {
					fprintf(stderr, "Fatal Error 0005.\n");
				}
				printGenealogyAndExit(gen,-1);
			}
			event_chains[gen].events[new_event].type = event_chains[gen].events[orig_event].type;
			event_chains[gen].events[new_event].node_id = event_chains[gen].events[orig_event].node_id;
			event_chains[gen].events[orig_event].type = DUMMY;
		} else {
			new_event = rubberband_migs[gen].new_events[i];
			event_chains[gen].events[orig_event].type = event_chains[gen].events[new_event].type;
			// if removing first event in pop, update incoming number of lineages
			if( event_chains[gen].first_event[pop] == new_event ) {
				event_chains[gen].events[ event_chains[gen].events[new_event].next ].num_lineages =
						event_chains[gen].events[new_event].num_lineages;
			}
			removeEvent(gen,new_event);
		}
	}

	for(pop=0; pop<dataSetup.popTree->numPops; pop++) {

		if(affected_pops[pop]) {
			delta_lnLd += recalcStats(gen,pop);
		}
	}

	if(!do_or_redo) {
		rubberband_migs[gen].num_moved_events = 0;
	}

	return delta_lnLd;

}



/*	traceLineage
	Traces a lineage up a genetree. 
	It reconnect == 0, then traces original edge, reduces number of lineages
	in events along the way by 1, and records the stat changes in genetree_stats_delta[0].
	If reconnect == 1, then re-coalesces lineage with rest of genetree according to 
	coalescence with migration process determined by model parameters.
	records stat changes in genetree_stats_delta[0].

	In both cases, sets the appropriate entries in mig_spr_stats.

	Returns -1 if too many migration events were sampled (otherwise 0).

 */
int traceLineage(int gen, int node, int reconnect) {

	int i, pop, event, node_id, mig_band, mig_source, proceed;
	int num_live_mig_bands, live_mig_bands[MAX_MIG_BANDS];
	double age=-1.0, t;
	// these variables are needed only for reconnecting or delta_lnLd computation
	int target, num_targets, targets[2*NS-1];
	double 	event_sample, rate, mig_rate, theta,
	heredity_factor = 1;

	/******** initialization start  **************/
	pop = nodePops[gen][node];

	if(node < dataSetup.numSamples) {
		event = event_chains[gen].first_event[pop];
		while(event_chains[gen].events[event].type != SAMPLES_START && event_chains[gen].events[event].type != END_CHAIN) {
			event = event_chains[gen].events[event].next;
		}
		if(event_chains[gen].events[event].type == END_CHAIN) {
			if(debug) {
				fprintf(stderr, "\nError: traceLineage: while tracing edge above node %d in gen %d.\n",node,gen);
				fprintf(stderr, "  Could not find start event in pop %d.\n",pop);
				fprintf(stderr, ".\n");
			} else {
				fprintf(stderr, "Fatal Error 0101.\n");
			}
		}
		event = event_chains[gen].events[event].next;
	} else {
		event = event_chains[gen].events[nodeEvents[gen][node]].next;
	}
	theta = dataSetup.popTree->pops[pop]->theta*heredity_factor;

	age = getNodeAge(dataState.lociData[gen], node);
	//age = nodes[node].age;
	mig_spr_stats.genetree_delta_lnLd[reconnect] = 0.0;
	if(!reconnect) {
		mig_spr_stats.num_old_migs = 0;
		if(node != getLocusRoot(dataState.lociData[gen])) {
			mig_spr_stats.father_event_old = nodeEvents[gen][ getNodeFather(dataState.lociData[gen], node) ];
		}
		// mig_spr_stats.father_event_old = nodes[nodes[node].father].event_id;
	} else {
		mig_spr_stats.num_new_migs = 0;
	}


	genetree_stats_delta[reconnect].num_changed_events = 0;
	// assume all migration bands and populations are affected
	genetree_stats_delta[reconnect].num_pops_changed = dataSetup.popTree->numPops;
	for(i=0; i<dataSetup.popTree->numPops; i++) {
		genetree_stats_delta[reconnect].pops_changed[i] = i;
		genetree_stats_delta[reconnect].coal_stats_delta[i] = 0.0;
	}
	genetree_stats_delta[reconnect].num_mig_bands_changed = dataSetup.popTree->numMigBands;
	for(i=0; i<dataSetup.popTree->numMigBands; i++) {
		genetree_stats_delta[reconnect].mig_bands_changed[i] = i;
		genetree_stats_delta[reconnect].mig_stats_delta[i] = 0.0;
	}

	mig_rate = 0.0;
	num_live_mig_bands = 0;
	for(mig_band = 0; mig_band < dataSetup.popTree->numMigBands; mig_band++) {
		if(dataSetup.popTree->migBands[mig_band].targetPop == pop &&
				dataSetup.popTree->migBands[mig_band].startTime < age &&
				dataSetup.popTree->migBands[mig_band].endTime   > age)
		{
			mig_rate += dataSetup.popTree->migBands[mig_band].migRate;
			live_mig_bands[num_live_mig_bands++] = mig_band;
		}
	}
	mig_source = -1;
	proceed = 1;
	/******** initialization end    **************/

	// loop will stop when reaching original father's event (if !reconnect)
	// or when creating new coalescent event for lineage (if reconnect)
	while(proceed) {


		// if last event in population, move to next pop
		if(event < 0) {
			if(dataSetup.popTree->pops[pop]->father == NULL) {
				// reject if trying to reconnect above OLDAGE
				if(reconnect == 1) {
					//          fprintf(stderr, "\nError: In traceLineage: while tracing edge above node %d in gen %d:",node,gen);
					//          fprintf(stderr, "  trying to reconnect above OLDAGE.\n");
					return -1;
				}
				// if detaching subtree, you can reach top only if detaching the root edge
				// this has to be considered in migration scenarios
				//				if(node != getLocusRoot(dataState.lociData[gen]) {
				if(debug) {
					fprintf(stderr, "\nError: traceLineage: while tracing edge above node %d in gen %d.\n",node,gen);
					fprintf(stderr, "  Reached top event at age %f (reconnect == %d).\n",age,reconnect);
					fprintf(stderr, "  Events visited:");
					for(i=0;i<genetree_stats_delta[reconnect].num_changed_events; i++) {
						fprintf(stderr, " %d",genetree_stats_delta[reconnect].changed_events[i]);
					}
					fprintf(stderr, ".\n");
				} else {
					fprintf(stderr, "Fatal Error 0006.\n");
				}
				printGenealogyAndExit(gen,-1);
				//				}
			}
			pop = dataSetup.popTree->pops[pop]->father->id;
			theta = dataSetup.popTree->pops[pop]->theta*heredity_factor;
			event = event_chains[gen].first_event[pop];

			if(fabs(mig_rate) > 0.00000001 || num_live_mig_bands > 0) {
				if(debug) {
					fprintf(stderr, "\nError: traceLineage: while tracing edge above node %d in gen %d.\n",node,gen);
					fprintf(stderr, "  %d living migration bands when moving to population %d (with total rate %f):\n",num_live_mig_bands,pop, mig_rate);
					for(i=0;i<num_live_mig_bands; i++) {
						fprintf(stderr, " %d (rate %f)",live_mig_bands[i], dataSetup.popTree->migBands[ live_mig_bands[i] ].migRate);
					}
					fprintf(stderr, ".\n");
				} else {
					fprintf(stderr, "Fatal Error 0007.\n");
				}
				printGenealogyAndExit(gen,-1);
			}
			mig_rate = 0.0;

			if(fabs(age/dataSetup.popTree->pops[pop]->age-1) > 0.01) {
				if(debug) {
					fprintf(stderr, "\nError traceLineage: while tracing edge above node %d in gen %d.\n",node,gen);
					fprintf(stderr, "  age recorded at start of pop %d is %g, while age of pop is %g (diff = %g).\n",pop, age, dataSetup.popTree->pops[pop]->age,age-dataSetup.popTree->pops[pop]->age);
					for(i=0;i<num_live_mig_bands; i++) {
						fprintf(stderr, " %d (rate %f)",live_mig_bands[i], dataSetup.popTree->migBands[ live_mig_bands[i] ].migRate);
					}
					fprintf(stderr, ".\n");
				} else {
					fprintf(stderr, "Fatal Error 0008.\n");
				}
				printGenealogyAndExit(gen,-1);
			}
			age = dataSetup.popTree->pops[pop]->age;
		}// end of moving to next pop

		//printf("\ntraceLineage(%d,%d,%d): event %d.",gen,node,reconnect,event);

		node_id = event_chains[gen].events[event].node_id;

		if(!reconnect) {
			// if(!reconnect) then reduce the number of lineages along pruned edge by 1
			event_chains[gen].events[event].num_lineages--;
			t = event_chains[gen].events[event].elapsed_time;
			age += t;
			proceed = (event != mig_spr_stats.father_event_old);

			// if event is an in-migration event of edge above node,
			// record it and follow it to source population
			if(event_chains[gen].events[event].type == IN_MIG){
				//				printf("\nIn traceLineage(%d,%d,%d) encountering in-mig event %d for mignode %d on branch %d.",
				//								gen,node,reconnect, event, node_id, genetree_migs[gen].mignodes[node_id].gtree_branch);
				if(genetree_migs[gen].mignodes[node_id].gtree_branch == node)	{
					//					printf("... following.");
					mig_band = genetree_migs[gen].mignodes[node_id].migration_band;
					mig_source = genetree_migs[gen].mignodes[node_id].source_event;
					mig_spr_stats.old_migs[mig_spr_stats.num_old_migs++] = node_id;
				}
			}
		} else {
			// if(reconnet) then sample new event within interval (coal or mig)
			// if coal - create new event, record point of coalescence and terminate process
			// if mig  - create 2 new events, and record source of migration for later processing
			rate = mig_rate + 2*event_chains[gen].events[event].num_lineages/theta;
			// rate can be zero, if no lineages and no incoming migration bands
			if(rate <= 0) {
				t = event_chains[gen].events[event].elapsed_time;
			} else {
				t=rndexp( 1/rate );
			}

			if(t >= event_chains[gen].events[event].elapsed_time) {
				// no event sampled in this interval
				t = event_chains[gen].events[event].elapsed_time;
				age += t;
			} else {
				// sample event in this interval
				age += t;
				event_sample = rate*rndu();

				if(event_sample < mig_rate) {
					// migration event - figure out where to migrate
					if(MAX_MIGS <= genetree_migs[gen].num_migs + mig_spr_stats.num_new_migs  - mig_spr_stats.num_old_migs) {
						misc_stats.not_enough_migs++;
						return -1;
					}
					if(event_sample == 0.0 && debug) {
						fprintf(stderr, "\ntraceLineage() event_sample=%g when rate=%g at pop=%d, age=%g.\n",event_sample,rate,pop,age);
						printPopulationTree(dataSetup.popTree, stderr, 1);
						printLocusGenTree(dataState.lociData[gen], stderr, nodePops[gen], nodeEvents[gen]);
						printEventChains(stderr, gen);
						fprintf(stderr, "\n\n********************************************************************\n\n");
					}
					for(i=0; event_sample >= 0 && i < num_live_mig_bands; i++) {
						event_sample -= dataSetup.popTree->migBands[ live_mig_bands[i] ].migRate;
					}
					if(event_sample >= 0.0) {	// at this point i = num_live_mig_bands
						fprintf(stderr, "\ntraceLineage() after iterating through all migration bands, event_sample=%g when rate=%g, mig_rate=%g at pop=%d, age=%g.\n",
								event_sample,rate,mig_rate,pop,age);
						if(num_live_mig_bands <= 0) {
							if(debug) {
								fprintf(stderr, "\nError: traceLineage: while re-coalescing edge above node %d in gen %d:\n",node,gen);
								fprintf(stderr, "  num_live_mig_bands=%d but migration rate is %g (event_sample = %g, pop=%d, age=%g).\n",
										num_live_mig_bands,mig_rate,event_sample,pop,age);
							} else {
								fprintf(stderr, "Fatal Error 0009.\n");
							}
							printGenealogyAndExit(gen,-1);
						}
						fprintf(stderr, "%d live mig bands:",num_live_mig_bands);
						for(i=0, event_sample=mig_rate; i < num_live_mig_bands; i++) {
							fprintf(stderr, " %d (rate=%g)", live_mig_bands[i], dataSetup.popTree->migBands[ live_mig_bands[i] ].migRate);
							event_sample -= dataSetup.popTree->migBands[ live_mig_bands[i] ].migRate;
						}
						fprintf(stderr, " mig_rate - all rates = %g.\n", event_sample);
						printPopulationTree(dataSetup.popTree, stderr, 1);
						printLocusGenTree(dataState.lociData[gen], stderr, nodePops[gen], nodeEvents[gen]);
						printEventChains(stderr, gen);
						fprintf(stderr, "\n\n********************************************************************\n\n");
					}

					if(i <= 0) {
						if(debug) {
							fprintf(stderr, "\nError: traceLineage: while re-coalescing edge above node %d in gen %d:\n",node,gen);
							fprintf(stderr, "  i=%d and event_sample=%g (migration rate = %g, pop=%d, age=%g).\n",
									i,event_sample,mig_rate,pop,age);
						} else {
							fprintf(stderr, "Fatal Error 0009A.\n");
						}
						printGenealogyAndExit(gen,-1);
					}
					// create in and out events and follow out event
					// do not mark the events yet (unless accepted)
					mig_spr_stats.new_migs_bands[mig_spr_stats.num_new_migs] = mig_band = live_mig_bands[i-1];
					if(dataSetup.popTree->migBands[mig_band].targetPop != pop) {
						if(debug) {
							fprintf(stderr, "\nError: traceLineage: while re-coalescing edge above node %d in gen %d:\n",node,gen);
							fprintf(stderr, "  selected mig_band=%d but target population %d does not fit current population %d (age=%g).\n",
									mig_band,dataSetup.popTree->migBands[mig_band].targetPop,pop,age);
							fprintf(stderr, "%d live mig bands:",num_live_mig_bands);
							for(i=0, event_sample=mig_rate; i < num_live_mig_bands; i++) {
								fprintf(stderr, " %d (rate=%g)", live_mig_bands[i], dataSetup.popTree->migBands[ live_mig_bands[i] ].migRate);
								event_sample -= dataSetup.popTree->migBands[ live_mig_bands[i] ].migRate;
							}
							fprintf(stderr, " mig_rate - all rates = %g.\n", event_sample);
						} else {
							fprintf(stderr, "Fatal Error 0009B.\n");
						}
						printGenealogyAndExit(gen,-1);
					}
					mig_spr_stats.new_migs_ages [mig_spr_stats.num_new_migs] = age;
					event = mig_spr_stats.new_migs_in[mig_spr_stats.num_new_migs]  = createEventBefore(gen,pop,event,t);
					// mark source event for migration
					mig_source = mig_spr_stats.new_migs_out[mig_spr_stats.num_new_migs] =
							createEvent(gen,dataSetup.popTree->migBands[mig_band].sourcePop,age);
					if( mig_source < 0 ) {
						if(debug) {
							fprintf(stderr, "Error: problem encountered when creating new migration event for locus %d, mig band %d pops %d-->%d, in trace lineage.\n",
									gen, mig_band,dataSetup.popTree->migBands[mig_band].sourcePop, dataSetup.popTree->migBands[mig_band].targetPop );
						} else {
							fprintf(stderr, "Fatal Error 0010.\n");
						}
						printGenealogyAndExit(gen,-1);
					}
					mig_spr_stats.num_new_migs++;
				} else {
					// coalescence event - figure out with whom to coalesce
					//getLineagesAtInterval(gen,event,pop,node,targets);  <-- UNUSED
					num_targets = getEdgesForTimePop(gen, (age-t)+event_chains[gen].events[event].elapsed_time/2, pop, node, targets);
					if(num_targets !=  event_chains[gen].events[event].num_lineages) {
						if(debug) {
							fprintf(stderr, "\nError: in traceLineage: %d targets and %d lineages (at gen %d, node %d, pop %d, , event %d, age %g)\n",
									num_targets, event_chains[gen].events[event].num_lineages, gen, node, pop, event,
									(age-t)+event_chains[gen].events[event].elapsed_time/2);
						} else {
							fprintf(stderr, "Fatal Error 0011.\n");
						}
						printGenealogyAndExit(gen,-1);
					}
					/*					printf("\nFound %d targets at gen %d for node %d, at event %d, pop %d, age %f:",
								num_targets, gen, node, event, pop, age);
                                for(i=0; i< num_targets; i++) {
                                printf(" %d", targets[i]);
                                }
                                printf("\n"); 							
					 */
					// this simulates choosing a lineage uniformly at random
					i = (int)((event_sample - mig_rate) * theta/2);
					target = targets[i];
					// attach lineage in genetree and create new event
					//printf("\nGrafting subtree below node %d to edge above node %d at time %g.",node, target, age);

					executeGenSPR(dataState.lociData[gen], node, target, age);
					mig_spr_stats.father_pop_new = pop;

					//GraftNode(node, target, age, pop);
					//printf("\nOriginal data log-likelihood was %g and after re-grafting it is %g.",data.lnpDi[gen], lnpD_gen(gen));
					mig_spr_stats.target = target;
					mig_spr_stats.father_event_new = event = createEventBefore(gen,pop,event,t);
					// signal to terminate
					proceed = 0;
				}// if event == mig/coal
			}// if event in interval
		}// if (reconnect)

		// record stat changes
		genetree_stats_delta[reconnect].coal_stats_delta[pop] += 2*event_chains[gen].events[event].num_lineages*t;
		for(i=0; i<num_live_mig_bands; i++) {
			genetree_stats_delta[reconnect].mig_stats_delta[live_mig_bands[i]] += t;
		}
		genetree_stats_delta[reconnect].changed_events[genetree_stats_delta[reconnect].num_changed_events++] = event;

		// add log-likelihood of no events during interval
		mig_spr_stats.genetree_delta_lnLd[reconnect] -= (mig_rate + 2*event_chains[gen].events[event].num_lineages/theta) * t;

		// if event is migration (original or new)
		// follow to source and recalculate living mig_bands
		// if event starts/ends migration band, change mig_rate
		if(mig_source >= 0) {
			// add log-likelihood of migration event
			mig_spr_stats.genetree_delta_lnLd[reconnect] += log(dataSetup.popTree->migBands[mig_band].migRate);
			event = mig_source;
			pop = dataSetup.popTree->migBands[mig_band].sourcePop;
			theta = dataSetup.popTree->pops[pop]->theta*heredity_factor;
			mig_source = -1;
			mig_rate = 0.0;
			num_live_mig_bands = 0;
			for(mig_band = 0; mig_band < dataSetup.popTree->numMigBands; mig_band++) {
				if(dataSetup.popTree->migBands[mig_band].targetPop == pop &&
						dataSetup.popTree->migBands[mig_band].startTime <= age &&
						dataSetup.popTree->migBands[mig_band].endTime   >  age)
				{
					mig_rate += dataSetup.popTree->migBands[mig_band].migRate;
					live_mig_bands[num_live_mig_bands++] = mig_band;
				}
			}

		} else if(event_chains[gen].events[event].type == MIG_BAND_START) {
			mig_rate += dataSetup.popTree->migBands[node_id].migRate;
			live_mig_bands[num_live_mig_bands] = node_id;
			num_live_mig_bands++;
		} else if(event_chains[gen].events[event].type == MIG_BAND_END) {
			mig_rate -= dataSetup.popTree->migBands[node_id].migRate;
			if(num_live_mig_bands == 1 && mig_rate > 0.0000001 && debug) {
				fprintf(stderr, "\ntraceLineage() mig_rate=%g when all mig bands have been emptied at pop=%d, age=%g.\n",mig_rate,pop,age);
				printPopulationTree(dataSetup.popTree, stderr, 1);
				printLocusGenTree(dataState.lociData[gen], stderr, nodePops[gen], nodeEvents[gen]);
				printEventChains(stderr, gen);
				fprintf(stderr, "\n\n********************************************************************\n\n");
			}
			if(num_live_mig_bands == 1) {
				mig_rate = 0.0;
			}
			for(i=0; i< num_live_mig_bands; i++) {
				if(live_mig_bands[i] == node_id) {
					live_mig_bands[i] = live_mig_bands[--num_live_mig_bands];
					break;
				}
			}
		}

		event = event_chains[gen].events[event].next;

	}// end of while

	// add log-likelihood of father coalescent
	mig_spr_stats.genetree_delta_lnLd[reconnect] += log(2/theta);


	return 0;
}

/* replaceMigNodes
   This function is called by UpdateGB_MigSPR, if SPR is accepted.
   - it removes all old migration events
   - it replaces old migration nodes with new ones
   - it removes excess old mignodes
   - it defines new mignodes, if necessary
   This function is the only one which modifies the gen mignodes.
   It uses data in mig_spr_stats.	
 */

int replaceMigNodes(int gen, int node) {

	int	i, j, mig, mig_band,
	max_num_migs = (mig_spr_stats.num_old_migs > mig_spr_stats.num_new_migs) ? mig_spr_stats.num_old_migs : mig_spr_stats.num_new_migs ;

	for(i=0; i<max_num_migs; i++) {
		if(i < mig_spr_stats.num_old_migs) {
			// replace old mignode
			mig = mig_spr_stats.old_migs[i];
			removeEvent(gen, genetree_migs[gen].mignodes[mig].source_event);
			removeEvent(gen, genetree_migs[gen].mignodes[mig].target_event);
			// update migration stats
			mig_band = genetree_migs[gen].mignodes[mig].migration_band;
			genetree_stats[gen].num_migs[mig_band]--;
			genetree_stats_total.num_migs[mig_band]--;

		} else {
			// locate free mignode
			for(mig=0; mig<MAX_MIGS; mig++) {
				if(genetree_migs[gen].mignodes[mig].migration_band < 0)   break;
			}
			if(mig == MAX_MIGS) {
				if(debug) {
					fprintf(stderr, "Error: replaceMigNodes: not enough free migs. Num old %d, num new %d, existing %d, max %d.\n",
							mig_spr_stats.num_old_migs, mig_spr_stats.num_new_migs, genetree_migs[gen].num_migs, MAX_MIGS);
				} else {
					fprintf(stderr, "Fatal Error 0012.\n");
				}
				printGenealogyAndExit(gen,-1);
			}
			genetree_migs[gen].living_mignodes[genetree_migs[gen].num_migs++] = mig;
			genetree_migs[gen].mignodes[mig].gtree_branch = node;
		}

		// at this point mig is an index of either a discarded old mignode or a free one

		if(i<mig_spr_stats.num_new_migs) {
			genetree_migs[gen].mignodes[mig].source_event 	= mig_spr_stats.new_migs_out[i];
			genetree_migs[gen].mignodes[mig].target_event 	= mig_spr_stats.new_migs_in[i];
			genetree_migs[gen].mignodes[mig].age 			= mig_spr_stats.new_migs_ages[i];
			mig_band 											= mig_spr_stats.new_migs_bands[i];
			genetree_migs[gen].mignodes[mig].migration_band = mig_band;
			genetree_migs[gen].mignodes[mig].source_pop 	= dataSetup.popTree->migBands[mig_band].sourcePop;
			genetree_migs[gen].mignodes[mig].target_pop 	= dataSetup.popTree->migBands[mig_band].targetPop;

			event_chains[gen].events[mig_spr_stats.new_migs_out[i]].type = OUT_MIG;
			event_chains[gen].events[mig_spr_stats.new_migs_out[i]].node_id = mig;

			event_chains[gen].events[mig_spr_stats.new_migs_in[i]].type = IN_MIG;
			event_chains[gen].events[mig_spr_stats.new_migs_in[i]].node_id = mig;

			// update migration stats
			genetree_stats[gen].num_migs[mig_band]++;
			genetree_stats_total.num_migs[mig_band]++;

			//         printf("\nCreating new mignode %d for gen %d: mig band %d, branch %d, events (in,out) %d,%d, age %g.",
			//         				mig, gen, mig_band, genetree_migs[gen].mignodes[mig].gtree_branch,
			//         				mig_spr_stats.new_migs_in[i], mig_spr_stats.new_migs_out[i], mig_spr_stats.new_migs_ages[i]);

		} else {
			// brand this mignode as "free" and remove from living list
			genetree_migs[gen].mignodes[mig].migration_band = -1;
			for(j=0; j<genetree_migs[gen].num_migs; j++) {
				if(genetree_migs[gen].living_mignodes[j] == mig) {
					genetree_migs[gen].living_mignodes[j] = genetree_migs[gen].living_mignodes[--genetree_migs[gen].num_migs];
					break;
				}
			}
		}
	}
	return 0;
}


/* considerEventMove
   computes the modifications required for changing a specific genetree
   by moving a given event to a new position (in populations ancestral or 
   descendant to source population). Event can be a coalescence or migration
   event. Computes change in event chains (by adding 
   new event and not removing the original one), and changes in statistics 
   required for fast computation of genetree likelihood. 
   This procedure	updates all fields of genetree_stats_delta and returns the 
   delta in log-likelihood of genetree due to this step.
 */

double considerEventMove(int gen, int instance, int event_id, int source_pop, double original_age, int target_pop, double new_age) {

	int  new_event, bottom_event, top_event, bottom_pop;
	double top_age, bottom_age,  delta_lnLd;


	// find new place
	new_event = createEvent(gen, target_pop, new_age);
	if( new_event < 0 ) {
		if(debug) {
			fprintf(stderr, "Error: Error in creating new event in consider event move.\n");
		} else {
			fprintf(stderr, "Fatal Error 0013.\n");
		}
		printGenealogyAndExit(gen,-1);
	}

#ifdef LOG_STEPS
	fprintf(ioSetup.debugFile, "considerEventMove: gen %d, event %d, pops %d--> %d, ages %g --> %g. New event %d.\n",
			gen, event_id, source_pop,target_pop, original_age, new_age, new_event);
#endif


	genetree_stats_delta[instance].original_event = event_id;
	genetree_stats_delta[instance].updated_event  = new_event;

	if(new_age > original_age) {
		genetree_stats_delta[instance].num_lin_delta = (event_chains[gen].events[event_id].type == OUT_MIG) ? (-1) : (1);
		bottom_event = event_chains[gen].events[event_id].next;
		top_event = new_event;
		bottom_pop = source_pop;
		top_age = new_age;
		bottom_age = original_age;
	} else {
		genetree_stats_delta[instance].num_lin_delta = (event_chains[gen].events[event_id].type == OUT_MIG) ? (1) : (-1);
		bottom_event = event_chains[gen].events[new_event].next;
		top_event = event_id;
		bottom_pop = target_pop;
		top_age = original_age;
		bottom_age = new_age;
	}

#ifdef DEBUG_NODE_CHANGE
	printf("Computing coal stats delta...\n");
#endif	 

	// compute changes in coalescence statistics
	computeCoalStatsDelta(instance, gen, bottom_event, bottom_pop, top_event, genetree_stats_delta[instance].num_lin_delta);

#ifdef DEBUG_NODE_CHANGE
	printf("Computing mig stats delta...\n");
#endif	 
	// compute changes in migration statistics across relevant migration bands
	computeMigStatsDelta(instance, bottom_age, bottom_pop, top_age, genetree_stats_delta[instance].num_lin_delta);

#ifdef DEBUG_NODE_CHANGE
	printf("Computing delta in log-likelihood...\n");
#endif	 

	delta_lnLd = computeDeltaLnLd(gen, instance);

#ifdef DEBUG_NODE_CHANGE
	printf("Done.\n");
#endif	 

	// if event moved was a coalescent event and it was moved to another population
	// then take consider density of coalescence in likelihood.
	// Note: migration events cannot be moved to other bands.
	if(event_chains[gen].events[event_id].type == COAL && source_pop != target_pop) {
		delta_lnLd += log(dataSetup.popTree->pops[source_pop]->theta/dataSetup.popTree->pops[target_pop]->theta);
	}

	return delta_lnLd;
}

/*	computeDeltaLnLd
	Computes the delta in log-likelihood of genetree when considering a sampling move.
	This is based on differences of statistics already computed in 'genetree_stats_delta[instance]'.
	Does not take into account other changes (such as moving a coalescence event from one
	population to the other).
 */

double computeDeltaLnLd(int gen, int instance) {
	int i;
	double delta_lnLd = 0 , heredity_factor = 1;

	for(i=0; i<genetree_stats_delta[instance].num_pops_changed; i++) {
		delta_lnLd -= genetree_stats_delta[instance].coal_stats_delta[i] / dataSetup.popTree->pops[ genetree_stats_delta[instance].pops_changed[i] ]->theta;
	}

	delta_lnLd /= heredity_factor;

	// migration - difference
	for(i=0; i<genetree_stats_delta[instance].num_mig_bands_changed; i++) {
		delta_lnLd -= genetree_stats_delta[instance].mig_stats_delta[i] * dataSetup.popTree->migBands[ genetree_stats_delta[instance].mig_bands_changed[i] ].migRate;
	}

	return delta_lnLd;
}


/* acceptEventChainChanges
   Removes the original event (and keeps the new one).
   Updates all genetree statistics according to delta.
   Updates num_lineages for all relevant events.
 */
int acceptEventChainChanges(int gen, int instance) {
	int i, pop, mig_band, res=0;
	double heredity_factor = 1;


#ifdef DEBUG_NODE_CHANGE
	printf("Accepting changes in gen %d. Will replace original event %d with updated %d (type %d).\n",
			gen, genetree_stats_delta[instance].original_event,
			genetree_stats_delta[instance].updated_event,
			event_chains[gen].events[genetree_stats_delta[instance].original_event].type);
#endif

	// change genetree statistics
	for(i=0; i<genetree_stats_delta[instance].num_pops_changed; i++) {
		pop = genetree_stats_delta[instance].pops_changed[i];
		genetree_stats[gen].coal_stats[pop] += genetree_stats_delta[instance].coal_stats_delta[i];
		genetree_stats_total.coal_stats[pop]  += genetree_stats_delta[instance].coal_stats_delta[i]/heredity_factor;
	}
	for(i=0; i<genetree_stats_delta[instance].num_mig_bands_changed; i++) {
		mig_band = genetree_stats_delta[instance].mig_bands_changed[i];
		genetree_stats[gen].mig_stats[mig_band] += genetree_stats_delta[instance].mig_stats_delta[i];
		genetree_stats_total.mig_stats[mig_band]  += genetree_stats_delta[instance].mig_stats_delta[i];
	}

	// change number of lineages in affected interval

	i = genetree_stats_delta[instance].num_changed_events-1;

	// if original event is in this list (meaning that it must be last), skip it.
	if(genetree_stats_delta[instance].changed_events[i] == genetree_stats_delta[instance].original_event)
		i--;

	for( ; i>=0 ; i--) {
		event_chains[gen].events[genetree_stats_delta[instance].changed_events[i]].num_lineages += genetree_stats_delta[instance].num_lin_delta;
	}

	if(genetree_stats_delta[instance].updated_event >= 0) {


		// update pointers between event and node
		event_chains[gen].events[genetree_stats_delta[instance].updated_event].node_id = event_chains[gen].events[genetree_stats_delta[instance].original_event].node_id;
		event_chains[gen].events[genetree_stats_delta[instance].updated_event].type = event_chains[gen].events[genetree_stats_delta[instance].original_event].type;

		switch(event_chains[gen].events[genetree_stats_delta[instance].updated_event].type) {
			case(COAL):
#ifdef DEBUG_NODE_CHANGE
    																																														  printf("Node %d in gen %d corresponds now to event %d (rather than event %d).\n",
    																																																  event_chains[gen].events[genetree_stats_delta[instance].updated_event].node_id,
																																																	  gen,
																																																	  genetree_stats_delta[instance].updated_event,
																																																	  genetree_stats_delta[instance].original_event);
#endif
			nodeEvents[gen][ event_chains[gen].events[genetree_stats_delta[instance].updated_event].node_id ] = genetree_stats_delta[instance].updated_event;
			//nodes[event_chains[gen].events[genetree_stats_delta[instance].updated_event].node_id].event_id = genetree_stats_delta[instance].updated_event;
			break;
			case(OUT_MIG):
    																																														  genetree_migs[gen].mignodes[event_chains[gen].events[genetree_stats_delta[instance].updated_event].node_id].source_event =
    																																																  genetree_stats_delta[instance].updated_event;
			break;
			case(IN_MIG):
    																																														  genetree_migs[gen].mignodes[event_chains[gen].events[genetree_stats_delta[instance].updated_event].node_id].target_event =
    																																																  genetree_stats_delta[instance].updated_event;
			break;
			default:
				if(debug) {
					fprintf(stderr, "\nError: acceptEventChainChanges: Illegal event type %d.\n", event_chains[gen].events[genetree_stats_delta[instance].updated_event].type);
				} else {
					fprintf(stderr, "Fatal Error 0014.\n");
				}
				printGenealogyAndExit(gen,-1);
		}

		res = removeEvent(gen, genetree_stats_delta[instance].original_event);
	}




	// initialize genetree_stats_delta
	genetree_stats_delta[instance].num_pops_changed = 0;
	genetree_stats_delta[instance].num_mig_bands_changed = 0;
	genetree_stats_delta[instance].num_changed_events = 0;
	genetree_stats_delta[instance].num_lin_delta = 0;
	genetree_stats_delta[instance].original_event = -1;
	genetree_stats_delta[instance].updated_event = -1;

	return res;
}


/* rejectEventChainChanges
   Simply removes the new event added.
 */
int rejectEventChainChanges(int gen, int instance) {

	int res = 1;


	if(genetree_stats_delta[instance].updated_event >= 0) {
		res = removeEvent(gen, genetree_stats_delta[instance].updated_event);
#ifdef DEBUG_NODE_CHANGE
		printf("Rejecting changes in gen %d. Removing updated event %d.\n",gen, genetree_stats_delta[instance].updated_event);
#endif
	}
	// initialize genetree_stats_delta
	genetree_stats_delta[instance].num_pops_changed = 0;
	genetree_stats_delta[instance].num_mig_bands_changed = 0;
	genetree_stats_delta[instance].num_changed_events = 0;
	genetree_stats_delta[instance].num_lin_delta = 0;
	genetree_stats_delta[instance].original_event = -1;
	genetree_stats_delta[instance].updated_event = -1;


	return res;

}

/* removeEvent
   Removes specified event (coalescence of migration) from chain and adds it to free events.
 */
int removeEvent(int gen, int event){
	int 	next_event = event_chains[gen].events[event].next,
			prev_event = event_chains[gen].events[event].prev;


	event_chains[gen].events[next_event].elapsed_time += event_chains[gen].events[event].elapsed_time;
	event_chains[gen].events[next_event].prev = prev_event;

	if(prev_event < 0) {
		// weird way for getting population of event
		// MAYBE CHANGE!!!
		for(prev_event=next_event; event_chains[gen].events[prev_event].type != END_CHAIN; prev_event = event_chains[gen].events[prev_event].next)
		{ ; }

		event_chains[gen].first_event[event_chains[gen].events[prev_event].node_id] = next_event;
	} else {
		event_chains[gen].events[prev_event].next = next_event;
	}

	// add event to chain of free events
	next_event = event_chains[gen].free_events;
	event_chains[gen].events[event].next = next_event;
	event_chains[gen].events[next_event].prev = event;
	event_chains[gen].free_events = event;

	event_chains[gen].events[event].elapsed_time = 0;
	event_chains[gen].events[event].num_lineages  = 0;
	event_chains[gen].events[event].node_id  = -1;

#ifdef DEBUG_EVENT_CHAIN
	printf("Removing event %d from gen %d. First free event is %d.\n", event, gen, event_chains[gen].free_events);
#endif

	return 0;
}

/* createEventBefore
   Creates a new event before specified event.
   Returns the id of the event.
 */

int createEventBefore(int gen, int pop, int event, double elapsed_time) {

	int prev_event = event_chains[gen].events[event].prev,
			new_event = event_chains[gen].free_events;

	event_chains[gen].free_events = event_chains[gen].events[new_event].next;

	if(event_chains[gen].free_events < 0) {
		if(debug) {
			fprintf(stderr, "\nError: Empty event pool in gen %d.\n", gen);
		} else {
			fprintf(stderr, "Fatal Error 0015.\n");
		}
		printGenealogyAndExit(gen,-1);
	}

	// make changes in elapsed time and pointers
	event_chains[gen].events[new_event].next = event;
	event_chains[gen].events[new_event].prev = prev_event;
	event_chains[gen].events[new_event].num_lineages = event_chains[gen].events[event].num_lineages;
	event_chains[gen].events[new_event].elapsed_time = elapsed_time;
	event_chains[gen].events[new_event].type = DUMMY;
	event_chains[gen].events[event].prev = new_event;
	event_chains[gen].events[event].elapsed_time -= elapsed_time;
	if(prev_event < 0) {
		event_chains[gen].first_event[pop] = new_event;
	} else {
		event_chains[gen].events[prev_event].next = new_event;
	}

#ifdef DEBUG_EVENT_CHAIN
	printf("Added new event %d to gen %d. Next free event %d.\n",new_event,gen, event_chains[gen].free_events);
#endif

	return new_event;
}


/* createEvent
   Creates a new event in specified population at given time.
   Changes only elapsed time of subsequent event. Makes no other
   changes to event chain. num_lineanges of new event is as the one 
   of subsequent event.
   Returns the id of the event.
 */

int createEvent(int gen, int pop, double age) {

	int event;
	double delta_time = age - dataSetup.popTree->pops[pop]->age;

	if(delta_time < 0) {
		if(debug) {
			fprintf(stderr, "\nError: createEvent: time specified %g is smaller than age of target population %d (%g).\n",
					age, pop, dataSetup.popTree->pops[pop]->age);
		} else {
			fprintf(stderr, "Fatal Error 0016.\n");
		}
		return (-1);
	}

	if(pop != dataSetup.popTree->rootPop && age > dataSetup.popTree->pops[pop]->father->age + 0.000001) {
		if(debug) {
			fprintf(stderr, "\nError: createEvent: time specified %g is greater than age of parent population %d (%g).\n",
					age, dataSetup.popTree->pops[pop]->father->id, dataSetup.popTree->pops[pop]->father->age);
		} else {
			fprintf(stderr, "Fatal Error 0017.\n");
		}
		return (-1);
	}

	// find spot for new event
	for(event = event_chains[gen].first_event[pop];
			event_chains[gen].events[event].type != END_CHAIN && event_chains[gen].events[event].elapsed_time < delta_time;
			//			(event>=0 && event_chains[gen].events[event].elapsed_time < delta_time);
			event = event_chains[gen].events[event].next) {

		delta_time -= event_chains[gen].events[event].elapsed_time;

	}

	if(event_chains[gen].events[event].elapsed_time < delta_time) {
		if(event_chains[gen].events[event].elapsed_time < delta_time - 0.000001) {
			if(debug) {
				fprintf(stderr, "\nError: createEvent: trying to insert new event in pop %d, gen %d at time %g, %g above END_CHAIN event.\n",
						pop, gen, age, delta_time - event_chains[gen].events[event].elapsed_time);
			} else {
				fprintf(stderr, "Fatal Error 0018.\n");
			}
			printGenealogyAndExit(gen,-1);
		}
		delta_time = event_chains[gen].events[event].elapsed_time;
	}

	return createEventBefore(gen,pop,event,delta_time);

}


/*	adjustRootEvents
	adjusts elapsed times of all root events after every UpdateTau() and mixing()	
 */
int adjustRootEvents() {
	int gen, event;
	double age;

	for(gen=0; gen<dataSetup.numLoci; gen++) {

		event = event_chains[gen].first_event[dataSetup.popTree->rootPop];
		age = dataSetup.popTree->pops[ dataSetup.popTree->rootPop ]->age;
		while(event_chains[gen].events[event].next >= 0) {
			age += event_chains[gen].events[event].elapsed_time;
			event = event_chains[gen].events[event].next;
		}

		event_chains[gen].events[event].elapsed_time = OLDAGE - age;
	}

	return 0;
}

/*	computeMigStatsDelta
	Computes the difference in migration statistics in a genetree caused by changing the number of 
	lineages in a consecutive interval in time by some constant number (typically +1 or -1).
	The interval is indicated by starting point - bottom_age at bottom_pop and end point
	top_age (at ancestral population of bottom_pop).
	'instance' is the index of genetree_stats_delta to be updated.
	This procedure updates the following fields of genetree_stats_delta[instance]:
	- num_mig_bands_changed
	- mig_bands_changed array
	- mig_stats_delta array
 */
int computeMigStatsDelta(int instance, double bottom_age, int bottom_pop, double top_age, int num_lins_delta){

	int mig_band;
	double delta_time;

	genetree_stats_delta[instance].num_mig_bands_changed = 0;
	for(mig_band = 0; mig_band<dataSetup.popTree->numMigBands; mig_band++) {
		// consider only bands coming into populations ancestral to bottom one
		if(!dataSetup.popTree->pops[ dataSetup.popTree->migBands[mig_band].targetPop ]->isAncestralTo[bottom_pop])
			continue;
		// find intersection interval between migration band and affected interval
		delta_time = min2(dataSetup.popTree->migBands[mig_band].endTime , top_age) - max2(dataSetup.popTree->migBands[mig_band].startTime , bottom_age);

		if(delta_time <= 0)
			continue;

		// mig_band is affected by change. Note that and compute delta
		// mig_stat = num_lins * elapsed_time
		// and so mig_stats_delta = num_lin_delta * delta_time
		genetree_stats_delta[instance].mig_bands_changed[genetree_stats_delta[instance].num_mig_bands_changed] = mig_band;
		genetree_stats_delta[instance].mig_stats_delta[genetree_stats_delta[instance].num_mig_bands_changed] = num_lins_delta * delta_time;
		genetree_stats_delta[instance].num_mig_bands_changed++;
	}

	return 0;

}
/*	computeCoalStatsDelta
	Computes the difference in coalescence statistics in a genetree caused by changing the number of 
	lineages in a series of consecutive events by some constant number (typically +1 or -1).
	bottom_event is assumed to be below top_event but not necessarily in the same population.
	'instance' is the index of genetree_stats_delta to be updated.
	This procedure updates the following fields of genetree_stats_delta[instance]:
	- num_changed_events
	- changed_events array
	- num_pops_changed
	- pops_changed array
	- coal_stats_delta array
 */

int computeCoalStatsDelta(int instance, int gen, int bottom_event, int bottom_pop, int top_event, int num_lins_delta){

	int 	pop = bottom_pop, event = bottom_event;

	// initialize genetree_stats_delta
	genetree_stats_delta[instance].num_pops_changed = 1;
	genetree_stats_delta[instance].pops_changed[0] = pop;
	genetree_stats_delta[instance].coal_stats_delta[0] = 0;

	genetree_stats_delta[instance].num_changed_events = 0;


	while(event >= 0) {

		// change coalescence stats for event
		// coal_stats = num_lins*(num_lins-1)*elapsed_time
		// so coal_stats_delta = ( num_lins_delta^2 - num_lins_delta + 2*num_lins_delta*num_lins ) * elapsed_time
		genetree_stats_delta[instance].coal_stats_delta[genetree_stats_delta[instance].num_pops_changed-1] +=
				num_lins_delta * (num_lins_delta - 1 + 2*event_chains[gen].events[event].num_lineages)
				* event_chains[gen].events[event].elapsed_time;

		genetree_stats_delta[instance].changed_events[genetree_stats_delta[instance].num_changed_events] = event;
		genetree_stats_delta[instance].num_changed_events++;

		// halting condition for loop
		if(event == top_event)
			break;

		event = event_chains[gen].events[event].next;
		// if you reach top of population, move to next one
		if(event < 0) {
			if(dataSetup.popTree->pops[pop]->father == NULL) {
				if(debug) {
					fprintf(stderr, "\nError: computeCoalStatsDelta: Could not find top event %d from bottom event %d (bottom pop %d) in gen %d.\n",top_event, bottom_event,bottom_pop,gen);
				} else {
					fprintf(stderr, "Fatal Error 0019.\n");
				}
				printGenealogyAndExit(gen,-1);
			}
			pop = dataSetup.popTree->pops[pop]->father->id;
			event = event_chains[gen].first_event[pop];
			genetree_stats_delta[instance].pops_changed[genetree_stats_delta[instance].num_pops_changed] = pop;
			genetree_stats_delta[instance].coal_stats_delta[genetree_stats_delta[instance].num_pops_changed] = 0;
			genetree_stats_delta[instance].num_pops_changed++;
		}

	}

	return 0;
}


/* populationPostOrder
   Computes pots-order for population subtreetree rooted at pop.
   Writes down the post order in specified array.
   Returns size of subtree.
   A recursive procedure.
 */
int populationPostOrder(int pop, int* ordered_pops) {
	int size;

	// halting condition
	if(pop < dataSetup.popTree->numCurPops) {
		ordered_pops[0] = pop;
		return 1;
	}

	// compute post-order for every subtree and add root
	size  = populationPostOrder(dataSetup.popTree->pops[pop]->sons[0]->id, ordered_pops);
	size += populationPostOrder(dataSetup.popTree->pops[pop]->sons[1]->id, ordered_pops+size);
	ordered_pops[size] = pop;

	return size+1;
}


/*	constructEventChain
	Constructs event-chains for a given gen.
	Constructs everything from scratch. 
	Typically used only for initial genetrees or for testing.
	Records number of lineages only for first events in leaf populations.
	The rest are recorded by computeGenetreeStats
 */
int constructEventChain(int gen) {

	int 	i, pop, mig, node, event;
	double age;


	// initialize event chains (with a single END_CHAIN event)
	for(pop=0; pop<dataSetup.popTree->numPops; pop++) {
		event_chains[gen].events[pop].type = END_CHAIN;
		event_chains[gen].events[pop].next = -1;
		event_chains[gen].events[pop].prev = -1;
		event_chains[gen].events[pop].node_id = pop;
		// it is important to initialize 0 incoming lineages per population
		event_chains[gen].events[pop].num_lineages = 0;
		if(pop == dataSetup.popTree->rootPop) {
			event_chains[gen].events[pop].elapsed_time = OLDAGE - dataSetup.popTree->pops[dataSetup.popTree->rootPop]->age;
		}
		else {
			event_chains[gen].events[pop].elapsed_time = dataSetup.popTree->pops[pop]->father->age - dataSetup.popTree->pops[pop]->age;
		}
		event_chains[gen].first_event[pop] = pop;
		// this should stay constant
		event_chains[gen].last_event[pop] = pop;
	}


	// initialize free event chain
	event_chains[gen].free_events = dataSetup.popTree->numPops;
	event_chains[gen].events[dataSetup.popTree->numPops].prev = -1;
	event_chains[gen].events[event_chains[gen].total_events-1].next = -1;
	for(event = dataSetup.popTree->numPops; event < event_chains[gen].total_events-1; event++) {
		event_chains[gen].events[event].next = event+1;
		event_chains[gen].events[event+1].prev = event;
	}

	// migration band events
	for(mig=0; mig<dataSetup.popTree->numMigBands; mig++) {
		// start band event
		pop = dataSetup.popTree->migBands[mig].targetPop;
		age = dataSetup.popTree->migBands[mig].startTime;
		event=createEvent(gen, pop, age);
		if( event < 0 ) {
			if(debug) {
				fprintf(stderr, "Error: Unable to create new migration band start event.\n");
			} else {
				fprintf(stderr, "Fatal Error 0020.\n");
			}
			printGenealogyAndExit(gen,-1);
		}
		event_chains[gen].events[event].type = MIG_BAND_START;
		event_chains[gen].events[event].node_id = mig;
		// end band event (pop, pop1 cannot be root)
		age = dataSetup.popTree->migBands[mig].endTime;
		event=createEvent(gen, pop, age);
		if( event < 0 ) {
			if(debug) {
				fprintf(stderr, "Error: Unable to create new migration band end event.\n");
			} else {
				fprintf(stderr, "Fatal Error 0021.\n");
			}
			printGenealogyAndExit(gen,-1);
		}
		event_chains[gen].events[event].type = MIG_BAND_END;
		event_chains[gen].events[event].node_id = mig;
	}


	// samples start (for ancient samples)
	for(pop=0; pop<dataSetup.popTree->numCurPops; pop++) {
		event=createEvent(gen, pop, dataSetup.popTree->pops[pop]->sampleAge);
		event_chains[gen].events[event].type = SAMPLES_START;
	}

	// migration node events
	for(i=0; i<genetree_migs[gen].num_migs; i++) {
		mig = genetree_migs[gen].living_mignodes[i];
#ifdef DEBUG_EVENT_CHAIN
		printf("living migration %d: %d in gen %d.\n",i,mig,gen);
#endif
		// incoming migration event
		pop = genetree_migs[gen].mignodes[mig].target_pop;
		age = genetree_migs[gen].mignodes[mig].age;
		event=createEvent(gen, pop, age);
		if( event < 0 ) {
			if(debug) {
				fprintf(stderr, "Error: Unable to create new in migration event.\n");
			} else {
				fprintf(stderr, "Fatal Error 0022.\n");
			}
			printGenealogyAndExit(gen,-1);
		}
		event_chains[gen].events[event].type = IN_MIG;
		event_chains[gen].events[event].node_id = mig;
		genetree_migs[gen].mignodes[mig].target_event = event;
		// outgoing migration event
		pop = genetree_migs[gen].mignodes[mig].source_pop;
		event=createEvent(gen, pop, age);
		if( event < 0 ) {
			if(debug) {
				fprintf(stderr, "Error: Unable to create new out migration event.\n");
			} else {
				fprintf(stderr, "Fatal Error 0023.\n");
			}
			printGenealogyAndExit(gen,-1);
		}
		event_chains[gen].events[event].type = OUT_MIG;
		event_chains[gen].events[event].node_id = mig;
		genetree_migs[gen].mignodes[mig].source_event = event;
	}


	// coalescent node events + initialize incoming lineages per leaf populations
	for(node=0; node<2*dataSetup.numSamples-1; node++) {
		pop = nodePops[gen][node];
		// leaf
		if(node < dataSetup.numSamples) {
#ifdef DEBUG_EVENT_CHAIN
			printf("L%d ",pop);
#endif
			/*
      event = event_chains[gen].first_event[pop]; 
	  while(event_chains[gen].events[event].type != SAMPLES_START && event_chains[gen].events[event].type != END_CHAIN) {
		event = event_chains[gen].events[event].next;
	  }
	  if(event_chains[gen].events[event].type == END_CHAIN) {
		if(debug) {
	        fprintf(stderr, "\nError: constructEventChain: could not find start event in pop %d.\n",pop);
        	fprintf(stderr, ".\n");
		} else {
			fprintf(stderr, "Fatal Error 0102.\n");
		}
	  }
	  event = event_chains[gen].events[event].next;
      // although first event changes throughout this process,
      // because createEvent copies num_lineages of subsequent event, 
      // the number of starting lineages per leaf population is recorded correctly.
      event_chains[gen].events[event].num_lineages++;
			 */
		}
		//coalescent node
		else {
			age = getNodeAge(dataState.lociData[gen],node);
			//			age = gnodes[gen][node].age;
#ifdef DEBUG_EVENT_CHAIN
			printf("C%d %d %f ", pop, node, age);
#endif
			event = createEvent(gen, pop, age);
			if( event < 0 ) {
				if(debug) {
					fprintf(stderr, "Error: Unable to create new coalescence event.\n");
				} else {
					fprintf(stderr, "Fatal Error 0024.\n");
				}
				printGenealogyAndExit(gen,-1);
			}
			event_chains[gen].events[event].type = COAL;
			event_chains[gen].events[event].node_id = node;
			nodeEvents[gen][node] = event;
		}
	}
#ifdef DEBUG_EVENT_CHAIN
	printf("\n");
#endif
	return 0;
}


/*	computeTotalStats
	Computes all statistics in genetree_stats_total using pre-calculated genetree
	statistics for every gen. Note that coal_stats takes gen heredity into 
	consideration in the total statistics.
 */

int computeTotalStats() {
	double heredity_factor = 1;
	int pop, mig_band,gen;

	//initialize values
	for(pop=0; pop<dataSetup.popTree->numPops; pop++) {
		genetree_stats_total.coal_stats[pop] = 0;
		genetree_stats_total.num_coals[pop]  = 0;
	}
	for(mig_band=0; mig_band<dataSetup.popTree->numMigBands; mig_band++) {
		genetree_stats_total.mig_stats[mig_band] = 0;
		genetree_stats_total.num_migs[mig_band]  = 0;
	}

	for(gen=0; gen<dataSetup.numLoci; gen++) {
		//		if(data.est_heredity)
		//			heredity_factor = data.heredity[gen];

		for(pop=0; pop<dataSetup.popTree->numPops; pop++) {
			genetree_stats_total.coal_stats[pop] += genetree_stats[gen].coal_stats[pop]/heredity_factor;
			genetree_stats_total.num_coals[pop]  += genetree_stats[gen].num_coals[pop];
		}
		for(mig_band=0; mig_band<dataSetup.popTree->numMigBands; mig_band++) {
			genetree_stats_total.mig_stats[mig_band] += genetree_stats[gen].mig_stats[mig_band];
			genetree_stats_total.num_migs[mig_band]  += genetree_stats[gen].num_migs[mig_band];
		}
	}


	return 0;
}
/*** end of computeNodeStats ***/


/*	computeNodeStats
	Computes all statistics for coalescent nodes in all genealogies to find model violations
 */

int computeNodeStats() {

	int res, gen, leaf1, leaf2, node, pop;


	/*** point to arrays/matrices in genetree_node_stats ***/
	double*** probCoalMatrix      = genetree_node_stats.probCoalMatrix;       // Prob[sample pair coalesce in pop]
	double*** probFirstCoalMatrix = genetree_node_stats.probFirstCoalMatrix;  // Prob[sample pair is the first coalescence in pop]
	double*** coalTimeMatrix      = genetree_node_stats.coalTimeMatrix;       // mean coal time for coalescence of pair (given they coalesce in pop)
	// auxilliary arrays
	int**     lcaMatrix           = genetree_node_stats.lcaMatrix;            // matrix of LCAs for all pairs of leaves
	int*      leafArray_aux       = genetree_node_stats.leafArray;            // array of leaves for LCA computation
	int*      firstNodesInPop     = genetree_node_stats.firstNodesInPop;      // array with the first coal node in each pop
	double*   firstNodesAges      = genetree_node_stats.firstNodesAges;       // array with the age of first coal node in each pop
	double*   nodeAges            = genetree_node_stats.nodeAges;             // array with the ages of all internal nodes

	/*** initialize matrices ***/
	for(leaf1=0; leaf1<dataSetup.numSamples; leaf1++) {
		for(leaf2=0; leaf2<dataSetup.numSamples; leaf2++) {
			for(pop=0; pop<dataSetup.popTree->numPops; pop++) {
				coalTimeMatrix[leaf1][leaf2][pop] = 0.0;
				probFirstCoalMatrix[leaf1][leaf2][pop] = 0.0;
				probCoalMatrix[leaf1][leaf2][pop] = 0.0;
			}// end of for(pop)
			lcaMatrix[leaf1][leaf2] = -1;
		}// end of for(leaf2)
	}// end of for(leaf1)

	/*** collect stats per gen - sum ***/
	for(gen=0; gen<dataSetup.numLoci; gen++) {
		res = computePairwiseLCAs (dataState.lociData[gen], lcaMatrix, leafArray_aux);
		if(!res) {
			if(debug) {
				fprintf(stderr, "\nError: computeNodeStats: error computing LCA matrix for gen %d.\n", gen);
			} else {
				fprintf(stderr, "Fatal Error 0202.\n");
			}
			printGenealogyAndExit(gen,-1);
		}
		/*** figure out node ages and first node in pop ***/
		for(pop=0; pop<dataSetup.popTree->numPops; pop++) {
			firstNodesInPop[pop] = -1;
			firstNodesAges[pop] = 0.0;
		}// end of for(pop)

		for(node=dataSetup.numSamples; node<2*dataSetup.numSamples-1; node++) {
			nodeAges[node] = getNodeAge(dataState.lociData[gen],node);
			pop = nodePops[gen][node];
			if(firstNodesInPop[pop] < 0 || nodeAges[node] < firstNodesAges[pop]) {
				firstNodesInPop[pop] = node;
				firstNodesAges[pop]  = nodeAges[node];
			}
		}// end of for(node)


		/*** convert LCAs to pops and figure first node in pop ***/
		for(leaf1=0; leaf1<dataSetup.numSamples; leaf1++) {
			for(leaf2=0; leaf2<dataSetup.numSamples; leaf2++) {
				node = lcaMatrix[leaf1][leaf2];
				if(node<0) {
					if(leaf1 != leaf2) {
						if(debug) {
							fprintf(stderr, "\nError: computeNodeStats: invalid LCA matrix for gen %d.\n", gen);
						} else {
							fprintf(stderr, "Fatal Error 0203.\n");
						}
					} else {
						continue;
					}
				}

				pop  = nodePops[gen][node];
				probCoalMatrix[leaf1][leaf2][pop] += 1.0;
				coalTimeMatrix[leaf1][leaf2][pop] += nodeAges[node];
				if(firstNodesInPop[pop] == node) {
					probFirstCoalMatrix[leaf1][leaf2][pop] += 1.0;
				}
			}// end of for(leaf2)
		}// end of for(leaf1)


	}// end of for(gen)

	/*** normalize stats ***/
	for(leaf1=0; leaf1<dataSetup.numSamples; leaf1++) {
		for(leaf2=0; leaf2<dataSetup.numSamples; leaf2++) {
			for(pop=0; pop<dataSetup.popTree->numPops; pop++) {
				if(probCoalMatrix[leaf1][leaf2][pop] > 0) {
					coalTimeMatrix[leaf1][leaf2][pop] /= probCoalMatrix[leaf1][leaf2][pop];
				}
				probFirstCoalMatrix[leaf1][leaf2][pop] /= dataSetup.numLoci;
				probCoalMatrix[leaf1][leaf2][pop] /= dataSetup.numLoci;
			}// end of for(pop)
		}// end of for(leaf2)
	}// end of for(leaf1)

	return 0;
}
/*** end of computeNodeStats ***/













void debug_print_cladestats(int clade, int gen, char* methodName){
	printf("=== %s, clade:%s, gen:%d, num_events:%d, num_coals_total:%d === >>\n",
			methodName, dataSetup.popTree->pops[clade]->name, gen, clade_stats[clade].num_events, clade_stats[clade].num_coals_total);
	for (int i = 0 ; i < clade_stats[clade].num_events ; i++){
		printf("age:%0.25f, elapsed_time:%0.25f, num_lin:%d, type:%s \n",
				clade_stats[clade].sorted_ages[i], clade_stats[clade].elapsed_times[i],
				clade_stats[clade].num_lineages[i],	getEventTypeName(clade_stats[clade].event_types[i]));
	}
	printf("====================================================================================\n");
}

#define TRUE 1 //TODO - where should these consts be?!
#define FALSE 0
int	isLeafPopulation(int pop){
	Population *population, *left_son, *right_son;

	population = dataSetup.popTree->pops[pop];

	left_son = population->sons[LEFT];
	right_son = population->sons[RIGHT];

	if (left_son || right_son){
		return FALSE;
	} else {
		return TRUE;
	}
}

void test_validateRootCladeVsFlatStats(){
	int rootPop = dataSetup.popTree->rootPop;
	double rootCladeCoalStats = clade_stats[rootPop].coal_stats_total;
	double flatCoalStats = genetree_stats_flat.coal_stats_flat;

	if (rootCladeCoalStats != flatCoalStats) {
		fprintf(stderr, "Error: got root_clade_stats=%0.35f and flat_stats=%0.35f but expected same result.\n",
				rootCladeCoalStats, flatCoalStats);
		exit(-1);
	}
}


double getCoalStats(double* elapsed_times, int* num_lineages, int size){
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

/***
 * This method is currently conjoined to "appendPopToClade". It must run inside it, so don't move it.
 * TODO - extract this test from "appendPopToClade"
 */
void test_compareGphocsVsCladePopCoalStats(int clade, int gen, int startingPoint){
	double gphocsValue = genetree_stats[gen].coal_stats[clade];
	double cladeValue = getCoalStats(clade_stats[clade].elapsed_times + startingPoint,
			clade_stats[clade].num_lineages + startingPoint, clade_stats[clade].num_events);

	if (gphocsValue != cladeValue) {
		printf("Error when comparing gphocs pop_coal_stats with clade pop_coal_stats:\n"
				"got gphocs value:%0.25f and clade value:%0.25f but expected same result.\n",
				gphocsValue, cladeValue);
				exit(-1);
	}
}


/***
 * cleans up all variables in clade_stats
 */
void initCladeStats(){
	for (int clade = 0; clade < dataSetup.popTree->numPops; clade++) {
		initSpecificCladeStats(clade);
	}
}
void initSpecificCladeStats(int clade){
	clade_stats[clade].coal_stats_total = 0.0;
	clade_stats[clade].mig_stats_total 	= 0.0;
	clade_stats[clade].num_coals_total 	= 0;
	clade_stats[clade].num_migs_total 	= 0;
	clade_stats[clade].num_events 		= 0;
}





/***
 * author - Ron visbord
 *
 * computeFlatStats
	Computes all statistics in genetree_stats_flat needed for likelihood of
	gene trees under null model with single population
 */

void computeFlatStats() {
	double deltaT, coalStat, migStat;
	int gen, pop, mig_band, i, numLins;

	// number of events (migration, coalescent)
	genetree_stats_flat.num_migs_total = 0;
	genetree_stats_flat.num_coals_total = 0;
	for(mig_band=0; mig_band<dataSetup.popTree->numMigBands; mig_band++) {
		genetree_stats_flat.num_migs_total += genetree_stats_total.num_migs[mig_band];
	}
	for(pop=0; pop<dataSetup.popTree->numPops; pop++) {
		genetree_stats_flat.num_coals_total += genetree_stats_total.num_coals[pop];
	}

	// coalescent and migration stats
	coalStat = 0.0;
	migStat = 0.0;
	for(gen=0; gen<dataSetup.numLoci; gen++) {

		if(0 == getSortedAges(dataState.lociData[gen], genetree_stats_flat.sortedAgesArray)) {
			if(debug) {
				fprintf(stderr, "\nError: computeFlatStats: error sorting node ages for gen %d.\n", gen);
			} else {
				fprintf(stderr, "Fatal Error 0201.\n");
			}
			printGenealogyAndExit(gen,-1);
		}
		for(i=0, numLins=dataSetup.numSamples ; numLins>1; numLins--, i++) {
			deltaT = (i==0) ? genetree_stats_flat.sortedAgesArray[0] : genetree_stats_flat.sortedAgesArray[i] - genetree_stats_flat.sortedAgesArray[i-1];
			migStat += deltaT*numLins;
			coalStat += deltaT*numLins*(numLins-1);
		}
	}
	genetree_stats_flat.coal_stats_flat = coalStat;
	genetree_stats_flat.mig_stats_flat  = migStat;
}

void computeCladeStats() {

	initCladeStats();
	//	computeCladeNumMigs();
	computeCladeNumCoals();
	computeCladeCoalStats();
}


void computeCladeNumCoals(){
	computeCladeNumCoals_rec(dataSetup.popTree->rootPop);
}

void computeCladeNumCoals_rec(int pop){
	int leftSon, rightSon;

	if (isLeafPopulation(pop)){
		clade_stats[pop].num_coals_total = genetree_stats_total.num_coals[pop];
	} else{
		leftSon = dataSetup.popTree->pops[pop]->sons[LEFT]->id;
		rightSon = dataSetup.popTree->pops[pop]->sons[RIGHT]->id;

		computeCladeNumCoals_rec(leftSon);
		computeCladeNumCoals_rec(rightSon);

		clade_stats[pop].num_coals_total = genetree_stats_total.num_coals[pop] +
				clade_stats[leftSon].num_coals_total + clade_stats[rightSon].num_coals_total;
	}
}


void computeCladeCoalStats(){
	for(int gen=0; gen<dataSetup.numLoci; gen++) {
		computeCladeCoalStats_rec(dataSetup.popTree->rootPop, gen);
	}
}

void computeCladeCoalStats_rec(int clade, int gen) {
	int leftSon, rightSon;

	if (isLeafPopulation(clade)){
		fillupLeafCladeStats(clade, gen);
	} else{
		leftSon = dataSetup.popTree->pops[clade]->sons[LEFT]->id;
		rightSon = dataSetup.popTree->pops[clade]->sons[RIGHT]->id;

		computeCladeCoalStats_rec(leftSon, gen);
		computeCladeCoalStats_rec(rightSon, gen);

		fillupCladeStats(clade, gen);
	}
//	debug_print_cladestats(clade, gen, "computeCladeCoalStats_rec");
}

void fillupLeafCladeStats(int clade, int gen){
	appendPopToClade(clade, gen, 0); // since this is a leaf, starting point of the clade_stats arrays should be 0
}

void appendPopToClade(int clade, int gen, int startingPoint){
	int i = startingPoint;
	int event = event_chains[gen].first_event[clade];
	double cladeStartTime = dataSetup.popTree->pops[clade]->age;
	double eventAge = cladeStartTime;

	for ( ; event >= 0 ; i++, event = event_chains[gen].events[event].next){
		eventAge += event_chains[gen].events[event].elapsed_time;
		clade_stats[clade].sorted_ages[i] = eventAge;
		clade_stats[clade].elapsed_times[i] = event_chains[gen].events[event].elapsed_time;
		clade_stats[clade].num_lineages[i] = event_chains[gen].events[event].num_lineages;
		clade_stats[clade].event_types[i] = event_chains[gen].events[event].type;
	}
	clade_stats[clade].num_events = i;
	clade_stats[clade].coal_stats_total += genetree_stats[gen].coal_stats[clade];

#ifdef CHECKCLADE
	test_compareGphocsVsCladePopCoalStats(clade, gen, startingPoint);
#endif
}


void fillupCladeStats(int clade, int gen){
	addChildenIntoCladeStats(clade, gen);
	addCurrentPopIntoCladeStats(clade, gen);
}

void addChildenIntoCladeStats(int clade, int gen){
	mergeChildern(clade, gen);
	addChildrenCladeStats(clade, gen);
}

void mergeChildern(int clade, int gen){
	int i, j = 0, k = 0;
	double leftAge, rightAge;
	int leftSon, rightSon;
	leftSon = dataSetup.popTree->pops[clade]->sons[LEFT]->id;
	rightSon = dataSetup.popTree->pops[clade]->sons[RIGHT]->id;
	int m = clade_stats[leftSon].num_events;
	int n = clade_stats[rightSon].num_events;

	for (i = 0 ; i < m + n; ) {
		if ( j < m && k < n) {
			leftAge = clade_stats[leftSon].sorted_ages[j];
			rightAge = clade_stats[rightSon].sorted_ages[k];
			if (leftAge < rightAge){
				clade_stats[clade].event_types[i] = clade_stats[leftSon].event_types[j];
				clade_stats[clade].sorted_ages[i] = clade_stats[leftSon].sorted_ages[j];
				clade_stats[clade].num_lineages[i] =
						clade_stats[leftSon].num_lineages[j] + clade_stats[rightSon].num_lineages[k];
				j++;
			}
			else{
				clade_stats[clade].event_types[i] = clade_stats[rightSon].event_types[k];
				clade_stats[clade].sorted_ages[i] = clade_stats[rightSon].sorted_ages[k];
				clade_stats[clade].num_lineages[i] =
						clade_stats[leftSon].num_lineages[j] + clade_stats[rightSon].num_lineages[k];
				k++;
			}
			i++;
		}
		else if (j == m) {
			for (; i < m + n ;) {
				clade_stats[clade].event_types[i] = clade_stats[rightSon].event_types[k];
				clade_stats[clade].sorted_ages[i] = clade_stats[rightSon].sorted_ages[k];
				clade_stats[clade].num_lineages[i] =
						clade_stats[leftSon].num_lineages[j-1] + clade_stats[rightSon].num_lineages[k];
				k++;
				i++;
			}
		}
		else {
			for (; i < m + n;) {
				clade_stats[clade].event_types[i] = clade_stats[leftSon].event_types[j];
				clade_stats[clade].sorted_ages[i] = clade_stats[leftSon].sorted_ages[j];
				clade_stats[clade].num_lineages[i] =
						clade_stats[leftSon].num_lineages[j] + clade_stats[rightSon].num_lineages[k-1];
				j++;
				i++;
			}
		}
	}
	for (i = 0 ; i < m + n; i++ ) {
		clade_stats[clade].elapsed_times[i] = (i==0) ? 0 :
				clade_stats[clade].sorted_ages[i] - clade_stats[clade].sorted_ages[i-1];
	}

	clade_stats[clade].num_events = m + n;
}

void addChildrenCladeStats(int clade, int gen){
	clade_stats[clade].coal_stats_total +=
			getCoalStats(clade_stats[clade].elapsed_times, clade_stats[clade].num_lineages, clade_stats[clade].num_events);
}

void addCurrentPopIntoCladeStats(int clade, int gen){
	appendPopToClade(clade, gen, clade_stats[clade].num_events); // start filling the clade_stats arrays from the last known event
}















/*	computeGenetreeStats
	Computes the statistics of a given gene tree. Assumes event chains are built, 
	but number of lineages is ONLY set for first events in the leaf populations.
	Sets number of lineages for each non-leaf event by traversing the 
	population tree post-order. In parallel, also records the statistics.
 */
int computeGenetreeStats(int gen) {

	int i, pop, pop_queue[2*NSPECIES-1];		// post-order queue of populations

	// go over all chains and compute num_lineages per each event
	// also update genetree statistics
	populationPostOrder(dataSetup.popTree->rootPop, pop_queue);

	for(i=0; i<dataSetup.popTree->numPops; i++) {
		pop = pop_queue[i];
		// if not leaf population get number of in-lineages from end-events of son populations.
		if(pop >= dataSetup.popTree->numCurPops) {
			event_chains[gen].events[ event_chains[gen].first_event[pop] ].num_lineages =
					event_chains[gen].events[ dataSetup.popTree->pops[pop]->sons[0]->id ].num_lineages +
					event_chains[gen].events[ dataSetup.popTree->pops[pop]->sons[1]->id ].num_lineages;
		} else {
			event_chains[gen].events[ event_chains[gen].first_event[pop] ].num_lineages = 0;
		}

		recalcStats(gen,pop);
	}
	return 0;
}


/* recalcStats
   Re-calculates stats for given population in given gen.
   Writes down stats in genetree_stats_check and then compares to prior stats to return the log-likelihood.
   The log-likelihood computation ignores changes in number of migs/coals!!
 */

double recalcStats(int gen, int pop) {
	int n, id, mig_band, event;
	int live_mig_bands[MAX_MIG_BANDS];
	int num_live_mig_bands = 0;

	double t, heredity_factor = 1;
	double delta_lnLd = 0.0;

	genetree_stats_check.coal_stats[pop] = 0.0;
	genetree_stats_check.num_coals[pop]  = 0;
	event = event_chains[gen].first_event[pop];
	n  = event_chains[gen].events[event].num_lineages;


	// follow event chain and set number of lineages per interval according to previous event
	// also update statistics
	for( ; event>=0; event = event_chains[gen].events[event].next) {

		event_chains[gen].events[event].num_lineages = n;
		id = event_chains[gen].events[event].node_id;
		t  = event_chains[gen].events[event].elapsed_time;

		genetree_stats_check.coal_stats[pop] += n*(n-1)*t;

		for(mig_band=0; mig_band<num_live_mig_bands; mig_band++) {
			genetree_stats_check.mig_stats[live_mig_bands[mig_band]] += n*t;
		}

		switch(event_chains[gen].events[event].type) {
			case(SAMPLES_START):
																																n += dataSetup.numSamplesPerPop[pop];
			break;
			case(COAL):
																																genetree_stats_check.num_coals[pop]++;
			n--;
			break;
			case(IN_MIG):
    																															// figure out migration band and update its statistics
    																															mig_band = genetree_migs[gen].mignodes[id].migration_band;
			genetree_stats_check.num_migs[mig_band]++;
			n--;
			break;
			case(OUT_MIG):
																																n++;
			break;
			case(MIG_BAND_START):
    																															live_mig_bands[num_live_mig_bands++] = id;
			// initialize statistics for this new migration band
			genetree_stats_check.num_migs[id] = 0;
			genetree_stats_check.mig_stats[id] = 0.0;
			break;
			case(MIG_BAND_END):
							// compare and copy stats for mig band
							delta_lnLd -= (genetree_stats_check.mig_stats[id] - genetree_stats[gen].mig_stats[id]) * dataSetup.popTree->migBands[id].migRate;
			genetree_stats_total.mig_stats[id] += genetree_stats_check.mig_stats[id] - genetree_stats[gen].mig_stats[id];
			genetree_stats_total.num_migs[id]  += genetree_stats_check.num_migs[id] - genetree_stats[gen].num_migs[id];
			genetree_stats[gen].mig_stats[id] = genetree_stats_check.mig_stats[id];
			genetree_stats[gen].num_migs[id] = genetree_stats_check.num_migs[id];
			// remove mig band from living list
			for(mig_band=0; mig_band<num_live_mig_bands; mig_band++) {
				if(live_mig_bands[mig_band] == id)
					break;
			}
			if(mig_band == num_live_mig_bands) {
				if(debug) {
					fprintf(stderr, "\nError: recalcStats: migration band %d not alive in population %d, gen %d.\n",
							event_chains[gen].events[event].node_id, pop, gen);
					fprintf(stderr, "Live migration bands are:");
					for(mig_band=0; mig_band<num_live_mig_bands; mig_band++) {
						fprintf(stderr, " %d",live_mig_bands[mig_band]);
					}
					fprintf(stderr, "\n");
				} else {
					fprintf(stderr, "Fatal Error 0025.\n");
				}
				printGenealogyAndExit(gen,-1);
			}
			live_mig_bands[mig_band] = live_mig_bands[--num_live_mig_bands];
			break;
			case(DUMMY):
			case(END_CHAIN):
			break;
			default:
				if(debug) {
					fprintf(stderr, "\nError: recalcStats: event of unknown type %d in population %d, gen %d.\n",
							event_chains[gen].events[event].type, pop, gen);
				} else {
					fprintf(stderr, "Fatal Error 0026.\n");
				}
				printGenealogyAndExit(gen,-1);
				break;
		}// end of switch

	}// end of for(event)

	if(num_live_mig_bands != 0) {
		if(debug) {
			fprintf(stderr, "\nError: recalcStats: number of live mig bands %d at end of population %d in gen %d.\n",
					num_live_mig_bands, pop, gen);
		} else {
			fprintf(stderr, "Fatal Error 0027.\n");
		}
		printGenealogyAndExit(gen,-1);
	}

	delta_lnLd -= (genetree_stats_check.coal_stats[pop] - genetree_stats[gen].coal_stats[pop]) /
			( dataSetup.popTree->pops[pop]->theta * heredity_factor);
	genetree_stats_total.coal_stats[pop] += (genetree_stats_check.coal_stats[pop] - genetree_stats[gen].coal_stats[pop])/heredity_factor;
	genetree_stats_total.num_coals[pop]  += genetree_stats_check.num_coals[pop] - genetree_stats[gen].num_coals[pop];
	genetree_stats[gen].coal_stats[pop] = genetree_stats_check.coal_stats[pop];
	genetree_stats[gen].num_coals[pop] = genetree_stats_check.num_coals[pop];


	return delta_lnLd;
}
/*** end of recalcStats ***/


/* recalcStats_partitioned
   Re-calculates stats for given population in given gen - partitions population into K equally sized chunks
   Writes down stats in genetree_stats_total_partitioned.
  returns 0   
 */

int recalcStats_partitioned(int gen, int pop) {
	int partition, numPartitions;
	double age, partitionTime, partitionEndTime;
	int n, event;
	//    --  do not record migration stats for now
	// int live_mig_bands[MAX_MIG_BANDS];
	// int mig_band, num_live_mig_bands = 0;

	double t;
	double PERCISION = 0.0000001;

	// do not initialize so that stats are accumulated
	// genetree_stats_check.coal_stats[pop] = 0.0;
	// genetree_stats_check.num_coals[pop]  = 0;
	event = event_chains[gen].first_event[pop];
	n  = event_chains[gen].events[event].num_lineages;
	age = dataSetup.popTree->pops[pop]->age;

	if(pop == dataSetup.popTree->rootPop) {
		numPartitions = 1;
		partitionTime = OLDAGE - age;
	} else {
		numPartitions = dataSetup.numPopPartitions;
		partitionTime = (dataSetup.popTree->pops[pop]->father->age - age) / numPartitions;
	}
	partition = 0;
	partitionEndTime = age + partitionTime;

	// follow event chain and set number of lineages per interval according to previous event
	// also update statistics
	for( ; event>=0; event = event_chains[gen].events[event].next) {

		event_chains[gen].events[event].num_lineages = n;
		t  = event_chains[gen].events[event].elapsed_time;
		age += t;

		// check to see if crosses to next partition
		while(age > partitionEndTime) {

			// calculate contribution for previous partition
			t -= (age-partitionEndTime);
			genetree_stats_total_partitioned[partition].coal_stats[pop] += n*(n-1)*t;

			// move on to next partition
			partition++;
			t = age - partitionEndTime;
			partitionEndTime += partitionTime;

			// check that partitions were not exhausted
			if(partition >= numPartitions) {
				if(t<=PERCISION) {
					t=0.0;
					partition--;
					partitionEndTime -= partitionTime;
					break;
				} else {
					if(debug) {
						fprintf(stderr, "\nError: recalcStats_partitioned: overflow time of %g in pop %d, gen %d.\n", t, pop, gen);
					} else {
						fprintf(stderr, "Fatal Error 9001.\n");
					}
					printGenealogyAndExit(gen,-1);
					break;
				}
			}
		}// end of while(age > partitionEndTime)

		genetree_stats_total_partitioned[partition].coal_stats[pop] += n*(n-1)*t;

		// do not record migration stats for now
		//
		// for(mig_band=0; mig_band<num_live_mig_bands; mig_band++) {
		//   genetree_stats_total_partitioned[partition].mig_stats[live_mig_bands[mig_band]] += n*t;
		// }

		switch(event_chains[gen].events[event].type) {
			case(SAMPLES_START):
						n += dataSetup.numSamplesPerPop[pop];
			break;
			case(COAL):
						genetree_stats_total_partitioned[partition].num_coals[pop]++;
			n--;
			break;
			case(IN_MIG):
						  // figure out migration band and update its statistics  --  do not record migration stats for now
						  //
						  // mig_band = genetree_migs[gen].mignodes[id].migration_band;
						  // genetree_stats_total_partitioned[partition].num_migs[mig_band]++;
						  n--;
			break;
			case(OUT_MIG):
						  n++;
			break;
			case(MIG_BAND_START):
						  // initialize statistics for this new migration band    --  do not record migration stats for now
						  //
						  // live_mig_bands[num_live_mig_bands++] = id;
						  // genetree_stats_total_partitioned[partition].num_migs[id] = 0;
						  // genetree_stats_total_partitioned[partition].mig_stats[id] = 0.0;
						  break;
			case(MIG_BAND_END):
						  // compare and copy stats for mig band    --   DO NOT DO FOR NOW !!
						  /*
delta_lnLd -= (genetree_stats_check.mig_stats[id] - genetree_stats[gen].mig_stats[id]) * dataSetup.popTree->migBands[id].migRate;
genetree_stats_total.mig_stats[id] += genetree_stats_check.mig_stats[id] - genetree_stats[gen].mig_stats[id];
genetree_stats_total.num_migs[id]  += genetree_stats_check.num_migs[id] - genetree_stats[gen].num_migs[id];
genetree_stats[gen].mig_stats[id] = genetree_stats_check.mig_stats[id];
genetree_stats[gen].num_migs[id] = genetree_stats_check.num_migs[id];
//remove mig band from living list
for(mig_band=0; mig_band<num_live_mig_bands; mig_band++) {
if(live_mig_bands[mig_band] == id)
break;
}
if(mig_band == num_live_mig_bands) {
if(debug) {
fprintf(stderr, "\nError: recalcStats: migration band %d not alive in population %d, gen %d.\n",
event_chains[gen].events[event].node_id, pop, gen);
fprintf(stderr, "Live migration bands are:");
for(mig_band=0; mig_band<num_live_mig_bands; mig_band++) {
fprintf(stderr, " %d",live_mig_bands[mig_band]);
}
fprintf(stderr, "\n");
} else {
fprintf(stderr, "Fatal Error 0025.\n");
}
printGenealogyAndExit(gen,-1);
}
live_mig_bands[mig_band] = live_mig_bands[--num_live_mig_bands];
						   */
						  break;
			case(DUMMY):
			case(END_CHAIN):
			break;
			default:
				if(debug) {
					fprintf(stderr, "\nError: recalcStats_partitioned: event of unknown type %d in population %d, gen %d.\n",
							event_chains[gen].events[event].type, pop, gen);
				} else {
					fprintf(stderr, "Fatal Error 9026.\n");
				}
				printGenealogyAndExit(gen,-1);
				break;
		}// end of switch

	}// end of for(event)
	if(partition != numPartitions-1 || partitionEndTime - age > PERCISION) {
		if(debug) {
			fprintf(stderr, "\nError: recalcStats: excess of time left, partition %d of %d, time %g of %g, at  population %d in gen %d.\n",
					partition+1, numPartitions, age, partitionEndTime, pop, gen);
		} else {
			fprintf(stderr, "Fatal Error 9027.\n");
		}
		printGenealogyAndExit(gen,-1);
	}

	/*   --  no migration stuff  !!
  if(num_live_mig_bands != 0) {
		if(debug) {
			fprintf(stderr, "\nError: recalcStats: number of live mig bands %d at end of population %d in gen %d.\n",
				num_live_mig_bands, pop, gen);
		} else {
			fprintf(stderr, "Fatal Error 9027.\n");
		}
    printGenealogyAndExit(gen,-1);
  }
	 */

	// updateing global structures   --   do we need to do this here ??

	return 0;
}
/*** end of recalcStats_partitioned ***/


/*	gtreeLnLikelihood
	Computes likelihood of genetree (with migration events) of a given gen.
	Assumes statistics are already computed.
 */

double gtreeLnLikelihood(int gen) {

	int pop, mig_band, sample, node;
	double lnLd, mig_rate, theta, heredity_factor=1;

	lnLd = 0;

	for(pop=0; pop < dataSetup.popTree->numPops; pop++) {
		theta = dataSetup.popTree->pops[pop]->theta * heredity_factor;
		lnLd +=
				genetree_stats[gen].num_coals[pop]  * log(2/theta)	-
				genetree_stats[gen].coal_stats[pop] / (theta);
	}

	for(mig_band=0; mig_band<dataSetup.popTree->numMigBands; mig_band++) {
		mig_rate = dataSetup.popTree->migBands[mig_band].migRate;
		if(mig_rate > 0.0) {
			lnLd +=
					genetree_stats[gen].num_migs[mig_band]  * log(mig_rate)	-
					genetree_stats[gen].mig_stats[mig_band] * mig_rate;
		}
	}

	// considering admixture
	for(sample=0; sample<admixed_samples.number; sample++) {
		node = admixed_samples.samples[sample];
		if(nodePops[gen][node] == admixed_samples.popPairs[sample][0]) {
			lnLd += log(1-admixture_status.admixtureCoefficients[sample]);
		} else {
			lnLd += log(admixture_status.admixtureCoefficients[sample]);
		}
	}// end of for(sample)


	return lnLd;

}


/* checkAll()
   invokes checkGtreeStructure on all loci and checks genetree_stats_total
   also checks data log likelihood
 */
int checkAll() {
	int gen, mig_band, pop, heredity_factor, res = 1;
	double PERCISION = 0.0000001;
	double lnLd_gen, genLnLd, dataLnLd;

	for(pop=0; pop<dataSetup.popTree->numPops; pop++) {
		genetree_stats_total_check.num_coals[pop]  = 0;
		genetree_stats_total_check.coal_stats[pop] = 0.0;
	}
	for(mig_band=0; mig_band<dataSetup.popTree->numMigBands; mig_band++) {
		genetree_stats_total_check.num_migs[mig_band]  = 0;
		genetree_stats_total_check.mig_stats[mig_band] = 0.0;
	}

	genLnLd = 0.0;
	dataLnLd = 0.0;

	for(gen=0; gen<dataSetup.numLoci; gen++) {
		if(!checkGtreeStructure(gen)) {
			fprintf(stderr, "\nError: Checking gene tree structure failed\n");
			printGenealogyAndExit(gen,0);
			return 0;
		}


		if(!checkLocusDataLikelihood(dataState.lociData[gen])) {
			fprintf(stderr, "\nError: checking recorded likelihood for gen %d!", gen);
			printGenealogyAndExit(gen,0);
			return 0;
		}

		lnLd_gen = gtreeLnLikelihood(gen);

		if(fabs(genLogLikelihood[gen] - lnLd_gen) > PERCISION && fabs(1-genLogLikelihood[gen]/lnLd_gen) > PERCISION) {
			fprintf(stderr, "\nError: in recorded genealogy log likelihood for gen %d: (recorded = %g, actual = %g, diff = %g, 1-ratio = %g).\n",
					gen, genLogLikelihood[gen], lnLd_gen, genLogLikelihood[gen] - lnLd_gen, 1-genLogLikelihood[gen]/lnLd_gen);
			printGenealogyAndExit(gen,0);
			return 0;
		}
		genLogLikelihood[gen] = lnLd_gen;
		dataLnLd += getLocusDataLikelihood(dataState.lociData[gen]);
		genLnLd  += lnLd_gen;

		heredity_factor = 1;
		for(pop=0; pop<dataSetup.popTree->numPops; pop++) {
			genetree_stats_total_check.num_coals[pop]  += genetree_stats[gen].num_coals[pop];
			genetree_stats_total_check.coal_stats[pop] += genetree_stats[gen].coal_stats[pop]/heredity_factor;
		}
		for(mig_band=0; mig_band<dataSetup.popTree->numMigBands; mig_band++) {
			genetree_stats_total_check.num_migs[mig_band]  += genetree_stats[gen].num_migs[mig_band];
			genetree_stats_total_check.mig_stats[mig_band] += genetree_stats[gen].mig_stats[mig_band];
		}


	}

	for(pop=0; pop<dataSetup.popTree->numPops; pop++) {
		if(fabs(genetree_stats_total_check.coal_stats[pop] - genetree_stats_total.coal_stats[pop]) > PERCISION && fabs(1-genetree_stats_total_check.coal_stats[pop]/genetree_stats_total.coal_stats[pop]) > PERCISION) {
			if(debug) {
				fprintf(stderr, "\nError: checking total coal stats for pop %d: %f (saved) %f (actual) %g (difference) %g (1-ratio)",
						pop, genetree_stats_total.coal_stats[pop], genetree_stats_total_check.coal_stats[pop],
						genetree_stats_total_check.coal_stats[pop] - genetree_stats_total.coal_stats[pop],
						1 - genetree_stats_total_check.coal_stats[pop] / genetree_stats_total.coal_stats[pop]);
			} else {
				fprintf(stderr, "Fatal Error 0028.\n");
			}
			res = 0;
		}
		genetree_stats_total.coal_stats[pop] = genetree_stats_total_check.coal_stats[pop];
		if(genetree_stats_total_check.num_coals[pop] != genetree_stats_total.num_coals[pop]) {
			if(debug) {
				fprintf(stderr, "\nError: checking total num coals for pop %d: %d (saved) %d (actual)",
						pop, genetree_stats_total.num_coals[pop], genetree_stats_total_check.num_coals[pop]);
			} else {
				fprintf(stderr, "Fatal Error 0029.\n");
			}
			res = 0;
		}

		genetree_stats_total.coal_stats[pop] = genetree_stats_total_check.coal_stats[pop];
	}
	for(mig_band=0; mig_band<dataSetup.popTree->numMigBands; mig_band++) {
		if(fabs(genetree_stats_total_check.mig_stats[mig_band] - genetree_stats_total.mig_stats[mig_band]) > PERCISION && fabs(1-genetree_stats_total_check.mig_stats[mig_band]/genetree_stats_total.mig_stats[mig_band]) > PERCISION) {
			if(debug) {
				fprintf(stderr, "\nError: checking total mig stats for mig band %d: %f (saved) %f (actual), %f (diff) %f (1-ratio)",
						mig_band, genetree_stats_total.mig_stats[mig_band], genetree_stats_check.mig_stats[mig_band],
						genetree_stats_total_check.mig_stats[mig_band] - genetree_stats_total.mig_stats[mig_band],
						1-genetree_stats_total_check.mig_stats[mig_band] / genetree_stats_total.mig_stats[mig_band]);
			} else {
				fprintf(stderr, "Fatal Error 0030.\n");
			}
			res = 0;
		}
		genetree_stats_total.mig_stats[mig_band] = genetree_stats_total_check.mig_stats[mig_band];
		if(genetree_stats_total_check.num_migs[mig_band] != genetree_stats_total.num_migs[mig_band]) {
			if(debug) {
				fprintf(stderr, "\nError: checking total num migs for mig band %d: %d (saved) %d (actual)",
						mig_band, genetree_stats_total.num_migs[mig_band], genetree_stats_total_check.num_migs[mig_band]);
			} else {
				fprintf(stderr, "Fatal Error 0031.\n");
			}
			res = 0;
		}

		genetree_stats_total.mig_stats[mig_band] = genetree_stats_total_check.mig_stats[mig_band];
	}
	if(fabs(dataState.dataLogLikelihood - dataLnLd) > PERCISION && fabs(1-dataState.dataLogLikelihood/dataLnLd) > PERCISION) {
		if(debug) {
			fprintf(stderr, "\nError: in recorded data log likelihood: (recorded = %g, actual = %g, diff = %g, 1-ratio = %g).\n", 
					dataState.dataLogLikelihood, dataLnLd, dataState.dataLogLikelihood - dataLnLd, 1-dataState.dataLogLikelihood/dataLnLd);
		} else {
			fprintf(stderr, "Fatal Error 0032.\n");
		}
		res = 0;
	}
	dataState.dataLogLikelihood = dataLnLd;
	// check recorded total log-likelihood
	dataLnLd = (genLnLd+dataLnLd)/dataSetup.numLoci;
	if(fabs(1-dataState.logLikelihood/dataLnLd) > PERCISION) {
		if(debug) {
			fprintf(stderr, "\nError: in recorded total log likelihood: (recorded = %g, actual = %g).\n", dataState.logLikelihood, dataLnLd);
		} else {
			fprintf(stderr, "Fatal Error 0033.\n");
		}
		res = 0;
	}
	dataState.logLikelihood = dataLnLd;

	if(!res){
		printf("\n\n");
		//		exit(-1);
	}

	return res;

}


/* checkGtreeStructure
   Checks consistency of genetree data - nodes, mignodes, events and stats
 */

int checkGtreeStructure(int gen) {

	int i, n, pop, mig_band, event, id, node, mig, num_migs;
	int num_living_mig_bands, living_mig_bands[MAX_MIG_BANDS], pop_queue[2*NSPECIES-1];
	int pop_lins_in[2*NSPECIES-1];
	int father;
	int res = 1;

	double age, delta_t, PERCISION = 0.0000000001;


	// essentially follow the same path as procedure computeGenetreeStats,
	// but validates with genetree nodes

	// initialize number of in-lineages for leaf pops and order pops
	for(pop=0; pop<dataSetup.popTree->numPops; pop++) {
		pop_lins_in[pop] = 0;
	}
	// MARK CHANGHES

	// for(node=0; node<2*dataSetup.numSamples-1; node++) {
	//  if(node < dataSetup.numSamples)
	//     pop_lins_in[ nodePops[gen][node] ]++;
	// }
	populationPostOrder(dataSetup.popTree->rootPop, pop_queue);

	for(i=0; i<dataSetup.popTree->numPops; i++) {
		pop = pop_queue[i];
		genetree_stats_check.coal_stats[pop] = 0.0;
		genetree_stats_check.num_coals[pop]  = 0;
		event = event_chains[gen].first_event[pop];
		n  = pop_lins_in[pop];
		age = dataSetup.popTree->pops[pop]->age;
		num_living_mig_bands = 0;


		// follow event chain and check all saved data
		for( ; event>=0; event = event_chains[gen].events[event].next) {


			if(event_chains[gen].events[event].num_lineages != n) {
				if(debug) {
					fprintf(stderr, "\nError: checking genetree for gen %d: ",gen);
					fprintf(stderr, "number of lineages of event %d doesn't match: %d, %d.\n", event, event_chains[gen].events[event].num_lineages, n);
				} else {
					fprintf(stderr, "Fatal Error 0034.\n");
				}
				res = 0;
			}

			if(event_chains[gen].events[event].next >=0 && event != event_chains[gen].events[event_chains[gen].events[event].next].prev) {
				if(debug) {
					fprintf(stderr, "\nError: checking genetree for gen %d: ",gen);
					fprintf(stderr, "prev event of next event %d doesn't match: %d, %d.\n",
							event_chains[gen].events[event].next, event, event_chains[gen].events[event_chains[gen].events[event].next].prev);
				} else {
					fprintf(stderr, "Fatal Error 0035.\n");
				}
				res = 0;
			}
			id = event_chains[gen].events[event].node_id;
			delta_t  = event_chains[gen].events[event].elapsed_time;
			age += delta_t;
			//			if(delta_t <= PERCISION) {
			//				printf("\nError checking genetree for gen %d: ",gen);
			//				printf("Too small time interval for event %d: %g.", event, delta_t);
			//				res = 0;
			//			}

			genetree_stats_check.coal_stats[pop] += n*(n-1)*delta_t;

			for(mig_band=0; mig_band<num_living_mig_bands; mig_band++) {
				genetree_stats_check.mig_stats[living_mig_bands[mig_band]] += n*delta_t;
			}

			switch(event_chains[gen].events[event].type) {
				case(SAMPLES_START):
																																																n += dataSetup.numSamplesPerPop[pop];
				if(fabs(dataSetup.popTree->pops[pop]->sampleAge - age) > PERCISION) {
					//					if(fabs(gnodes[gen][id].age - age) > PERCISION) {
					if(debug) {
						fprintf(stderr, "\nError: checking genetree for gen %d: ",gen);
						fprintf(stderr, "sample age of pop %d (event %d) doesn't match: %f, %f.\n", pop, event, age, dataSetup.popTree->pops[pop]->sampleAge);
					} else {
						fprintf(stderr, "Fatal Error 0036.\n");
					}
					res = 0;
				}
				break;
				case(COAL):
        																																														genetree_stats_check.num_coals[pop]++;
				n--;
				if(fabs(getNodeAge(dataState.lociData[gen],id) - age) > PERCISION) {
					//					if(fabs(gnodes[gen][id].age - age) > PERCISION) {
					if(debug) {
						fprintf(stderr, "\nError: checking genetree for gen %d: ",gen);
						fprintf(stderr, "age of node %d (event %d) doesn't match: %g, %g.\n", id, event, age, getNodeAge(dataState.lociData[gen],id));
					} else {
						fprintf(stderr, "Fatal Error 0036.\n");
					}
					res = 0;
				}

				if(nodePops[gen][id] != pop) {
					if(debug) {
						fprintf(stderr, "\nError: checking genetree for gen %d: ",gen);
						fprintf(stderr, "population of node %d (event %d) doesn't match: %d, %d.", id, event, pop, nodePops[gen][id]);
					} else {
						fprintf(stderr, "Fatal Error 0037.\n");
					}
					res = 0;
				}
				if(nodeEvents[gen][id] != event) {
					if(debug) {
						fprintf(stderr, "\nError: checking genetree for gen %d: ",gen);
						fprintf(stderr, "event of node %d doesn't match: %d, %d.", id, event, nodeEvents[gen][id]);
					} else {
						fprintf(stderr, "Fatal Error 0038.\n");
					}
					res = 0;
				}

				break;
				case(IN_MIG):
        																																														// figure out migration band and update its statistics
        																																														mig_band = genetree_migs[gen].mignodes[id].migration_band;
				if(mig_band<0 || mig_band> dataSetup.popTree->numMigBands) {
					if(debug) {
						fprintf(stderr, "\nError: checking genetree for gen %d: ",gen);
						fprintf(stderr, "mig band %d invalid for mig node %d .", mig_band, id);
					} else {
						fprintf(stderr, "Fatal Error 0039a.\n");
					}
					res = 0;
				}
				genetree_stats_check.num_migs[mig_band]++;
				n--;
				if(fabs(genetree_migs[gen].mignodes[id].age - age) > PERCISION) {
					if(debug) {
						fprintf(stderr, "\nError: checking genetree for gen %d: ",gen);
						fprintf(stderr, "age of mignode %d (target event %d) doesn't match: %g, %g.", id, event, age, genetree_migs[gen].mignodes[id].age);
					} else {
						fprintf(stderr, "Fatal Error 0039.\n");
					}
					res = 0;
				}
				if(dataSetup.popTree->migBands[mig_band].targetPop != pop) {
					if(debug) {
						fprintf(stderr, "\nError: checking genetree for gen %d: ",gen);
						fprintf(stderr, "target population of migration band %d of mignode %d (target event %d) doesn't match: %d, %d.",
								mig_band, id, event, pop, dataSetup.popTree->migBands[mig_band].targetPop);
					} else {
						fprintf(stderr, "Fatal Error 0040.\n");
					}
					res = 0;
				}
				break;
				case(OUT_MIG):
        																																														n++;
				mig_band = genetree_migs[gen].mignodes[id].migration_band;
				if(fabs(genetree_migs[gen].mignodes[id].age - age) > PERCISION) {
					if(debug) {
						fprintf(stderr, "\nError: checking genetree for gen %d: ",gen);
						fprintf(stderr, "age of mignode %d (source event %d) doesn't match: %g, %g.", id, event, age, genetree_migs[gen].mignodes[id].age);
					} else {
						fprintf(stderr, "Fatal Error 0041.\n");
					}
					res = 0;
				}
				if(dataSetup.popTree->migBands[mig_band].sourcePop != pop) {
					if(debug) {
						fprintf(stderr, "\nError: checking genetree for gen %d: ",gen);
						fprintf(stderr, "source population of migration band %d of mignode %d (source event %d) doesn't match: %d, %d.",
								mig_band, id, event, pop, dataSetup.popTree->migBands[mig_band].sourcePop);
					} else {
						fprintf(stderr, "Fatal Error 0042.\n");
					}
					res = 0;
				}
				if(genetree_migs[gen].mignodes[id].source_event != event) {
					if(debug) {
						fprintf(stderr, "\nError: checking genetree for gen %d: ",gen);
						fprintf(stderr, "source event of node %d doesn't match: %d, %d.", id, event, genetree_migs[gen].mignodes[id].source_event);
					} else {
						fprintf(stderr, "Fatal Error 0043.\n");
					}
					res = 0;
				}
				break;
				case(MIG_BAND_START):
        																																														living_mig_bands[num_living_mig_bands] = id;
				num_living_mig_bands++;
				// initialize statistics for this new migration band
				genetree_stats_check.num_migs[id] = 0;
				genetree_stats_check.mig_stats[id] = 0.0;
				if(fabs(dataSetup.popTree->migBands[id].startTime - age) > PERCISION) {
					if(debug) {
						fprintf(stderr, "\nError: checking genetree for gen %d: ",gen);
						fprintf(stderr, "start time of migration band %d (start event %d) doesn't match: %g, %g.", id, event, age, dataSetup.popTree->migBands[id].startTime);
					} else {
						fprintf(stderr, "Fatal Error 0044.\n");
					}
					res = 0;
				}
				if(dataSetup.popTree->migBands[id].targetPop != pop) {
					if(debug) {
						fprintf(stderr, "\nError: checking genetree for gen %d: ",gen);
						fprintf(stderr, "target population of migration band %d (start event %d) doesn't match: %d, %d.", id, event, pop, dataSetup.popTree->migBands[id].targetPop);
					} else {
						fprintf(stderr, "Fatal Error 0045.\n");
					}
					res = 0;
				}
				if(fabs(age - max2(dataSetup.popTree->pops[pop]->age, dataSetup.popTree->pops[dataSetup.popTree->migBands[id].sourcePop]->age)) > PERCISION) {
					if(debug) {
						fprintf(stderr, "\nError: checking genetree for gen %d: ",gen);
						fprintf(stderr, "start time of migration band %d (start event %d) doesn't match max of start times for populations %d,%d: %g, max(%g,%g).",
								id, event, pop, dataSetup.popTree->migBands[id].sourcePop, age, dataSetup.popTree->pops[pop]->age, dataSetup.popTree->pops[dataSetup.popTree->migBands[id].sourcePop]->age);
					} else {
						fprintf(stderr, "Fatal Error 0046.\n");
					}
					res = 0;
				}
				break;
				case(MIG_BAND_END):
        																																														for(mig_band=0; mig_band<num_living_mig_bands; mig_band++) {
        																																															if(living_mig_bands[mig_band] == id)
        																																																break;
        																																														}
				if(mig_band == num_living_mig_bands) {
					if(debug) {
						fprintf(stderr, "\nError: checking genetree for gen %d: ",gen);
						fprintf(stderr, "migration band %d (start event %d) ended without starting.", id, event);
					} else {
						fprintf(stderr, "Fatal Error 0047.\n");
					}
					res = 0;
				} else{
					living_mig_bands[mig_band] = living_mig_bands[--num_living_mig_bands];
				}
				if(fabs(dataSetup.popTree->migBands[id].endTime - age) > PERCISION) {
					if(debug) {
						fprintf(stderr, "\nError: checking genetree for gen %d: ",gen);
						fprintf(stderr, "end time of migration band %d (end event %d) doesn't match: %g, %g.", id, event, age, dataSetup.popTree->migBands[id].endTime);
					} else {
						fprintf(stderr, "Fatal Error 0048.\n");
					}
					res = 0;
				}
				if(dataSetup.popTree->migBands[id].targetPop != pop) {
					if(debug) {
						fprintf(stderr, "\nError: checking genetree for gen %d: ",gen);
						fprintf(stderr, "target population of migration band %d (end event %d) doesn't match: %d, %d.", id, event, pop, dataSetup.popTree->migBands[id].targetPop);
					} else {
						fprintf(stderr, "Fatal Error 0049.\n");
					}
					res = 0;
				}
				if(fabs( age - min2(	dataSetup.popTree->pops[pop]->father->age,
						dataSetup.popTree->pops[dataSetup.popTree->migBands[id].sourcePop]->father->age)   ) > PERCISION) {
					if(debug) {
						fprintf(stderr, "\nError: checking genetree for gen %d: ",gen);
						fprintf(stderr, "end time of migration band %d (end event %d) doesn't match min of start times for parent populations %d,%d: %g, min(%g,%g).",
								id, event, dataSetup.popTree->pops[pop]->father->id,
								dataSetup.popTree->pops[dataSetup.popTree->migBands[id].sourcePop]->father->id,
								age, dataSetup.popTree->pops[pop]->father->age,
								dataSetup.popTree->pops[dataSetup.popTree->migBands[id].sourcePop]->father->age);
					} else {
						fprintf(stderr, "Fatal Error 0050.\n");
					}
					res = 0;
				}

				break;
				case(END_CHAIN):
        																																														if(id != pop) {
        																																															if(debug) {
        																																																fprintf(stderr, "\nError: checking genetree for gen %d: ",gen);
        																																																fprintf(stderr, "id of end event %d for population %d doesn't match: %d.", event, pop, id);
        																																															} else {
        																																																fprintf(stderr, "Fatal Error 0051.\n");
        																																															}
        																																															res = 0;
        																																														}
				if(num_living_mig_bands != 0) {
					if(debug) {
						fprintf(stderr, "\nError: checking genetree for gen %d: ",gen);
						fprintf(stderr, "there are %d living migration bands remaining in population %d:", num_living_mig_bands,pop);
						for(i=0; i<num_living_mig_bands; i++) {
							fprintf(stderr," %d",living_mig_bands[i]);
						}
					} else {
						fprintf(stderr, "Fatal Error 0052.\n");
					}
					res = 0;
				}
				if(event_chains[gen].events[event].next >= 0) {
					if(debug) {
						fprintf(stderr, "\nError: checking genetree for gen %d: ",gen);
						fprintf(stderr, "remaining event %d after end event %d for pop %d.", event_chains[gen].events[event].next, event, pop);
					} else {
						fprintf(stderr, "Fatal Error 0053.\n");
					}
					res = 0;
				}

				//					if(pop == dataSetup.popTree->rootPop) {
				//						if(fabs(OLDAGE - age) > PERCISION) {
				//							printf("\nError checking genetree for gen %d: ",gen);
				//							printf("end time of population %d (end event %d) doesn't match OLDAGE: %f, %f.", pop, event, age, OLDAGE);
				//							res = 0;
				//						}
				//					}
				if(pop != dataSetup.popTree->rootPop) {
					// add to incoming lineages of parent population
					pop_lins_in[dataSetup.popTree->pops[pop]->father->id] += n;
					if(fabs(dataSetup.popTree->pops[pop]->father->age - age) > PERCISION) {
						if(debug) {
							fprintf(stderr, "\nError: checking genetree for gen %d: ",gen);
							fprintf(stderr, "end time of population %d (end event %d) doesn't match start of father pop %d: %f, %f.",
									pop, event, dataSetup.popTree->pops[pop]->father->id, age, dataSetup.popTree->pops[pop]->father->age);
						} else {
							fprintf(stderr, "Fatal Error 0054.\n");
						}
						res = 0;
					}
				}
				break;
				default:
					if(debug) {
						fprintf(stderr, "\nError: checking genetree for gen %d: ",gen);
						fprintf(stderr, "Event %d of illegal type: %d.", event, event_chains[gen].events[event].type);
					} else {
						fprintf(stderr, "Fatal Error 0055.\n");
					}
					res = 0;
					break;
			}// end of switch

		}// end of for(event)

		if(num_living_mig_bands != 0) {
			if(debug) {
				fprintf(stderr, "\nError: checking genetree for gen %d: number of live mig bands %d at end of population %d.\n",
						gen, num_living_mig_bands, pop);
			} else {
				fprintf(stderr, "Fatal Error 0056.\n");
			}
			printGenealogyAndExit(gen,-1);
		}

	}// end of for(pop)


	// check computed stats
	for(pop=0; pop<dataSetup.popTree->numPops; pop++) {
		if(fabs(genetree_stats_check.coal_stats[pop] - genetree_stats[gen].coal_stats[pop]) > PERCISION) {
			if(debug) {
				fprintf(stderr, "\nError: checking genetree for gen %d: ",gen);
				fprintf(stderr, "coal stats for pop %d don't match: %g, %g (diff = %g)", pop, genetree_stats[gen].coal_stats[pop], genetree_stats_check.coal_stats[pop],
						genetree_stats[gen].coal_stats[pop]- genetree_stats_check.coal_stats[pop]);
			} else {
				fprintf(stderr, "Fatal Error 0057.\n");
			}
			res = 0;
		}
		genetree_stats[gen].coal_stats[pop] = genetree_stats_check.coal_stats[pop];
		if(genetree_stats_check.num_coals[pop] != genetree_stats[gen].num_coals[pop]) {
			if(debug) {
				fprintf(stderr, "\nError: checking genetree for gen %d: ",gen);
				fprintf(stderr, "num coals for pop %d don't match: %d, %d.", pop, genetree_stats[gen].num_coals[pop], genetree_stats_check.num_coals[pop]);
			} else {
				fprintf(stderr, "Fatal Error 0058.\n");
			}
			res = 0;
		}
	}
	for(mig_band=0; mig_band<dataSetup.popTree->numMigBands; mig_band++) {
		if(fabs(genetree_stats_check.mig_stats[mig_band] - genetree_stats[gen].mig_stats[mig_band]) > PERCISION) {
			if(debug) {
				fprintf(stderr, "\nError: checking genetree for gen %d: ",gen);
				fprintf(stderr, "mig stats for migration band %d don't match: %f, %f.", mig_band, genetree_stats[gen].mig_stats[mig_band], genetree_stats_check.mig_stats[mig_band]);
			} else {
				fprintf(stderr, "Fatal Error 0059.\n");
			}
			res = 0;
		}
		genetree_stats[gen].mig_stats[mig_band] = genetree_stats_check.mig_stats[mig_band];
		if(genetree_stats_check.num_migs[mig_band] != genetree_stats[gen].num_migs[mig_band]) {
			if(debug) {
				fprintf(stderr, "\nError: checking genetree for gen %d: ",gen);
				fprintf(stderr, "num migs for migration band %d don't match: %d, %d.", mig_band, genetree_stats[gen].num_migs[mig_band], genetree_stats_check.num_migs[mig_band]);
			} else {
				fprintf(stderr, "Fatal Error 0060.\n");
			}
			res = 0;
		}
	}

	// check tree structure
	num_migs = 0;
	for(node=0; node<2*dataSetup.numSamples-1; node++) {

		pop = nodePops[gen][node];
		mig = findFirstMig(gen,node,-1);
		age = getNodeAge(dataState.lociData[gen],node);
		//age = gnodes[gen][node].age;
		while(mig >= 0) {
			num_migs++;
			if(genetree_migs[gen].mignodes[mig].gtree_branch != node) {
				if(debug) {
					fprintf(stderr, "\nError: checking genetree for gen %d: ",gen);
					fprintf(stderr, "migration node %d for node %d points to wrong edge: %d.", mig, node, genetree_migs[gen].mignodes[mig].gtree_branch);
				} else {
					fprintf(stderr, "Fatal Error 0061.\n");
				}
				res = 0;
			}
			if(genetree_migs[gen].mignodes[mig].age <= age) {
				if(debug) {
					fprintf(stderr, "\nError: checking genetree for gen %d: ",gen);
					fprintf(stderr, "migration node %d for node %d has conflicting age: %f, %f.", mig, node, genetree_migs[gen].mignodes[mig].age, age);
				} else {
					fprintf(stderr, "Fatal Error 0062.\n");
				}
				res = 0;
			}
			for(i=0; i<genetree_migs[gen].num_migs; i++) {
				if(genetree_migs[gen].living_mignodes[i] == mig)   break;
			}
			if(i >= genetree_migs[gen].num_migs) {
				if(debug) {
					fprintf(stderr, "\nError: checking genetree for gen %d: ",gen);
					fprintf(stderr, "mignode %d doesn't appear in list of living mignodes.", mig);
				} else {
					fprintf(stderr, "Fatal Error 0063.\n");
				}
				res = 0;
			}
			age = genetree_migs[gen].mignodes[mig].age;
			mig_band = genetree_migs[gen].mignodes[mig].migration_band;

			if(dataSetup.popTree->migBands[mig_band].sourcePop != genetree_migs[gen].mignodes[mig].source_pop) {
				if(debug) {
					fprintf(stderr, "\nError: checking genetree for gen %d: ",gen);
					fprintf(stderr, "source population of mignode %d doesn't match mig band %d: %d, %d.",
							mig, mig_band, genetree_migs[gen].mignodes[mig].source_pop, dataSetup.popTree->migBands[mig_band].sourcePop);
				} else {
					fprintf(stderr, "Fatal Error 0064.\n");
				}
				res = 0;
			}
			if(dataSetup.popTree->migBands[mig_band].targetPop != genetree_migs[gen].mignodes[mig].target_pop) {
				if(debug) {
					fprintf(stderr, "\nError: checking genetree for gen %d: ",gen);
					fprintf(stderr, "target population of mignode %d doesn't match mig band %d: %d, %d.",
							mig, mig_band, genetree_migs[gen].mignodes[mig].target_pop, dataSetup.popTree->migBands[mig_band].targetPop);
				} else {
					fprintf(stderr, "Fatal Error 0065.\n");
				}
				res = 0;
			}
			if(age <= dataSetup.popTree->migBands[mig_band].startTime || age >= dataSetup.popTree->migBands[mig_band].endTime) {
				if(debug) {
					fprintf(stderr, "\nError: checking genetree for gen %d: ",gen);
					fprintf(stderr, "age of mignode %d conflicts with times of migration band %d: %f, [%f,%f].",
							mig, mig_band, age, dataSetup.popTree->migBands[mig_band].startTime, dataSetup.popTree->migBands[mig_band].endTime);
				} else {
					fprintf(stderr, "Fatal Error 0066.\n");
				}
				res = 0;
			}
			if(!dataSetup.popTree->pops[ genetree_migs[gen].mignodes[mig].target_pop ]->isAncestralTo[pop]) {
				if(debug) {
					fprintf(stderr, "\nError: checking genetree for gen %d: ",gen);
					fprintf(stderr, "target pop %d of mignode %d is not ancestral to previous pop %d.", genetree_migs[gen].mignodes[mig].target_pop, mig, pop);
				} else {
					fprintf(stderr, "Fatal Error 0067.\n");
				}
				res = 0;
			}

			pop = genetree_migs[gen].mignodes[mig].source_pop;
			mig = findFirstMig(gen,node,age);
		}// end of while(mig>=0)

		father = getNodeFather(dataState.lociData[gen], node);
		if(father >= 0) {
			//		if(gnodes[gen][node].father >= 0) {
			/*			if(gnodes[gen][gnodes[gen][node].father].nson != 2) {
                    printf("\nError checking genetree for gen %d: ",gen);
                    printf("father %d of node %d has number of sons: %d.", gnodes[gen][node].father, node, gnodes[gen][gnodes[gen][node].father].nson);
                    res = 0;
                    } else
			 */
			if(getNodeSon(dataState.lociData[gen], father, 0) != node && getNodeSon(dataState.lociData[gen], father, 1) != node) {
				//			if(gnodes[gen][gnodes[gen][node].father].sons[0] != node  && gnodes[gen][gnodes[gen][node].father].sons[1] != node) {
				if(debug) {
					fprintf(stderr, "\nError: checking genetree for gen %d: ",gen);
					fprintf(stderr, "father %d of node %d has no matching son: %d, %d.", father, node,
							getNodeSon(dataState.lociData[gen], father, 0), getNodeSon(dataState.lociData[gen], father, 1));
				} else {
					fprintf(stderr, "Fatal Error 0068.\n");
				}
				res = 0;
			}
			if(getNodeAge(dataState.lociData[gen], father) <= age) {
				//			if(gnodes[gen][gnodes[gen][node].father].age <= age) {
				if(debug) {
					fprintf(stderr, "\nError: checking genetree for gen %d: ",gen);
					fprintf(stderr, "age of father %d is smaller than previous: %g, %g (diff %g).", father, getNodeAge(dataState.lociData[gen], father), age,getNodeAge(dataState.lociData[gen], father)- age);
				} else {
					fprintf(stderr, "Fatal Error 0069.\n");
				}
				res = 0;
			}
			if(!dataSetup.popTree->pops[ nodePops[gen][father] ]->isAncestralTo[pop]) {
				//			if(!dataSetup.popTree->pops[ gnodes[gen][gnodes[gen][node].father].ipop ]->isAncestralTo[pop]) {
				if(debug) {
					fprintf(stderr, "\nError: checking genetree for gen %d: ",gen);
					fprintf(stderr, "pop %d of father %d is not ancestral to previous pop %d.", nodePops[gen][father], father, pop);
				} else {
					fprintf(stderr, "Fatal Error 0070.\n");
				}
				res = 0;
			}
		} else if(node != getLocusRoot(dataState.lociData[gen])) {
			if(debug) {
				fprintf(stderr, "\nError: checking genetree for gen %d: ",gen);
				fprintf(stderr, "node %d has no father but is not root %d.", node, getLocusRoot(dataState.lociData[gen]));
			} else {
				fprintf(stderr, "Fatal Error 0071.\n");
			}
			res = 0;
		}


	}// end of for(node)

	if(num_migs != genetree_migs[gen].num_migs) {
		if(debug) {
			fprintf(stderr, "\nError: checking genetree for gen %d: ",gen);
			fprintf(stderr, "number of recorded migration nodes doesn't match: %d, %d.", genetree_migs[gen].num_migs, num_migs);
		} else {
			fprintf(stderr, "Fatal Error 0072.\n");
		}
		res = 0;
	}
	for(mig_band=0; mig_band < dataSetup.popTree->numMigBands; mig_band++) {
		num_migs -= genetree_stats[gen].num_migs[mig_band];
	}
	if(num_migs != 0) {
		if(debug) {
			fprintf(stderr, "\nError checking genetree for gen %d: ",gen);
			fprintf(stderr, "number of recorded migration nodes doesn't match sum of recorded per mig band: %d, %d.", genetree_migs[gen].num_migs, genetree_migs[gen].num_migs-num_migs);
		} else {
			fprintf(stderr, "Fatal Error 0073.\n");
		}
		res = 0;
	}



	return res;
}

/* synchronizeEvents
   Synchronizes times recorded by event chains to ones recorded by elemenets (gene tree nodes, migration nodes, population divergences, etc.)
 */

int synchronizeEvents(int gen) {

	int i, pop, event, id;
	int pop_queue[2*NSPECIES-1];
	int res = 1;

	double realAge=0.0, age, PERCISION = 0.0000001;


	populationPostOrder(dataSetup.popTree->rootPop, pop_queue);

	for(i=0; i<dataSetup.popTree->numPops; i++) {
		pop = pop_queue[i];
		event = event_chains[gen].first_event[pop];
		age = dataSetup.popTree->pops[pop]->age;


		// follow event chain and check all saved data
		for( ; event>=0; event = event_chains[gen].events[event].next) {

			id = event_chains[gen].events[event].node_id;
			age += event_chains[gen].events[event].elapsed_time;

			switch(event_chains[gen].events[event].type) {
				case(SAMPLES_START):
																																																			realAge = dataSetup.popTree->pops[pop]->sampleAge;
				break;
				case(COAL):
																																																			realAge = getNodeAge(dataState.lociData[gen],id);
				break;
				case(IN_MIG):
				case(OUT_MIG):
				realAge = genetree_migs[gen].mignodes[id].age;
				break;
				case(MIG_BAND_START):
																																																			realAge = dataSetup.popTree->migBands[id].startTime;
				break;
				case(MIG_BAND_END):
																																																			realAge = dataSetup.popTree->migBands[id].endTime;
				break;
				case(END_CHAIN):
																																																			if(pop != dataSetup.popTree->rootPop) {
																																																				realAge = dataSetup.popTree->pops[pop]->father->age;
																																																			} else {
																																																				realAge = age;
																																																			}
				break;
				default:
					realAge = age;
					break;
			}// end of switch

			if(fabs(realAge - age) > PERCISION) {
				if(debug) {
					fprintf(stderr, "\nError: synchronizing events in genetree for gen %d: ",gen);
					fprintf(stderr, "ellapsed time of event %d (type %d) doesn't match age of element it represents: el %g, ev %g (diff = %g> precision of %g).\n", 
							event, event_chains[gen].events[event].type, realAge, age, realAge - age, PERCISION);
				} else {
					fprintf(stderr, "Fatal Error 0075.\n");
				}
				res = 0;
			}

			event_chains[gen].events[event].elapsed_time += realAge-age;

			if(event_chains[gen].events[event].elapsed_time < -PERCISION) {
				if(debug) {
					fprintf(stderr, "\nError: synchronizing events in genetree for gen %d: ",gen);
					fprintf(stderr, "after correction, ellapsed time of event %d (type %d) becomes negative: ET=%g, diff = %g.\n", 
							event, event_chains[gen].events[event].type, event_chains[gen].events[event].elapsed_time - (realAge-age), realAge - age);
				} else {
					fprintf(stderr, "Fatal Error 0076.\n");
				}
				res = 0;
			} else if(event_chains[gen].events[event].elapsed_time < 0.0) {
				event_chains[gen].events[event].elapsed_time = 0.0;
			}

			age=realAge;

		}// end of for(event)

	}// end of for(pop)

	return res;
}

/* printEventChains
 */

int printEventChains(FILE* stream, int gen) {
	int pop, event, i, mig;

	fprintf(stream, "Event chains for gen %d.\n\nType codes: 0 - COAL, 1- IN_MIG, 2 - OUT_MIG, 3 - MIG_BAND_START, 4 - MIG_BAND_END, 5 - END_CHAIN.\n\n",gen);

	for(pop=0;pop<dataSetup.popTree->numPops;pop++) {
		for(event = event_chains[gen].first_event[pop]; event>=0; event = event_chains[gen].events[event].next) {

			fprintf(stream, "Pop %2d: Event %3d of type %d for node %2d: num-lins = %2d, elapsed-time=%8g, prev = %3d, next=%3d.\n",
					pop, event,
					event_chains[gen].events[event].type,
					event_chains[gen].events[event].node_id,
					event_chains[gen].events[event].num_lineages,
					event_chains[gen].events[event].elapsed_time,
					event_chains[gen].events[event].prev,
					event_chains[gen].events[event].next);

		}
	}

	fprintf(stream, "---------------------------------------------------------------\n");
	if(genetree_migs[gen].num_migs > 0) {
		for(i=0; i< genetree_migs[gen].num_migs; i++) {
			mig = genetree_migs[gen].living_mignodes[i];
			fprintf(stream, "Mignode %d in gen %d: mig band %d, branch %d, events (in,out) %d,%d, age %g.\n",
					mig, gen, genetree_migs[gen].mignodes[mig].migration_band, genetree_migs[gen].mignodes[mig].gtree_branch,
					genetree_migs[gen].mignodes[mig].target_event, genetree_migs[gen].mignodes[mig].source_event, genetree_migs[gen].mignodes[mig].age);
		}
		fprintf(stream, "---------------------------------------------------------------\n");
	}


	return 0;
}

/*************************************  EVENT CHAIN PROCEDURES   ****************************************/
