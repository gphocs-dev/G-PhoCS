#pragma once


#define MAX_MIG_BANDS	40			// max migration bands in the population tree
#define MAX_MIGS		10			// max migration events per genealogy
#define NSPECIES		20			// max # of species
#define NS				200			// max # of sequences
#define OLDAGE			999			// upper bound on age (can be extended...)
#define MAX_EVENTS   (NS + 2*NSPECIES + 3*MAX_MIGS)
#define NUM_DELTA_STATS_INSTANCES 2

#define DEBUG_NODE_CHANGE_NOT
#define DEBUG_RUBBERBAND_NOT

//============================ END OF FILE ====================================
