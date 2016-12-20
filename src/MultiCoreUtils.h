#ifndef MultiCoreUtils
#define MultiCoreUtils

#include <omp.h>


#define THREAD_COUNT_GPHOCS  100 // This indicates the maximum thread count possible. System will always select MIN(Envirnoment setting , THREAD_COUNT_GPHOCS)
#define THREAD_MT_ON //Comment out to turn off the MT
#define THREAD_SCHEDULING_STRATEGY static //use static, dynamic or guided. see openmp help for difference. you can also set: dynamic,100 if you wish to specify the chunk size. DEFAULT = static

/* flags to disable or enable MT on specific methods
 * comment out a DEFINE to disable.
 */

#define THREAD_UpdateGB_InternalNode
#define THREAD_UpdateGB_MigrationNode
#define THREAD_UpdateGB_MigSPR
#define THREAD_UpdateTau
#define THREAD_UpdateMigRates
#define THREAD_mixing
//#define THREAD_UpdateTheta
#define THREAD_UpdateSampleAge



//#define RECORD_METHOD_TIMES








#endif
