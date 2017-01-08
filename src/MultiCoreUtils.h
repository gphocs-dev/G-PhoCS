#ifndef MultiCoreUtils
#define MultiCoreUtils

#ifdef ENABLE_OMP_THREADS

	#include <omp.h>

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

#else
	extern void omp_set_num_threads(int n);
	extern int omp_get_max_threads();
	extern int omp_get_thread_num();
#endif


//#define RECORD_METHOD_TIMES








#endif
