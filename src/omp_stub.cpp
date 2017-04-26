#ifndef ENABLE_OMP_THREADS
//_OPENMP
extern "C"
{
  void omp_set_dynamic(int n){}
  void omp_set_num_threads(int n) {}
  int omp_get_num_threads(){return 1;}
  int omp_get_thread_num(){return 0;}
  int omp_get_max_threads(){return 1;}
}
#endif



