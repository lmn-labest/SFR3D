#include<HccaTime.h>

double getTimeC(void){
/**********************************************************************
 * GETTIMEC : Obter tempo                                             *
 * -------------------------------------------------------------------*
 * Parametros de Entrada:                                             *
 * -------------------------------------------------------------------*
 * -------------------------------------------------------------------*
 * Parametros de Saida:                                               *
 * -------------------------------------------------------------------*
 * -------------------------------------------------------------------*
*********************************************************************/

/*funcao MPI */
#if _MPICH_
  double t=0.0;
  t = MPI_Wtime();
  return t;
/*funcao do openmp*/
#elif _OPENMP 
  #if _MPICH_
    double t=0.0;
    t = MPI_Wtime();
    return t;
  #else
    double t=0.0;
    t = omp_get_wtime();
    return t;
  #endif
/*funcao do nativa do C (baixa precisao)*/
#elif _WIN32
   clock_t start;
   double  t;

   start = clock();
   t = (double)(start) / CLOCKS_PER_SEC;
   return t;
#else
  double t =0.0;
  struct timeval start;

  gettimeofday(&start,NULL);
  t = (double) (start.tv_sec + start.tv_usec/1000000.e0);
  
  return t;
#endif

}

/*********************************************************************/
